template<class Count> class ReplayPlan {
    using Base=fullrange::Kernel<Count>;
    using Choose=typename Base::Combinations;
    const Index& index_;
    const Choose& choose_;
    orderreplay::Metrics& metrics_;
    std::vector<Count> value_;
    Vertices top_,ranks_,ordered_,first_,before_,active_rows_,active_order_;
    bool valid_=false;

    Count coefficient(int q,int r) {++metrics_.coefficient_calls;return choose_(q,r);}
    Count product(Count value,Vertex z) {
        const typename Base::Wider full=typename Base::Wider(value)*z;
        if(full>=Base::infinity)throw std::overflow_error("replay group product overflow");
        return static_cast<Count>(full);
    }
    void add(Vertex v,Count value) {
        if(value){Base::checked_add(value_[v],value);++metrics_.nonzero_additions;}
    }
    Count common(Vertex p,Vertex q,Vertex z,int r) {
        Count value=product(coefficient(q,r-1),z);
        if(index_.zero_choice[index_.rows[p].group])Base::checked_add(value,coefficient(q,r));
        return value;
    }
    void memory() {
        metrics_.scratch_bytes=value_.capacity()*sizeof(Count)
            +(top_.capacity()+ranks_.capacity()+ordered_.capacity()+first_.capacity()+before_.capacity()
              +active_rows_.capacity()+active_order_.capacity())*sizeof(Vertex);
    }
    void compile(const Vertices& order) {
        const auto begin=Clock::now();
        for(Vertex i=0;i<order.size();++i)ranks_[order[i]]=i;
        ordered_.resize(index_.members.size());first_.resize(index_.rows.size());before_.resize(index_.rows.size());
        for(Vertex p:active_rows_) {
            const auto& row=index_.rows[p];Vertex first=index_.members[row.begin];
            for(Offset i=row.begin;i<row.hold_end;++i) {
                ++metrics_.compile_reads;
                if(ranks_[index_.members[i]]<ranks_[first])first=index_.members[i];
            }
            first_[p]=first;
            std::copy(index_.members.begin()+row.hold_end,index_.members.begin()+row.end,ordered_.begin()+row.hold_end);
            auto less=[&](Vertex a,Vertex b){++metrics_.sort_comparisons;return ranks_[a]<ranks_[b];};
            std::sort(ordered_.begin()+row.hold_end,ordered_.begin()+row.pivot_end,less);
            std::sort(ordered_.begin()+row.pivot_end,ordered_.begin()+row.end,less);
            Offset at=row.hold_end;
            while(at<row.pivot_end && ranks_[ordered_[at]]<ranks_[first]){++at;++metrics_.compile_reads;}
            before_[p]=static_cast<Vertex>(at-row.hold_end);
        }
        valid_=true;++metrics_.compilations;metrics_.compile_ms+=ms(begin);memory();
    }
    void forward_counts(int s) {
        for(Vertex p:active_rows_) {
            const auto& row=index_.rows[p];if(!row.valid(s))continue;
            const int r=s-static_cast<int>(row.holds());
            if(row.group==absent) {
                add(first_[p],coefficient(row.pivots()-before_[p],r));++metrics_.forward_reads;
                for(Vertex j=0;j<before_[p];++j) {
                    add(ordered_[row.hold_end+j],coefficient(static_cast<int>(row.pivots()-j)-1,r-1));
                    ++metrics_.forward_reads;
                }
                continue;
            }
            Offset qi=row.hold_end,xi=row.pivot_end;
            const Vertex limit=ranks_[first_[p]];
            while(true) {
                const Vertex qr=qi<row.pivot_end?ranks_[ordered_[qi]]:absent;
                const Vertex xr=xi<row.end?ranks_[ordered_[xi]]:absent;
                if(std::min(qr,xr)>=limit)break;
                if(qr<xr) {
                    const Vertex v=ordered_[qi++];
                    Count value=product(coefficient(row.pivot_end-qi,r-2),row.end-xi);
                    if(index_.zero_choice[row.group])Base::checked_add(value,coefficient(row.pivot_end-qi,r-1));
                    add(v,value);
                }else {
                    const Vertex v=ordered_[xi++];add(v,coefficient(row.pivot_end-qi,r-1));
                }
                ++metrics_.forward_reads;
            }
            add(first_[p],common(p,row.pivot_end-qi,row.end-xi,r));++metrics_.forward_reads;
        }
    }
    void certificate(int s,std::span<const Count> upper) {
        for(Vertex p:active_rows_) {
            const auto& row=index_.rows[p];if(!row.valid(s))continue;
            const int r=s-static_cast<int>(row.holds());const Count level=upper[first_[p]];
            if(row.group==absent) {
                Offset at=row.hold_end,at_level=row.pivot_end;
                while(at<row.pivot_end) {
                    const Count key=upper[ordered_[at]];
                    if(key>=level && at_level==row.pivot_end)at_level=at;
                    if(key>level)break;
                    Offset next=at+1;
                    while(next<row.pivot_end && upper[ordered_[next]]==key)++next;
                    const Count weight=key?coefficient(static_cast<int>(row.pivot_end-at)-1,r-1):Count{0};
                    for(Vertex i=at;i<next;++i){add(ordered_[i],weight);++metrics_.certificate_reads;}
                    at=next;
                }
                const Count weight=level?coefficient(row.pivot_end-at_level,r):Count{0};
                for(Offset i=row.begin;i<row.hold_end;++i) {
                    const Vertex v=index_.members[i];if(upper[v]==level)add(v,weight);++metrics_.certificate_reads;
                }
                continue;
            }
            Offset qi=row.hold_end,xi=row.pivot_end;
            auto next_key=[&]() -> Count {
                return std::min(qi<row.pivot_end?upper[ordered_[qi]]:Base::infinity,
                                xi<row.end?upper[ordered_[xi]]:Base::infinity);
            };
            auto group=[&](Count key) {
                Offset qe=qi,xe=xi;
                while(qe<row.pivot_end && upper[ordered_[qe]]==key)++qe;
                while(xe<row.end && upper[ordered_[xe]]==key)++xe;
                if(qe>qi) {
                    Count weight=0;
                    if(key) {
                        const int q=static_cast<int>(row.pivot_end-qi)-1;
                        weight=product(coefficient(q,r-2),row.end-xi);
                        if(index_.zero_choice[row.group])Base::checked_add(weight,coefficient(q,r-1));
                    }
                    for(Vertex i=qi;i<qe;++i){add(ordered_[i],weight);++metrics_.certificate_reads;}
                }
                if(xe>xi) {
                    const Count weight=key?coefficient(row.pivot_end-qi,r-1):Count{0};
                    for(Vertex i=xi;i<xe;++i){add(ordered_[i],weight);++metrics_.certificate_reads;}
                }
                qi=qe;xi=xe;
            };
            while(next_key()<level)group(next_key());
            const Count weight=level?common(p,row.pivot_end-qi,row.end-xi,r):Count{0};
            for(Offset i=row.begin;i<row.hold_end;++i) {
                const Vertex v=index_.members[i];if(upper[v]==level)add(v,weight);++metrics_.certificate_reads;
            }
            if(next_key()==level)group(level);
        }
    }
public:
    ReplayPlan(Vertex n,const Index& index,const Choose& choose,bool enabled,orderreplay::Metrics& metrics)
        :index_(index),choose_(choose),metrics_(metrics) {
        if(enabled) {
            value_.resize(n);top_.assign(n,1);ranks_.resize(n);active_order_.reserve(n);
            active_rows_.resize(index.rows.size());std::iota(active_rows_.begin(),active_rows_.end(),0);
            for(const auto& row:index.rows)for(Offset i=row.begin;i<row.end;++i)
                top_[index.members[i]]=std::max(top_[index.members[i]],row.hi);
        }
        memory();
    }
    void invalidate(){valid_=false;}
    size_t bytes() const{return metrics_.scratch_bytes;}
    const std::vector<Count>& counts() const{return value_;}
    void forward(int s,const Vertices& order,std::span<Count> upper) {
        std::erase_if(active_rows_,[&](Vertex p){return index_.rows[p].hi<static_cast<Vertex>(s);});
        if(!valid_) {
            active_order_.clear();std::fill(value_.begin(),value_.end(),0);
            for(Vertex v:order)if(top_[v]>=static_cast<Vertex>(s))active_order_.push_back(v);
            compile(active_order_);
        }else std::erase_if(active_order_,[&](Vertex v) {
            if(top_[v]>=static_cast<Vertex>(s))return false;
            value_[v]=0;return true;
        });
        for(Vertex v:active_order_)value_[v]=0;
        forward_counts(s);std::fill(upper.begin(),upper.end(),0);Count level=0;
        for(Vertex v:active_order_){level=std::max(level,value_[v]);upper[v]=level;}
    }
    bool trial(int s,const Vertices& order,std::span<Count> upper) {
        const auto start=Clock::now();const uint64_t compilations=metrics_.compilations;
        forward(s,order,upper);for(Vertex v:active_order_)value_[v]=0;certificate(s,upper);
        orderreplay::Layer layer;layer.s=s;layer.compiled=metrics_.compilations!=compilations;
        for(Vertex v:active_order_){layer.failures+=value_[v]<upper[v];layer.positive+=upper[v]!=0;}
        layer.accepted=layer.failures==0;
        if(layer.accepted)++metrics_.accepted;else ++metrics_.rejected;
        metrics_.layers.push_back(layer);memory();metrics_.replay_ms+=ms(start);return layer.accepted;
    }
};
