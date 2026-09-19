#pragma once

#include "active.hpp"
#include <type_traits>

namespace terminal {
using namespace allsize;

struct BuildWork {
    uint64_t states=0, degree_tests=0, minimum_writes=0;
    uint64_t full_groups=0, partial_groups=0, choices=0, saved_members=0;
    size_t minimum_scratch_bytes=0;
};

using Offset = uint64_t;   // position in `members`; a graph's clique tree can hold more than 2^32 incidences (com-lj, hollywood)
using RowId = uint64_t;    // row (leaf) index; com-lj has more than 2^30 rows.  Reverse codes pack (row << 2) | role in 64 bits.
struct Row {
    Offset begin, hold_end, pivot_end, end;
    Vertex group, lo, hi;
    Vertex holds() const { return static_cast<Vertex>(hold_end-begin); }
    Vertex pivots() const { return static_cast<Vertex>(pivot_end-hold_end); }
    bool valid(int s) const { return lo<=static_cast<Vertex>(s) && static_cast<Vertex>(s)<=hi; }
};

struct Index {
    std::vector<Row> rows;
    Vertices members;
    std::vector<uint64_t> reverse;        // per vertex: (row << 2) | role, role 0 hold, 1 pivot, 2 choice
    std::vector<RowId> group_row;
    std::vector<size_t> reverse_off;
    std::vector<uint8_t> zero_choice;
    BuildWork work;
    int maximum;
    explicit Index(int s): maximum(s) {}

    void append(const Vertices& h,const Vertices& q,const Vertices& x={},bool zero=false) {
        require(!h.empty(),"terminal requires a hold");
        require(group_row.size()<absent,"group ID overflow");
        const Offset begin=members.size(),mid=begin+h.size(),last=mid+q.size();
        const Vertex low=std::max<size_t>(2,h.size()+(!x.empty() && !zero));
        const Vertex high=std::min<size_t>(maximum,h.size()+q.size()+!x.empty());
        require(low<=high,"invalid terminal size interval");
        Vertex group=absent;
        if(!x.empty()) {
            require(x.size()>=2,"singleton group");
            group=group_row.size();group_row.push_back(rows.size());zero_choice.push_back(zero);
            work.choices+=x.size();work.saved_members+=(x.size()-1)*(h.size()+q.size());
            if(zero)++work.full_groups;else ++work.partial_groups;
        }
        members.insert(members.end(),h.begin(),h.end());
        members.insert(members.end(),q.begin(),q.end());
        members.insert(members.end(),x.begin(),x.end());
        rows.push_back({begin,mid,last,static_cast<Offset>(members.size()),group,low,high});
    }
    void prepare(Vertex n) {
        reverse_off.assign(static_cast<size_t>(n)+1,0);
        for(Vertex v:members)++reverse_off[v+1];
        std::partial_sum(reverse_off.begin(),reverse_off.end(),reverse_off.begin());
        reverse.resize(members.size());
        auto cursor=reverse_off;
        for(RowId p=0;p<rows.size();++p) {
            const auto& row=rows[p];
            for(Offset i=row.begin;i<row.end;++i) {
                const uint64_t role=i<row.hold_end?0:(i<row.pivot_end?1:2);
                reverse[cursor[members[i]]++]=(p<<2)|role;
            }
        }
    }
    std::span<const uint64_t> touching(Vertex v) const {
        return std::span<const uint64_t>(reverse).subspan(reverse_off[v],reverse_off[v+1]-reverse_off[v]);
    }
    size_t bytes() const {
        return rows.capacity()*sizeof(Row)
            +members.capacity()*sizeof(Vertex)+(reverse.capacity()+group_row.capacity())*sizeof(uint64_t)
            +reverse_off.capacity()*sizeof(size_t)+zero_choice.capacity();
    }
};

// Mode 0 includes the old direct-empty-child shortcut, but no merging.
template<int Mode> class Builder {
    const Graph& graph_;
    Index& index_;
    Vertices holds_,pivots_,minimum_;
    void visit(Vertices candidates) {
        ++index_.work.states;
        const size_t maximum=index_.maximum;
        if(holds_.size()>maximum || holds_.size()+pivots_.size()+candidates.size()<2)return;
        if(holds_.size()==maximum){index_.append(holds_,{});return;}
        auto close=[&] {
            const size_t saved=pivots_.size();
            pivots_.insert(pivots_.end(),candidates.begin(),candidates.end());
            index_.append(holds_,pivots_);pivots_.resize(saved);
        };
        if(candidates.empty() || holds_.size()+1==maximum){close();return;}
        Vertices universal;
        Vertex pivot=absent;
        size_t best=0,lowest=std::numeric_limits<size_t>::max();
        if constexpr(Mode==2)minimum_.clear();
        for(Vertex u:candidates) {
            ++index_.work.degree_tests;
            const size_t degree=intersection_size(candidates,graph_.row(u));
            if(degree+1==candidates.size())universal.push_back(u);
            else {
                if(pivot==absent || degree>best){pivot=u;best=degree;}
                if constexpr(Mode==2) {
                    if(degree<lowest){lowest=degree;minimum_.clear();}
                    if(degree==lowest){minimum_.push_back(u);++index_.work.minimum_writes;}
                }
            }
        }
        if constexpr(Mode==2)index_.work.minimum_scratch_bytes=
            std::max(index_.work.minimum_scratch_bytes,minimum_.capacity()*sizeof(Vertex));
        if(universal.size()==candidates.size()){close();return;}
        const size_t saved=pivots_.size();
        Vertices remaining;
        std::set_difference(candidates.begin(),candidates.end(),universal.begin(),universal.end(),
                            std::back_inserter(remaining));
        candidates.swap(remaining);
        pivots_.insert(pivots_.end(),universal.begin(),universal.end());
        if(best==universal.size()) {
            if constexpr(Mode>0)index_.append(holds_,pivots_,candidates,true);
            else {
                pivots_.push_back(pivot);index_.append(holds_,pivots_);pivots_.pop_back();
                for(Vertex u:candidates)if(u!=pivot) {
                    holds_.push_back(u);index_.append(holds_,pivots_);holds_.pop_back();
                }
            }
            pivots_.resize(saved);return;
        }
        if constexpr(Mode==2) {
            if(lowest==universal.size() && minimum_.size()>=2) {
                index_.append(holds_,pivots_,minimum_,false);
                remaining.clear();
                std::set_difference(candidates.begin(),candidates.end(),minimum_.begin(),minimum_.end(),
                                    std::back_inserter(remaining));
                candidates.swap(remaining);
            }
        }
        Vertices branches{pivot};
        for(Vertex u:candidates)if(u!=pivot && !graph_.adjacent(pivot,u))branches.push_back(u);
        for(Vertex u:branches) {
            const auto at=std::lower_bound(candidates.begin(),candidates.end(),u);
            require(at!=candidates.end() && *at==u,"missing terminal branch");
            candidates.erase(at);
            auto& selected=u==pivot?pivots_:holds_;
            selected.push_back(u);visit(intersect(candidates,graph_.row(u)));selected.pop_back();
        }
        pivots_.resize(saved);
    }
public:
    Builder(const Graph& graph,Index& index): graph_(graph),index_(index) {}
    void run() {
        for(Vertex v=0;v<graph_.n;++v) {
            Vertices later;
            for(Vertex u:graph_.row(v))if(u>v)later.push_back(u);
            holds_.assign(1,v);pivots_.clear();visit(std::move(later));
        }
    }
};

inline void build(const Graph& graph,Index& index,int mode) {
    if(mode==0)Builder<0>(graph,index).run();
    else if(mode==1)Builder<1>(graph,index).run();
    else {require(mode==2,"unknown builder mode");Builder<2>(graph,index).run();}
}

struct Metrics {
    uint64_t positive_writes=0,affected_rows=0,group_updates=0;
    uint64_t lazy_discards=0,lazy_writes=0,coefficient_calls=0;
};

#include "terminal_replay.inc"

template<class Count> struct Solver {
    using Base=fullrange::Kernel<Count>;
    using Choose=typename Base::Combinations;
    using Queue=typename Base::Queue;
    struct Value {Count h=0,q=0,x=0;};
    struct Result {typename Base::Common common;orderdp::Extra extra;Metrics metrics;orderreplay::Metrics replay;};

    static Count multiply(Count a,Vertex b) {
        const typename Base::Wider value=typename Base::Wider(a)*b;
        if(value>=Base::infinity)throw std::overflow_error("factored product overflow");
        return static_cast<Count>(value);
    }
    static Value coefficients(const Index& index,RowId p,int s,Vertex q,Vertex z,const Choose& choose,
                              Metrics& metrics) {
        const Row& row=index.rows[p];const int r=s-static_cast<int>(row.holds());
        Value value;
        if(row.group==absent) {
            value.h=choose(q,r);value.q=q?choose(static_cast<int>(q)-1,r-1):Count{0};
            metrics.coefficient_calls+=2;return value;
        }
        value.x=choose(q,r-1);value.h=multiply(value.x,z);
        value.q=q?multiply(choose(static_cast<int>(q)-1,r-2),z):Count{0};
        metrics.coefficient_calls+=2;
        if(index.zero_choice[row.group]) {
            Base::checked_add(value.h,choose(q,r));
            if(q)Base::checked_add(value.q,choose(static_cast<int>(q)-1,r-1));
            metrics.coefficient_calls+=2;
        }
        return value;
    }

    // Sink: optional row consumer `sink(int s, std::span<const Count> row)`.  With a sink the core matrix is a two-row
    // window (rows s-1 and s), every finished row is delivered once, in increasing s, and result.common.data.core holds
    // only the window afterwards.  Without a sink (default) the full matrix is returned as before.
    template<bool Audit=false,bool Replay=false,class Sink=std::nullptr_t> static Result solve(const Graph& graph,const Index& index,const Choose& choose,
                  const Vertices& ordinary,const std::vector<Count>* expected=nullptr,Sink sink=nullptr) {
        constexpr bool Stream=!std::is_same_v<Sink,std::nullptr_t>;
        const auto start=Clock::now();const Vertex n=graph.n;
        Result result;auto& data=result.common.data;auto& stats=result.common.stats;
        auto& work=data.work;auto& extra=result.extra;auto& metrics=result.metrics;
        ReplayPlan<Count> plan(n,index,choose,Replay,result.replay);
        auto slot=[&](int s){return static_cast<size_t>(Stream?(s&1):s)*n;};
        auto emit=[&](int s,std::span<const Count> row){if constexpr(Stream)sink(s,row);};
        data.core.assign((Stream?size_t{2}:static_cast<size_t>(index.maximum+1))*n,0);
        std::copy(ordinary.begin(),ordinary.end(),data.core.begin()+slot(2));
        emit(2,std::span<const Count>(data.core).subspan(slot(2),n));
        Vertices order(n),next_order;std::iota(order.begin(),order.end(),0);
        if(!std::is_sorted(ordinary.begin(),ordinary.end()))
            std::sort(order.begin(),order.end(),[&](Vertex a,Vertex b){return ordinary[a]<ordinary[b] || (ordinary[a]==ordinary[b] && a<b);});
        next_order.reserve(n);
        const Vertex maximum=ordinary.empty()?0:*std::max_element(ordinary.begin(),ordinary.end());
        struct Cache {Count key=Base::infinity,value=0;};
        std::vector<Cache> cache(4096);
        for(int s=3;s<=index.maximum;++s) {
            auto upper=std::span<Count>(data.core).subspan(slot(s),n);
            const auto previous=std::span<const Count>(data.core).subspan(slot(s-1),n);
            if constexpr(Stream)std::fill(upper.begin(),upper.end(),0);   // the slot held row s-2
            if constexpr(Replay) {
                const bool accepted=plan.trial(s,order,upper);
                data.state_bytes=std::max(data.state_bytes,plan.bytes()+(order.capacity()+next_order.capacity())*sizeof(Vertex)+cache.capacity()*sizeof(Cache));
                if constexpr(Audit) {
                    require(expected,"replay needs audit reference");
                    for(Vertex v=0;v<n;++v) {
                        require(upper[v]>=(*expected)[static_cast<size_t>(s)*n+v],"replay upper below true core");
                        if(accepted)require(upper[v]==(*expected)[static_cast<size_t>(s)*n+v],"false factored replay acceptance");
                    }
                }
                if(accepted) {
                    emit(s,upper);
                    if(std::all_of(upper.begin(),upper.end(),[](const Count& value){return value==0;}))break;
                    continue;
                }
                std::fill(upper.begin(),upper.end(),0);
            }
            std::vector<Count> support(n,0),wh(index.rows.size(),0),wp(index.rows.size(),0),wx(index.group_row.size(),0);
            Vertices count(index.rows.size(),0),choices(index.group_row.size(),0);
            std::vector<Offset> choice_off(index.group_row.size(),0);Vertices choice_size(index.group_row.size(),0),scratch;
            std::vector<uint8_t> live(n,0),touched(index.rows.size(),0),dead(index.rows.size(),1),dirty(n,0);
            for(RowId p=0;p<index.rows.size();++p) {
                const auto& row=index.rows[p];if(!row.valid(s))continue;
                dead[p]=0;count[p]=row.pivots();
                const Vertex z=row.end-row.pivot_end;
                const auto value=coefficients(index,p,s,count[p],z,choose,metrics);
                wh[p]=value.h;wp[p]=value.q;
                auto initialize=[&](Offset begin,Offset end,Count weight) {
                    if(!weight)return;
                    work.count_reads+=end-begin;
                    for(Offset i=begin;i<end;++i)Base::checked_add(support[index.members[i]],weight);
                };
                initialize(row.begin,row.hold_end,value.h);initialize(row.hold_end,row.pivot_end,value.q);
                if(row.group!=absent) {
                    const Vertex g=row.group;wx[g]=value.x;choices[g]=z;
                    initialize(row.pivot_end,row.end,value.x);
                    choice_off[g]=scratch.size();
                    if(value.x) {
                        choice_size[g]=z;
                        scratch.insert(scratch.end(),index.members.begin()+row.pivot_end,index.members.begin()+row.end);
                    }
                }
            }
            Vertex remaining=0;next_order.clear();std::fill(cache.begin(),cache.end(),Cache{});
            const auto bounds=Clock::now();
            for(Vertex v=0;v<n;++v) {
                if(!support[v]){next_order.push_back(v);continue;}
                live[v]=1;++remaining;const Count a=previous[v];
                auto& entry=cache[static_cast<size_t>((a^(a>>17)^(a>>37))&(cache.size()-1))];
                if(entry.key!=a){entry={a,Base::integer_upper(a,s-2,maximum,stats)};++stats.cache_misses;}
                else ++stats.cache_hits;
                upper[v]=entry.value;stats.initial_capped+=upper[v]<support[v];
                require(upper[v]>0,"zero bound for a clique member");
            }
            stats.bounds_ms+=ms(bounds);
            Queue heap(support,upper,false,extra);
            Vertices batch,changed;std::vector<RowId> affected;size_t cursor=0;
            auto stream_key=[&]() -> Count {
                while(cursor<order.size() && !live[order[cursor]]){++cursor;++extra.order_reads;}
                return cursor<order.size()?upper[order[cursor]]:Base::infinity;
            };
            std::vector<uint64_t> cliques;
            if constexpr(Audit)cliques=bottomup::clique_masks(graph,s);
            auto audit=[&] {
                if constexpr(Audit) {
                    require(n<=63 && expected,"missing small-graph audit reference");
                    uint64_t mask=0;for(Vertex v=0;v<n;++v)if(live[v])mask|=uint64_t{1}<<v;
                    std::vector<uint64_t> actual(n,0);
                    for(auto clique:cliques)if((clique&mask)==clique)
                        for(Vertex v=0;v<n;++v)if((clique>>v)&1)++actual[v];
                    Count minimum=Base::infinity,last=0;
                    for(Vertex v:order)if(live[v]) {
                        require(support[v]==actual[v],"factored residual degree mismatch");
                        require(upper[v]>=last,"unsorted inherited upper stream");last=upper[v];
                        require(heap.contains(v)==(support[v]<upper[v]),"factored queue membership mismatch");
                        if(heap.contains(v))require(heap.key(v)==std::min(support[v],upper[v]),"factored queue key mismatch");
                        minimum=std::min(minimum,std::min(upper[v],support[v]));
                    }
                    for(Vertex v=0;v<n;++v)if(!live[v])require(upper[v]==(*expected)[static_cast<size_t>(s)*n+v],"factored completed label mismatch");
                    require(minimum==std::min(stream_key(),heap.first_key()),"factored minimum mismatch");heap.audit();
                }
            };
            auto memory=[&] {
                data.state_bytes=std::max(data.state_bytes,
                    (support.capacity()+wh.capacity()+wp.capacity()+wx.capacity())*sizeof(Count)
                    +(count.capacity()+choices.capacity()+2*choice_off.capacity()+choice_size.capacity()+scratch.capacity()
                      +order.capacity()+next_order.capacity()+batch.capacity()+affected.capacity()+changed.capacity())*sizeof(Vertex)
                    +live.capacity()+touched.capacity()+dead.capacity()+dirty.capacity()+heap.bytes()+cache.capacity()*sizeof(Cache)+plan.bytes());
            };
            audit();memory();
            if(!remaining){extra.skipped_layers=index.maximum-s;emit(s,upper);break;}
            Count level=0;
            while(remaining) {
                ++stats.batches;level=std::max(level,std::min(stream_key(),heap.first_key()));
                require(level!=Base::infinity,"missing factored minimum");batch.clear();
                while(std::min(stream_key(),heap.first_key())<=level) {
                    Vertex v;
                    if(heap.first_key()<=level)v=heap.pop();
                    else {v=order[cursor++];++extra.order_reads;++extra.implicit_pops;require(!heap.contains(v),"implicit heap duplicate");}
                    require(live[v],"duplicate factored removal");live[v]=0;--remaining;
                    stats.early_removals+=support[v]>level;upper[v]=level;next_order.push_back(v);batch.push_back(v);++work.events;
                }
                if(!remaining){audit();memory();break;}
                affected.clear();changed.clear();
                for(Vertex v:batch)for(uint64_t code:index.touching(v)) {
                    ++work.source_reads;const RowId p=code>>2;const unsigned role=code&3;
                    if(dead[p])continue;
                    if(!touched[p]){touched[p]=1;affected.push_back(p);}
                    if(role==0)dead[p]=1;
                    else if(role==1){require(count[p]>0,"pivot counter underflow");--count[p];}
                    else {auto& z=choices[index.rows[p].group];require(z>0,"choice counter underflow");--z;}
                }
                auto subtract=[&](Vertex v,Count loss) {
                    if(!live[v])return;
                    require(support[v]>=loss,"factored degree underflow");support[v]-=loss;++metrics.positive_writes;
                    if(!dirty[v]){dirty[v]=1;changed.push_back(v);}
                };
                auto scan=[&](Offset begin,Offset end,Count loss) {
                    if(!loss)return;
                    work.target_reads+=end-begin;
                    for(Offset i=begin;i<end;++i)subtract(index.members[i],loss);
                };
                for(RowId p:affected) {
                    ++metrics.affected_rows;const auto& row=index.rows[p];const Vertex g=row.group;
                    const Value value=dead[p]?Value{}:coefficients(index,p,s,count[p],g==absent?0:choices[g],choose,metrics);
                    require(value.h<=wh[p] && value.q<=wp[p],"negative factored common loss");
                    const Count lh=wh[p]-value.h,lq=wp[p]-value.q;
                    wh[p]=value.h;wp[p]=value.q;touched[p]=0;if(!value.h)dead[p]=1;
                    scan(row.begin,row.hold_end,lh);scan(row.hold_end,row.pivot_end,lq);
                    if(g==absent)continue;
                    ++metrics.group_updates;require(value.x<=wx[g],"negative choice loss");
                    const Count lx=wx[g]-value.x;wx[g]=value.x;
                    if(!lx)continue;
                    const Offset begin=choice_off[g];Vertex length=choice_size[g],at=0;
                    while(at<length) {
                        Vertex v=scratch[begin+at];++work.target_reads;bool replacement=false;
                        while(!live[v]) {
                            ++metrics.lazy_discards;--length;
                            if(at==length)break;
                            v=scratch[begin+length];++work.target_reads;replacement=true;
                        }
                        if(at==length)break;
                        if(replacement){scratch[begin+at]=v;++metrics.lazy_writes;}
                        subtract(v,lx);++at;
                    }
                    choice_size[g]=length;
                }
                for(Vertex v:changed) {
                    dirty[v]=0;const Count next=std::min(upper[v],support[v]);
                    if(heap.contains(v)) {
                        if(next<heap.key(v)){heap.decrease(v,next);++work.updates;}
                        else ++stats.unchanged_keys;
                    }else if(support[v]<upper[v]){heap.insert(v,next);++work.updates;}
                    else ++stats.unchanged_keys;
                }
                audit();memory();
            }
            require(next_order.size()==n,"incomplete factored order");order.swap(next_order);plan.invalidate();
            emit(s,upper);
        }
        data.peel_ms=ms(start);return result;
    }
};
}
