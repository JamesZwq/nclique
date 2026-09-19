#pragma once

#include "generated.hpp"

namespace orderreplay {
using namespace fullrange;

struct Layer {
    int s=0;
    uint64_t failures=0,positive=0;
    bool accepted=false,compiled=false;
};

struct Metrics {
    uint64_t forward_reads=0,certificate_reads=0,counter_updates=0;
    uint64_t coefficient_calls=0,nonzero_additions=0,compile_reads=0,sort_comparisons=0;
    uint64_t compilations=0,accepted=0,rejected=0;
    size_t scratch_bytes=0;
    double compile_ms=0,replay_ms=0;
    std::vector<Layer> layers;
};

template<class Count> class Plan {
    using Base=Kernel<Count>;
    using Choose=typename Base::Combinations;
    const Graph& graph_;
    const Layout& index_;
    const Choose& choose_;
    int mode_;
    Metrics& work_;
    std::vector<Count> value_;
    Vertices ranks_,ordered_,first_,before_,count_;
    std::vector<uint8_t> alive_;
    bool valid_=false;

    Count coefficient(int q,int r) {
        ++work_.coefficient_calls;
        return choose_(q,r);
    }
    void add(Vertex v,Count value) {
        if(value) {Base::checked_add(value_[v],value);++work_.nonzero_additions;}
    }
    void memory() {
        work_.scratch_bytes=value_.capacity()*sizeof(Count)
            +(ranks_.capacity()+ordered_.capacity()+first_.capacity()+before_.capacity()+count_.capacity())*sizeof(Vertex)
            +alive_.capacity();
    }
    void compile(const Vertices& order) {
        const auto start=Clock::now();
        ranks_.resize(graph_.n);
        for(Vertex i=0;i<graph_.n;++i)ranks_[order[i]]=i;
        ordered_=index_.paths.vertices;
        first_.resize(index_.paths.size());before_.resize(index_.paths.size());
        for(size_t p=0;p<index_.paths.size();++p) {
            const size_t begin=index_.paths.off[p],middle=begin+index_.paths.holds[p],end=index_.paths.off[p+1];
            require(begin<middle,"empty mandatory set");
            Vertex first=ordered_[begin];
            for(size_t i=begin;i<middle;++i) {
                ++work_.compile_reads;
                if(ranks_[ordered_[i]]<ranks_[first])first=ordered_[i];
            }
            std::sort(ordered_.begin()+middle,ordered_.begin()+end,[&](Vertex u,Vertex v) {
                ++work_.sort_comparisons;return ranks_[u]<ranks_[v];
            });
            size_t at=middle;
            while(at<end && ranks_[ordered_[at]]<ranks_[first]){++at;++work_.compile_reads;}
            first_[p]=first;before_[p]=at-middle;
        }
        valid_=true;++work_.compilations;work_.compile_ms+=ms(start);memory();
    }
    void compiled_forward(int s) {
        for(size_t p=0;p<index_.paths.size();++p)if(index_.valid(p,s)) {
            const int h=index_.paths.holds[p],q=index_.paths.row(p).size()-h,r=s-h;
            ++work_.forward_reads;
            add(first_[p],coefficient(q-before_[p],r));
            const size_t middle=index_.paths.off[p]+h;
            for(Vertex j=0;j<before_[p];++j) {
                ++work_.forward_reads;
                add(ordered_[middle+j],coefficient(q-j-1,r-1));
            }
        }
    }
    void compiled_certificate(int s,std::span<const Count> upper) {
        for(size_t p=0;p<index_.paths.size();++p)if(index_.valid(p,s)) {
            const size_t begin=index_.paths.off[p],middle=begin+index_.paths.holds[p],end=index_.paths.off[p+1];
            const int r=s-static_cast<int>(index_.paths.holds[p]);
            const Count level=upper[first_[p]];
            size_t at=middle,at_level=end;
            while(at<end) {
                const Count key=upper[ordered_[at]];
                if(key>=level && at_level==end)at_level=at;
                if(key>level)break;
                size_t next=at+1;
                while(next<end && upper[ordered_[next]]==key)++next;
                const Count weight=key ? coefficient(end-at-1,r-1) : Count{0};
                for(size_t i=at;i<next;++i){++work_.certificate_reads;add(ordered_[i],weight);}
                at=next;
            }
            const Count weight=level ? coefficient(end-at_level,r) : Count{0};
            for(size_t i=begin;i<middle;++i) {
                ++work_.certificate_reads;
                if(upper[ordered_[i]]==level)add(ordered_[i],weight);
            }
        }
    }
    void reset_counters(int s) {
        count_.resize(index_.paths.size());alive_.resize(index_.paths.size());
        for(size_t p=0;p<index_.paths.size();++p) {
            count_[p]=index_.paths.row(p).size()-index_.paths.holds[p];
            alive_[p]=index_.valid(p,s);
        }
    }
    void erase(Vertex v) {
        for(Vertex occurrence:index_.touching(v)) {
            ++work_.counter_updates;const size_t p=index_.owner[occurrence];
            if(!alive_[p])continue;
            if(index_.pivot(occurrence)){require(count_[p]>0,"replay counter underflow");--count_[p];}
            else alive_[p]=0;
        }
    }
    void count_vertex(Vertex v,int s,bool forward) {
        for(Vertex occurrence:index_.touching(v)) {
            if(forward)++work_.forward_reads;else ++work_.certificate_reads;
            const size_t p=index_.owner[occurrence];
            if(!alive_[p])continue;
            const int role=index_.pivot(occurrence);
            add(v,coefficient(static_cast<int>(count_[p])-role,s-static_cast<int>(index_.paths.holds[p])-role));
        }
    }
    void counter_forward(int s,const Vertices& order) {
        reset_counters(s);
        for(Vertex v:order){count_vertex(v,s,true);erase(v);}
    }
    void counter_certificate(int s,const Vertices& order,std::span<const Count> upper) {
        reset_counters(s);size_t at=0;
        while(at<order.size()) {
            size_t end=at+1;
            while(end<order.size() && upper[order[end]]==upper[order[at]])++end;
            for(size_t i=at;i<end;++i)count_vertex(order[i],s,false);
            for(size_t i=at;i<end;++i)erase(order[i]);
            at=end;
        }
    }
public:
    Plan(const Graph& graph,const Layout& index,const Choose& choose,int mode,Metrics& work)
        :graph_(graph),index_(index),choose_(choose),mode_(mode),work_(work) {
        if(mode_)value_.resize(graph.n);
        memory();
    }
    void invalidate(){valid_=false;}
    size_t bytes()const{return work_.scratch_bytes;}
    const std::vector<Count>& certificate_counts()const{return value_;}
    void forward(int s,const Vertices& order,std::span<Count> upper) {
        if(mode_==2 && !valid_)compile(order);
        std::fill(value_.begin(),value_.end(),0);
        if(mode_==2)compiled_forward(s);else counter_forward(s,order);
        Count level=0;
        for(Vertex v:order){level=std::max(level,value_[v]);upper[v]=level;}
    }
    bool trial(int s,const Vertices& order,std::span<Count> upper) {
        const auto start=Clock::now();const uint64_t compilations=work_.compilations;
        forward(s,order,upper);
        std::fill(value_.begin(),value_.end(),0);
        if(mode_==2)compiled_certificate(s,upper);else counter_certificate(s,order,upper);
        Layer layer;layer.s=s;layer.compiled=work_.compilations!=compilations;
        for(Vertex v=0;v<graph_.n;++v){layer.failures+=value_[v]<upper[v];layer.positive+=upper[v]!=0;}
        layer.accepted=layer.failures==0;
        if(layer.accepted)++work_.accepted;else ++work_.rejected;
        work_.layers.push_back(layer);memory();work_.replay_ms+=ms(start);
        return layer.accepted;
    }
};

template<class Count> struct Solver : Kernel<Count> {
    using Base=Kernel<Count>;
    using Combinations=typename Base::Combinations;
    using Queue=typename Base::Queue;
    using Wider=typename Base::Wider;
    using Base::infinity;
    using Base::checked_add;
    using Base::contribution;
    using Base::capped_choose;
    using Base::integer_upper;
    struct Result {typename Base::Common common;orderdp::Extra extra;Metrics replay;};
#include "order_replay_solve.inc"
};
}
