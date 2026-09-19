#pragma once

#include "bottomup.hpp"
#include <array>

namespace tworoads {
using namespace allsize;

enum class Bound { None, Ordinary, Ratio, Integer };

struct Statistics {
    uint64_t batches=0,initial_capped=0,tighter_bounds=0,early_removals=0;
    uint64_t unchanged_keys=0,cache_hits=0,cache_misses=0,binomial_steps=0,audits=0;
    uint64_t components=0,largest_component=0,component_vertices=0,component_scans=0;
    uint64_t splits=0,shared_events=0,shared_assignments=0,scalar_target_checks=0;
    uint64_t copied_index_bytes=0;
    double bounds_ms=0,components_ms=0;
};

struct Result {
    allsize::Result data;
    Statistics stats;
};

inline Count capped_choose(uint64_t n,int r,Count cap,Statistics& stats) {
    if (!cap || r<0 || static_cast<uint64_t>(r)>n) return 0;
    r=static_cast<int>(std::min<uint64_t>(r,n-r));
    unsigned __int128 value=1;
    for (int j=1; j<=r; ++j) {
        ++stats.binomial_steps;
        value=value*(n-j+1)/j;
        if (value>=cap) return cap;
    }
    return static_cast<Count>(std::min<unsigned __int128>(cap,value));
}

inline Count integer_upper(Count a,int r,Vertex maximum,Statistics& stats) {
    Count answer=0;
    uint64_t upper=maximum;
    for (int j=r; j>=1 && a; --j) {
        require(upper>=static_cast<uint64_t>(j),"invalid canonical upper argument");
        uint64_t low=j-1,high=upper+1;
        while (high-low>1) {
            const uint64_t middle=low+(high-low)/2;
            if (capped_choose(middle,j,a+1,stats)<=a) low=middle;
            else high=middle;
        }
        const Count term=capped_choose(low,j,a+1,stats);
        require(term<=a && term>0,"invalid canonical term");
        a-=term;
        const Count next=capped_choose(low,j+1,infinity-1-answer,stats);
        answer+=next;
        if (answer==infinity-1) return answer;
        upper=low-1;
    }
    require(a==0,"incomplete canonical expansion");
    return answer;
}

template<Bound Kind,bool Audit=false>
Result capped(const Graph& graph,const Layout& index,const Combinations& choose,
              const Vertices& ordinary,const std::vector<Count>* expected=nullptr) {
    const auto start=Clock::now();
    const Vertex n=graph.n;
    Result result;
    auto& data=result.data;
    auto& work=data.work;
    auto& stats=result.stats;
    data.core.assign(static_cast<size_t>(index.maximum+1)*n,0);
    std::copy(ordinary.begin(),ordinary.end(),data.core.begin()+2*static_cast<size_t>(n));
    const Vertex maximum=ordinary.empty() ? 0 : *std::max_element(ordinary.begin(),ordinary.end());
    struct Cache { Count key=infinity,value=0; };
    std::vector<Cache> cache;
    if constexpr (Kind==Bound::Integer) cache.resize(4096);
    for (int s=3; s<=index.maximum; ++s) {
        auto keys=std::span<Count>(data.core).subspan(static_cast<size_t>(s)*n,n);
        const auto previous=std::span<const Count>(data.core).subspan(static_cast<size_t>(s-1)*n,n);
        std::vector<Count> support(n,0),upper;
        std::vector<Count> wh(index.paths.size(),0),wp(index.paths.size(),0);
        Vertices count(index.paths.size(),0);
        std::vector<uint8_t> touched(index.paths.size(),0),dead(index.paths.size(),1),dirty(n,0);
        for (size_t p=0; p<index.paths.size(); ++p) if (index.valid(p,s)) {
            dead[p]=0;
            count[p]=index.paths.row(p).size()-index.paths.holds[p];
            wh[p]=contribution(index,p,false,s,count[p],choose);
            wp[p]=contribution(index,p,true,s,count[p],choose);
            for (size_t i=index.paths.off[p]; i<index.paths.off[p+1]; ++i) {
                ++work.count_reads;
                checked_add(support[index.paths.vertices[i]],index.pivot(i) ? wp[p] : wh[p]);
            }
        }
        if constexpr (Kind==Bound::None) std::copy(support.begin(),support.end(),keys.begin());
        else {
            const auto bounding=Clock::now();
            upper.assign(n,0);
            if constexpr (Kind==Bound::Integer) std::fill(cache.begin(),cache.end(),Cache{});
            for (Vertex v=0; v<n; ++v) {
                const Count original=capped_choose(ordinary[v],s-1,infinity-1,stats);
                Count value=original;
                if constexpr (Kind==Bound::Ratio || Kind==Bound::Integer) {
                    const unsigned __int128 ratio=ordinary[v]<static_cast<Vertex>(s-1) ? 0 :
                        static_cast<unsigned __int128>(previous[v])*(ordinary[v]-s+2)/(s-1);
                    value=static_cast<Count>(std::min<unsigned __int128>(value,ratio));
                }
                if constexpr (Kind==Bound::Integer) {
                    if (value) {
                        const Count a=previous[v];
                        Cache& entry=cache[(a^(a>>17)^(a>>37))&(cache.size()-1)];
                        if (entry.key!=a) {
                            ++stats.cache_misses;
                            entry={a,integer_upper(a,s-2,maximum,stats)};
                        } else ++stats.cache_hits;
                        value=std::min(value,entry.value);
                    }
                }
                upper[v]=value;
                keys[v]=std::min(value,support[v]);
                stats.initial_capped+=keys[v]<support[v];
                stats.tighter_bounds+=value<original;
                if constexpr (Audit) {
                    require(expected && value>=(*expected)[static_cast<size_t>(s)*n+v],"invalid transferred upper bound");
                    require(!support[v] || keys[v],"zero bound on positive initial degree");
                }
            }
            stats.bounds_ms+=ms(bounding);
        }
        sizesplit::LayerQueue heap(keys);
        Vertices batch,affected,changed;
        std::vector<uint64_t> cliques;
        if constexpr (Audit) cliques=bottomup::clique_masks(graph,s);
        auto audit=[&] {
            if constexpr (Audit) {
                ++stats.audits;
                require(expected!=nullptr,"missing oracle");
                uint64_t live=0;
                for (Vertex v=0; v<n; ++v) if (heap.contains(v)) live|=uint64_t{1}<<v;
                std::vector<Count> actual(n,0);
                for (uint64_t clique : cliques) if ((clique&live)==clique)
                    for (Vertex v=0; v<n; ++v) if ((clique>>v)&1) ++actual[v];
                for (Vertex v=0; v<n; ++v) {
                    if (heap.contains(v)) {
                        require(support[v]==actual[v],"raw residual degree differs");
                        Count key=support[v];
                        if constexpr (Kind!=Bound::None) key=std::min(key,upper[v]);
                        require(keys[v]==key,"capped priority differs");
                    } else require(keys[v]==(*expected)[static_cast<size_t>(s)*n+v],"completed capped output differs");
                }
                heap.audit();
            }
        };
        auto memory=[&] {
            const size_t bytes=(support.capacity()+upper.capacity()+wh.capacity()+wp.capacity())*sizeof(Count)
                +(count.capacity()+batch.capacity()+affected.capacity()+changed.capacity())*sizeof(Vertex)
                +touched.capacity()+dead.capacity()+dirty.capacity()+heap.bytes()+cache.capacity()*sizeof(Cache);
            data.state_bytes=std::max(data.state_bytes,bytes);
        };
        audit();
        Count level=0;
        while (!heap.empty()) {
            ++stats.batches;
            level=std::max(level,heap.key(heap.first()));
            batch.clear();
            while (!heap.empty() && heap.key(heap.first())<=level) {
                const Vertex v=heap.first();
                heap.erase(v);
                stats.early_removals+=support[v]>level;
                keys[v]=level;
                batch.push_back(v);
                ++work.events;
            }
            if (heap.empty()) { memory(); audit(); break; }
            affected.clear(); changed.clear();
            for (Vertex v : batch) for (Vertex occurrence : index.touching(v)) {
                ++work.source_reads;
                const size_t p=index.owner[occurrence];
                if (dead[p]) continue;
                if (!touched[p]) { touched[p]=1; affected.push_back(p); }
                if (index.pivot(occurrence)) {
                    require(count[p]>0,"negative optional count");
                    --count[p];
                } else dead[p]=1;
            }
            for (Vertex p : affected) {
                if (count[p]<static_cast<Vertex>(s)-index.paths.holds[p]) dead[p]=1;
                const Count next_h=dead[p] ? 0 : contribution(index,p,false,s,count[p],choose);
                const Count next_p=dead[p] ? 0 : contribution(index,p,true,s,count[p],choose);
                require(next_h<=wh[p] && next_p<=wp[p],"negative batch loss");
                const Count loss_h=wh[p]-next_h,loss_p=wp[p]-next_p;
                wh[p]=next_h; wp[p]=next_p; touched[p]=0;
                if (!loss_h && !loss_p) continue;
                for (size_t i=index.paths.off[p]; i<index.paths.off[p+1]; ++i) {
                    ++work.target_reads;
                    const Count delta=index.pivot(i) ? loss_p : loss_h;
                    const Vertex v=index.paths.vertices[i];
                    if (!delta || !heap.contains(v)) continue;
                    require(support[v]>=delta,"raw degree underflow");
                    support[v]-=delta;
                    if (!dirty[v]) { dirty[v]=1; changed.push_back(v); }
                }
            }
            for (Vertex v : changed) {
                dirty[v]=0;
                Count next=support[v];
                if constexpr (Kind!=Bound::None) next=std::min(next,upper[v]);
                if (next<keys[v]) { heap.decrease(v,next); ++work.updates; }
                else ++stats.unchanged_keys;
            }
            memory();
            audit();
        }
        memory();
    }
    data.peel_ms=ms(start);
    return result;
}
}
