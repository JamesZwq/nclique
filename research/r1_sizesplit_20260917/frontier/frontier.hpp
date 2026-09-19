#pragma once

#include "../split.hpp"
#include <bit>

namespace commonfront {
using namespace sizesplit;

struct Report : sizesplit::Report {
    uint64_t exposures=0, ready_seed_reads=0, split_vertex_reads=0, ready_notifications=0;
};

template<bool Audit=false>
Report peel(const Layout& index, Vertex n, const Combinations& choose,
            const Vertices& ordinary, bool share_rows=true, const std::vector<Count>* expected=nullptr) {
    const auto start=Clock::now();
    require(index.maximum>=2 && index.maximum<=32,"unsupported maximum size");
    Report report;
    auto& r=report.data;
    r.core.assign(static_cast<size_t>(index.maximum+1)*n,0);
    std::copy(ordinary.begin(),ordinary.end(),r.core.begin()+2*static_cast<size_t>(n));
    report.output_bytes=r.core.capacity()*sizeof(Count);
    if (index.maximum==2) { r.peel_ms=ms(start); return report; }
    auto cell=[&](Vertex s,Vertex v)->Count& { return r.core[static_cast<size_t>(s)*n+v]; };
    for (size_t p=0; p<index.paths.size(); ++p) {
        const int q=index.paths.row(p).size()-index.paths.holds[p];
        for (int s=std::max<int>(3,index.lo[p]); s<=static_cast<int>(index.hi[p]); ++s) {
            const Count h=contribution(index,p,false,s,q,choose), z=contribution(index,p,true,s,q,choose);
            for (size_t i=index.paths.off[p]; i<index.paths.off[p+1]; ++i) {
                ++r.work.count_reads;
                checked_add(cell(s,index.paths.vertices[i]),index.pivot(i) ? z : h);
            }
        }
    }
    std::vector<LayerQueue> layers;
    layers.reserve(index.maximum-2);
    std::vector<uint64_t> eligible(n,0),pending(n,0);
    std::array<Vertex,33> waiting{};
    Vertices sizes;
    for (int s=3; s<=index.maximum; ++s) {
        const uint64_t bit=uint64_t{1}<<s;
        for (Vertex v=0; v<n; ++v) {
            if (cell(s,v)) pending[v]|=bit;
            else eligible[v]|=bit;
        }
        layers.emplace_back(std::span<Count>(r.core).subspan(static_cast<size_t>(s)*n,n));
        if (!layers.back().empty()) sizes.push_back(s);
    }
    size_t base_bytes=layers.capacity()*sizeof(LayerQueue)+(eligible.capacity()+pending.capacity())*sizeof(uint64_t);
    for (const auto& layer : layers) base_bytes+=layer.bytes();
    State initial;
    initial.pivots.resize(index.paths.size());
    initial.holds_live.assign(index.paths.size(),1);
    for (size_t p=0; p<index.paths.size(); ++p)
        initial.pivots[p]=index.paths.row(p).size()-index.paths.holds[p];
    if constexpr (Audit) initial.removed.assign(n,0);
    size_t live_bytes=initial.bytes(),live_states=1;
    r.state_bytes=base_bytes+live_bytes;
    report.peak_group_states=1;
    auto mask_for=[](const Vertices& active) {
        uint64_t mask=0;
        for (Vertex s : active) mask|=uint64_t{1}<<s;
        return mask;
    };
    auto is_ready=[&](Vertex v,uint64_t mask) {
        return (pending[v]&mask) && (eligible[v]&mask)==mask;
    };
    auto expose=[&](Vertex s,uint64_t mask,Vertices* ready) {
        auto& layer=layers[s-3];
        if (!waiting[s] && !layer.empty()) layer.level=layer.key(layer.first());
        while (!layer.empty() && layer.key(layer.first())==layer.level) {
            const Vertex u=layer.first();
            const uint64_t bit=uint64_t{1}<<s;
            require((pending[u]&bit) && !(eligible[u]&bit),"duplicate minimum exposure");
            layer.erase(u);
            eligible[u]|=bit;
            ++waiting[s];
            ++report.exposures;
            if (ready) {
                ++report.ready_notifications;
                if (is_ready(u,mask)) ready->push_back(u);
            }
        }
    };
    auto audit=[&](const Vertices& active,const State& state) {
        if constexpr (Audit) {
            require(expected!=nullptr,"audit requires oracle");
            ++report.audited_states;
            for (size_t p=0; p<index.paths.size(); ++p) {
                bool alive=true;
                Vertex q=0;
                for (size_t i=index.paths.off[p]; i<index.paths.off[p+1]; ++i) {
                    const bool live=!state.removed[index.paths.vertices[i]];
                    if (index.pivot(i)) q+=live;
                    else alive &= live;
                }
                require(alive==bool(state.holds_live[p]),"hold state differs");
                if (alive) require(q==state.pivots[p],"pivot state differs");
            }
            for (Vertex s : active) {
                const auto& layer=layers[s-3];
                const uint64_t bit=uint64_t{1}<<s;
                layer.audit();
                Vertex exposed=0;
                for (Vertex v=0; v<n; ++v) {
                    if (pending[v]&bit) {
                        require(!state.removed[v],"pending vertex physically absent");
                        Count raw=0;
                        for (Vertex i : index.touching(v)) {
                            const size_t p=index.owner[i];
                            if (index.valid(p,s) && state.holds_live[p])
                                checked_add(raw,contribution(index,p,index.pivot(i),s,state.pivots[p],choose));
                        }
                        require(cell(s,v)==std::max(layer.level,raw),"clipped count differs");
                        require((*expected)[static_cast<size_t>(s)*n+v]>=layer.level,"floor exceeds answer");
                        if (eligible[v]&bit) {
                            ++exposed;
                            require(!layer.contains(v) && cell(s,v)==layer.level,"invalid exposed minimum");
                        } else require(layer.contains(v) && cell(s,v)>layer.level,"incomplete exposure");
                    } else {
                        require(!layer.contains(v),"completed vertex still queued");
                        require(cell(s,v)==(*expected)[static_cast<size_t>(s)*n+v],"completed answer differs");
                        if (!state.removed[v]) require(cell(s,v)==0,"nonzero ghost");
                    }
                }
                require(exposed==waiting[s],"minimum pending count differs");
                require(waiting[s] || layer.empty(),"unfinished layer has no exposed minimum");
            }
        }
    };
    for (Vertex s : sizes) expose(s,0,nullptr);
    struct Loss { Vertex s; Count h,z,k; };
    auto mix=[&](uint64_t x) { report.trace_hash^=x; report.trace_hash*=1099511628211ULL; };
    std::function<void(Vertices,State&)> run;
    run=[&](Vertices active,State& state) {
        std::erase_if(active,[&](Vertex s){ return !waiting[s] && layers[s-3].empty(); });
        if (active.empty()) return;
        uint64_t mask=mask_for(active);
        Vertices ready;
        ready.reserve(n);
        live_bytes+=ready.capacity()*sizeof(Vertex);
        r.state_bytes=std::max(r.state_bytes,base_bytes+live_bytes);
        for (Vertex v=0; v<n; ++v) {
            ++report.ready_seed_reads;
            if (is_ready(v,mask)) ready.push_back(v);
        }
        size_t cursor=0;
        audit(active,state);
        while (true) {
            std::erase_if(active,[&](Vertex s){ return !waiting[s] && layers[s-3].empty(); });
            mask=mask_for(active);
            if (active.empty()) break;
            if (cursor==ready.size()) {
                if constexpr (Audit) for (Vertex u=0; u<n; ++u)
                    require(!is_ready(u,mask),"missed a common minimum");
                Vertex selected=absent;
                const uint64_t bit=uint64_t{1}<<active.front();
                for (Vertex v=0; v<n; ++v) {
                    ++report.split_vertex_reads;
                    if ((pending[v]&bit) && (eligible[v]&bit)) { selected=v; break; }
                }
                require(selected!=absent,"no split vertex in unfinished layer");
                Vertices yes,no;
                for (Vertex s : active)
                    ((eligible[selected]&(uint64_t{1}<<s)) ? yes : no).push_back(s);
                require(!yes.empty() && !no.empty(),"split without a conflict");
                ++report.splits;
                live_bytes-=ready.capacity()*sizeof(Vertex);
                Vertices().swap(ready);
                State other=state;
                report.copied_bytes+=other.bytes();
                live_bytes+=other.bytes();
                ++live_states;
                report.peak_group_states=std::max(report.peak_group_states,live_states);
                r.state_bytes=std::max(r.state_bytes,base_bytes+live_bytes);
                run(std::move(yes),state);
                run(std::move(no),other);
                live_bytes-=other.bytes();
                --live_states;
                return;
            }
            const Vertex v=ready[cursor++];
            require(is_ready(v,mask),"stale ready candidate");
            const uint64_t assigned=pending[v]&mask;
            const uint64_t count=std::popcount(assigned);
            pending[v]&=~mask;
            ++r.work.events;
            report.assignments+=count;
            if (count>1) { ++report.shared_events; report.shared_assignments+=count; }
            mix(v); mix(assigned);
            for (Vertex s : active) {
                if (!(assigned&(uint64_t{1}<<s))) { ++report.ghost_steps; continue; }
                require(waiting[s]>0,"negative minimum bucket");
                --waiting[s];
                if constexpr (Audit)
                    require(cell(s,v)==(*expected)[static_cast<size_t>(s)*n+v],"frontier assignment differs");
            }
            for (Vertex occurrence : index.touching(v)) {
                ++r.work.source_reads;
                const size_t p=index.owner[occurrence];
                if (!state.holds_live[p]) continue;
                const bool pivot=index.pivot(occurrence);
                const int q=state.pivots[p];
                std::array<Loss,33> losses;
                size_t used=0;
                for (Vertex s : active) if (index.valid(p,s)) {
                    Count h=contribution(index,p,false,s,q,choose),z=contribution(index,p,true,s,q,choose);
                    if (pivot) {
                        h-=contribution(index,p,false,s,q-1,choose);
                        z-=contribution(index,p,true,s,q-1,choose);
                    }
                    if (h || z) losses[used++]={s,h,z,layers[s-3].level};
                    if constexpr (Audit) if (!(assigned&(uint64_t{1}<<s)))
                        require(!h && !z,"ghost has positive loss");
                }
                if (pivot) {
                    require(state.pivots[p]>0,"negative pivot count");
                    --state.pivots[p];
                } else state.holds_live[p]=0;
                ++report.state_updates;
                if (!used) continue;
                auto update=[&](size_t i,const Loss& loss) {
                    ++report.scalar_target_checks;
                    const Vertex u=index.paths.vertices[i];
                    auto& layer=layers[loss.s-3];
                    const Count amount=index.pivot(i) ? loss.z : loss.h;
                    if (amount && layer.contains(u) && layer.key(u)>loss.k) {
                        layer.decrease(u,amount>=layer.key(u)-loss.k ? loss.k : layer.key(u)-amount);
                        ++r.work.updates;
                    }
                };
                if (share_rows) {
                    for (size_t i=index.paths.off[p]; i<index.paths.off[p+1]; ++i) {
                        ++r.work.target_reads;
                        for (size_t j=0; j<used; ++j) update(i,losses[j]);
                    }
                } else {
                    for (size_t j=0; j<used; ++j)
                        for (size_t i=index.paths.off[p]; i<index.paths.off[p+1]; ++i) {
                            ++r.work.target_reads;
                            update(i,losses[j]);
                        }
                }
            }
            if constexpr (Audit) state.removed[v]=1;
            for (Vertex s : active) expose(s,mask,&ready);
            require(ready.size()<=n,"ready queue contains duplicates");
            audit(active,state);
        }
        live_bytes-=ready.capacity()*sizeof(Vertex);
    };
    run(std::move(sizes),initial);
    r.peel_ms=ms(start);
    return report;
}
}
