#pragma once

#include "floor.hpp"
#include <array>

namespace sizesplit {
using namespace envelope;

class LayerQueue {
    std::span<Count> keys_;
    Vertices heap_, position_;
    bool less(Vertex a, Vertex b) const {
        return keys_[a] < keys_[b] || (keys_[a] == keys_[b] && a < b);
    }
    void swap_at(size_t a, size_t b) {
        std::swap(heap_[a], heap_[b]);
        position_[heap_[a]] = a;
        position_[heap_[b]] = b;
    }
    void up(size_t at) {
        while (at && less(heap_[at], heap_[(at-1)/2])) {
            const size_t parent = (at-1)/2;
            swap_at(at, parent);
            at = parent;
        }
    }
    void down(size_t at) {
        while (2*at+1 < heap_.size()) {
            size_t child = 2*at+1;
            if (child+1 < heap_.size() && less(heap_[child+1], heap_[child])) ++child;
            if (!less(heap_[child], heap_[at])) break;
            swap_at(at, child);
            at = child;
        }
    }
public:
    Count level = 0;
    explicit LayerQueue(std::span<Count> keys): keys_(keys), position_(keys.size(), absent) {
        heap_.reserve(std::count_if(keys.begin(), keys.end(), [](Count x){ return x != 0; }));
        for (Vertex v=0; v<keys.size(); ++v) if (keys[v]) {
            position_[v] = heap_.size();
            heap_.push_back(v);
        }
        for (size_t i=heap_.size()/2; i; --i) down(i-1);
    }
    bool empty() const { return heap_.empty(); }
    bool contains(Vertex v) const { return position_[v] != absent; }
    Vertex first() const { return heap_.front(); }
    Count key(Vertex v) const { return keys_[v]; }
    size_t bytes() const { return (heap_.capacity()+position_.capacity())*sizeof(Vertex); }
    void erase(Vertex v) {
        require(contains(v), "erase missing layer handle");
        const size_t at = position_[v];
        swap_at(at, heap_.size()-1);
        heap_.pop_back();
        position_[v] = absent;
        if (at < heap_.size()) {
            if (at && less(heap_[at], heap_[(at-1)/2])) up(at);
            else down(at);
        }
    }
    void decrease(Vertex v, Count value) {
        require(contains(v) && value<=keys_[v], "invalid layer decrease");
        keys_[v] = value;
        up(position_[v]);
    }
    void audit() const {
        for (size_t i=0; i<heap_.size(); ++i) {
            require(position_[heap_[i]]==i, "layer handle corrupt");
            require(!i || !less(heap_[i], heap_[(i-1)/2]), "layer heap corrupt");
            require(keys_[heap_[i]]>=level, "layer key below floor");
        }
    }
};

struct Report {
    allsize::Result data;
    uint64_t assignments=0, shared_events=0, shared_assignments=0, splits=0;
    uint64_t eligibility_checks=0, scalar_target_checks=0, copied_bytes=0;
    uint64_t ghost_steps=0, state_updates=0, audited_states=0;
    uint64_t trace_hash=1469598103934665603ULL;
    size_t peak_group_states=0, output_bytes=0;
};

struct Config {
    bool share_rows=true;
    bool force_singletons=false;
};

struct State {
    Vertices pivots;
    std::vector<uint8_t> holds_live, removed;
    size_t bytes() const {
        return pivots.capacity()*sizeof(Vertex)+holds_live.capacity()+removed.capacity();
    }
};

template<bool Audit=false>
Report peel(const Layout& index, Vertex n, const Combinations& choose,
            const Vertices& ordinary, Config config={}, const std::vector<Count>* expected=nullptr) {
    const auto start = Clock::now();
    require(index.maximum>=2 && index.maximum<=32, "unsupported maximum size");
    Report report;
    auto& r = report.data;
    r.core.assign(static_cast<size_t>(index.maximum+1)*n, 0);
    std::copy(ordinary.begin(), ordinary.end(), r.core.begin()+2*static_cast<size_t>(n));
    report.output_bytes = r.core.capacity()*sizeof(Count);
    if (index.maximum==2) { r.peel_ms=ms(start); return report; }
    auto cell = [&](int s, Vertex v) -> Count& { return r.core[static_cast<size_t>(s)*n+v]; };
    for (size_t p=0; p<index.paths.size(); ++p) {
        const int q = index.paths.row(p).size()-index.paths.holds[p];
        for (int s=std::max<int>(3,index.lo[p]); s<=static_cast<int>(index.hi[p]); ++s) {
            const Count h=contribution(index,p,false,s,q,choose), z=contribution(index,p,true,s,q,choose);
            for (size_t i=index.paths.off[p]; i<index.paths.off[p+1]; ++i) {
                ++r.work.count_reads;
                checked_add(cell(s,index.paths.vertices[i]), index.pivot(i) ? z : h);
            }
        }
    }
    std::vector<LayerQueue> layers;
    layers.reserve(index.maximum-2);
    Vertices sizes;
    for (int s=3; s<=index.maximum; ++s) {
        layers.emplace_back(std::span<Count>(r.core).subspan(static_cast<size_t>(s)*n,n));
        if (!layers.back().empty()) sizes.push_back(s);
    }
    size_t heap_bytes=layers.capacity()*sizeof(LayerQueue);
    for (const auto& layer : layers) heap_bytes+=layer.bytes();
    State initial;
    initial.pivots.resize(index.paths.size());
    initial.holds_live.assign(index.paths.size(),1);
    for (size_t p=0; p<index.paths.size(); ++p)
        initial.pivots[p]=index.paths.row(p).size()-index.paths.holds[p];
    if constexpr (Audit) initial.removed.assign(n,0);
    size_t live_bytes=initial.bytes(), live_states=1;
    r.state_bytes=heap_bytes+live_bytes;
    report.peak_group_states=1;
    auto mix = [&](uint64_t value) { report.trace_hash^=value; report.trace_hash*=1099511628211ULL; };
    auto audit = [&](const Vertices& active, const State& state) {
        if constexpr (Audit) {
            require(expected!=nullptr, "audit needs oracle");
            ++report.audited_states;
            for (size_t p=0; p<index.paths.size(); ++p) {
                bool alive=true;
                Vertex q=0;
                for (size_t i=index.paths.off[p]; i<index.paths.off[p+1]; ++i) {
                    const bool live=!state.removed[index.paths.vertices[i]];
                    if (index.pivot(i)) q+=live;
                    else alive &= live;
                }
                require(alive==bool(state.holds_live[p]), "hold flag differs from residual set");
                if (alive) require(q==state.pivots[p], "pivot count differs from residual set");
            }
            for (Vertex s : active) {
                const auto& layer=layers[s-3];
                layer.audit();
                for (Vertex v=0; v<n; ++v) {
                    if (layer.contains(v)) {
                        require(!state.removed[v], "queued vertex physically absent");
                        Count raw=0;
                        for (Vertex i : index.touching(v)) {
                            const size_t p=index.owner[i];
                            if (index.valid(p,s) && state.holds_live[p])
                                checked_add(raw,contribution(index,p,index.pivot(i),s,state.pivots[p],choose));
                        }
                        require(cell(s,v)==std::max(layer.level,raw), "clipped count differs");
                        require((*expected)[static_cast<size_t>(s)*n+v]>=layer.level, "floor exceeds pending answer");
                    } else {
                        require(cell(s,v)==(*expected)[static_cast<size_t>(s)*n+v], "completed answer differs");
                        if (!state.removed[v]) require(cell(s,v)==0, "nonzero ghost vertex");
                    }
                }
            }
        }
    };
    struct Loss { Vertex s; Count h,z,k; };
    std::function<void(Vertices,State&)> run;
    run = [&](Vertices active, State& state) {
        audit(active,state);
        while (!active.empty()) {
            std::erase_if(active,[&](Vertex s){ return layers[s-3].empty(); });
            if (active.empty()) break;
            const Vertex v=layers[active.front()-3].first();
            std::array<Vertex,33> yes,no;
            size_t yes_size=0,no_size=0;
            for (Vertex s : active) {
                ++report.eligibility_checks;
                const auto& layer=layers[s-3];
                if (!layer.contains(v) || layer.key(v)==layer.key(layer.first())) yes[yes_size++]=s;
                else no[no_size++]=s;
            }
            if (no_size) {
                ++report.splits;
                State other=state;
                report.copied_bytes+=other.bytes();
                live_bytes+=other.bytes();
                ++live_states;
                r.state_bytes=std::max(r.state_bytes,heap_bytes+live_bytes);
                report.peak_group_states=std::max(report.peak_group_states,live_states);
                run(Vertices(yes.begin(),yes.begin()+yes_size),state);
                run(Vertices(no.begin(),no.begin()+no_size),other);
                live_bytes-=other.bytes();
                --live_states;
                return;
            }
            ++r.work.events;
            uint64_t count=0, mask=0;
            for (Vertex s : active) {
                auto& layer=layers[s-3];
                if (!layer.contains(v)) { ++report.ghost_steps; continue; }
                require(layer.key(v)>=layer.level, "decreasing assigned layer level");
                layer.level=layer.key(v);
                if constexpr (Audit)
                    require(cell(s,v)==(*expected)[static_cast<size_t>(s)*n+v], "shared assignment differs");
                layer.erase(v);
                ++count;
                mask|=uint64_t{1}<<s;
            }
            report.assignments+=count;
            if (count>1) { ++report.shared_events; report.shared_assignments+=count; }
            mix(v); mix(mask);
            for (Vertex occurrence : index.touching(v)) {
                ++r.work.source_reads;
                const size_t p=index.owner[occurrence];
                if (!state.holds_live[p]) continue;
                const bool pivot=index.pivot(occurrence);
                const int q=state.pivots[p];
                std::array<Loss,33> losses;
                size_t used=0;
                for (Vertex s : active) if (index.valid(p,s)) {
                    Count h=contribution(index,p,false,s,q,choose), z=contribution(index,p,true,s,q,choose);
                    if (pivot) {
                        h-=contribution(index,p,false,s,q-1,choose);
                        z-=contribution(index,p,true,s,q-1,choose);
                    }
                    if (h || z) losses[used++]={s,h,z,layers[s-3].level};
                    if constexpr (Audit) if (!(mask & (uint64_t{1}<<s)))
                        require(!h && !z, "ghost has nonzero target loss");
                }
                if (pivot) {
                    require(state.pivots[p]>0, "negative pivot counter");
                    --state.pivots[p];
                } else state.holds_live[p]=0;
                ++report.state_updates;
                if (!used) continue;
                auto update = [&](size_t i,const Loss& loss) {
                    ++report.scalar_target_checks;
                    const Vertex u=index.paths.vertices[i];
                    auto& layer=layers[loss.s-3];
                    const Count amount=index.pivot(i) ? loss.z : loss.h;
                    if (amount && layer.contains(u) && layer.key(u)>loss.k) {
                        const Count next=amount>=layer.key(u)-loss.k ? loss.k : layer.key(u)-amount;
                        layer.decrease(u,next);
                        ++r.work.updates;
                    }
                };
                if (config.share_rows) {
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
            audit(active,state);
        }
    };
    if (config.force_singletons) {
        for (Vertex s : sizes) {
            State state=initial;
            report.copied_bytes+=state.bytes();
            live_bytes+=state.bytes();
            r.state_bytes=std::max(r.state_bytes,heap_bytes+live_bytes);
            report.peak_group_states=std::max<size_t>(report.peak_group_states,2);
            run(Vertices{s},state);
            live_bytes-=state.bytes();
        }
    } else run(std::move(sizes),initial);
    r.peel_ms=ms(start);
    return report;
}
}
