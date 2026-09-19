#pragma once

#include "split.hpp"
#include <bit>

namespace bottomup {
using namespace partial;
inline double ms(Clock::time_point start) { return 1000*elapsed(start); }

struct Work {
    uint64_t expanded=0, degree_queries=0, intersections=0, root_reads=0;
    uint64_t frontier_reads=0, resumed=0, reused=0, discarded=0, bound_tests=0;
    uint64_t count_reads=0, source_reads=0, target_reads=0, updates=0, events=0;
};

struct State {
    Vertices holds, optional, candidates;
    size_t list_bytes() const {
        return (holds.capacity()+optional.capacity()+candidates.capacity())*sizeof(Vertex);
    }
};

inline size_t frontier_bytes(const std::vector<State>& states) {
    size_t bytes=states.capacity()*sizeof(State);
    for (const auto& state : states) bytes+=state.list_bytes();
    return bytes;
}

class Search {
    const Graph& graph_;
    int size_;
    const Vertices& bound_;
    Work& work_;
    std::vector<State> next_;
    Vertices holds_,optional_;

    void save(Vertices candidates) {
        next_.push_back({holds_,optional_,std::move(candidates)});
    }
    void expand(Vertices candidates) {
        ++work_.expanded;
        if (holds_.size()+optional_.size()+candidates.size()<static_cast<size_t>(size_)) {
            ++work_.discarded;
            return;
        }
        if (candidates.size()<=1) {
            const size_t saved=optional_.size();
            optional_.insert(optional_.end(),candidates.begin(),candidates.end());
            save({});
            optional_.resize(saved);
            return;
        }
        if (holds_.size()+1==static_cast<size_t>(size_)) { save(std::move(candidates)); return; }
        require(holds_.size()+1<static_cast<size_t>(size_),"invalid pause boundary");
        Vertices universal;
        Vertex pivot=absent;
        size_t best=0;
        for (Vertex u : candidates) {
            ++work_.degree_queries;
            const size_t degree=intersection_size(candidates,graph_.row(u));
            if (degree+1==candidates.size()) universal.push_back(u);
            else if (pivot==absent || degree>best) { pivot=u; best=degree; }
        }
        const size_t saved=optional_.size();
        if (universal.size()==candidates.size()) {
            optional_.insert(optional_.end(),candidates.begin(),candidates.end());
            save({});
            optional_.resize(saved);
            return;
        }
        Vertices remaining;
        std::set_difference(candidates.begin(),candidates.end(),universal.begin(),universal.end(),
            std::back_inserter(remaining));
        candidates.swap(remaining);
        optional_.insert(optional_.end(),universal.begin(),universal.end());
        Vertices branches{pivot};
        for (Vertex u : candidates) if (u!=pivot && !graph_.adjacent(pivot,u)) branches.push_back(u);
        for (Vertex u : branches) {
            const auto at=std::lower_bound(candidates.begin(),candidates.end(),u);
            require(at!=candidates.end() && *at==u,"missing pivot branch");
            candidates.erase(at);
            auto& selected=u==pivot ? optional_ : holds_;
            selected.push_back(u);
            ++work_.intersections;
            expand(intersect(candidates,graph_.row(u)));
            selected.pop_back();
        }
        optional_.resize(saved);
    }
public:
    Search(const Graph& graph,int s,const Vertices& bound,Work& work)
        : graph_(graph),size_(s),bound_(bound),work_(work) {}

    std::vector<State> fresh() {
        for (Vertex v=0; v<graph_.n; ++v) if (bound_[v]>=static_cast<Vertex>(size_)) {
            Vertices candidates;
            const auto row=graph_.row(v);
            for (auto it=std::upper_bound(row.begin(),row.end(),v); it!=row.end(); ++it) {
                ++work_.root_reads;
                if (bound_[*it]>=static_cast<Vertex>(size_)) candidates.push_back(*it);
            }
            if (candidates.size()+1<static_cast<size_t>(size_)) continue;
            holds_.assign(1,v);
            optional_.clear();
            expand(std::move(candidates));
        }
        return std::move(next_);
    }

    std::vector<State> resume(std::vector<State> old) {
        next_.reserve(old.size());
        for (auto& state : old) {
            bool valid=true;
            for (Vertex v : state.holds) {
                ++work_.frontier_reads;
                if (bound_[v]<static_cast<Vertex>(size_)) { valid=false; break; }
            }
            if (!valid) { ++work_.discarded; state={}; continue; }
            auto remove=[&](Vertex v) {
                ++work_.frontier_reads;
                return bound_[v]<static_cast<Vertex>(size_);
            };
            std::erase_if(state.optional,remove);
            std::erase_if(state.candidates,remove);
            if (state.holds.size()+state.optional.size()+state.candidates.size()<static_cast<size_t>(size_)) {
                ++work_.discarded;
                state={};
                continue;
            }
            if (state.candidates.empty()) {
                ++work_.reused;
                next_.push_back(std::move(state));
            } else {
                ++work_.resumed;
                holds_=std::move(state.holds);
                optional_=std::move(state.optional);
                expand(std::move(state.candidates));
            }
        }
        return std::move(next_);
    }
};

inline Index materialize(const std::vector<State>& states,Vertex n,int s) {
    Index index;
    index.holds.reserve(states.size());
    index.off.reserve(states.size()+1);
    size_t members=0;
    for (const auto& state : states) members+=state.holds.size()+state.optional.size()+state.candidates.size();
    index.vertices.reserve(members);
    for (const auto& state : states) {
        require(state.candidates.empty() || state.holds.size()+1==static_cast<size_t>(s),"unresolved state at wrong size");
        index.append(state.holds,state.optional);
        index.vertices.insert(index.vertices.end(),state.candidates.begin(),state.candidates.end());
        index.off.back()=index.vertices.size();
    }
    index.transpose(n);
    return index;
}

inline void update_bounds(Vertices& bound,std::span<const Count> labels,int s,Work& work) {
    for (Vertex v=0; v<bound.size(); ++v) if (bound[v]>=static_cast<Vertex>(s)) {
        Vertex low=s-1,high=bound[v]+1;
        while (high-low>1) {
            const Vertex middle=low+(high-low)/2;
            ++work.bound_tests;
            if (allsize::choose_at_most(middle-1,s-1,labels[v])) low=middle;
            else high=middle;
        }
        bound[v]=low;
    }
}

inline std::vector<uint64_t> clique_masks(const Graph& graph,int s) {
    require(graph.n<20,"audit graph too large");
    std::vector<uint64_t> masks;
    for (uint64_t mask=1; mask<(uint64_t{1}<<graph.n); ++mask) {
        if (std::popcount(mask)!=s) continue;
        bool clique=true;
        for (Vertex u=0; u<graph.n && clique; ++u) if ((mask>>u)&1)
            for (Vertex v=u+1; v<graph.n; ++v)
                if (((mask>>v)&1) && !graph.adjacent(u,v)) { clique=false; break; }
        if (clique) masks.push_back(mask);
    }
    return masks;
}

inline void audit_index(const Graph& graph,const Index& index,int s) {
    std::vector<uint64_t> actual;
    for (size_t p=0; p<index.size(); ++p) {
        const auto row=index.row(p);
        if (index.holds[p]>static_cast<Vertex>(s) || row.size()<static_cast<size_t>(s)) continue;
        const size_t h=index.holds[p],q=row.size()-h;
        uint64_t required=0,all=0;
        for (size_t i=0; i<row.size(); ++i) {
            require(row[i]<graph.n && !(all&(uint64_t{1}<<row[i])),"invalid duplicate path member");
            all|=uint64_t{1}<<row[i];
            if (i<h) required|=uint64_t{1}<<row[i];
        }
        for (uint64_t mask=0; mask<(uint64_t{1}<<q); ++mask) if (std::popcount(mask)==s-static_cast<int>(h)) {
            uint64_t clique=required;
            for (size_t i=0; i<q; ++i) if ((mask>>i)&1) clique|=uint64_t{1}<<row[h+i];
            actual.push_back(clique);
        }
    }
    std::sort(actual.begin(),actual.end());
    require(std::adjacent_find(actual.begin(),actual.end())==actual.end(),"duplicate represented clique");
    require(actual==clique_masks(graph,s),"current frontier clique coverage differs");
}

struct PeelResult {
    std::vector<Count> core;
    size_t state_bytes=0;
    uint64_t audits=0;
};

template<bool Audit=false>
PeelResult peel(const Graph& graph,const Index& index,int s,const Combinations& choose,Work& work) {
    const Vertex n=graph.n;
    std::vector<Count> support(n,0),keys(n),wh(index.size(),0),wp(index.size(),0);
    Vertices count(index.size(),0);
    std::vector<uint8_t> touched(index.size(),0),dead(index.size(),1),dirty(n,0);
    for (size_t p=0; p<index.size(); ++p) {
        const auto row=index.row(p);
        if (index.holds[p]>static_cast<Vertex>(s) || row.size()<static_cast<size_t>(s)) continue;
        dead[p]=0;
        count[p]=row.size()-index.holds[p];
        const int needed=s-index.holds[p];
        wh[p]=choose(count[p],needed);
        wp[p]=choose(static_cast<int>(count[p])-1,needed-1);
        for (size_t i=0; i<row.size(); ++i) {
            ++work.count_reads;
            checked_add(support[row[i]],i<index.holds[p] ? wh[p] : wp[p]);
        }
    }
    keys=support;
    sizesplit::LayerQueue heap(keys);
    PeelResult result;
    result.core.assign(n,0);
    Vertices batch,affected,changed;
    std::vector<uint64_t> masks;
    if constexpr (Audit) masks=clique_masks(graph,s);
    auto audit=[&] {
        if constexpr (Audit) {
            ++result.audits;
            uint64_t live=0;
            for (Vertex v=0; v<n; ++v) if (heap.contains(v)) live|=uint64_t{1}<<v;
            std::vector<Count> actual(n,0);
            for (uint64_t mask : masks) if ((mask&live)==mask)
                for (Vertex v=0; v<n; ++v) if ((mask>>v)&1) ++actual[v];
            for (Vertex v=0; v<n; ++v) if (heap.contains(v))
                require(actual[v]==support[v] && keys[v]==support[v],"residual support differs");
            heap.audit();
        }
    };
    auto memory=[&] {
        const size_t bytes=(support.capacity()+keys.capacity()+wh.capacity()+wp.capacity())*sizeof(Count)
            +(count.capacity()+batch.capacity()+affected.capacity()+changed.capacity())*sizeof(Vertex)
            +touched.capacity()+dead.capacity()+dirty.capacity()+heap.bytes()+result.core.capacity()*sizeof(Count);
        result.state_bytes=std::max(result.state_bytes,bytes);
    };
    audit();
    Count level=0;
    while (!heap.empty()) {
        level=std::max(level,heap.key(heap.first()));
        batch.clear();
        while (!heap.empty() && heap.key(heap.first())<=level) {
            const Vertex v=heap.first();
            heap.erase(v);
            result.core[v]=level;
            batch.push_back(v);
            ++work.events;
        }
        if (heap.empty()) { memory(); break; }
        affected.clear(); changed.clear();
        for (Vertex v : batch) for (Vertex code : index.touching(v)) {
            ++work.source_reads;
            const Vertex p=code/2;
            if (dead[p]) continue;
            if (!touched[p]) { touched[p]=1; affected.push_back(p); }
            if (code&1) { require(count[p]>0,"pivot count underflow"); --count[p]; }
            else dead[p]=1;
        }
        for (Vertex p : affected) {
            const int needed=s-index.holds[p];
            if (count[p]<static_cast<Vertex>(needed)) dead[p]=1;
            const Count next_h=dead[p] ? 0 : choose(count[p],needed);
            const Count next_p=dead[p] ? 0 : choose(static_cast<int>(count[p])-1,needed-1);
            require(next_h<=wh[p] && next_p<=wp[p],"negative path loss");
            const Count loss_h=wh[p]-next_h,loss_p=wp[p]-next_p;
            wh[p]=next_h; wp[p]=next_p; touched[p]=0;
            const auto row=index.row(p);
            for (size_t i=0; i<row.size(); ++i) {
                ++work.target_reads;
                const Count delta=i<index.holds[p] ? loss_h : loss_p;
                const Vertex v=row[i];
                if (!delta || !heap.contains(v)) continue;
                require(support[v]>=delta,"support underflow");
                support[v]-=delta;
                if (!dirty[v]) { dirty[v]=1; changed.push_back(v); }
            }
        }
        for (Vertex v : changed) {
            dirty[v]=0;
            heap.decrease(v,support[v]);
            ++work.updates;
        }
        memory();
        audit();
    }
    memory();
    return result;
}

struct Layer {
    int s=0;
    size_t vertices=0,ordinary_vertices=0,paths=0,members=0,unresolved=0,frontier_bytes=0,index_bytes=0;
    double search_ms=0,index_ms=0,peel_ms=0,bounds_ms=0;
};

struct Result {
    std::vector<Count> core;
    std::vector<Layer> layers;
    Work work;
    double search_ms=0,index_ms=0,peel_ms=0,bounds_ms=0,total_ms=0;
    size_t frontier_bytes=0,index_bytes=0,state_bytes=0,simultaneous_bytes=0;
    uint64_t audits=0;
};

template<bool Audit=false>
Result run(const Graph& graph,const Vertices& ordinary,int maximum,const Combinations& choose,const std::string& mode) {
    const auto start=Clock::now();
    require(maximum>=2 && maximum<=32,"unsupported size range");
    require(mode=="eager" || mode=="fixed" || mode=="rebuild" || mode=="resume" || mode=="resume-ordinary","unknown bottom-up mode");
    Result result;
    result.core.assign(static_cast<size_t>(maximum+1)*graph.n,0);
    std::copy(ordinary.begin(),ordinary.end(),result.core.begin()+2*static_cast<size_t>(graph.n));
    Vertices bound(graph.n);
    for (Vertex v=0; v<graph.n; ++v) bound[v]=std::min<Vertex>(maximum,ordinary[v]+1);
    std::vector<State> frontier;
    Index full;
    if (mode=="eager" && maximum>=3) {
        const auto building=Clock::now();
        allsize::Search(graph,maximum,full).run();
        result.search_ms=ms(building);
        const auto indexing=Clock::now();
        full.transpose(graph.n);
        result.index_ms=ms(indexing);
    }
    for (int s=3; s<=maximum; ++s) {
        Layer layer;
        layer.s=s;
        for (Vertex v=0; v<graph.n; ++v) {
            layer.vertices+=bound[v]>=static_cast<Vertex>(s);
            layer.ordinary_vertices+=ordinary[v]+1>=static_cast<Vertex>(s);
        }
        if (!layer.vertices) {
            result.layers.push_back(layer);
            continue;
        }
        Index current;
        if (mode!="eager") {
            const auto building=Clock::now();
            Search search(graph,s,bound,result.work);
            if ((mode=="resume" || mode=="resume-ordinary") && s>3) frontier=search.resume(std::move(frontier));
            else frontier=search.fresh();
            layer.search_ms=ms(building);
            layer.frontier_bytes=bottomup::frontier_bytes(frontier);
            for (const auto& state : frontier) layer.unresolved+=!state.candidates.empty();
            const auto indexing=Clock::now();
            current=materialize(frontier,graph.n,s);
            layer.index_ms=ms(indexing);
        }
        const Index& index=mode=="eager" ? full : current;
        layer.index_bytes=index.capacity_bytes();
        layer.paths=index.size();
        layer.members=index.vertices.size();
        if constexpr (Audit) audit_index(graph,index,s);
        const auto peeling=Clock::now();
        const auto labels=peel<Audit>(graph,index,s,choose,result.work);
        layer.peel_ms=ms(peeling);
        result.audits+=labels.audits;
        std::copy(labels.core.begin(),labels.core.end(),result.core.begin()+static_cast<size_t>(s)*graph.n);
        if (mode=="resume" || mode=="rebuild") {
            const auto updating=Clock::now();
            update_bounds(bound,labels.core,s,result.work);
            layer.bounds_ms=ms(updating);
        }
        if constexpr (Audit) {
            for (int t=s+1; t<=maximum; ++t) for (uint64_t mask : clique_masks(graph,t))
                for (Vertex v=0; v<graph.n; ++v) if ((mask>>v)&1)
                    require(bound[v]>=static_cast<Vertex>(t),"future clique vertex incorrectly excluded");
        }
        result.search_ms+=layer.search_ms;
        result.index_ms+=layer.index_ms;
        result.peel_ms+=layer.peel_ms;
        result.bounds_ms+=layer.bounds_ms;
        result.frontier_bytes=std::max(result.frontier_bytes,layer.frontier_bytes);
        result.index_bytes=std::max(result.index_bytes,layer.index_bytes);
        result.state_bytes=std::max(result.state_bytes,labels.state_bytes);
        result.simultaneous_bytes=std::max(result.simultaneous_bytes,
            layer.frontier_bytes+layer.index_bytes+labels.state_bytes+bound.capacity()*sizeof(Vertex));
        result.layers.push_back(layer);
    }
    result.total_ms=ms(start);
    return result;
}
}
