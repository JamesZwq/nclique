#pragma once

#include "capped.hpp"
#include "frontier.hpp"

namespace tworoads {
inline void add_work(allsize::Work& to,const allsize::Work& from) {
    to.events+=from.events; to.source_reads+=from.source_reads; to.target_reads+=from.target_reads;
    to.count_reads+=from.count_reads; to.updates+=from.updates; to.rank_calls+=from.rank_calls;
    to.list_moves+=from.list_moves;
}

inline Result shared(const Layout& index,Vertex n,const Combinations& choose,const Vertices& ordinary,bool rows=true) {
    auto old=commonfront::peel(index,n,choose,ordinary,rows);
    Result result;
    result.data=std::move(old.data);
    result.stats.splits=old.splits;
    result.stats.shared_events=old.shared_events;
    result.stats.shared_assignments=old.shared_assignments;
    result.stats.scalar_target_checks=old.scalar_target_checks;
    return result;
}

template<bool Audit=false>
Result components(const Graph& graph,const Layout& index,const Combinations& choose,
                  const Vertices& ordinary,bool rows=true,const std::vector<Count>* expected=nullptr) {
    const auto start=Clock::now();
    Result result;
    auto& data=result.data;
    auto& stats=result.stats;
    data.core.assign(static_cast<size_t>(index.maximum+1)*graph.n,0);
    std::copy(ordinary.begin(),ordinary.end(),data.core.begin()+2*static_cast<size_t>(graph.n));
    if (index.maximum==2) { data.peel_ms=ms(start); return result; }
    DSU dsu(graph.n);
    for (size_t p=0; p<index.paths.size(); ++p) if (index.valid(p,3)) {
        const auto row=index.paths.row(p);
        require(index.paths.holds[p]>0,"component discovery needs rooted paths");
        const size_t used=index.paths.holds[p]==3 ? 3 : row.size();
        for (size_t i=0; i<used; ++i) {
            ++stats.component_scans;
            dsu.join(row[0],row[i]);
        }
    }
    Vertices root_group(graph.n,absent),group(graph.n,absent),local(graph.n,absent);
    std::vector<Vertices> vertices,paths;
    for (Vertex v=0; v<graph.n; ++v) {
        const Vertex root=dsu.find(v);
        if (dsu.size[root]==1) continue;
        if (root_group[root]==absent) {
            root_group[root]=vertices.size();
            vertices.emplace_back(); paths.emplace_back();
        }
        const Vertex c=root_group[root];
        group[v]=c; local[v]=vertices[c].size(); vertices[c].push_back(v);
    }
    for (size_t p=0; p<index.paths.size(); ++p) if (index.hi[p]>=3) {
        const auto row=index.paths.row(p);
        const Vertex c=group[row[0]];
        require(c!=absent,"higher clique has no triangle component");
        paths[c].push_back(p);
        if constexpr (Audit) {
            const size_t used=index.paths.holds[p]==static_cast<Vertex>(index.maximum) ? index.paths.holds[p] : row.size();
            for (size_t i=0; i<used; ++i) require(group[row[i]]==c,"path crosses components");
        }
    }
    stats.components=vertices.size();
    size_t metadata=(dsu.parent.capacity()+dsu.size.capacity()+root_group.capacity()+group.capacity()+local.capacity())*sizeof(Vertex)
        +(vertices.capacity()+paths.capacity())*sizeof(Vertices);
    for (size_t c=0; c<vertices.size(); ++c) {
        stats.component_vertices+=vertices[c].size();
        stats.largest_component=std::max<uint64_t>(stats.largest_component,vertices[c].size());
        metadata+=(vertices[c].capacity()+paths[c].capacity())*sizeof(Vertex);
    }
    stats.components_ms=ms(start);
    for (size_t c=0; c<vertices.size(); ++c) {
        const auto copying=Clock::now();
        const Vertex n=vertices[c].size();
        Layout piece(Graph::from_edges(0,{}),index.maximum);
        Vertices holds,pivots,local_ordinary(n);
        for (Vertex p : paths[c]) {
            holds.clear(); pivots.clear();
            const auto row=index.paths.row(p);
            for (size_t i=0; i<index.paths.holds[p]; ++i) holds.push_back(local[row[i]]);
            if (index.paths.holds[p]<static_cast<Vertex>(index.maximum))
                for (size_t i=index.paths.holds[p]; i<row.size(); ++i) pivots.push_back(local[row[i]]);
            piece.paths.append(holds,pivots);
        }
        piece.prepare(n);
        for (Vertex v=0; v<n; ++v) local_ordinary[v]=ordinary[vertices[c][v]];
        stats.copied_index_bytes+=piece.bytes();
        stats.components_ms+=ms(copying);
        std::vector<Count> local_expected;
        if constexpr (Audit) {
            require(expected!=nullptr,"component audit needs oracle");
            local_expected.assign(static_cast<size_t>(index.maximum+1)*n,0);
            for (int s=2; s<=index.maximum; ++s) for (Vertex v=0; v<n; ++v)
                local_expected[static_cast<size_t>(s)*n+v]=(*expected)[static_cast<size_t>(s)*graph.n+vertices[c][v]];
            std::vector<std::pair<Vertex,Vertex>> edges;
            for (Vertex u : vertices[c]) for (Vertex v : graph.row(u))
                if (u<v && group[v]==c) edges.emplace_back(local[u],local[v]);
            const auto subgraph=Graph::from_edges(n,std::move(edges));
            for (int s=3; s<=index.maximum; ++s) bottomup::audit_index(subgraph,piece.paths,s);
        }
        auto part=commonfront::peel<Audit>(piece,n,choose,local_ordinary,rows,Audit ? &local_expected : nullptr);
        for (int s=3; s<=index.maximum; ++s) for (Vertex v=0; v<n; ++v)
            data.core[static_cast<size_t>(s)*graph.n+vertices[c][v]]=part.data.core[static_cast<size_t>(s)*n+v];
        add_work(data.work,part.data.work);
        stats.splits+=part.splits;
        stats.shared_events+=part.shared_events;
        stats.shared_assignments+=part.shared_assignments;
        stats.scalar_target_checks+=part.scalar_target_checks;
        stats.audits+=part.audited_states;
        data.state_bytes=std::max(data.state_bytes,metadata+piece.bytes()+part.data.state_bytes
            +part.data.core.capacity()*sizeof(Count)+(holds.capacity()+pivots.capacity()+local_ordinary.capacity())*sizeof(Vertex));
    }
    data.state_bytes=std::max(data.state_bytes,metadata);
    data.peel_ms=ms(start);
    return result;
}
}
