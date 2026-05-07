// Implementation of the post-peel (1,s)-nucleus hierarchy build.
//
// Memory note: V3 peel itself runs at the §7.3 ~10Σ peel-state budget;
// nothing here is allocated until peel completes.  We then take O(Σ)
// scratch to build the leaf-first vertex lists (one bucket-CSR pass
// over vtxLeafIds), run union-find with leaves as connector nodes, and
// free the scratch on return.

#include "BuildHierarchyR1.h"
#include "NCliqueCoreDecomposition.h"  // ST_V2_Data

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

namespace nucleus_hier {

namespace {

struct HierNode {
    int     id;
    double  k_birth;
    double  k_death;
    int     parent;
    int     size_birth;
    int     size_death;
};

// Path-compressed find. parent_uf[r] == r marks a root.
inline daf::Size dsu_find(std::vector<daf::Size>& parent_uf, daf::Size x) {
    while (parent_uf[x] != x) {
        parent_uf[x] = parent_uf[parent_uf[x]];
        x = parent_uf[x];
    }
    return x;
}

}  // namespace

void buildAndDumpHierarchyFromCSR(
    const ST_V2_Data& v3,
    const double*     coreV,
    daf::Size         numVertices,
    daf::Size         s,
    int               min_size,
    const char*       out_path)
{
    if (!coreV || !out_path) return;
    const size_t numLeaves = v3.numLeaves;
    if (numLeaves == 0) {
        std::cerr << "PIVOTER_DUMP_HIER: no leaves\n";
        return;
    }

    // ---- Step 1: bucket-CSR pass.  Compute kP_L and per-leaf vertex
    //              counts in one scan over vtxLeafIds.  Then build
    //              leaf-first CSR (leafFirstOff / leafFirstIds) so we
    //              can list every leaf's vertices in O(|leaf|).
    std::vector<double>     kP(numLeaves, std::numeric_limits<double>::infinity());
    std::vector<daf::Size>  leafCount(numLeaves, 0);

    auto vtxLeafIdAt = [&](size_t i) -> daf::Size { return v3.vtxLeafIds[i]; };

    for (daf::Size v = 0; v < numVertices; ++v) {
        const double cv = coreV[v];
        const size_t off0 = v3.vtxLeafOff[v];
        const size_t off1 = v3.vtxLeafOff[v + 1];
        for (size_t i = off0; i < off1; ++i) {
            daf::Size L = vtxLeafIdAt(i);
            ++leafCount[L];
            if (cv >= 0.0 && cv < kP[L]) kP[L] = cv;
        }
    }

    // Prefix-sum to leafFirstOff; reuse leafCount as fill pointers.
    std::vector<size_t> leafFirstOff(numLeaves + 1, 0);
    for (size_t L = 0; L < numLeaves; ++L) leafFirstOff[L + 1] = leafFirstOff[L] + leafCount[L];
    const size_t totalIncidences = leafFirstOff[numLeaves];
    std::vector<daf::Size> leafFirstIds(totalIncidences);
    {
        std::vector<size_t> fillPos(leafFirstOff.begin(), leafFirstOff.end() - 1);
        for (daf::Size v = 0; v < numVertices; ++v) {
            const size_t off0 = v3.vtxLeafOff[v];
            const size_t off1 = v3.vtxLeafOff[v + 1];
            for (size_t i = off0; i < off1; ++i) {
                daf::Size L = vtxLeafIdAt(i);
                leafFirstIds[fillPos[L]++] = v;
            }
        }
    }

    // ---- Step 2: filter leaves & sort by descending kP ----
    // Skip leaves with size < s (can't form an s-clique) or kP <= 0
    // (no admissible level).  Stable order-by-leafId tie-break.
    struct LeafEv { double kP; daf::Size L; };
    std::vector<LeafEv> events;
    events.reserve(numLeaves);
    for (size_t L = 0; L < numLeaves; ++L) {
        if (leafCount[L] < (daf::Size)s) continue;
        if (!std::isfinite(kP[L]) || kP[L] <= 0) continue;
        events.push_back({kP[L], (daf::Size)L});
    }
    std::sort(events.begin(), events.end(),
              [](const LeafEv& a, const LeafEv& b) {
                  if (a.kP != b.kP) return a.kP > b.kP;
                  return a.L < b.L;
              });

    // ---- Step 3: DSU over (vertices ∪ leaves).  Leaves act as
    //              connector nodes — union(v, leaf_id) means "v is in
    //              leaf L".  Two vertices share a DSU root iff they
    //              are s-connected through admissible CPI cliques.
    const daf::Size dsuSize = numVertices + (daf::Size)numLeaves;
    auto leafSlot = [&](daf::Size L) -> daf::Size { return numVertices + L; };

    std::vector<daf::Size> parent_uf(dsuSize);
    std::iota(parent_uf.begin(), parent_uf.end(), (daf::Size)0);
    std::vector<int>       size_uf(dsuSize, 0);   // 0 = not yet active
    std::vector<int>       state_node(dsuSize, -1);
    std::vector<HierNode>  nodes;
    nodes.reserve(events.size());

    // Reusable scratch for distinct-roots collection.
    std::vector<daf::Size> distinct_reps;
    distinct_reps.reserve(64);
    std::vector<uint8_t>   seen_mark(dsuSize, 0);
    std::vector<daf::Size> seen_clear;
    seen_clear.reserve(64);

    auto count_vertex_in_component = [&](daf::Size root) {
        // size_uf[root] tracks DSU size including leaves.  For the
        // size_birth / size_death output we want VERTEX count — compute
        // by walking the DSU at the end (cheaper than re-tracking
        // separately on every union).  Here we just defer.
        (void)root;
        return 0;  // placeholder; populated in a final pass
    };
    (void)count_vertex_in_component;

    // ---- Step 4: process leaves in descending birth order ----
    for (const auto& e : events) {
        const daf::Size L = e.L;
        const double    kP_L = e.kP;
        const size_t off0 = leafFirstOff[L];
        const size_t off1 = leafFirstOff[L + 1];

        // Collect distinct DSU roots among leaf members + the leaf
        // slot itself.  Activate any new singletons (size_uf == 0).
        distinct_reps.clear();
        for (daf::Size x : seen_clear) seen_mark[x] = 0;
        seen_clear.clear();

        auto add_rep = [&](daf::Size x) {
            if (size_uf[x] == 0) size_uf[x] = 1;
            daf::Size r = dsu_find(parent_uf, x);
            if (!seen_mark[r]) {
                seen_mark[r] = 1;
                seen_clear.push_back(r);
                distinct_reps.push_back(r);
            }
        };
        add_rep(leafSlot(L));
        for (size_t i = off0; i < off1; ++i) add_rep(leafFirstIds[i]);
        if (distinct_reps.empty()) continue;

        // Emit a hierarchy node for any rep that doesn't yet have one.
        struct Branch { int nid; double kb; int sz; daf::Size rep; };
        std::vector<Branch> branches;
        branches.reserve(distinct_reps.size());
        for (daf::Size r : distinct_reps) {
            if (state_node[r] == -1) {
                HierNode n{};
                n.id        = (int)nodes.size();
                n.k_birth   = kP_L;
                n.k_death   = 0;
                n.parent    = -1;
                n.size_birth = size_uf[r];
                n.size_death = size_uf[r];
                nodes.push_back(n);
                state_node[r] = n.id;
            }
            branches.push_back({state_node[r], nodes[state_node[r]].k_birth,
                                size_uf[r], r});
        }

        // Union all reps.
        daf::Size new_root = distinct_reps[0];
        for (size_t i = 1; i < distinct_reps.size(); ++i) {
            daf::Size rA = dsu_find(parent_uf, new_root);
            daf::Size rB = dsu_find(parent_uf, distinct_reps[i]);
            if (rA == rB) continue;
            if (size_uf[rA] >= size_uf[rB]) {
                parent_uf[rB] = rA;
                size_uf[rA] += size_uf[rB];
                new_root = rA;
            } else {
                parent_uf[rA] = rB;
                size_uf[rB] += size_uf[rA];
                new_root = rB;
            }
        }
        new_root = dsu_find(parent_uf, new_root);

        // Elder rule among the merged branches.
        if (branches.size() >= 2) {
            std::sort(branches.begin(), branches.end(),
                      [](const Branch& a, const Branch& b) {
                          if (a.kb != b.kb) return a.kb > b.kb;
                          return a.sz > b.sz;
                      });
            int elder = branches[0].nid;
            for (size_t i = 1; i < branches.size(); ++i) {
                int c = branches[i].nid;
                nodes[c].k_death    = kP_L;
                nodes[c].size_death = branches[i].sz;
                nodes[c].parent     = elder;
            }
            state_node[new_root] = elder;
        } else {
            state_node[new_root] = branches[0].nid;
        }
    }

    // ---- Final pass: roots get size_death = vertex count of their DSU
    //                  component (excluding leaf nodes, which are
    //                  bookkeeping artifacts not part of paper Def 4).
    std::vector<int> rootVertexCount(dsuSize, 0);
    for (daf::Size v = 0; v < numVertices; ++v) {
        if (size_uf[v] == 0) continue;          // never touched a leaf
        daf::Size r = dsu_find(parent_uf, v);
        ++rootVertexCount[r];
    }
    for (auto& n : nodes) {
        if (n.parent != -1) continue;            // not a root
        // find any rep mapping to this root via state_node lookup
        // (linear scan is fine — #nodes is small relative to Σ)
    }
    // simpler: walk roots once, assign root size to its hierarchy node
    for (daf::Size r = 0; r < dsuSize; ++r) {
        if (parent_uf[r] != r) continue;
        int nid = state_node[r];
        if (nid < 0) continue;
        if (nodes[nid].parent != -1) continue;
        nodes[nid].size_death = rootVertexCount[r];
    }

    // ---- Write CSV ----
    std::FILE* f = std::fopen(out_path, "w");
    if (!f) {
        std::cerr << "PIVOTER_DUMP_HIER: failed to open " << out_path << "\n";
        return;
    }
    std::fprintf(f, "id,k_birth,k_death,parent,size_birth,size_death,persistence\n");
    int written = 0;
    for (const auto& n : nodes) {
        if (n.size_birth < min_size && n.size_death < min_size) continue;
        const double persistence = n.k_birth - n.k_death;
        std::fprintf(f, "%d,%.0f,%.0f,%d,%d,%d,%.0f\n",
                     n.id, n.k_birth, n.k_death, n.parent,
                     n.size_birth, n.size_death, persistence);
        ++written;
    }
    std::fclose(f);
    std::cerr << "PIVOTER_DUMP_HIER: wrote " << written
              << " hierarchy nodes (" << nodes.size() << " total) to "
              << out_path << "\n";
}

}  // namespace nucleus_hier
