/*
 * SDCT_Par7: Stateless Independent Parallel BK
 *
 * Strategy: Precompute each vertex's P(v) member ordering to match Par6's
 * cumulative state, then solve each subproblem independently in parallel.
 *
 * Phase 1 (serial, fast): Simulate Par6's cumulative vertexSets swaps.
 *   For each vertex v in order 0..N-1, record which P(v) members are at
 *   positions [newBeginP..newBeginR) and their position-based ordering.
 *   Store per-vertex data: the ordered list of P(v) members.
 *
 * Phase 2 (parallel): For each vertex v, build local vertexSets/neighborsInP
 *   using the precomputed ordering, then run recurse7.
 */

#include <omp.h>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <cstdio>

#include "degeneracy_algorithm_cliques_V.h"
#include "misc.h"
#include "tree/MultiBranchTree.h"
#include "graph/DynamicGraph.h"
#include "graph/Graph.h"
#include "dataStruct/CliqueCSR.hpp"

extern double nCr[1001][401];

struct Slab7 {
    int* base = nullptr;
    int  top  = 0;
    int  cap  = 0;
    int* alloc(int n) { if (n <= 0) n = 1; int* p = base + top; top += n; return p; }
    int  save()        { return top; }
    void restore(int t){ top = t; }
};

struct Out7 {
    std::vector<int> buf;
    std::vector<int> off;
    void reserve(size_t n){ buf.reserve(n*8); off.reserve(n); }
    size_t size() const { return off.size(); }
    size_t total_vertices() const {
        size_t tv = 0;
        for (size_t i = 0; i < off.size(); i++) tv += (size_t)buf[off[i]];
        return tv;
    }
    void add(const int* keepV, int keepSz, const int* dropV, int dropSz, const int* l2g) {
        int total = keepSz + dropSz;
        off.push_back((int)buf.size());
        buf.push_back(total);
        int tmp[MAX_CSIZE];
        int ci = 0;
        for(int i = 0; i < keepSz; i++) tmp[ci++] = l2g[keepV[i]];
        if (dropV) for(int i = 0; i < dropSz; i++) tmp[ci++] = l2g[dropV[i]];
        // insertion sort by global ID
        for(int i = 1; i < ci; i++){
            int gv = tmp[i], j = i-1;
            while(j >= 0 && tmp[j] > gv){ tmp[j+1] = tmp[j]; j--; }
            tmp[j+1] = gv;
        }
        for(int i = 0; i < ci; i++) buf.push_back(tmp[i]);
    }
};

// findPivot7: uses MY_MIN + break (valid after correct ordering setup)
static int findPivot7(
    int** pivotCands, int* nPivotCands,
    int* vertexSets, int* vertexLookup,
    int** neighborsInP, int* numNeighbors,
    int beginP, int beginR, Slab7& slab)
{
    int pivot = -1, maxInP = -1;
    int pSz = beginR - beginP;

    for (int j = beginP; j < beginR; j++) {
        int u = vertexSets[j];
        int numPotential = std::min(pSz, numNeighbors[u]);
        int inP = 0;
        for (int k = 0; k < numPotential; k++) {
            int nl = vertexLookup[neighborsInP[u][k]];
            if (nl >= beginP && nl < beginR) inP++;
            else break;
        }
        if (inP > maxInP) { maxInP = inP; pivot = u; }
    }

    int* buf = slab.alloc(pSz);
    for (int i = 0; i < pSz; i++) buf[i] = vertexSets[beginP + i];
    *nPivotCands = pSz;

    int numPivotPotential = std::min(pSz, numNeighbors[pivot]);
    for (int j = 0; j < numPivotPotential; j++) {
        int nb = neighborsInP[pivot][j];
        int nl = vertexLookup[nb];
        if (nl >= beginP && nl < beginR) buf[nl - beginP] = -1;
        else break;
    }
    for (int i = 0; i < *nPivotCands; ) {
        if (buf[i] == -1) { (*nPivotCands)--; buf[i] = buf[*nPivotCands]; }
        else i++;
    }
    *pivotCands = buf;
    return pivot;
}

// moveToR7: exact match of moveToRDegeneracyCliques
static void moveToR7(
    int vertex, int* vertexSets, int* vertexLookup,
    int** neighborsInP, int* numNeighbors,
    int* pBeginP, int* pBeginR, int* pNewBeginP, int* pNewBeginR)
{
    int vl = vertexLookup[vertex];
    (*pBeginR)--;
    vertexSets[vl] = vertexSets[*pBeginR]; vertexLookup[vertexSets[*pBeginR]] = vl;
    vertexSets[*pBeginR] = vertex; vertexLookup[vertex] = *pBeginR;

    *pNewBeginP = *pBeginP;
    *pNewBeginR = *pBeginP;
    int sizeOfP = *pBeginR - *pBeginP;

    for (int j = *pBeginP; j < *pBeginR; j++) {
        int nb = vertexSets[j];
        int numPotential = std::min(sizeOfP, numNeighbors[nb]);
        for (int k = 0; k < numPotential; k++) {
            if (neighborsInP[nb][k] == vertex) {
                vertexSets[j] = vertexSets[*pNewBeginR]; vertexLookup[vertexSets[*pNewBeginR]] = j;
                vertexSets[*pNewBeginR] = nb; vertexLookup[nb] = *pNewBeginR;
                (*pNewBeginR)++;
                break;
            }
        }
    }

    for (int j = *pNewBeginP; j < *pNewBeginR; j++) {
        int tv = vertexSets[j];
        int numPotential = std::min(sizeOfP, numNeighbors[tv]);
        int cnt = 0;
        for (int k = 0; k < numPotential; k++) {
            int nb = neighborsInP[tv][k];
            int nl = vertexLookup[nb];
            if (nl >= *pNewBeginP && nl < *pNewBeginR) {
                neighborsInP[tv][k] = neighborsInP[tv][cnt];
                neighborsInP[tv][cnt++] = nb;
            }
        }
    }
}

// recurse7
static void recurse7(
    int* vertexSets, int* vertexLookup,
    int** neighborsInP, int* numNeighbors,
    int beginP, int beginR,
    int* keepV, int keepSz, int* dropV, int dropSz,
    int max_k, int min_k,
    const int* l2g,
    Slab7& slab, Out7& out)
{
    if (keepSz > max_k) return;
    if (beginP >= beginR) {
        int cSize = keepSz + dropSz;
        if (cSize < min_k) return;
        if (keepSz == max_k)
            out.add(keepV, keepSz, nullptr, 0, l2g);
        else
            out.add(keepV, keepSz, dropV, dropSz, l2g);
        return;
    }

    int* cands; int nCands;
    int arenaTop = slab.save();
    int pivot = findPivot7(&cands, &nCands, vertexSets, vertexLookup,
                           neighborsInP, numNeighbors, beginP, beginR, slab);

    for (int ci = 0; ci < nCands; ci++) {
        int vertex = cands[ci], newBeginP, newBeginR;
        moveToR7(vertex, vertexSets, vertexLookup, neighborsInP, numNeighbors,
                 &beginP, &beginR, &newBeginP, &newBeginR);
        if (vertex == pivot) {
            dropV[dropSz] = vertex;
            recurse7(vertexSets, vertexLookup, neighborsInP, numNeighbors,
                     newBeginP, newBeginR, keepV, keepSz, dropV, dropSz + 1,
                     max_k, min_k, l2g, slab, out);
        } else {
            keepV[keepSz] = vertex;
            recurse7(vertexSets, vertexLookup, neighborsInP, numNeighbors,
                     newBeginP, newBeginR, keepV, keepSz + 1, dropV, dropSz,
                     max_k, min_k, l2g, slab, out);
        }
    }
    slab.restore(arenaTop);
}

CliqueCSR<int> SDCT_Par7(Graph& G, int max_k, int min_k)
{
    int N  = (int)G.getGraphNodeSize();
    int nT = omp_get_max_threads();
    printf("SDCT_Par7: %d threads, N=%d\n", nT, N);

    long long E = (long long)G.getGraphEdgeSize();
    int arena_cap = (int)std::min(128LL*1024*1024,
                                   std::max(32LL*E/nT, (long long)N*64));
    printf("SDCT_Par7: arena %.0f MB/thread\n", arena_cap*4.0/(1024.0*1024.0));

    // ============================================================
    // Phase 1 (serial): Simulate Par6's cumulative vertexSets swaps
    // to precompute, for each vertex v, the ordered P(v) members
    // (in the exact position order that Par6's fillInP6 would produce).
    // ============================================================
    double t_pre0 = omp_get_wtime();

    // Simulate Par6's cumulative state
    std::vector<int> sim_vs(N), sim_vl(N);
    for (int i = 0; i < N; i++) { sim_vs[i] = i; sim_vl[i] = i; }

    // Per-vertex: store the ordered P(v) members (global IDs) in position order
    // pv_members[v] = list of P(v) members in the order they appear in
    //                 vertexSets[newBeginP..newBeginR) after fillInP swap
    std::vector<std::vector<int>> pv_members(N);

    int sim_beginX = 0, sim_beginP = 0, sim_beginR = N;

    for (int vertex = 0; vertex < N; vertex++) {
        // Simulate fillInP6: swap vertex into R position
        int vl = sim_vl[vertex];
        sim_beginR--;
        sim_vs[vl] = sim_vs[sim_beginR]; sim_vl[sim_vs[sim_beginR]] = vl;
        sim_vs[sim_beginR] = vertex; sim_vl[vertex] = sim_beginR;

        int newBeginR = sim_beginR;
        int newBeginP = sim_beginR;

        // Swap later neighbors of vertex into P section
        auto [nb, ne] = G.getNbr(vertex);
        for (int idx = nb; idx < ne; idx++) {
            int neighbor = G.adj_list[idx];
            if (neighbor <= vertex) continue;
            int nl = sim_vl[neighbor];
            newBeginP--;
            sim_vs[nl] = sim_vs[newBeginP]; sim_vl[sim_vs[newBeginP]] = nl;
            sim_vs[newBeginP] = neighbor; sim_vl[neighbor] = newBeginP;
        }

        // Record P(v) members in their position order
        int p = newBeginR - newBeginP;
        pv_members[vertex].resize(p);
        for (int i = 0; i < p; i++)
            pv_members[vertex][i] = sim_vs[newBeginP + i];

        // Note: we do NOT undo the swaps — cumulative state carries forward
        // (same as Par6). beginR already decremented. beginP stays at 0.
    }

    printf("SDCT_Par7: precompute %.1f ms\n", (omp_get_wtime()-t_pre0)*1000.0);

    // ============================================================
    // Phase 2 (parallel): For each vertex v, build local state
    // using precomputed pv_members ordering and run recurse7.
    // ============================================================
    std::vector<int*> slab_mem((size_t)nT);
    for(int t = 0; t < nT; t++)
        slab_mem[t] = (int*)std::malloc((size_t)arena_cap * sizeof(int));

    std::vector<Out7> thread_out((size_t)nT);
    double t0 = omp_get_wtime();

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        Slab7 slab;
        slab.base = slab_mem[tid];
        slab.cap  = arena_cap;
        slab.top  = 0;

        Out7& out = thread_out[(size_t)tid];
        out.reserve((size_t)(N/nT+1)*20);

        std::vector<int> gmark(N, -1);
        int keepV[MAX_CSIZE], dropV[MAX_CSIZE];

        #pragma omp for schedule(dynamic,16) nowait
        for (int v = 0; v < N; v++) {
            const auto& pv = pv_members[v];
            int p = (int)pv.size();

            if (p == 0) {
                if (1 >= min_k) {
                    out.off.push_back((int)out.buf.size());
                    out.buf.push_back(1);
                    out.buf.push_back(v);
                }
                continue;
            }

            int slabTop = slab.save();

            // l2g: local ID -> global ID
            // Local ID i corresponds to vertex pv[i] (in Par6's position order)
            int* l2g = slab.alloc(p + 1);
            for (int i = 0; i < p; i++) {
                l2g[i] = pv[i];
                gmark[pv[i]] = i;
            }
            l2g[p] = v;

            // vertexSets[i] = i, vertexLookup[i] = i
            int* lsets = slab.alloc(p);
            int* llook = slab.alloc(p);
            for (int i = 0; i < p; i++) { lsets[i] = i; llook[i] = i; }

            // Build neighborsInP: iterate j=0..p-1 (matching Par6's vertexSets position order)
            // For each local j, scan its adjacency for laterNeighbor > l2g[j] (global ID).
            // This matches fillInPandX's iteration exactly.

            // Count
            int* lnbCnt = slab.alloc(p);
            for (int i = 0; i < p; i++) lnbCnt[i] = 0;

            for (int j = 0; j < p; j++) {
                int gj = l2g[j];
                auto [nb2, ne2] = G.getNbr(gj);
                for (int idx2 = nb2; idx2 < ne2; idx2++) {
                    int w = G.adj_list[idx2];
                    if (w <= gj) continue;
                    int lw = gmark[w];
                    if (lw < 0) continue;
                    lnbCnt[j]++;
                    lnbCnt[lw]++;
                }
            }

            // Allocate
            std::vector<int*> lnbr_ptrs((size_t)p);
            for (int i = 0; i < p; i++)
                lnbr_ptrs[i] = slab.alloc(lnbCnt[i] > 0 ? lnbCnt[i] : 1);

            // Fill (same order as fillInPandX)
            for (int i = 0; i < p; i++) lnbCnt[i] = 0;
            for (int j = 0; j < p; j++) {
                int gj = l2g[j];
                auto [nb2, ne2] = G.getNbr(gj);
                for (int idx2 = nb2; idx2 < ne2; idx2++) {
                    int w = G.adj_list[idx2];
                    if (w <= gj) continue;
                    int lw = gmark[w];
                    if (lw < 0) continue;
                    lnbr_ptrs[j][lnbCnt[j]++] = lw;
                    lnbr_ptrs[lw][lnbCnt[lw]++] = j;
                }
            }

            keepV[0] = p;
            recurse7(lsets, llook, lnbr_ptrs.data(), lnbCnt,
                     0, p,
                     keepV, 1, dropV, 0,
                     max_k - 1, min_k - 1,
                     l2g, slab, out);

            for (int i = 0; i < p; i++) gmark[pv[i]] = -1;
            slab.restore(slabTop);
        }
    }

    for(int t = 0; t < nT; t++) std::free(slab_mem[t]);
    double t1 = omp_get_wtime();

    std::vector<size_t> coff((size_t)(nT+1), 0), voff((size_t)(nT+1), 0);
    for(int t = 0; t < nT; t++){
        coff[(size_t)(t+1)] = coff[(size_t)t] + thread_out[(size_t)t].size();
        voff[(size_t)(t+1)] = voff[(size_t)t] + thread_out[(size_t)t].total_vertices();
    }
    size_t TC = coff[(size_t)nT], TV = voff[(size_t)nT];

    std::vector<int> csrOff((size_t)(TC+1));
    std::vector<int> csrDat((size_t)TV);

    #pragma omp parallel for schedule(static) num_threads(nT)
    for(int t = 0; t < nT; t++){
        size_t base = coff[(size_t)t];
        const auto& o = thread_out[(size_t)t];
        for(size_t oi = 0; oi < o.size(); oi++)
            csrOff[base+oi] = o.buf[o.off[oi]];
    }
    size_t run = 0;
    for(size_t i = 0; i < TC; i++){ size_t s=(size_t)csrOff[i]; csrOff[i]=(int)run; run+=s; }
    csrOff[TC] = (int)run;

    #pragma omp parallel for schedule(static) num_threads(nT)
    for(int t = 0; t < nT; t++){
        size_t base = coff[(size_t)t];
        const auto& o = thread_out[(size_t)t];
        for(size_t oi = 0; oi < o.size(); oi++){
            size_t dst = (size_t)csrOff[base+oi];
            int pos = o.off[oi]; int total_sz = o.buf[pos++];
            for(int i = 0; i < total_sz; i++) csrDat[dst++] = o.buf[pos++];
        }
    }

    CliqueCSR<int> res;
    res.init_from_flat(std::move(csrOff), std::move(csrDat));
    printf("SDCT_Par7: precompute=%.1fms par=%.1fms merge=%.1fms total=%zu\n",
           (t0-t_pre0)*1000.0, (t1-t0)*1000.0, (omp_get_wtime()-t1)*1000.0, TC);
    return res;
}
