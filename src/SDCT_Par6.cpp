/*
 * SDCT_Par6: Correct high-efficiency parallel SDCT
 *
 * Correctness: uses EXACTLY the same algorithm as SDCT/SDCT_Parallel.
 * The only changes vs SDCT_Parallel:
 *   1. neighborsInP allocated from per-thread slab arena (no Calloc/Free in hot path)
 *   2. findBestPivot result stored in arena slab (no Calloc)
 *   3. Parallel two-pass flush instead of serial merge
 *   4. All per-thread memory pre-allocated serially
 */

#include <omp.h>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <iostream>

#include "degeneracy_algorithm_cliques_V.h"
#include "degeneracy_helper.h"
#include "LinkedList.h"
#include "MemoryManager.h"
#include "misc.h"
#include "tree/MultiBranchTree.h"
#include "graph/DynamicGraph.h"
#include "graph/Graph.h"
#include "dataStruct/CliqueCSR.hpp"

extern double nCr[1001][401];

// ---- Slab arena ----
struct SlabArena6 {
    int* base = nullptr;
    int  top  = 0;
    int  cap  = 0;
    int* alloc(int n)     { if (n <= 0) n = 1; int* p = base + top; top += n; return p; }
    int  save()           { return top; }
    void restore(int t)   { top = t; }
};

// ---- Leaf storage ----
struct LeafArena6 {
    std::vector<int> buf;
    std::vector<int> offsets;

    void reserve(size_t n) { buf.reserve(n * 10); offsets.reserve(n); }

    // Direct add from raw int arrays — no intermediate vector<TreeGraphNode>
    void add(const int* keepV, int keepSz, const int* dropV, int dropSz) {
        offsets.push_back((int)buf.size());
        int total = keepSz + dropSz;
        buf.push_back(total);
        // Interleave as (v, isPivot) pairs, sort inline
        int tmp[MAX_CSIZE * 2];  // v, isPivot pairs
        int cnt = 0;
        for (int i = 0; i < keepSz; i++) { tmp[cnt++] = keepV[i]; tmp[cnt++] = 0; }
        for (int i = 0; i < dropSz; i++) { tmp[cnt++] = dropV[i]; tmp[cnt++] = 1; }
        // insertion sort pairs by v value (total is small in practice)
        for (int i = 2; i < cnt; i += 2) {
            int kv = tmp[i], kp = tmp[i+1], j = i - 2;
            while (j >= 0 && tmp[j] > kv) { tmp[j+2] = tmp[j]; tmp[j+3] = tmp[j+1]; j -= 2; }
            tmp[j+2] = kv; tmp[j+3] = kp;
        }
        for (int i = 0; i < cnt; i++) buf.push_back(tmp[i]);
    }

    size_t size() const { return offsets.size(); }

    // Flush to DynamicGraph (kept for correctness testing with small graphs)
    void flush_range(DynamicGraph<TreeGraphNode>& out, size_t base_idx) {
        for (int oi = 0; oi < (int)offsets.size(); oi++) {
            int pos = offsets[oi], sz = buf[pos++];
            std::vector<TreeGraphNode> leaf; leaf.reserve((size_t)sz);
            for (int i = 0; i < sz; i++) {
                leaf.emplace_back((daf::Size)buf[pos], (bool)buf[pos+1]); pos += 2;
            }
            out.adj_list[base_idx + (size_t)oi] = std::move(leaf);
        }
    }

    // Fast flush to CSR — no per-clique heap allocation; parallel memcpy-friendly
    // out_offsets[base_clique .. base_clique+size()] and out_data[base_data..]
    // must be pre-allocated by caller (two-pass: first collect sizes, then fill)
    void fill_csr_offsets(size_t base_clique, std::vector<size_t>& out_offsets) const {
        for (int oi = 0; oi < (int)offsets.size(); oi++) {
            int pos = offsets[oi], sz = buf[pos];  // buf[pos] = total
            out_offsets[base_clique + (size_t)oi] = (size_t)sz;
        }
    }

    // After prefix-sum, fill data region starting at out_data + base_data
    void fill_csr_data(size_t base_data, std::vector<int>& out_data) const {
        size_t dst = base_data;
        for (int oi = 0; oi < (int)offsets.size(); oi++) {
            int pos = offsets[oi], sz = buf[pos++];
            // Copy only vertex IDs (every other int: v, isPivot, v, isPivot, ...)
            for (int i = 0; i < sz; i++) { out_data[dst++] = buf[pos]; pos += 2; }
        }
    }

    size_t total_vertices() const {
        // total (v,isPivot) pairs stored = (buf.size() - offsets.size()) / 2
        return (buf.size() - offsets.size()) / 2;
    }
};

// ---- findBestPivot: exact copy of findBestPivotNonNeighborsDegeneracyCliques
//      but allocates from slab instead of Calloc ----
static int findPivot6(
    int** pivotNonNeighbors, int* numNonNeighbors,
    int* vertexSets, int* vertexLookup,
    int** neighborsInP, int* numNeighbors,
    int beginP, int beginR,
    SlabArena6& arena)
{
    int pivot = -1, maxInP = -1;

    // Find pivot
    for (int j = beginP; j < beginR; j++) {
        int vertex = vertexSets[j];
        int numPotential = MY_MIN(beginR - beginP, numNeighbors[vertex]);
        int inP = 0;
        for (int k = 0; k < numPotential; k++) {
            int nl = vertexLookup[neighborsInP[vertex][k]];
            if (nl >= beginP && nl < beginR) inP++;
            else break;
        }
        if (inP > maxInP) { maxInP = inP; pivot = vertex; }
    }

    // Collect non-neighbors of pivot using slab (not Calloc)
    int sz = beginR - beginP;
    int* buf = arena.alloc(sz);
    // copy P into buf
    for (int i = 0; i < sz; i++) buf[i] = vertexSets[beginP + i];
    *numNonNeighbors = sz;

    int numPivotNbr = MY_MIN(sz, numNeighbors[pivot]);
    // mark pivot's neighbors by setting to -1
    for (int j = 0; j < numPivotNbr; j++) {
        int nb = neighborsInP[pivot][j];
        int nl = vertexLookup[nb];
        if (nl >= beginP && nl < beginR) {
            buf[nl - beginP] = -1;
        } else break;
    }
    // compact: remove -1 entries
    int j = 0;
    while (j < *numNonNeighbors) {
        if (buf[j] == -1) {
            (*numNonNeighbors)--;
            buf[j] = buf[*numNonNeighbors];
        } else {
            j++;
        }
    }

    *pivotNonNeighbors = buf;
    return pivot;
}

// ---- moveToR: exact copy of moveToRDegeneracyCliques ----
static void moveToR6(
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

    // Find neighbors of vertex in P
    for (int j = *pBeginP; j < *pBeginR; j++) {
        int neighbor = vertexSets[j];
        int numPotential = MY_MIN(sizeOfP, numNeighbors[neighbor]);
        for (int k = 0; k < numPotential; k++) {
            if (neighborsInP[neighbor][k] == vertex) {
                vertexSets[j] = vertexSets[*pNewBeginR];
                vertexLookup[vertexSets[*pNewBeginR]] = j;
                vertexSets[*pNewBeginR] = neighbor;
                vertexLookup[neighbor] = *pNewBeginR;
                (*pNewBeginR)++;
                break;
            }
        }
    }

    // Compact neighborsInP for new P
    for (int j = *pNewBeginP; j < *pNewBeginR; j++) {
        int tv = vertexSets[j];
        int numPotential = MY_MIN(sizeOfP, numNeighbors[tv]);
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

// ---- fillInP: exact copy of fillInPandXForRecursiveCallDegeneracyCliquesEdgeGraph
//      but uses slab instead of Calloc ----
static void fillInP6(
    int vertex, int* vertexSets, int* vertexLookup,
    Graph& edgeGraph, int** neighborsInP, int* numNeighbors,
    int* pBeginX, int* pBeginP, int* pBeginR,
    int* pNewBeginX, int* pNewBeginP, int* pNewBeginR,
    SlabArena6& arena)
{
    int vl = vertexLookup[vertex];
    (*pBeginR)--;
    vertexSets[vl] = vertexSets[*pBeginR]; vertexLookup[vertexSets[*pBeginR]] = vl;
    vertexSets[*pBeginR] = vertex; vertexLookup[vertex] = *pBeginR;
    *pNewBeginR = *pBeginR; *pNewBeginP = *pBeginR;

    auto [nb, ne] = edgeGraph.getNbr(vertex);
    for (int idx = nb; idx < ne; ++idx) {
        int neighbor = edgeGraph.adj_list[idx];
        if (neighbor <= vertex) continue;
        int nl = vertexLookup[neighbor]; (*pNewBeginP)--;
        vertexSets[nl] = vertexSets[*pNewBeginP]; vertexLookup[vertexSets[*pNewBeginP]] = nl;
        vertexSets[*pNewBeginP] = neighbor; vertexLookup[neighbor] = *pNewBeginP;
    }
    *pNewBeginX = *pNewBeginP;
    int pSize = *pNewBeginR - *pNewBeginP;

    // Allocate neighborsInP from slab instead of Calloc
    for (int j = *pNewBeginP; j < *pNewBeginR; j++) {
        int v = vertexSets[j]; numNeighbors[v] = 0;
        int as = MY_MIN(pSize, (int)edgeGraph.getNbrCount(v));
        neighborsInP[v] = arena.alloc(as > 0 ? as : 1);
    }
    for (int j = *pNewBeginP; j < *pNewBeginR; j++) {
        int v = vertexSets[j];
        auto [nb2, ne2] = edgeGraph.getNbr(v);
        for (int idx = nb2; idx < ne2; ++idx) {
            int w = edgeGraph.adj_list[idx];
            if (w <= v) continue;
            int wl = vertexLookup[w];
            if (wl >= *pNewBeginP && wl < *pNewBeginR) {
                neighborsInP[v][numNeighbors[v]++] = w;
                neighborsInP[w][numNeighbors[w]++] = v;
            }
        }
    }
}

// ---- recurse6: exact port of listAllCliquesDegeneracyRecursive_VedgeGraphSDCT
//      using slab for pivot candidates ----
static void recurse6(
    int* vertexSets, int* vertexLookup,
    int** neighborsInP, int* numNeighbors,
    int beginP, int beginR,
    int* keepV, int keepSz, int* dropV, int dropSz,
    int max_k, int min_k,
    SlabArena6& arena, LeafArena6& output)
{
    if (keepSz > max_k) return;
    if (beginP >= beginR) {
        int cSize = keepSz + dropSz;
        if (cSize < min_k) return;
        if (keepSz == max_k)
            output.add(keepV, keepSz, nullptr, 0);
        else
            output.add(keepV, keepSz, dropV, dropSz);
        return;
    }

    int* cands; int nCands;
    int arenaTop = arena.save();
    int pivot = findPivot6(&cands, &nCands, vertexSets, vertexLookup,
                           neighborsInP, numNeighbors, beginP, beginR, arena);

    for (int ci = 0; ci < nCands; ci++) {
        int vertex = cands[ci], newBeginP, newBeginR;
        moveToR6(vertex, vertexSets, vertexLookup, neighborsInP, numNeighbors,
                 &beginP, &beginR, &newBeginP, &newBeginR);
        if (vertex == pivot) {
            dropV[dropSz] = vertex;
            recurse6(vertexSets, vertexLookup, neighborsInP, numNeighbors,
                     newBeginP, newBeginR, keepV, keepSz, dropV, dropSz + 1,
                     max_k, min_k, arena, output);
        } else {
            keepV[keepSz] = vertex;
            recurse6(vertexSets, vertexLookup, neighborsInP, numNeighbors,
                     newBeginP, newBeginR, keepV, keepSz + 1, dropV, dropSz,
                     max_k, min_k, arena, output);
        }
    }
    arena.restore(arenaTop);
}

// ============================================================
// Public entry point
// ============================================================
DynamicGraph<TreeGraphNode> SDCT_Par6(Graph& edgeGraph, int max_k, int min_k) {
    int size     = (int)edgeGraph.getGraphNodeSize();
    int nthreads = omp_get_max_threads();
    printf("SDCT_Par6: %d threads\n", nthreads);

    long long total_edges = (long long)edgeGraph.getGraphEdgeSize();
    int arena_cap = (int)std::min(
        512LL * 1024 * 1024,
        std::max(8LL * total_edges, (long long)size * 16));
    printf("SDCT_Par6: arena_cap=%d ints/thread (%.0f MB/thread)\n",
           arena_cap, (double)arena_cap * 4.0 / (1024.0*1024.0));

    struct ThreadBufs {
        int*  vertexSets   = nullptr;
        int*  vertexLookup = nullptr;
        int*  numNeighbors = nullptr;
        int** neighborsInP = nullptr;
        int*  arena_mem    = nullptr;
    };
    std::vector<ThreadBufs> bufs((size_t)nthreads);
    double t0 = omp_get_wtime();
    for (int t = 0; t < nthreads; t++) {
        bufs[t].vertexSets   = (int*) std::malloc((size_t)size * sizeof(int));
        bufs[t].vertexLookup = (int*) std::malloc((size_t)size * sizeof(int));
        bufs[t].numNeighbors = (int*) std::malloc((size_t)size * sizeof(int));
        bufs[t].neighborsInP = (int**)std::malloc((size_t)size * sizeof(int*));
        bufs[t].arena_mem    = (int*) std::malloc((size_t)arena_cap * sizeof(int));
    }
    printf("SDCT_Par6: pre-alloc %.1f ms\n", (omp_get_wtime()-t0)*1000.0);

    std::vector<LeafArena6> thread_leaves((size_t)nthreads);

    double t_par0 = omp_get_wtime();
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();

        SlabArena6 arena;
        arena.base = bufs[tid].arena_mem;
        arena.cap  = arena_cap;
        arena.top  = 0;

        int*  vertexSets   = bufs[tid].vertexSets;
        int*  vertexLookup = bufs[tid].vertexLookup;
        int*  numNeighbors = bufs[tid].numNeighbors;
        int** neighborsInP = bufs[tid].neighborsInP;

        // First-touch init: same as SDCT_Parallel
        for (int i = 0; i < size; i++) {
            vertexSets[i]   = i;
            vertexLookup[i] = i;
            numNeighbors[i] = 1;
            neighborsInP[i] = arena.alloc(1);  // dummy 1-int slot (never read before overwrite)
        }
        // arena.top == size now

        LeafArena6& output = thread_leaves[(size_t)tid];
        output.reserve((size_t)std::max(1, size / nthreads) * 20);

        #pragma omp barrier

        // Each thread maintains cumulative state (same as SDCT_Parallel).
        int beginX = 0, beginP = 0, beginR = size;
        int keepV[MAX_CSIZE], dropV[MAX_CSIZE];
        double t_work = 0.0, t_ps = omp_get_wtime();

        #pragma omp for schedule(guided) nowait
        for (int vertex = 0; vertex < size; vertex++) {
            double tw0 = omp_get_wtime();

            // Save arena top: restore after each vertex to reclaim slab memory.
            // vertexSets/vertexLookup retain modifications (cumulative state).
            int arenaBase = arena.save();

            int newBeginX, newBeginP, newBeginR;
            fillInP6(vertex, vertexSets, vertexLookup, edgeGraph,
                     neighborsInP, numNeighbors,
                     &beginX, &beginP, &beginR,
                     &newBeginX, &newBeginP, &newBeginR, arena);

            keepV[0] = vertex;
            recurse6(vertexSets, vertexLookup, neighborsInP, numNeighbors,
                     newBeginP, newBeginR,
                     keepV, 1, dropV, 0,
                     max_k, min_k, arena, output);

            arena.restore(arenaBase);
            beginR = beginR + 1;  // same as serial SDCT

            t_work += omp_get_wtime() - tw0;
        }

        double t_pe = omp_get_wtime(), dur = t_pe - t_ps;
        #pragma omp critical
        printf("Thread %2d: work=%.1fms total=%.1fms eff=%.0f%% cliques=%zu\n",
               tid, t_work*1000.0, dur*1000.0,
               dur > 0 ? 100.0*t_work/dur : 0.0, output.size());
    }
    printf("SDCT_Par6: parallel %.1f ms\n", (omp_get_wtime()-t_par0)*1000.0);

    for (int t = 0; t < nthreads; t++) {
        std::free(bufs[t].vertexSets);   std::free(bufs[t].vertexLookup);
        std::free(bufs[t].numNeighbors); std::free(bufs[t].neighborsInP);
        std::free(bufs[t].arena_mem);
    }

    // --- Parallel CSR merge (no per-clique heap allocation) ---
    double t_m0 = omp_get_wtime();

    // Pass 1: count cliques and vertices per thread
    std::vector<size_t> clique_off((size_t)(nthreads+1), 0);
    std::vector<size_t> vert_off((size_t)(nthreads+1), 0);
    for (int t = 0; t < nthreads; t++) {
        clique_off[(size_t)(t+1)] = clique_off[(size_t)t] + thread_leaves[(size_t)t].size();
        vert_off[(size_t)(t+1)]   = vert_off[(size_t)t]   + thread_leaves[(size_t)t].total_vertices();
    }
    size_t total_cliques  = clique_off[(size_t)nthreads];
    size_t total_vertices = vert_off[(size_t)nthreads];

    // Allocate flat CSR arrays
    // offset_[i] = start in data_ for clique i; offset_[total_cliques] = total_vertices
    std::vector<size_t> csr_offsets(total_cliques + 1);
    std::vector<int>    csr_data(total_vertices);

    // Pass 2 (parallel): fill per-clique sizes into csr_offsets[base..base+n]
    #pragma omp parallel for schedule(static) num_threads(nthreads)
    for (int t = 0; t < nthreads; t++) {
        size_t base = clique_off[(size_t)t];
        const auto& la = thread_leaves[(size_t)t];
        for (size_t oi = 0; oi < la.size(); oi++) {
            int pos = la.offsets[oi];
            int sz  = la.buf[pos];  // stored as first int per entry
            csr_offsets[base + oi] = (size_t)sz;
        }
    }

    // Serial prefix sum over csr_offsets (fast: one pass over total_cliques+1 ints)
    size_t running = 0;
    for (size_t i = 0; i <= total_cliques; i++) {
        size_t sz = csr_offsets[i];
        csr_offsets[i] = running;
        running += sz;
    }

    // Pass 3 (parallel): fill vertex data using pre-computed offsets
    #pragma omp parallel for schedule(static) num_threads(nthreads)
    for (int t = 0; t < nthreads; t++) {
        size_t base_clique = clique_off[(size_t)t];
        const auto& la = thread_leaves[(size_t)t];
        for (size_t oi = 0; oi < la.size(); oi++) {
            size_t dst = csr_offsets[base_clique + oi];
            int pos    = la.offsets[oi];
            int sz     = la.buf[pos++];
            // buf layout: sz, v0, p0, v1, p1, ...  — only copy vertex IDs
            for (int i = 0; i < sz; i++) { csr_data[dst++] = la.buf[pos]; pos += 2; }
        }
    }

    double t_merge_ms = (omp_get_wtime() - t_m0) * 1000.0;
    printf("SDCT_Par6: CSR merge %.1f ms | total cliques: %zu | total verts: %zu | wall %.1f ms\n",
           t_merge_ms, total_cliques, total_vertices, (omp_get_wtime()-t_par0)*1000.0);

    // Build DynamicGraph from CSR (needed for callers that consume DynamicGraph)
    double t_dg0 = omp_get_wtime();
    DynamicGraph<TreeGraphNode> treeGraph(size);
    treeGraph.adj_list.resize(total_cliques);
    #pragma omp parallel for schedule(static) num_threads(nthreads)
    for (int t = 0; t < nthreads; t++) {
        size_t base = clique_off[(size_t)t];
        const auto& la = thread_leaves[(size_t)t];
        for (size_t oi = 0; oi < la.size(); oi++) {
            size_t ci  = base + oi;
            size_t src = csr_offsets[ci];
            size_t sz  = csr_offsets[ci+1] - src;
            // Reconstruct with isPivot info from original buf
            int pos = la.offsets[oi]; pos++;  // skip sz field
            std::vector<TreeGraphNode> leaf; leaf.reserve(sz);
            for (size_t i = 0; i < sz; i++) {
                leaf.emplace_back((daf::Size)la.buf[pos], (bool)la.buf[pos+1]); pos += 2;
            }
            treeGraph.adj_list[ci] = std::move(leaf);
        }
    }
    printf("SDCT_Par6: DynamicGraph build %.1f ms\n", (omp_get_wtime()-t_dg0)*1000.0);

    return treeGraph;
}

// ============================================================
// CSR-only entry point — returns CliqueCSR<int>, no DynamicGraph
// Avoids ALL per-clique vector<> heap allocation in merge.
// ============================================================
CliqueCSR<int> SDCT_Par6_CSR(Graph& edgeGraph, int max_k, int min_k) {
    int size     = (int)edgeGraph.getGraphNodeSize();
    int nthreads = omp_get_max_threads();
    printf("SDCT_Par6_CSR: %d threads\n", nthreads);

    long long total_edges = (long long)edgeGraph.getGraphEdgeSize();
    int arena_cap = (int)std::min(
        512LL * 1024 * 1024,
        std::max(8LL * total_edges, (long long)size * 16));

    struct ThreadBufs {
        int*  vertexSets   = nullptr;
        int*  vertexLookup = nullptr;
        int*  numNeighbors = nullptr;
        int** neighborsInP = nullptr;
        int*  arena_mem    = nullptr;
    };
    std::vector<ThreadBufs> bufs((size_t)nthreads);
    for (int t = 0; t < nthreads; t++) {
        bufs[t].vertexSets   = (int*) std::malloc((size_t)size * sizeof(int));
        bufs[t].vertexLookup = (int*) std::malloc((size_t)size * sizeof(int));
        bufs[t].numNeighbors = (int*) std::malloc((size_t)size * sizeof(int));
        bufs[t].neighborsInP = (int**)std::malloc((size_t)size * sizeof(int*));
        bufs[t].arena_mem    = (int*) std::malloc((size_t)arena_cap * sizeof(int));
    }

    std::vector<LeafArena6> thread_leaves((size_t)nthreads);

    double t_par0 = omp_get_wtime();
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        SlabArena6 arena;
        arena.base = bufs[tid].arena_mem;
        arena.cap  = arena_cap;
        arena.top  = 0;

        int*  vertexSets   = bufs[tid].vertexSets;
        int*  vertexLookup = bufs[tid].vertexLookup;
        int*  numNeighbors = bufs[tid].numNeighbors;
        int** neighborsInP = bufs[tid].neighborsInP;

        for (int i = 0; i < size; i++) {
            vertexSets[i]   = i;
            vertexLookup[i] = i;
            numNeighbors[i] = 1;
            neighborsInP[i] = arena.alloc(1);
        }

        LeafArena6& output = thread_leaves[(size_t)tid];
        output.reserve((size_t)std::max(1, size / nthreads) * 20);

        #pragma omp barrier

        int beginX = 0, beginP = 0, beginR = size;
        int keepV[MAX_CSIZE], dropV[MAX_CSIZE];
        double t_work = 0.0, t_ps = omp_get_wtime();

        #pragma omp for schedule(guided) nowait
        for (int vertex = 0; vertex < size; vertex++) {
            double tw0 = omp_get_wtime();
            int arenaBase = arena.save();
            int newBeginX, newBeginP, newBeginR;
            fillInP6(vertex, vertexSets, vertexLookup, edgeGraph,
                     neighborsInP, numNeighbors,
                     &beginX, &beginP, &beginR,
                     &newBeginX, &newBeginP, &newBeginR, arena);
            keepV[0] = vertex;
            recurse6(vertexSets, vertexLookup, neighborsInP, numNeighbors,
                     newBeginP, newBeginR,
                     keepV, 1, dropV, 0,
                     max_k, min_k, arena, output);
            arena.restore(arenaBase);
            beginR = beginR + 1;
            t_work += omp_get_wtime() - tw0;
        }

        double dur = omp_get_wtime() - t_ps;
        #pragma omp critical
        printf("Thread %2d: work=%.1fms total=%.1fms eff=%.0f%% cliques=%zu\n",
               tid, t_work*1000.0, dur*1000.0,
               dur > 0 ? 100.0*t_work/dur : 0.0, output.size());
    }
    printf("SDCT_Par6_CSR: parallel %.1f ms\n", (omp_get_wtime()-t_par0)*1000.0);

    for (int t = 0; t < nthreads; t++) {
        std::free(bufs[t].vertexSets);   std::free(bufs[t].vertexLookup);
        std::free(bufs[t].numNeighbors); std::free(bufs[t].neighborsInP);
        std::free(bufs[t].arena_mem);
    }

    // Parallel three-pass CSR merge
    double t_m0 = omp_get_wtime();
    std::vector<size_t> clique_off((size_t)(nthreads+1), 0);
    std::vector<size_t> vert_off((size_t)(nthreads+1), 0);
    for (int t = 0; t < nthreads; t++) {
        clique_off[(size_t)(t+1)] = clique_off[(size_t)t] + thread_leaves[(size_t)t].size();
        vert_off[(size_t)(t+1)]   = vert_off[(size_t)t]   + thread_leaves[(size_t)t].total_vertices();
    }
    size_t total_cliques  = clique_off[(size_t)nthreads];
    size_t total_vertices = vert_off[(size_t)nthreads];

    CliqueCSR<int> result;
    // Pre-allocate: offset_ needs total_cliques+1 entries, data_ needs total_vertices
    result.reserve_cliques(total_cliques);
    result.reserve_vertices(total_vertices);

    // We build CSR directly: fill offsets array then data array in parallel
    // Use internal vectors via a two-pass: accumulate sizes, prefix-sum, fill data
    std::vector<size_t> csr_offsets(total_cliques + 1);

    #pragma omp parallel for schedule(static) num_threads(nthreads)
    for (int t = 0; t < nthreads; t++) {
        size_t base = clique_off[(size_t)t];
        const auto& la = thread_leaves[(size_t)t];
        for (size_t oi = 0; oi < la.size(); oi++) {
            csr_offsets[base + oi] = (size_t)la.buf[la.offsets[oi]];  // sz field
        }
    }
    size_t running = 0;
    for (size_t i = 0; i < total_cliques; i++) {
        size_t sz = csr_offsets[i];
        csr_offsets[i] = running;
        running += sz;
    }
    csr_offsets[total_cliques] = running;

    std::vector<int> csr_data(total_vertices);
    #pragma omp parallel for schedule(static) num_threads(nthreads)
    for (int t = 0; t < nthreads; t++) {
        size_t base = clique_off[(size_t)t];
        const auto& la = thread_leaves[(size_t)t];
        for (size_t oi = 0; oi < la.size(); oi++) {
            size_t dst = csr_offsets[base + oi];
            int pos    = la.offsets[oi];
            int sz     = la.buf[pos++];
            for (int i = 0; i < sz; i++) { csr_data[dst++] = la.buf[pos]; pos += 2; }
        }
    }

    // Build CliqueCSR from the flat arrays (single memcpy equivalent)
    for (size_t i = 0; i < total_cliques; i++) {
        result.add_clique(csr_data.data() + csr_offsets[i],
                          csr_offsets[i+1] - csr_offsets[i]);
    }

    printf("SDCT_Par6_CSR: merge %.1f ms | total cliques: %zu | wall %.1f ms\n",
           (omp_get_wtime()-t_m0)*1000.0, total_cliques, (omp_get_wtime()-t_par0)*1000.0);

    return result;
}
