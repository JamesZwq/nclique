/*
 * SDCT_Par7: Stateless Independent Parallel BK
 *
 * Each vertex v gets a fully independent subproblem:
 *   P(v) = {u > v : u ∈ N(v)}  (global IDs, degeneracy order)
 *   Root v is always in the partial clique.
 *   No shared state between vertices → perfect parallelism.
 *
 * Correctness: cliqueCount (combinatorial k-clique counts via nCr) is
 * pivot-invariant, so any valid BK pivot ordering is correct.
 * No serial precompute phase needed.
 *
 * Output stores isPivot flag per vertex for cliqueCount computation.
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

// Per-thread output buffer storing vertex IDs and isPivot flags
struct Out7 {
    std::vector<int> buf;       // [total, v0, v1, ..., total, v0, ...]
    std::vector<uint8_t> pbuf;  // isPivot flags parallel to vertex entries in buf
    std::vector<int> off;       // offset[i] = position in buf where clique i starts
    void reserve(size_t n){ buf.reserve(n*8); pbuf.reserve(n*8); off.reserve(n); }
    size_t size() const { return off.size(); }
    size_t total_vertices() const {
        size_t tv = 0;
        for (size_t i = 0; i < off.size(); i++) tv += (size_t)buf[off[i]];
        return tv;
    }
    // Add a clique leaf with keepV (non-pivot) and dropV (pivot) vertices
    // l2g maps local IDs to global IDs
    void add(const int* keepV, int keepSz, const int* dropV, int dropSz, const int* l2g) {
        int total = keepSz + dropSz;
        off.push_back((int)buf.size());
        buf.push_back(total);  // size header (no isPivot for this entry)
        // Build sorted (globalID, isPivot) pairs
        struct VP { int gid; uint8_t pivot; };
        VP tmp[MAX_CSIZE];
        int ci = 0;
        for(int i = 0; i < keepSz; i++) { tmp[ci++] = {l2g[keepV[i]], 0}; }
        if (dropV) for(int i = 0; i < dropSz; i++) { tmp[ci++] = {l2g[dropV[i]], 1}; }
        // insertion sort by global ID
        for(int i = 1; i < ci; i++){
            VP v = tmp[i]; int j = i-1;
            while(j >= 0 && tmp[j].gid > v.gid){ tmp[j+1] = tmp[j]; j--; }
            tmp[j+1] = v;
        }
        for(int i = 0; i < ci; i++) {
            buf.push_back(tmp[i].gid);
            pbuf.push_back(tmp[i].pivot);
        }
    }
};

// findPivot7: full scan (Par7 has bidirectional neighbor lists, no packing invariant)
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
        int inP = 0;
        for (int k = 0; k < numNeighbors[u]; k++) {
            int nl = vertexLookup[neighborsInP[u][k]];
            if (nl >= beginP && nl < beginR) inP++;
        }
        if (inP > maxInP) { maxInP = inP; pivot = u; }
    }

    int* buf = slab.alloc(pSz);
    for (int i = 0; i < pSz; i++) buf[i] = vertexSets[beginP + i];
    *nPivotCands = pSz;

    for (int j = 0; j < numNeighbors[pivot]; j++) {
        int nb = neighborsInP[pivot][j];
        int nl = vertexLookup[nb];
        if (nl >= beginP && nl < beginR) buf[nl - beginP] = -1;
    }
    for (int i = 0; i < *nPivotCands; ) {
        if (buf[i] == -1) { (*nPivotCands)--; buf[i] = buf[*nPivotCands]; }
        else i++;
    }
    *pivotCands = buf;
    return pivot;
}

// moveToR7: full scan (no early break, no MIN bound)
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

    for (int j = *pBeginP; j < *pBeginR; j++) {
        int nb = vertexSets[j];
        for (int k = 0; k < numNeighbors[nb]; k++) {
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
        int cnt = 0;
        for (int k = 0; k < numNeighbors[tv]; k++) {
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
            auto [nb, ne] = G.getNbr(v);

            // P(v) = {u > v : u ∈ N(v)}
            int p = 0;
            for (int idx = nb; idx < ne; idx++)
                if (G.adj_list[idx] > v) p++;

            if (p == 0) {
                if (1 >= min_k) {
                    out.off.push_back((int)out.buf.size());
                    out.buf.push_back(1);
                    out.buf.push_back(v);
                    out.pbuf.push_back(0); // non-pivot singleton
                }
                continue;
            }

            int slabTop = slab.save();

            // l2g[i] = global ID of i-th P(v) member (sorted by global ID)
            int* l2g = slab.alloc(p + 1);
            int li = 0;
            for (int idx = nb; idx < ne; idx++) {
                int u = G.adj_list[idx];
                if (u > v) { gmark[u] = li; l2g[li] = u; li++; }
            }
            l2g[p] = v;

            // vertexSets[i] = i, vertexLookup[i] = i
            int* lsets = slab.alloc(p);
            int* llook = slab.alloc(p);
            for (int i = 0; i < p; i++) { lsets[i] = i; llook[i] = i; }

            // Build neighborsInP: iterate j=0..p-1 in order.
            // For each local j, scan its later neighbors w (l2g[w] > l2g[j]).
            // numNeighbors[u] = TOTAL P-degree.
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

            std::vector<int*> lnbr_ptrs((size_t)p);
            for (int i = 0; i < p; i++)
                lnbr_ptrs[i] = slab.alloc(lnbCnt[i] > 0 ? lnbCnt[i] : 1);

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
                     max_k, min_k,
                     l2g, slab, out);

            for (int idx = nb; idx < ne; idx++) gmark[G.adj_list[idx]] = -1;
            slab.restore(slabTop);
        }
    }

    for(int t = 0; t < nT; t++) std::free(slab_mem[t]);
    double t1 = omp_get_wtime();

    // --- CSR merge with pivot flags ---
    std::vector<size_t> coff((size_t)(nT+1), 0), voff((size_t)(nT+1), 0);
    for(int t = 0; t < nT; t++){
        coff[(size_t)(t+1)] = coff[(size_t)t] + thread_out[(size_t)t].size();
        voff[(size_t)(t+1)] = voff[(size_t)t] + thread_out[(size_t)t].total_vertices();
    }
    size_t TC = coff[(size_t)nT], TV = voff[(size_t)nT];

    std::vector<int> csrOff((size_t)(TC+1));
    std::vector<int> csrDat((size_t)TV);
    std::vector<uint8_t> csrPiv((size_t)TV);

    // Pass 1: fill per-clique sizes
    #pragma omp parallel for schedule(static) num_threads(nT)
    for(int t = 0; t < nT; t++){
        size_t base = coff[(size_t)t];
        const auto& o = thread_out[(size_t)t];
        for(size_t oi = 0; oi < o.size(); oi++)
            csrOff[base+oi] = o.buf[o.off[oi]];
    }
    // Prefix sum
    size_t run = 0;
    for(size_t i = 0; i < TC; i++){ size_t s=(size_t)csrOff[i]; csrOff[i]=(int)run; run+=s; }
    csrOff[TC] = (int)run;

    // Pass 2: fill vertex data and pivot flags
    #pragma omp parallel for schedule(static) num_threads(nT)
    for(int t = 0; t < nT; t++){
        size_t base = coff[(size_t)t];
        const auto& o = thread_out[(size_t)t];
        // Track pivot buffer offset: pbuf has entries only for vertex data
        // (not for the size header). We need to compute the pivot offset
        // for each clique in this thread's output.
        size_t pivOff = 0;
        for(size_t oi = 0; oi < o.size(); oi++){
            size_t dst = (size_t)csrOff[base+oi];
            int pos = o.off[oi]; int total_sz = o.buf[pos++];
            for(int i = 0; i < total_sz; i++) {
                csrDat[dst+i] = o.buf[pos+i];
                csrPiv[dst+i] = o.pbuf[pivOff+i];
            }
            pos += total_sz;
            pivOff += (size_t)total_sz;
        }
    }

    CliqueCSR<int> res;
    res.init_from_flat(std::move(csrOff), std::move(csrDat), std::move(csrPiv));
    printf("SDCT_Par7: par=%.1fms merge=%.1fms total=%zu\n",
           (t1-t0)*1000.0, (omp_get_wtime()-t1)*1000.0, TC);
    return res;
}