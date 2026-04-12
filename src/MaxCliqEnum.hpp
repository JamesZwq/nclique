#pragma once
//
// MaxCliqEnum: maximal clique enumerator using the quick-cliques algorithm
// (Eppstein, Löffler, Strash — SEA 2011).
// Degeneracy-ordered BK with P∪X pivoting, vertexSets/vertexLookup data structures.
// Returns vector<vector<daf::Size>> of maximal cliques with size ≥ minSize.
//

#include "Global/Global.h"
#include "misc.h"
#include <algorithm>
#include <vector>
#include <cstring>

namespace maxcliq {

// ─── Pivot selection: P∪X ───
// Finds pivot from P∪X maximizing |N(u)∩P|, returns P \ N(pivot).
inline int findPivotPX(int **pivotNonNeighbors, int *numNonNeighbors,
                       int *vertexSets, int *vertexLookup,
                       int **neighborsInP, int *numNeighbors,
                       int beginX, int beginP, int beginR)
{
    int pivot = -1, maxCnt = -1;
    // Search P∪X (from beginX to beginR)
    for (int j = beginX; j < beginR; ++j) {
        int v = vertexSets[j];
        int numPot = std::min(beginR - beginP, numNeighbors[v]);
        int cnt = 0;
        for (int k = 0; k < numPot; ++k) {
            int nb = neighborsInP[v][k];
            int loc = vertexLookup[nb];
            if (loc >= beginP && loc < beginR) cnt++;
            else break;
        }
        if (cnt > maxCnt) { maxCnt = cnt; pivot = v; }
    }

    int sizeOfP = beginR - beginP;
    *pivotNonNeighbors = (int *)Calloc(sizeOfP, sizeof(int));
    memcpy(*pivotNonNeighbors, &vertexSets[beginP], sizeOfP * sizeof(int));
    *numNonNeighbors = sizeOfP;

    int numPivotNbr = std::min(sizeOfP, numNeighbors[pivot]);
    for (int k = 0; k < numPivotNbr; ++k) {
        int nb = neighborsInP[pivot][k];
        int loc = vertexLookup[nb];
        if (loc >= beginP && loc < beginR)
            (*pivotNonNeighbors)[loc - beginP] = -1;
        else break;
    }

    // Compact
    int w = 0;
    for (int j = 0; j < *numNonNeighbors; ++j) {
        if ((*pivotNonNeighbors)[j] != -1)
            (*pivotNonNeighbors)[w++] = (*pivotNonNeighbors)[j];
    }
    *numNonNeighbors = w;
    return pivot;
}

// ─── moveToR with X support (ported from quick-cliques) ───
inline void moveToR(int vertex,
                    int *vertexSets, int *vertexLookup,
                    int **neighborsInP, int *numNeighbors,
                    int *pBeginX, int *pBeginP, int *pBeginR,
                    int *pNewBeginX, int *pNewBeginP, int *pNewBeginR)
{
    int vertexLocation = vertexLookup[vertex];
    (*pBeginR)--;
    vertexSets[vertexLocation] = vertexSets[*pBeginR];
    vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
    vertexSets[*pBeginR] = vertex;
    vertexLookup[vertex] = *pBeginR;

    *pNewBeginX = *pBeginP;
    *pNewBeginP = *pBeginP;
    *pNewBeginR = *pBeginP;

    int sizeOfP = *pBeginR - *pBeginP;

    // Move X neighbors of vertex into new X
    int j = *pBeginX;
    while (j < *pNewBeginX) {
        int neighbor = vertexSets[j];
        int neighborLocation = j;
        int incrementJ = 1;
        int numPot = std::min(sizeOfP, numNeighbors[neighbor]);
        for (int k = 0; k < numPot; ++k) {
            if (neighborsInP[neighbor][k] == vertex) {
                (*pNewBeginX)--;
                vertexSets[neighborLocation] = vertexSets[(*pNewBeginX)];
                vertexLookup[vertexSets[(*pNewBeginX)]] = neighborLocation;
                vertexSets[(*pNewBeginX)] = neighbor;
                vertexLookup[neighbor] = (*pNewBeginX);
                incrementJ = 0;
            }
        }
        if (incrementJ) j++;
    }

    // Move P neighbors of vertex into new P
    j = *pBeginP;
    while (j < *pBeginR) {
        int neighbor = vertexSets[j];
        int neighborLocation = j;
        int numPot = std::min(sizeOfP, numNeighbors[neighbor]);
        for (int k = 0; k < numPot; ++k) {
            if (neighborsInP[neighbor][k] == vertex) {
                vertexSets[neighborLocation] = vertexSets[(*pNewBeginR)];
                vertexLookup[vertexSets[(*pNewBeginR)]] = neighborLocation;
                vertexSets[(*pNewBeginR)] = neighbor;
                vertexLookup[neighbor] = (*pNewBeginR);
                (*pNewBeginR)++;
            }
        }
        j++;
    }

    // Update neighborsInP for new X∪P
    j = *pNewBeginX;
    while (j < *pNewBeginR) {
        int thisVertex = vertexSets[j];
        int numPot = std::min(sizeOfP, numNeighbors[thisVertex]);
        int numNeighborsInP = 0;
        for (int k = 0; k < numPot; ++k) {
            int nb = neighborsInP[thisVertex][k];
            int loc = vertexLookup[nb];
            if (loc >= *pNewBeginP && loc < *pNewBeginR) {
                neighborsInP[thisVertex][k] = neighborsInP[thisVertex][numNeighborsInP];
                neighborsInP[thisVertex][numNeighborsInP] = nb;
                numNeighborsInP++;
            }
        }
        j++;
    }
}

// ─── moveFromRToX ───
inline void moveFromRToX(int vertex,
                         int *vertexSets, int *vertexLookup,
                         int *pBeginX, int *pBeginP, int *pBeginR)
{
    int vertexLocation = vertexLookup[vertex];
    vertexSets[vertexLocation] = vertexSets[*pBeginP];
    vertexLookup[vertexSets[*pBeginP]] = vertexLocation;
    vertexSets[*pBeginP] = vertex;
    vertexLookup[vertex] = *pBeginP;
    (*pBeginP)++;
    (*pBeginR)++;
}

// ─── BK recursion with P∪X pivot ───
static void bkRecurse(
    int *vertexSets, int *vertexLookup,
    int **neighborsInP, int *numNeighbors,
    int beginX, int beginP, int beginR,
    std::vector<int> &clique,
    std::vector<std::vector<daf::Size>> &output,
    int minSize)
{
    int sizeOfP = beginR - beginP;

    // Look-ahead pruning
    if (minSize > 0 && (int)clique.size() + sizeOfP < minSize) return;

    // P empty: maximal iff X also empty
    if (beginX >= beginP && beginP >= beginR) {
        if ((int)clique.size() >= minSize) {
            std::vector<daf::Size> c(clique.begin(), clique.end());
            std::sort(c.begin(), c.end());
            output.push_back(std::move(c));
        }
        return;
    }
    if (beginP >= beginR) return; // P empty, X non-empty → not maximal

    int *candidates;
    int numCandidates;
    findPivotPX(&candidates, &numCandidates,
                vertexSets, vertexLookup,
                neighborsInP, numNeighbors,
                beginX, beginP, beginR);

    if (numCandidates != 0) {
        for (int i = 0; i < numCandidates; ++i) {
            int vertex = candidates[i];
            clique.push_back(vertex);

            int newBeginX, newBeginP, newBeginR;
            moveToR(vertex, vertexSets, vertexLookup,
                    neighborsInP, numNeighbors,
                    &beginX, &beginP, &beginR,
                    &newBeginX, &newBeginP, &newBeginR);

            bkRecurse(vertexSets, vertexLookup,
                      neighborsInP, numNeighbors,
                      newBeginX, newBeginP, newBeginR,
                      clique, output, minSize);

            clique.pop_back();

            moveFromRToX(vertex, vertexSets, vertexLookup,
                         &beginX, &beginP, &beginR);
        }

        // Restore: move candidates from X back to P
        for (int i = 0; i < numCandidates; ++i) {
            int vertex = candidates[i];
            int vertexLocation = vertexLookup[vertex];
            beginP--;
            vertexSets[vertexLocation] = vertexSets[beginP];
            vertexLookup[vertexSets[beginP]] = vertexLocation;
            vertexSets[beginP] = vertex;
            vertexLookup[vertex] = beginP;
        }
    }

    Free(candidates);
}

} // namespace maxcliq

// Enumerate maximal cliques in graph g (degeneracy-ordered)
// Uses efficient vertexSets/vertexLookup data structures with P∪X pivoting
// minSize: only return cliques with size >= minSize (0 = all)
static std::vector<std::vector<daf::Size>> enumerateMaximalCliques(Graph &g, int minSize = 0) {
    int n = g.getGraphNodeSize();
    std::vector<std::vector<daf::Size>> result;

    int *vertexSets = (int *)Calloc(n, sizeof(int));
    int *vertexLookup = (int *)Calloc(n, sizeof(int));
    int **neighborsInP = (int **)Calloc(n, sizeof(int *));
    int *numNeighbors = (int *)Calloc(n, sizeof(int));

    for (int i = 0; i < n; ++i) {
        vertexLookup[i] = i;
        vertexSets[i] = i;
        neighborsInP[i] = (int *)Calloc(1, sizeof(int));
        numNeighbors[i] = 1;
    }

    int beginX = 0, beginP = 0, beginR = n;

    for (int v = 0; v < n; ++v) {
        int newBeginX, newBeginP, newBeginR;

        // This sets up P (later neighbors of v) but NOT X
        fillInPandXForRecursiveCallDegeneracyCliquesEdgeGraph(v,
            vertexSets, vertexLookup, g,
            neighborsInP, numNeighbors,
            &beginX, &beginP, &beginR,
            &newBeginX, &newBeginP, &newBeginR);

        // Set up X: swap earlier neighbors of v into X region before P
        // (fillInPandX only sets up P; we need X for maximality check)
        {
            auto [nb, ne] = g.getNbr(v);
            for (int idx = nb; idx < ne; ++idx) {
                int neighbor = g.adj_list[idx];
                if (neighbor >= v) continue; // only earlier neighbors
                int neighborLocation = vertexLookup[neighbor];

                newBeginX--;
                vertexSets[neighborLocation] = vertexSets[newBeginX];
                vertexLookup[vertexSets[newBeginX]] = neighborLocation;
                vertexSets[newBeginX] = neighbor;
                vertexLookup[neighbor] = newBeginX;

                // Rebuild neighborsInP for this X vertex
                Free(neighborsInP[neighbor]);
                neighborsInP[neighbor] = (int *)Calloc(
                    std::min(newBeginR - newBeginP, (int)g.getNbrCount(neighbor)),
                    sizeof(int));
                numNeighbors[neighbor] = 0;

                // Fill in neighbors of this X vertex that are in P
                auto [nb2, ne2] = g.getNbr(neighbor);
                for (int idx2 = nb2; idx2 < ne2; ++idx2) {
                    int laterNeighbor = g.adj_list[idx2];
                    int laterLoc = vertexLookup[laterNeighbor];
                    if (laterLoc >= newBeginP && laterLoc < newBeginR) {
                        neighborsInP[neighbor][numNeighbors[neighbor]] = laterNeighbor;
                        numNeighbors[neighbor]++;
                    }
                }
            }
        }

        std::vector<int> clique = {v};
        maxcliq::bkRecurse(vertexSets, vertexLookup,
                           neighborsInP, numNeighbors,
                           newBeginX, newBeginP, newBeginR,
                           clique, result, minSize);

        beginR = beginR + 1;
    }

    Free(vertexSets);
    Free(vertexLookup);
    for (int i = 0; i < n; ++i) Free(neighborsInP[i]);
    Free(neighborsInP);
    Free(numNeighbors);

    return result;
}
