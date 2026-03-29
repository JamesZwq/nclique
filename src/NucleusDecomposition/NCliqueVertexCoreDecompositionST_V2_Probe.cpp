//
// NCliqueVertexCoreDecompositionST_V2_Probe.cpp
//
// Feasibility probe for interleaved construction-decomposition.
// Measures what fraction of vertices could be peeled during SDCT
// construction using the finalization property: after vertex v's
// subtree completes in the degeneracy-ordered BK, countingV[v] is final.
//
// This does NOT actually peel — it just simulates to measure the potential.
//

#include "NCliqueCoreDecomposition.h"
#include "SDCT_Augmented.h"
#include <chrono>
#include <vector>
#include <algorithm>
#include <cstring>

extern double nCr[1001][401];

void NCliqueVertexCoreDecomposition_ST_V2_InterleavedProbe(
    Graph &edgeGraph, daf::CliqueSize k)
{
    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numVertices = edgeGraph.getGraphNodeSize();

    // Track per-vertex support accumulation
    auto *countingV = new double[numVertices];
    memset(countingV, 0, numVertices * sizeof(double));

    // Track which vertices are finalized (all their leaves have been emitted)
    std::vector<uint8_t> finalized(numVertices, 0);
    size_t numFinalized = 0;
    size_t numPeelableAtFinalization = 0;
    double runningMinCore = 0;

    // Track per-vertex: the vertex that "owns" this top-level BK call
    // so we know when a vertex's support is complete
    int currentTopVertex = -1;

    // Use a modified NoTree that has a per-vertex callback
    // We need the top-level loop to call onVertexDone after each vertex's subtree

    // Phase 1: Run SDCT and collect support + finalization order

    // We can't easily hook into SDCT_Augmented_NoTree's top-level loop
    // without modifying it. Instead, let's simulate: run the full SDCT,
    // then replay the finalization order.

    // Since SDCT processes vertices 0, 1, 2, ... in degeneracy order,
    // after vertex v's subtree completes, v will never appear in future leaves.
    // So countingV[v] is final after the last leaf containing v is emitted.

    // Strategy: track last-seen leaf for each vertex, then determine
    // finalization order.

    struct VertexInfo {
        daf::Size lastLeafId = 0;  // last leaf containing this vertex
        double finalSupport = 0;
    };
    std::vector<VertexInfo> vinfo(numVertices);

    size_t totalLeaves = 0;

    size_t numLeaves = SDCT_Augmented_NoTree(edgeGraph, k, /*min_k=*/1,
        [&](daf::Size leafId, const daf::StaticVector<int> &keepV,
            const daf::StaticVector<int> &dropV)
        {
            int pivotC = (int)dropV.size();
            int keepC = (int)keepV.size();
            int needP = (int)k - keepC;

            double wKeep = (needP >= 0 && needP <= pivotC) ? nCr[pivotC][needP] : 0.0;
            double wPivot = (needP >= 1 && pivotC >= 1) ? nCr[pivotC - 1][needP - 1] : 0.0;

            for (int i = 0; i < keepC; ++i) {
                daf::Size v = keepV[i];
                countingV[v] += wKeep;
                vinfo[v].lastLeafId = leafId;
            }
            for (int i = 0; i < pivotC; ++i) {
                daf::Size v = dropV[i];
                countingV[v] += wPivot;
                vinfo[v].lastLeafId = leafId;
            }
            totalLeaves++;
        });

    // Now countingV[v] = final support for all v.
    // Store final supports.
    for (daf::Size v = 0; v < numVertices; ++v) {
        vinfo[v].finalSupport = countingV[v];
    }

    // Simulate interleaved peeling:
    // Process vertices in degeneracy order (0, 1, 2, ...).
    // After each vertex v is "finalized", check if countingV[v] <= runningMinCore.
    // If so, it could be peeled immediately.
    //
    // Key insight: in the real peeling, earlier-peeled vertices reduce
    // neighbors' support. Here we use FINAL support (no early peeling deltas),
    // which gives a LOWER BOUND on peelability.
    // The actual interleaved version would peel more because early peeling
    // would reduce other vertices' support, enabling more peeling.

    // Build finalization order: vertex v is finalized when its lastLeafId's
    // vertex-owner finishes. In degeneracy ordering, vertex 0's subtree is
    // first, so all vertices whose lastLeafId falls within vertex 0's subtree
    // range become finalizable after vertex 0.
    //
    // Simpler approach: vertex v's support is final after ALL leaves containing v
    // have been emitted. The last leaf containing v is vinfo[v].lastLeafId.
    // Since leaves are emitted sequentially, vertex v is finalized when
    // leafCounter > vinfo[v].lastLeafId.

    // Sort vertices by their lastLeafId to determine finalization order
    std::vector<daf::Size> finOrder(numVertices);
    std::iota(finOrder.begin(), finOrder.end(), 0);
    std::sort(finOrder.begin(), finOrder.end(), [&](daf::Size a, daf::Size b) {
        return vinfo[a].lastLeafId < vinfo[b].lastLeafId;
    });

    // Simulate a greedy peeling from the finalization order
    double minCoreSoFar = 0;
    size_t peelableEarly = 0;     // could be peeled during construction
    size_t peelableLater = 0;     // final support > 0 but finalized late
    size_t totalWithSupport = 0;

    for (daf::Size v = 0; v < numVertices; ++v) {
        if (vinfo[v].finalSupport > 0) totalWithSupport++;
    }

    // Track: for each vertex in finalization order, what fraction of
    // total leaves have been emitted when it finalizes?
    size_t earlyHalf = 0;   // finalized in first half of leaf emission
    size_t lateHalf = 0;
    size_t halfLeaves = numLeaves / 2;

    for (auto v : finOrder) {
        if (vinfo[v].finalSupport <= 0) continue;
        if (vinfo[v].lastLeafId < halfLeaves) earlyHalf++;
        else lateHalf++;
    }

    // Count vertices peelable at minimum core level
    // Sort by final support to simulate peeling
    std::vector<std::pair<double, daf::Size>> supportOrder;
    supportOrder.reserve(numVertices);
    for (daf::Size v = 0; v < numVertices; ++v) {
        if (vinfo[v].finalSupport > 0)
            supportOrder.push_back({vinfo[v].finalSupport, v});
    }
    std::sort(supportOrder.begin(), supportOrder.end());

    // Find how many vertices are peeled at each core level
    // and what fraction of them are finalized by the time they're peelable
    double curCore = 0;
    size_t peeled = 0;
    size_t peelableBeforeAllDone = 0;

    for (auto &[support, v] : supportOrder) {
        curCore = std::max(curCore, support);
        peeled++;
        // This vertex is peelable when curCore >= support.
        // It's finalized when lastLeafId < numLeaves.
        // It could be peeled during construction if lastLeafId < numLeaves * (peeled / totalWithSupport)
        // (rough heuristic: its finalization fraction < its peel fraction)
        double finFraction = (double)(vinfo[v].lastLeafId + 1) / numLeaves;
        double peelFraction = (double)peeled / totalWithSupport;
        if (finFraction <= peelFraction * 0.8) { // 80% margin
            peelableBeforeAllDone++;
        }
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::cout << "=== Interleaved Peeling Feasibility Probe ===" << std::endl;
    std::cout << "Vertices: " << numVertices << ", Leaves: " << numLeaves << std::endl;
    std::cout << "Vertices with support > 0: " << totalWithSupport << std::endl;
    std::cout << "Finalized in first half of leaf emission: " << earlyHalf
              << " (" << (100.0 * earlyHalf / totalWithSupport) << "%)" << std::endl;
    std::cout << "Finalized in second half: " << lateHalf
              << " (" << (100.0 * lateHalf / totalWithSupport) << "%)" << std::endl;
    std::cout << "Peelable before all leaves done (conservative): " << peelableBeforeAllDone
              << " (" << (100.0 * peelableBeforeAllDone / totalWithSupport) << "%)" << std::endl;
    std::cout << "Probe took: " << elapsed << " ms" << std::endl;

    delete[] countingV;
}
