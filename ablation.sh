#!/bin/bash
set -e

# Ablation experiment script — builds 4 variants and benchmarks each
# Variants:
#   V0: Baseline (HEAD, heap-based peeling)
#   V1: +Bucket peeling only
#   V2: +Bucket +Batch Phase C (decouple tree mutations from support updates)
#   V3: +Bucket +Batch Phase C +Adaptive buildParallel (full)

PROJECT_DIR="$HOME/nclique_tmp"
SRC="$PROJECT_DIR/src/NucleusDecomposition/NucleusCoreDecompositionRemoveSclique.cpp"
HDR="$PROJECT_DIR/src/dataStruct/CliqueHashMap.h"
BUILD_DIR="$PROJECT_DIR/build"
BIN="$BUILD_DIR/bin/degeneracy_cliques"

GRAPHS=("/data/wenqianz/com-dblp.edges" "/data/wenqianz/web-Google.edges" "/data/wenqianz/com-youtube.edges")
GRAPH_NAMES=("com-dblp" "web-Google" "com-youtube")
PARAMS=("3 4" "4 5")

export OMP_NUM_THREADS=16

# Save current (V3 full) source
cp "$SRC" "$SRC.v3"
cp "$HDR" "$HDR.v3"

echo "============================================================"
echo "ABLATION EXPERIMENT"
echo "Server: $(hostname), Cores: $(nproc), Threads: $OMP_NUM_THREADS"
echo "Date: $(date)"
echo "============================================================"
echo ""

run_benchmark() {
    local variant="$1"
    echo "=== Variant: $variant ==="

    # Build
    cd "$BUILD_DIR"
    make degeneracy_cliques -j16 2>&1 | tail -1
    cd "$PROJECT_DIR"

    for gi in "${!GRAPHS[@]}"; do
        local graph="${GRAPHS[$gi]}"
        local gname="${GRAPH_NAMES[$gi]}"
        for param in "${PARAMS[@]}"; do
            local r=$(echo $param | cut -d' ' -f1)
            local s=$(echo $param | cut -d' ' -f2)
            echo "--- $variant | $gname r=$r s=$s ---"
            PIVOTER_COMPARE=1 OMP_NUM_THREADS=16 timeout 300 "$BIN" "$graph" $r $s 2>&1 | \
                grep -E "(clique Index build|countingPerRClique|time:|Time Breakdown|Init:|Pop:|Structure:|Support:|Heap:|Comparison|Reference|Optimized)"
            echo ""
        done
    done
}

########################################
# V0: Baseline (HEAD — heap-based)
########################################
echo ">>> Building V0: Baseline (heap) ..."
cd "$PROJECT_DIR"
git checkout -- "$SRC" "$HDR" 2>/dev/null || true
run_benchmark "V0_Baseline"

########################################
# V1: +Bucket peeling only
# Replace heap init/pop/update with bucket, but keep Phase C interleaved
########################################
echo ">>> Building V1: +Bucket ..."
# Restore V3 full source, then revert Phase C and buildParallel changes
cp "$SRC.v3" "$SRC"
cp "$HDR.v3" "$HDR"

# Revert adaptive buildParallel back to serial build
python3 -c "
import re
with open('$SRC', 'r') as f:
    src = f.read()
# Replace adaptive buildParallel block with simple serial build
old = '''    daf::timeCount(\"clique Index build\", [&]() {
#ifdef _OPENMP
        // Adaptive: use parallel build only when tree is large enough to amortize mutex overhead
        if (tree.adj_list.size() > 100000) {
            cliqueIndex.buildParallel(tree, edgeGraph.adj_list.size());
        } else {
            cliqueIndex.build(tree, edgeGraph.adj_list.size());
        }
#else
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
#endif
    });'''
new = '''    daf::timeCount(\"clique Index build\", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });'''
src = src.replace(old, new)
with open('$SRC', 'w') as f:
    f.write(src)
"

# Revert batch Phase C back to interleaved Phase C
python3 -c "
with open('$SRC', 'r') as f:
    src = f.read()

old_phase_c = '''        // ============ Phase C1: Serial tree mutations only ============
        {
            auto t_struct_B = std::chrono::high_resolution_clock::now();
            for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
                auto leafId = changedLeaf[idx];
                const auto leaf = tree.adj_list[leafId];
                LeafResult &res = leafResults[idx];

                for (auto leafV : leaf) {
                    if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                    else treeGraphV.removeNbr(leafV.v, {leafId, false});
                }
                for (auto &newLeaf : res.newLeaves) {
                    auto newId = tree.addNode(newLeaf);
                    const auto &stored = tree.adj_list[newId];
                    for (auto i : stored) {
                        if (i.isPivot) treeGraphV.addNbr(i.v, {newId, true});
                        else treeGraphV.addNbr(i.v, {newId, false});
                    }
                    if (newId >= changedLeafIndex.size())
                        changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                }
                tree.removeNode(leafId);
            }
            duration_structure += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_struct_B).count();
        }

        // ============ Phase C2: Batch support update + single bucket move ============
        {
            auto t_supp = std::chrono::high_resolution_clock::now();
            // Apply all deltas, collect unique affected IDs via flag
            std::vector<daf::Size> affectedIds;
            for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
                LeafResult &res = leafResults[idx];
                for (const auto &p : res.incr) {
                    countingRClique[p.first] += p.second;
                }
                for (const auto &p : res.decr) {
                    if (!rCliqueInHeap[p.first]) continue;
                    countingRClique[p.first] -= p.second;
                    int newB = std::max(0, (int)countingRClique[p.first]);
                    int oldB = bucket_of[p.first];
                    if (newB != oldB) {
                        auto& oldVec = buckets[oldB];
                        daf::Size myPos = pos_in_bucket[p.first];
                        if (myPos < oldVec.size() - 1) {
                            daf::Size last = oldVec.back();
                            oldVec[myPos] = last;
                            pos_in_bucket[last] = myPos;
                        }
                        oldVec.pop_back();
                        if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
                        bucket_of[p.first] = newB;
                        pos_in_bucket[p.first] = buckets[newB].size();
                        buckets[newB].push_back(p.first);
                        if (newB < curBucket) curBucket = newB;
                    }
                }
            }
            duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_supp).count();
        }'''

new_phase_c = '''        // ============ Phase C: Per-leaf serial mutations + support + bucket ============
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            auto leafId = changedLeaf[idx];
            const auto leaf = tree.adj_list[leafId];
            LeafResult &res = leafResults[idx];

            auto t_struct_B = std::chrono::high_resolution_clock::now();
            for (auto leafV : leaf) {
                if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                else treeGraphV.removeNbr(leafV.v, {leafId, false});
            }
            for (auto &newLeaf : res.newLeaves) {
                auto newId = tree.addNode(newLeaf);
                const auto &stored = tree.adj_list[newId];
                for (auto i : stored) {
                    if (i.isPivot) treeGraphV.addNbr(i.v, {newId, true});
                    else treeGraphV.addNbr(i.v, {newId, false});
                }
                if (newId >= changedLeafIndex.size())
                    changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
            }
            duration_structure += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_struct_B).count();

            auto t_supp = std::chrono::high_resolution_clock::now();
            for (const auto &p : res.incr) countingRClique[p.first] += p.second;
            duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_supp).count();

            auto t_heap = std::chrono::high_resolution_clock::now();
            for (const auto &p : res.decr) {
                countingRClique[p.first] -= p.second;
                int newB = std::max(0, (int)countingRClique[p.first]);
                int oldB = bucket_of[p.first];
                if (newB != oldB) {
                    auto& oldVec = buckets[oldB];
                    daf::Size myPos = pos_in_bucket[p.first];
                    if (myPos < oldVec.size() - 1) {
                        daf::Size last = oldVec.back();
                        oldVec[myPos] = last;
                        pos_in_bucket[last] = myPos;
                    }
                    oldVec.pop_back();
                    if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
                    bucket_of[p.first] = newB;
                    pos_in_bucket[p.first] = buckets[newB].size();
                    buckets[newB].push_back(p.first);
                    if (newB < curBucket) curBucket = newB;
                }
            }
            duration_heap += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_heap).count();

            auto t_struct_C = std::chrono::high_resolution_clock::now();
            tree.removeNode(leafId);
            duration_structure += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_struct_C).count();
        }'''

src = src.replace(old_phase_c, new_phase_c)
with open('$SRC', 'w') as f:
    f.write(src)
"
run_benchmark "V1_Bucket"

########################################
# V2: +Bucket +Batch Phase C
########################################
echo ">>> Building V2: +Bucket +BatchC ..."
cp "$SRC.v3" "$SRC"
cp "$HDR.v3" "$HDR"
# Revert adaptive buildParallel back to serial build
python3 -c "
with open('$SRC', 'r') as f:
    src = f.read()
old = '''    daf::timeCount(\"clique Index build\", [&]() {
#ifdef _OPENMP
        // Adaptive: use parallel build only when tree is large enough to amortize mutex overhead
        if (tree.adj_list.size() > 100000) {
            cliqueIndex.buildParallel(tree, edgeGraph.adj_list.size());
        } else {
            cliqueIndex.build(tree, edgeGraph.adj_list.size());
        }
#else
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
#endif
    });'''
new = '''    daf::timeCount(\"clique Index build\", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });'''
src = src.replace(old, new)
with open('$SRC', 'w') as f:
    f.write(src)
"
run_benchmark "V2_Bucket_BatchC"

########################################
# V3: +Bucket +Batch Phase C +Adaptive buildParallel (full)
########################################
echo ">>> Building V3: Full (all optimizations) ..."
cp "$SRC.v3" "$SRC"
cp "$HDR.v3" "$HDR"
run_benchmark "V3_Full"

# Cleanup
rm -f "$SRC.v3" "$HDR.v3"

echo "============================================================"
echo "ABLATION EXPERIMENT COMPLETE"
echo "============================================================"
