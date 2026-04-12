#define R2_LAB_NAMESPACE DCLPHardLeafHybridLab
#define R2_LAB_ENTRYPOINT PlusNucleusEdgeCoreDecompositionSet_HardLeafHybridLab_AutoInternal
#define R2_LAB_DCLP_LABEL "HardLeaf Hybrid Lab"
#define R2_LAB_HYBRID_LABEL "HardLeaf Hybrid Lab"
#define R2_LAB_ENABLE_CASEB_LEAF_REUSE 1
#define R2_LAB_ENABLE_CASEB_ZEROEDGE_CSR 1

#include "NCliqueEdgeCoreDecompositionHybridLab.cpp"

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_HardLeafHybridLab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    return DCLPHardLeafHybridLab::PlusNucleusEdgeCoreDecompositionSet_DCLP(tree, edgeGraph, treeGraphV, k);
}
