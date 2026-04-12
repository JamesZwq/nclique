#define R2_LAB_NAMESPACE DCLPHardLeafLab
#define R2_LAB_ENTRYPOINT PlusNucleusEdgeCoreDecompositionSet_HardLeafLab_AutoInternal
#define R2_LAB_DCLP_LABEL "HardLeaf Lab"
#define R2_LAB_HYBRID_LABEL "HardLeaf Hybrid Auto"
#define HARD_LEAF_LAB_PROFILE 1

#include "NCliqueEdgeCoreDecompositionHybridLab.cpp"

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_HardLeafLab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    return DCLPHardLeafLab::PlusNucleusEdgeCoreDecompositionSet_DCLP(tree, edgeGraph, treeGraphV, k);
}
