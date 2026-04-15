#define R2_LAB_NAMESPACE DCLPDefectLab
#define R2_LAB_ENTRYPOINT PlusNucleusEdgeCoreDecompositionSet_DefectLab_AutoInternal
#define R2_LAB_DCLP_LABEL "Defect Lab"
#define R2_LAB_HYBRID_LABEL "Defect Lab"
#define R2_LAB_ENABLE_CASEB_DEFECT_GRAPH 1

#include "NCliqueEdgeCoreDecompositionHybridLab.cpp"

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_DefectLab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    return DCLPDefectLab::PlusNucleusEdgeCoreDecompositionSet_DCLP(tree, edgeGraph, treeGraphV, k);
}
