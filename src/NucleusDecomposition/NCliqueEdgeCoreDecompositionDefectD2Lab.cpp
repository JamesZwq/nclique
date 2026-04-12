#define R2_LAB_NAMESPACE DCLPDefectD2Lab
#define R2_LAB_ENTRYPOINT PlusNucleusEdgeCoreDecompositionSet_DefectD2Lab_AutoInternal
#define R2_LAB_DCLP_LABEL "Defect D2 Lab"
#define R2_LAB_HYBRID_LABEL "Defect D2 Lab"
#define R2_LAB_ENABLE_CASEB_DEFECT_GRAPH 1
#define R2_LAB_ENABLE_CASEB_LEAF_REUSE 1
#define R2_LAB_ENABLE_CASEB_ZEROEDGE_CSR 1
#define R2_LAB_ENABLE_DEFECT_D2_IMMEDIATE 1

#include "NCliqueEdgeCoreDecompositionHybridLab.cpp"

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_DefectD2Lab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    int minLeafVerts = 20;
    int minPivots = 16;
    int minActiveEdges = 2;
    if (const char *env = std::getenv("PIVOTER_R2_DEFECT_D2_LAB_MIN_LEAF_VERTS"))
        minLeafVerts = std::max(0, std::atoi(env));
    if (const char *env = std::getenv("PIVOTER_R2_DEFECT_D2_LAB_MIN_PIVOTS"))
        minPivots = std::max(0, std::atoi(env));
    if (const char *env = std::getenv("PIVOTER_R2_DEFECT_D2_LAB_MIN_ACTIVE_EDGES"))
        minActiveEdges = std::max(0, std::atoi(env));

    const DCLPDefectD2Lab::DCLP::DCLPOptions options{
        "Defect D2 Lab",
        false,
        0,
        minLeafVerts,
        minPivots,
        minActiveEdges
    };
    return DCLPDefectD2Lab::PlusNucleusEdgeCoreDecompositionSet_DCLP_Impl(
        tree, edgeGraph, treeGraphV, k, options);
}
