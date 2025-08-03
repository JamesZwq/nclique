//
// Created by 张文谦 on 25-3-25.
//

#ifndef NUCLEUSCOREDECOMPOSITION_H
#define NUCLEUSCOREDECOMPOSITION_H


#include <Global/Global.h>

#include "../tree/MultiBranchTree.h"

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > NucleusCoreDecomposition(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > NucleusCoreDecompositionHierarchy(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

std::vector<std::pair<std::vector<daf::Size>, int> > NucleusCoreDecompositionRClique(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);


std::vector<std::pair<std::vector<daf::Size>, int> > NucleusCoreDecompositionRCliqueOpt(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);
#endif //NUCLEUSCOREDECOMPOSITION_H
