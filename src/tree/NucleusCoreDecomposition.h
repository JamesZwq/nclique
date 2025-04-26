//
// Created by 张文谦 on 25-3-25.
//

#ifndef NUCLEUSCOREDECOMPOSITION_H
#define NUCLEUSCOREDECOMPOSITION_H


#include <Global/Global.h>

#include "MultiBranchTree.h"

void NucleusCoreDecomposition(MultiBranchTree &tree,
                              daf::StaticVector<TreeNode *> leafList,
                              Graph treeGraphV,
                              Graph leafGraph,
                              daf::CliqueSize r, daf::CliqueSize s);


#endif //NUCLEUSCOREDECOMPOSITION_H
