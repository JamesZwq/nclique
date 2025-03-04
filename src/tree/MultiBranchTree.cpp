//
// Created by 张文谦 on 25-3-3.
//

#include "MultiBranchTree.h"
#include "misc.h"

//{
//    std::cout << "Tree:\n";
//    root->prettyPrint(std::cout, 0);
//}
extern double nCr[1001][401];
daf::CliqueSize TreeNode::maxCliques = 0;

void MultiBranchTree::printTree() const {
    std::cout << "Tree:\n";
    root->prettyPrint(std::cout, 0);
}

void MultiBranchTree::cliqueCountHelper(TreeNode* node, daf::CliqueSize pivotCount, daf::CliqueSize nonPivotCount, std::vector<daf::Size> &cliqueCounts) {
    if (node->children.empty()) {
        int rsize = pivotCount + nonPivotCount;
        for (int i = 0; i <= pivotCount; i++) {
            int k = rsize - i;
            cliqueCounts[k] += nCr[pivotCount][i];
        }
        return;
    }

    // 对于每个子节点，若子节点是 pivot 则 pivotCount+1，否则 nonPivotCount+1
    for (const auto child : node->children) {
        if (child->isPivot) {
            cliqueCountHelper(child, pivotCount + 1, nonPivotCount, cliqueCounts);
        } else {
            cliqueCountHelper(child, pivotCount, nonPivotCount + 1, cliqueCounts);
        }
    }
}

void MultiBranchTree::cliqueCount() {
    // 如果根节点没有子节点，则返回空计数（这通常不会发生）
    std::vector<daf::Size> counts(TreeNode::maxCliques + 1, 0);
    for (auto child : root->children) {
        if (child->isPivot) {
            cliqueCountHelper(child, 1, 0, counts);
        } else {
            cliqueCountHelper(child, 0, 1, counts);
        }
    }
    std::cout <<counts << std::endl;
}