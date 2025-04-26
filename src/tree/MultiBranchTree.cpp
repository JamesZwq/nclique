//
// Created by 张文谦 on 25-3-3.
//

#include "MultiBranchTree.h"
#include "misc.h"

#include <functional>
#include <graph/Graph.h>
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

daf::StaticVector<double> MultiBranchTree::cliqueCount() const {
    daf::StaticVector<double> counts(root->MaxDeep + 1);
    counts.c_size = root->MaxDeep + 1;
    std::function<void(TreeNode *, daf::CliqueSize, daf::CliqueSize)> helper =
            [&](TreeNode *node, daf::CliqueSize pivotCount, daf::CliqueSize nonPivotCount) {
        if (node->children.empty()) {
            const daf::CliqueSize rsize = pivotCount + nonPivotCount;
            for (daf::CliqueSize i = 0; i <= pivotCount; i++) {
                const daf::Size k = rsize - i;
                counts[k] += nCr[pivotCount][i];
            }
            return;
        }
        for (const auto child: node->children) {
            if (child->isPivot) helper(child, pivotCount + 1, nonPivotCount);
            else helper(child, pivotCount, nonPivotCount + 1);
        }
    };

    helper(root, 0, 0);
    // std::cout << counts << std::endl;
    return counts;
}

double MultiBranchTree::unionCliqueCount(const Graph &leafGraph, const daf::CliqueSize k) const {
    double count = 0;

    std::function<void(TreeNode *, daf::CliqueSize, daf::CliqueSize)> helper =
            [&](const TreeNode *node, const daf::CliqueSize pivotCount, const daf::CliqueSize nonPivotCount) {
        if (nonPivotCount > k) {
            return;
        }
        if (node->children.empty()) {
            if (pivotCount + nonPivotCount < k) {
                return;
            }
            count += nCr[pivotCount][k - nonPivotCount];
            return;
        }
        for (const auto child: node->children) {
            if (child->isPivot) {
                if (leafGraph.getNbrCount(child->v) == 1) {
                    helper(child, pivotCount, nonPivotCount);
                } else {
                    helper(child, pivotCount + 1, nonPivotCount);
                }
            } else {
                if (leafGraph.getNbrCount(child->v) != 1) {
                    helper(child, pivotCount, nonPivotCount + 1);
                }
            }
        }
    };

    helper(root, 0, 0);
    std::cout << count << std::endl;
    return count;
}


void MultiBranchTree::unionCliqueList(const Graph &leafGraph, const daf::CliqueSize k) const {
    double count = 0;
    daf::StaticVector<daf::Size> povit;
    daf::StaticVector<daf::Size> keep;
    daf::Size numList = 0;
    std::function<void(const TreeNode *, daf::StaticVector<daf::Size>, daf::StaticVector<daf::Size>)> helper =
            [&](const TreeNode *node, daf::StaticVector<daf::Size> povitC, daf::StaticVector<daf::Size> keepC) {
        if (keepC.size() > k) {
            return;
        }
        if (node->children.empty()) {
            if (keepC.size() + povitC.size() < k) {
                return;
            }
            count += nCr[povitC.size()][k - keepC.size()];
            numList++;
            // std::cout << "povitC: " << povitC << " keepC: " << keepC << std::endl;
            // daf::generateCombinations<daf::Size>(keepC, povitC, k);
            return;
        }
        for (const auto child: node->children) {
            if (child->isPivot) {
                if (leafGraph.getNbrCount(child->v) == 1) {
                    helper(child, povitC, keepC);
                } else {
                    povitC.push_back(child->v);
                    helper(child, povitC, keepC);
                    povitC.pop_back();
                }
            } else {
                if (leafGraph.getNbrCount(child->v) != 1) {
                    keepC.push_back(child->v);
                    helper(child, povitC, keepC);
                    keepC.pop_back();
                }
            }
        }
    };
    helper(root, povit, keep);
    std::cout << "numList: " << numList << std::endl;
    std::cout << count << std::endl;
    povit.free();
    keep.free();
}


void MultiBranchTree::CliqueList(const daf::CliqueSize k) const {
    double count = 0;
    daf::StaticVector<daf::Size> povit;
    daf::StaticVector<daf::Size> keep;
    std::function<void(const TreeNode *, daf::StaticVector<daf::Size>, daf::StaticVector<daf::Size>)> helper =
            [&](const TreeNode *node, daf::StaticVector<daf::Size> povitC, daf::StaticVector<daf::Size> keepC) {
        if (keepC.size() > k) {
            return;
        }
        if (node->children.empty()) {
            if (keepC.size() + povitC.size() < k) {
                return;
            }
            count += nCr[povitC.size()][k - keepC.size()];
            std::cout << "povitC: " << povitC << " keepC: " << keepC << std::endl;
            return;
        }
        for (const auto child: node->children) {
            if (child->isPivot) {
                povitC.push_back(child->v);
                helper(child, povitC, keepC);
                povitC.pop_back();
            } else {
                keepC.push_back(child->v);
                helper(child, povitC, keepC);
                keepC.pop_back();
            }
        }
    };

    helper(root, povit, keep);
    std::cout << count << std::endl;
}

// daf::CliqueSize MultiBranchTree::initMaxDeepHelper(TreeNode *node, daf::CliqueSize deep) {
//     node->MaxDeep = deep;
//     for (const auto &child : node->children) {
//         node->MaxDeep = std::max(node->MaxDeep,initMaxDeepHelper(child, deep + 1));
//     }
//     return node->MaxDeep;
// }


void MultiBranchTree::initMaxDeep() const {
    std::function<daf::CliqueSize(TreeNode *, daf::CliqueSize)> dfs =
            [&](TreeNode *node, daf::CliqueSize deep) -> daf::CliqueSize {
        node->MaxDeep = deep;
        for (const auto &child: node->children) {
            node->MaxDeep = std::max(node->MaxDeep, dfs(child, deep + 1));
        }
        return node->MaxDeep;
    };
    dfs(root, 0);

    TreeNode::maxCliques = root->MaxDeep;
}

// static void printCliqueCountByLeafsListRow(TreeNode **leafs, const daf::Size maxId) {
//     std::vector<double> counts(TreeNode::maxCliques + 1, 0);
//     for (daf::Size i = 0; i < maxId; i++) {
//         auto node = leafs[i];
//         daf::Size pivotCount = 0;
//         daf::Size nonPivotCount = 0;
//         while (node != nullptr) {
//             if (node->isPivot) pivotCount++;
//             else nonPivotCount++;
//             node = node->parent;
//             if (node->v == ROOTID) {
//                 for (daf::CliqueSize j = 0; j <= pivotCount; j++) {
//                     daf::Size k = pivotCount + nonPivotCount - j;
//                     counts[k] += nCr[pivotCount][j];
//                 }
//                 break;
//             }
//         }
//     }
//
//     std::cout << counts << std::endl;
// }


[[nodiscard]] daf::StaticVector<TreeNode *> MultiBranchTree::getLeafsList(const daf::Size maxId,
                                                                          const daf::CliqueSize minK) const {
    // auto *leafs = new TreeNode*[maxId];
    daf::StaticVector<TreeNode *> leafs(maxId);
    leafs.c_size = maxId;
    std::function<void(TreeNode *)> rec = [&](TreeNode *node) {
        if (node->children.empty()) leafs[node->leafId] = node;
        for (const auto &child: node->children) {
            if (child->MaxDeep < minK) {
                continue;
            }
            rec(child);
        }
    };

    rec(root);

    return leafs;
}

[[nodiscard]] daf::Size MultiBranchTree::initLeafsParentAndId(const daf::CliqueSize minK) const {
    std::function<daf::Size(TreeNode *, daf::Size)> rec =
            [&](TreeNode *node, daf::Size maxId) -> daf::Size {
        if (node->children.empty()) {
            node->leafId = maxId;
            return maxId + 1;
        }
        for (const auto &child: node->children) {
            if (child->MaxDeep < minK) {
                continue;
            }
            maxId = rec(child, maxId);
            child->parent = node;
        }
        return maxId;
    };

    return rec(root, 0);
}

[[nodiscard]] double MultiBranchTree::numLeafs() const {
    // 定义一个 lambda 用于递归计算叶子数
    std::function<double(TreeNode *)> rec = [&](TreeNode *node) -> double {
        if (node->children.empty()) {
            return 1;
        }
        double sum = 0;
        for (const auto &child: node->children) {
            sum += rec(child);
        }
        return sum;
    };
    return rec(root);
}

void TreeNode::prettyPrint(std::ostream &os, int indent) const {
    // 打印缩进
    for (int i = 0; i < indent; i++) {
        os << "  "; // 两个空格缩进
    }
    // 如果是 pivot 节点，在前面添加标记
    if (isPivot) {
        os << "[Pivot] ";
    }
    os << v << " deep-" << MaxDeep;
    if (leafId != std::numeric_limits<daf::Size>::max()) {
        os << " leafId-" << leafId;
    }
    os << std::endl;
    // 递归打印每个子节点，缩进加 1
    // auto cCount = 0;
    for (const auto &child: children) {
        // std::cout << " v: " << v << " child: " << ++cCount << std::endl;
        child->prettyPrint(os, indent + 1);
    }
}

// void removeNodes(const TreeNode *leaf, daf::StaticVector<TreeNode *> &nodes) {
//     daf::Size index = 0;
//     auto currentNode = leaf;
//     auto previousNode = nullptr;
//     while (currentNode != nullptr) {
//         if (nodes[index] == currentNode) {
//             ++index;
//         }
//         currentNode = currentNode->parent;
//     }
// }


// Efficient in-place removal from leaf to its first multi-child ancestor
// Assumes 'nodes' has reserved capacity >= depth of tree
TreeNode *MultiBranchTree::removeNodes(TreeNode *leaf, daf::StaticVector<daf::Size> &nodes) {
    // Efficiently remove a chain of nodes (nodes) from leaf→root path,
    // reparenting their child to the next ancestor.
    TreeNode *curr = leaf;
    size_t idx = 0, n = nodes.size();
    bool leafRemoved = (n > 0 && nodes[0] == leaf->v);
    TreeNode *lastAncestor = nullptr;

    while (curr && idx < n) {
        if (curr->v == nodes[idx]) {
            TreeNode *parent = curr->parent;
            // Reparent all children of curr to its parent
            if (parent) {
                auto &parentKids = parent->children;
                for (TreeNode *child: curr->children) {
                    child->parent = parent;
                    parentKids.push_back(child);
                }
            }
            // Remove curr from its parent's children
            if (parent) {
                auto &siblings = parent->children;
                for (size_t i = 0, sz = siblings.size(); i < sz; ++i) {
                    if (siblings[i] == curr) {
                        siblings[i] = siblings.back();
                        siblings.pop_back();
                        break;
                    }
                }
            }
            // Clear curr's children and sever link
            curr->children.clear();
            curr->parent = nullptr;

            lastAncestor = parent;
            curr = parent;
            ++idx;
        } else {
            curr = curr->parent;
        }
    }
    // If the original leaf was removed, return its first surviving ancestor
    if (leafRemoved) {
        // lastAncestor->leafId = leaf->leafId;
        if (lastAncestor) {
            lastAncestor->leafId = leaf->leafId;
        }
        return lastAncestor ? lastAncestor : getRoot();
    }
    // Otherwise return original leaf
    return leaf;
}


// remove all node that only belong to the current leaf;
void MultiBranchTree::removeLeaf(TreeNode *leaf) {
    TreeNode *curr = leaf;
    while (curr->parent
           && !curr->parent->isRoot()
           && curr->parent->children.size() == 1) {
        curr = curr->parent;
    }

    for (size_t i = 0; i < curr->parent->children.size(); ++i) {
        if (curr->parent->children[i] == curr) {
            curr->parent->children[i] = curr->parent->children.back();
            curr->parent->children.pop_back();
            break;
        }
    }
    // delete curr;
}
