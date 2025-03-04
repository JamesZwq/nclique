#ifndef TREE_H
#define TREE_H

#include <vector>
#include <iostream>
#include <limits>
#include <sstream>
#include "Global.h"

class TreeNode {
public:
    daf::Size v;
    bool isPivot;  // 注意将原来的 isPiovt 改为 isPivot
    std::vector<TreeNode *> children;
    daf::CliqueSize deep = 1;

    static daf::CliqueSize maxCliques;

    TreeNode(int data, bool isPivot) : v(data), isPivot(isPivot) {
        children.reserve(4);
    }

    TreeNode (const TreeNode&) = delete;
    TreeNode& operator=(const TreeNode&) = delete;

    TreeNode* addChild(int data, bool isPivot) {
        auto *child = new TreeNode(data, isPivot);
        child->deep = deep + 1;
        maxCliques = std::max(maxCliques, child->deep);
        children.push_back(child);
        return child;
    }

    ~TreeNode() {
        for (auto child : children) {
            delete child;
        }
    }

    // 递归打印节点，indent 表示缩进级别
    void prettyPrint(std::ostream &os, int indent = 0) const {
        // 打印缩进
        for (int i = 0; i < indent; i++) {
            os << "  "; // 两个空格缩进
        }
        // 如果是 pivot 节点，在前面添加标记
        if (isPivot) {
            os << "[Pivot] ";
        }
        os << v << "\n";
        // 递归打印每个子节点，缩进加 1
        auto cCount = 0;
        for (const auto& child : children) {
            // std::cout << " v: " << v << " child: " << ++cCount << std::endl;
            child->prettyPrint(os, indent + 1);
        }
    }
};

class MultiBranchTree {
public:

    TreeNode *root;

    MultiBranchTree() {
        root = new TreeNode(std::numeric_limits<daf::Size>::max(), false);
    }

    MultiBranchTree(const MultiBranchTree&) = delete;
    MultiBranchTree& operator=(const MultiBranchTree&) = delete;

    TreeNode* getRoot() const {
        return root;
    }

    // 重载输出运算符调用 prettyPrint
    friend std::ostream &operator<<(std::ostream &os, const MultiBranchTree &tree) {
        tree.root->prettyPrint(os, 0);
        return os;
    }

    void printTree() const;


    void cliqueCount();

private:
    void cliqueCountHelper(TreeNode *node, daf::CliqueSize pivotCount, daf::CliqueSize nonPivotCount,
                           std::vector<daf::Size> &cliqueCounts);
    // void cliqueCountHelper(TreeNode* node, int pivotCount, int nonPivotCount, std::vector<daf::CliqueSize> &cliqueCounts) const;
};

#endif //TREE_H