#ifndef TREE_H
#define TREE_H

#include <fstream>
#include <vector>
#include <iostream>
#include <limits>
#include <sstream>

#include <vector>
#include <iostream>
#include <limits>
#include <sstream>

#include <iostream>
#include "Global.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>

class TreeNode {
public:
    daf::Size v;
    bool isPivot; // 注意将原来的 isPiovt 改为 isPivot
    std::vector<TreeNode *> children;
    daf::CliqueSize MaxDeep = 1;
    // daf::Size leafId = std::numeric_limits<daf::Size>::max();
    static daf::CliqueSize maxCliques;
    TreeNode *parent;

    TreeNode() : v(0), isPivot(false) {
        children.reserve(4);
    }


    TreeNode(daf::Size data, bool isPivot) : v(data), isPivot(isPivot) {
        children.reserve(4);
    }

    TreeNode(const TreeNode &) = delete;

    TreeNode &operator=(const TreeNode &) = delete;

    TreeNode *addChild(daf::Size data, bool isPivot) {
        auto *child = new TreeNode(data, isPivot);
        // child->deep = deep + 1;
        // maxCliques = std::max(maxCliques, child->deep);
        children.push_back(child);
        return child;
    }

    virtual ~TreeNode() {
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
        os << v << " deep-" << MaxDeep << std::endl;
        // 递归打印每个子节点，缩进加 1
        // auto cCount = 0;
        for (const auto &child: children) {
            // std::cout << " v: " << v << " child: " << ++cCount << std::endl;
            child->prettyPrint(os, indent + 1);
        }
    }

private:
    // Boost.Serialization 需要访问私有成员
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int) {
        ar & v;
        ar & isPivot;
        ar & children;
        ar & MaxDeep;
        ar & parent;
    }

};

class MultiBranchTree {
public:
    TreeNode *root;

    MultiBranchTree() {
        root = new TreeNode(std::numeric_limits<daf::Size>::max(), false);
    }

    MultiBranchTree(const MultiBranchTree &) = delete;

    MultiBranchTree &operator=(const MultiBranchTree &) = delete;

    TreeNode *getRoot() const {
        return root;
    }

    // 重载输出运算符调用 prettyPrint
    friend std::ostream &operator<<(std::ostream &os, const MultiBranchTree &tree) {
        tree.root->prettyPrint(os, 0);
        return os;
    }

    void initMaxDeep() const;

    void printTree() const;

    void cliqueCount();

    void serialize(const std::string &filename) const {
        // auto filename = "/data/wenqianz/test.bin";
        std::ofstream ofs(filename, std::ios::binary);
        if (!ofs) {
            std::cerr << "Failed to serialize Tree, failed to open " << filename << std::endl;
            return;
        }
        std::cout << "Tree serialized to " << filename << std::endl;

        try {
            boost::archive::binary_oarchive oa(ofs);
            oa << *this;
        } catch(const boost::archive::archive_exception &e) {
            std::cerr << "Archive exception: " << e.what() << std::endl;
            return;
        }
        std::cout << "serialize success" << std::endl;
    }

    static MultiBranchTree *deserialize(const std::string &filename) {
        auto *newTree = new MultiBranchTree();
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs) {
            std::cerr << "Failed to deserialize Tree, failed to open " << filename << std::endl;
            return nullptr;
        }
        boost::archive::binary_iarchive ia(ifs);
        ia >> *newTree;
        return newTree;
    }

    /**
     * @brief 计算树中叶子节点的数量。
     *
     * 使用内部 lambda 表达式进行递归遍历，从根节点开始统计所有叶子节点（即没有子节点的节点）。
     *
     * @return 叶子节点总数。
     */
    [[nodiscard]] double numLeafs() const {
        // 定义一个 lambda 用于递归计算叶子数
        std::function<double(TreeNode*)> rec = [&](TreeNode* node) -> double {
            if (node->children.empty()) {
                return 1;
            }
            double sum = 0;
            for (const auto &child : node->children) {
                sum += rec(child);
            }
            return sum;
        };
        return rec(root);
    }

    /**
     * 初始化叶子节点的parent指针
     */
    void initLeafsParent() const {
        // 定义一个 lambda 用于递归初始化叶子节点的 ID
        std::function<void(TreeNode*)> rec = [&](TreeNode* node) -> void {
            if (node->children.empty()) return;
            for (const auto &child : node->children) {
                rec(child);
                child->parent = node;
            }
        };
        rec(root);
    }

private:

    void cliqueCountHelper(TreeNode *node, daf::CliqueSize pivotCount, daf::CliqueSize nonPivotCount,
                           std::vector<double> &cliqueCounts);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int) {
        ar & root;
    }
    // void cliqueCountHelper(TreeNode* node, int pivotCount, int nonPivotCount, std::vector<daf::CliqueSize> &cliqueCounts) const;
};

#endif //TREE_H
