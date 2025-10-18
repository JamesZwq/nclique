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
#include "../Global/Global.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>
#include <graph/Graph.h>

static constexpr daf::Size ROOTID = std::numeric_limits<daf::Size>::max();
extern double nCr[1001][401];


class TreeNode final {
public:
    daf::Size v;
    bool isPivot; // 注意将原来的 isPiovt 改为 isPivot
    std::vector<TreeNode *> children;
    daf::CliqueSize MaxDeep = 1;
    daf::Size leafId = std::numeric_limits<daf::Size>::max();
    static daf::CliqueSize maxCliques;
    TreeNode *parent = nullptr;

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
        children.push_back(child);
        return child;
    }

    bool operator==(const TreeNode &other) const {
        return
                v == other.v &&
                isPivot == other.isPivot &&
                MaxDeep == other.MaxDeep &&
                leafId == other.leafId;
    }

    ~TreeNode() {
        // std::cout << "delete node: " << v << std::endl;
        for (const auto child: children) delete child;
    }

    // 递归打印节点，indent 表示缩进级别
    void prettyPrint(std::ostream &os, int indent = 0) const;

    void prettyPrint() const { prettyPrint(std::cout, 0); }

    friend std::ostream &operator<<(std::ostream &os, const TreeNode &node) {
        node.prettyPrint(os, 0);
        return os;
    }

    bool isRoot() const {
        return v == ROOTID;
    }

    // TreeNode &operator=(daf::StaticVector<TreeNode>::iterator node);


private:
    // Boost.Serialization 需要访问私有成员
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int) {
        ar & v;
        ar & isPivot;
        ar & children;
        ar & MaxDeep;
        ar & leafId;
        ar & parent;
    }
};

class MultiBranchTree {
public:
    TreeNode *root;

    MultiBranchTree() {
        root = nullptr;
    }


    MultiBranchTree(const MultiBranchTree &) = delete;

    MultiBranchTree &operator=(const MultiBranchTree &) = delete;

    MultiBranchTree(MultiBranchTree &&oldT) noexcept {
        if (this != &oldT) {
            root = oldT.root;
            oldT.root = nullptr;
        }
    };

    MultiBranchTree &operator=(MultiBranchTree &&oldT) = delete;

    [[nodiscard]] TreeNode *getRoot() {
        if (root == nullptr) {
            root = new TreeNode(ROOTID, false);
        }
        return root;
    }

    [[nodiscard]] TreeNode *getRoot() const {
        if (root == nullptr) {
            std::cerr << "root is null" << std::endl;
        }
        return root;
    }

    // 重载输出运算符调用 prettyPrint
    friend std::ostream &operator<<(std::ostream &os, const MultiBranchTree &tree) {
        tree.root->prettyPrint(os, 0);
        return os;
    }

    void initMaxDeep() const;

    void printTree() const;

    [[nodiscard]] daf::StaticVector<double> cliqueCount() const;

    [[nodiscard]] double unionCliqueCount(const Graph &leafGraph, daf::CliqueSize k) const;

    void unionCliqueList(const Graph &leafGraph, daf::CliqueSize k) const;

    void CliqueList(daf::CliqueSize k) const;

    // remove the nodes from the tree in the leaf, assume node must not in other leaf;
    // the order of nodes must be the same as the order from leaf to root
    // return new leaf after removal
    TreeNode *removeNodes(TreeNode *leaf, daf::StaticVector<daf::Size> &nodes);

    // remove all node that only belong to the current leaf;

    void removeLeaf(TreeNode *leaf);

    void serialize(const std::string &filename) const {
        // auto filename = "/data/_/test.bin";
        std::ofstream ofs(filename, std::ios::binary);
        if (!ofs) {
            std::cerr << "Failed to serialize Tree, failed to open " << filename << std::endl;
            return;
        }
        std::cout << "Tree serialized to " << filename << std::endl;

        try {
            boost::archive::binary_oarchive oa(ofs);
            oa << *this;
        } catch (const boost::archive::archive_exception &e) {
            std::cerr << "Archive exception: " << e.what() << std::endl;
            return;
        }
        std::cout << "serialize success" << std::endl;
    }

    static MultiBranchTree deserialize(const std::string &filename) {
        MultiBranchTree newTree;
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs) {
            std::cerr << "Failed to serialize Tree, failed to open " << filename << std::endl;
        }
        boost::archive::binary_iarchive ia(ifs);
        ia >> newTree;
        TreeNode::maxCliques = newTree.root->MaxDeep;
        return newTree;
    }

    /**
     * @brief 计算树中叶子节点的数量。
     *
     * 使用内部 lambda 表达式进行递归遍历，从根节点开始统计所有叶子节点（即没有子节点的节点）。
     *
     * @return 叶子节点总数。
     */
    [[nodiscard]] double numLeafs() const;

    /**
     * 初始化叶子节点的parent指针, 和叶子节点的ID
     */
    [[nodiscard]] daf::Size initLeafsParentAndId(daf::CliqueSize minK = 0) const;

    [[nodiscard]] daf::StaticVector<TreeNode *> getLeafsList(daf::Size maxId, daf::CliqueSize minK = 0) const;

    //
    //
    // [[nodiscard]] TreeNode **getLeafsListRowData(const daf::Size maxId) const {
    //     auto *leafs = new TreeNode *[maxId];
    //     // daf::StaticVector<TreeNode*> leafs(maxId);
    //     // leafs.c_size = maxId;
    //     std::function<void(TreeNode *)> rec = [&](TreeNode *node) {
    //         if (node->children.empty()) leafs[node->leafId] = node;
    //         for (const auto &child: node->children) rec(child);
    //     };
    //
    //     rec(root);
    //
    //     return leafs;
    // }

    static void printCliqueCountByLeafsList(daf::StaticVector<TreeNode *> leafs, const daf::Size maxId) {
        std::vector<double> counts(TreeNode::maxCliques + 1, 0);
        for (daf::Size i = 0; i < maxId; i++) {
            auto node = leafs[i];
            daf::Size pivotCount = 0;
            daf::Size nonPivotCount = 0;
            while (node != nullptr) {
                if (node->isPivot) pivotCount++;
                else nonPivotCount++;
                node = node->parent;
                if (node->v == ROOTID) {
                    for (daf::Size j = 0; j <= pivotCount; j++) {
                        daf::Size k = pivotCount + nonPivotCount - j;
                        counts[k] += nCr[pivotCount][j];
                    }
                    break;
                }
            }
        }

        std::cout << counts << std::endl;
    }


    static void printCliques(daf::StaticVector<TreeNode *> leafs, const daf::Size minK = 0) {
        daf::StaticVector<daf::Size> cliques;
        for (auto node: leafs) {
            if (node->MaxDeep < minK) {
                continue;
            }
            cliques.push_back(node->v);
            auto parent = node->parent;
            while (!parent->isRoot()) {
                cliques.push_back(parent->v);
                parent = parent->parent;
            }
            std::ranges::sort(cliques);
            std::cout << "leaf: " << node->leafId << " cliques: " << cliques << std::endl;
            cliques.clear();
        }
        cliques.free();
    }


    ~MultiBranchTree() {
        // std::cout << "delete tree" << std::endl;
        delete root;
    }

    // static void printCliqueCountByLeafsListRow(TreeNode **leafs, daf::Size maxId);

private:
    void cliqueCountHelper(TreeNode *node, daf::CliqueSize pivotCount, daf::CliqueSize nonPivotCount,
                           std::vector<double> &cliqueCounts) const;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int) {
        ar & root;
    }

    // void cliqueCountHelper(TreeNode* node, int pivotCount, int nonPivotCount, std::vector<daf::CliqueSize> &cliqueCounts) const;
};

#endif //TREE_H
