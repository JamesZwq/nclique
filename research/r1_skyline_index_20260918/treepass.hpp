#pragma once
// The canonical tree of one size and the own node of every vertex (2026-09-23).
//
// Operation for operation the procedure of make_tree_row in count.cpp: the vertices of positive value sorted by
// decreasing value with the same std::sort call, one level per value, the same union-find calls in the same order
// (union by size, ties to the first root, path halving), the same child lists and the same node ids (a level's nodes in
// increasing order of their final root).  The chain index built from it is therefore byte-identical; chain_index_tool
// checks it against make_tree_row under CHAIN_TREE_CHECK=1, and tail_check / the .cx comparison check the files.
// What changes is the cost: the per-size arrays over the vertices and the rows are allocated once and reset through
// the entries a size touched, the pending child lists are linked lists over node ids instead of one vector per vertex,
// node values stay in the count type, and everything a row visit reads or writes (size interval, holds, active counts,
// representative) sits in one 20-byte record per row; the row itself is read only when it becomes live.
#include "../r1_terminal_20260918/terminal.hpp"
#include <span>

namespace treepass {
using namespace allsize;

template<class T> class Pass {
public:
    // the last size's tree: nodes in creation order, parent (-1 for a root), top value, children (CSR, creation order)
    std::vector<int> parent;
    std::vector<T> top;
    std::vector<uint32_t> child_off;   // children of node x: child[child_off[x] .. child_off[x + 1])
    std::vector<int> child;
    std::vector<int> leaf;             // own node of every vertex (-1 when its value is 0)

    Pass(const terminal::Index& index, Vertex n)
        : leaf(n, -1), index_(index), n_(n), dsu_(n), size_(n, 1), active_(n, 0), cur_(n, -1), head_(n, -1), tail_(n, -1),
          state_(index.rows.size()) {
        std::iota(dsu_.begin(), dsu_.end(), 0);
        for (size_t p = 0; p < index.rows.size(); ++p) {
            const auto& row = index.rows[p];
            require(row.hi < 65536 && row.lo < 65536, "size interval exceeds 16 bits");
            state_[p].lo = static_cast<uint16_t>(row.lo); state_[p].hi = static_cast<uint16_t>(row.hi); state_[p].holds = row.holds();
        }
    }

    void run(std::span<const T> core, int s) {
        require(core.size() == n_, "row size");
        for (Vertex v : order_) leaf[v] = -1;                      // the previous size's own nodes
        parent.clear(); top.clear(); child.clear(); child_off.assign(1, 0); next_.clear();
        order_.clear();
        for (Vertex v = 0; v < n_; ++v) if (core[v] > 0) order_.push_back(v);
        std::sort(order_.begin(), order_.end(), [&](Vertex a, Vertex b) { return core[a] > core[b]; });
        const uint32_t sv = static_cast<uint32_t>(s);
        const Vertex* const members = index_.members.data();
        size_t at = 0;
        while (at < order_.size()) {
            const T value = core[order_[at]];
            size_t end = at;
            while (end < order_.size() && core[order_[end]] == value) ++end;
            ++stamp_; touched_.clear();
            for (size_t z = at; z < end; ++z) {
                const Vertex v = order_[z];
                active_[v] = 1; touched_.push_back(find(v));
                for (uint64_t code : index_.touching(v)) {
                    const size_t p = code >> 2; const unsigned role = code & 3;
                    RowState& st = state_[p];
                    if (st.lo > sv || st.hi < sv) continue;                            // row not valid at s
                    if (!st.ah && !st.aq && st.rep < 0) rows_.push_back(p);
                    if (role == 0) ++st.ah; else if (role == 1) ++st.aq;
                    if (st.rep < 0) {
                        if (st.ah == st.holds && st.holds + st.aq >= sv) {             // the row becomes live
                            st.rep = static_cast<int32_t>(v);
                            const auto& row = index_.rows[p];
                            for (terminal::Offset i = row.begin; i < row.pivot_end; ++i) if (active_[members[i]]) unite(v, members[i]);
                        }
                    } else unite(v, static_cast<Vertex>(st.rep));
                }
            }
            std::sort(touched_.begin(), touched_.end());
            touched_.erase(std::unique(touched_.begin(), touched_.end()), touched_.end());
            for (Vertex old : touched_) {
                const Vertex r = find(old);
                if (r != old) continue;
                const int id = static_cast<int>(parent.size());
                parent.push_back(-1); top.push_back(value);
                for (int c = head_[r]; c >= 0; c = next_[c]) { child.push_back(c); parent[c] = id; }
                child_off.push_back(static_cast<uint32_t>(child.size()));
                head_[r] = tail_[r] = -1; cur_[r] = id;
            }
            for (size_t z = at; z < end; ++z) leaf[order_[z]] = cur_[find(order_[z])];
            at = end;
        }
        for (Vertex v : order_) require(leaf[v] >= 0 && top[leaf[v]] == core[v], "own node level");
        for (size_t x = 0; x < parent.size(); ++x) if (parent[x] >= 0) require(top[parent[x]] < top[x], "parent order");
        // reset what this size touched
        for (Vertex v : order_) { dsu_[v] = v; size_[v] = 1; active_[v] = 0; cur_[v] = -1; head_[v] = tail_[v] = -1; }
        for (size_t p : rows_) { state_[p].ah = 0; state_[p].aq = 0; state_[p].rep = -1; }
        rows_.clear();
    }
    size_t nodes() const { return parent.size(); }

private:
    Vertex find(Vertex x) { while (dsu_[x] != x) { dsu_[x] = dsu_[dsu_[x]]; x = dsu_[x]; } return x; }
    Vertex join(Vertex a, Vertex b) {
        a = find(a); b = find(b); if (a == b) return a;
        if (size_[a] < size_[b]) std::swap(a, b);
        dsu_[b] = a; size_[a] += size_[b]; return a;
    }
    void add_child(Vertex root, int node) {                        // once per node and level, in call order
        if (node < 0) return;
        if (marked_.size() <= static_cast<size_t>(node)) marked_.resize(static_cast<size_t>(node) + 1, 0);
        if (marked_[node] == stamp_) return;
        marked_[node] = stamp_;
        if (next_.size() <= static_cast<size_t>(node)) next_.resize(static_cast<size_t>(node) + 1, -1);
        next_[node] = -1;
        if (head_[root] < 0) head_[root] = node; else next_[tail_[root]] = node;
        tail_[root] = node;
    }
    Vertex unite(Vertex a, Vertex b) {
        a = find(a); b = find(b); if (a == b) return a;
        add_child(a, cur_[a]); add_child(b, cur_[b]);
        const Vertex r = join(a, b), other = r == a ? b : a;
        if (head_[other] >= 0) {                                   // append the other root's pending children
            if (head_[r] < 0) head_[r] = head_[other]; else next_[tail_[r]] = head_[other];
            tail_[r] = tail_[other]; head_[other] = tail_[other] = -1;
        }
        cur_[r] = -1; touched_.push_back(find(r)); return r;
    }

    const terminal::Index& index_;
    Vertex n_;
    Vertices dsu_, size_;
    std::vector<uint8_t> active_;
    std::vector<int> cur_, head_, tail_, next_;
    struct RowState { uint16_t lo = 0, hi = 0; uint32_t holds = 0, ah = 0, aq = 0; int32_t rep = -1; };   // rep >= 0: live
    std::vector<RowState> state_;
    std::vector<uint64_t> marked_;
    uint64_t stamp_ = 0;
    Vertices order_, touched_;
    std::vector<size_t> rows_;
};
}
