#pragma once
// hierarchy.hpp: per-size canonical trees over the all-size solver (make_tree_row), the brute-force helpers used by the
// selftest (CountDSU, complete, split_graph) and the input preparation of the solver harness.  Taken from the research
// line r1_skyline_index_20260918 (count.cpp) without its stage-1 counting program.
#include "solver/terminal.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sys/resource.h>

using namespace orderreplay;
#include "solver/shared_harness.inc"

struct CountDSU {
    Vertices p, sz;
    explicit CountDSU(Vertex n):p(n),sz(n,1){std::iota(p.begin(),p.end(),0);}
    Vertex find(Vertex x){while(p[x]!=x){p[x]=p[p[x]];x=p[x];}return x;}
    Vertex join(Vertex a,Vertex b){a=find(a);b=find(b);if(a==b)return a;if(sz[a]<sz[b])std::swap(a,b);p[b]=a;sz[a]+=sz[b];return a;}
};

static cpp_int choose_int(unsigned n,unsigned r) { return binomial(n,r); }
static cpp_int shadow(int s, cpp_int k) {
    if(!k)return 0; cpp_int ans=0;
    for(int r=s;r>=1 && k>0;--r) {
        unsigned lo=r,hi=r;
        while(choose_int(hi,r)<=k) { if(hi>1000000) throw std::overflow_error("cascade search"); hi*=2; }
        while(lo+1<hi) { unsigned mid=lo+(hi-lo)/2; if(choose_int(mid,r)<=k)lo=mid;else hi=mid; }
        k-=choose_int(lo,r); ans+=choose_int(lo,r-1);
    }
    return ans;
}

struct Node { cpp_int hi; int parent=-1, creator=-1; std::vector<int> children; };
template<class T> struct Tree {
    std::vector<Node> nodes; std::vector<int> leaf; uint64_t l1=0;
    uint64_t chains=0; uint64_t depth=0;
    // Only selftest consumes these snapshots; they are component labels after a core level.
    std::map<cpp_int,std::vector<int>,std::greater<cpp_int>> snapshot;
};

template<class T> static Tree<T> make_tree_row(const Graph& g,const terminal::Index& index,
        std::span<const T> core,int s,bool save_snapshots=false) {
    const Vertex n=g.n; Tree<T> out; out.leaf.assign(n,-1); require(core.size()==n,"row size");
    std::vector<Vertex> order; for(Vertex v=0;v<n;++v)if(core[v]>0)order.push_back(v);
    std::sort(order.begin(),order.end(),[&](Vertex a,Vertex b){return core[a]>core[b];});
    CountDSU uf(n); std::vector<uint8_t> active(n),live(index.rows.size());
    std::vector<uint32_t> ah(index.rows.size()),aq(index.rows.size()); std::vector<int> rep(index.rows.size(),-1),cur(n,-1);
    std::vector<std::vector<int>> pending(n); std::vector<uint64_t> marked; uint64_t stamp=0;
    auto add_child=[&](Vertex root,int child){ if(child<0)return; if(marked.size()<=static_cast<size_t>(child))marked.resize(child+1); if(marked[child]!=stamp){marked[child]=stamp;pending[root].push_back(child);} };
    size_t at=0;
    while(at<order.size()) {
        const T value=core[order[at]]; size_t end=at;while(end<order.size()&&core[order[end]]==value)++end;
        ++stamp; std::vector<Vertex> touched;
        auto touch=[&](Vertex r){r=uf.find(r);touched.push_back(r);};
        auto unite=[&](Vertex a,Vertex b) {
            a=uf.find(a);b=uf.find(b);if(a==b)return a;
            add_child(a,cur[a]);add_child(b,cur[b]);
            Vertex r=uf.join(a,b),other=(r==a?b:a);
            if(!pending[other].empty()){pending[r].insert(pending[r].end(),pending[other].begin(),pending[other].end());pending[other].clear();}
            cur[r]=-1; touch(r); return r;
        };
        for(size_t z=at;z<end;++z) {
            Vertex v=order[z];active[v]=1;touch(v);
            for(uint64_t code:index.touching(v)) {
                const size_t p=code>>2;const unsigned role=code&3;const auto& row=index.rows[p];if(!row.valid(s))continue;
                if(role==0)++ah[p]; else if(role==1)++aq[p];
                if(!live[p] && ah[p]==row.holds() && row.holds()+aq[p]>=static_cast<Vertex>(s)) {
                    live[p]=1;rep[p]=v;
                    for(terminal::Offset i=row.begin;i<row.pivot_end;++i)if(active[index.members[i]])unite(v,index.members[i]);
                } else if(live[p]) unite(v,static_cast<Vertex>(rep[p]));
            }
        }
        std::sort(touched.begin(),touched.end());touched.erase(std::unique(touched.begin(),touched.end()),touched.end());
        std::map<Vertex,Vertex> firstv;for(size_t z=at;z<end;++z){Vertex r=uf.find(order[z]);if(!firstv.contains(r))firstv[r]=order[z];}
        for(Vertex old:touched) { Vertex r=uf.find(old); if(r!=old)continue; require(firstv.contains(r),"touched root without a level vertex"); Node node;node.hi=cpp_int(value);node.creator=firstv[r];node.children=std::move(pending[r]);
            int id=out.nodes.size();for(int c:node.children)out.nodes[c].parent=id;out.nodes.push_back(std::move(node));cur[r]=id; }
        for(size_t z=at;z<end;++z)out.leaf[order[z]]=cur[uf.find(order[z])];
        if(save_snapshots) { std::vector<int> labels(n,-1);std::map<Vertex,int> ids;int next=0;for(Vertex v:order)if(active[v]){Vertex r=uf.find(v);if(!ids.contains(r))ids[r]=next++;labels[v]=ids[r];}out.snapshot.emplace(cpp_int(value),std::move(labels)); }
        at=end;
    }
    for(Vertex v=0;v<n;++v)if(core[v]>0){require(out.leaf[v]>=0,"missing leaf");require(out.nodes[out.leaf[v]].hi==cpp_int(core[v]),"L1 leaf level");++out.l1;}
    for(const auto& node:out.nodes)if(!node.children.empty())for(int c:node.children)require(node.hi<out.nodes[c].hi,"L1 parent order");
    for(const auto& node:out.nodes)if(node.children.size()==1)++out.chains;
    for(size_t x=0;x<out.nodes.size();++x){uint64_t d=1;for(int p=out.nodes[x].parent;p>=0;p=out.nodes[p].parent)++d;out.depth=std::max(out.depth,d);}
    return out;
}

template<class T> static Tree<T> make_tree(const Graph& g,const terminal::Index& index,
        const std::vector<T>& core,int s,bool save_snapshots=false) {
    return make_tree_row<T>(g,index,std::span<const T>(core).subspan(static_cast<size_t>(s)*g.n,g.n),s,save_snapshots);
}

static std::vector<int> twins(const Graph& g,std::vector<std::vector<Vertex>>& groups) {
    std::map<std::vector<Vertex>,int> ids;std::vector<int> c(g.n);
    for(Vertex v=0;v<g.n;++v){std::vector<Vertex> key(g.row(v).begin(),g.row(v).end());key.insert(std::lower_bound(key.begin(),key.end(),v),v);auto [it,newone]=ids.emplace(std::move(key),groups.size());if(newone)groups.emplace_back();c[v]=it->second;groups[c[v]].push_back(v);}return c;
}

static Graph complete(Vertex n){std::vector<std::pair<Vertex,Vertex>> e;for(Vertex a=0;a<n;++a)for(Vertex b=a+1;b<n;++b)e.emplace_back(a,b);return Graph::from_edges(n,std::move(e));}
static Graph split_graph(Vertex h,Vertex c,bool extra=false){Vertex n=h+c+(extra?2:0);std::vector<std::pair<Vertex,Vertex>> e;for(Vertex a=0;a<h;++a)for(Vertex b=a+1;b<n;++b)e.emplace_back(a,b);if(extra)e.emplace_back(n-2,n-1);return Graph::from_edges(n,std::move(e));}
