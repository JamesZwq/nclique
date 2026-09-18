#include "common.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sys/resource.h>

using namespace orderreplay;
#include "../r1_orderreplay_20260917/shared_harness.inc"

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

template<class T> static Tree<T> make_tree(const Graph& g,const terminal::Index& index,
        const std::vector<T>& core,int s,bool save_snapshots=false) {
    const Vertex n=g.n; Tree<T> out; out.leaf.assign(n,-1);
    std::vector<Vertex> order; for(Vertex v=0;v<n;++v)if(core[static_cast<size_t>(s)*n+v]>0)order.push_back(v);
    std::sort(order.begin(),order.end(),[&](Vertex a,Vertex b){return core[static_cast<size_t>(s)*n+a]>core[static_cast<size_t>(s)*n+b];});
    CountDSU uf(n); std::vector<uint8_t> active(n),live(index.rows.size());
    std::vector<uint32_t> ah(index.rows.size()),aq(index.rows.size()); std::vector<int> rep(index.rows.size(),-1),cur(n,-1);
    std::vector<std::vector<int>> pending(n); std::vector<uint64_t> marked; uint64_t stamp=0;
    auto add_child=[&](Vertex root,int child){ if(child<0)return; if(marked.size()<=static_cast<size_t>(child))marked.resize(child+1); if(marked[child]!=stamp){marked[child]=stamp;pending[root].push_back(child);} };
    size_t at=0;
    while(at<order.size()) {
        const T value=core[static_cast<size_t>(s)*n+order[at]]; size_t end=at;while(end<order.size()&&core[static_cast<size_t>(s)*n+order[end]]==value)++end;
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
            for(Vertex code:index.touching(v)) {
                Vertex p=code>>2,role=code&3;const auto& row=index.rows[p];if(!row.valid(s))continue;
                if(role==0)++ah[p]; else if(role==1)++aq[p];
                if(!live[p] && ah[p]==row.holds() && row.holds()+aq[p]>=static_cast<Vertex>(s)) {
                    live[p]=1;rep[p]=v;
                    for(Vertex i=row.begin;i<row.pivot_end;++i)if(active[index.members[i]])unite(v,index.members[i]);
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
    for(Vertex v=0;v<n;++v)if(core[static_cast<size_t>(s)*n+v]>0){require(out.leaf[v]>=0,"missing leaf");require(out.nodes[out.leaf[v]].hi==cpp_int(core[static_cast<size_t>(s)*n+v]),"L1 leaf level");++out.l1;}
    for(const auto& node:out.nodes)if(!node.children.empty())for(int c:node.children)require(node.hi<out.nodes[c].hi,"L1 parent order");
    for(const auto& node:out.nodes)if(node.children.size()==1)++out.chains;
    for(size_t x=0;x<out.nodes.size();++x){uint64_t d=1;for(int p=out.nodes[x].parent;p>=0;p=out.nodes[p].parent)++d;out.depth=std::max(out.depth,d);}
    return out;
}

static std::vector<int> twins(const Graph& g,std::vector<std::vector<Vertex>>& groups) {
    std::map<std::vector<Vertex>,int> ids;std::vector<int> c(g.n);
    for(Vertex v=0;v<g.n;++v){std::vector<Vertex> key(g.row(v).begin(),g.row(v).end());key.insert(std::lower_bound(key.begin(),key.end(),v),v);auto [it,newone]=ids.emplace(std::move(key),groups.size());if(newone)groups.emplace_back();c[v]=it->second;groups[c[v]].push_back(v);}return c;
}

struct CountStats { uint64_t n=0,m=0,d=0,smax=0,bits=0,ncls=0,maxcls=0,apv=0,apc=0,ev=0,ec=0,nt=0,ntmax=0,empty=0,chain=0,depth=0,residue=0,reszero=0,cert=0,f2=0,f2bad=0,l1=0,l2=0; double av=0,mu=0,avc=0,muc=0; };

template<class T> static CountStats measure(const Input& in,unsigned bits,bool test=false,uint64_t* partition_checks=nullptr) {
    const Graph& g=in.graph;const int maximum=std::max(2,static_cast<int>(in.d)+1);Layout layout(g,maximum);layout.prepare(g.n);
    typename Kernel<T>::Combinations choose(in.d+1,maximum); terminal::Index index(maximum);terminal::build(g,index,0);index.prepare(g.n);
    auto out=terminal::Solver<T>::solve(g,index,choose,in.ordinary);auto fixed=Kernel<T>{}.template fixed_sparse<true>(layout,g.n,choose,in.ordinary);require(out.common.data.core==fixed.core,"terminal core differs from frozen fixed_sparse control");const auto& core=out.common.data.core;
    std::vector<std::vector<Vertex>> groups;auto cls=twins(g,groups);CountStats q;q.n=g.n;q.m=g.m;q.d=in.d;q.smax=maximum;q.bits=bits;q.ncls=groups.size();for(auto& x:groups)q.maxcls=std::max<uint64_t>(q.maxcls,x.size());
    for(auto& x:groups)for(int s=2;s<=maximum;++s)for(Vertex u:x)require(core[static_cast<size_t>(s)*g.n+u]==core[static_cast<size_t>(s)*g.n+x[0]],"F5 twin core equality");
    std::vector<int> omega(g.n),sigma(g.n);std::vector<std::vector<uint8_t>> sky(g.n,std::vector<uint8_t>(maximum+1));std::map<std::pair<int,cpp_int>,cpp_int> cache;
    auto sig=[&](int s,const T& v)->cpp_int{auto key=std::make_pair(s,cpp_int(v));auto it=cache.find(key);if(it!=cache.end())return it->second;return cache.emplace(key,shadow(s,key.second)).first->second;};
    for(Vertex v=0;v<g.n;++v){for(int s=2;s<=maximum;++s)if(core[static_cast<size_t>(s)*g.n+v]>0)omega[v]=s;for(int s=2;s<=omega[v];++s)require(core[static_cast<size_t>(s)*g.n+v]>0,"F1 support nesting");sigma[v]=omega[v]+1;for(int s=2;s<=omega[v];++s)if(cpp_int(core[static_cast<size_t>(s)*g.n+v])==choose_int(omega[v]-1,s-1)&&sigma[v]==omega[v]+1)sigma[v]=s;for(int s=sigma[v];s<=omega[v];++s)require(cpp_int(core[static_cast<size_t>(s)*g.n+v])==choose_int(omega[v]-1,s-1),"F8 certified tail");
        for(int s=2;s<omega[v];++s){auto a=cpp_int(core[static_cast<size_t>(s)*g.n+v]),b=sig(s,core[static_cast<size_t>(s+1)*g.n+v]);++q.f2;if(a<b){++q.f2bad;require(false,"F2 shadow bound");}bool zero=a==b;sky[v][s]=!zero;if(s<sigma[v]){++q.residue;if(zero)++q.reszero;}}if(omega[v]>=2)sky[v][omega[v]]=1; if(omega[v]>=2&&sigma[v]==omega[v]+1)++q.residue;q.apv+=omega[v]>=2?omega[v]-1:0;q.ev+=std::count(sky[v].begin(),sky[v].end(),uint8_t{1});}
    for(size_t c=0;c<groups.size();++c){Vertex v=groups[c][0];q.apc+=omega[v]>=2?omega[v]-1:0;q.ec+=std::count(sky[v].begin(),sky[v].end(),uint8_t{1});}
    uint64_t active=0,activec=0;for(Vertex v=0;v<g.n;++v)active+=omega[v]>=2;for(auto& z:groups)activec+=omega[z[0]]>=2;q.av=active?double(q.apv)/active:0;q.mu=active?double(q.ev)/active:0;q.avc=activec?double(q.apc)/activec:0;q.muc=activec?double(q.ec)/activec:0;q.cert=q.apv-q.residue;
    std::vector<Tree<T>> trees;trees.reserve(maximum-1);for(int s=2;s<=maximum;++s){trees.push_back(make_tree(g,index,core,s,test));auto& tr=trees.back();q.nt+=tr.nodes.size();q.ntmax=std::max<uint64_t>(q.ntmax,tr.nodes.size());q.chain+=tr.chains;q.depth=std::max(q.depth,tr.depth);q.l1+=tr.l1;}
    for(int s=2;s<=maximum;++s){auto& tr=trees[s-2];std::vector<uint8_t> used(tr.nodes.size());for(Vertex v=0;v<g.n;++v)if(sky[v][s]&&tr.leaf[v]>=0)used[tr.leaf[v]]=1;for(auto x:used)q.empty+=!x;}
    for(int s=2;s<maximum;++s){auto& a=trees[s-2];auto& b=trees[s-1];for(Vertex v=0;v<g.n;++v)if(s<omega[v]&&cpp_int(core[static_cast<size_t>(s)*g.n+v])==sig(s,core[static_cast<size_t>(s+1)*g.n+v])){int node=b.leaf[v];Vertex u=static_cast<Vertex>(b.nodes[node].creator);cpp_int level=shadow(s,b.nodes[node].hi);int target=a.leaf[u];while(a.nodes[target].parent>=0&&a.nodes[a.nodes[target].parent].hi>=level)target=a.nodes[target].parent;require(target==a.leaf[v],"L2 chain");++q.l2;}}
    if(test){ for(int s=2;s<=maximum;++s){auto& tr=trees[s-2];auto cl=bottomup::clique_masks(g,s);for(const auto& [level,got]:tr.snapshot){std::vector<uint8_t> inside(g.n);for(Vertex v=0;v<g.n;++v)inside[v]=core[static_cast<size_t>(s)*g.n+v]>=static_cast<T>(level);CountDSU brute(g.n);for(auto mask:cl) {bool ok=true;for(Vertex v=0;v<g.n;++v)if((mask>>v)&1)ok&=inside[v];if(ok){Vertex first=absent;for(Vertex v=0;v<g.n;++v)if((mask>>v)&1){if(first==absent)first=v;else brute.join(first,v);}}} std::map<Vertex,int> ids;std::map<int,int> canon;int next=0,nextc=0;for(Vertex v=0;v<g.n;++v)if(inside[v]){Vertex r=brute.find(v);if(!ids.contains(r))ids[r]=next++;require(got[v]>=0,"inactive vertex inside core");if(!canon.contains(got[v]))canon[got[v]]=nextc++;require(canon[got[v]]==ids[r],"connectivity partition mismatch");}else require(got[v]<0,"active vertex outside core");if(partition_checks)++*partition_checks;}
        // Canonical interval check: every constructed node has exactly the brute component's contiguous level interval.
        std::map<std::vector<Vertex>,std::pair<uint64_t,uint64_t>> bruteints,ours;
        for(const auto& [lv,unused]:tr.snapshot){const uint64_t k=lv.template convert_to<uint64_t>();std::vector<uint8_t> inside(g.n);for(Vertex v=0;v<g.n;++v)inside[v]=core[static_cast<size_t>(s)*g.n+v]>=k;CountDSU dsu(g.n);for(auto mask:cl){bool ok=true;for(Vertex v=0;v<g.n;++v)if((mask>>v)&1)ok&=inside[v];if(ok){Vertex first=absent;for(Vertex v=0;v<g.n;++v)if((mask>>v)&1){if(first==absent)first=v;else dsu.join(first,v);}}}std::map<Vertex,std::vector<Vertex>> parts;for(Vertex v=0;v<g.n;++v)if(inside[v])parts[dsu.find(v)].push_back(v);for(auto& [r,set]:parts){auto& p=bruteints[set];if(!p.first||k<p.first)p.first=k;p.second=std::max(p.second,k);}}
        for(const auto& [level,labels]:tr.snapshot){std::map<int,std::vector<Vertex>> parts;for(Vertex v=0;v<g.n;++v)if(labels[v]>=0)parts[labels[v]].push_back(v);for(auto& [r,set]:parts){auto& p=ours[set];uint64_t k=level.template convert_to<uint64_t>();if(!p.first||k<p.first)p.first=k;p.second=std::max(p.second,k);}}require(bruteints==ours,"canonical interval mismatch"); }
    }
    return q;
}

static Graph complete(Vertex n){std::vector<std::pair<Vertex,Vertex>> e;for(Vertex a=0;a<n;++a)for(Vertex b=a+1;b<n;++b)e.emplace_back(a,b);return Graph::from_edges(n,std::move(e));}
static Graph split_graph(Vertex h,Vertex c,bool extra=false){Vertex n=h+c+(extra?2:0);std::vector<std::pair<Vertex,Vertex>> e;for(Vertex a=0;a<h;++a)for(Vertex b=a+1;b<n;++b)e.emplace_back(a,b);if(extra)e.emplace_back(n-2,n-1);return Graph::from_edges(n,std::move(e));}
template<class T> static void test_graph_inner(const Graph& g,uint64_t& cells,uint64_t& parts);
template<class T> static void test_graph(const Graph& g,uint64_t& cells,uint64_t& parts){
    try{test_graph_inner<T>(g,cells,parts);}catch(const std::exception& e){std::cerr<<"selftest failure on n="<<g.n<<" edges:";for(Vertex u=0;u<g.n;++u)for(Vertex w:g.row(u))if(u<w)std::cerr<<' '<<u<<'-'<<w;std::cerr<<'\n';throw;}}
template<class T> static void test_graph_inner(const Graph& g,uint64_t& cells,uint64_t& parts){Seeds z(g);Input in{g,z.ordinary,z.maximum};int mx=std::max(2,static_cast<int>(z.maximum)+1);const unsigned b=64;auto r=measure<T>(in,b,true,&parts);cells+=static_cast<uint64_t>(mx+1)*g.n;require(r.f2bad==0,"selftest F2 violation");}
// Kruskal-Katona is tight on colex initial segments: the first k s-subsets of {0..11}
// in colex order have exactly shadow(s,k) distinct (s-1)-subsets.
static void colex_tests(){for(int s=2;s<=4;++s){std::vector<std::vector<int>> subs;std::vector<int> idx(s);std::iota(idx.begin(),idx.end(),0);
    while(true){subs.push_back(idx);int i=s-1;while(i>=0&&idx[i]==12-s+i)--i;if(i<0)break;++idx[i];for(int j=i+1;j<s;++j)idx[j]=idx[j-1]+1;}
    std::sort(subs.begin(),subs.end(),[](const std::vector<int>&a,const std::vector<int>&b){for(int i=a.size()-1;i>=0;--i)if(a[i]!=b[i])return a[i]<b[i];return false;});
    for(int k=1;k<=60;++k){std::set<std::vector<int>> sh;for(int t=0;t<k;++t)for(int drop=0;drop<s;++drop){std::vector<int> sub;for(int j=0;j<s;++j)if(j!=drop)sub.push_back(subs[t][j]);sh.insert(sub);}
        require(cpp_int(sh.size())==shadow(s,k),"colex shadow brute force");}}}
static void shadow_tests(){for(int s=2;s<=6;++s)for(int k=0;k<=200;++k)require(shadow(s,k)<=shadow(s,k+1),"shadow monotonicity");for(int s=2;s<=6;++s)for(int a=s;a<=40;++a)require(shadow(s,choose_int(a,s))==choose_int(a,s-1),"single cascade shadow");require(shadow(2,3)==3&&shadow(2,4)==4&&shadow(2,21)==7&&shadow(3,35)==21&&shadow(3,21)==17&&shadow(3,22)==18,"shadow examples");colex_tests();}
static void selftest(){shadow_tests();uint64_t graphs=0,cells=0,parts=0;std::mt19937_64 rng(20260918);for(Vertex n=0;n<=6;++n){std::vector<std::pair<Vertex,Vertex>> p;for(Vertex a=0;a<n;++a)for(Vertex b=a+1;b<n;++b)p.emplace_back(a,b);for(uint64_t mask=0;mask<(uint64_t{1}<<p.size());++mask){std::vector<std::pair<Vertex,Vertex>> e;for(size_t i=0;i<p.size();++i)if(mask>>i&1)e.push_back(p[i]);test_graph<uint64_t>(Graph::from_edges(n,std::move(e)),cells,parts);++graphs;}}for(int z=0;z<200;++z){Vertex n=7+rng()%4;std::vector<std::pair<Vertex,Vertex>> e;for(Vertex a=0;a<n;++a)for(Vertex b=a+1;b<n;++b)if(rng()%2)e.emplace_back(a,b);test_graph<uint64_t>(Graph::from_edges(n,std::move(e)),cells,parts);++graphs;}for(bool x:{false,true})for(Vertex h:{1,2,4}){test_graph<uint64_t>(split_graph(h,4,x),cells,parts);++graphs;}test_graph<uint64_t>(complete(8),cells,parts);++graphs;std::cout<<"{\"passed\":true,\"graphs\":"<<graphs<<",\"cells\":"<<cells<<",\"partitions_compared\":"<<parts<<"}\n";}
static void print(const CountStats& x){uint64_t ws=x.apc+4*x.nt,wi=x.ec+x.ec/4+6*x.nt;std::cout<<std::fixed<<std::setprecision(6)<<"{\"passed\":true,\"n\":"<<x.n<<",\"m\":"<<x.m<<",\"d\":"<<x.d<<",\"s_max\":"<<x.smax<<",\"count_bits\":"<<x.bits<<",\"n_cls\":"<<x.ncls<<",\"max_class_size\":"<<x.maxcls<<",\"active_pairs_vertex\":"<<x.apv<<",\"active_pairs_class\":"<<x.apc<<",\"E_vertex\":"<<x.ev<<",\"E_class\":"<<x.ec<<",\"avg_traj\":"<<x.av<<",\"mu\":"<<x.mu<<",\"avg_traj_class\":"<<x.avc<<",\"mu_class\":"<<x.muc<<",\"N_T\":"<<x.nt<<",\"N_T_max\":"<<x.ntmax<<",\"empty_nodes\":"<<x.empty<<",\"chain_nodes\":"<<x.chain<<",\"max_depth\":"<<x.depth<<",\"residue_cells\":"<<x.residue<<",\"residue_zero_delta\":"<<x.reszero<<",\"certified_cells\":"<<x.cert<<",\"F2_checked\":"<<x.f2<<",\"F2_violations\":"<<x.f2bad<<",\"L1_checks\":"<<x.l1<<",\"L2_checks\":"<<x.l2<<",\"words_S_trees\":"<<ws<<",\"words_skyline\":"<<wi<<",\"lhs\":"<<x.ncls*(x.avc-x.muc)<<",\"rhs\":"<<2*x.nt+x.muc*x.ncls/4<<",\"verdict\":\""<<(wi<ws?"skyline smaller":"S trees")<<"\"}\n";}
int main(int argc,char** argv){try{if(argc==2&&std::string(argv[1])=="--selftest"){selftest();return 0;}require(argc==3&&std::string(argv[1])=="--graph","usage: count --selftest | --graph path");Input in=prepare(argv[2]);int maximum=std::max(2,static_cast<int>(in.d)+1);Layout l(in.graph,maximum);l.prepare(in.graph.n);unsigned b=width(count_bound(in.graph,l,in.d));if(b==64)print(measure<uint64_t>(in,b));else if(b==128)print(measure<unsigned __int128>(in,b));else if(b==256)print(measure<boost::multiprecision::uint256_t>(in,b));else print(measure<boost::multiprecision::uint512_t>(in,b));}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
