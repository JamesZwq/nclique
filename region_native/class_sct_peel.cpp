// class_sct_peel.cpp — PEEL-SUITABILITY probe for the weighted class SCT.
//
// Question: build the s-clique SCT as disjoint CCPath leaves. For a chosen
// r-clique PATTERN P (a multiplicity vector, sum=r, support is a clique),
// insert P as a forbidden dead-box into EVERY leaf, then
//   sum_leaves support_count(leaf, b=0)
// should equal the number of weighted s-cliques that DO NOT "contain" P,
// where "contain" must be defined at the actual-vertex level to be a faithful
// nucleus peel. We compare against an EXPANDED simple-graph oracle.
//
// Build (from region_native/):
//   g++ -O3 -std=c++17 -I../src/NucleusDecomposition -o class_sct_peel class_sct_peel.cpp

#include "../src/NucleusDecomposition/CCPathCore.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <random>
#include <vector>

using ccpath::CCPath;
using ccpath::Vec;

static double PASCAL[129][129];
static void initPascal(){
    for(int n=0;n<=128;++n){PASCAL[n][0]=1.0;
        for(int k=1;k<=n;++k)PASCAL[n][k]=PASCAL[n-1][k-1]+(k<=n-1?PASCAL[n-1][k]:0.0);}
}
static double nCr_fn(int n,int k){ if(k<0||n<0||k>n)return 0.0; return PASCAL[n][k]; }

struct ClassG { int C=0; std::vector<int> w; std::vector<std::vector<char>> A;
    bool adj(int i,int j)const{return A[i][j];} };
struct OpenC { int c; int wres; };

// ---- builder (verbatim from class_sct.cpp) ----
static void emitLeaf(const ClassG& G, const std::vector<std::pair<int,int>>& spine,
                     const std::vector<OpenC>& pool, int k, std::vector<CCPath>& out){
    const int C=G.C; Vec h=ccpath::zeros_vec(C),n=ccpath::zeros_vec(C),
        ell=ccpath::zeros_vec(C),u=ccpath::zeros_vec(C);
    int forced=0;
    for(auto&sp:spine){int c=sp.first,mult=sp.second;n[c]=(int16_t)G.w[c];
        ell[c]=(int16_t)mult;u[c]=(int16_t)mult;forced+=mult;}
    int poolCap=0;
    for(auto&pc:pool){int c=pc.c;n[c]=(int16_t)G.w[c];u[c]=(int16_t)G.w[c];poolCap+=(int)G.w[c];}
    CCPath p; p.h=h;p.n=n;p.ell=ell;p.u=u;p.T=k;
    if(forced>k)return; if(forced+poolCap<k)return; out.push_back(std::move(p));
}
static void gen(const ClassG& G, std::vector<std::pair<int,int>>& spine,int spineSum,
                std::vector<OpenC>& pool,int poolSum,std::vector<OpenC> P,int k,
                std::vector<CCPath>& out){
    if(spineSum>k)return;
    if(P.empty()){emitLeaf(G,spine,pool,k,out);return;}
    int bestIdx=0;long bestScore=-1;
    for(size_t i=0;i<P.size();++i){long sc=0;
        for(size_t j=0;j<P.size();++j) if(j!=i&&G.adj(P[i].c,P[j].c)) sc+=P[j].wres;
        if(sc>bestScore){bestScore=sc;bestIdx=(int)i;}}
    int pc=P[bestIdx].c;
    bool universal=true;
    for(auto&q:P) if(q.c!=pc&&!G.adj(pc,q.c)){universal=false;break;}
    if(universal){std::vector<OpenC> Pr;Pr.reserve(P.size()-1);
        for(auto&q:P) if(q.c!=pc)Pr.push_back(q);
        pool.push_back(OpenC{pc,G.w[pc]});
        gen(G,spine,spineSum,pool,poolSum+G.w[pc],std::move(Pr),k,out);
        pool.pop_back();return;}
    std::vector<int> nonNb;
    for(auto&q:P) if(q.c!=pc&&!G.adj(pc,q.c))nonNb.push_back(q.c);
    std::sort(nonNb.begin(),nonNb.end());
    {std::vector<OpenC> Pp;Pp.reserve(P.size());
        for(auto&q:P) if(q.c==pc||G.adj(pc,q.c))Pp.push_back(q);
        gen(G,spine,spineSum,pool,poolSum,std::move(Pp),k,out);}
    for(size_t t=0;t<nonNb.size();++t){int v=nonNb[t];int wv=0;
        for(auto&q:P) if(q.c==v){wv=q.wres;break;}
        std::vector<char> drop(G.C,0); for(size_t q=0;q<=t;++q)drop[nonNb[q]]=1;
        std::vector<OpenC> base;base.reserve(P.size());
        for(auto&q:P){if(drop[q.c])continue; if(!G.adj(v,q.c))continue; base.push_back(q);}
        for(int mtt=1;mtt<=wv;++mtt){ if(spineSum+mtt>k)break;
            spine.push_back({v,mtt}); std::vector<OpenC> child=base;
            gen(G,spine,spineSum+mtt,pool,poolSum,std::move(child),k,out); spine.pop_back();}}
}
static std::vector<CCPath> buildClassSCT(const ClassG& G,int k){
    std::vector<CCPath> out; std::vector<OpenC> P;
    for(int c=0;c<G.C;++c) if(G.w[c]>0)P.push_back(OpenC{c,G.w[c]});
    std::sort(P.begin(),P.end(),[](const OpenC&a,const OpenC&b){return a.c<b.c;});
    std::vector<std::pair<int,int>> spine; std::vector<OpenC> pool;
    gen(G,spine,0,pool,0,std::move(P),k,out); return out;
}
static double sctTotalForbidden(const std::vector<CCPath>& leaves){
    double t=0.0; for(const auto&L:leaves){Vec b=ccpath::zeros_vec(L.m());
        t+=ccpath::support_count(L,b,nCr_fn);} return t;
}

// ---- expanded simple graph for an honest oracle ----
struct Expanded {
    int N=0; std::vector<int> vclass; std::vector<std::vector<char>> Adj;
    std::vector<std::vector<int>> classVerts; // verts of each class
};
static Expanded expand(const ClassG& G){
    Expanded E; E.classVerts.assign(G.C,{});
    for(int c=0;c<G.C;++c) for(int t=0;t<G.w[c];++t){ E.classVerts[c].push_back(E.N); E.vclass.push_back(c); E.N++; }
    E.Adj.assign(E.N,std::vector<char>(E.N,0));
    for(int i=0;i<E.N;++i)for(int j=0;j<E.N;++j){ if(i==j)continue;
        if(E.vclass[i]==E.vclass[j])E.Adj[i][j]=1; else E.Adj[i][j]=G.A[E.vclass[i]][E.vclass[j]]; }
    return E;
}
// count s-cliques in expanded graph that DO NOT contain a fixed vertex set R (an r-clique).
// "contain R" = the s-clique's vertex set is a superset of R.
static long long countSCliquesAvoidingSet(const Expanded& E,int s,const std::vector<int>& R){
    // R must itself be a clique; we count s-subsets that are cliques and NOT superset of R.
    std::vector<char> inR(E.N,0); for(int v:R)inR[v]=1;
    long long cnt=0; std::vector<int> stk;
    std::function<void(int)> rec=[&](int start){
        if((int)stk.size()==s){
            // is R subset of stk?
            bool sup=true; for(int v:R){bool f=false; for(int u:stk) if(u==v){f=true;break;} if(!f){sup=false;break;}}
            if(!sup)++cnt; return;
        }
        if(E.N-start < s-(int)stk.size())return;
        for(int v=start;v<E.N;++v){bool ok=true; for(int u:stk) if(!E.Adj[u][v]){ok=false;break;}
            if(ok){stk.push_back(v);rec(v+1);stk.pop_back();}}
    };
    rec(0); return cnt;
}
static long long countSCliquesTotal(const Expanded& E,int s){
    long long cnt=0; std::vector<int> stk;
    std::function<void(int)> rec=[&](int start){
        if((int)stk.size()==s){++cnt;return;}
        if(E.N-start < s-(int)stk.size())return;
        for(int v=start;v<E.N;++v){bool ok=true; for(int u:stk) if(!E.Adj[u][v]){ok=false;break;}
            if(ok){stk.push_back(v);rec(v+1);stk.pop_back();}}
    };
    rec(0); return cnt;
}

// Pick a concrete representative r-clique of vertices for a pattern P (mult vector).
// For class c with P[c]=m, take the FIRST m vertices of that class.
static std::vector<int> patternToVertexSet(const Expanded& E,const std::vector<int>& P){
    std::vector<int> R;
    for(int c=0;c<(int)P.size();++c) for(int t=0;t<P[c];++t) R.push_back(E.classVerts[c][t]);
    return R;
}

int main(){
    initPascal();
    std::mt19937 rng(424242u);
    printf("==== PEEL SUITABILITY: delete one r-clique pattern, recount s-cliques ====\n");
    printf("Compare (SCT leaves with P dead-boxed) vs expanded-oracle 'avoid superset of R'.\n\n");

    long long cases=0, fails=0, shown=0;
    // also separately track: does SCT total before peel match oracle total?
    for(int trial=0; trial<4000; ++trial){
        int C = 3 + (int)(rng()%5);            // 3..7
        ClassG G; G.C=C; G.w.resize(C);
        for(int c=0;c<C;++c) G.w[c]=1+(int)(rng()%3); // 1..3
        G.A.assign(C,std::vector<char>(C,0));
        double pEdge=(double)(rng()%101)/100.0;
        for(int i=0;i<C;++i)for(int j=i+1;j<C;++j){
            char e=((double)(rng()%1000)/1000.0<pEdge)?1:0; G.A[i][j]=G.A[j][i]=e; }
        int r = 1 + (int)(rng()%2);            // r in 1..2
        int s = r + 1 + (int)(rng()%3);        // s in r+1..r+3

        Expanded E = expand(G);
        long long totS = countSCliquesTotal(E,s);
        // sanity: SCT (no forbidden) equals total s-cliques
        auto leaves = buildClassSCT(G,s);
        double sct0 = sctTotalForbidden(leaves);
        if(std::abs(sct0-(double)totS)>0.5){
            ++fails; if(shown<10){++shown;
                printf("PRE-PEEL MISMATCH trial=%d C=%d s=%d sct=%.0f oracle=%lld\n",trial,C,s,sct0,totS);}
            continue;
        }
        if(totS==0) continue;

        // enumerate all r-clique patterns P (mult vec, sum=r, clique support) and peel each.
        std::vector<int> Pm(C,0); std::vector<int> chosen;
        std::function<void(int,int)> rec=[&](int c,int rem){
            if(rem==0){
                // build vertex set R for this pattern, must be realizable (each m_c <= w_c)
                bool realizable=true; for(int cc=0;cc<C;++cc) if(Pm[cc]>G.w[cc])realizable=false;
                if(realizable){
                    Vec Pv=ccpath::zeros_vec(C); for(int cc=0;cc<C;++cc)Pv[cc]=(int16_t)Pm[cc];
                    // dead-box P into every leaf, recount
                    double peeled=0.0;
                    for(auto L:leaves){ // copy
                        ccpath::insert_antichain(L.forbidden, Pv);
                        Vec b=ccpath::zeros_vec(L.m());
                        peeled+=ccpath::support_count(L,b,nCr_fn);
                    }
                    // oracle: number of s-cliques NOT containing R, summed over ALL realizations
                    // of P? No — dead-boxing pattern P removes s-cliques whose multiplicity
                    // vector y >= P, i.e. that USE at least P[c] vertices of class c FOR EVERY c.
                    // The honest oracle: count weighted s-cliques whose mult vector y does NOT
                    // dominate P. We compute that directly from the expanded graph by counting
                    // s-cliques and checking the per-class usage.
                    long long oracleKept = 0;
                    {
                        std::vector<int> stk;
                        std::function<void(int)> r2=[&](int start){
                            if((int)stk.size()==s){
                                // per-class usage
                                std::vector<int> use(C,0); for(int v:stk)use[E.vclass[v]]++;
                                bool dom=true; for(int cc=0;cc<C;++cc) if(use[cc]<Pm[cc]){dom=false;break;}
                                if(!dom)++oracleKept; return;
                            }
                            if(E.N-start < s-(int)stk.size())return;
                            for(int v=start;v<E.N;++v){bool ok=true;
                                for(int u:stk) if(!E.Adj[u][v]){ok=false;break;}
                                if(ok){stk.push_back(v);r2(v+1);stk.pop_back();}}
                        };
                        r2(0);
                    }
                    ++cases;
                    if(std::abs(peeled-(double)oracleKept)>0.5){
                        ++fails;
                        if(shown<20){++shown;
                            printf("PEEL-FAIL trial=%d C=%d s=%d r=%d  P=[",trial,C,s,r);
                            for(int cc=0;cc<C;++cc)printf("%d ",Pm[cc]);
                            printf("] peeled=%.1f oracleKept=%lld  totS=%lld\n",peeled,oracleKept,totS);
                        }
                    }
                }
                return;
            }
            if(c>=C)return;
            rec(c+1,rem);
            bool compat=true; for(int x:chosen) if(!G.adj(c,x)){compat=false;break;}
            if(compat){int hi=std::min(rem,G.w[c]); chosen.push_back(c);
                for(int j=1;j<=hi;++j){Pm[c]=j;rec(c+1,rem-j);} Pm[c]=0; chosen.pop_back();}
        };
        rec(0,r);
    }
    printf("\n[peel] pattern-delete cases=%lld  FAIL=%lld  => %s\n",
           cases, fails, fails==0?"PASS (dead-box removes exactly the dominating s-cliques)":"FAIL");
    return fails==0?0:1;
}
