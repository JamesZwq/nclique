// class_sct_stress.cpp — boundary + worst-case leaf-count stress for the builder.
#include "ClassSCT.h"
#include <cstdio>
#include <chrono>
#include <vector>
#include <functional>
#include <cmath>

static double PASCAL[2050][2050];
static bool pinit=false;
static void initPascal(){ // big table so large weights/k don't OOB
  for(int n=0;n<=2048;++n){PASCAL[n][0]=1.0;
    for(int k=1;k<=n;++k)PASCAL[n][k]=PASCAL[n-1][k-1]+(k<=n-1?PASCAL[n-1][k]:0.0);}
  pinit=true;
}
static double nCr_fn(int n,int k){if(k<0||n<0||k>n||n>2048)return 0.0;return PASCAL[n][k];}

static double sctTot(const std::vector<CCPath>&L){double t=0;
  for(auto&l:L){Vec b=ccpath::zeros_vec(l.m());t+=ccpath::support_count(l,b,nCr_fn);}return t;}

// expanded-graph oracle for big-weight boundary correctness
static long long oracle(const ClassG&G,int k){
  std::vector<int> vc; for(int c=0;c<G.C;++c)for(int t=0;t<G.w[c];++t)vc.push_back(c);
  int N=(int)vc.size(); if(k==0)return 1; if(k>N)return 0;
  std::vector<std::vector<char>> A(N,std::vector<char>(N,0));
  for(int i=0;i<N;++i)for(int j=0;j<N;++j){if(i==j)continue;
    A[i][j]= (vc[i]==vc[j])?1:G.A[vc[i]][vc[j]];}
  long long c=0; std::vector<int> stk;
  std::function<void(int)> rec=[&](int s){if((int)stk.size()==k){++c;return;}
    if(N-s<k-(int)stk.size())return;
    for(int v=s;v<N;++v){bool ok=true;for(int u:stk)if(!A[u][v]){ok=false;break;}
      if(ok){stk.push_back(v);rec(v+1);stk.pop_back();}}};
  rec(0); return c;
}

int main(){
  initPascal();
  printf("==== BOUNDARY: large weight / large k (>64) correctness ====\n");
  // 1 class w=70 k=35: C(70,35) is astronomically large but representable as double.
  { ClassG G; G.C=1; G.w={70}; G.A.assign(1,std::vector<char>(1,0));
    auto L=buildClassSCT(G,35); double sc=sctTot(L);
    printf("  1cls w=70 k=35: sct=%.6g  C(70,35)=%.6g  %s\n",
           sc, nCr_fn(70,35), std::abs(sc-nCr_fn(70,35))<1e-3*nCr_fn(70,35)?"OK":"FAIL");}
  // small enough to oracle-check with big-ish weights:
  { ClassG G; G.C=2; G.w={9,9}; G.A.assign(2,std::vector<char>(2,1));G.A[0][0]=G.A[1][1]=0;
    int k=10; auto L=buildClassSCT(G,k); double sc=sctTot(L); long long ex=oracle(G,k);
    printf("  K2 w=9 k=10: sct=%.0f oracle=%lld %s\n",sc,ex,std::abs(sc-(double)ex)<0.5?"OK":"FAIL");}
  { ClassG G; G.C=3; G.w={8,8,8}; G.A.assign(3,std::vector<char>(3,1));for(int i=0;i<3;++i)G.A[i][i]=0;
    int k=12; auto L=buildClassSCT(G,k); double sc=sctTot(L); long long ex=oracle(G,k);
    printf("  K3 w=8 k=12: sct=%.0f oracle=%lld %s\n",sc,ex,std::abs(sc-(double)ex)<0.5?"OK":"FAIL");}

  printf("\n==== WORST-CASE LEAF COUNT (Moon-Moser & complement-of-matching) ====\n");
  // Moon-Moser graph: complete multipartite K_{3,3,...,3} (t groups of 3 NON-adjacent).
  // = complement of disjoint triangles. Max # maximal cliques = 3^t. Each maximal
  // clique picks one vertex per group. With weights, the SCT should still be small
  // per maximal clique. We test t groups, each group = 3 classes mutually non-adj,
  // all cross-group pairs adjacent.
  printf("  -- Moon-Moser K_{3,3,..} (t groups of 3), w=1, k=t --\n");
  for(int t=1;t<=6;++t){
    int C=3*t; ClassG G; G.C=C; G.w.assign(C,1); G.A.assign(C,std::vector<char>(C,1));
    for(int i=0;i<C;++i)G.A[i][i]=0;
    for(int gp=0;gp<t;++gp)for(int a=0;a<3;++a)for(int b=a+1;b<3;++b){
      int x=gp*3+a,y=gp*3+b; G.A[x][y]=G.A[y][x]=0;}
    int k=t;
    auto t0=std::chrono::high_resolution_clock::now();
    auto L=buildClassSCT(G,k);
    auto t1=std::chrono::high_resolution_clock::now();
    double ms=std::chrono::duration<double,std::milli>(t1-t0).count();
    long long ex = (C<=24)? oracle(G,k) : -1;
    double sc=sctTot(L);
    printf("  t=%d (C=%2d) k=%d: leaves=%-7zu  3^t=%-6.0f  sct=%.0f oracle=%lld  %.2fms %s\n",
           t,C,k,L.size(),std::pow(3.0,t),sc,ex,ms,
           (ex<0||std::abs(sc-(double)ex)<0.5)?"":"<<FAIL");
  }
  printf("\n  -- complement of perfect matching (cocktail party) w=1, sweep k --\n");
  for(int C=10;C<=24;C+=2){
    ClassG G; G.C=C; G.w.assign(C,1); G.A.assign(C,std::vector<char>(C,1));
    for(int i=0;i<C;++i)G.A[i][i]=0;
    for(int i=0;i+1<C;i+=2)G.A[i][i+1]=G.A[i+1][i]=0;
    int k=C/2; // max clique size = C/2
    auto t0=std::chrono::high_resolution_clock::now();
    auto L=buildClassSCT(G,k);
    auto t1=std::chrono::high_resolution_clock::now();
    double ms=std::chrono::duration<double,std::milli>(t1-t0).count();
    printf("  C=%2d k=%2d: leaves=%-8zu  2^(C/2)=%-8.0f  %.2fms\n",
           C,k,L.size(),std::pow(2.0,C/2),ms);
  }

  printf("\n  -- WEIGHTED Moon-Moser (w=3) to see weight x structure interaction --\n");
  for(int t=1;t<=5;++t){
    int C=3*t; ClassG G; G.C=C; G.w.assign(C,3); G.A.assign(C,std::vector<char>(C,1));
    for(int i=0;i<C;++i)G.A[i][i]=0;
    for(int gp=0;gp<t;++gp)for(int a=0;a<3;++a)for(int b=a+1;b<3;++b){
      int x=gp*3+a,y=gp*3+b; G.A[x][y]=G.A[y][x]=0;}
    int k=std::min(C,6);
    auto t0=std::chrono::high_resolution_clock::now();
    auto L=buildClassSCT(G,k);
    auto t1=std::chrono::high_resolution_clock::now();
    double ms=std::chrono::duration<double,std::milli>(t1-t0).count();
    printf("  t=%d (C=%2d) w=3 k=%d: leaves=%-8zu  %.2fms\n",t,C,k,L.size(),ms);
  }
  return 0;
}
