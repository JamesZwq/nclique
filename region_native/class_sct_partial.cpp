// class_sct_partial.cpp — the decisive peel-suitability question:
//
// When peeling proceeds INSTANCE by INSTANCE (the real nucleus algorithm peels
// the single min-support r-clique, not a whole symmetric orbit at once), can the
// pattern-level dead-box represent the state "some instances of pattern P are
// dead, others alive"?
//
// Concretely: take classes {0,1} weight (w0,w1), s-clique leaf. r=2, pattern
// P=(1,1) has w0*w1 distinct instances (one per cross pair). Peel them ONE at a
// time. After peeling exactly ONE instance R0=(vert a of class0, vert b of
// class1), the s-cliques that should lose 1 support are those whose vertex set
// contains {a,b}. We compare:
//   (i) the EXACT alive s-clique count after removing the single R0 (expanded);
//   (ii) what the dead-box machinery can express: it can only forbid the whole
//        PATTERN (a,b) at the class level, i.e. forbid usage>= (1,1), which kills
//        ALL s-cliques touching class0&class1, not just those through {a,b}.
//
// If (i) != (ii) for a single-instance peel, then the class SCT dead-box is a
// PATTERN/ORBIT peel, and single-instance peeling needs additional machinery.

#include "../src/NucleusDecomposition/CCPathCore.h"
#include <cstdio>
#include <functional>
#include <vector>
using ccpath::CCPath; using ccpath::Vec;

static double PASCAL[129][129];
static void initPascal(){for(int n=0;n<=128;++n){PASCAL[n][0]=1.0;
  for(int k=1;k<=n;++k)PASCAL[n][k]=PASCAL[n-1][k-1]+(k<=n-1?PASCAL[n-1][k]:0.0);}}
static double nCr_fn(int n,int k){if(k<0||n<0||k>n)return 0.0;return PASCAL[n][k];}

int main(){
  initPascal();
  printf("==== PARTIAL (instance-level) PEEL probe ====\n\n");

  // expanded simple graph for K(w0)+K(w1) fully joined = K_{w0+w1}? no: two
  // classes fully adjacent => the whole thing is one clique of size w0+w1.
  // s-clique = any s-subset. Peel ONE cross r-clique R0={a,b}.
  auto trial=[&](int w0,int w1,int s){
    int N=w0+w1;
    // vertices 0..w0-1 class0, w0..N-1 class1; complete graph (both classes adj).
    // total s-cliques = C(N,s).
    double total=nCr_fn(N,s);
    // pick R0 = {vertex 0 (class0), vertex w0 (class1)}  (one specific (1,1) inst)
    // exact: s-cliques containing both 0 and w0 = C(N-2, s-2).
    double containR0=nCr_fn(N-2,s-2);
    double aliveExact = total - containR0;  // peel exactly that one instance

    // dead-box of PATTERN (1,1): forbids usage>=(1,1). support after:
    int C=2; CCPath L; L.h=Vec(C,0); L.n=Vec{(int16_t)w0,(int16_t)w1};
    L.ell=Vec(C,0); L.u=L.n; L.T=s;
    Vec b=ccpath::zeros_vec(C);
    CCPath L2=L; ccpath::insert_antichain(L2.forbidden, Vec{(int16_t)1,(int16_t)1});
    double aliveDeadbox=ccpath::support_count(L2,b,nCr_fn);

    printf("w=(%d,%d) s=%d : total=%.0f  alive after ONE-instance peel(exact)=%.0f  "
           "alive after PATTERN(1,1) dead-box=%.0f   %s\n",
           w0,w1,s,total,aliveExact,aliveDeadbox,
           (std::abs(aliveExact-aliveDeadbox)<0.5?"[match]":"[DIFFER -> orbit peel only]"));
  };
  trial(2,2,3);
  trial(3,3,3);
  trial(3,2,3);
  trial(4,3,4);
  trial(2,2,2);

  printf("\nNote: a single class with w>1 also has intra-class r-cliques. The class\n");
  printf("dead-box can only express a PATTERN (multiplicity) threshold, so it peels\n");
  printf("symmetric orbits atomically. This is correct IF the nucleus algorithm is\n");
  printf("defined on ORBITS (all symmetric r-cliques share a core number, peeled\n");
  printf("together). It is NOT a drop-in for classic per-instance peeling.\n");
  return 0;
}
