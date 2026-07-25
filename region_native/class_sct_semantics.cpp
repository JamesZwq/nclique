// class_sct_semantics.cpp — clarify what a per-leaf dead-box of pattern P means
// for nucleus peeling, vs. peeling a SPECIFIC r-clique instance (vertex set).
//
// Two distinct quantities for an s-clique pattern Y (usage vector, sum=s) and an
// r-clique pattern P (usage vector, sum=r, P<=Y componentwise on a clique):
//
//   (A) #(s-clique instances with usage Y that CONTAIN a fixed instance R of P)
//       = prod_c C(Y_c - P_c, ...)? Actually for fixed R (specific P_c verts of
//       class c chosen), an s-instance with usage Y contains R iff its chosen
//       Y_c-subset of class c is a superset of R's P_c verts.
//       count over Y-instances containing R = prod_c C(w_c - P_c, Y_c - P_c).
//
//   (B) support_count semantics: dead-box pattern P removes Y-instances with
//       Y>=P entirely (all prod_c C(w_c,Y_c) of them). It does NOT subtract a
//       per-instance contribution.
//
// This program quantifies the difference and states which peel model the SCT
// dead-box implements, with concrete numbers.

#include "../src/NucleusDecomposition/CCPathCore.h"
#include <cstdio>
#include <vector>
using ccpath::CCPath; using ccpath::Vec;

static double PASCAL[129][129];
static void initPascal(){for(int n=0;n<=128;++n){PASCAL[n][0]=1.0;
  for(int k=1;k<=n;++k)PASCAL[n][k]=PASCAL[n-1][k-1]+(k<=n-1?PASCAL[n-1][k]:0.0);}}
static double nCr_fn(int n,int k){if(k<0||n<0||k>n)return 0.0;return PASCAL[n][k];}

int main(){
  initPascal();
  printf("==== SEMANTICS OF A PER-LEAF DEAD-BOX (pattern P) ====\n\n");

  // Single clique leaf: classes {0,1} weight {w0,w1}, fully adjacent, s-clique.
  // Build the trivial 1-leaf SCT by hand: n=(w0,w1), ell=0, u=n, T=s.
  auto trial=[&](int w0,int w1,int s,int p0,int p1){
    int C=2;
    CCPath L; L.h=Vec(C,0); L.n=Vec{(int16_t)w0,(int16_t)w1};
    L.ell=Vec(C,0); L.u=L.n; L.T=s;
    Vec b=ccpath::zeros_vec(C);
    double before=ccpath::support_count(L,b,nCr_fn);
    Vec P=Vec{(int16_t)p0,(int16_t)p1};
    CCPath L2=L; ccpath::insert_antichain(L2.forbidden,P);
    double after=ccpath::support_count(L2,b,nCr_fn);
    double removed=before-after;

    // model (B): instances with usage Y>=P removed entirely.
    // sum over Y (Y0+Y1=s, Y0>=p0,Y1>=p1, Y<=w) of C(w0,Y0)C(w1,Y1)
    double modelB=0.0;
    for(int Y0=p0;Y0<=w0;++Y0){int Y1=s-Y0; if(Y1<p1||Y1>w1)continue;
      modelB+=nCr_fn(w0,Y0)*nCr_fn(w1,Y1);}

    // model (A): peel ONE specific instance R of pattern P; count (s-instance, R)
    // incidences where R subset of s-instance. For each s-instance usage Y>=P,
    // the number of its sub-instances equal to R is: it contains R iff Y's chosen
    // verts include R's. # s-instances containing the FIXED R = prod C(w_c-p_c, Y_c-p_c).
    // Summed over all R-instances would double count; classic nucleus peels ONE R
    // at a time. Here we report, for a SINGLE fixed R, how many s-instances contain it.
    double modelA_oneR=0.0;
    for(int Y0=p0;Y0<=w0;++Y0){int Y1=s-Y0; if(Y1<p1||Y1>w1)continue;
      modelA_oneR+=nCr_fn(w0-p0,Y0-p0)*nCr_fn(w1-p1,Y1-p1);}

    printf("w=(%d,%d) s=%d P=(%d,%d): support before=%.0f after=%.0f removed=%.0f | "
           "modelB(all Y>=P insts)=%.0f  modelA(s-insts containing ONE fixed R)=%.0f\n",
           w0,w1,s,p0,p1,before,after,removed,modelB,modelA_oneR);
  };

  trial(3,3,3, 1,1);
  trial(3,3,4, 1,1);
  trial(4,2,3, 2,0);
  trial(4,4,4, 2,1);
  trial(5,5,5, 1,1);
  trial(3,3,2, 1,1);

  printf("\nInterpretation:\n");
  printf(" * The dead-box matches model (B): it removes EVERY s-clique instance\n");
  printf("   whose usage vector dominates P, all at once. removed==modelB.\n");
  printf(" * Model (A) (peel ONE specific r-clique vertex set R) gives a DIFFERENT,\n");
  printf("   smaller number. A pattern dead-box is NOT a single-instance peel.\n");
  return 0;
}
