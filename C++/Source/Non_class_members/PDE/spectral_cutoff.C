
// Header Lorene:
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"
#include "math.h"
#include "metric.h"
#include "param.h"
#include "param_elliptic.h"
#include "vector.h"
#include "scalar.h"
#include "spheroid.h"
#include "diff.h"
#include "proto.h"
#include "unites.h"
#include "tensor.h"
#include "sym_tensor.h"

// Spectral cutoff used in tensor elliptic solvers, and solving for stationary black hole spacetimes


namespace Lorene {
void coupe_l_tous( Sym_tensor& hij,Sym_tensor& aa, Scalar& nn,Scalar& ppsi,Vector& bb, int ntt, int cutoff){

  nn.annule_l(2*(ntt-1) - cutoff, 2*(ntt-1));
  ppsi.annule_l(2*(ntt-1) - cutoff, 2*(ntt-1)); // Warning! only true for SYMMETRY in theta and phi.
  Scalar bb1 = bb.set(1);
  bb1.annule_l(2*(ntt-1) - cutoff, 2*(ntt-1));
  Scalar mmu = bb.mu();
  Scalar etta = bb.eta();
 
  mmu.annule_l(2*(ntt-1) - cutoff, 2*(ntt-1));
  etta.annule_l(2*(ntt-1) - cutoff, 2*(ntt-1));
  
  bb.set_vr_eta_mu(bb1, etta, mmu);
  
  tensor_coupe_l(aa, ntt, cutoff);
  tensor_coupe_l(hij, ntt, cutoff);

  // hij_new.set_auxiliary(hrrBC, tilde_etaBC, mmuAsr, wwBC, xxA, hh -hrrBC);


  return;

}

void tensor_coupe_l( Sym_tensor& ten, int ntt, int cutoff){

  Scalar ten1 = ten.set(1,1); 
  ten1.annule_l(2*(ntt-1) - cutoff, 2*(ntt-1));  // Warning! only true for SYMMETRY in theta and phi.
  int dzp = ten1.get_dzpuis();

  Scalar eta = ten.eta();
  Scalar mmu = ten.mu();
  Scalar xxx = ten.xxx();
  Scalar www = ten.www();
  Scalar smalltrace = ten.set(2,2) + ten.set(3,3);
 
  eta.annule_l(2*(ntt-1) - cutoff, 2*(ntt-1));
  mmu.annule_l(2*(ntt-1) - cutoff, 2*(ntt-1));
  xxx.annule_l(2*(ntt-1) - cutoff, 2*(ntt-1));
  www.annule_l(2*(ntt-1) - cutoff, 2*(ntt-1));
  smalltrace.annule_l(2*(ntt-1) - cutoff, 2*(ntt-1));

  eta.div_r_dzpuis(dzp);
  mmu.div_r_dzpuis(dzp); // set_auxiliary needs quantities rescaled over r


  
  ten.set_auxiliary(ten1, eta, mmu, www, xxx, smalltrace);
  
  return;
}
  
}
