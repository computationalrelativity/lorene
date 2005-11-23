#include "scalar.h"
#include "monopole.h"
#include "param_elliptic.h"
#include "graphique.h"

int Monopole::solve_config(double precis, int itemax, double relax) {

  double erreur ;
  int ite = 1 ;
  
  init_big_W() ;
  init_big_H() ;
  do_small_h() ;
  do_small_w() ;

  Scalar big_W_old(mp) ;
  Scalar big_H_old(mp) ;

  Scalar source_W (compute_source_W()) ;
  Scalar source_H (compute_source_H()) ;

  Param_elliptic param_W (source_W) ;

  Scalar F_W (mp) ;
  F_W = 1 ;
  F_W.std_spectral_base() ;
  F_W.annule(1, nz-1) ;
  param_W.set_variable_F(F_W) ;

  Scalar G_W (mp) ;
  G_W = 1 ;
  G_W.std_spectral_base() ;
  Scalar copie(G_W) ;
  G_W.mult_r() ;
  G_W.annule(1, nz-1) ;
  copie.annule_domain(0) ;
  G_W += copie ;
  param_W.set_variable_G(G_W) ;

  param_W.inc_l_quant(0) ;
  for (int i=1 ; i<nz ; i++)
    param_W.set_helmholtz_minus (i, 1., source_W) ;
  
  Param_elliptic param_H (source_H) ;
  
  Scalar F_H (mp) ;
  F_H = 1 ;
  F_H.std_spectral_base() ;
  F_H.annule(0, nz-2) ;
  param_H.set_variable_F (F_H) ;
  param_H.inc_l_quant(0) ;
  for (int i=1 ; i<nz ; i++)
    param_H.set_helmholtz_minus (i, beta, source_H) ;

  do {
    big_W_old = big_W ;
    big_H_old = big_H ;
    
    do_small_w() ;
    do_small_h() ;
    
    source_W = compute_source_W() ;   
    small_w = source_W.sol_elliptic(param_W) ;
    do_big_W() ;
    big_W = relax*big_W + (1-relax)*big_W_old ;

    source_H = compute_source_H() ;
    small_h = source_H.sol_elliptic(param_H) ;
    do_big_H() ;
    big_H = relax*big_H + (1-relax)*big_H_old ;

    erreur = max(diffrelmax(big_W, big_W_old)) ;
    cout << "Iteration " << ite << " : " << erreur << endl ;
    ite ++ ;
  }
  while ((erreur > precis) && (ite < itemax)) ;
  return ite ;
}
