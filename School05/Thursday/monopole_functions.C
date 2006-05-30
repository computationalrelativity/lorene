#include <math.h>

#include "scalar.h"
#include "monopole.h"
#include "proto.h"
#include "graphique.h"

void Monopole::init_big_W() {
  double rlim = mp.val_r (0,1,0,0) ;
 
  big_W.allocate_all() ;
  //for (int i=0 ; i<nz-1 ; i++)
       big_W.set_domain(0) =  (1-0.5*radius*radius/rlim/rlim).domain(0) ;
  for (int i=1 ; i<nz ; i++) big_W.set_domain(i) = (0.5*exp(-(radius-rlim)/rlim)).domain(i) ;
  
  big_W.std_spectral_base() ;
}

void Monopole::init_big_H() {
 cout << "Monopole::init_big_H is not implemented" << endl ;
}

void Monopole::do_small_w () {
  small_w = big_W ;
  
  Scalar auxi ((big_W-1)/radius) ;
  auxi.set_inner_boundary(0,0) ;
  small_w.set_domain(0) = auxi.domain(0) ;

  // Basis :
  small_w.std_spectral_base_odd() ;
}

void Monopole::do_big_W () {
 cout << "Monopole::do_big_W is not implemented" << endl ;
}

void Monopole::do_small_h() {  
 cout << "Monopole::do_small_h is not implemented" << endl ;
}

void Monopole::do_big_H() {
 cout << "Monopole::do_big_H is not implemented" << endl ;
}

Scalar Monopole::compute_source_W() const {
    cout << "Monopole::compute_source_w is not implemented" << endl ;
    return Scalar(mp) ; /// To avoid compilation error
}

Scalar Monopole::compute_source_H() const {

 Scalar source (mp) ;
  source.allocate_all() ;
  
  // near r = 0 ;
  Scalar source_noyau (mp) ;
  source_noyau = 2*big_H*
    (small_w*small_w+2*small_w/radius) + beta*beta/2.*big_H*(big_H*big_H-1) ;
  source_noyau.set_inner_boundary(0,0) ;
  source.set_domain(0) = source_noyau.domain(0) ;
  
  // CED :
  Scalar source_zec (mp) ;
  source_zec = 2*big_W*big_W*(small_h+1) ;
  source_zec.set_dzpuis(2) ;
  Scalar auxi (mp) ;
  auxi = beta*beta/2.*small_h*small_h*(3.+small_h) ;
  auxi.std_spectral_base() ;
  auxi.inc_dzpuis(2) ;
  source_zec += auxi ;
  source.set_domain(nz-1) = source_zec.domain(nz-1) ;
  source.set_dzpuis(source_zec.get_dzpuis()) ;
  
  // Coquilles :
  Scalar source_H (mp) ;
  source_H = 2*big_W*big_W*big_H/radius/radius + beta*beta/2.*big_H*(big_H*big_H-3) ;
  for (int i=1 ; i<nz-1 ; i++)
      source.set_domain(i) = source_H.domain(i) ;

  source.std_spectral_base_odd() ;
  return source ;
}
  
double Monopole::give_a() const {
  cout << "Monopole::give_a is not implemented" << endl ;
  return 1 ; // To avoid compilation errors
}

double Monopole::give_b() const {
  cout << "Monopole::give_b is not implemented" << endl ;
  return 1 ; /// To avoid compilatin error
}

