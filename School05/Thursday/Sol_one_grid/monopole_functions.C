#include <math.h>

#include "scalar.h"
#include "monopole.h"
#include "proto.h"
#include "graphique.h"

void Monopole::init_big_W() {

  double rlim = mp.val_r (0,1,0,0) ;
 
  // Polynomial near r=0
  Scalar auxi_un (mp) ;
  auxi_un = 1-0.5*radius*radius/rlim/rlim ;
  auxi_un.annule (1, nz-1) ;

  // Exponential at infinity :
  Scalar auxi_deux (mp) ;
  auxi_deux = 0.5*exp(-(radius-rlim)/rlim) ;
  auxi_deux.set_outer_boundary(nz-1,0) ;
  auxi_deux.annule_domain(0) ;

  big_W = auxi_un+auxi_deux ;
  big_W.std_spectral_base() ;
}

void Monopole::init_big_H() {

  double rlim = mp.val_r (0,1,0,0) ;
  Scalar auxi_un (radius/rlim) ;
  auxi_un.annule(1,nz-1) ;

  Scalar auxi_deux (mp) ;
  auxi_deux = 1 ;
  auxi_deux.annule_domain(0) ;

  big_H = auxi_un + auxi_deux ;
  big_H.std_spectral_base_odd() ;
}

void Monopole::do_small_w () {

  Scalar auxi_un (mp) ;
  auxi_un = (big_W-1)/radius ;
  // Put zero at the origin
  auxi_un.set_inner_boundary(0,0) ;
  auxi_un.annule(1, nz-1) ;

  Scalar auxi_deux (big_W) ;
  auxi_deux.annule_domain(0) ;
  small_w = auxi_un+auxi_deux ;

  // Basis :
  small_w.std_spectral_base_odd() ;
}

void Monopole::do_big_W () {

  Scalar auxi_un (mp) ;
  auxi_un = small_w*radius ;
  // Put zero at origin
  auxi_un.set_inner_boundary(0,0) ;
  auxi_un = auxi_un+1 ;
  auxi_un.annule(1, nz-1) ;

  Scalar auxi_deux (small_w) ;
  auxi_deux.annule_domain(0) ;

  big_W = auxi_un+auxi_deux ;

  // Basis :
  big_W.std_spectral_base() ;
}

void Monopole::do_small_h() {  
  Scalar auxi (big_H-1) ;
  auxi.annule(0, nz-2) ;

  small_h = big_H ;
  small_h.annule_domain(nz-1) ;
  small_h += auxi ;
  small_h.std_spectral_base_odd() ;
}

void Monopole::do_big_H() {
  Scalar auxi (small_h+1) ;
  auxi.annule(0, nz-2) ;

  big_H = small_h;
  big_H.annule_domain(nz-1) ;
  big_H += auxi ;
  big_H.std_spectral_base_odd() ;
}

Scalar Monopole::compute_source_W() const {

  // Near r = 0
  Scalar source_w (mp) ;
  source_w = small_w*small_w*small_w + 3*small_w*small_w/radius +
    (1+radius*small_w)*big_H*big_H/radius ;
  source_w.set_inner_boundary(0,0) ;
  source_w.annule(1, nz-1) ;

  // Near infinity
  Scalar source_zec (mp) ;
  source_zec = big_W*small_h*(small_h+2) ;
  source_zec.std_spectral_base() ;
  source_zec.inc_dzpuis(2) ;
  Scalar auxi(big_W*(big_W*big_W-1)) ;
  auxi.set_dzpuis(2) ;
  source_zec += auxi ;
  source_zec += 2*big_W.dsdr()/radius ;
  source_zec.set_outer_boundary(nz-1, 0) ;
  source_zec.annule(0, nz-2) ;
  
  // In the shells :
  Scalar source_W (mp) ;
  source_W = big_W*((big_W*big_W-1)/radius/radius + 
		    big_H*big_H-1.) ;
  source_W.annule_domain(nz-1) ;
  source_W += 2*big_W.dsdr()/radius ;
  source_W.annule_domain(nz-1) ;
  source_W.annule_domain(0) ;
  source_W.std_spectral_base() ;

  Scalar source (source_W+source_w+source_zec) ;
  source.set_spectral_base(small_w.get_spectral_va().base) ;
  
  return source ;
}

Scalar Monopole::compute_source_H() const {

  // near r = 0 ;
  Scalar source_noyau (mp) ;
  source_noyau = 2*big_H*
    (small_w*small_w+2*small_w/radius) + beta*beta/2.*big_H*(big_H*big_H-1) ;
  source_noyau.annule(1,nz-1) ;
  source_noyau.set_inner_boundary(0,0) ;
  
  // ZEC :
  Scalar source_zec (mp) ;
  source_zec = 2*big_W*big_W*(small_h+1) ;
  source_zec.set_dzpuis(2) ;
  Scalar auxi (mp) ;
  auxi = beta*beta/2.*small_h*small_h*(3.+small_h) ;
  auxi.std_spectral_base() ;
  auxi.inc_dzpuis(2) ;
  source_zec += auxi ;
  source_zec.std_spectral_base() ;
  source_zec.annule(0, nz-2) ;
  
  // Coquilles :
  Scalar source_H (mp) ;
  source_H = 2*big_W*big_W*big_H/radius/radius + beta*beta/2.*big_H*(big_H*big_H-3) ;
  source_H.annule_domain(0) ;
  source_H.annule_domain(nz-1) ;

  Scalar source (source_noyau+source_zec+source_H) ;
  source.set_spectral_base(big_H.get_spectral_va().base) ;
  return source ;
}
  
double Monopole::give_a() const {
  double res = big_H.dsdr().val_grid_point(0,0,0,0) ;
  return res ;
}

double Monopole::give_b() const {
  Scalar auxi (big_W.dsdr()) ;
  auxi.dec_dzpuis(2) ;
  double res = auxi.dsdr().val_grid_point(0,0,0,0)/2. ;
  return -res ;
}
