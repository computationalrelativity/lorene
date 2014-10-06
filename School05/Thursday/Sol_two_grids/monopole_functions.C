#include <math.h>

#include "scalar.h"
#include "monopole.h"
#include "proto.h"
#include "graphique.h"

void Monopole::init_big_W() {

  double rlim = mp_W.val_r (0,1,0,0) ;
 
  // Polynomial near r=0
  Scalar auxi_un (mp_W) ;
  auxi_un = 1-0.5*radius_on_W*radius_on_W/rlim/rlim ;
  auxi_un.annule (1, nz_W-1) ;

  // Exponential at infinity :
  Scalar auxi_deux (mp_W) ;
  auxi_deux = 0.5*exp(-(radius_on_W-rlim)/rlim) ;
  auxi_deux.set_outer_boundary(nz_W-1,0) ;
  auxi_deux.annule_domain(0) ;

  big_W = auxi_un+auxi_deux ;
  big_W.std_spectral_base() ;
  
  big_W_on_H.set_etat_qcq() ;
  big_W_on_H.import(big_W) ;
  big_W_on_H.std_spectral_base() ;
}

void Monopole::init_big_H() {
  double rlim = mp_H.val_r (0,1,0,0) ;
  Scalar auxi_un (radius_on_H/rlim) ;
  auxi_un.annule(1,nz_H-1) ;

  Scalar auxi_deux (mp_H) ;
  auxi_deux = 1 ;
  auxi_deux.annule_domain(0) ;

  big_H = auxi_un + auxi_deux ;
  big_H.std_spectral_base_odd() ;
  
  big_H_on_W.import(big_H) ;
  big_H_on_W.std_spectral_base_odd() ;
}

void Monopole::do_small_w () {

  Scalar auxi_un (mp_W) ;
  auxi_un = (big_W-1)/radius_on_W ;
  // Put zero at the origin
  auxi_un.set_inner_boundary(0,0) ;
  auxi_un.annule(1, nz_W-1) ;

  Scalar auxi_deux (big_W) ;
  auxi_deux.annule_domain(0) ;
  small_w = auxi_un+auxi_deux ;

  // Basis :
  small_w.std_spectral_base_odd() ;
 
  Scalar auxi_un_on_H (mp_H) ;
  auxi_un_on_H = (big_W_on_H-1)/radius_on_H ;
  // Put zero at the origin
  auxi_un_on_H.set_inner_boundary(0,0) ;
  auxi_un_on_H.annule(1, nz_H-1) ;

  Scalar auxi_deux_on_H (big_W_on_H) ;
  auxi_deux_on_H.annule_domain(0) ;
  small_w_on_H = auxi_un_on_H+auxi_deux_on_H ;

  // Basis :
  small_w_on_H.std_spectral_base_odd() ;
}

void Monopole::do_big_W () {

  Scalar auxi_un (mp_W) ;
  auxi_un = small_w*radius_on_W ;
  // Put zero at origin
  auxi_un.set_inner_boundary(0,0) ;
  auxi_un = auxi_un+1 ;
  auxi_un.annule(1, nz_W-1) ;

  Scalar auxi_deux (small_w) ;
  auxi_deux.annule_domain(0) ;

  big_W = auxi_un+auxi_deux ;

  // Basis :
  big_W.std_spectral_base() ;
  big_W_on_H.import(big_W) ;
  big_W_on_H.std_spectral_base() ;
}

void Monopole::do_small_h() {  

  Scalar auxi (big_H-1) ;
  auxi.annule(0, nz_H-2) ;

  small_h = big_H ;
  small_h.annule_domain(nz_H-1) ;
  small_h += auxi ;
  small_h.std_spectral_base_odd() ;

  Scalar auxi_on_W (big_H_on_W-1) ;
  auxi_on_W.annule(0, nz_W-2) ;

  small_h_on_W = big_H_on_W ;
  small_h_on_W.annule_domain(nz_W-1) ;
  small_h_on_W += auxi_on_W ;
  small_h_on_W.std_spectral_base_odd() ;
}

void Monopole::do_big_H() {
  Scalar auxi (small_h+1) ;
  auxi.annule(0, nz_H-2) ;

  big_H = small_h;
  big_H.annule_domain(nz_H-1) ;
  big_H += auxi ;
  big_H.std_spectral_base_odd() ;

  big_H_on_W.import(big_H) ; 
  big_H_on_W.std_spectral_base_odd() ;
}

Scalar Monopole::compute_source_W() const {

  // Near r = 0
  Scalar source_w (mp_W) ;
  source_w = small_w*small_w*small_w + 3*small_w*small_w/radius_on_W +
    (1+radius_on_W*small_w)*big_H_on_W*big_H_on_W/radius_on_W ;
  source_w.set_inner_boundary(0,0) ;
  source_w.annule(1, nz_W-1) ;

  // Near infinity
  Scalar source_zec (mp_W) ;
  source_zec = big_W*small_h_on_W*(small_h_on_W+2) ;
  source_zec.std_spectral_base() ;
  source_zec.inc_dzpuis(2) ;
  Scalar auxi(big_W*(big_W*big_W-1)) ;
  auxi.set_dzpuis(2) ;
  source_zec += auxi ;
  source_zec += 2*big_W.dsdr()/radius_on_W ;
  source_zec.set_outer_boundary(nz_W-1, 0) ;
  source_zec.annule(0, nz_W-2) ;
  
  // In the shells :
  Scalar source_W (mp_W) ;
  source_W = big_W*((big_W*big_W-1)/radius_on_W/radius_on_W + 
		    big_H_on_W*big_H_on_W-1.) ;
  source_W.annule_domain(nz_W-1) ;
  source_W += 2*big_W.dsdr()/radius_on_W ;
  source_W.annule_domain(nz_W-1) ;
  source_W.annule_domain(0) ;
  source_W.std_spectral_base() ;

  Scalar source (source_W+source_w+source_zec) ;
  source.set_spectral_base(small_w.get_spectral_va().base) ;

  return source ;
}

Scalar Monopole::compute_source_H() const {

  // near r = 0 ;
  Scalar source_noyau (mp_H) ;
  source_noyau = 2*big_H*
    (small_w_on_H*small_w_on_H+2*small_w_on_H/radius_on_H) + 
    beta*beta/2.*big_H*(big_H*big_H-1) ;
  source_noyau.annule(1,nz_H-1) ;
  source_noyau.set_inner_boundary(0,0) ;
  
  // ZEC :
  Scalar source_zec (mp_H) ;
  source_zec = 2*big_W_on_H*big_W_on_H*(small_h+1) ;
  source_zec.set_dzpuis(2) ;
  Scalar auxi (mp_H) ;
  auxi = beta*beta/2.*small_h*small_h*(3.+small_h) ;
  auxi.std_spectral_base() ;
  auxi.inc_dzpuis(2) ;
  source_zec += auxi ;
  source_zec.std_spectral_base() ;
  source_zec.annule(0, nz_H-2) ;
  
  // Coquilles :
  Scalar source_H (mp_H) ;
  source_H = 2*big_W_on_H*big_W_on_H*big_H/radius_on_H/radius_on_H + 
    beta*beta/2.*big_H*(big_H*big_H-3) ;
  source_H.annule_domain(0) ;
  source_H.annule_domain(nz_H-1) ;

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

Scalar Monopole::density_on_W() const {

	double par_b = give_b() ;

	// Part W' :
	Scalar part_w_prime (big_W.dsdr()*big_W.dsdr()) ;
	part_w_prime.div_r() ;
	part_w_prime.div_r() ;
	part_w_prime.annule_domain(nz_W-1) ;
	
	Scalar part_zec (big_W.dsdr()) ;
	part_zec.dec_dzpuis(2) ;
	part_zec.set_dzpuis(2) ;
	part_zec = part_zec*big_W.dsdr() ;
	part_zec.annule (0, nz_W-2) ;
	part_w_prime = part_w_prime + part_zec ;
	
	// Part W2-1
	Scalar part_w2 ((big_W*big_W-1)*(big_W*big_W-1)/2)  ;
	for (int i=0 ; i<4 ; i++)
	    part_w2.div_r() ;
	part_w2.annule_domain(nz_W-1) ;
	
	part_zec = (big_W*big_W-1)*(big_W*big_W-1)/2 ;
	part_zec.set_dzpuis(4) ;
	part_zec.annule (0, nz_W-2) ;
	part_w2 = part_w2 + part_zec ;
	
	Scalar res (part_w_prime+part_w2) ;
	res.std_spectral_base() ;
	res = res/4./M_PI ;
	
	return res ;	
}

Scalar Monopole::density_on_H() const {

	// Part en Hprime :
	Scalar part_h_prime (big_H.dsdr()*big_H.dsdr()/2) ;
	part_h_prime.annule_domain(nz_H-1) ;
	
	Scalar part_zec (big_H.dsdr()*big_H.dsdr()/2) ;
	part_zec.annule (0, nz_H-2) ;
	part_h_prime = part_h_prime + part_zec ;
	
	
	//Part en h2
	Scalar part_h2 ((big_H*big_H-1)*(big_H*big_H-1)*beta*beta/8) ;
	part_h2.inc_dzpuis(4) ; 
	
	// Part WH :
	Scalar part_wh (big_W_on_H*big_W_on_H*big_H*big_H) ;
	part_wh.div_r() ;
	part_wh.div_r() ;
	part_wh.annule_domain(nz_H-1) ;
	
	part_zec = big_W_on_H*big_W_on_H*big_H*big_H ;
	part_zec.set_dzpuis(2) ;
	part_zec.inc_dzpuis(2) ;
	part_zec.annule (0, nz_H-2) ;
	part_wh = part_wh + part_zec ;
	
	
	Scalar res (part_h_prime+part_h2+part_wh) ;
	res.std_spectral_base() ;
	res = res/4./M_PI ;
	
	return res ;	
}

double Monopole::energy() const {

    Scalar part_W (density_on_W()) ;
    Scalar part_H (density_on_H()) ;
    
    return part_W.integrale() + part_H.integrale() ;
}


  
