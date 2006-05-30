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
// Tout ce qui suit est implémenté lors du TP//
void Monopole::init_big_H() {
    double rlim = mp.val_r (0,1,0,0) ;// taille radiale caractéristique; il importe de fixer la grille et le mapping par la suite. 
    big_H.allocate_all() ; // donne de la place physique à W; ce n'est pas fait par défaut, pour gagner de la mémoire...
     big_H.set_domain(0) =  (1+10*radius/rlim).domain(0) ;
  for (int i=1 ; i<nz ; i++) big_H.set_domain(i) = (10*exp(-(radius-rlim)/rlim)).domain(i) ;
  
  big_H.std_spectral_base() ; // correspond aux comportements asymptotiques; les deux parties de la courbe ont un raccordement (au moins) C0 en r=1 du mapping 
}

void Monopole::do_small_w () {
  small_w = big_W ;
  
  Scalar auxi ((big_W-1)/radius) ;
  auxi.set_inner_boundary(0,0) ;
  small_w.set_domain(0) = auxi.domain(0) ;

  // Basis :
  small_w.std_spectral_base_odd() ; // Pourquoi???? Il est juste dit que w est impair au niveau de l'origine... 
 }

void Monopole::do_big_W () {
    big_W = small_w;
   Scalar auxi (radius*small_w -1.) ;
  auxi.set_inner_boundary(0,1.) ;
  big_W.set_domain(0) = auxi.domain(0) ;

  // w et W coincident ailleurs que dans le noyau!
}

void Monopole::do_small_h() {  // un peu une redite de précédemment... 
  small_h = big_H ;
  
  Scalar auxi (big_H -1.) ;
  auxi.set_outer_boundary(2,0) ;
  small_h.set_domain(2) = auxi.domain(2) ;

  // Basis :
  small_h.std_spectral_base() ; // pas pris la odd, mais la complete...
}

void Monopole::do_big_H() {
    big_H=small_h;
    Scalar auxi (small_h+1.);
    auxi.set_outer_boundary (2,1.) ; //! la limite de H en l'infini est 1!
    big_H.set_domain(2) = auxi.domain(2);
  
    big_H.std_spectral_base();
}

Scalar Monopole::compute_source_W() const {

 Scalar source (mp) ;
  source.allocate_all() ;
  
  // near r = 0 ;
  Scalar source_noyau (mp) ;
  source_noyau = small_w * small_w * small_w + 3* small_w*small_w/radius + ( 1 + radius*small_w)* big_H*big_H/radius;
  source_noyau.set_inner_boundary(0,0) ;
  source.set_domain(0) = source_noyau.domain(0) ;
  
  // CED :
  Scalar source_zec (mp) ;
  Scalar big_Wp = big_W.dsdr() ;
  source_zec =small_h*big_W*(small_h+2.)+big_W*(big_W*big_W-1)/(radius*radius)+2*big_Wp/radius ;
  source_zec.set_dzpuis(2) ;// Il est généralement fixé a deux dans les codes... vérifier sa véritable importance... 
  source.set_domain(nz-1) = source_zec.domain(nz-1) ;
  source.set_dzpuis(source_zec.get_dzpuis()) ;
  
  // Coquilles :
  Scalar source_H (mp) ;
 Scalar big_Wp = big_W.dsdr() ;
  source_H = big_W*(big_H*big_H-1)+big_W*(big_W*big_W-1)/ (radius*radius) + 2*big_Wp / radius ;
  for (int i=1 ; i<nz-1 ; i++)
      source.set_domain(i) = source_H.domain(i) ;

  source.std_spectral_base() ;
  return source ;

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
  source_zec = 2*big_W*big_W*(small_h+1)/(radius*radius) ;// erreur??? j'ai rajouté le r²)
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

