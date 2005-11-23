#include "map.h"
#include "scalar.h"
#include "monopole.h"
#include "utilitaires.h"

// Standard constructor :
Monopole::Monopole (const Map_af& mpi, double bet)  : mp(mpi), radius (mp), beta(bet), 
  big_W (mp), small_w (mp), big_H (mp), small_h(mp) {
  
  radius.set_etat_qcq() ;
  radius = mp.r ;
  nz = mp.get_mg()->get_nzone() ;
}
  

// Constructor by copy 
Monopole::Monopole (const Monopole& so) : mp(so.mp), radius(so.radius), nz(so.nz),
  beta (so.beta), big_W (so.big_W),small_w (so.small_w), 
  big_H(so.big_H),  small_h(so.small_h) {} 

// Constructor from a file :
Monopole::Monopole (const Map_af& mpi, FILE* fiche) : 
  mp(mpi), radius(mp), big_W (mp), small_w(mp), big_H (mp), small_h(mp) {
  
  nz = mp.get_mg()->get_nzone() ;
  fread_be (&beta, sizeof(double), 1, fiche) ;
  
  radius.set_etat_qcq() ;
  radius = mp.r ;

  Scalar big_W_file (mp, *mp.get_mg(), fiche) ;
  big_W = big_W_file ;
  
  Scalar small_w_file (mp, *mp.get_mg(), fiche) ;
  small_w = small_w_file ;

  Scalar big_H_file (mp, *mp.get_mg(), fiche) ;
  big_H = big_H_file ;
  
  Scalar small_h_file (mp, *mp.get_mg(), fiche) ;
  small_h = small_h_file ;
}

// Destructor
Monopole::~Monopole() {} ;

// Assignement :
void Monopole::operator= (const Monopole& so) {
  assert (mp == so.mp) ;
  beta = so.beta ;

  radius = so.radius ;
  nz = so.nz ;
  
  big_W = so.big_W ;
  small_w = so.small_w ;
  big_H = so.big_H ;
  small_h = so.small_h ;
}

// Write to a file :
void Monopole::sauve (FILE* fiche) const {

  fwrite_be (&beta, sizeof(double), 1, fiche) ;

  big_W.sauve(fiche) ;
  small_w.sauve (fiche) ;
  big_H.sauve(fiche) ;
  small_h.sauve (fiche) ;
}

  
