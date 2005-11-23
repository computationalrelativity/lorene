#include "map.h"
#include "scalar.h"
#include "monopole.h"
#include "utilitaires.h"

// Standard constructor :
Monopole::Monopole (const Map_af& mpi_W, const Map_af& mpi_H, 
		    double bet)  : mp_W(mpi_W), mp_H(mpi_H),  
  radius_on_W(mp_W), radius_on_H(mp_H), beta(bet), 
  big_W (mp_W), small_w (mp_W), big_H (mp_H), small_h(mp_H), 
  big_W_on_H (mp_H), small_w_on_H (mp_H), 
  big_H_on_W (mp_W), small_h_on_W(mp_W) {
  
  radius_on_W.set_etat_qcq() ;
  radius_on_W = mp_W.r ;
  nz_W = mp_W.get_mg()->get_nzone() ;
  
  radius_on_H.set_etat_qcq() ;
  radius_on_H = mp_H.r ;
  nz_H = mp_H.get_mg()->get_nzone() ;
}
  

// Constructor by copy 
Monopole::Monopole (const Monopole& so) : mp_W(so.mp_W), mp_H(so.mp_H), 
  radius_on_W(so.radius_on_W), radius_on_H(so.radius_on_H),
  nz_W(so.nz_W), nz_H(so.nz_H),
  beta (so.beta),  big_W (so.big_W),small_w (so.small_w), 
  big_H(so.big_H),  small_h(so.small_h), big_W_on_H (so.big_W_on_H), 
  small_w_on_H (so.small_w_on_H), big_H_on_W(so.big_H_on_W), 
  small_h_on_W (so.small_h_on_W) {} ;

// Constructor from a file :
Monopole::Monopole (const Map_af& mpi_W, const Map_af& mpi_H, FILE* fiche) : 
  mp_W(mpi_W), mp_H(mpi_H), radius_on_W(mp_W), radius_on_H(mp_H),
  big_W (mp_W), small_w(mp_W), big_H (mp_H), small_h(mp_H), 
  big_W_on_H (mp_H), small_w_on_H(mp_H), big_H_on_W(mp_W), small_h_on_W (mp_W) {
  
  fread_be (&beta, sizeof(double), 1, fiche) ;
  
  radius_on_W.set_etat_qcq() ;
  radius_on_W = mp_W.r ;
  radius_on_H.set_etat_qcq() ;
  radius_on_H = mp_H.r ;

  nz_W = mp_W.get_mg()->get_nzone() ;
  nz_H = mp_H.get_mg()->get_nzone() ;
  
  Scalar big_W_file (mp_W, *mp_W.get_mg(), fiche) ;
  big_W = big_W_file ;
  
  Scalar small_w_file (mp_W, *mp_W.get_mg(), fiche) ;
  small_w = small_w_file ;

  Scalar big_H_file (mp_H, *mp_H.get_mg(), fiche) ;
  big_H = big_H_file ;
  
  Scalar small_h_file (mp_H, *mp_H.get_mg(), fiche) ;
  small_h = small_h_file ;

   Scalar big_W_on_H_file (mp_H, *mp_H.get_mg(), fiche) ;
  big_W_on_H = big_W_on_H_file ;
  
  Scalar small_w_on_H_file (mp_H, *mp_H.get_mg(), fiche) ;
  small_w_on_H = small_w_on_H_file ;

  Scalar big_H_on_W_file (mp_W, *mp_W.get_mg(), fiche) ;
  big_H_on_W = big_H_on_W_file ;
  
  Scalar small_h_on_W_file (mp_W, *mp_W.get_mg(), fiche) ;
  small_h_on_W = small_h_on_W_file ;
}

// Destructor
Monopole::~Monopole() {} ;

// Assignement :
void Monopole::operator= (const Monopole& so) {
  assert (mp_W == so.mp_W) ;
  assert (mp_H == so.mp_H) ;

  beta = so.beta ;

  radius_on_W = so.radius_on_W ;
  radius_on_H = so.radius_on_H ;
  
  nz_W = so.nz_W ;
  nz_H = so.nz_H ;
  
  big_W = so.big_W ;
  small_w = so.small_w ;
  big_H = so.big_H ;
  small_h = so.small_h ;
  
  big_W_on_H = so.big_W_on_H ;
  small_w_on_H = so.small_w_on_H ;
  big_H_on_W = so.big_H_on_W ;
  small_h_on_W = so.small_h_on_W ;
}

// Write to a file :
void Monopole::sauve (FILE* fiche) const {

  fwrite_be (&beta, sizeof(double), 1, fiche) ;

  big_W.sauve(fiche) ;
  small_w.sauve (fiche) ;
  big_H.sauve(fiche) ;
  small_h.sauve (fiche) ;

  big_W_on_H.sauve(fiche) ;
  small_w_on_H.sauve (fiche) ;
  big_H_on_W.sauve(fiche) ;
  small_h_on_W.sauve (fiche) ;
}

  
