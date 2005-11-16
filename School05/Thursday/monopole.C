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

// Destructor
Monopole::~Monopole() {}

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


  
