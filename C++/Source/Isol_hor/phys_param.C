/*
 *  Method of class Isol_hor to compute physical parameters of the horizon
 *
 *    (see file isol_hor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Jose Luis Jaramillo
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

char phys_param_C[] = "$Header$" ;

/*
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>
#include <assert.h>

// Lorene headers
#include "time_slice.h"
#include "isol_hor.h"
#include "metric.h"
#include "evolution.h"
#include "unites.h"
#include "scalar.h"
#include "vector.h"
#include "graphique.h"
#include "utilitaires.h"



Vector Isol_hor::radial_vect_hor()  {

  Vector get_radial_vect (ff.get_mp(), CON, *(ff.get_triad()) ) ;
       
  get_radial_vect.set(1) = gam_uu()(1,1) ;
 
  get_radial_vect.set(2) = gam_uu()(1,2) ;

  get_radial_vect.set(3) = gam_uu()(1,3) ;

  get_radial_vect = get_radial_vect / (psi() * psi() * sqrt(gam_uu()(1,1))) ;

  get_radial_vect.std_spectral_base() ;


  return get_radial_vect ;

}


Vector Isol_hor::beta_bound_cart() {


  const Map& map = ff.get_mp() ; 

  Vector tmp_vect = nn() * radial_vect_hor() ;

  tmp_vect.set_etat_zero() ;

  Scalar vel_ang (map) ;
  
  vel_ang = omega_hor() ;

  vel_ang.std_spectral_base() ;

  vel_ang.mult_rsint() ;

  tmp_vect.set(3) = tmp_vect(3)  -  vel_ang ;

  

  tmp_vect.change_triad(map.get_bvect_cart() ) ;
 
  /*
     for (int i=1 ; i<4 ; i++) {
	cout<<"component: ("<<i<<") of beta_cart"<<endl ;      
	des_profile(tmp_vect(i), 1.000001, 10., M_PI/2., 0., "vector comp beta_cart(tmp) ", "radial distance") ;
      }
  */
  return tmp_vect ;

}








Scalar Isol_hor::darea_hor() {
  
  Scalar tmp = sqrt( gam_dd()(2,2) * gam_dd()(3,3) - gam_dd()(2,3) * gam_dd()(2,3)) ;
  
  tmp.std_spectral_base() ;
  
  return tmp ;
  
}



double Isol_hor::radius_hor()  {

  //  if( !(radius_evol.is_known(jtime)) {

    //    assert(   .isknown(jtime) ) ;

  //  Scalar tmp = sqrt( gam_dd()(2,2) * gam_dd()(3,3) - gam_dd()(2,3) * gam_dd()(2,3)) ;

  //  tmp.std_spectral_base() ;



  Map_af map_affine (ff.get_mp()) ; 

  double resu =  map_affine.integrale_surface(darea_hor(), 1.0000000001) / (4. * M_PI);

    //   radius_evol.update(resu, jtime, thetime(jtime)) ;

    //  }

    //  return radius_evol[jtime] ;

  return resu ;

}




double Isol_hor::ang_mom_hor() {


  // Vector \partial_phi
  Vector phi (ff.get_mp(), CON, *(ff.get_triad()) ) ;

  Scalar tmp (ff.get_mp() ) ;
  tmp = 1 ;
  tmp.std_spectral_base() ;

  tmp.mult_rsint() ;

  phi.set(1) = 0. ;
  phi.set(2) = 0. ;
  phi.set(3) = tmp ;
  
  phi.std_spectral_base() ;
 
  
  Scalar k_rphi = contract(contract( radial_vect_hor(), 0, k_dd(), 0), 0, phi, 0) / (8. * M_PI) ;

  Scalar integrand = k_rphi * darea_hor() ;   // we correct with the curved element of area 

  Map_af map_affine (ff.get_mp() ) ;

  double tmp_ang_mom = map_affine.integrale_surface(integrand, 1.0000000001) ;

  return tmp_ang_mom ;

}


// Mass  (fundamental constants made 1)
double Isol_hor::mass_hor() {
  
  double rr = radius_hor() ;

  double  tmp = sqrt( pow( rr, 4) + 4 * pow( ang_mom_hor(), 2) ) / ( 2 * rr ) ;
											
  return tmp ;

}





// Surface gravity
double Isol_hor::kappa_hor() {
  
  double rr = radius_hor() ;

  double jj = ang_mom_hor() ;

  double tmp = (pow( rr, 4) - 4 * pow( jj, 2)) / ( 2 * pow( rr, 3) *  sqrt( pow( rr, 4) + 4 * pow( jj, 2) ) ) ;
  

  return tmp ;


}



// Orbital velocity
double Isol_hor::omega_hor() {
  
  double rr = radius_hor() ;

  double jj = ang_mom_hor() ;

  double tmp = 2 * jj / ( rr *  sqrt( pow( rr, 4) + 4 * pow( jj, 2) ) ) ;
  

  return tmp ;


}

