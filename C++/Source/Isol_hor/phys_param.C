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
 * $Id$
 * $Log$
 * Revision 1.10  2005/04/15 11:54:21  jl_jaramillo
 * function to compute the expansion of spherical surfaces
 *
 * Revision 1.9  2005/03/22 13:25:36  f_limousin
 * Small changes. The angular velocity and A^{ij} are computed
 * with a differnet sign.
 *
 * Revision 1.8  2005/03/03 10:10:14  f_limousin
 * Add the function area_hor().
 *
 * Revision 1.7  2005/02/07 10:35:42  f_limousin
 * Minor changes.
 *
 * Revision 1.6  2004/12/22 18:16:16  f_limousin
 * Mny different changes.
 *
 * Revision 1.5  2004/11/18 12:30:01  jl_jaramillo
 * Definition of b_tilde
 *
 * Revision 1.4  2004/10/29 15:44:13  jl_jaramillo
 * ADM angular momentum added.
 *
 * Revision 1.3  2004/09/17 13:37:21  f_limousin
 * Correction of an error in calculation of the radius
 *
 * Revision 1.2  2004/09/09 16:54:53  f_limousin
 * Add the 2 lines $Id$Log: for CVS
 *
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
#include "isol_hor.h"
#include "metric.h"
#include "evolution.h"
#include "unites.h"
#include "scalar.h"
#include "vector.h"
#include "graphique.h"
#include "utilitaires.h"


const Vector Isol_hor::radial_vect_hor() const {

  Vector get_radial_vect (ff.get_mp(), CON, *(ff.get_triad()) ) ;
       
  get_radial_vect.set(1) = gam_uu()(1,1) ;
 
  get_radial_vect.set(2) = gam_uu()(1,2) ;

  get_radial_vect.set(3) = gam_uu()(1,3) ;

  get_radial_vect = get_radial_vect / sqrt(gam_uu()(1,1)) ;

  get_radial_vect.std_spectral_base() ;


  return get_radial_vect ;

}


// Think of defining this as a pointer
const Vector Isol_hor::tradial_vect_hor() const {

  Vector get_radial_vect (ff.get_mp(), CON, *(ff.get_triad()) ) ;
       
  get_radial_vect.set(1) = (met_gamt.con())(1,1) ;
 
  get_radial_vect.set(2) = (met_gamt.con())(1,2) ;

  get_radial_vect.set(3) = (met_gamt.con())(1,3) ;

  get_radial_vect = get_radial_vect / sqrt((met_gamt.con())(1,1)) ;

  get_radial_vect.std_spectral_base() ;


  return get_radial_vect ;

}


const Scalar Isol_hor::b_tilde()const {

  Scalar tmp = contract( beta(), 0, met_gamt.radial_vect()
			 .down(0, met_gamt), 0) ;
  
  return tmp ;

}


const Scalar Isol_hor::darea_hor() const {
  
  Scalar tmp = sqrt( gam_dd()(2,2) * gam_dd()(3,3) - gam_dd()(2,3) 
		     * gam_dd()(2,3)) ;
  
  tmp.std_spectral_base() ;
  
  return tmp ;
  
}

double Isol_hor::area_hor() const {
    
    Scalar integrand (darea_hor()) ;
    integrand.raccord(1) ;

    return mp.integrale_surface(integrand, radius + 1e-15) ;

}


double Isol_hor::radius_hor() const {

  double resu =  area_hor() / (4. * M_PI);

  resu = pow(resu, 1./2.) ;

  return resu ;

}


double Isol_hor::ang_mom_hor()const {

  // Vector \partial_phi
  Vector phi (ff.get_mp(), CON, *(ff.get_triad()) ) ;

  Scalar tmp (ff.get_mp() ) ;
  tmp = 1 ;
  tmp.std_spectral_base() ;
  tmp.mult_rsint() ;

  phi.set(1) = 0. ;
  phi.set(2) = 0. ;
  phi.set(3) = tmp ; 
  
  Scalar k_rphi = contract(contract( radial_vect_hor(), 0, k_dd(), 0), 0, 
			   phi, 0) / (8. * M_PI) ;

  Scalar integrand = k_rphi * darea_hor() ;   // we correct with the curved 
                                              // element of area 

  double ang_mom = mp.integrale_surface(integrand, 1.0000000001) ;

  return ang_mom ;

}


// Mass  (fundamental constants made 1)
double Isol_hor::mass_hor()const {
  
  double rr = radius_hor() ;

  double  tmp = sqrt( pow( rr, 4) + 4 * pow( ang_mom_hor(), 2) ) / ( 2 * rr ) ;
									      
  return tmp ;

}


// Surface gravity
double Isol_hor::kappa_hor() const{
  
  double rr = radius_hor() ;

  double jj = ang_mom_hor() ;

  double tmp = (pow( rr, 4) - 4 * pow( jj, 2)) / ( 2 * pow( rr, 3) 
			 *  sqrt( pow( rr, 4) + 4 * pow( jj, 2) ) ) ;
  
  return tmp ;

}


// Orbital velocity
double Isol_hor::omega_hor()const {
  
  double rr = radius_hor() ;

  double jj = ang_mom_hor() ;

  double tmp = 2 * jj / ( rr *  sqrt( pow( rr, 4) + 4 * pow( jj, 2) ) ) ;
  
  return tmp ;

}


// ADM angular momentum

double Isol_hor::ang_mom_adm()const {

  Scalar integrand =  (k_dd()(1,3) - gam_dd()(1,3) * trk()) / (8. * M_PI)  ;

  integrand.mult_rsint() ;  // in order to pass from the triad 
                            // component to the coordinate basis

  double tmp = mp.integrale_surface_infini(integrand) ;

  return  tmp ;

}

// Expansion

Scalar Isol_hor::expansion() const {

  Scalar expa = contract(gam().radial_vect().derive_cov(gam()), 0,1) 
    + contract(contract(k_dd(), 0, gam().radial_vect(), 0), 
	       0, gam().radial_vect(), 0) - trk() ; 

  return expa ;
}
