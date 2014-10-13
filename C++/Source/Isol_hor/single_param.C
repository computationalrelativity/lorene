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

char single_param_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:01  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:11  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2007/04/13 15:28:35  f_limousin
 * Lots of improvements, generalisation to an arbitrary state of
 * rotation, implementation of the spatial metric given by Samaya.
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>
#include <cassert>

// Lorene headers
#include "isol_hor.h"
#include "metric.h"
#include "evolution.h"
#include "unites.h"
#include "scalar.h"
#include "vector.h"
#include "graphique.h"
#include "utilitaires.h"



namespace Lorene {
const Scalar Single_hor::b_tilde()const {

  Scalar tmp = contract( beta, 0, tgam.radial_vect()
			 .down(0, tgam), 0) ;
  
  return tmp ;

}

const Scalar Single_hor::darea_hor() const {
  
  Scalar tmp = sqrt( get_gam().cov()(2,2) * get_gam().cov()(3,3) - 
		     get_gam().cov()(2,3) * get_gam().cov()(2,3)) ;
  
  tmp.std_spectral_base() ;
  
  return tmp ;
  
}

double Single_hor::area_hor() const {
    
    Scalar integrand (darea_hor()) ;
    integrand.raccord(1) ;

    return mp.integrale_surface(integrand, radius + 1e-15) ;

}

double Single_hor::radius_hor() const {

  double resu =  area_hor() / (4. * M_PI);

  resu = pow(resu, 1./2.) ;

  return resu ;

}

double Single_hor::ang_mom_hor()const {

  // Vector \partial_phi
  Vector phi (ff.get_mp(), CON, *(ff.get_triad()) ) ;

  Scalar tmp (ff.get_mp() ) ;
  tmp = 1 ;
  tmp.std_spectral_base() ;
  tmp.mult_rsint() ;

  phi.set(1) = 0. ;
  phi.set(2) = 0. ;
  phi.set(3) = tmp ; 
  
  Scalar k_rphi = contract(contract( get_gam().radial_vect(), 0, 
				     get_k_dd(), 0), 0, 
			   phi, 0) / (8. * M_PI) ;

  Scalar integrand = k_rphi * darea_hor() ;   // we correct with the curved 
                                              // element of area 

  double ang_mom = mp.integrale_surface(integrand, radius + 1e-15) ;

  return ang_mom ;

}

// Mass  (fundamental constants made 1)
double Single_hor::mass_hor()const {
  
  double rr = radius_hor() ;

  double  tmp = sqrt( pow( rr, 4) + 4 * pow( ang_mom_hor(), 2) ) / ( 2 * rr ) ;
									      
  return tmp ;

}

// Surface gravity
double Single_hor::kappa_hor() const{
  
  double rr = radius_hor() ;

  double jj = ang_mom_hor() ;

  double tmp = (pow( rr, 4) - 4 * pow( jj, 2)) / ( 2 * pow( rr, 3) 
			 *  sqrt( pow( rr, 4) + 4 * pow( jj, 2) ) ) ;
  
  return tmp ;

}

// Orbital velocity
double Single_hor::omega_hor()const {
  
  double rr = radius_hor() ;

  double jj = ang_mom_hor() ;

  double tmp = 2 * jj / ( rr *  sqrt( pow( rr, 4) + 4 * pow( jj, 2) ) ) ;
  
  return tmp ;

}

// ADM angular momentum

double Single_hor::ang_mom_adm()const {

  Scalar integrand =  (get_k_dd()(1,3) - get_gam().cov()(1,3) * trK) / 
    (8. * M_PI)  ;

  integrand.mult_rsint() ;  // in order to pass from the triad 
                            // component to the coordinate basis

  double tmp = mp.integrale_surface_infini(integrand) ;

  return  tmp ;

}

// Expansion

Scalar Single_hor::expansion() const {

  Scalar expa = contract(get_gam().radial_vect().derive_cov(get_gam()), 0,1) 
    + contract(contract(get_k_dd(), 0, get_gam().radial_vect(), 0), 
	       0, get_gam().radial_vect(), 0) - trK ; 

  return expa ;
}
}
