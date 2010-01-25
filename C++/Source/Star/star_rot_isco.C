/*
 * Method of class Star_rot to compute the location of the ISCO
 *
 * (see file star_rot.h for documentation).
 *
 */

/*
 *   Copyright (c) 2010 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 J. Leszek Zdunik
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


char star_rot_isco_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2010/01/25 22:33:04  e_gourgoulhon
 * First implementation.
 *
 * Revision 1.1  2010/01/25 18:15:52  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <math.h>

// Headers Lorene
#include "star_rot.h"
#include "param.h"
#include "utilitaires.h"

double funct_star_rot_isco(double, const Param& ) ; 

//=============================================================================
//		r_isco()
//=============================================================================

double Star_rot::r_isco(ostream* ost) const {

    if (p_r_isco == 0x0) {    // a new computation is required

    // First and second derivatives of metric functions
    // ------------------------------------------------

    int nzm1 = mp.get_mg()->get_nzone() - 1 ;
    Scalar dnphi = nphi.dsdr() ;
    dnphi.annule_domain(nzm1) ;
    Scalar ddnphi = dnphi.dsdr() ;    	// d^2/dr^2 N^phi

    Scalar tmp = nn.dsdr() ;
    tmp.annule_domain(nzm1) ;
    Scalar ddnnn = tmp.dsdr() ; 		// d^2/dr^2 N

    tmp = bbb.dsdr() ;
    tmp.annule_domain(nzm1) ;
    Scalar ddbbb = tmp.dsdr() ; 		// d^2/dr^2 B

    // Constructing the velocity of a particle corotating with the star
    // ----------------------------------------------------------------

    Scalar bdlog = bbb.dsdr() / bbb ;
    Scalar ndlog = nn.dsdr() / nn ;
    Scalar bsn = bbb / nn ;

    Scalar r(mp) ;
    r = mp.r ;

    Scalar r2= r*r ;

    bdlog.annule_domain(nzm1) ;
    ndlog.annule_domain(nzm1) ;
    bsn.annule_domain(nzm1) ;
    r2.annule_domain(nzm1) ;

    // ucor_plus - the velocity of corotating particle on the circular orbit
    Scalar ucor_plus = (r2 * bsn * dnphi +
        sqrt ( r2 * r2 * bsn *bsn * dnphi * dnphi +
		4 * r2 * bdlog * ndlog + 4 * r * ndlog) ) /
		2 / (1 + r * bdlog ) ;

    Scalar factor_u2 = r2 * (2 * ddbbb / bbb - 2 * bdlog * bdlog +
    					  4 * bdlog * ndlog ) +
       2 * r2 * r2 * bsn * bsn * dnphi * dnphi +
       4 * r * ( ndlog - bdlog ) - 6 ;

    Scalar factor_u1 = 2 * r * r2 * bsn * ( 2 * ( ndlog - bdlog ) *
       				dnphi - ddnphi ) ;

    Scalar factor_u0 = - r2 * ( 2 * ddnnn / nn - 2 * ndlog * ndlog +
					 4 * bdlog * ndlog ) ;

    // Scalar field the zero of which will give the marginally stable orbit
    Scalar orbit = factor_u2 * ucor_plus * ucor_plus +
				factor_u1 * ucor_plus + factor_u0 ;

    // Search for the zero
    // -------------------

    int l_ms = nzet ;  // index of the domain containing the MS orbit


    Param par_ms ;
    par_ms.add_int(l_ms) ;
    par_ms.add_scalar(orbit) ;

    // Preliminary location of the zero
    // of the orbit function with an error = 0.01
    // The orbit closest to the star
    double theta_ms = M_PI / 2. ; // pi/2
    double phi_ms = 0. ;

    double xi_min = -1. ;
    double xi_max = 1. ;

    double resloc_old ;
    double xi_f = xi_min;

    orbit.std_spectral_base() ;
    const Valeur& vorbit = orbit.get_spectral_va() ;

    double resloc = vorbit.val_point(l_ms, xi_min, theta_ms, phi_ms) ;
	
	for (int iloc=0; iloc<200; iloc++) {
		xi_f = xi_f + 0.01 ;
     	resloc_old = resloc ;
     	resloc = vorbit.val_point(l_ms, xi_f, theta_ms, phi_ms) ;
     	if ((resloc * resloc_old) < double(0) ) {
			xi_min = xi_f - 0.01 ;
			xi_max = xi_f ;
			break ;
		}
  	}
  	
  	
  	if (ost != 0x0) {
		*ost <<
		"Star_rot::isco : preliminary location of zero of MS function :"
		 << endl ;
		*ost << "    xi_min = " << xi_min << "  f(xi_min) = " <<
     		 vorbit.val_point(l_ms, xi_min, theta_ms, phi_ms) << endl ;
 		*ost << "    xi_max = " << xi_max << "  f(xi_max) = " <<
     		vorbit.val_point(l_ms, xi_max, theta_ms, phi_ms) << endl ;
    }
     	
    double xi_ms, r_ms ;  	
	
	if ( vorbit.val_point(l_ms, xi_min, theta_ms, phi_ms) *
 	     vorbit.val_point(l_ms, xi_max, theta_ms, phi_ms) < double(0) ) {

     	double precis_ms = 1.e-12 ;    // precision in the determination of xi_ms
     	int nitermax_ms = 100 ;	       // max number of iterations

     	int niter ;
     	xi_ms = zerosec(funct_star_rot_isco, par_ms, xi_min, xi_max,
     					precis_ms, nitermax_ms, niter) ;

  		if (ost != 0x0) {
     		* ost <<
     		"    number of iterations used in zerosec to locate the ISCO : "
	  		 << niter << endl ;
     		*ost << "    zero found at xi = " << xi_ms << endl ;
        }

      	r_ms = mp.val_r(l_ms, xi_ms, theta_ms, phi_ms) ;

     }
     else {
	    xi_ms = -1 ;
     	r_ms = ray_eq() ;
     }
  	
 	p_r_isco = new double (
 		(bbb.get_spectral_va()).val_point(l_ms, xi_ms, theta_ms, phi_ms) * r_ms
    					  ) ;

	// Determination of the frequency at the marginally stable orbit
  	// -------------------------------------------------------------				

	ucor_plus.std_spectral_base() ;
	double ucor_msplus = (ucor_plus.get_spectral_va()).val_point(l_ms, xi_ms, theta_ms,
												  phi_ms) ;
	double nobrs = (bsn.get_spectral_va()).val_point(l_ms, xi_ms, theta_ms, phi_ms) ;
	double nphirs = (nphi.get_spectral_va()).val_point(l_ms, xi_ms, theta_ms, phi_ms) ;
	
 	p_f_isco = new double ( ( ucor_msplus / nobrs / r_ms + nphirs ) /
                       			(double(2) * M_PI) ) ;

	// Specific angular momentum on ms orbit
	// -------------------------------------				
	p_lspec_isco=new double (ucor_msplus/sqrt(1.-ucor_msplus*ucor_msplus)*
	((bbb.get_spectral_va()).val_point(l_ms, xi_ms, theta_ms, phi_ms)) * r_ms );

	// Specific energy on ms orbit
	// ---------------------------			
	p_espec_isco=new double (( 1./nobrs / r_ms / ucor_msplus  + nphirs) *
	(ucor_msplus/sqrt(1.-ucor_msplus*ucor_msplus)*
	((bbb.get_spectral_va()).val_point(l_ms, xi_ms, theta_ms, phi_ms)) * r_ms ));

        // Determination of the Keplerian frequency at the equator
	// -------------------------------------------------------------


	double ucor_eqplus = (ucor_plus.get_spectral_va()).val_point(l_ms, -1, theta_ms,phi_ms)
	  ;
	double nobeq = (bsn.get_spectral_va()).val_point(l_ms, -1, theta_ms, phi_ms) ;
	double nphieq = (nphi.get_spectral_va()).val_point(l_ms, -1, theta_ms, phi_ms) ;

	p_f_eq = new double ( ( ucor_eqplus / nobeq / ray_eq() + nphieq ) /
			      (double(2) * M_PI) ) ;

	

    }  // End of computation

    return *p_r_isco ;

}



//=============================================================================
//		f_isco()
//=============================================================================

double Star_rot::f_isco() const {

    if (p_f_isco == 0x0) {    // a new computation is required

    	r_isco() ; 		// f_isco is computed in the method r_isco()

    	assert(p_f_isco != 0x0) ;
    }

    return *p_f_isco ;

}

//=============================================================================
//		lspec_isco()
//=============================================================================

double Star_rot::lspec_isco() const {

    if (p_lspec_isco == 0x0) {    // a new computation is required

    	r_isco() ; 	// lspec_isco is computed in the method r_isco()

    	assert(p_lspec_isco != 0x0) ;
    }

    return *p_lspec_isco ;

}

//=============================================================================
//		espec_isco()
//=============================================================================

double Star_rot::espec_isco() const {

    if (p_espec_isco == 0x0) {    // a new computation is required

    	r_isco() ; 	// espec_isco is computed in the method r_isco()

    	assert(p_espec_isco != 0x0) ;
    }

    return *p_espec_isco ;

}


//=============================================================================
//              f_eq()
//=============================================================================

double Star_rot::f_eq() const {
  
  if (p_f_eq == 0x0) {    // a new computation is required

    r_isco() ;              // f_eq is computed in the method r_isco()

    assert(p_f_eq != 0x0) ;
  }

  return *p_f_eq ;

}


//=============================================================================
//	Function used to locate the MS orbit
//=============================================================================


double funct_star_rot_isco(double xi, const Param& par){

    // Retrieval of the information:
    int l_ms = par.get_int() ;
    const Scalar& orbit = par.get_scalar() ;
    const Valeur& vorbit = orbit.get_spectral_va() ;

    // Value et the desired point
    double theta = M_PI / 2. ;
    double phi = 0 ;
    return vorbit.val_point(l_ms, xi, theta, phi) ;

}






