/*
 * Method of class Etoile_rot to compute eccentric orbits
 *
 * (see file etoile.h for documentation).
 *
 */

/*
 *   Copyright (c) 2001 Dorota Gondek-Rosinska
 *   Copyright (c) 2001 Eric Gourgoulhon
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


char et_rot_f_eccentric_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  2001/02/08  15:13:24  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Headers C
#include <math.h>

// Headers Lorene
#include "etoile.h"
#include "param.h"
#include "utilitaires.h"

//=============================================================================
//		r_isco()
//=============================================================================

double Etoile_rot::f_eccentric(double ecc, double periast, ostream* ost) const {

    cout << "Etoile_rot::f_eccentric not ready yet !" << endl ; 
    abort() ; 

    // First and second derivatives of metric functions
    // ------------------------------------------------

    int nzm1 = mp.get_mg()->get_nzone() - 1 ;
    Cmp dnphi = nphi().dsdr() ;
    dnphi.annule(nzm1) ;
    Cmp ddnphi = dnphi.dsdr() ;    	// d^2/dr^2 N^phi

    Cmp tmp = nnn().dsdr() ;
    tmp.annule(nzm1) ;
    Cmp ddnnn = tmp.dsdr() ; 		// d^2/dr^2 N

    tmp = bbb().dsdr() ;
    tmp.annule(nzm1) ;
    Cmp ddbbb = tmp.dsdr() ; 		// d^2/dr^2 B

    // Constructing the velocity of a particle corotating with the star
    // ----------------------------------------------------------------

    Cmp bdlog = bbb().dsdr() / bbb() ;
    Cmp ndlog = nnn().dsdr() / nnn() ;
    Cmp bsn = bbb() / nnn() ;

    Cmp r(mp) ;
    r = mp.r ;

    Cmp r2= r*r ;

    bdlog.annule(nzm1) ;
    ndlog.annule(nzm1) ;
    bsn.annule(nzm1) ;
    r2.annule(nzm1) ;

    // ucor_plus - the velocity of corotating particle on the circular orbit
    Cmp ucor_plus = (r2 * bsn * dnphi +
        sqrt ( r2 * r2 * bsn *bsn * dnphi * dnphi +
		4 * r2 * bdlog * ndlog + 4 * r * ndlog) ) /
		2 / (1 + r * bdlog ) ;

    Cmp factor_u2 = r2 * (2 * ddbbb / bbb() - 2 * bdlog * bdlog +
    					  4 * bdlog * ndlog ) +
       2 * r2 * r2 * bsn * bsn * dnphi * dnphi +
       4 * r * ( ndlog - bdlog ) - 6 ;

    Cmp factor_u1 = 2 * r * r2 * bsn * ( 2 * ( ndlog - bdlog ) *
       				dnphi - ddnphi ) ;

    Cmp factor_u0 = - r2 * ( 2 * ddnnn / nnn() - 2 * ndlog * ndlog +
					 4 * bdlog * ndlog ) ;

    // Scalar field the zero of which will give the marginally stable orbit
    Cmp orbit = factor_u2 * ucor_plus * ucor_plus +
				factor_u1 * ucor_plus + factor_u0 ;

    // Search for the zero
    // -------------------

    int l_ms = nzet ;  // index of the domain containing the MS orbit


    Param par_ms ;
    par_ms.add_int(l_ms) ;
    par_ms.add_cmp(orbit) ;

    // Preliminary location of the zero
    // of the orbit function with an error = 0.01
    // The orbit closest to the star
    double theta_ms = M_PI / 2. ; // pi/2
    double phi_ms = 0. ;

    double xi_min = -1. ;
    double xi_max = 1. ;

    double resloc_old ;
    double xi_f = xi_min;

    orbit.std_base_scal() ;
    const Valeur& vorbit = orbit.va ;

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
		"Etoile_rot::isco : preliminary location of zero of MS function :"
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
 		(bbb().va).val_point(l_ms, xi_ms, theta_ms, phi_ms) * r_ms
    					  ) ;

	// Determination of the frequency at the marginally stable orbit
  	// -------------------------------------------------------------				

    ucor_plus.std_base_scal() ;
	double ucor_msplus = (ucor_plus.va).val_point(l_ms, xi_ms, theta_ms,
												  phi_ms) ;
	double nobrs = (bsn.va).val_point(l_ms, xi_ms, theta_ms, phi_ms) ;
	double nphirs = (nphi().va).val_point(l_ms, xi_ms, theta_ms, phi_ms) ;
	
 	p_f_isco = new double ( ( ucor_msplus / nobrs / r_ms + nphirs ) /
                       			(double(2) * M_PI) ) ;

    return 0 ;

}









