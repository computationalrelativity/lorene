/*
 * Method of class Etoile_rot to compute the location of the ISCO
 *
 * (see file etoile.h for documentation).
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char et_rot_isco_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/06 15:13:09  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2014/07/04 12:09:06  j_novak
 * New argument in zerosec(): a boolean (false by default) for aborting if the number of iteration is greater than the max.
 *
 * Revision 1.3  2011/01/07 18:20:08  m_bejger
 * Correcting for the case of stiff EOS, in which ISCO may be farther than the first domain outside the star - now searching all non-compactified domains
 *
 * Revision 1.2  2001/12/06 15:11:43  jl_zdunik
 * Introduction of the new function f_eq() in the class Etoile_rot
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2001/03/26  09:31:13  jlz
 * New functions: espec_isco() and lspec_isco().
 *
 * Revision 2.1  2000/11/18  23:18:08  eric
 * Ajout de l'argument ost a Etoile_rot::r_isco. Ajout de l'argument ost a Etoile_rot::r_isco.
 *
 * Revision 2.0  2000/11/18  21:10:41  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene
#include "etoile.h"
#include "param.h"
#include "utilitaires.h"

double fonct_etoile_rot_isco(double, const Param&) ;


//=============================================================================
//		r_isco()
//=============================================================================

double Etoile_rot::r_isco(ostream* ost) const {

    if (p_r_isco == 0x0) {    // a new computation is required

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
    orbit.std_base_scal() ;
    
    // Search for the zero
    // -------------------

    double r_ms, theta_ms, phi_ms, xi_ms, 
    	   xi_min = -1, xi_max = 1;  
    int l_ms = nzet, l ;
    bool find_status = false ; 

    const Valeur& vorbit = orbit.va ;
    
	for(l = nzet; l <= nzm1; l++) { 

    // Preliminary location of the zero
    // of the orbit function with an error = 0.01
    theta_ms = M_PI / 2. ; // pi/2
    phi_ms = 0. ;

    xi_min = -1. ;
    xi_max = 1. ;

    double resloc_old ;
    double xi_f = xi_min;
    
    double resloc = vorbit.val_point(l, xi_min, theta_ms, phi_ms) ;
	
	for (int iloc=0; iloc<200; iloc++) {
		xi_f = xi_f + 0.01 ;
     	resloc_old = resloc ;
     	resloc = vorbit.val_point(l, xi_f, theta_ms, phi_ms) ;
     	if ( resloc * resloc_old < double(0) ) {
			xi_min = xi_f - 0.01 ;
			xi_max = xi_f ;
			l_ms = l ; 
			find_status = true ;  
			break ;
		} 
		
  	  }
  	
    } 
  	
  	Param par_ms ;
    par_ms.add_int(l_ms) ;
    par_ms.add_cmp(orbit) ;

    if(find_status) { 
      
      double precis_ms = 1.e-13 ;    // precision in the determination of xi_ms
      int nitermax_ms = 200 ;	       // max number of iterations
      
      int niter ;
      xi_ms = zerosec(fonct_etoile_rot_isco, par_ms, xi_min, xi_max,
		      precis_ms, nitermax_ms, niter, false) ;
      
      if (niter > nitermax_ms) {
	cerr << "Etoile_rot::r_isco : " << endl ;
	cerr << "result may be wrong ... " << endl ;
	cerr << "Continuing." << endl ;
      }

      if (ost != 0x0) {
	* ost <<
	  "    number of iterations used in zerosec to locate the ISCO : "
	      << niter << endl ;
	*ost << "    zero found at xi = " << xi_ms << endl ;
      }
      
      r_ms = mp.val_r(l_ms, xi_ms, theta_ms, phi_ms) ;
      
    } else { 
      
      // assuming the ISCO is under the surface of a star 
      r_ms = ray_eq() ; 
      xi_ms = -1 ; 
      l_ms = nzet ; 
      
    } 
    
    p_r_isco = new double ((bbb().va).val_point(l_ms, xi_ms, theta_ms, phi_ms)
			   * r_ms ) ;
    
	// Determination of the frequency at the marginally stable orbit
  	// -------------------------------------------------------------				

    ucor_plus.std_base_scal() ;
	double ucor_msplus = (ucor_plus.va).val_point(l_ms, xi_ms, theta_ms,
												  phi_ms) ;
	double nobrs = (bsn.va).val_point(l_ms, xi_ms, theta_ms, phi_ms) ;
	double nphirs = (nphi().va).val_point(l_ms, xi_ms, theta_ms, phi_ms) ;
	
 	p_f_isco = new double ( ( ucor_msplus / nobrs / r_ms + nphirs ) /
                       			(double(2) * M_PI) ) ;

	// Specific angular momentum on ms orbit
	// -------------------------------------				
	p_lspec_isco=new double (ucor_msplus/sqrt(1.-ucor_msplus*ucor_msplus)*
	((bbb().va).val_point(l_ms, xi_ms, theta_ms, phi_ms)) * r_ms );

	// Specific energy on ms orbit
	// ---------------------------			
	p_espec_isco=new double (( 1./nobrs / r_ms / ucor_msplus  + nphirs) *
	(ucor_msplus/sqrt(1.-ucor_msplus*ucor_msplus)*
	((bbb().va).val_point(l_ms, xi_ms, theta_ms, phi_ms)) * r_ms ));

    // Determination of the Keplerian frequency at the equator
	// -------------------------------------------------------------

	double ucor_eqplus = (ucor_plus.va).val_point(l_ms, -1, theta_ms, phi_ms) ;
	double nobeq = (bsn.va).val_point(l_ms, -1, theta_ms, phi_ms) ;
	double nphieq = (nphi().va).val_point(l_ms, -1, theta_ms, phi_ms) ;

	p_f_eq = new double ( ( ucor_eqplus / nobeq / ray_eq() + nphieq ) /
			      (double(2) * M_PI) ) ;

    }  // End of computation

    return *p_r_isco ;

}

//=============================================================================
//		f_isco()
//=============================================================================

double Etoile_rot::f_isco() const {

    if (p_f_isco == 0x0) {    // a new computation is required

    	r_isco() ; 		// f_isco is computed in the method r_isco()

    	assert(p_f_isco != 0x0) ;
    }

    return *p_f_isco ;

}

//=============================================================================
//		lspec_isco()
//=============================================================================

double Etoile_rot::lspec_isco() const {

    if (p_lspec_isco == 0x0) {    // a new computation is required

    	r_isco() ; 	// lspec_isco is computed in the method r_isco()

    	assert(p_lspec_isco != 0x0) ;
    }

    return *p_lspec_isco ;

}

//=============================================================================
//		espec_isco()
//=============================================================================

double Etoile_rot::espec_isco() const {

    if (p_espec_isco == 0x0) {    // a new computation is required

    	r_isco() ; 	// espec_isco is computed in the method r_isco()

    	assert(p_espec_isco != 0x0) ;
    }

    return *p_espec_isco ;

}


//=============================================================================
//              f_eq()
//=============================================================================

double Etoile_rot::f_eq() const {
  
  if (p_f_eq == 0x0) {    // a new computation is required

    r_isco() ;              // f_eq is computed in the method r_isco()

    assert(p_f_eq != 0x0) ;
  }

  return *p_f_eq ;

}


//=============================================================================
//	Function used to locate the MS orbit
//=============================================================================


double fonct_etoile_rot_isco(double xi, const Param& par){

    // Retrieval of the information:
    int l_ms = par.get_int() ;
    const Cmp& orbit = par.get_cmp() ;
    const Valeur& vorbit = orbit.va ;

    // Value et the desired point
    double theta = M_PI / 2. ;
    double phi = 0 ;
    return vorbit.val_point(l_ms, xi, theta, phi) ;

}

