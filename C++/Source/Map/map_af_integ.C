/*
 *  Method of the class Map_af for computing the integral of a Cmp over
 *  all space.
 */

/*
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
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


char map_af_integ_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2008/08/27 08:48:05  jl_cornou
 * Added R_JACO02 case
 *
 * Revision 1.5  2007/10/05 15:56:19  j_novak
 * Addition of new bases for the non-symmetric case in theta.
 *
 * Revision 1.4  2003/12/19 16:21:43  j_novak
 * Shadow hunt
 *
 * Revision 1.3  2003/10/19 19:58:15  e_gourgoulhon
 * Access to Base_val::b now via Base_val::get_b().
 *
 * Revision 1.2  2003/10/15 10:35:27  e_gourgoulhon
 * Changed cast (double *) into static_cast<double*>.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  2000/01/28  16:09:37  eric
 * Remplacement du ci.get_dzpuis() == 4 par ci.check_dzpuis(4).
 *
 * Revision 1.1  1999/12/09  10:45:43  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <math.h>


// Headers Lorene
#include "map.h"
#include "cmp.h"

Tbl* Map_af::integrale(const Cmp& ci) const {

    static double* cx_tcp = 0 ;	    // Coefficients theta, dev. en cos(2l theta)
    static double* cx_rrp = 0 ;	    // Coefficients r rare, dev. en T_{2i}
    static double* cx_rf_x2 = 0 ;	// Coefficients r fin, int. en r^2
    static double* cx_rf_x = 0 ;	// Coefficients r fin, int en r
    static double* cx_rf = 0 ;		// Coefficients r fin, int en cst.

    static int nt_cp_pre = 0 ;
    static int nr_p_pre = 0 ;
    static int nr_f_pre = 0 ;

    assert(ci.get_etat() != ETATNONDEF) ; 

    int nz = mg->get_nzone() ; 
    
    Tbl* resu = new Tbl(nz) ;
    
    if (ci.get_etat() == ETATZERO) {
	resu->annule_hard() ; 
	return resu ; 
    }
    
    assert( ci.get_etat() == ETATQCQ ) ; 
    
    (ci.va).coef() ;	// The spectral coefficients are required
    
    const Mtbl_cf* p_mti = (ci.va).c_cf ; 

    assert( ci.check_dzpuis(4) ) ;	    // dzpuis must be equal to 4  
    
    assert(p_mti->get_etat() == ETATQCQ) ; 

    resu->set_etat_qcq() ;  // Allocates the memory for the array of double

    // Loop of the various domains
    // ---------------------------
    for (int l=0 ; l<nz ; l++) {
	
	const Tbl* p_tbi = p_mti->t[l] ; 
	
	if ( p_tbi->get_etat() == ETATZERO ) {
	    resu->t[l] = 0 ; 
	}
	else {	    // Case where the computation must be done
	
	    assert( p_tbi->get_etat() == ETATQCQ ) ; 
	
	    int nt = mg->get_nt(l) ; 
	    int nr = mg->get_nr(l) ;
	    
	    int base = (p_mti->base).get_b(l) ; 
	    int base_r = base & MSQ_R ;
	    int base_t = base & MSQ_T ;
	    int base_p = base & MSQ_P ;
	    
	    // ----------------------------------
	    // Integration on theta -> array in r
	    // ----------------------------------

	    double* s_tr = new double[nr] ;	// Partial integral on theta 
	    double* x_spec = p_tbi->t ;	// Pointer on the spectral coefficients
	    
	    switch (base_t) {
	    
		case T_COS_P: case T_COSSIN_CP: {
		    if (nt > nt_cp_pre) {  // Initialization of factors for summation
			nt_cp_pre = nt ;
			cx_tcp 
			    = static_cast<double*>(realloc(cx_tcp, nt*sizeof(double))) ;
			for (int j=0 ; j<nt ; j++) {
			    cx_tcp[j] = 2./(1. - 4.*j*j) ;  // Factor 2 symmetry
			}
		    }
		    
		    // Summation : 
		    for (int i=0 ; i<nr ; i++) s_tr[i] = 0 ;
		    for (int j=0 ; j<nt ; j++) {
			for (int i=0 ; i<nr ; i++) {
			    s_tr[i] += cx_tcp[j] * x_spec[i] ;
			}
			x_spec += nr ;	// Next theta
		    }
		    break ;
		}
		case T_COSSIN_C: {
		    // Summation : 
		    for (int i=0 ; i<nr ; i++) s_tr[i] = 0 ;
		    for (int j=0 ; j<nt ; j++) {
			if ((j%2)==0)
			    for (int i=0 ; i<nr ; i++) {
				s_tr[i] += (2. / (1.-j*j)) * x_spec[i] ;
			    }
			x_spec += nr ;	// Next theta
		    }
		    break ;		    
		}
		default: {
		    cout << "Map_af::integrale: unknown theta basis ! " << endl ;
		    abort () ;
		    break ;
		}
		    
	    }	// End of the various cases on base_t

	    // ----------------
	    // Integration on r
	    // ----------------

	    double som = 0 ;
	    double som_x2 = 0;
	    double som_x = 0 ;
	    double som_c = 0 ;
	
	    switch(base_r) {
	    case R_CHEBP: case R_CHEBPIM_P: case R_CHEBPI_P :{
		assert(beta[l] == 0) ;
		if (nr > nr_p_pre) {  // Initialization of factors for summation
		    nr_p_pre = nr ;
		    cx_rrp = static_cast<double*>(realloc(cx_rrp, nr*sizeof(double))) ;
		    for (int i=0 ; i<nr ; i++) {
			cx_rrp[i] = (3. - 4.*i*i) / 
				    (9. - 40. * i*i + 16. * i*i*i*i) ;
		    }
		}
	    
		for (int i=0 ; i<nr ; i++) {
		    som += cx_rrp[i] * s_tr[i] ;
		}
		double rmax = alpha[l] ;
		som *= rmax*rmax*rmax ;
		break ;
	    }
	    
	    case R_CHEB: {
		if (nr > nr_f_pre) {  // Initialization of factors for summation
		    nr_f_pre = nr ;
		    cx_rf_x2 = static_cast<double*>(realloc(cx_rf_x2, nr*sizeof(double))) ;
		    cx_rf_x  = static_cast<double*>(realloc(cx_rf_x, nr*sizeof(double))) ;
		    cx_rf    = static_cast<double*>(realloc(cx_rf, nr*sizeof(double))) ;
		    for (int i=0 ; i<nr ; i +=2 ) {
			cx_rf_x2[i] = 2.*(3. - i*i)/(9. - 10. * i*i + i*i*i*i) ;
			cx_rf_x[i] = 0 ;
			cx_rf[i] = 2./(1. - i*i) ;
		    }
		    for (int i=1 ; i<nr ; i +=2 ) {
			cx_rf_x2[i] = 0 ;
			cx_rf_x[i] = 2./(4. - i*i) ;
			cx_rf[i] = 0 ;
		    }
		}
	    
		for (int i=0 ; i<nr ; i +=2 ) {
		    som_x2 += cx_rf_x2[i] * s_tr[i] ;
		}
		for (int i=1 ; i<nr ; i +=2 ) {
		    som_x += cx_rf_x[i] * s_tr[i] ;
		}
		for (int i=0 ; i<nr ; i +=2 ) {
		    som_c += cx_rf[i] * s_tr[i] ;
		}
		double a = alpha[l] ;
		double b = beta[l] ;
		som = a*a*a * som_x2 + 2.*a*a*b * som_x + a*b*b * som_c ;
		break ;
	    }

	    case R_JACO02: {
		if (nr > nr_f_pre) {  // Initialization of factors for summation
		    nr_f_pre = nr ;
		    cx_rf_x2 = static_cast<double*>(realloc(cx_rf_x2, nr*sizeof(double))) ;
		    cx_rf_x  = static_cast<double*>(realloc(cx_rf_x, nr*sizeof(double))) ;
		    cx_rf    = static_cast<double*>(realloc(cx_rf, nr*sizeof(double))) ;
		    double signe = 1 ;
		    for (int i=0 ; i<nr ; i +=1 ) {
			cx_rf_x2[i] = 0 ;
			cx_rf_x[i] = 2*signe/double(i+1)/double(i+2);
			cx_rf[i] = 2*signe ;
			signe = - signe ;
		    }
			cx_rf_x2[0] = double(8)/double(3) ;
		}
	    
		for (int i=0 ; i<nr ; i +=1 ) {
		    som_x2 += cx_rf_x2[i] * s_tr[i] ;
		}
		for (int i=1 ; i<nr ; i +=1 ) {
		    som_x += cx_rf_x[i] * s_tr[i] ;
		}
		for (int i=0 ; i<nr ; i +=1 ) {
		    som_c += cx_rf[i] * s_tr[i] ;
		}
		double a = alpha[l] ;
		double b = beta[l] ;
		assert(b == a) ;
		som = a*a*a * som_x2 + 2.*a*a*(b-a) * som_x + a*(b-a)*(b-a) * som_c ;
		break ;
	    }
	    
	    case R_CHEBU: {
		assert(beta[l] == - alpha[l]) ;
		if (nr > nr_f_pre) {  // Initialization of factors for summation
		    nr_f_pre = nr ;
		    cx_rf    = static_cast<double*>(realloc(cx_rf, nr*sizeof(double))) ;
		    for (int i=0 ; i<nr ; i +=2 ) {
			cx_rf[i] = 2./(1. - i*i) ;
		    }
		    for (int i=1 ; i<nr ; i +=2 ) {
			cx_rf[i] = 0 ;
		    }
		}
	    
		for (int i=0 ; i<nr ; i +=2 ) {
		    som_c += cx_rf[i] * s_tr[i] ;
		}
		som = - alpha[l] * som_c ;
		break ;
	    }


	    default: {
		cout << "Map_af::integrale: unknown r basis ! " << endl ;
		abort () ;
		break ;
	    }
	    
	    }	// End of the various cases on base_r
	    
	    // ------------------
	    // Integration on phi
	    // ------------------

	    switch (base_p) {

	    case P_COSSIN: {
		som *= 2. * M_PI ;
		break ;
	    }
	    
	    case P_COSSIN_P: {
		som *= 2. * M_PI ;
		break ;
	    }
	    
	    default: {
		cout << "Map_af::integrale: unknown phi basis ! " << endl ;
		abort () ;
		break ;
	    }
	    
	    }	// End of the various cases on base_p

	    // Final result for this domain:
	    // ----------------------------
	       
	    resu->t[l] = som ; 
	    delete [] s_tr ;

	}   // End of the case where the computation must be done
	
	
    }	// End of the loop onto the domains
        
    return resu ;
    
}
