/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char map_radial_reevaluate_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:13:14  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2007/05/15 12:43:57  p_grandclement
 * Scalar version of reevaluate
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/01/04  15:24:00  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>

// Headers Lorene
#include "map.h"
#include "cmp.h"
#include "param.h"
#include "scalar.h"

void Map_radial::reevaluate(const Map* mp_prev0, int nzet, Cmp& uu) const { 
    
    const Map_radial* mp_prev = dynamic_cast<const Map_radial*>(mp_prev0) ; 

    if (mp_prev == 0x0) {
	cout << 
	    "Map_radial::reevaluate : the mapping mp_prev does not belong"
	    << endl ; 
	cout << " to the class Map_radial !" << endl ; 
	abort() ;  
    }

    int nz = mg->get_nzone() ; 

    // Protections
    // -----------
    assert(uu.get_mp() == this) ; 
    assert(uu.get_dzpuis() == 0) ; 
    assert(uu.get_etat() != ETATNONDEF) ; 
    assert(mp_prev->mg == mg) ; 
    assert(nzet <= nz) ; 
    
    
    // Maybe nothing to do ?
    if ( uu.get_etat() == ETATZERO ) {
	return ; 
    }
    
    assert(uu.get_etat() == ETATQCQ) ; 
    (uu.va).coef() ;	// the coefficients of uu are required 

    // Copy of the coefficients
    Mtbl_cf va_cf = *((uu.va).c_cf) ;
    
    // Preparation of uu.va for the storage of the new values.
    // ------------------------------------------------------
    
    (uu.va).set_etat_c_qcq() ;	// Allocates the memory for the Mtbl uu.va.c
				//  if it does not exist already
				// Destroys the coefficients

    Mtbl& va_c = *((uu.va).c) ; // Reference to the Mtbl uu.va.c
				
    va_c.set_etat_qcq() ;	// Allocates the memory for the Tbl's in each
				//  domain if they do not exist already


    // Values of the Coord r
    // ---------------------

    if ( r.c == 0x0 ) r.fait() ; 
    const Mtbl& rc = *(r.c) ; 
    
    // Precision for val_lx_jk :
    // ------------------------
    int nitermax = 100 ; // Maximum number of iterations in the secant method
    int niter ;		 // Number of iterations effectively used 
    double precis = 1.e-15 ; // Absolute precision in the secant method
    Param par ; 
    par.add_int(nitermax) ; 
    par.add_int_mod(niter) ; 
    par.add_double(precis) ; 
    
    // Loop on the domains where the evaluation is to be performed
    // -----------------------------------------------------------
    
    for (int l=0; l<nzet; l++) {
	int nr = mg->get_nr(l) ; 
	int nt = mg->get_nt(l) ; 
	int np = mg->get_np(l) ; 

	va_c.t[l]->set_etat_qcq() ;	// Allocates the array of double to 
					//  store the result 

	double* ptx = (va_c.t[l])->t ;	// Pointer on the allocated array
	
	double* pr = (rc.t[l])->t ;	// Pointer on the values of r

	// Loop on all the grid points in the considered domain

	for (int k=0; k<np; k++) {
	    for (int j=0; j<nt; j++) {
		for (int i=0; i<nr; i++) {

		    // domain l0 and value of the coordinate xi0 of the 
		    //  considered point in the previous mapping
		
		    int l0 ; 
		    double xi0 ; 
		    mp_prev->val_lx_jk(*pr, j, k, par, l0, xi0) ;
		
		    // Value of uu at this point
		
		    *ptx = va_cf.val_point_jk(l0, xi0, j, k) ; 
		    
		    // next point
		    pr++ ; 
		    ptx++ ; 
		}
	    }
	}

	
    }  // End of the loop on the domains where the evaluation had to be performed

    // In the remaining domains, uu is set to zero:
    // -------------------------------------------
    
    uu.annule(nzet, nz - 1) ; 
    
    
    
    
}

void Map_radial::reevaluate(const Map* mp_prev0, int nzet, Scalar& uu) const { 
    
    const Map_radial* mp_prev = dynamic_cast<const Map_radial*>(mp_prev0) ; 

    if (mp_prev == 0x0) {
	cout << 
	    "Map_radial::reevaluate : the mapping mp_prev does not belong"
	    << endl ; 
	cout << " to the class Map_radial !" << endl ; 
	abort() ;  
    }

    int nz = mg->get_nzone() ; 

    // Protections
    // -----------
    assert(uu.get_mp() == *this) ; 
    assert(uu.get_dzpuis() == 0) ; 
    assert(uu.get_etat() != ETATNONDEF) ; 
    assert(mp_prev->mg == mg) ; 
    assert(nzet <= nz) ; 
    
    
    // Maybe nothing to do ?
    if ( uu.get_etat() == ETATZERO ) {
	return ; 
    }
    
    assert(uu.get_etat() == ETATQCQ) ; 
    uu.set_spectral_va().coef() ;	// the coefficients of uu are required 

    // Copy of the coefficients
    Mtbl_cf va_cf = *(uu.set_spectral_va().c_cf) ;
    
    // Preparation of uu.va for the storage of the new values.
    // ------------------------------------------------------
    
    uu.set_spectral_va().set_etat_c_qcq() ;	// Allocates the memory for the Mtbl uu.va.c
				//  if it does not exist already
				// Destroys the coefficients

    Mtbl& va_c = *(uu.set_spectral_va().c) ; // Reference to the Mtbl uu.va.c
				
    va_c.set_etat_qcq() ;	// Allocates the memory for the Tbl's in each
				//  domain if they do not exist already


    // Values of the Coord r
    // ---------------------

    if ( r.c == 0x0 ) r.fait() ; 
    const Mtbl& rc = *(r.c) ; 
    
    // Precision for val_lx_jk :
    // ------------------------
    int nitermax = 100 ; // Maximum number of iterations in the secant method
    int niter ;		 // Number of iterations effectively used 
    double precis = 1.e-15 ; // Absolute precision in the secant method
    Param par ; 
    par.add_int(nitermax) ; 
    par.add_int_mod(niter) ; 
    par.add_double(precis) ; 
    
    // Loop on the domains where the evaluation is to be performed
    // -----------------------------------------------------------
    
    for (int l=0; l<nzet; l++) {
	int nr = mg->get_nr(l) ; 
	int nt = mg->get_nt(l) ; 
	int np = mg->get_np(l) ; 

	va_c.t[l]->set_etat_qcq() ;	// Allocates the array of double to 
					//  store the result 

	double* ptx = (va_c.t[l])->t ;	// Pointer on the allocated array
	
	double* pr = (rc.t[l])->t ;	// Pointer on the values of r

	// Loop on all the grid points in the considered domain

	for (int k=0; k<np; k++) {
	    for (int j=0; j<nt; j++) {
		for (int i=0; i<nr; i++) {

		    // domain l0 and value of the coordinate xi0 of the 
		    //  considered point in the previous mapping
		
		    int l0 ; 
		    double xi0 ; 
		    mp_prev->val_lx_jk(*pr, j, k, par, l0, xi0) ;
		
		    // Value of uu at this point
		
		    *ptx = va_cf.val_point_jk(l0, xi0, j, k) ; 
		    
		    // next point
		    pr++ ; 
		    ptx++ ; 
		}
	    }
	}

	
    }  // End of the loop on the domains where the evaluation had to be performed

    // In the remaining domains, uu is set to zero:
    // -------------------------------------------
    
    uu.annule(nzet, nz - 1) ; 
    
    
    
    
}
