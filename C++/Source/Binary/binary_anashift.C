/*
 * Method of class Binary to set some analytical form to the shift vector.
 *
 * (see file binary.h for documentation).
 */

/*
 *   Copyright (c) 2004 Francois Limousin
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


char binary_anashift_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/01/20 15:21:23  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "binary.h"

void Binary::analytical_shift(){
    
/*

    #include "unites.h"
    if (this == 0x0) {	// To avoid any compilation warning 
	cout << f_unit << msol << km << mevpfm3 << qpig ;
    }
        
    for (int i=0; i<2; i++) {

	// Radius of the star:
	double a0 = et[i]->ray_eq() ; 
    
	// Mass ratio
	double p_mass = et[i]->mass_g() / et[1-i]->mass_g() ; 
    
	// G M Omega R / (1+p) 
	double www = ggrav * et[i]->mass_g() * omega 
		    * separation() / (1. + p_mass) ;  
    
	const Map& mp = et[i]->get_mp() ; 
	Scalar tmp(mp) ;  
	Scalar tmp_ext(mp) ;  
	int nzet = et[i]->get_nzet() ; 
	int nzm1 = mp.get_mg()->get_nzone() - 1 ; 
    
	// Computation of w_shift 
	// ----------------------
	// X component
	// -----------

	et[i]->set_w_shift().set_triad(mp.get_bvect_cart()) ;
	et[i]->set_w_shift().set(0) = 0 ; 

	// Y component
	// -----------

// For the incompressible case :
	tmp = - 6  * www / a0 * ( 1 - (mp.r)*(mp.r) / (3*a0*a0) ) ; 

// For the compressible (n=1) case : 
//	Mtbl xi = M_PI * mp.r / a0 ; 
//	Mtbl sinc = sin(xi) / xi ; 	
//	 The value of sinc is set to 1 at the origin
//	for (int k=0; k<mp.get_mg()->get_np(0); k++) {
//	    for (int j=0; j<mp.get_mg()->get_nt(0); j++) {
//		sinc.set(0, k, j, 0) = 1 ; 
//	    }
//	}
//	tmp = - 4 * www / a0 * ( 1 + sinc ) ; 

	tmp.annule(nzet, nzm1) ; 
	tmp_ext = - 4 * www / mp.r ;
	tmp_ext.annule(0, nzet-1) ; 
    
	et[i]->set_w_shift().set(1) = tmp + tmp_ext ; 

	// Z component
	// -----------
	et[i]->set_w_shift().set(2) = 0 ; 

	et[i]->set_w_shift().set_triad(mp.get_bvect_spher()) ;
	et[i]->set_w_shift().std_spectral_base() ; 
	    
	// Computation of khi_shift
	// ------------------------

	tmp = 2 * www / a0 * (mp.y) * ( 1 - 3 * (mp.r)*(mp.r) / (5*a0*a0) ) ;
	tmp.annule(nzet, nzm1) ; 
	tmp_ext = 0.8 * www * a0*a0 * (mp.sint) * (mp.sinp) 
					    / (mp.r * mp.r) ;   
	tmp_ext.annule(0, nzet-1) ; 

	et[i]->set_khi_shift() = tmp + tmp_ext ; 

	// Sets the standard spectral bases for a scalar field
	et[i]->set_khi_shift().std_spectral_base() ; 	    
    
    }
*/ 
}
