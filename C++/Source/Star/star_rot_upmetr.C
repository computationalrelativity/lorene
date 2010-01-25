/*
 * Methods Star_rot::update_metric and Star_rot::extrinsic_curvature
 *
 * (see file star_rot.h for documentation)
 *
 */

/*
 *   Copyright (c) 2010 Eric Gourgoulhon
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


char star_rot_upmetr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2010/01/25 18:15:52  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "star_rot.h"


void Star_rot::update_metric() {
 
    // Lapse function N
    // ----------------
    
    nn = exp( unsurc2 * logn ) ; 

    nn.std_spectral_base() ;   // set the bases for spectral expansions
    
    
    // Metric factor A^2
    // -----------------
    
    a_car = exp( 2*unsurc2*( dzeta - logn ) ) ; 

    a_car.std_spectral_base() ;   // set the bases for spectral expansions

    // Metric factor B
    // ---------------
    
    Scalar tmp = tggg ; 
    tmp.div_rsint() ;	        //... Division of tG by r sin(theta)

    bbb = (1 + tmp) / nn ; 

    bbb.std_spectral_base() ;   // set the bases for spectral expansions
        
    b_car = bbb * bbb ; 
    
    // Full 3-metric
    // -------------

    Sym_tensor gam(mp, COV, mp.get_bvect_spher()) ; 
    gam.set(1,1) = a_car ; 
    gam.set(1,2) = 0 ; 
    gam.set(1,3) = 0 ; 
    gam.set(2,2) = a_car ; 
    gam.set(2,3) = 0 ; 
    gam.set(3,3) = b_car ;

    gamma = gam ;

    // Tensor B^{-2} K_{ij} and Scalar A^2 K_{ij} K^{ij}
    // -------------------------------------------------
    
    extrinsic_curvature() ; 
    
  
    // The derived quantities are no longer up to date : 
    // -----------------------------------------------

    del_deriv() ;  

}


/*************************************************************************************/


void Star_rot::extrinsic_curvature(){
    

	// ---------------------------------------
	// Special treatment for axisymmetric case
	// ---------------------------------------
	
 	if ( (mp.get_mg())->get_np(0) == 1) {
 	
 		tkij.set_etat_zero() ;		// initialisation
				
		// Computation of K_xy
		// -------------------
		
		Scalar dnpdr = nphi.dsdr() ; 		// d/dr (N^phi)
 		Scalar dnpdt = nphi.srdsdt() ; 		// 1/r d/dtheta (N^phi)
 		
 		// What follows is valid only for a mapping of class Map_radial :	
		assert( dynamic_cast<const Map_radial*>(&mp) != 0x0 ) ;
		
		if (dnpdr.get_etat() == ETATQCQ) {
		    // multiplication by sin(theta)
		    dnpdr.set_spectral_va() = (dnpdr.get_spectral_va()).mult_st() ;	
 		}
		
		if (dnpdt.get_etat() == ETATQCQ) {
		    // multiplication by cos(theta)
		    dnpdt.set_spectral_va() = (dnpdt.get_spectral_va()).mult_ct() ;	
 		}
 	
		Scalar tmp = dnpdr + dnpdt ;
 	
		tmp.mult_rsint() ;	// multiplication by r sin(theta)
 	
		tkij.set(1,2) = - 0.5 * tmp / nn ; 	// component (x,y)
 	
 	
		// Computation of K_yz
		// -------------------
 	
		dnpdr = nphi.dsdr() ; 		// d/dr (N^phi)
 		dnpdt = nphi.srdsdt() ; 		// 1/r d/dtheta (N^phi)
 		
		if (dnpdr.get_etat() == ETATQCQ) {
		    // multiplication by cos(theta)
		    dnpdr.set_spectral_va() = (dnpdr.get_spectral_va()).mult_ct() ;	
 		}
		
		if (dnpdt.get_etat() == ETATQCQ) {
		    // multiplication by sin(theta)
		    dnpdt.set_spectral_va() = (dnpdt.get_spectral_va()).mult_st() ;	
 		}
 	
		tmp = dnpdr - dnpdt ;
		
		tmp.mult_rsint() ;	// multiplication by r sin(theta)
 		
		tkij.set(2,3) = - 0.5 * tmp / nn ; 	// component (y,z)
 	
		// The other components are set to zero
		// ------------------------------------
		tkij.set(1,1) = 0 ;	// component (x,x)
		tkij.set(1,3) = 0 ;     // component (x,z)
		tkij.set(2,2) = 0 ;    	// component (y,y)
		tkij.set(3,3) = 0 ;     // component (z,z)
 	
 	}
    else {

    // ------------
    // General case
    // ------------

    	// Gradient (Cartesian components) of the shift
    	// D_j N^i
    
    	Tensor dn = - beta.derive_cov( mp.flat_met_cart() ) ;
    
    	// Trace of D_j N^i = divergence of N^i :
    	Scalar divn = contract(dn, 0, 1) ;
    
    	if (divn.get_etat() == ETATQCQ) {
    
		// Computation of B^{-2} K_{ij}
		// ----------------------------
		tkij.set_etat_qcq() ;
		for (int i=1; i<=3; i++) {
		    for (int j=i; j<=3; j++) {
			  tkij.set(i, j) = dn(i, j) + dn(j, i)  ;
		    }
		    tkij.set(i, i) -= double(2) /double(3) * divn ;
		}
    
		tkij = - 0.5 * tkij / nn ;
	
    	}
    	else{
		assert( divn.get_etat() == ETATZERO ) ;
		tkij.set_etat_zero() ;
    	}
   }
    
    // Computation of A^2 K_{ij} K^{ij}
    // --------------------------------
        
    ak_car = 0 ;
    
    for (int i=1; i<=3; i++) {
	for (int j=1; j<=3; j++) {
	
	    ak_car += tkij(i, j) * tkij(i, j) ;
	
	}
    }
    
    ak_car = b_car * ak_car ;
    
}

