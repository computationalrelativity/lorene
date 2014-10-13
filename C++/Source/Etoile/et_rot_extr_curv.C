/*
 * Member function of class Etoile_rot to compute the extrinsic curvature
 */

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


char et_rot_extr_curv_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/13 08:52:57  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2000/11/18  17:14:35  eric
 * Traitement du cas np=1 (axisymetrie).
 *
 * Revision 2.1  2000/10/06  15:07:10  eric
 * Traitement des cas ETATZERO.
 *
 * Revision 2.0  2000/09/18  16:15:45  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "etoile.h"

namespace Lorene {
void Etoile_rot::extrinsic_curvature(){
    

	// ---------------------------------------
	// Special treatment for axisymmetric case
	// ---------------------------------------
	
 	if ( (mp.get_mg())->get_np(0) == 1) {
 	
 		tkij.set_etat_zero() ;		// initialisation
		
		
		// Computation of K_xy
		// -------------------
		
		Cmp dnpdr = nphi().dsdr() ; 		// d/dr (N^phi)
 		Cmp dnpdt = nphi().srdsdt() ; 		// 1/r d/dtheta (N^phi)
 		
 		// What follows is valid only for a mapping of class Map_radial :	
		assert( dynamic_cast<const Map_radial*>(&mp) != 0x0 ) ;
		
		if (dnpdr.get_etat() == ETATQCQ) {
		    dnpdr.va = (dnpdr.va).mult_st() ;	// multiplication by sin(theta)
 		}
		
		if (dnpdt.get_etat() == ETATQCQ) {
		    dnpdt.va = (dnpdt.va).mult_ct() ;	// multiplication by cos(theta)
 		}
 	
 	    Cmp tmp = dnpdr + dnpdt ;
 	
 	    tmp.mult_rsint() ;			// multiplication by r sin(theta)
 	
 		if (tmp.get_etat() != ETATZERO) {
 			tkij.set_etat_qcq() ;
 			tkij.set(0,1) = - 0.5 * tmp / nnn() ; 	// component (x,y)
 		}
 	
		// Computation of K_yz
		// -------------------
 	
		dnpdr = nphi().dsdr() ; 		// d/dr (N^phi)
 		dnpdt = nphi().srdsdt() ; 		// 1/r d/dtheta (N^phi)
 		
		if (dnpdr.get_etat() == ETATQCQ) {
		    dnpdr.va = (dnpdr.va).mult_ct() ;	// multiplication by cos(theta)
 		}
		
		if (dnpdt.get_etat() == ETATQCQ) {
		    dnpdt.va = (dnpdt.va).mult_st() ;	// multiplication by sin(theta)
 		}
 	
 	    tmp = dnpdr - dnpdt ;
		
 	    tmp.mult_rsint() ;			// multiplication by r sin(theta)
 		
 		if (tmp.get_etat() != ETATZERO) {
 			if (tkij.get_etat() != ETATQCQ) {
 				tkij.set_etat_qcq() ;
 			}
   			tkij.set(1,2) = - 0.5 * tmp / nnn() ; 	// component (y,z)
 	    }
 	
 	    // The other components are set to zero
 	    // ------------------------------------
  		if (tkij.get_etat() == ETATQCQ) {
			tkij.set(0,0) = 0 ;			// component (x,x)
			tkij.set(0,2) = 0 ;         // component (x,z)
			tkij.set(1,1) = 0 ;         // component (y,y)
			tkij.set(2,2) = 0 ;         // component (z,z)
		}	
 	
 	
 	
 	}
    else {

    // ------------
    // General case
    // ------------

    	// Gradient (Cartesian components) of the shift
    	// D_j N^i
    
    	Tenseur dn = shift.gradient() ;
    
    	// Trace of D_j N^i = divergence of N^i :
    	Tenseur divn = contract(dn, 0, 1) ;
    
    	if (divn.get_etat() == ETATQCQ) {
    
			// Computation of B^{-2} K_{ij}
			// ----------------------------
			tkij.set_etat_qcq() ;
			for (int i=0; i<3; i++) {
	    		for (int j=i; j<3; j++) {
					tkij.set(i, j) = dn(i, j) + dn(j, i)  ;
	    		}
	    		tkij.set(i, i) -= double(2) /double(3) * divn() ;
			}
    
			tkij = - 0.5 * tkij / nnn ;
	
    	}
    	else{
			assert( divn.get_etat() == ETATZERO ) ;
			tkij.set_etat_zero() ;
    	}
   	}
    
    // Computation of A^2 K_{ij} K^{ij}
    // --------------------------------
    
    if (tkij.get_etat() == ETATZERO) {
		ak_car = 0 ;
    }
    else {
		ak_car.set_etat_qcq() ;
    
		ak_car.set() = 0 ;
    
		for (int i=0; i<3; i++) {
	    	for (int j=0; j<3; j++) {
	
				ak_car.set() += tkij(i, j) * tkij(i, j) ;
	
	    	}
		}
    
		ak_car = b_car * ak_car ;
    }
    
}

}
