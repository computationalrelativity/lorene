/*
 * Method of class Etoile_bin to compute the extrinsic curvature tensor
 *
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


char et_bin_extr_curv_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  2000/03/07  14:51:49  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "etoile.h"

void Etoile_bin::extrinsic_curvature(){
    
    // Components of shift_auto with respect to the Cartesian triad
    //  (d/dx, d/dy, d/dz) of the mapping : 
    Tenseur shift_auto_local = shift_auto ; 
    shift_auto_local.change_triad( mp.get_bvect_cart() ) ; 
    
    // Gradient (partial derivatives with respect to the Cartesian coordinates
    //           of the mapping)
    // D_j N^i 
    
    Tenseur dn = shift_auto_local.gradient() ; 
    
    // Return to the absolute reference frame
    dn.change_triad(ref_triad) ; 
    
    // Trace of D_j N^i = divergence of N^i : 
    Tenseur divn = contract(dn, 0, 1) ; 
    
    // Computation of A^2 K^{ij}
    // -------------------------
    tkij_auto.set_etat_qcq() ; 
    for (int i=0; i<3; i++) {
	for (int j=i; j<3; j++) {
	    tkij_auto.set(i, j) = dn(i, j) + dn(j, i)  ; 
	}
	tkij_auto.set(i, i) -= double(2) /double(3) * divn() ; 
    }
    
    tkij_auto = - 0.5 * tkij_auto / nnn ; 
    
    // Computation of A^2 K_{ij} K^{ij}
    // --------------------------------
    
    akcar_auto.set_etat_qcq() ; 
    
    akcar_auto.set() = 0 ; 
    
    for (int i=0; i<3; i++) {
	for (int j=0; j<3; j++) {
	
	    akcar_auto.set() += tkij_auto(i, j) * tkij_auto(i, j) ; 
	
	}
    }
    
    akcar_auto = a_car * akcar_auto ; 
    
    
}
