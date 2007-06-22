/*
 *  Method of class Star_bhns to compute the extrinsic curvature tensor
 *
 *    (see file star_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005-2006 Keisuke Taniguchi
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

char star_bhns_extr_curv_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:31:05  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
#include <math.h>

// Lorene headers
#include "star_bhns.h"

void Star_bhns::extr_curv_bhns() {

    // Computation of \tilde{A}_{NS}^{ij}
    // ----------------------------------

    Scalar divshift(mp) ;
    divshift = d_shift_auto(1,1) + d_shift_auto(2,2)
      + d_shift_auto(3,3) ;
    divshift.std_spectral_base() ;

    Sym_tensor flat_taij(mp, CON, mp.get_bvect_cart()) ;
    flat_taij.set_etat_qcq() ;

    for (int i=1; i<=3; i++) {
        for (int j=1; j<=3; j++) {
	    flat_taij.set(i,j) = d_shift_auto(i,j)
	      + d_shift_auto(j,i)
	      - 2. * divshift * flat.con()(i,j) / 3. ;
	}
    }
    flat_taij.std_spectral_base() ;

    taij_auto = 0.5 * pow(confo_auto+0.5, 6.) * flat_taij
      / (lapse_auto+0.5) ;
    taij_auto.std_spectral_base() ;


    // Computation of \tilde{A}_{NS}^{ij} \tilde{A}^{NS}_{ij}
    // ------------------------------------------------------

    Sym_tensor flat_dshift(mp, COV, mp.get_bvect_cart()) ;
    flat_dshift.set_etat_qcq() ;

    for (int i=1; i<=3; i++) {
        for (int j=1; j<=3; j++) {
	    flat_dshift.set(i,j) =
	      flat.cov()(j,1) * d_shift_auto(i,1)
	      + flat.cov()(j,2) * d_shift_auto(i,2)
	      + flat.cov()(j,3) * d_shift_auto(i,3)
	      + flat.cov()(i,1) * d_shift_auto(j,1)
	      + flat.cov()(i,2) * d_shift_auto(j,2)
	      + flat.cov()(i,3) * d_shift_auto(j,3)
	      - 2. * divshift * flat.cov()(i,j) / 3. ;
	}
    }
    flat_dshift.std_spectral_base() ;

    Sym_tensor taij_down(mp, COV, mp.get_bvect_cart()) ;
    taij_down.set_etat_qcq() ;

    taij_down = 0.5 * pow(confo_auto+0.5, 6.) * flat_dshift
      / (lapse_auto+0.5) ;
    taij_down.std_spectral_base() ;

    taij_quad_auto = 0. ;

    for (int i=1; i<=3; i++) {
        for (int j=1; j<=3; j++) {
	    taij_quad_auto += taij_down(i,j) * taij_auto(i,j) ;
	}
    }

    taij_quad_auto.std_spectral_base() ;

}
