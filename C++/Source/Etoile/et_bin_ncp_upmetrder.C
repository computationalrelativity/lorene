/*
 * Methods Et_bin_ncp::update_metric_der_comp
 *
 * (see file et_bin_ncp.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003 Francois Limousin
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


char et_bin_ncp_upmetrder_C[] = "$Header$" ;

/*
 * $Header$
 *
 */

// Headers Lorene
#include "et_bin_ncp.h"

void Et_bin_ncp::update_metric_der_comp(const Et_bin_ncp& comp) {

    // Computation of d_logn_comp
    // --------------------------

    if ( (comp.d_logn_auto).get_etat() == ETATZERO ) {
	d_logn_comp.set_etat_zero() ;
    }
    else{
      d_logn_comp = logn_comp.gradient() ;
    }

    d_logn_comp.change_triad(mp.get_bvect_cart()) ;

 
    // Computation of tkij_comp
    // ------------------------

    if ( (comp.tkij_auto).get_etat() == ETATZERO ) {
	tkij_comp.set_etat_zero() ;
    }
    else{

      // Components of shift_comp with respect to the Cartesian triad
      //  (d/dx, d/dy, d/dz) of the mapping :
      Tenseur shift_comp_local = shift_comp ;
 
      // Gradient tilde (partial derivatives with respect to
      //           the Cartesian coordinates of the mapping)
      // D~_j beta^i

      Tenseur dshift_comp = shift_comp_local.gradient() ;

      // Trace of D~_j beta^i  :
      Tenseur divshift_comp = contract(shift_comp_local.derive_cov(gtilde), 0, 1) ;

      // Computation of K^{ij}
      // -------------------------
      tkij_comp.set_etat_qcq() ;

      for (int i=0; i<3; i++) {
	for (int j=i; j<3; j++) {
	    tkij_comp.set(i, j) = dshift_comp(i, j) + dshift_comp(j, i) - double(2) /double(3) * divshift_comp() * (gtilde.con())(i,j) ; 

      }

      tkij_comp = - 0.5 * tkij_comp / nnn ;

    }

    tkij_comp.set_std_base() ;

    if (relativistic) {
	// Computation of kcar_comp
	// -------------------------
    
	kcar_comp.set_etat_qcq() ;
    
	kcar_comp.set() = 0 ;

	for (int i=0; i<3; i++) {
	    for (int j=0; j<3; j++) {

	    kcar_comp.set() +=  contract(met_gamma.cov() * contract(met_gamma.cov() * tkij_comp, 1, 3), 1, 3)() ; 

	    }
	}

	kcar_comp.set_std_base() ;

    }

    // The derived quantities are obsolete
    // -----------------------------------

    del_deriv() ;


    }
}
