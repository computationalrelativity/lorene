/*
 *  Methods Et_bin_nsbh::update_metric_der_comp
 *
 *    (see file et_bin_nsbh.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003 Philippe Grandclement
 *                 2003 Keisuke Taniguchi
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

char et_bin_nsbh_upmetr_der_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2014/10/13 08:52:56  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2007/04/24 20:14:45  f_limousin
 * Implementation of Dirichlet and Neumann BC for the lapse
 *
 * Revision 1.4  2005/10/18 13:12:33  p_grandclement
 * update of the mixted binary codes
 *
 * Revision 1.3  2005/08/29 15:10:17  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.2  2004/06/07 11:08:31  k_taniguchi
 * A minor change.
 *
 * Revision 1.1  2003/10/24 12:29:21  k_taniguchi
 * Method of update metric for the BH companion
 *
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "et_bin_nsbh.h"
#include "bhole.h"

        //---------------------------------------//
	//      Version with BH companion       //
        //---------------------------------------//


namespace Lorene {
void Et_bin_nsbh::update_metric_der_comp(const Bhole& comp) {

    // Computation of Grad(N) ---> stored in d_n_comp
    // -------------------------------------------------

    Tenseur dncomp = ( comp.get_n_auto() ).gradient() ;

    if ( dncomp.get_etat() == ETATZERO ) {
	d_n_comp.set_etat_zero() ;
    }
    else{

	// 1/ Division by r^2 of comp.d_n_auto in the ZEC
	dncomp.dec2_dzpuis() ;

	// 2/ Interpolation of the result

	d_n_comp.set_etat_qcq() ;
	(d_n_comp.set(0)).import( dncomp(0) ) ;  // d/dx sym.
	(d_n_comp.set(1)).import( dncomp(1) ) ; // d/dy antisym.
	(d_n_comp.set(2)).import( dncomp(2) ) ;  // d/dz sym.

    }
    d_n_comp.set_std_base() ;
    d_n_comp.inc2_dzpuis() ;
    d_n_comp.set_triad( *(dncomp.get_triad()) ) ;
    d_n_comp.change_triad(ref_triad) ;
    

    // Computation of Grad(Psi) ---> stored in d_confpsi_comp
    // ------------------------------------------------------

    Tenseur dpsicomp = ( comp.get_psi_auto() ).gradient() ;

    if ( dpsicomp.get_etat() == ETATZERO ) {
	d_confpsi_comp.set_etat_zero() ;
    }
    else {
	// 1/ Division by r^2 of comp.d_confpsi_auto in the ZEC
        dpsicomp.dec2_dzpuis() ;

	// 2/ Interpolation of the result

	d_confpsi_comp.set_etat_qcq() ;

	(d_confpsi_comp.set(0)).import(dpsicomp(0) ) ;  // d/dx sym.
	(d_confpsi_comp.set(1)).import(dpsicomp(1) ) ; // d/dy antisym.
	(d_confpsi_comp.set(2)).import(dpsicomp(2) ) ;  // d/dz sym.

    }

    d_confpsi_comp.set_std_base() ;
    d_confpsi_comp.inc2_dzpuis() ;
    d_confpsi_comp.set_triad( *(dpsicomp.get_triad()) ) ;
    d_confpsi_comp.change_triad(ref_triad) ;

    // The derived quantities are obsolete
    // -----------------------------------
    del_deriv() ;
}
}
