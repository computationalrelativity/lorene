/*
 *  Methods Et_bin_nsbh::update_metric
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:17:53  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:52:56  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2007/04/24 20:14:45  f_limousin
 * Implementation of Dirichlet and Neumann BC for the lapse
 *
 * Revision 1.2  2005/08/29 15:10:17  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.1  2003/10/24 12:28:39  k_taniguchi
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

		    //----------------------------------//
		    //	 Version a BH companion 	//
		    //----------------------------------//

namespace Lorene {
void Et_bin_nsbh::update_metric(const Bhole& comp) {

    // Computation of quantities coming from the companion
    // ---------------------------------------------------

    // Computes N_comp  ---> stored in logn_comp
    if ( (comp.get_n_auto()).get_etat() == ETATZERO ) {
	n_comp.set_etat_zero() ;
    }
    else{
	n_comp.set_etat_qcq() ;
	(n_comp.set()).import( comp.get_n_auto()() ) ;
	n_comp.set_std_base() ;   // set the bases for spectral expansions
    }


    // Computes Psi_comp  ---> stored in beta_comp
    if ( (comp.get_psi_auto()).get_etat() == ETATZERO ) {
	confpsi_comp.set_etat_zero() ;
    }
    else{
	confpsi_comp.set_etat_qcq() ;
	(confpsi_comp.set()).import( comp.get_psi_auto()() ) ;
	confpsi_comp.set_std_base() ; // set the bases for spectral expansions
    }

 
    // Computes N^i_comp  ---> stored in shift_comp
    if ( (comp.get_shift_auto()).get_etat() == ETATZERO ) {
	shift_comp.set_etat_zero() ;
    }
    else{
	shift_comp.set_etat_qcq() ;

	(shift_comp.set(0)).import( comp.get_shift_auto()(0) ) ;  // N^x antisym
	(shift_comp.set(1)).import( comp.get_shift_auto()(1) ) ;   // N^y sym.
	(shift_comp.set(2)).import( comp.get_shift_auto()(2) ) ;  // N^z anisym

	shift_comp.set_std_base() ;   // set the bases for spectral expansions
    }

    // Lapse function N
    // ----------------

    nnn = n_auto + n_comp ;

    // Conformal factor confpsi
    // ------------------------

    confpsi = confpsi_auto + confpsi_comp ;

    confpsi.set_std_base() ;   // set the bases for spectral expansions
    
    // Conformal factor A^2
    // ---------------------
    
    a_car =  pow(confpsi, 4.);
    a_car.set_std_base() ;   // set the bases for spectral expansions
    
    // Shift vector N^i
    // ----------------

    shift = shift_auto + shift_comp ;

    // Derivatives of metric coefficients
    // ----------------------------------

    // ... (d/dX,d/dY,d/dZ)(n_auto) :
    d_n_auto = n_auto.gradient() ;    // (d/dx, d/dy, d/dz)
    d_n_auto.change_triad(ref_triad) ;   // --> (d/dX, d/dY, d/dZ)

    // ... (d/dX,d/dY,d/dZ)(confpsi_auto) :
    d_confpsi_auto = confpsi_auto.gradient() ;    // (d/dx, d/dy, d/dz)
    d_confpsi_auto.change_triad(ref_triad) ;   // --> (d/dX, d/dY, d/dZ)

    // The derived quantities are obsolete
    // -----------------------------------

    del_deriv() ;
}
}
