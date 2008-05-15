/*
 *  Method of class Star_bhns to compute metric quantities from
 *   the companion black-hole and total metric quantities
 *
 *    (see file star_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005,2007 Keisuke Taniguchi
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

char star_bhns_upmetr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2008/05/15 19:17:39  k_taniguchi
 * Change of some parameters.
 *
 * Revision 1.1  2007/06/22 01:32:37  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
//#include <>

// Lorene headers
#include "star_bhns.h"
#include "hole_bhns.h"
#include "utilitaires.h"

void Star_bhns::update_metric_bhns(const Hole_bhns& hole,
				   const Star_bhns& star_prev,
				   double relax_met_comp) {

    //-----------------------------------------------------
    // Computation of quantities coming from the companion
    //-----------------------------------------------------

    const Map& mp_bh (hole.get_mp()) ;

    // Lapconf function
    // ----------------

    if ( (hole.get_lapconf_auto()).get_etat() == ETATZERO ) {
        lapconf_comp.set_etat_zero() ;
    }
    else {
        lapconf_comp.set_etat_qcq() ;
	lapconf_comp.import( hole.get_lapconf_auto() ) ;
	lapconf_comp.std_spectral_base() ;
    }

    // Shift vector
    // ------------

    if ( (hole.get_shift_auto())(2).get_etat() == ETATZERO ) {
        assert( (hole.get_shift_auto())(1).get_etat() == ETATZERO ) ;
	assert( (hole.get_shift_auto())(3).get_etat() == ETATZERO ) ;

	shift_comp.set_etat_zero() ;
    }
    else {
        shift_comp.set_etat_qcq() ;
	shift_comp.set_triad(mp.get_bvect_cart()) ;

	Vector comp_shift(hole.get_shift_auto()) ;
	comp_shift.change_triad(mp_bh.get_bvect_cart()) ;
	comp_shift.change_triad(mp.get_bvect_cart()) ;

	assert( *(shift_comp.get_triad()) == *(comp_shift.get_triad()) ) ;

	(shift_comp.set(1)).import( comp_shift(1) ) ;
	(shift_comp.set(2)).import( comp_shift(2) ) ;
	(shift_comp.set(3)).import( comp_shift(3) ) ;

	shift_comp.std_spectral_base() ;
    }

    // Conformal factor
    // ----------------

    if ( (hole.get_confo_auto()).get_etat() == ETATZERO ) {
        confo_comp.set_etat_zero() ;
    }
    else {
        confo_comp.set_etat_qcq() ;
	confo_comp.import( hole.get_confo_auto() ) ;
	confo_comp.std_spectral_base() ;
    }

    //----------------------------------------------------
    // Relaxation on lapconf_comp, shift_comp, confo_comp
    //----------------------------------------------------

    double relax_met_comp_jm1 = 1. - relax_met_comp ;

    lapconf_comp = relax_met_comp * lapconf_comp
      + relax_met_comp_jm1 * (star_prev.lapconf_comp) ;

    shift_comp = relax_met_comp * shift_comp
      + relax_met_comp_jm1 * (star_prev.shift_comp) ;

    confo_comp = relax_met_comp * confo_comp
      + relax_met_comp_jm1 * (star_prev.confo_comp) ;

    //-------------------------
    // Total metric quantities
    //-------------------------

    lapconf_tot = lapconf_auto + lapconf_comp ;
    lapconf_tot.std_spectral_base() ;

    shift_tot = shift_auto + shift_comp ;
    shift_tot.std_spectral_base() ;

    confo_tot = confo_auto + confo_comp ;
    confo_tot.std_spectral_base() ;

    psi4 = pow(confo_tot, 4.) ;
    psi4.std_spectral_base() ;

    lapse_auto = lapconf_auto / confo_tot ;
    lapse_auto.std_spectral_base() ;

    lapse_tot = lapconf_tot / confo_tot ;
    lapse_tot.std_spectral_base() ;

    //---------------------------------
    // Derivative of metric quantities
    //---------------------------------

    d_lapconf_auto.set_etat_qcq() ;
    for (int i=1; i<=3; i++)
        d_lapconf_auto.set(i) = lapconf_auto.deriv(i) ;

    d_lapconf_auto.std_spectral_base() ;

    d_shift_auto.set_etat_qcq() ;
    for (int i=1; i<=3; i++) {
        for (int j=1; j<=3; j++) {
	    d_shift_auto.set(i,j) = shift_auto(j).deriv(i) ;
        }
    }

    d_shift_auto.std_spectral_base() ;

    d_confo_auto.set_etat_qcq() ;
    for (int i=1; i<=3; i++)
        d_confo_auto.set(i) = confo_auto.deriv(i) ;

    d_confo_auto.std_spectral_base() ;

    // ... extrinsic curvature (taij_auto and taij_quad_auto)
    //    extr_curv_bhns() ;

    // The derived quantities are obsolete
    // -----------------------------------

    del_deriv() ;

}
