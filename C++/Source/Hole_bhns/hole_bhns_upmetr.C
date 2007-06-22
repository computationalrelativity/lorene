/*
 *  Method of class Hole_bhns to compute metric quantities from
 *   the companion neutron-star and total metric quantities
 *
 *    (see file hole_bhns.h for documentation).
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

char hole_bhns_upmetr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:25:31  k_taniguchi
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
#include "hole_bhns.h"
#include "star_bhns.h"
#include "utilitaires.h"
#include "unites.h"

void Hole_bhns::update_metric_bhns(const Star_bhns& star,
				   const Hole_bhns& hole_prev,
				   double relax) {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    //-----------------------------------------------------
    // Computation of quantities coming from the companion
    //-----------------------------------------------------

    const Map& mp_ns (star.get_mp()) ;

    double mass = ggrav * mass_bh ;

    // Lapse function
    // --------------

    if ( (star.get_lapse_auto()).get_etat() == ETATZERO ) {
        lapse_comp.set_etat_zero() ;
    }
    else {
        lapse_comp.set_etat_qcq() ;
	lapse_comp.import( star.get_lapse_auto() ) ;
	lapse_comp.std_spectral_base() ;
    }

    // Shift vector
    // ------------

    if ( (star.get_shift_auto())(2).get_etat() == ETATZERO ) {
        assert( (star.get_shift_auto())(1).get_etat() == ETATZERO ) ;
	assert( (star.get_shift_auto())(3).get_etat() == ETATZERO ) ;

        shift_comp.set_etat_zero() ;
    }
    else {
        shift_comp.set_etat_qcq() ;
	shift_comp.set_triad(mp.get_bvect_cart()) ;

	Vector comp_shift(star.get_shift_auto()) ;
	comp_shift.change_triad(mp_ns.get_bvect_cart()) ;
	comp_shift.change_triad(mp.get_bvect_cart()) ;

	assert( *(shift_comp.get_triad()) == *(comp_shift.get_triad()) ) ;

	(shift_comp.set(1)).import( comp_shift(1) ) ;
	(shift_comp.set(2)).import( comp_shift(2) ) ;
	(shift_comp.set(3)).import( comp_shift(3) ) ;

	shift_comp.std_spectral_base() ;
    }

    // Conformal factor
    // ----------------

    if ( (star.get_confo_auto()).get_etat() == ETATZERO ) {
        confo_comp.set_etat_zero() ;
    }
    else {
        confo_comp.set_etat_qcq() ;
	confo_comp.import( star.get_confo_auto() ) ;
	confo_comp.std_spectral_base() ;
    }

    //--------------------------------------------------
    // Relaxation on lapse_comp, shift_comp, confo_comp
    //--------------------------------------------------

    double relax_jm1 = 1. - relax ;

    lapse_comp = relax * lapse_comp + relax_jm1 * (hole_prev.lapse_comp) ;

    shift_comp = relax * shift_comp + relax_jm1 * (hole_prev.shift_comp) ;

    confo_comp = relax * confo_comp + relax_jm1 * (hole_prev.confo_comp) ;


    // Coordinates
    // -----------
    Scalar rr(mp) ;
    rr = mp.r ;
    rr.std_spectral_base() ;
    Scalar st(mp) ;
    st = mp.sint ;
    st.std_spectral_base() ;
    Scalar ct(mp) ;
    ct = mp.cost ;
    ct.std_spectral_base() ;
    Scalar sp(mp) ;
    sp = mp.sinp ;
    sp.std_spectral_base() ;
    Scalar cp(mp) ;
    cp = mp.cosp ;
    cp.std_spectral_base() ;

    Vector ll(mp, CON, mp.get_bvect_cart()) ;
    ll.set_etat_qcq() ;
    ll.set(1) = st * cp ;
    ll.set(2) = st * sp ;
    ll.set(3) = ct ;
    ll.std_spectral_base() ;


    if (kerrschild) {

        //----------------------------------------------
        // Metric quantities from the analytic solution
        //----------------------------------------------

	lapse_auto_bh = 1. / sqrt(1.+2.*mass/rr) ;
	lapse_auto_bh.std_spectral_base() ;
	lapse_auto_bh.annule_domain(0) ;
	lapse_auto_bh.raccord(1) ;

	confo_auto_bh = 1. ;
	confo_auto_bh.std_spectral_base() ;

	shift_auto_bh = 2. * lapse_auto_bh * lapse_auto_bh * mass * ll / rr ;
	shift_auto_bh.std_spectral_base() ;
	shift_auto_bh.annule_domain(0) ;

	//---------------------------------
	// Derivative of metric quantities
	//---------------------------------

	d_lapse_auto_rs.set_etat_qcq() ;
	for (int i=1; i<=3; i++)
	    d_lapse_auto_rs.set(i) = lapse_auto_rs.deriv(i) ;

	d_lapse_auto_rs.std_spectral_base() ;

	d_lapse_auto_bh.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	    d_lapse_auto_bh.set(i) = pow(lapse_auto_bh,3.) * mass * ll(i)
	      / rr / rr ;
	}
	d_lapse_auto_bh.std_spectral_base() ;
	d_lapse_auto_bh.annule_domain(0) ;
	d_lapse_auto_bh.inc_dzpuis(2) ;

	d_lapse_auto = d_lapse_auto_rs + d_lapse_auto_bh ;
	d_lapse_auto.std_spectral_base() ;

	d_shift_auto_rs.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        d_shift_auto_rs.set(i,j) = shift_auto_rs(j).deriv(i) ;
	    }
	}

	d_shift_auto_rs.std_spectral_base() ;

	d_shift_auto_bh.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        d_shift_auto_bh.set(i,j) = 2.*lapse_auto_bh*lapse_auto_bh*mass
		  * (flat.con()(i,j)
		     - 2.*lapse_auto_bh*lapse_auto_bh*(1.+mass/rr)
		     * ll(i) * ll(j))
		  / rr / rr ;
	    }
	}
	d_shift_auto_bh.std_spectral_base() ;
	d_shift_auto_bh.annule_domain(0) ;
	d_shift_auto_bh.inc_dzpuis(2) ;

	d_shift_auto = d_shift_auto_rs + d_shift_auto_bh ;
	d_shift_auto.std_spectral_base() ;

	d_confo_auto_rs.set_etat_qcq() ;
	for (int i=1; i<=3; i++)
	    d_confo_auto_rs.set(i) = confo_auto_rs.deriv(i) ;

	d_confo_auto_rs.std_spectral_base() ;

	d_confo_auto_bh.set_etat_qcq() ;
	for (int i=1; i<=3; i++)
	    d_confo_auto_bh.set(i) = 0. ;

	d_confo_auto_bh.std_spectral_base() ;

	d_confo_auto = d_confo_auto_rs + d_confo_auto_bh ;
	d_confo_auto.std_spectral_base() ;

    }
    else {  // Isotropic coordinates with the maximal slicing

        // Sets C/M^2 for each case of the lapse boundary condition
        // --------------------------------------------------------
        double cc ;

	if (bc_lapse_nd) {  // Neumann boundary condition
	    if (bc_lapse_fs) {  // First condition
	        // d\alpha/dr = 0
	        // --------------
	        cc = 2. ;
	    }
	    else {  // Second condition
	        // d\alpha/dr = \alpha/(2 rah)
	        // ---------------------------
	        cc = 0.5 * (sqrt(17.) - 1.) ;
	    }
	}
	else {  // Dirichlet boundary condition
	    if (bc_lapse_fs) {  // First condition
	        // \alpha = 1/2
	        // ------------
	        cc = 2. ;
	    }
	    else {  // Second condition
	        // \alpha = 1/sqrt(2.)
	        // -------------------
	        cc = 2. * sqrt(2.) ;
	    }
	}

        //----------------------------------------------
        // Metric quantities from the analytic solution
        //----------------------------------------------

	Scalar r_are(mp) ;
	r_are = r_coord(bc_lapse_nd, bc_lapse_fs) ;
	r_are.std_spectral_base() ;

	lapse_auto_bh = sqrt(1. - 2.*mass/r_are/rr
			     + cc*cc*pow(mass/r_are/rr, 4.)) ;
	lapse_auto_bh.std_spectral_base() ;
	lapse_auto_bh.annule_domain(0) ;
	lapse_auto_bh.raccord(1) ;

	confo_auto_bh = sqrt(r_are) ;
	confo_auto_bh.std_spectral_base() ;
	confo_auto_bh.annule_domain(0) ;
	confo_auto_bh.raccord(1) ;

	shift_auto_bh = mass * mass * cc * ll / rr / rr / pow(r_are, 3.) ;
	shift_auto_bh.std_spectral_base() ;
	shift_auto_bh.annule_domain(0) ;
	for (int i=1; i<=3; i++)
	    shift_auto_bh.set(i).raccord(1) ;

	//---------------------------------
	// Derivative of metric quantities
	//---------------------------------

	d_lapse_auto_rs.set_etat_qcq() ;
	for (int i=1; i<=3; i++)
	    d_lapse_auto_rs.set(i) = lapse_auto_rs.deriv(i) ;

	d_lapse_auto_rs.std_spectral_base() ;

	d_lapse_auto_bh.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	    d_lapse_auto_bh.set(i) =
	      (mass/r_are/rr - 2.*cc*cc*pow(mass/r_are/rr,4.)) * ll(i) / rr ;
	}
	d_lapse_auto_bh.std_spectral_base() ;
	d_lapse_auto_bh.annule_domain(0) ;
	d_lapse_auto_bh.inc_dzpuis(2) ;

	d_lapse_auto = d_lapse_auto_rs + d_lapse_auto_bh ;
	d_lapse_auto.std_spectral_base() ;
	d_lapse_auto.annule_domain(0) ;
	for (int i=1; i<=3; i++)
	    d_lapse_auto.set(i).raccord(1) ;

	d_shift_auto_rs.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        d_shift_auto_rs.set(i,j) = shift_auto_rs(j).deriv(i) ;
	    }
	}

	d_shift_auto_rs.std_spectral_base() ;

	d_shift_auto_bh.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        d_shift_auto_bh.set(i,j) =
		  mass*mass*cc*(flat.con()(i,j)
				-3.*sqrt(1. - 2.*mass/r_are/rr
					 +cc*cc*pow(mass/r_are/rr,4.))
				*ll(i)*ll(j))
		  / pow(r_are*rr,3.) ;
	    }
	}
	d_shift_auto_bh.std_spectral_base() ;
	d_shift_auto_bh.annule_domain(0) ;
	d_shift_auto_bh.inc_dzpuis(2) ;

	d_shift_auto = d_shift_auto_rs + d_shift_auto_bh ;
	d_shift_auto.std_spectral_base() ;
	d_shift_auto.annule_domain(0) ;
	for (int i=1; i<=3; i++) {
	    for (int j=1; j<=3; j++) {
	        d_shift_auto.set(i,j).raccord(1) ;
	    }
	}

	d_confo_auto_rs.set_etat_qcq() ;
	for (int i=1; i<=3; i++)
	    d_confo_auto_rs.set(i) = confo_auto_rs.deriv(i) ;

	d_confo_auto_rs.std_spectral_base() ;

	d_confo_auto_bh.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	    d_confo_auto_bh.set(i) = 0.5*sqrt(r_are)
	      *(sqrt(1. - 2.*mass/r_are/rr +cc*cc*pow(mass/r_are/rr,4.)) - 1.)
	      * ll(i) / rr ;
	}
	d_confo_auto_bh.std_spectral_base() ;
	d_confo_auto_bh.annule_domain(0) ;
	d_confo_auto_bh.inc_dzpuis(2) ;

	d_confo_auto = d_confo_auto_rs + d_confo_auto_bh ;
	d_confo_auto.std_spectral_base() ;
	d_confo_auto.annule_domain(0) ;
	for (int i=1; i<=3; i++)
	    d_confo_auto.set(i).raccord(1) ;

    }

    //-------------------------
    // Total metric quantities
    //-------------------------

    lapse_auto = lapse_auto_rs + lapse_auto_bh ;
    lapse_auto.std_spectral_base() ;
    lapse_auto.annule_domain(0) ;
    lapse_auto.raccord(1) ;

    lapse_tot = lapse_auto_rs + lapse_auto_bh + lapse_comp ;
    lapse_tot.std_spectral_base() ;
    lapse_tot.annule_domain(0) ;
    lapse_tot.raccord(1) ;

    shift_auto = shift_auto_rs + shift_auto_bh ;
    shift_auto.std_spectral_base() ;
    shift_auto.annule_domain(0) ;
    for (int i=1; i<=3; i++) {
        shift_auto.set(i).raccord(1) ;
    }

    shift_tot = shift_auto_rs + shift_auto_bh + shift_comp ;
    shift_tot.std_spectral_base() ;
    shift_tot.annule_domain(0) ;
    for (int i=1; i<=3; i++) {
        shift_tot.set(i).raccord(1) ;
    }

    confo_auto = confo_auto_rs + confo_auto_bh ;
    confo_auto.std_spectral_base() ;
    confo_auto.annule_domain(0) ;
    confo_auto.raccord(1) ;

    confo_tot = confo_auto_rs + confo_auto_bh + confo_comp ;
    confo_tot.std_spectral_base() ;
    confo_tot.annule_domain(0) ;
    confo_tot.raccord(1) ;

    // The derived quantities are obsolete
    // -----------------------------------

    del_deriv() ;

}
