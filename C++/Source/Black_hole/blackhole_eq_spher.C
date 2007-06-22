/*
 *  Method of class Black_hole to compute a single black hole
 *
 *    (see file blackhole.h for documentation).
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

char blackhole_eq_spher_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:19:11  k_taniguchi
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
#include "blackhole.h"
#include "cmp.h"
#include "tenseur.h"
#include "param.h"
#include "unites.h"
#include "proto.h"
#include "utilitaires.h"
#include "graphique.h"

void Black_hole::equilibrium_spher(bool neumann, bool first, double precis) {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    // Initializations
    // ---------------

    const Mg3d* mg = mp.get_mg() ;
    int nz = mg->get_nzone() ;          // total number of domains

    double mass = ggrav * mass_bh ;

    // Inner boundary condition
    // ------------------------

    Valeur bc_laps(mg->get_angu()) ;
    Valeur bc_conf(mg->get_angu()) ;

    Valeur bc_shif_x(mg->get_angu()) ;
    Valeur bc_shif_y(mg->get_angu()) ;
    Valeur bc_shif_z(mg->get_angu()) ;

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

    // Sets C/M^2 for each case of the lapse boundary condition
    // --------------------------------------------------------
    double cc ;

    if (neumann) {  // Neumann boundary condition
        if (first) {  // First condition
	  // d\alpha/dr = 0
	  // --------------
	  cc = 2. ;
	}
	else {  // Second condition
	  // d\alpha/dr = \alpha/(2 rah)
	  // ---------------------------
	  cc = 0.5 * (sqrt(17.) - 1.) ;

	  //	  cc = 0.5 * (sqrt(17.) + 1.) ;
	}
    }
    else {  // Dirichlet boundary condition
       if (first) {  // First condition
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


    // Ininital metric
    if (kerrschild) {

        update_metric() ;

        lapse = lapse_bh ;
	lapse.std_spectral_base() ;
	lapse_rs = 0. ;
	lapse_rs.std_spectral_base() ;

	confo = 1. ;
	confo.std_spectral_base() ;

	shift = shift_bh ;
	shift.std_spectral_base() ;
	shift_rs.set_etat_zero() ;

    }
    else {  // Isotropic coordinates

        Scalar r_are(mp) ;
	r_are = r_coord(neumann, first) ;
	r_are.std_spectral_base() ;
	//	r_are.annule_domain(0) ;
	//	r_are.raccord(1) ;

	cout << "r_are:" << endl ;
	for (int l=0; l<nz; l++) {
	  cout << r_are.val_grid_point(l,0,0,0) << endl ;
	}

        lapse = sqrt(1. - 2.*mass/r_are/rr
		     + cc*cc*pow(mass/r_are/rr,4.)) ;
	lapse.std_spectral_base() ;
	lapse.annule_domain(0) ;
	lapse.raccord(1) ;

	confo = sqrt(r_are) ;
	confo.std_spectral_base() ;
	//	confo.annule_domain(0) ;
	//	confo.raccord(1) ;

	for (int i=1; i<=3; i++) {
	    shift.set(i) = mass * mass * cc * ll(i) / rr / rr
	      / pow(r_are,3.) ;
	}

	shift.std_spectral_base() ;

	for (int i=1; i<=3; i++) {
	    shift.set(i).annule_domain(0) ;
	    shift.set(i).raccord(1) ;
	}

	/*
	des_profile( r_are, 0., 20, M_PI/2., 0.,
		     "Areal coordinate",
		     "Areal (theta=pi/2, phi=0)" ) ;

	des_profile( lapse, 0., 20, M_PI/2., 0.,
		     "Initial lapse function of BH",
		     "Lapse (theta=pi/2, phi=0)" ) ;

	des_profile( confo, 0., 20, M_PI/2., 0.,
		     "Initial conformal factor of BH",
		     "Confo (theta=pi/2, phi=0)" ) ;

	des_profile( shift(1), 0., 20, M_PI/2., 0.,
		     "Initial shift vector (X) of BH",
		     "Shift (theta=pi/2, phi=0)" ) ;

	des_coupe_vect_z( shift, 0., -3., 0.5, 4,
		       "Shift vector of BH") ;
	*/
    }

    // Auxiliary quantities
    // --------------------

    Scalar source_lapse(mp) ;
    Scalar source_confo(mp) ;
    Vector source_shift(mp, CON, mp.get_bvect_cart()) ;

    Scalar lapse_m1(mp) ;  // = lapse - 1 (only for the isotropic case)
    Scalar confo_m1(mp) ;  // = confo - 1

    Scalar dlapse(mp) ;
    Scalar dlconfo(mp) ;
    Scalar lconfo(mp) ;
    Scalar lappsi(mp) ;
    Scalar confo6(mp) ;

    Scalar lapse_jm1(mp) ;
    Scalar confo_jm1(mp) ;
    Vector shift_jm1(mp, CON, mp.get_bvect_cart()) ;

    double diff_lp = 1. ;
    double diff_cf = 1. ;
    double diff_sx = 1. ;
    double diff_sy = 1. ;
    double diff_sz = 1. ;

    int mermax = 200 ;          // max number of iterations

    //======================================//
    //          Start of iteration          //
    //======================================//
    /*
    for (int mer=0;
	 (diff_lp > precis) || (diff_cf > precis) && (mer < mermax); mer++) {

    for (int mer=0;
	 (diff_sx > precis) || (diff_sy > precis) || (diff_sz > precis)
	   && (mer < mermax); mer++) {
    */
    for (int mer=0;
	 (diff_lp > precis) && (mer < mermax); mer++) {

        cout << "--------------------------------------------------" << endl ;
	cout << "step: " << mer << endl ;
	cout << "diff_lapse = " << diff_lp << endl ;
	cout << "diff_confo = " << diff_cf << endl ;
	cout << "diff_shift : x = " << diff_sx
	     << "  y = " << diff_sy << "  z = " << diff_sz << endl ;

	if (kerrschild) {
	    lapse_jm1 = lapse_rs ;
	    confo_jm1 = confo ;
	    shift_jm1 = shift_rs ;
	}
	else {
	    lapse_jm1 = lapse ;
	    confo_jm1 = confo ;
	    shift_jm1 = shift ;
	}

	//------------------------------------------//
	//  Computation of the extrinsic curvature  //
	//------------------------------------------//

	if (kerrschild) {
	    update_metric() ;
	}

	extr_curv_bh() ;

	//-------------------------------------------------------------//
	//  Resolution of the Poisson equation for the lapse function  //
	//-------------------------------------------------------------//

	// Source term
	// -----------

	if (kerrschild) {

	    lconfo = log(confo_jm1) ;
	    lconfo.std_spectral_base() ;

	    Scalar dlcon = lconfo.dsdr() ;

	    Scalar dlap = lapse_jm1.dsdr() ;

	    Scalar ddlap = dlap.dsdr() ;

	    Vector dlps(mp, COV, mp.get_bvect_cart()) ;
	    dlps.set_etat_qcq() ;
	    for (int i=1; i<=3; i++)
	        dlps.set(i) = lapse_jm1.deriv(i) ;

	    dlps.std_spectral_base() ;

	    Vector dlco(mp, CON, mp.get_bvect_cart()) ;
	    dlco.set_etat_qcq() ;
	    for (int i=1; i<=3; i++) {
	        dlco.set(i) = flat.con()(i,1)%(lconfo.deriv(1))
		  + flat.con()(i,2)%(lconfo.deriv(2))
		  + flat.con()(i,3)%(lconfo.deriv(3))
		  - 2. * mass * lapse_bh % lapse_bh % ll(i) % dlcon / rr ;
	    }
	    dlco.std_spectral_base() ;

	    Vector dtrk(mp, COV, mp.get_bvect_cart()) ;
	    dtrk.set_etat_qcq() ;
	    for (int i=1; i<=3; i++)
	        dtrk.set(i) = -2.*pow(lapse_bh,5.)*mass
		  *(2.+10.*mass/rr+9.*mass*mass/rr/rr)/pow(rr,3.)*ll(i) ;

	    dtrk.std_spectral_base() ;

	    Scalar tmp1 = lapse % taij_quad_rs / pow(confo_jm1, 8.) ;
	    tmp1.std_spectral_base() ;
	    tmp1.inc_dzpuis(4-tmp1.get_dzpuis()) ;

	    Scalar tmp2 = -2. * (dlps(1)%dlco(1) + dlps(2)%dlco(2)
				 + dlps(3)%dlco(3)) ;
	    tmp2.std_spectral_base() ;
	    tmp2.inc_dzpuis(4-tmp2.get_dzpuis()) ;

	    Scalar confo4(mp) ;
	    confo4 = pow(confo_jm1,4.) ;
	    confo4.std_spectral_base() ;

	    Scalar tmp3 = 4.*lapse_jm1%confo4*pow(lapse_bh,6.)
	      *mass*mass*(1.+3.*mass/rr)*(1.+3.*mass/rr)/3./pow(rr,4.) ;
	    tmp3.std_spectral_base() ;
	    tmp3.inc_dzpuis(4-tmp3.get_dzpuis()) ;

	    Scalar tmp4 = confo4 % (shift_jm1(1)%dtrk(1)
				    + shift_jm1(2)%dtrk(2)
				    + shift_jm1(3)%dtrk(3)) ;
	    tmp4.std_spectral_base() ;
	    tmp4.inc_dzpuis(4-tmp4.get_dzpuis()) ;

	    Scalar tmp5 = 2. * mass * lapse_bh % lapse_bh % ddlap / rr ;
	    tmp5.annule_domain(0) ;
	    tmp5.std_spectral_base() ;
	    tmp5.inc_dzpuis(4-tmp5.get_dzpuis()) ;

	    Scalar tmp6 = mass * pow(lapse_bh,4.) * dlap
	      * (3.+8.*mass/rr) / rr / rr ;
	    tmp6.annule_domain(0) ;
	    tmp6.std_spectral_base() ;
	    tmp6.inc_dzpuis(4-tmp6.get_dzpuis()) ;

	    Scalar tmp7 = -2. * pow(lapse_bh,5.) * mass * dlcon / rr / rr ;
	    tmp7.annule_domain(0) ;
	    tmp7.std_spectral_base() ;
	    tmp7.inc_dzpuis(4-tmp7.get_dzpuis()) ;

	    Scalar tmp8 = 4.*pow(lapse_bh,7.)*mass*mass/3./pow(rr,4.)
	      *(2.*(lapse_bh/lapse - 1.)*(4.+12.*mass/rr+9.*mass*mass/rr/rr)
		*confo4 + 3.*(confo4 - 1.)) ;
	    tmp8.annule_domain(0) ;
	    tmp8.std_spectral_base() ;
	    tmp8.inc_dzpuis(4-tmp8.get_dzpuis()) ;

	    source_lapse = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7
	      + tmp8 ;

	    /*
	    Scalar lap_bh7(mp) ;
	    lap_bh7 = 1./pow(1.+2.*mass/rr,3.5) ;
	    lap_bh7.std_spectral_base() ;

	    Scalar tmp1_exact(mp) ;
	    tmp1_exact = 8.*lap_bh7*mass*mass*(2.+3.*mass/rr)*(2.+3.*mass/rr)
	      /3./pow(rr,4.) ;
	    tmp1_exact.annule_domain(0) ;
	    tmp1_exact.std_spectral_base() ;
	    tmp1_exact.inc_dzpuis(4) ;

	    Scalar tmp3_exact(mp) ;
	    tmp3_exact = 4.*lap_bh7*mass*mass*(1.+3.*mass/rr)*(1.+3.*mass/rr)
	      /3./pow(rr,4.) ;
	    tmp3_exact.annule_domain(0) ;
	    tmp3_exact.std_spectral_base() ;
	    tmp3_exact.inc_dzpuis(4) ;

	    Scalar tmp4_exact(mp) ;
	    tmp4_exact = -4.*lap_bh7*mass*mass
	      *(2.+10*mass/rr+9.*mass*mass/rr/rr) / pow(rr,4.) ;
	    tmp4_exact.annule_domain(0) ;
	    tmp4_exact.std_spectral_base() ;
	    tmp4_exact.inc_dzpuis(4) ;

	    Scalar tmp5_exact(mp) ;
	    tmp5_exact = -2.*lap_bh7*mass*mass*(2.+mass/rr)/pow(rr,4.) ;
	    tmp5_exact.annule_domain(0) ;
	    tmp5_exact.std_spectral_base() ;
	    tmp5_exact.inc_dzpuis(4) ;

	    Scalar tmp6_exact(mp) ;
	    tmp6_exact = lap_bh7*mass*mass*(3.+8.*mass/rr)/pow(rr,4.) ;
	    tmp6_exact.annule_domain(0) ;
	    tmp6_exact.std_spectral_base() ;
	    tmp6_exact.inc_dzpuis(4) ;
	    */

	    /*
	    cout << "source_lapse - exact" << endl ;
	    cout << source_lapse - tmp1_exact - tmp3_exact - tmp4_exact
	      - tmp5_exact - tmp6_exact << endl ;
	    arrete() ;
	    */

	    /*
	    Scalar sou_lap_exact(mp) ;
	    sou_lap_exact = 3. * pow(lap_bh2,2.5) * mass * mass
	      / pow(rr, 4.) ;
	    sou_lap_exact.annule_domain(0) ;
	    sou_lap_exact.std_spectral_base() ;
	    sou_lap_exact.inc_dzpuis(4-sou_lap_exact.get_dzpuis()) ;

	    source_lapse = sou_lap_exact ;
	    */
	    /*
	    cout << "sou_lap_exact - exact" << endl ;
	    cout << sou_lap_exact - tmp1_exact - tmp3_exact - tmp4_exact
	      - tmp5_exact - tmp6_exact << endl ;
	    arrete() ;
	    */
	    /*
	    cout << "resi_lap" << endl ;
	    cout << source_lapse - sou_lap_exact << endl ;
	    arrete() ;
	    cout << "source_lapse" << endl ;
	    cout << source_lapse << endl ;
	    arrete() ;
	    cout << "sou_lap_exact" << endl ;
	    cout << sou_lap_exact << endl ;
	    arrete() ;
	    */
	}
	else {  // Isotropic coordinates with the maximal slicing

	    lconfo = log(confo_jm1) ;
	    lconfo.std_spectral_base() ;
	    dlapse = lapse_jm1.dsdr() ;
	    dlconfo = lconfo.dsdr() ;

	    Scalar tmp1 = lapse_jm1 % taij_quad / pow(confo_jm1, 8.) ;
	    tmp1.std_spectral_base() ;

	    Scalar tmp2 = - 2. * dlapse % dlconfo ;
	    tmp2.std_spectral_base() ;
	    tmp2.inc_dzpuis(4-tmp2.get_dzpuis()) ;

	    /*
	    Scalar tmp_taij_quad = taij_quad ;
	    tmp_taij_quad.dec_dzpuis(4) ;
	    tmp_taij_quad.std_spectral_base() ;

	    des_profile( tmp_taij_quad, 0., 20, M_PI/2., 0.,
		     "A_ij A^ij",
		     "A_ij A^ij (theta=pi/2, phi=0)" ) ;
	    */

	    source_lapse = tmp1 + tmp2 ;

	}

	source_lapse.annule_domain(0) ;
	source_lapse.set_dzpuis(4) ;
	source_lapse.std_spectral_base() ;

	/*
	Scalar tmp_source = source_lapse ;
	tmp_source.dec_dzpuis(4) ;
	tmp_source.std_spectral_base() ;

	des_profile( tmp_source, 0., 20, M_PI/2., 0.,
		     "Source term of lapse",
		     "source_lapse (theta=pi/2, phi=0)" ) ;
	*/

	bc_laps = bc_lapse(neumann, first) ;


	if (kerrschild) {

	    lapse_rs.set_etat_qcq() ;

	    if (neumann) {
	        lapse_rs = source_lapse.poisson_neumann(bc_laps, 0) ;
	    }
	    else {
	        lapse_rs = source_lapse.poisson_dirichlet(bc_laps, 0) ;
	    }

	    // Re-construction of the lapse function
	    // -------------------------------------
	    lapse_rs.annule_domain(0) ;
	    lapse_rs.raccord(1) ;

	    lapse = lapse_rs + lapse_bh ;
	    lapse.annule_domain(0) ;
	    lapse.raccord(1) ;

	}
	else {  // Isotropic coordinates with the maximal slicing

	    lapse_m1.set_etat_qcq() ;

	    if (neumann) {
	        lapse_m1 = source_lapse.poisson_neumann(bc_laps, 0) ;
	    }
	    else {
	        lapse_m1 = source_lapse.poisson_dirichlet(bc_laps, 0) ;
	    }

	    // Re-construction of the lapse function
	    // -------------------------------------
	    lapse = lapse_m1 + 1. ;
	    lapse.annule_domain(0) ;
	    lapse.raccord(1) ;
	    /*
	    des_profile( lapse, 0., 20, M_PI/2., 0.,
			 "Lapse function of BH",
			 "Lapse (theta=pi/2, phi=0)" ) ;
	    */
	}

	//---------------------------------------------------------------//
	//  Resolution of the Poisson equation for the conformal factor  //
	//---------------------------------------------------------------//

	// Source term
	// -----------

	if (kerrschild) {

	    Scalar dcnf = confo_jm1.dsdr() ;
	    Scalar ddcnf = dcnf.dsdr() ;

	    Scalar tmp1 = - 0.125 * taij_quad / pow(confo_jm1, 7.) ;
	    tmp1.std_spectral_base() ;
	    //	    cout << tmp1.get_dzpuis() << endl ;
	    tmp1.inc_dzpuis(4-tmp1.get_dzpuis()) ;

	    Scalar tmp2 = mass*mass*pow(lapse_bh,4.)*confo_jm1/pow(rr, 4) ;
	    tmp2.std_spectral_base() ;
	    //	    cout << tmp2.get_dzpuis() << endl ;
	    tmp2.inc_dzpuis(4-tmp2.get_dzpuis()) ;

	    Scalar confo5(mp) ;
	    confo5 = pow(confo_jm1, 5.) ;
	    confo5.std_spectral_base() ;

	    //	    Scalar tmp3 = confo5 % trace_k % trace_k / 12. ;
	    Scalar tmp3 = confo5 * pow(lapse_bh,6.)*mass*mass
	      *(1.+3.*mass/rr)*(1.+3.*mass/rr)/3./pow(rr,4.) ;
	    tmp3.std_spectral_base() ;
	    //	    cout << tmp3.get_dzpuis() << endl ;
	    tmp3.inc_dzpuis(4-tmp3.get_dzpuis()) ;

	    Scalar tmp4 = 2. * mass * lapse_bh % lapse_bh % ddcnf / rr ;
	    tmp4.std_spectral_base() ;
	    //	    cout << tmp4.get_dzpuis() << endl ;
	    tmp4.inc_dzpuis(4-tmp4.get_dzpuis()) ;

	    Scalar tmp5 = pow(lapse_bh,4.)*dcnf*mass*(3.+8.*mass/rr)/rr/rr ;
	    tmp5.std_spectral_base() ;
	    //	    cout << tmp5.get_dzpuis() << endl ;
	    tmp5.inc_dzpuis(4-tmp5.get_dzpuis()) ;

	    source_confo = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 ;

	    //	    source_confo = 0. ;

	}
	else {  // Isotropic coordinates with the maximal slicing

	    Scalar tmp1 = - 0.125 * taij_quad / pow(confo_jm1, 7.) ;
	    tmp1.std_spectral_base() ;
	    tmp1.inc_dzpuis(4-tmp1.get_dzpuis()) ;

	    source_confo = tmp1 ;

	}

	source_confo.annule_domain(0) ;
	source_confo.set_dzpuis(4) ;
	source_confo.std_spectral_base() ;

	bc_conf = bc_confo() ;

	confo_m1.set_etat_qcq() ;

	confo_m1 = source_confo.poisson_neumann(bc_conf, 0) ;

	// Re-construction of the conformal factor
	// ---------------------------------------

	confo = confo_m1 + 1. ;
	confo.annule_domain(0) ;
	confo.raccord(1) ;

	//-----------------------------------------------------------//
	//  Resolution of the Poisson equation for the shift vector  //
	//-----------------------------------------------------------//

	// Source term
	// -----------

	confo6 = pow(confo_jm1, 6.) ;
	confo6.std_spectral_base() ;

	Vector dlappsi(mp, COV, mp.get_bvect_cart()) ;
	for (int i=1; i<=3; i++)
	    dlappsi.set(i) = (lapse_jm1.deriv(i)
			      - 6.*lapse*confo_jm1.deriv(i)/confo_jm1)
	      / confo6 ;

	dlappsi.std_spectral_base() ;
	dlappsi.annule_domain(0) ;
	/*
	for (int i=1; i<=3; i++)
	    (dlappsi.set(i)).inc_dzpuis(4-dlappsi(i).get_dzpuis()) ;
	*/
	/*
	Vector dlappsi_exact(mp, COV, mp.get_bvect_cart()) ;
	for (int i=1; i<=3; i++)
	    dlappsi_exact.set(i) = mass * ll(i) / pow(1.+2.*mass/rr,1.5)
	      /pow(rr,2.) ;

	dlappsi_exact.annule_domain(0) ;
	dlappsi_exact.std_spectral_base() ;
	for (int i=1; i<=3; i++)
	    (dlappsi_exact.set(i)).dec_dzpuis(dlappsi_exact(i).get_dzpuis()) ;
	*/

	if (kerrschild) {

	    Vector llcomb(mp, COV, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        llcomb.set(i) = pow(lapse_bh,3.)*ll(i)*mass/rr/rr/confo6 ;

	    llcomb.annule_domain(0) ;
	    llcomb.std_spectral_base() ;

	    Vector dtrk(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        dtrk.set(i) = -2.*pow(lapse_bh,7.)*mass
		  *(2.+10.*mass/rr+9.*mass*mass/rr/rr)*ll(i)/pow(rr,3.) ;

	    dtrk.annule_domain(0) ;
	    dtrk.std_spectral_base() ;
	    /*
	    for (int i=1; i<=3; i++)
	        (dtrk.set(i)).inc_dzpuis(4-dtrk(i).get_dzpuis()) ;
	    */
	    /*
	    Vector dtrk_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        dtrk_exact.set(i) = -2.*pow(lap_bh2,3.5)*mass
		  *(2.+10.*mass/rr+9.*mass*mass/rr/rr)*ll(i)/pow(rr,3.) ;

	    dtrk_exact.annule_domain(0) ;
	    dtrk_exact.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (dtrk_exact.set(i)).dec_dzpuis(dtrk_exact(i).get_dzpuis()) ;
	    */

	    Vector lldsh(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        lldsh.set(i) = shift_jm1(i).dsdr() ;

	    lldsh.std_spectral_base() ;
	    /*
	    for (int i=1; i<=3; i++)
	        (lldsh.set(i)).inc_dzpuis(4-lldsh(i).get_dzpuis()) ;
	    */
	    /*
	    Vector lldsh_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        lldsh_exact.set(i) = -2.*lap_bh2*lap_bh2*mass*ll(i)/rr/rr ;
	    lldsh_exact.annule_domain(0) ;
	    lldsh_exact.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (lldsh_exact.set(i)).dec_dzpuis(lldsh_exact(i).get_dzpuis()) ;
	    */

	    Vector lldlldsh(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        lldlldsh.set(i) = lldsh(i).dsdr() ;

	    lldlldsh.std_spectral_base() ;
	    /*
	    for (int i=1; i<=3; i++)
	        (lldlldsh.set(isource_shift = source_shift - sou_shif_exact ;)).inc_dzpuis(4-(lldlldsh(i).get_dzpuis())) ;
	    */

	    /*
	    Vector lldlldsh_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        lldlldsh_exact.set(i) = 4.*pow(lap_bh2,3.)*mass*ll(i)
		  / pow(rr,3.) ;

	    lldlldsh_exact.annule_domain(0) ;
	    lldlldsh_exact.std_spectral_base() ;

	    lldlldsh.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (lldlldsh_exact.set(i)).dec_dzpuis(lldlldsh_exact(i).get_dzpuis()) ;
	    */

	    Scalar divshf(mp) ;
	    divshf = shift_jm1(1).deriv(1) + shift_jm1(2).deriv(2)
	      + shift_jm1(3).deriv(3) ;
	    divshf.std_spectral_base() ;
	    //	    divshf.inc_dzpuis(4-divshf.get_dzpuis()) ;

	    /*
	    Scalar divshf_exact(mp) ;
	    divshf_exact = 2.*lap_bh2*lap_bh2*mass*(1.+4.*mass/rr)/rr/rr ;
	    divshf_exact.std_spectral_base() ;
	    divshf_exact.dec_dzpuis(divshf_exact.get_dzpuis()) ;
	    */

	    Scalar llddivsh = divshf.dsdr() ;
	    llddivsh.std_spectral_base() ;
	    //	    llddivsh.inc_dzpuis(4-llddivsh.get_dzpuis()) ;

	    /*
	    Scalar llddivsh_exact(mp) ;
	    llddivsh_exact = -4.*pow(lap_bh2,3.)*mass
	      *(1.+6.*mass/rr+4.*mass*mass/rr/rr)/pow(rr,3.) ;
	    llddivsh_exact.std_spectral_base() ;
	    llddivsh_exact.dec_dzpuis(llddivsh_exact.get_dzpuis()) ;
	    */

	    Scalar llshift(mp) ;
	    llshift = ll(1)%shift_jm1(1) + ll(2)%shift_jm1(2)
	      + ll(3)%shift_jm1(3) ;
	    llshift.std_spectral_base() ;
	    //	    llshift.inc_dzpuis(4-llshift.get_dzpuis()) ;

	    /*
	    Scalar llshift_exact(mp) ;
	    llshift_exact = 2.*lap_bh2*mass/rr ;
	    llshift_exact.std_spectral_base() ;
	    llshift_exact.dec_dzpuis(llshift_exact.get_dzpuis()) ;
	    */

	    Vector dllsh(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	      dllsh.set(i) = llshift.deriv(i) ;

	    dllsh.std_spectral_base() ;
	    /*
	    for (int i=1; i<=3; i++)
	        (dllsh.set(i)).inc_dzpuis(4-dllsh(i).get_dzpuis()) ;
	    */
	    /*
	    Vector dllsh_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        dllsh_exact.set(i) = -2.*lap_bh2*lap_bh2*mass*ll(i)/rr/rr ;

	    dllsh_exact.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (dllsh_exact.set(i)).dec_dzpuis(dllsh_exact(i).get_dzpuis()) ;
	    */

	    Scalar lldllsh = llshift.dsdr() ;
	    lldllsh.std_spectral_base() ;
	    //	    lldllsh.inc_dzpuis(4-lldllsh.get_dzpuis()) ;

	    /*
	    Scalar lldllsh_exact(mp) ;
	    lldllsh_exact = -2.*lap_bh2*lap_bh2*mass/rr/rr ;
	    lldllsh_exact.std_spectral_base() ;
	    lldllsh_exact.dec_dzpuis(lldllsh_exact.get_dzpuis()) ;
	    */

	    Vector tmp0 = 2. * contract(taij_rs, 1, llcomb, 0) ;
	    tmp0.annule_domain(0) ;
	    tmp0.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp0.set(i)).inc_dzpuis(4-tmp0(i).get_dzpuis()) ;

	    Vector tmp1 = 2. * contract(taij, 1, dlappsi, 0) ;
	    tmp1.annule_domain(0) ;
	    tmp1.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp1.set(i)).inc_dzpuis(4-tmp1(i).get_dzpuis()) ;

	    /*
	    Vector tmp1_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        tmp1_exact.set(i) = -8.*pow(lap_bh2,4.)*mass*mass
		  *(2.+3.*mass/rr)*ll(i) / 3. / pow(rr,4.) ;

	    tmp1_exact.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp1_exact.set(i)).dec_dzpuis(tmp1_exact(i).get_dzpuis()) ;
	    */

	    Vector tmp2 = 4. * lapse_jm1 * dtrk / 3. ;
	    tmp2.annule_domain(0) ;
	    tmp2.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmp2.set(i)).inc_dzpuis(4-(tmp2(i).get_dzpuis())) ;

	    /*
	    Vector tmp2_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        tmp2_exact.set(i) = -8.*pow(lap_bh2,4.)*mass
		  *(2.+10.*mass/rr+9.*mass*mass/rr/rr)*ll(i)/3./pow(rr,3.) ;
	    tmp2_exact.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmp2_exact.set(i)).dec_dzpuis((tmp2_exact(i).get_dzpuis())) ;
	    */

	    Vector tmp3 = pow(lapse_bh,4.)*lldsh*mass*(3.+8.*mass/rr)/rr/rr ;
	    tmp3.annule_domain(0) ;
	    tmp3.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp3.set(i)).inc_dzpuis(4-(tmp3(i).get_dzpuis())) ;

	    /*
	    Vector tmp3_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        tmp3_exact.set(i) = -2.*pow(lap_bh2,4.)*mass*mass
		  *(3.+8.*mass/rr)*ll(i)/pow(rr,4.) ;
	    tmp3_exact.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp3_exact.set(i)).dec_dzpuis((tmp3_exact(i).get_dzpuis())) ;
	    */

	    Vector tmp4 = -pow(lapse_bh,4.)*shift_jm1*mass*(3.+8.*mass/rr)
	      /pow(rr,3.) ;
	    tmp4.annule_domain(0) ;
	    tmp4.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp4.set(i)).inc_dzpuis(4-(tmp4(i).get_dzpuis())) ;

	    /*
	    Vector tmp4_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        tmp4_exact.set(i) = -2.*pow(lap_bh2,3.)*mass*mass
		  *(3.+8.*mass/rr)*ll(i)/pow(rr,4.) ;
	    tmp4_exact.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp4_exact.set(i)).dec_dzpuis((tmp4_exact(i).get_dzpuis())) ;
	    */

	    Vector tmp5 = 2.*mass*lapse_bh*lapse_bh
	      *(lldlldsh+llddivsh*ll/3.)/rr ;
	    tmp5.annule_domain(0) ;
	    tmp5.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp5.set(i)).inc_dzpuis(4-(tmp5(i).get_dzpuis())) ;

	    /*
	    Vector tmp5_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        tmp5_exact.set(i) = 16.*pow(lap_bh2,4.)*mass*mass
		  *(1.-3.*mass/rr-2.*mass*mass/rr/rr)*ll(i)/3./pow(rr,4.) ;
	    tmp5_exact.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp5_exact.set(i)).dec_dzpuis((tmp5_exact(i).get_dzpuis())) ;
	    */

	    Vector tmp6 = mass*lapse_bh*lapse_bh
	      *(dllsh/3.-4.*divshf*ll)/rr/rr ;
	    tmp6.annule_domain(0) ;
	    tmp6.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp6.set(i)).inc_dzpuis(4-(tmp6(i).get_dzpuis())) ;

	    /*
	    Vector tmp6_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        tmp6_exact.set(i) = -2.*pow(lap_bh2,3.)*mass*mass
		  *(13.+48.*mass/rr)*ll(i)/3./pow(rr,4.) ;
	    tmp6_exact.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp6_exact.set(i)).dec_dzpuis((tmp6_exact(i).get_dzpuis())) ;
	    */

	    Vector tmp7 = 2.*pow(lapse_bh,4.)*lldllsh*mass*ll*(9.+11.*mass/rr)
	      /3./rr/rr ;
	    tmp7.annule_domain(0) ;
	    tmp7.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp7.set(i)).inc_dzpuis(4-(tmp7(i).get_dzpuis())) ;

	    /*
	    Vector tmp7_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        tmp7_exact.set(i) = -4.*pow(lap_bh2,4.)*mass*mass
		  *(9.+11.*mass/rr)*ll(i)/3./pow(rr,4.) ;
	    tmp7_exact.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp7_exact.set(i)).dec_dzpuis((tmp7_exact(i).get_dzpuis())) ;
	    */

	    Vector tmp8 = pow(lapse_bh,6.)*llshift*mass*ll
	      *(25.+106.*mass/rr+96.*mass*mass/rr/rr)/3./pow(rr,3.) ;
	    tmp8.annule_domain(0) ;
	    tmp8.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp8.set(i)).inc_dzpuis(4-(tmp8(i).get_dzpuis())) ;

	    /*
	    Vector tmp8_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        tmp8_exact.set(i) = 2.*pow(lap_bh2,4.)*mass*mass
		  *(25.+106.*mass/rr+96.*mass*mass/rr/rr)*ll(i)/3./pow(rr,4.) ;
	    tmp8_exact.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp8_exact.set(i)).dec_dzpuis((tmp8_exact(i).get_dzpuis())) ;
	    */

	    Vector tmp9 = 8.*pow(lapse_bh,8.)*mass*mass*(2.+3.*mass/rr)
	      *(1. - lapse_bh/lapse) * ll / 3. / pow(rr,4.) ;
	    tmp9.annule_domain(0) ;
	    tmp9.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (tmp9.set(i)).inc_dzpuis(4-(tmp9(i).get_dzpuis())) ;

	    source_shift = tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6
	        + tmp7 + tmp8 + tmp9 ;

	    /*
	    Vector sou_shif_exact(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++)
	        sou_shif_exact.set(i) = -16.*pow(lap_bh2,3.)*mass
		  *(1.+6.*mass/rr+4.*mass*mass/rr/rr)*ll(i)/3./pow(rr,3.) ;
	    sou_shif_exact.annule_domain(0) ;
	    sou_shif_exact.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	      (sou_shif_exact.set(i)).inc_dzpuis(4-sou_shif_exact(i).get_dzpuis()) ;

	    source_shift = source_shift - sou_shif_exact ;
	    */
	    /*
	    Vector tmp9 = tmp1_exact + tmp2_exact + tmp3_exact + tmp4_exact
	      + tmp5_exact + tmp6_exact + tmp7_exact + tmp8_exact ;
	    cout << "source_shift - tmp9" << endl ;
	    cout << source_shift(1) - tmp9(1) << endl ;
	    arrete() ;
	    */
	    /*
	    cout << "source_shift - sou_shif_exact" << endl ;
	    cout << source_shift(1) - sou_shif_exact(1) << endl ;
	    arrete() ;
	    */

	}
	else {  // Isotropic coordinates with the maximal slicing

	    source_shift = 2. * contract(taij, 1, dlappsi, 0) ;

	}

	source_shift.annule_domain(0) ;
	/*
	for (int i=1; i<=3; i++)
	    (source_shift.set(i)).raccord(1) ;
	*/
	/*
	for (int i=1; i<=3; i++) {
	    (source_shift.set(i)).set_dzpuis(4) ;
	}
	*/
	source_shift.std_spectral_base() ;

	for (int i=1; i<=3; i++) {
	    if (source_shift(i).get_etat() != ETATZERO)
	        source_shift.set(i).filtre(4) ;
	}

	/*
	for (int i=1; i<=3; i++) {
	    if (source_shift(i).dz_nonzero()) {
	        assert( source_shift(i).get_dzpuis() == 4 ) ;
	    }
	    else {
	        (source_shift.set(i)).set_dzpuis(4) ;
	    }
	}
	*/

	Tenseur source_p(mp, 1, CON, mp.get_bvect_cart()) ;
	source_p.set_etat_qcq() ;
	for (int i=0; i<3; i++) {
	    source_p.set(i) = Cmp(source_shift(i+1)) ;
	}

	Tenseur resu_p(mp, 1, CON, mp.get_bvect_cart()) ;
	resu_p.set_etat_qcq() ;

	for (int i=0; i<3; i++) {
	    resu_p.set(i) = shift_jm1(i+1) ;
	}

	bc_shif_x = bc_shift_x(0.) ;  // Non-rotating BH
	bc_shif_y = bc_shift_y(0.) ;  // Non-rotating BH
	bc_shif_z = bc_shift_z() ;
	/*
	cout << bc_shif_x << endl ;
	arrete() ;
	cout << bc_shif_y << endl ;
	arrete() ;
	cout << bc_shif_z << endl ;
	arrete() ;
	*/
	poisson_vect_frontiere(1./3., source_p, resu_p,
			       bc_shif_x, bc_shif_y, bc_shif_z,
			       0, precis, 20) ;


	if (kerrschild) {
	    for (int i=1; i<=3; i++)
	        shift_rs.set(i) = resu_p(i-1) ;

	    for (int i=1; i<=3; i++)
	        shift.set(i) = shift_rs(i) + shift_bh(i) ;

	    shift_rs.annule_domain(0) ;
	}
	else {  // Isotropic coordinates with the maximal slicing
	    for (int i=1; i<=3; i++)
	        shift.set(i) = resu_p(i-1) ;
	}

	shift.annule_domain(0) ;

	for (int i=1; i<=3; i++)
	    shift.set(i).raccord(1) ;


	/*
	Tbl diff_shftx = diffrel(shift(1), shift_ex(1)) ;
	double diff_shfx = diff_shftx(1) ;
	for (int l=2; l<nz; l++) {
	    diff_shfx += diff_shftx(l) ;
	}
	diff_shfx /= nz ;

	cout << "diff_shfx : " << diff_shfx << endl ;
	*/

	//------------------------------------------------//
	//  Relative difference in the metric quantities  //
	//------------------------------------------------//

	if (kerrschild) {

	    cout << "Mass_bh : " << mass_bh / msol << " [M_sol]" << endl ;
	    double rad_apphor = rad_ah() ;
	    cout << "        : " <<  0.5 * rad_apphor / ggrav / msol
		 << " [M_sol]" << endl ;

	}
	else {  // Isotropic coordinates with the maximal slicing

	    cout << "Mass_bh : " << mass_bh / msol << " [M_sol]" << endl ;

	}
	/*
	des_profile( lapse, 0., 20, M_PI/2., 0.,
		     "Lapse function of BH",
		     "Lapse (theta=pi/2, phi=0)" ) ;

	des_profile( confo, 0., 20, M_PI/2., 0.,
		     "Conformal factor of BH",
		     "Confo (theta=pi/2, phi=0)" ) ;

	des_coupe_vect_z( shift, 0., -3., 0.5, 4,
		       "Shift vector of BH") ;
	*/
	// Difference is calculated only outside the inner boundary.

	// Lapse function
	// --------------
	Tbl diff_lapse(nz) ;

	if (kerrschild) {

	    diff_lapse = diffrel(lapse_rs, lapse_jm1) ;

	}
	else {  // Isotropic coordinates with the maximal slicing

	    diff_lapse = diffrel(lapse, lapse_jm1) ;

	}
	/*
	cout << "lapse: " << endl ;
	for (int l=0; l<nz; l++) {
	  cout << lapse.val_grid_point(l,0,0,0) << endl ;
	}

	cout << "lapse_jm1: " << endl ;
	for (int l=0; l<nz; l++) {
	  cout << lapse_jm1.val_grid_point(l,0,0,0) << endl ;
	}
	*/

	cout << "Relative difference in the lapse function   : " << endl ;
	for (int l=0; l<nz; l++) {
	    cout << diff_lapse(l) << "  " ;
	}
	cout << endl ;

	diff_lp = diff_lapse(1) ;
	for (int l=2; l<nz; l++) {
	    diff_lp += diff_lapse(l) ;
	}
	diff_lp /= nz ;

	// Conformal factor
	// ----------------
	Tbl diff_confo = diffrel(confo, confo_jm1) ;

	cout << "Relative difference in the conformal factor : " << endl ;
	for (int l=0; l<nz; l++) {
	    cout << diff_confo(l) << "  " ;
	}
	cout << endl ;

	diff_cf = diff_confo(1) ;
	for (int l=2; l<nz; l++) {
	    diff_cf += diff_confo(l) ;
	}
	diff_cf /= nz ;

	// Shift vector
	// ------------
	Tbl diff_shift_x(nz) ;
	Tbl diff_shift_y(nz) ;
	Tbl diff_shift_z(nz) ;

	if (kerrschild) {

	    diff_shift_x = diffrel(shift_rs(1), shift_jm1(1)) ;

	}
	else { // Isotropic coordinates with the maximal slicing

	    diff_shift_x = diffrel(shift(1), shift_jm1(1)) ;

	}

	cout << "Relative difference in the shift vector (x) : " << endl ;
	for (int l=0; l<nz; l++) {
	    cout << diff_shift_x(l) << "  " ;
	}
	cout << endl ;

	diff_sx = diff_shift_x(1) ;
	for (int l=2; l<nz; l++) {
	    diff_sx += diff_shift_x(l) ;
	}
	diff_sx /= nz ;

	if (kerrschild) {

	    diff_shift_y = diffrel(shift_rs(2), shift_jm1(2)) ;

	}
	else { // Isotropic coordinates with the maximal slicing

	    diff_shift_y = diffrel(shift(2), shift_jm1(2)) ;

	}

	cout << "Relative difference in the shift vector (y) : " << endl ;
	for (int l=0; l<nz; l++) {
	    cout << diff_shift_y(l) << "  " ;
	}
	cout << endl ;

	diff_sy = diff_shift_y(1) ;
	for (int l=2; l<nz; l++) {
	    diff_sy += diff_shift_y(l) ;
	}
	diff_sy /= nz ;

	if (kerrschild) {

	    diff_shift_z = diffrel(shift_rs(3), shift_jm1(3)) ;

	}
	else { // Isotropic coordinates with the maximal slicing

	    diff_shift_z = diffrel(shift(3), shift_jm1(3)) ;

	}

	cout << "Relative difference in the shift vector (z) : " << endl ;
	for (int l=0; l<nz; l++) {
	    cout << diff_shift_z(l) << "  " ;
	}
	cout << endl ;

	diff_sz = diff_shift_z(1) ;
	for (int l=2; l<nz; l++) {
	    diff_sz += diff_shift_z(l) ;
	}
	diff_sz /= nz ;

	// Next step
	// ---------

	double irr_gm, adm_gm, kom_gm ;
	irr_gm = mass_irr() / mass_bh - 1. ;
	adm_gm = mass_adm() / mass_bh - 1. ;
	kom_gm = mass_kom() / mass_bh - 1. ;
	cout << "Irreducible mass :       " << mass_irr() / msol << endl ;
	cout << "Gravitaitonal mass :     " << mass_bh / msol << endl ;
	cout << "ADM mass :               " << mass_adm() / msol << endl ;
	cout << "Komar mass :             " << mass_kom() / msol << endl ;
	cout << "Diff. (Madm-Mirr)/Mirr : " << mass_adm()/mass_irr() - 1.
	     << endl ;
	cout << "Diff. (Mkom-Mirr)/Mirr : " << mass_kom()/mass_irr() - 1.
	     << endl ;
	cout << "Diff. (Madm-Mkom)/Madm : " << 1. - mass_kom()/mass_adm()
	     << endl ;
	cout << "Diff. (Mirr-Mg)/Mg :     " << irr_gm << endl ;
	cout << "Diff. (Madm-Mg)/Mg :     " << adm_gm << endl ;
	cout << "Diff. (Mkom-Mg)/Mg :     " << kom_gm << endl ;

	cout << endl ;

	del_deriv() ;

	// Relaxation
	/*
	lapse = 0.75 * lapse + 0.25 * lapse_jm1 ;
	confo = 0.75 * confo + 0.25 * confo_jm1 ;
	shift = 0.75 * shift + 0.25 * shift_jm1 ;
	*/
	/*
	des_profile( lapse, 0., 20, M_PI/2., 0.,
		     "Lapse function of BH",
		     "Lapse (theta=pi/2, phi=0)" ) ;

	des_profile( confo, 0., 20, M_PI/2., 0.,
		     "Conformal factor of BH",
		     "Confo (theta=pi/2, phi=0)" ) ;

	des_profile( shift(1), 0., 20, M_PI/2., 0.,
		     "Shift vector (X) of BH",
		     "Shift (theta=pi/2, phi=0)" ) ;
	*/
    }  // End of iteration loop

    //====================================//
    //          End of iteration          //
    //====================================//

    // Exact solution
    // --------------
    Scalar lapse_exact(mp) ;
    Scalar confo_exact(mp) ;
    Vector shift_exact(mp, CON, mp.get_bvect_cart()) ;

    if (kerrschild) {

        lapse_exact = 1./sqrt(1.+2.*mass/rr) ;

        confo_exact = 1. ;

        for (int i=1; i<=3; i++)
	    shift_exact.set(i) = 2.*mass*lapse_exact%lapse_exact%ll(i)/rr ;

    }
    else {

        Scalar rare(mp) ;
	rare = r_coord(neumann, first) ;
	rare.std_spectral_base() ;

        lapse_exact = sqrt(1. - 2.*mass/rare/rr
			   + cc*cc*pow(mass/rare/rr,4.)) ;

	confo_exact = sqrt(rare) ;

	for (int i=1; i<=3; i++) {
	    shift_exact.set(i) = mass * mass * cc * ll(i) / rr / rr
	      / pow(rare,3.) ;
	}

    }

    lapse_exact.annule_domain(0) ;
    lapse_exact.std_spectral_base() ;
    lapse_exact.raccord(1) ;

    confo_exact.annule_domain(0) ;
    confo_exact.std_spectral_base() ;
    confo_exact.raccord(1) ;

    shift_exact.annule_domain(0) ;
    shift_exact.std_spectral_base() ;
    for (int i=1; i<=3; i++)
      shift_exact.set(i).raccord(1) ;

    Scalar lapse_resi = lapse - lapse_exact ;
    Scalar confo_resi = confo - confo_exact ;
    Vector shift_resi = shift - shift_exact ;
    /*
    des_profile( lapse, 0., 20, M_PI/2., 0.,
		 "Lapse function",
		 "Lapse (theta=pi/2, phi=0)" ) ;

    des_profile( lapse_exact, 0., 20, M_PI/2., 0.,
		 "Exact lapse function",
		 "Exact lapse (theta=pi/2, phi=0)" ) ;

    des_profile( lapse_resi, 0., 20, M_PI/2., 0.,
		 "Residual of the lapse function",
		 "Delta Lapse (theta=pi/2, phi=0)" ) ;

    des_profile( confo, 0., 20, M_PI/2., 0.,
		 "Conformal factor",
		 "Confo (theta=pi/2, phi=0)" ) ;

    des_profile( confo_exact, 0., 20, M_PI/2., 0.,
		 "Exact conformal factor",
		 "Exact confo (theta=pi/2, phi=0)" ) ;

    des_profile( confo_resi, 0., 20, M_PI/2., 0.,
		 "Residual of the conformal factor",
		 "Delta Confo (theta=pi/2, phi=0)" ) ;

    des_profile( shift(1), 0., 20, M_PI/2., 0.,
		 "Shift vector (X)",
		 "Shift (X) (theta=pi/2, phi=0)" ) ;

    des_profile( shift_exact(1), 0., 20, M_PI/2., 0.,
		 "Exact shift vector (X)",
		 "Exact shift (X) (theta=pi/2, phi=0)" ) ;

    des_profile( shift_resi(1), 0., 20, M_PI/2., 0.,
		 "Residual of the shift vector X",
		 "Delta shift (X) (theta=pi/2, phi=0)" ) ;
    */
    /*
    des_coupe_vect_z( shift_resi, 0., -3., 0.5, 4,
		      "Delta Shift vector of BH") ;
    */

    // Relative difference in the lapse function
    Tbl diff_lapse_exact = diffrel(lapse, lapse_exact) ;
    diff_lapse_exact.set(0) = 0. ;
    cout << "Relative difference in the lapse function   : " << endl ;
    for (int l=0; l<nz; l++) {
        cout << diff_lapse_exact(l) << "  " ;
    }
    cout << endl ;

    // Relative difference in the conformal factor
    Tbl diff_confo_exact = diffrel(confo, confo_exact) ;
    diff_confo_exact.set(0) = 0. ;
    cout << "Relative difference in the conformal factor : " << endl ;
    for (int l=0; l<nz; l++) {
        cout << diff_confo_exact(l) << "  " ;
    }
    cout << endl ;

    // Relative difference in the shift vector
    Tbl diff_shift_exact_x = diffrel(shift(1), shift_exact(1)) ;
    Tbl diff_shift_exact_y = diffrel(shift(2), shift_exact(2)) ;
    Tbl diff_shift_exact_z = diffrel(shift(3), shift_exact(3)) ;

    cout << "Relative difference in the shift vector (x) : " << endl ;
    for (int l=0; l<nz; l++) {
        cout << diff_shift_exact_x(l) << "  " ;
    }
    cout << endl ;
    cout << "Relative difference in the shift vector (y) : " << endl ;
    for (int l=0; l<nz; l++) {
        cout << diff_shift_exact_y(l) << "  " ;
    }
    cout << endl ;
    cout << "Relative difference in the shift vector (z) : " << endl ;
    for (int l=0; l<nz; l++) {
        cout << diff_shift_exact_z(l) << "  " ;
    }
    cout << endl ;

    //---------------------------------//
    //          Info printing          //
    //---------------------------------//


}
