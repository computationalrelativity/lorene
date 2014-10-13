/*
 *  Method of class Black_hole to compute a single black hole
 *
 *    (see file blackhole.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005-2007 Keisuke Taniguchi
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
 * Revision 1.5  2014/10/13 08:52:46  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:02  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2008/07/02 20:42:53  k_taniguchi
 * Modification of the argument and so on.
 *
 * Revision 1.2  2008/05/15 19:26:30  k_taniguchi
 * Change of some parameters.
 *
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
#include <cmath>

// Lorene headers
#include "blackhole.h"
#include "cmp.h"
#include "tenseur.h"
#include "param.h"
#include "unites.h"
#include "proto.h"
#include "utilitaires.h"
#include "graphique.h"

namespace Lorene {
void Black_hole::equilibrium_spher(bool neumann, bool first,
				   double spin_omega, double precis,
				   double precis_shift) {

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

    Valeur bc_lpcnf(mg->get_angu()) ;
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
	  // d(\alpha \psi)/dr = 0
	  // ---------------------
	  cc = 2. * (sqrt(13.) - 1.) / 3. ;
	}
	else {  // Second condition
	  // d(\alpha \psi)/dr = (\alpha \psi)/(2 rah)
	  // -----------------------------------------
	  cc = 4. / 3. ;
	}
    }
    else {  // Dirichlet boundary condition
       if (first) {  // First condition
	 // (\alpha \psi) = 1/2
	 // -------------------
	 cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
	 abort() ;
       }
       else {  // Second condition
	 // (\alpha \psi) = 1/sqrt(2.) \psi_KS
	 // ----------------------------------
	 cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
	 abort() ;
	 //	 cc = 2. * sqrt(2.) ;
       }
    }


    // Ininital metric
    if (kerrschild) {

        lapconf_bh = 1. / sqrt(1. + 2. * mass / rr) ;
	lapconf_bh.std_spectral_base() ;
	lapconf_bh.annule_domain(0) ;
	lapconf_bh.raccord(1) ;

	lapconf = lapconf_bh ;
	lapconf.std_spectral_base() ;
	lapconf_rs = 0. ;
	lapconf_rs.std_spectral_base() ;

        lapse = lapconf ;
	lapse.std_spectral_base() ;

	confo = 1. ;
	confo.std_spectral_base() ;

	Scalar lapse_bh(mp) ;
	lapse_bh = 1. / sqrt(1. + 2. * mass / rr) ;
	lapse_bh.std_spectral_base() ;
	lapse_bh.annule_domain(0) ;
	lapse_bh.raccord(1) ;

	for (int i=1; i<=3; i++) {
	    shift_bh.set(i) = 2. * lapse_bh * lapse_bh * mass * ll(i) / rr ;
	}
	shift_bh.std_spectral_base() ;

	shift = shift_bh ;
	shift.std_spectral_base() ;
	shift_rs.set_etat_zero() ;

    }
    else {  // Isotropic coordinates

        Scalar r_are(mp) ;
	r_are = r_coord(neumann, first) ;
	r_are.std_spectral_base() ;
	r_are.annule_domain(0) ;
	r_are.raccord(1) ;
	/*
	cout << "r_are:" << endl ;
	for (int l=0; l<nz; l++) {
	  cout << r_are.val_grid_point(l,0,0,0) << endl ;
	}
	*/

	// Exact, non-spinning case
	/*
	lapconf = sqrt(1. - 2.*mass/r_are/rr
		       + cc*cc*pow(mass/r_are/rr,4.))
	  * sqrt(r_are) ;
	*/
	lapconf = sqrt(1. - 1.9*mass/r_are/rr
		       + cc*cc*pow(mass/r_are/rr,4.))
	  * sqrt(r_are) ;
	lapconf.std_spectral_base() ;
	lapconf.annule_domain(0) ;
	lapconf.raccord(1) ;

	/*
	lapse = sqrt(1. - 2.*mass/r_are/rr
		     + cc*cc*pow(mass/r_are/rr,4.)) ;
	*/
	lapse = sqrt(1. - 1.9*mass/r_are/rr
		     + cc*cc*pow(mass/r_are/rr,4.)) ;
	lapse.std_spectral_base() ;
	lapse.annule_domain(0) ;
	lapse.raccord(1) ;

	//	confo = sqrt(r_are) ;
	confo = sqrt(0.9*r_are) ;
	confo.std_spectral_base() ;
	confo.annule_domain(0) ;
	confo.raccord(1) ;

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

    Scalar source_lapconf(mp) ;
    Scalar source_confo(mp) ;
    Vector source_shift(mp, CON, mp.get_bvect_cart()) ;

    Scalar lapconf_m1(mp) ;  // = lapconf - 1 (only for the isotropic case)
    Scalar confo_m1(mp) ;  // = confo - 1

    Scalar lapconf_jm1(mp) ;
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
    for (int mer=0; (diff_lp > precis) && (diff_cf > precis) && (diff_sx > precis) && (diff_sy > precis) && (diff_sz > precis) && (mer < mermax); mer++) {

        cout << "--------------------------------------------------" << endl ;
	cout << "step: " << mer << endl ;
	cout << "diff_lapconf = " << diff_lp << endl ;
	cout << "diff_confo =   " << diff_cf << endl ;
	cout << "diff_shift : x = " << diff_sx
	     << "  y = " << diff_sy << "  z = " << diff_sz << endl ;

	if (kerrschild) {
	    lapconf_jm1 = lapconf_rs ;
	    confo_jm1 = confo ;
	    shift_jm1 = shift_rs ;
	}
	else {
	    lapconf_jm1 = lapconf ;
	    confo_jm1 = confo ;
	    shift_jm1 = shift ;
	}

	//------------------------------------------//
	//  Computation of the extrinsic curvature  //
	//------------------------------------------//

	extr_curv_bh() ;

	//---------------------------------------------------------------//
	//  Resolution of the Poisson equation for the lapconf function  //
	//---------------------------------------------------------------//

	// Source term
	// -----------

	if (kerrschild) {

	  cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	  abort() ;

	}
	else {  // Isotropic coordinates with the maximal slicing

	    source_lapconf = 7. * lapconf_jm1 % taij_quad
	      / pow(confo_jm1, 8.) / 8. ;

	}

	source_lapconf.annule_domain(0) ;
	source_lapconf.set_dzpuis(4) ;
	source_lapconf.std_spectral_base() ;

	/*
	Scalar tmp_source = source_lapse ;
	tmp_source.dec_dzpuis(4) ;
	tmp_source.std_spectral_base() ;

	des_profile( tmp_source, 0., 20, M_PI/2., 0.,
		     "Source term of lapse",
		     "source_lapse (theta=pi/2, phi=0)" ) ;
	*/

	bc_lpcnf = bc_lapconf(neumann, first) ;


	if (kerrschild) {

	    cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	    abort() ;
	    /*
	    lapconf_rs.set_etat_qcq() ;

	    if (neumann) {
	        lapconf_rs = source_lapconf.poisson_neumann(bc_lpcnf, 0) ;
	    }
	    else {
	        lapconf_rs = source_lapconf.poisson_dirichlet(bc_lpcnf, 0) ;
	    }

	    // Re-construction of the lapse function
	    // -------------------------------------
	    lapconf_rs.annule_domain(0) ;
	    lapconf_rs.raccord(1) ;

	    lapconf = lapconf_rs + lapconf_bh ;
	    lapconf.annule_domain(0) ;
	    lapconf.raccord(1) ;
	    */
	}
	else {  // Isotropic coordinates with the maximal slicing

	    lapconf_m1.set_etat_qcq() ;

	    if (neumann) {
	        lapconf_m1 = source_lapconf.poisson_neumann(bc_lpcnf, 0) ;
	    }
	    else {
	        lapconf_m1 = source_lapconf.poisson_dirichlet(bc_lpcnf, 0) ;
	    }

	    // Re-construction of the lapse function
	    // -------------------------------------
	    lapconf = lapconf_m1 + 1. ;
	    lapconf.annule_domain(0) ;
	    lapconf.raccord(1) ;
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

	    cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	    abort() ;

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

	Scalar confo7(mp) ;
	confo7 = pow(confo_jm1, 7.) ;
	confo7.std_spectral_base() ;

	Vector dlappsi(mp, COV, mp.get_bvect_cart()) ;
	for (int i=1; i<=3; i++)
	    dlappsi.set(i) = (lapconf_jm1.deriv(i)
			      - 7.*lapconf*confo_jm1.deriv(i)/confo_jm1)
	      / confo7 ;

	dlappsi.std_spectral_base() ;
	dlappsi.annule_domain(0) ;


	if (kerrschild) {

	    cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	    abort() ;

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
	    resu_p.set(i) = Cmp(shift_jm1(i+1)) ;
	}

	bc_shif_x = bc_shift_x(spin_omega) ;  // Rotating BH
	bc_shif_y = bc_shift_y(spin_omega) ;  // Rotating BH
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
			       0, precis_shift, 14) ;


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

	// Lapconf function
	// ----------------
	Tbl diff_lapconf(nz) ;

	if (kerrschild) {

	    diff_lapconf = diffrel(lapconf_rs, lapconf_jm1) ;

	}
	else {  // Isotropic coordinates with the maximal slicing

	    diff_lapconf = diffrel(lapconf, lapconf_jm1) ;

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

	cout << "Relative difference in the lapconf function   : " << endl ;
	for (int l=0; l<nz; l++) {
	    cout << diff_lapconf(l) << "  " ;
	}
	cout << endl ;

	diff_lp = diff_lapconf(1) ;
	for (int l=2; l<nz; l++) {
	    diff_lp += diff_lapconf(l) ;
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

	// Mass
	if (kerrschild) {

	    cout << "Mass_bh : " << mass_bh / msol << " [M_sol]" << endl ;
	    double rad_apphor = rad_ah() ;
	    cout << "        : " <<  0.5 * rad_apphor / ggrav / msol
		 << " [M_sol]" << endl ;

	}
	else {  // Isotropic coordinates with the maximal slicing

	    cout << "Mass_bh :                " << mass_bh / msol
		 << " [M_sol]" << endl ;

	}

	// ADM mass, Komar mass
	// --------------------

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
    Scalar lapconf_exact(mp) ;
    Scalar confo_exact(mp) ;
    Vector shift_exact(mp, CON, mp.get_bvect_cart()) ;

    if (kerrschild) {

        lapconf_exact = 1./sqrt(1.+2.*mass/rr) ;

        confo_exact = 1. ;

        for (int i=1; i<=3; i++)
	    shift_exact.set(i) =
	      2.*mass*lapconf_exact%lapconf_exact%ll(i)/rr ;

    }
    else {

        Scalar rare(mp) ;
	rare = r_coord(neumann, first) ;
	rare.std_spectral_base() ;

        lapconf_exact = sqrt(1. - 2.*mass/rare/rr
			   + cc*cc*pow(mass/rare/rr,4.)) * sqrt(rare) ;

	confo_exact = sqrt(rare) ;

	for (int i=1; i<=3; i++) {
	    shift_exact.set(i) = mass * mass * cc * ll(i) / rr / rr
	      / pow(rare,3.) ;
	}

    }

    lapconf_exact.annule_domain(0) ;
    lapconf_exact.std_spectral_base() ;
    lapconf_exact.raccord(1) ;

    confo_exact.annule_domain(0) ;
    confo_exact.std_spectral_base() ;
    confo_exact.raccord(1) ;

    shift_exact.annule_domain(0) ;
    shift_exact.std_spectral_base() ;
    for (int i=1; i<=3; i++)
      shift_exact.set(i).raccord(1) ;

    Scalar lapconf_resi = lapconf - lapconf_exact ;
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

    // Relative difference in the lapconf function
    Tbl diff_lapconf_exact = diffrel(lapconf, lapconf_exact) ;
    diff_lapconf_exact.set(0) = 0. ;
    cout << "Relative difference in the lapconf function   : " << endl ;
    for (int l=0; l<nz; l++) {
        cout << diff_lapconf_exact(l) << "  " ;
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
}
