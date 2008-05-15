/*
 *  Method of class Hole_bhns to compute black-hole metric quantities
 *   in a black hole-neutron star binary
 *
 *    (see file hole_bhns.h for documentation).
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

char hole_bhns_equilibrium_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2008/05/15 19:05:12  k_taniguchi
 * Change of some parameters.
 *
 * Revision 1.1  2007/06/22 01:24:36  k_taniguchi
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
#include "hole_bhns.h"
#include "cmp.h"
#include "tenseur.h"
#include "param.h"
#include "eos.h"
#include "unites.h"
#include "proto.h"
#include "utilitaires.h"
//#include "graphique.h"

void Hole_bhns::equilibrium_bhns(int mer, int mermax_bh,
				 int filter_r, int filter_r_s, int filter_p_s,
				 double x_rot, double y_rot, double precis,
				 double omega_orb, double resize_bh,
				 const Tbl& fact_resize, Tbl& diff) {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    // Initializations
    // ---------------

    const Mg3d* mg = mp.get_mg() ;
    int nz = mg->get_nzone() ;          // total number of domains

    // Re-adjustment of the boundary of domains
    // ----------------------------------------

    double rr_in_1 = mp.val_r(1, -1., M_PI/2, 0.) ;

    /*
    // Three shells outside the shell including NS
    // -------------------------------------------

    // Resize of the outer boundary of the shell including the NS
    double rr_out_nm5 = mp.val_r(nz-5, 1., M_PI/2., 0.) ;
    mp.resize(nz-5, rr_in_1/rr_out_nm5 * fact_resize(1)) ;

    // Resize of the innner boundary of the shell including the NS
    double rr_out_nm6 = mp.val_r(nz-6, 1., M_PI/2., 0.) ;
    mp.resize(nz-6, rr_in_1/rr_out_nm6 * fact_resize(0)) ;

    if (mer % 2 == 0) {

      // Resize of the domain N-2
      double rr_out_nm2 = mp.val_r(nz-2, 1., M_PI/2., 0.) ;
      mp.resize(nz-2, 8. * rr_in_1 * fact_resize(1) / rr_out_nm2) ;

      // Resize of the domain N-3
      double rr_out_nm3 = mp.val_r(nz-3, 1., M_PI/2., 0.) ;
      mp.resize(nz-3, 4. * rr_in_1 * fact_resize(1) / rr_out_nm3) ;

      // Resize of the domain N-4
      double rr_out_nm4 = mp.val_r(nz-4, 1., M_PI/2., 0.) ;
      mp.resize(nz-4, 2. * rr_in_1 * fact_resize(1) / rr_out_nm4) ;

      if (nz > 7) {

	// Resize of the domain 1
	double rr_out_1 = mp.val_r(1, 1., M_PI/2., 0.) ;
	mp.resize(1, rr_in_1/rr_out_1 * resize_bh) ;

	if (nz > 8) {

	  // Resize of the domain from 2 to N-7
	  double rr_out_1_new = mp.val_r(1, 1., M_PI/2., 0.) ;
	  double rr_out_nm6_new = mp.val_r(nz-6, 1., M_PI/2., 0.) ;
	  double dr = (rr_out_nm6_new - rr_out_1_new) / double(nz - 7) ;

	  for (int i=1; i<nz-7; i++) {

	    double rr = rr_out_1_new + i * dr ;
	    double rr_out_ip1 = mp.val_r(i+1, 1., M_PI/2., 0.) ;
	    mp.resize(i+1, rr/rr_out_ip1) ;

	  }

	}

      }

    }
    */

    /*
    // Two shells outside the shell including NS
    // -----------------------------------------

    // Resize of the outer boundary of the shell including the NS
    double rr_out_nm4 = mp.val_r(nz-4, 1., M_PI/2., 0.) ;
    mp.resize(nz-4, rr_in_1/rr_out_nm4 * fact_resize(1)) ;

    // Resize of the innner boundary of the shell including the NS
    double rr_out_nm5 = mp.val_r(nz-5, 1., M_PI/2., 0.) ;
    mp.resize(nz-5, rr_in_1/rr_out_nm5 * fact_resize(0)) ;

    //    if (mer % 2 == 0) {

    // Resize of the domain N-2
    double rr_out_nm2 = mp.val_r(nz-2, 1., M_PI/2., 0.) ;
    mp.resize(nz-2, 3. * rr_in_1 * fact_resize(1) / rr_out_nm2) ;

    // Resize of the domain N-3
    double rr_out_nm3 = mp.val_r(nz-3, 1., M_PI/2., 0.) ;
    mp.resize(nz-3, 1.5 * rr_in_1 * fact_resize(1) / rr_out_nm3) ;

    if (nz > 6) {

      // Resize of the domain 1
      double rr_out_1 = mp.val_r(1, 1., M_PI/2., 0.) ;
      mp.resize(1, rr_in_1/rr_out_1 * resize_bh) ;

      if (nz > 7) {

	// Resize of the domain from 2 to N-6
	double rr_out_nm5_new = mp.val_r(nz-5, 1., M_PI/2., 0.) ;

	for (int i=1; i<nz-6; i++) {

	  double rr_out_i = mp.val_r(i, 1., M_PI/2., 0.) ;

	  double rr_mid = rr_out_i
	    + (rr_out_nm5_new - rr_out_i) / double(nz - 5 - i) ;

	  double rr_2timesi = 2. * rr_out_i ;

	  if (rr_2timesi < rr_mid) {

	    double rr_out_ip1 = mp.val_r(i+1, 1., M_PI/2., 0.) ;
	    mp.resize(i+1, rr_2timesi / rr_out_ip1) ;

	  }
	  else {

	    double rr_out_ip1 = mp.val_r(i+1, 1., M_PI/2., 0.) ;
	    mp.resize(i+1, rr_mid / rr_out_ip1) ;

	  }  // End of else

	}  // End of i loop

      }  // End of (nz > 7) loop

    }  // End of (nz > 6) loop

    //    }  // End of (mer % 2) loop
    */

    // One shell outside the shell including NS
    // ----------------------------------------

    // Resize of the outer boundary of the shell including the NS
    double rr_out_nm3 = mp.val_r(nz-3, 1., M_PI/2., 0.) ;
    mp.resize(nz-3, rr_in_1/rr_out_nm3 * fact_resize(1)) ;

    // Resize of the innner boundary of the shell including the NS
    double rr_out_nm4 = mp.val_r(nz-4, 1., M_PI/2., 0.) ;
    mp.resize(nz-4, rr_in_1/rr_out_nm4 * fact_resize(0)) ;

    //    if (mer % 2 == 0) {

    // Resize of the domain N-2
    double rr_out_nm2 = mp.val_r(nz-2, 1., M_PI/2., 0.) ;
    mp.resize(nz-2, 2. * rr_in_1 * fact_resize(1) / rr_out_nm2) ;

    if (nz > 5) {

      // Resize of the domain 1
      double rr_out_1 = mp.val_r(1, 1., M_PI/2., 0.) ;
      mp.resize(1, rr_in_1/rr_out_1 * resize_bh) ;

      if (nz > 6) {

	// Resize of the domain from 2 to N-5
	double rr_out_nm4_new = mp.val_r(nz-4, 1., M_PI/2., 0.) ;

	for (int i=1; i<nz-5; i++) {

	  double rr_out_i = mp.val_r(i, 1., M_PI/2., 0.) ;

	  double rr_mid = rr_out_i
	    + (rr_out_nm4_new - rr_out_i) / double(nz - 4 - i) ;

	  double rr_2timesi = 2. * rr_out_i ;

	  if (rr_2timesi < rr_mid) {

	    double rr_out_ip1 = mp.val_r(i+1, 1., M_PI/2., 0.) ;
	    mp.resize(i+1, rr_2timesi / rr_out_ip1) ;

	  }
	  else {

	    double rr_out_ip1 = mp.val_r(i+1, 1., M_PI/2., 0.) ;
	    mp.resize(i+1, rr_mid / rr_out_ip1) ;

	  } // End of else

	}  // End of i loop

      }  // End of (nz > 6) loop

    }  // End of (nz > 5) loop

      //    } // End of (mer % 2) loop


    // Inner boundary condition
    // ------------------------

    Valeur bc_lapc(mg->get_angu()) ;
    Valeur bc_conf(mg->get_angu()) ;

    Valeur bc_shif_x(mg->get_angu()) ;
    Valeur bc_shif_y(mg->get_angu()) ;
    Valeur bc_shif_z(mg->get_angu()) ;

    // Error indicators
    // ----------------

    double& diff_lapconf = diff.set(0) ;
    double& diff_confo = diff.set(1) ;
    double& diff_shift_x = diff.set(2) ;
    double& diff_shift_y = diff.set(3) ;
    double& diff_shift_z = diff.set(4) ;

    Scalar lapconf_jm1 = lapconf_auto_rs ; // Lapconf function at previous step
    Scalar confo_jm1 = confo_auto_rs ;  // Conformal factor at preious step
    Vector shift_jm1 = shift_auto_rs ;  // Shift vector at previous step

    // Auxiliary quantities
    // --------------------

    Scalar source_lapconf(mp) ;
    Scalar source_confo(mp) ;
    Vector source_shift(mp, CON, mp.get_bvect_cart()) ;

    Scalar lapconf_m1(mp) ;  // = lapconf_auto_rs + 0.5
    Scalar confo_m1(mp) ;  // = confo_auto_rs + 0.5

    double mass = ggrav * mass_bh ;

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
    ll.set(1) = st % cp ;
    ll.set(2) = st % sp ;
    ll.set(3) = ct ;
    ll.std_spectral_base() ;

    Vector dlappsi(mp, COV, mp.get_bvect_cart()) ;
    for (int i=1; i<=3; i++) {
        dlappsi.set(i) = lapconf_auto_rs.deriv(i) + d_lapconf_comp(i)
	  - 7. * lapconf_tot % (confo_auto_rs.deriv(i) + d_confo_comp(i))
	  / confo_tot ;
    }

    dlappsi.std_spectral_base() ;

    //======================================//
    //          Start of iteration          //
    //======================================//

    for (int mer_bh=0; mer_bh<mermax_bh; mer_bh++) {

        cout << "--------------------------------------------------" << endl ;
	cout << "step: " << mer_bh << endl ;
	cout << "diff_lapconf = " << diff_lapconf << endl ;
	cout << "diff_confo = " << diff_confo << endl ;
	cout << "diff_shift : x = " << diff_shift_x
	     << "  y = " << diff_shift_y << "  z = " << diff_shift_z << endl ;

	if (kerrschild) {

	    cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	    abort() ;

	}  // End of Kerr-Schild
	else { // Isotropic coordinates with the maximal slicing

	    // Sets C/M^2 for each case of the lapse boundary condition
	    // --------------------------------------------------------
	    double cc ;

	    if (bc_lapconf_nd) {  // Neumann boundary condition
	      if (bc_lapconf_fs) {  // First condition
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
	      if (bc_lapconf_fs) {  // First condition
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
		//		cc = 2. * sqrt(2.) ;
	      }
	    }

	    Scalar r_are(mp) ;
	    r_are = r_coord(bc_lapconf_nd, bc_lapconf_fs) ;
	    r_are.std_spectral_base() ;

	    Scalar lapbh_iso(mp) ;
	    lapbh_iso = sqrt(1. - 2.*mass/r_are/rr
			     + cc*cc*pow(mass/r_are/rr,4.)) ;
	    lapbh_iso.std_spectral_base() ;
	    lapbh_iso.annule_domain(0) ;
	    lapbh_iso.raccord(1) ;

	    Scalar psibh_iso(mp) ;
	    psibh_iso = sqrt(r_are) ;
	    psibh_iso.std_spectral_base() ;
	    psibh_iso.annule_domain(0) ;
	    psibh_iso.raccord(1) ;

	    Scalar dlapbh_iso(mp) ;
	    dlapbh_iso = mass/r_are/rr - 2.*cc*cc*pow(mass/r_are/rr,4.) ;
	    dlapbh_iso.std_spectral_base() ;
	    dlapbh_iso.annule_domain(0) ;
	    dlapbh_iso.raccord(1) ;

	    //---------------------------------------------------------------//
	    //  Resolution of the Poisson equation for the lapconf function  //
	    //---------------------------------------------------------------//

	    // Source term
	    // -----------

	    Scalar tmpl1(mp) ;
	    tmpl1 = 0.875 * lapconf_tot % taij_quad_auto
	      / pow(confo_tot, 8.) ;
	    tmpl1.std_spectral_base() ;  // dzpuis = 4
	    tmpl1.annule_domain(0) ;
	    tmpl1.raccord(1) ;

	    Scalar tmpl2(mp) ;
	    tmpl2 = 0.875 * (lapconf_comp+0.5) * taij_quad_comp
	      * (pow(confo_tot/(confo_comp+0.5),6.)*(lapconf_comp+0.5)
		 /lapconf_tot - 1.)
	      / pow(confo_comp+0.5,8.) ;
	    tmpl2.std_spectral_base() ;
	    tmpl2.annule_domain(0) ;  // dzpuis = 4
	    tmpl2.raccord(1) ;

	    Scalar tmpl3(mp) ;
	    tmpl3 = 5.25 * cc * cc * pow(mass,4.) * lapbh_iso
	      * (pow(confo_tot,6.)*lapbh_iso/lapconf_tot - pow(psibh_iso,5.))
	      / pow(r_are*rr,6.) ;
	    tmpl3.std_spectral_base() ;
	    tmpl3.annule_domain(0) ;
	    tmpl3.raccord(1) ;

	    tmpl3.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    source_lapconf = tmpl1 + tmpl2 + tmpl3 ;
	    source_lapconf.std_spectral_base() ;

	    source_lapconf.annule_domain(0) ;
	    source_lapconf.raccord(1) ;
	    /*
	    if (source_lapconf.get_dzpuis() != 4) {
	      source_lapconf.set_dzpuis(4) ;
	    }
	    source_lapconf.std_spectral_base() ;
	    */
	    if (filter_r != 0) {
	      if (source_lapconf.get_etat() != ETATZERO) {
	        source_lapconf.filtre(filter_r) ;
	      }
	    }

	    bc_lapc = bc_lapconf() ;

	    lapconf_m1.set_etat_qcq() ;

	    if (bc_lapconf_nd) {
	      lapconf_m1 = source_lapconf.poisson_neumann(bc_lapc, 0) ;
	    }
	    else {
	      lapconf_m1 = source_lapconf.poisson_dirichlet(bc_lapc, 0) ;
	    }

	    // Re-construction of the lapconf function
	    // ---------------------------------------

	    lapconf_auto_rs = lapconf_m1 - 0.5 ;
	    lapconf_auto_rs.annule_domain(0) ;
	    lapconf_auto_rs.raccord(1) ;

	    lapconf_auto = lapconf_auto_rs + lapconf_auto_bh ;
	    lapconf_auto.annule_domain(0) ; // lapconf_auto,_comp->0.5 (r->inf)
	    lapconf_auto.raccord(1) ;       // lapconf_tot -> 1 (r->inf)


	    //---------------------------------------------------------------//
	    //  Resolution of the Poisson equation for the conformal factor  //
	    //---------------------------------------------------------------//

	    // Source term
	    // -----------

	    Scalar tmpc1 = - 0.125 * taij_quad_auto / pow(confo_tot, 7.) ;
	    tmpc1.std_spectral_base() ; // dzpuis = 4
	    tmpc1.annule_domain(0) ;
	    tmpc1.raccord(1) ;

	    Scalar tmpc2 = 0.75 * cc * cc * pow(mass,4.)
	      * (pow(psibh_iso,5.)
		 - pow(confo_tot,7.)*lapbh_iso*lapbh_iso
		 /lapconf_tot/lapconf_tot)
	      / pow(r_are*rr,6.) ;
	    tmpc2.std_spectral_base() ;
	    tmpc2.annule_domain(0) ;
	    tmpc2.raccord(1) ;

	    tmpc2.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    Scalar tmpc3 = 0.125 * taij_quad_comp
	      * (1. - pow(confo_tot/(confo_comp+0.5),7.)
		 *pow((lapconf_comp+0.5)/lapconf_tot,2.))
	      / pow(confo_comp+0.5, 7.) ;
	    tmpc3.std_spectral_base() ; // dzpuis = 4
	    tmpc3.annule_domain(0) ;
	    tmpc3.raccord(1) ;

	    source_confo = tmpc1 + tmpc2 + tmpc3 ;
	    source_confo.std_spectral_base() ;

	    source_confo.annule_domain(0) ;
	    source_confo.raccord(1) ;
	    /*
	    if (source_confo.get_dzpuis() != 4) {
	      source_confo.set_dzpuis(4) ;
	    }
	    source_confo.std_spectral_base() ;
	    */
	    if (filter_r != 0) {
	      if (source_confo.get_etat() != ETATZERO) {
	        source_confo.filtre(filter_r) ;
	      }
	    }

	    bc_conf = bc_confo(omega_orb, x_rot, y_rot) ;

	    confo_m1.set_etat_qcq() ;

	    confo_m1 = source_confo.poisson_neumann(bc_conf, 0) ;

	    // Re-construction of the conformal factor
	    // ---------------------------------------
	    confo_auto_rs = confo_m1 - 0.5 ;
	    confo_auto_rs.annule_domain(0) ;
	    confo_auto_rs.raccord(1) ;

	    confo_auto = confo_auto_rs + confo_auto_bh ;
	    confo_auto.annule_domain(0) ;  // confo_auto,_comp->0.5 (r->inf)
	    confo_auto.raccord(1) ;        // confo_tot -> 1 (r->inf)


	    //-----------------------------------------------------------//
	    //  Resolution of the Poisson equation for the shift vector  //
	    //-----------------------------------------------------------//

	    // Source term
	    // -----------

	    Vector dlapconf(mp, COV, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	      dlapconf.set(i) = lapconf_auto_rs.deriv(i)
		- 7. * lapconf_tot % confo_auto_rs.deriv(i)
		/ confo_tot ;
	    }

	    dlapconf.std_spectral_base() ;

	    Vector tmps1 = 2. * contract(taij_auto_rs, 1, dlappsi, 0)
	      / pow(confo_tot, 7.)
	      + 2. * contract(taij_comp, 1, dlapconf, 0)
	      * (lapconf_comp+0.5) / lapconf_tot / pow(confo_comp+0.5, 7.) ;
	    tmps1.std_spectral_base() ; // dzpuis = 4
	    tmps1.annule_domain(0) ;
	    for (int i=1; i<=3; i++) {
	        tmps1.set(i).raccord(1) ;
	    }

	    Vector tmps2(mp, CON, mp.get_bvect_cart()) ;
	    tmps2.set_etat_qcq() ;
	    for (int i=1; i<=3; i++) {
	      tmps2.set(i) = 2. * psibh_iso
		* (dlapbh_iso + 0.5*(lapbh_iso - 1.)
		   *(lapbh_iso - 7.*lapconf_tot/confo_tot))
		* (taij_tot_rs(i,1)%ll(1) + taij_tot_rs(i,2)%ll(2)
		   + taij_tot_rs(i,3)%ll(3)) / pow(confo_tot,7.) / rr ;
	    }
	    tmps2.std_spectral_base() ;
	    tmps2.annule_domain(0) ;
	    for (int i=1; i<=3; i++) {
	        tmps2.set(i).raccord(1) ;
	    }
	    for (int i=1; i<=3; i++) {
	        (tmps2.set(i)).inc_dzpuis(2) ;  // dzpuis : 2 -> 4
	    }

	    Vector tmps3(mp, CON, mp.get_bvect_cart()) ;
	    tmps3.set_etat_qcq() ;
	    for (int i=1; i<=3; i++) {
	        tmps3.set(i) = 2. * cc * mass * mass * lapbh_iso
		  * (dlappsi(i) - 3.*ll(i)*(ll(1)%dlappsi(1)
					    + ll(2)%dlappsi(2)
					    + ll(3)%dlappsi(3)))
		  / lapconf_tot / pow(r_are*rr,3.) ;
	    }
	    tmps3.std_spectral_base() ;
	    tmps3.annule_domain(0) ;
	    for (int i=1; i<=3; i++) {
	        tmps3.set(i).raccord(1) ;
	    }
	    for (int i=1; i<=3; i++) {
	        (tmps3.set(i)).inc_dzpuis(2) ;  // dzpuis : 2 -> 4
	    }

	    Vector tmps4 = 4. * cc * mass * mass
	      * (dlapbh_iso * (1. - psibh_iso*lapbh_iso/lapconf_tot)
		 + 0.5 * lapbh_iso * (lapbh_iso - 1.)
		 * (6.*(psibh_iso/confo_tot - 1.)
		    + psibh_iso*(1./confo_tot - lapbh_iso/lapconf_tot)))
	      * ll / rr / pow(r_are*rr,3.) ;
	    tmps4.std_spectral_base() ;
	    tmps4.annule_domain(0) ;
	    for (int i=1; i<=3; i++) {
	        tmps4.set(i).raccord(1) ;
	    }
	    for (int i=1; i<=3; i++) {
	        (tmps4.set(i)).inc_dzpuis(4) ;  // dzpuis : 0 -> 4
	    }

	    Vector dlappsi_comp(mp, COV, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	      dlappsi_comp.set(i) = ((lapconf_comp+0.5)/lapconf_tot - 1.)
		* d_lapconf_comp(i)
		- 7. * (lapconf_comp+0.5) * ((confo_comp+0.5)/confo_tot - 1.)
		* d_confo_comp(i) / (confo_comp+0.5) ;
	    }

	    dlappsi_comp.std_spectral_base() ;

	    Vector tmps5 = 2. * contract(taij_comp, 1, dlappsi_comp, 0)
	      / pow(confo_comp+0.5, 7.) ;
	    tmps5.std_spectral_base() ;
	    tmps5.annule_domain(0) ;
	    for (int i=1; i<=3; i++) {
	        tmps5.set(i).raccord(1) ;
	    }

	    source_shift = tmps1 + tmps2 + tmps3 + tmps4 + tmps5 ;
	    source_shift.std_spectral_base() ;
	    source_shift.annule_domain(0) ;

	    for (int i=1; i<=3; i++) {
	        source_shift.set(i).raccord(1) ;
	    }

	    if (filter_r_s != 0) {
	      for (int i=1; i<=3; i++) {
	        if (source_shift(i).get_etat() != ETATZERO)
		  source_shift.set(i).filtre(filter_r_s) ;
	      }
	    }

	    if (filter_p_s != 0) {
	      for (int i=1; i<=3; i++) {
	        if (source_shift(i).get_etat() != ETATZERO) {
		  source_shift.set(i).filtre_phi(filter_p_s, nz-1) ;
		  /*
		  for (int l=1; l<nz; l++) {
		    source_shift.set(i).filtre_phi(filter_p_s, l) ;
		  }
		  */
		}
	      }
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
	      resu_p.set(i) = Cmp(shift_auto_rs(i+1)) ;
	    }

	    // Boundary condition
	    bc_shif_x = bc_shift_x(omega_orb, y_rot) ;
	    bc_shif_y = bc_shift_y(omega_orb, x_rot) ;
	    bc_shif_z = bc_shift_z() ;

	    poisson_vect_frontiere(1./3., source_p, resu_p,
				   bc_shif_x, bc_shif_y, bc_shif_z,
				   0, precis, 7) ;

	    for (int i=1; i<=3; i++) {
	      shift_auto_rs.set(i) = resu_p(i-1) ;
	    }

	    shift_auto_rs.std_spectral_base() ;
	    shift_auto_rs.annule_domain(0) ;
	    for (int i=1; i<=3; i++) {
	        shift_auto_rs.set(i).raccord(1) ;
	    }

	    shift_auto = shift_auto_rs + shift_auto_bh ;
	    shift_auto.std_spectral_base() ;
	    shift_auto.annule_domain(0) ;
	    for (int i=1; i<=3; i++) {
	        shift_auto.set(i).raccord(1) ;
	    }

	}  // End of isotropic

	//------------------------------------------------//
	//  Relative difference in the metric quantities  //
	//------------------------------------------------//

	// Difference is calculated only outside the inner boundary.

	Tbl tdiff_lapconf = diffrel(lapconf_auto_rs, lapconf_jm1) ;
	tdiff_lapconf.set(0) = 0. ;
	cout << "Relative difference in the lapconf function : " << endl ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_lapconf(l) << "  " ;
	}
	cout << endl ;

	diff_lapconf = tdiff_lapconf(1) ;
	for (int l=2; l<nz; l++) {
	    diff_lapconf += tdiff_lapconf(l) ;
	}
	diff_lapconf /= nz ;

	Tbl tdiff_confo = diffrel(confo_auto_rs, confo_jm1) ;
	tdiff_confo.set(0) = 0. ;
	cout << "Relative difference in the conformal factor : " << endl ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_confo(l) << "  " ;
	}
	cout << endl ;

	diff_confo = tdiff_confo(1) ;
	for (int l=2; l<nz; l++) {
	    diff_confo += tdiff_confo(l) ;
	}
	diff_confo /= nz ;

	Tbl tdiff_shift_x = diffrel(shift_auto_rs(1), shift_jm1(1)) ;
	tdiff_shift_x.set(0) = 0. ;
	cout << "Relative difference in the shift vector (x) : " << endl ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_shift_x(l) << "  " ;
	}
	cout << endl ;

	diff_shift_x = tdiff_shift_x(1) ;
	for (int l=2; l<nz; l++) {
	    diff_shift_x += tdiff_shift_x(l) ;
	}
	diff_shift_x /= nz ;

	Tbl tdiff_shift_y = diffrel(shift_auto_rs(2), shift_jm1(2)) ;
	tdiff_shift_y.set(0) = 0. ;
	cout << "Relative difference in the shift vector (y) : " << endl ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_shift_y(l) << "  " ;
	}
	cout << endl ;

	diff_shift_y = tdiff_shift_y(1) ;
	for (int l=2; l<nz; l++) {
	    diff_shift_y += tdiff_shift_y(l) ;
	}
	diff_shift_y /= nz ;

	Tbl tdiff_shift_z = diffrel(shift_auto_rs(3), shift_jm1(3)) ;
	tdiff_shift_z.set(0) = 0. ;
	cout << "Relative difference in the shift vector (z) : " << endl ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_shift_z(l) << "  " ;
	}
	cout << endl ;

	diff_shift_z = tdiff_shift_z(1) ;
	for (int l=2; l<nz; l++) {
	    diff_shift_z += tdiff_shift_z(l) ;
	}
	diff_shift_z /= nz ;

	/*
	des_profile( lapconf_auto_rs, 0., 10.,
		     M_PI/2., 0., "Residual lapconf function of BH",
		     "Lapconf (theta=pi/2, phi=0)" ) ;

	des_profile( lapconf_auto_bh, 0., 10.,
		     M_PI/2., 0., "Analytic lapconf function of BH",
		     "Lapconf (theta=pi/2, phi=0)" ) ;

	des_profile( lapconf_auto, 0., 10.,
		     M_PI/2., 0., "Self lapconf function of BH",
		     "Lapconf (theta=pi/2, phi=0)" ) ;

	des_profile( lapconf_tot, 0., 10.,
		     M_PI/2., 0., "Total lapconf function of BH",
		     "Lapconf (theta=pi/2, phi=0)" ) ;

	des_profile( confo_auto_rs, 0., 10.,
		     M_PI/2., 0., "Residual conformal factor of BH",
		     "Confo (theta=pi/2, phi=0)" ) ;

	des_profile( confo_auto_bh, 0., 10.,
		     M_PI/2., 0., "Analytic conformal factor of BH",
		     "Confo (theta=pi/2, phi=0)" ) ;

	des_profile( confo_auto, 0., 10.,
		     M_PI/2., 0., "Self conformal factor of BH",
		     "Confo (theta=pi/2, phi=0)" ) ;

	des_profile( confo_tot, 0., 10.,
		     M_PI/2., 0., "Total conformal factor of BH",
		     "Confo (theta=pi/2, phi=0)" ) ;

	des_coupe_vect_z( shift_auto_rs, 0., -3., 0.5, 3,
			  "Residual shift vector of NS") ;

	des_coupe_vect_z( shift_auto_bh, 0., -3., 0.5, 3,
			  "Analytic shift vector of NS") ;

	des_coupe_vect_z( shift_auto, 0., -3., 0.5, 3,
			  "Self shift vector of NS") ;

	des_coupe_vect_z( shift_tot, 0., -3., 0.5, 3,
			  "Total Shift vector seen by NS") ;
	*/
    } // End of main loop

    //====================================//
    //          End of iteration          //
    //====================================//

}
