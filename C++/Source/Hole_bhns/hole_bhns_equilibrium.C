/*
 *  Method of class Hole_bhns to compute black-hole metric quantities
 *   in a black hole-neutron star binary
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

char hole_bhns_equilibrium_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
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
				 double omega_orb, double omega_spin,
				 double resize_bh,
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

    /*
    // One shell outside the shell including NS
    // ----------------------------------------

    // Resize of the outer boundary of the shell including the NS
    double rr_out_nm3 = mp.val_r(nz-3, 1., M_PI/2., 0.) ;
    mp.resize(nz-3, rr_in_1/rr_out_nm3 * fact_resize(1)) ;

    // Resize of the innner boundary of the shell including the NS
    double rr_out_nm4 = mp.val_r(nz-4, 1., M_PI/2., 0.) ;
    mp.resize(nz-4, rr_in_1/rr_out_nm4 * fact_resize(0)) ;

    if (mer % 2 == 0) {

      // Resize of the domain N-2
      double rr_out_nm2 = mp.val_r(nz-2, 1., M_PI/2., 0.) ;
      mp.resize(nz-2, 2. * rr_in_1 * fact_resize(1) / rr_out_nm2) ;

      if (nz > 5) {

	// Resize of the domain 1
	double rr_out_1 = mp.val_r(1, 1., M_PI/2., 0.) ;
	mp.resize(1, rr_in_1/rr_out_1 * resize_bh) ;

	if (nz > 6) {

	  // Resize of the domain from 2 to N-5
	  double rr_out_1_new = mp.val_r(1, 1., M_PI/2., 0.) ;
	  double rr_out_nm4_new = mp.val_r(nz-4, 1., M_PI/2., 0.) ;
	  double dr = (rr_out_nm4_new - rr_out_1_new) / double(nz - 5) ;

	  for (int i=1; i<nz-5; i++) {

	    double rr = rr_out_1_new + i * dr ;
	    double rr_out_ip1 = mp.val_r(i+1, 1., M_PI/2., 0.) ;
	    mp.resize(i+1, rr/rr_out_ip1) ;

	  }

	}

      }

    }
    */

    // Inner boundary condition
    // ------------------------

    Valeur bc_laps(mg->get_angu()) ;
    Valeur bc_conf(mg->get_angu()) ;

    Valeur bc_shif_x(mg->get_angu()) ;
    Valeur bc_shif_y(mg->get_angu()) ;
    Valeur bc_shif_z(mg->get_angu()) ;

    // Error indicators
    // ----------------

    double& diff_lapse = diff.set(0) ;
    double& diff_confo = diff.set(1) ;
    double& diff_shift_x = diff.set(2) ;
    double& diff_shift_y = diff.set(3) ;
    double& diff_shift_z = diff.set(4) ;

    Scalar lapse_jm1 = lapse_auto_rs ;  // Lapse function at previous step
    Scalar confo_jm1 = confo_auto_rs ;  // Conformal factor at preious step
    Vector shift_jm1 = shift_auto_rs ;  // Shift vector at previous step

    // Auxiliary quantities
    // --------------------

    Scalar source_lapse(mp) ;
    Scalar source_confo(mp) ;
    Vector source_shift(mp, CON, mp.get_bvect_cart()) ;

    Scalar lapse_m1(mp) ;  // = lapse_auto_rs + 0.5 for Kerr-Schild
                           // = lapse_auto_rs + 0.5 for Isotropic
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
    ll.set(1) = st * cp ;
    ll.set(2) = st * sp ;
    ll.set(3) = ct ;
    ll.std_spectral_base() ;

    Scalar dlap(mp) ;
    dlap = lapse_auto_rs.dsdr() + ll(1) * d_lapse_comp(1)
      + ll(2) * d_lapse_comp(2) + ll(3) * d_lapse_comp(3) ;

    dlap.std_spectral_base() ;

    Scalar dlconfo(mp) ;
    dlconfo = (confo_auto_rs.dsdr() + ll(1) * d_confo_comp(1)
	       + ll(2) * d_confo_comp(2)
	       + ll(3) * d_confo_comp(3)) / confo_tot ;

    dlconfo.std_spectral_base() ;

    Scalar dconf = confo_auto_rs.dsdr() + ll(1) * d_confo_comp(1)
      + ll(2) * d_confo_comp(2) + ll(3) * d_confo_comp(3) ;
    dconf.std_spectral_base() ;

    Scalar llshift(mp) ;
    llshift = ll(1) * (shift_auto_rs(1) + shift_comp(1))
      + ll(2) * (shift_auto_rs(2) + shift_comp(2))
      + ll(3) * (shift_auto_rs(3) + shift_comp(3)) ;
    llshift.std_spectral_base() ;

    Vector dlappsi(mp, COV, mp.get_bvect_cart()) ;
    for (int i=1; i<=3; i++) {
        dlappsi.set(i) = lapse_auto_rs.deriv(i) + d_lapse_comp(i)
	  - 6. * lapse_tot * (confo_auto_rs.deriv(i) + d_confo_comp(i))
	  / confo_tot ;
    }

    dlappsi.std_spectral_base() ;

    Vector dshift(mp, CON, mp.get_bvect_cart()) ;
    for (int i=1; i<=3; i++) {
        dshift.set(i) = shift_auto_rs(i).dsdr()
	  + ll(1)*d_shift_comp(1,i) + ll(2)*d_shift_comp(2,i)
	  + ll(3)*d_shift_comp(3,i) ;
    }
    dshift.std_spectral_base() ;

    Scalar divshift = shift_auto_rs(1).deriv(1) + shift_auto_rs(2).deriv(2)
      + shift_auto_rs(3).deriv(3)
      + d_shift_comp(1,1) + d_shift_comp(2,2) + d_shift_comp(3,3) ;
    divshift.std_spectral_base() ;

    Scalar ddivshif = divshift.dsdr() ;
    ddivshif.std_spectral_base() ;

    Scalar llshift_auto(mp) ;
    llshift_auto = ll(1)%shift_auto_rs(1) + ll(2)%shift_auto_rs(2)
      + ll(3)%shift_auto_rs(3) ;
    llshift_auto.std_spectral_base() ;

    Scalar dllshift = llshift_auto.dsdr()
      + ll(1) * ( ll(1)*d_shift_comp(1,1) + ll(2)*d_shift_comp(1,2)
		  + ll(3)*d_shift_comp(1,3) )
      + ll(2) * ( ll(1)*d_shift_comp(2,1) + ll(2)*d_shift_comp(2,2)
		  + ll(3)*d_shift_comp(2,3) )
      + ll(3) * ( ll(1)*d_shift_comp(3,1) + ll(2)*d_shift_comp(3,2)
		  + ll(3)*d_shift_comp(3,3) ) ;
    dllshift.std_spectral_base() ;

    Scalar orb_rot_x(mp) ;
    orb_rot_x = omega_orb * (mp.get_ori_x() - x_rot) ;
    orb_rot_x.std_spectral_base() ;

    Scalar orb_rot_y(mp) ;
    orb_rot_y = omega_orb * (mp.get_ori_y() - y_rot) ;
    orb_rot_y.std_spectral_base() ;


    //======================================//
    //          Start of iteration          //
    //======================================//

    for (int mer_bh=0; mer_bh<mermax_bh; mer_bh++) {

        cout << "--------------------------------------------------" << endl ;
	cout << "step: " << mer_bh << endl ;
	cout << "diff_lapse = " << diff_lapse << endl ;
	cout << "diff_confo = " << diff_confo << endl ;
	cout << "diff_shift : x = " << diff_shift_x
	     << "  y = " << diff_shift_y << "  z = " << diff_shift_z << endl ;

	if (kerrschild) {

	    //-------------------------------------------------------------//
	    //  Resolution of the Poisson equation for the lapse function  //
	    //-------------------------------------------------------------//

	    // Source term
	    // -----------

	    Scalar tmpl1 = lapse_tot * (taij_quad_tot_rs + taij_quad_tot_rot)
	      / pow(confo_tot, 8.) ;
	    tmpl1.std_spectral_base() ;
	    tmpl1.annule_domain(0) ;
	    // dzpuis = 4

	    Scalar tmpl2 = -2.*( (lapse_auto_rs.deriv(1)+d_lapse_comp(1))
				 *(confo_auto.deriv(1)+d_confo_comp(1))
				 + (lapse_auto_rs.deriv(2)+d_lapse_comp(2))
				 *(confo_auto.deriv(2)+d_confo_comp(2))
				 + (lapse_auto_rs.deriv(3)+d_lapse_comp(3))
				 *(confo_auto.deriv(3)+d_confo_comp(3)) )
	      / confo_tot
	      + 4. * lapse_auto_bh * lapse_auto_bh * mass * dlap * dlconfo
	      / rr ;
	    tmpl2.std_spectral_base() ;
	    tmpl2.annule_domain(0) ;
	    // dzpuis = 4

	    Scalar tmpl3 = -2. * mass * pow(lapse_auto_bh,5.) * dlconfo
	      / rr / rr ;
	    tmpl3.std_spectral_base() ;
	    tmpl3.annule_domain(0) ;
	    tmpl3.inc_dzpuis(2) ;  // dzpuis : 2 -> 4

	    Scalar tmpl4 = 4. * pow(lapse_auto_bh,6.) * mass * mass
	      * (1. + 3.*mass/rr) * (1. + 3.*mass/rr)
	      * (lapse_auto_rs + lapse_comp)
	      * pow(confo_tot,4.) / 3. / pow(rr,4.) ;
	    tmpl4.std_spectral_base() ;
	    tmpl4.annule_domain(0) ;
	    tmpl4.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    Scalar tmpl5 = -2. * pow(lapse_auto_bh,5.) * mass
	      * (2. + 10.*mass/rr + 9.*mass*mass/rr/rr)
	      * pow(confo_tot, 4.)
	      * (llshift + orb_rot_x*ll(2) - orb_rot_y*ll(1)) / pow(rr, 3.) ;
	    tmpl5.std_spectral_base() ;
	    tmpl5.annule_domain(0) ;
	    tmpl5.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    Scalar tmpl6 = 2. * lapse_auto_bh * lapse_auto_bh * mass
	      * dlap.dsdr() / rr ;
	    tmpl6.std_spectral_base() ;
	    tmpl6.annule_domain(0) ;
	    tmpl6.inc_dzpuis(1) ;  // dzpuis : 3 -> 4

	    Scalar tmpl7 = pow(lapse_auto_bh,4.) * mass * (3. + 8.*mass/rr)
	      * dlap / rr / rr ;
	    tmpl7.std_spectral_base() ;
	    tmpl7.annule_domain(0) ;
	    tmpl7.inc_dzpuis(2) ;  // dzpuis : 2 -> 4

	    Scalar tmpl8 = 4. * pow(lapse_auto_bh,7.) * mass * mass
	      * (2.*(lapse_auto_bh/lapse_tot - 1.) * pow(confo_tot,4.)
		 * (4. + 12.*mass/rr + 9.*mass*mass/rr/rr)
		 + 3.*(pow(confo_tot,4.) - 1.))
	      / 3. / pow(rr,4.) ;
	    tmpl8.std_spectral_base() ;
	    tmpl8.annule_domain(0) ;
	    tmpl8.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    source_lapse = tmpl1 + tmpl2 + tmpl3 + tmpl4 + tmpl5 + tmpl6
	      + tmpl7 + tmpl8 ;

	    source_lapse.annule_domain(0) ;
	    if (source_lapse.get_dzpuis() != 4) {
	      source_lapse.set_dzpuis(4) ;
	    }
	    source_lapse.std_spectral_base() ;
	    if (filter_r != 0) {
	      if (source_lapse.get_etat() != ETATZERO) {
	        source_lapse.filtre(filter_r) ;
	      }
	    }

	    bc_laps = bc_lapse() ;

	    lapse_m1.set_etat_qcq() ;

	    if (bc_lapse_nd) {
	      lapse_m1 = source_lapse.poisson_neumann(bc_laps, 0) ;
	    }
	    else {
	      lapse_m1 = source_lapse.poisson_dirichlet(bc_laps, 0) ;
	    }

	    // Re-construction of the lapse function
	    // -------------------------------------

	    lapse_auto_rs = lapse_m1 - 0.5 ;
	    lapse_auto_rs.annule_domain(0) ;
	    lapse_auto_rs.raccord(1) ;

	    lapse_auto = lapse_auto_rs + lapse_auto_bh ;
	    lapse_auto.annule_domain(0) ; // lapse_auto,_comp->0.5 (r->inf)
	    lapse_auto.raccord(1) ;       // lapse_tot -> 1 (r->inf)


	    //---------------------------------------------------------------//
	    //  Resolution of the Poisson equation for the conformal factor  //
	    //---------------------------------------------------------------//

	    // Source term
	    // -----------

	    Scalar tmpc1 = - 0.125 * (taij_quad_tot_rs + taij_quad_tot_rot)
	      / pow(confo_tot, 7.) ;
	    tmpc1.std_spectral_base() ; // dzpuis : 4

	    Scalar tmpc2 = 2. * lapse_auto_bh * lapse_auto_bh * mass
	      * dconf.dsdr() / rr ;
	    tmpc2.std_spectral_base() ;
	    tmpc2.inc_dzpuis(1) ;  // dzpuis : 3 -> 4

	    Scalar tmpc3 = pow(lapse_auto_bh,4.) * mass * (3.+8.*mass/rr)
	      * dconf / rr / rr ;
	    tmpc3.std_spectral_base() ;
	    tmpc3.inc_dzpuis(2) ;  // dzpuis : 2 -> 4

	    Scalar tmpc4 = pow(lapse_auto_bh,6.) * mass * mass * confo_tot
	      * ((1. - lapse_auto_bh*lapse_auto_bh/lapse_tot/lapse_tot)
		 * (4.+12.*mass/rr+9.*mass*mass/rr/rr)*pow(confo_tot,4.)
		 + 3.*(1.+2.*mass/rr)*(1.-pow(confo_tot,4.)))
	      / 3. / pow(rr,4.) ;
	    tmpc4.std_spectral_base() ;
	    tmpc4.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    source_confo = tmpc1 + tmpc2 + tmpc3 + tmpc4 ;

	    source_confo.annule_domain(0) ;
	    if (source_confo.get_dzpuis() != 4) {
	      source_confo.set_dzpuis(4) ;
	    }
	    source_confo.std_spectral_base() ;
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

	    Vector tmps1 = 2.
	      * contract(taij_tot_rs+taij_tot_rot, 1, dlappsi, 0)
	      / pow(confo_tot,6.) ;
	    tmps1.annule_domain(0) ;
	    tmps1.std_spectral_base() ;
	    // dzpuis = 4

	    Vector tmps2(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	        tmps2.set(i) = 2. * pow(lapse_auto_bh,3.) * mass
		  * (taij_tot_rs(i,1)*ll(1) + taij_tot_rs(i,2)*ll(2)
		     + taij_tot_rs(i,3)*ll(3))
		  / pow(confo_tot,6.) / rr / rr ;
	    }
	    tmps2.annule_domain(0) ;
	    tmps2.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps2.set(i)).inc_dzpuis(2) ;  // dzpuis : 2 -> 4

	    Vector tmps3(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	        tmps3.set(i) = 4. * pow(lapse_auto_bh,4.) * mass
		  * (2.+3.*mass/rr) * dlappsi(i) / 3. / lapse_tot / rr / rr ;
	    }
	    tmps3.annule_domain(0) ;
	    tmps3.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps3.set(i)).inc_dzpuis(2) ;  // dzpuis : 2 -> 4

	    Vector tmps4 = -4. * pow(lapse_auto_bh,6.) * mass
	      * (2.+3.*mass/rr) * (3.+2.*mass/rr) * ll
	      * (ll(1)*dlappsi(1) + ll(2)*dlappsi(2) + ll(3)*dlappsi(3))
	      / 3. / lapse_tot / rr / rr ;
	    tmps4.annule_domain(0) ;
	    tmps4.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps4.set(i)).inc_dzpuis(2) ;  // dzpuis : 2 -> 4

	    Vector tmps5(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	        tmps5.set(i) = 2. * lapse_auto_bh * lapse_auto_bh * mass
		  * dshift(i).dsdr() / rr ;
	    }
	    tmps5.annule_domain(0) ;
	    tmps5.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps5.set(i)).inc_dzpuis(1) ;  // dzpuis : 3 -> 4

	    Vector tmps6(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	        tmps6.set(i) = 2.*lapse_auto_bh*lapse_auto_bh*mass*ll(i)
		  *ddivshif/3./rr ;
	    }
	    tmps6.annule_domain(0) ;
	    tmps6.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps6.set(i)).inc_dzpuis(1) ;  // dzpuis : 3 -> 4

	    Vector tmps7(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	        tmps7.set(i) = - pow(lapse_auto_bh,4.) * mass * (3.+8.*mass/rr)
		  * (shift_auto_rs(i) + shift_comp(i)) / pow(rr, 3.) ;
	    }
	    tmps7.annule_domain(0) ;
	    tmps7.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps7.set(i)).inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    Vector tmps8(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	        tmps8.set(i) = pow(lapse_auto_bh,6.) * mass * ll(i) * llshift
		  * (25.+106.*mass/rr+96.*mass*mass/rr/rr) / 3./pow(rr, 3.) ;
	    }
	    tmps8.annule_domain(0) ;
	    tmps8.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps8.set(i)).inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    Vector tmps9(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	        tmps9.set(i) = - 8.* pow(lapse_auto_bh,7.) * mass * ll(i)
		  * (lapse_auto_rs + lapse_comp)
		  * (2.+10.*mass/rr+9.*mass*mass/rr/rr) / 3. / pow(rr, 3.) ;
	    }
	    tmps9.annule_domain(0) ;
	    tmps9.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps9.set(i)).inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    Vector tmps10(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	        tmps10.set(i) = - 2. * lapse_auto_bh * lapse_auto_bh
		  * mass * ll(i)
		  * (2. * divshift - lapse_auto_bh * lapse_auto_bh
		     * (3.+4.*mass/rr) * dllshift)
		  / rr / rr ;
	    }
	    tmps10.annule_domain(0) ;
	    tmps10.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps10.set(i)).inc_dzpuis(2) ;  // dzpuis : 2 -> 4

	    Vector tmps11(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	        tmps11.set(i) = pow(lapse_auto_bh,4.) * mass * (3.+8.*mass/rr)
		  * dshift(i) / rr / rr ;
	    }
	    tmps11.annule_domain(0) ;
	    tmps11.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps11.set(i)).inc_dzpuis(2) ;  // dzpuis : 2 -> 4

	    Vector tmps12(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	        tmps12.set(i) = lapse_auto_bh * lapse_auto_bh * mass
		  * llshift.deriv(i) / 3. / rr / rr ;
	    }
	    tmps12.annule_domain(0) ;
	    tmps12.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps12.set(i)).inc_dzpuis(2) ;  // dzpuis : 2 -> 4

	    Vector tmps13(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	        tmps13.set(i) = - 2. * pow(lapse_auto_bh,4.) * mass * mass
		  * ll(i) * dllshift / 3. / pow(rr, 3.) ;
	    }
	    tmps13.annule_domain(0) ;
	    tmps13.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps13.set(i)).inc_dzpuis(2) ;  // dzpuis : 2 -> 4

	    Vector tmps14 = 8. * pow(lapse_auto_bh,8.) * mass * mass
	      * (2.+3.*mass/rr) * (1. - lapse_auto_bh/lapse_tot) * ll
	      / 3. / pow(rr,4.) ;
	    tmps14.annule_domain(0) ;
	    tmps14.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps14.set(i)).inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    Vector tmps15(mp, CON, mp.get_bvect_cart()) ;
	    tmps15.set(1) = 2. * pow(lapse_auto_bh,4.) * mass
	      * ( 4.*(1+2.*mass/rr)
		  + 3.*mass*(1.-lapse_auto_bh/lapse_tot)/rr )
	      * orb_rot_y / 3. / pow(rr, 3.) ;
	    tmps15.set(2) = - 2. * pow(lapse_auto_bh,4.) * mass
	      * ( 4.*(1+2.*mass/rr)
		  + 3.*mass*(1.-lapse_auto_bh/lapse_tot)/rr )
	      * orb_rot_x / 3. / pow(rr, 3.) ;
	    tmps15.set(3) = 0. ;

	    tmps15.annule_domain(0) ;
	    tmps15.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps15.set(i)).inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    Vector tmps16(mp, CON, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	        tmps16.set(i) = 2. * pow(lapse_auto_bh,6.) * mass
		  * (orb_rot_x * ll(2) - orb_rot_y * ll(1)) * ll(i)
		  * ( 2.*(3.+4.*mass/rr)*(2.+5.*mass/rr)
		      +(5.+6.*mass/rr)*mass*(1.-lapse_auto_bh/lapse_tot)/rr)
		  / 3. /pow(rr, 3.) ;
	    }
	    tmps16.annule_domain(0) ;
	    tmps16.std_spectral_base() ;
	    for (int i=1; i<=3; i++)
	        (tmps16.set(i)).inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    source_shift = tmps1 + tmps2 + tmps3 + tmps4 + tmps5 + tmps6
	      + tmps7 + tmps8 + tmps9 + tmps10 + tmps11 + tmps12 + tmps13
	      + tmps14 + tmps15 + tmps16 ;

	    source_shift.annule_domain(0) ;
	    source_shift.std_spectral_base() ;

	    if (filter_r_s!= 0) {
	      for (int i=1; i<=3; i++) {
	        if (source_shift(i).get_etat() != ETATZERO)
		  source_shift.set(i).filtre(filter_r_s) ;
	      }
	    }

	    if (filter_p_s != 0) {
	      for (int i=1; i<=3; i++) {
	        if (source_shift(i).get_etat() != ETATZERO) {
		  for (int l=1; l<nz; l++) {
		    source_shift.set(i).filtre_phi(filter_p_s, l) ;
		  }
		}
	      }
	    }

	    for (int i=1; i<=3; i++) {
	      if (source_shift(i).dz_nonzero()) {
	        assert( source_shift(i).get_dzpuis() == 4 ) ;
	      }
	      else {
	        (source_shift.set(i)).set_dzpuis(4) ;
	      }
	    }

	    Tenseur source_p(mp, 1, CON, mp.get_bvect_cart()) ;
	    source_p.set_etat_qcq() ;
	    for (int i=0; i<3; i++) {
	      source_p.set(i) = Cmp(source_shift(i+1)) ;
	    }

	    Tenseur resu_p(mp, 1, CON, mp.get_bvect_cart()) ;
	    resu_p.set_etat_qcq() ;

	    for (int i=0; i<3; i++) {
	      resu_p.set(i) = shift_auto_rs(i+1) ;
	    }

	    // Boundary condition
	    bc_shif_x = bc_shift_x(omega_orb, omega_spin, y_rot) ;
	    bc_shif_y = bc_shift_y(omega_orb, omega_spin, x_rot) ;
	    bc_shif_z = bc_shift_z() ;

	    poisson_vect_frontiere(1./3., source_p, resu_p,
				   bc_shif_x, bc_shif_y, bc_shif_z,
				   0, precis, 20) ;

	    for (int i=1; i<=3; i++) {
	      shift_auto_rs.set(i) = resu_p(i-1) ;
	    }

	    shift_auto_rs.std_spectral_base() ;
	    shift_auto_rs.annule_domain(0) ;

	    shift_auto = shift_auto_rs + shift_auto_bh ;
	    shift_auto.std_spectral_base() ;
	    shift_auto.annule_domain(0) ;

	}  // End of Kerr-Schild
	else { // Isotropic coordinates with the maximal slicing

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

	    Scalar r_are(mp) ;
	    r_are = r_coord(bc_lapse_nd, bc_lapse_fs) ;
	    r_are.std_spectral_base() ;

	    Scalar lapbh_iso(mp) ;
	    lapbh_iso = sqrt(1. - 2.*mass/r_are/rr
			     + cc*cc*pow(mass/r_are/rr,4.)) ;
	    lapbh_iso.std_spectral_base() ;
	    lapbh_iso.annule_domain(0) ;
	    lapbh_iso.raccord(1) ;

	    Scalar dlapbh_iso(mp) ;
	    dlapbh_iso = mass/r_are/rr - 2.*cc*cc*pow(mass/r_are/rr,4.) ;
	    dlapbh_iso.std_spectral_base() ;
	    dlapbh_iso.annule_domain(0) ;
	    dlapbh_iso.raccord(1) ;

	    //-------------------------------------------------------------//
	    //  Resolution of the Poisson equation for the lapse function  //
	    //-------------------------------------------------------------//

	    // Source term
	    // -----------

	    Scalar tmpl1 = lapse_tot * taij_quad_auto / pow(confo_tot, 8.) ;
	    tmpl1.std_spectral_base() ;  // dzpuis = 4

	    Scalar tmpl2 = -2.*( lapse_auto_rs.deriv(1)
				 *(confo_auto_rs.deriv(1)+d_confo_comp(1))
				 + lapse_auto_rs.deriv(2)
				 *(confo_auto_rs.deriv(2)+d_confo_comp(2))
				 + lapse_auto_rs.deriv(3)
				 *(confo_auto_rs.deriv(3)+d_confo_comp(3))
				 + d_lapse_comp(1) * confo_auto_rs.deriv(1)
				 + d_lapse_comp(2) * confo_auto_rs.deriv(2)
				 + d_lapse_comp(3) * confo_auto_rs.deriv(3)
				 )
	      / confo_tot ;

	    tmpl2.std_spectral_base() ;  // dzpuis = 4

	    Scalar tmpl3 = - sqrt(r_are) * dlap * (lapbh_iso - 1.)
	      / rr / confo_tot
	      - 2. * dlconfo * dlapbh_iso / rr ;
	    tmpl3.std_spectral_base() ;
	    tmpl3.annule_domain(0) ;

	    tmpl3.inc_dzpuis(2) ;  // dzpuis : 2 -> 4

	    Scalar tmpl4 = dlapbh_iso * (lapbh_iso - 1.)
	      * (1. - sqrt(r_are)/confo_tot) / rr / rr
	      - 6. * cc * cc * pow(mass,4.) * lapbh_iso
	      * (1. - pow(confo_tot,4.)*lapbh_iso/r_are/r_are/lapse_tot)
	      / pow(r_are,4.) / pow(rr,6.) ;
	    tmpl4.std_spectral_base() ;
	    tmpl4.annule_domain(0) ;

	    tmpl4.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    Scalar tmpl5 = (lapse_comp+0.5) * taij_quad_comp
	      * (pow(confo_tot/(confo_comp+0.5),4.)*(lapse_comp+0.5)
		 /lapse_tot - 1.)
	      / pow(confo_comp+0.5,8.) ;
	    tmpl5.std_spectral_base() ;
	    tmpl5.annule_domain(0) ;

	    Scalar tmpl6 = 2. * (1. - (confo_comp+0.5)/confo_tot)
	      * (d_lapse_comp(1)*d_confo_comp(1)
		 + d_lapse_comp(2)*d_confo_comp(2)
		 + d_lapse_comp(3)*d_confo_comp(3)) / (confo_comp+0.5) ;
	    tmpl6.std_spectral_base() ;
	    tmpl6.annule_domain(0) ;  // dzpuis = 4

	    source_lapse = tmpl1 + tmpl2 + tmpl3 + tmpl4 + tmpl5 + tmpl6 ;

	    source_lapse.annule_domain(0) ;
	    if (source_lapse.get_dzpuis() != 4) {
	      source_lapse.set_dzpuis(4) ;
	    }
	    source_lapse.std_spectral_base() ;
	    if (filter_r != 0) {
	      if (source_lapse.get_etat() != ETATZERO) {
	        source_lapse.filtre(filter_r) ;
	      }
	    }

	    bc_laps = bc_lapse() ;

	    lapse_m1.set_etat_qcq() ;

	    if (bc_lapse_nd) {
	      lapse_m1 = source_lapse.poisson_neumann(bc_laps, 0) ;
	    }
	    else {
	      lapse_m1 = source_lapse.poisson_dirichlet(bc_laps, 0) ;
	    }

	    // Re-construction of the lapse function
	    // -------------------------------------

	    lapse_auto_rs = lapse_m1 - 0.5 ;
	    lapse_auto_rs.annule_domain(0) ;
	    lapse_auto_rs.raccord(1) ;

	    lapse_auto = lapse_auto_rs + lapse_auto_bh ;
	    lapse_auto.annule_domain(0) ; // lapse_auto,_comp->0.5 (r->inf)
	    lapse_auto.raccord(1) ;       // lapse_tot -> 1 (r->inf)


	    //---------------------------------------------------------------//
	    //  Resolution of the Poisson equation for the conformal factor  //
	    //---------------------------------------------------------------//

	    // Source term
	    // -----------

	    Scalar tmpc1 = - 0.125 * taij_quad_auto / pow(confo_tot, 7.) ;
	    tmpc1.std_spectral_base() ;
	    // dzpuis = 4

	    Scalar tmpc2 = 0.75 * cc * cc * pow(mass,4.)
	      * (pow(r_are,2.5)
		 - pow(confo_tot,5.)*lapbh_iso*lapbh_iso/lapse_tot/lapse_tot)
	      / pow(r_are*rr,6.) ;
	    tmpc2.std_spectral_base() ;
	    tmpc2.annule_domain(0) ;

	    tmpc2.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	    Scalar tmpc3 = 0.125 * taij_quad_comp
	      * (1. - pow(confo_tot/(confo_comp+0.5),5.)
		 *pow((lapse_comp+0.5)/lapse_tot,2.))
	      / pow(confo_comp+0.5, 7.) ;
	    tmpc3.std_spectral_base() ;
	    // dzpuis = 4

	    source_confo = tmpc1 + tmpc2 + tmpc3 ;

	    source_confo.annule_domain(0) ;
	    if (source_confo.get_dzpuis() != 4) {
	      source_confo.set_dzpuis(4) ;
	    }
	    source_confo.std_spectral_base() ;
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
	      dlapconf.set(i) = lapse_auto_rs.deriv(i)
		- 6. * lapse_tot * confo_auto_rs.deriv(i)
		/ confo_tot ;
	    }

	    dlapconf.std_spectral_base() ;

	    Vector tmps1 = 2. * contract(taij_auto_rs, 1, dlappsi, 0)
	      / pow(confo_tot, 6.)
	      + 2. * contract(taij_comp, 1, dlapconf, 0)
	      * (lapse_comp+0.5) / lapse_tot / pow(confo_comp+0.5, 6.) ;
	    tmps1.annule_domain(0) ;
	    tmps1.std_spectral_base() ;
	    // dzpuis = 4

	    Vector tmps2(mp, CON, mp.get_bvect_cart()) ;
	    tmps2.set_etat_qcq() ;
	    for (int i=1; i<=3; i++) {
	      tmps2.set(i) = 2. * (dlapbh_iso
				   - 3.*lapse_tot*sqrt(r_are)
				   * (lapbh_iso - 1.) / confo_tot)
		* (taij_tot_rs(i,1)*ll(1) + taij_tot_rs(i,2)*ll(2)
		   + taij_tot_rs(i,3)*ll(3)) / pow(confo_tot,6.) / rr ;
	    }
	    tmps2.annule_domain(0) ;
	    tmps2.std_spectral_base() ;
	    for (int i=1; i<=3; i++) {
	        (tmps2.set(i)).inc_dzpuis(2) ;  // dzpuis : 2 -> 4
	    }

	    Vector tmps3(mp, CON, mp.get_bvect_cart()) ;
	    tmps3.set_etat_qcq() ;
	    for (int i=1; i<=3; i++) {
	        tmps3.set(i) = 2. * cc * mass * mass * lapbh_iso
		  * (dlappsi(i) - 3.*ll(i)*(ll(1)*dlappsi(1)
					    + ll(2)*dlappsi(2)
					    + ll(3)*dlappsi(3)))
		  / lapse_tot / pow(r_are*rr,3.) ;
	    }
	    tmps3.annule_domain(0) ;
	    tmps3.std_spectral_base() ;
	    for (int i=1; i<=3; i++) {
	        (tmps3.set(i)).inc_dzpuis(2) ;  // dzpuis : 2 -> 4
	    }

	    Vector tmps4 = - 4. * cc * mass * mass
	      * (dlapbh_iso * (lapbh_iso/lapse_tot - 1.)
		 + 3. * lapbh_iso * (lapbh_iso - 1.)
		 * (1. - sqrt(r_are)/confo_tot))
	      * ll / rr / pow(r_are*rr,3.) ;
	    tmps4.annule_domain(0) ;
	    tmps4.std_spectral_base() ;
	    for (int i=1; i<=3; i++) {
	        (tmps4.set(i)).inc_dzpuis(4) ;  // dzpuis : 0 -> 4
	    }

	    Vector dlappsi_comp(mp, COV, mp.get_bvect_cart()) ;
	    for (int i=1; i<=3; i++) {
	      dlappsi_comp.set(i) = ((lapse_comp+0.5)/lapse_tot - 1.)
		* d_lapse_comp(i)
		- 6. * (lapse_comp+0.5) * ((confo_comp+0.5)/confo_tot - 1.)
		* d_confo_comp(i) / (confo_comp+0.5) ;
	    }

	    dlappsi_comp.std_spectral_base() ;

	    Vector tmps5 = 2. * contract(taij_comp, 1, dlappsi_comp, 0)
	      / pow(confo_comp+0.5, 6.) ;
	    tmps5.std_spectral_base() ;

	    source_shift = tmps1 + tmps2 + tmps3 + tmps4 + tmps5 ;

	    source_shift.annule_domain(0) ;
	    source_shift.std_spectral_base() ;

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

	    for (int i=1; i<=3; i++) {
	      if (source_shift(i).dz_nonzero()) {
	        assert( source_shift(i).get_dzpuis() == 4 ) ;
	      }
	      else {
	        (source_shift.set(i)).set_dzpuis(4) ;
	      }
	    }

	    Tenseur source_p(mp, 1, CON, mp.get_bvect_cart()) ;
	    source_p.set_etat_qcq() ;
	    for (int i=0; i<3; i++) {
	      source_p.set(i) = Cmp(source_shift(i+1)) ;
	    }

	    Tenseur resu_p(mp, 1, CON, mp.get_bvect_cart()) ;
	    resu_p.set_etat_qcq() ;

	    for (int i=0; i<3; i++) {
	      resu_p.set(i) = shift_auto_rs(i+1) ;
	    }

	    // Boundary condition
	    bc_shif_x = bc_shift_x(omega_orb, omega_spin, y_rot) ;
	    bc_shif_y = bc_shift_y(omega_orb, omega_spin, x_rot) ;
	    bc_shif_z = bc_shift_z() ;

	    poisson_vect_frontiere(1./3., source_p, resu_p,
				   bc_shif_x, bc_shif_y, bc_shif_z,
				   0, precis, 20) ;

	    for (int i=1; i<=3; i++) {
	      shift_auto_rs.set(i) = resu_p(i-1) ;
	    }

	    shift_auto_rs.std_spectral_base() ;
	    shift_auto_rs.annule_domain(0) ;

	    shift_auto = shift_auto_rs + shift_auto_bh ;
	    shift_auto.std_spectral_base() ;
	    shift_auto.annule_domain(0) ;

	}  // End of isotropic

	//------------------------------------------------//
	//  Relative difference in the metric quantities  //
	//------------------------------------------------//

	// Difference is calculated only outside the inner boundary.

	Tbl tdiff_lapse = diffrel(lapse_auto_rs, lapse_jm1) ;
	tdiff_lapse.set(0) = 0. ;
	cout << "Relative difference in the lapse function   : " << endl ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_lapse(l) << "  " ;
	}
	cout << endl ;

	diff_lapse = tdiff_lapse(1) ;
	for (int l=2; l<nz; l++) {
	    diff_lapse += tdiff_lapse(l) ;
	}
	diff_lapse /= nz ;

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
	des_profile( lapse_auto_rs, 0., 10.,
		     M_PI/2., 0., "Residual lapse function of BH",
		     "Lapse (theta=pi/2, phi=0)" ) ;

	des_profile( lapse_auto_bh, 0., 10.,
		     M_PI/2., 0., "Analytic lapse function of BH",
		     "Lapse (theta=pi/2, phi=0)" ) ;

	des_profile( lapse_auto, 0., 10.,
		     M_PI/2., 0., "Self lapse function of BH",
		     "Lapse (theta=pi/2, phi=0)" ) ;

	des_profile( confo_auto_rs, 0., 10.,
		     M_PI/2., 0., "Residual conformal factor of BH",
		     "Confo (theta=pi/2, phi=0)" ) ;

	des_profile( confo_auto_bh, 0., 10.,
		     M_PI/2., 0., "Analytic conformal factor of BH",
		     "Confo (theta=pi/2, phi=0)" ) ;

	des_profile( confo_auto, 0., 10.,
		     M_PI/2., 0., "Self conformal factor of BH",
		     "Confo (theta=pi/2, phi=0)" ) ;
	*/

    } // End of main loop

    //====================================//
    //          End of iteration          //
    //====================================//

}
