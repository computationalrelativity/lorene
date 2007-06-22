/*
 *  Method of class Star_bhns to compute the velocity scalar potential $\psi$
 *  by solving the continuity equation.
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

char star_bhns_vel_pot_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:33:14  k_taniguchi
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
#include "eos.h"
#include "param.h"
#include "cmp.h"
#include "tenseur.h"
#include "utilitaires.h"
#include "unites.h"

// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ;

double Star_bhns::velo_pot_bhns(const double& mass_bh, const double& sepa,
				bool kerrschild, int mermax, double precis,
				double relax) {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    int nzm1 = mp.get_mg()->get_nzone() - 1 ;

    //--------------------------------
    // Specific relativistic enthalpy      ---> hhh
    //--------------------------------

    Scalar hhh = exp(ent) ;
    hhh.std_spectral_base() ;

    //---------------------------------------------------
    // Computation of W^i = psi4 h Gamma_n B^i/lapse_tot
    // See Eq (62) from Gourgoulhon et al. (2001)
    //---------------------------------------------------

    Vector www = hhh * gam_euler * bsn * psi4 ;

    www.change_triad( mp.get_bvect_cart() ) ;  // components on the mapping
                                               // Cartesian basis

    //-------------------------------------------------
    // Constant value of W^i at the center of the star
    //-------------------------------------------------

    Vector v_orb(mp, CON, mp.get_bvect_cart()) ;
    v_orb.set_etat_qcq() ;

    for (int i=1; i<=3; i++) {
        v_orb.set(i) = www(i).val_grid_point(0, 0, 0, 0) ;
    }

    v_orb.set_triad( *(www.get_triad()) ) ;
    v_orb.std_spectral_base() ;

    //-------------------------------------------------
    // Source and coefficients a,b for poisson_compact (independent from psi0)
    //-------------------------------------------------

    Scalar dndh_log = eos.der_nbar_ent(ent, nzet) ;

    // In order to avoid any division by zero in the computation of zeta_h
    //  the value of dndh_log is set to 1 in the external domains:
    for (int l=nzet; l<=nzm1; l++) {
        dndh_log.set_domain(l) = 1. ;
    }

    Scalar zeta_h( ent / dndh_log ) ;
    zeta_h.std_spectral_base() ;

    double erreur ;

    Scalar source(mp) ;
    Vector bb(mp, CON, mp.get_bvect_cart()) ;
    Metric_flat flat_spher( mp.flat_met_spher() ) ;

    if (kerrschild) {

        double mass = ggrav * mass_bh ;

	Scalar xx(mp) ;
	xx = mp.x ;
	xx.std_spectral_base() ;
	Scalar yy(mp) ;
	yy = mp.y ;
	yy.std_spectral_base() ;
	Scalar zz(mp) ;
	zz = mp.z ;
	zz.std_spectral_base() ;

	double yns = mp.get_ori_y() ;

	Scalar rbh(mp) ;
	rbh = sqrt( (xx+sepa)*(xx+sepa) + (yy+yns)*(yy+yns) + zz*zz ) ;
	rbh.std_spectral_base() ;

	Vector ll(mp, CON, mp.get_bvect_cart()) ;
	ll.set_etat_qcq() ;
	ll.set(1) = (xx+sepa) / rbh ;
	ll.set(2) = (yy+yns) / rbh ;
	ll.set(3) = zz / rbh ;
	ll.std_spectral_base() ;

	Scalar msr(mp) ;
	msr = mass / rbh ;
	msr.std_spectral_base() ;

	Scalar lap_bh(mp) ;
	lap_bh = 1./sqrt(1.+2.*msr) ;
	lap_bh.std_spectral_base() ;

	Scalar lnlappsi(mp) ;
	lnlappsi = log( lapse_tot * confo_tot * confo_tot ) ;
	lnlappsi.std_spectral_base() ;

	bb = (1. - zeta_h) * ent.derive_con(flat_spher)
	  + zeta_h * lnlappsi.derive_con(flat_spher) ;

	Vector dentmb(mp, COV, mp.get_bvect_cart()) ;
	dentmb.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	    dentmb.set(i) = ent.deriv(i)
	      - (d_lapse_auto(i) + d_lapse_comp(i)) / lapse_tot
	      - 2. * (d_confo_auto(i) + d_confo_comp(i)) / confo_tot ;
	}
	dentmb.std_spectral_base() ;

	Scalar lldpsi0(mp) ;
	lldpsi0 = ll(1)*psi0.deriv(1) + ll(2)*psi0.deriv(2)
	  + ll(3)*psi0.deriv(3) ;
	lldpsi0.std_spectral_base() ;

	Scalar lldlldpsi0(mp) ;
	lldlldpsi0 = ll(1)*lldpsi0.deriv(1) + ll(2)*lldpsi0.deriv(2)
	  + ll(3)*lldpsi0.deriv(3) ;
	lldlldpsi0.std_spectral_base() ;

	Scalar llvorb(mp) ;
	llvorb = ll(1)*v_orb(1) + ll(2)*v_orb(2) + ll(3)*v_orb(3) ;
	llvorb.std_spectral_base() ;

	Scalar lldent(mp) ;
	lldent = ll(1)*ent.deriv(1) + ll(2)*ent.deriv(2)
	  + ll(3)*ent.deriv(3) ;
	lldent.std_spectral_base() ;

	Scalar lldlnlappsi(mp) ;
	lldlnlappsi = ll(1)*((d_lapse_auto(1)+d_lapse_comp(1))/lapse_tot
			     +2.*(d_confo_auto(1)+d_confo_comp(1))/confo_tot)
	  + ll(2)*((d_lapse_auto(2)+d_lapse_comp(2))/lapse_tot
		   +2.*(d_confo_auto(2)+d_confo_comp(2))/confo_tot)
	  + ll(3)*((d_lapse_auto(3)+d_lapse_comp(3))/lapse_tot
		   +2.*(d_confo_auto(3)+d_confo_comp(3))/confo_tot) ;

	lldlnlappsi.std_spectral_base() ;

	source = contract(www - v_orb, 0, ent.derive_cov(flat), 0)
	  + zeta_h*(contract(v_orb, 0, dentmb, 0)
		    +contract(www/gam_euler, 0, gam_euler.derive_cov(flat), 0)
		    )
	  + 2. * zeta_h * lap_bh * lap_bh * msr * lldlldpsi0
	  + (2.*lap_bh*lap_bh*msr*( (1.-zeta_h)*lldent + zeta_h*lldlnlappsi )
	     + zeta_h*pow(lap_bh,4.)*msr*(3.+8.*msr)/rbh)*(lldpsi0 + llvorb)
	  + 2.*zeta_h*psi4*hhh*gam_euler*pow(lap_bh,3.)
	  *msr*(1.+3.*msr)/rbh ;

	source.annule(nzet, nzm1) ;

    }  // End of Kerr-Schild
    else {  // Isotropic coordinates with the maximal slicing

        Scalar lnlappsi(mp) ;
	lnlappsi = log( lapse_tot * confo_tot * confo_tot ) ;
	lnlappsi.std_spectral_base() ;

	bb = (1. - zeta_h) * ent.derive_con(flat_spher)
	  + zeta_h * lnlappsi.derive_con(flat_spher) ;

	Vector dentmb(mp, COV, mp.get_bvect_cart()) ;
	dentmb.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	    dentmb.set(i) = ent.deriv(i)
	      - (d_lapse_auto(i) + d_lapse_comp(i)) / lapse_tot
	      - 2. * (d_confo_auto(i) + d_confo_comp(i)) / confo_tot ;
	}
	dentmb.std_spectral_base() ;

	source = contract(www - v_orb, 0, ent.derive_cov(flat), 0)
	  +zeta_h*(contract(v_orb, 0, dentmb, 0)
		   + contract(www/gam_euler, 0, gam_euler.derive_cov(flat), 0)
		   ) ;

	source.annule(nzet, nzm1) ;

    }  // End of Isotropic coordinates


    //----------------------------------------------------
    // Resolution by means of Map_radial::poisson_compact
    //----------------------------------------------------

    Param par ;
    int niter ;
    par.add_int(mermax) ;
    par.add_double(precis, 0) ;
    par.add_double(relax, 1) ;
    par.add_int_mod(niter) ;

    if (psi0.get_etat() == ETATZERO) {
        psi0 = 0. ;
    }

    Cmp source_cmp(source) ;
    Cmp zeta_h_cmp(zeta_h) ;
    Cmp psi0_cmp(psi0) ;

    Tenseur bb_cmp(mp, 1, CON, mp.get_bvect_spher()) ;
    bb_cmp.set_etat_qcq() ;
    Cmp bb_cmp1(bb(1)) ;
    Cmp bb_cmp2(bb(2)) ;
    Cmp bb_cmp3(bb(3)) ;
    bb_cmp.set(0) = bb_cmp1 ;
    bb_cmp.set(1) = bb_cmp2 ;
    bb_cmp.set(2) = bb_cmp3 ;

    source_cmp.va.ylm() ;

    mp.poisson_compact(source_cmp, zeta_h_cmp, bb_cmp, par, psi0_cmp) ;

    psi0 = psi0_cmp ;

    //-----------------------
    // Check of the solution
    //-----------------------

    Scalar bb_dpsi0 = contract(bb, 0, psi0.derive_cov(flat_spher), 0) ;

    Scalar oper = zeta_h * psi0.laplacian() + bb_dpsi0 ;

    source.set_spectral_va().ylm_i() ;

    erreur = diffrel(oper, source)(0) ;

    cout << "Check of the resolution of the continuity equation : "
	 << endl ;
    cout << "            norme(source) : " << norme(source)(0)
	 << "    diff oper/source : " << erreur << endl ;

    //--------------------------
    // Computation of grad(psi)
    //--------------------------

    for (int i=1; i<=3; i++)
        d_psi.set(i) = (psi0.derive_cov(flat))(i) + v_orb(i) ;


    // C^1 continuation of d_psi outside the star
    //  (to ensure a smooth enthalpy field accross the stellar surface)
    // ----------------------------------------------------------------

    d_psi.annule(nzet, nzm1) ;

    for (int i=1; i<=3; i++) {
        Cmp d_psi_i( d_psi(i) ) ;
	d_psi_i = raccord_c1(d_psi_i, nzet) ;
	d_psi.set(i) = d_psi_i ;
    }

    return erreur ;

}
