/*
 *  Method of class Et_bin_bhns_extr to compute the velocity scalar
 *  potential $\psi$ by solving the continuity equation
 *  in the Kerr-Schild background metric or in the conformally flat one
 *  with extreme mass ratio
 *
 *    (see file et_bin_bhns_extr.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004-2005 Keisuke Taniguchi
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

char et_bin_bhns_extr_velpot_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:52:55  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2005/02/28 23:17:18  k_taniguchi
 * Modification to include the case of the conformally flat background metric
 *
 * Revision 1.1  2004/11/30 20:51:56  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
//#include <math.h>

// Lorene headers
#include "et_bin_bhns_extr.h"
#include "scalar.h"
#include "metrique.h"
#include "etoile.h"
#include "eos.h"
#include "param.h"
#include "coord.h"
#include "unites.h"

// Local prototype
namespace Lorene {
Cmp raccord_c1(const Cmp& uu, int l1) ;

double Et_bin_bhns_extr::velocity_pot_extr(const double& mass,
					   const double& sepa,
					   int mermax, double precis,
					   double relax) {

  using namespace Unites ;

    if (kerrschild) {

        if (eos.identify() == 5 || eos.identify() == 4 ||
	    eos.identify() == 3) {

	    cout << "Etoile_bin::velocity_pot_extr" << endl ;
	    cout << "The code has not prepared for this kind of EOS!!!"
		 << endl;
	    abort() ;

	}  // End of strange stars case
	else {

	    int nzm1 = mp.get_mg()->get_nzone() - 1 ;

	    //--------------------------------
	    // Specific relativistic enthalpy		    ---> hhh
	    //--------------------------------

	    Tenseur hhh = exp(unsurc2 * ent) ;  // = 1 at the Newtonian limit
	    hhh.set_std_base() ;

	    //--------------------------------------------
	    // Computation of W^i = - A^2 h Gamma_n B^i/N
	    //--------------------------------------------

	    Tenseur www = - a_car * hhh * gam_euler * bsn ;

	    www.change_triad( mp.get_bvect_cart() ) ;
	                        // components on the mapping Cartesian basis

	    //-------------------------------------------------
	    // Constant value of W^i at the center of the star
	    //-------------------------------------------------

	    Tenseur v_orb(mp, 1, CON, mp.get_bvect_cart()) ;

	    v_orb.set_etat_qcq() ;
	    for (int i=0; i<3; i++) {
	        v_orb.set(i) = www(i)(0, 0, 0, 0) ;
	    }

	    v_orb.set_triad( *(www.get_triad()) ) ;
	    v_orb.set_std_base() ;

	    // Some auxiliary terms
	    // --------------------

	    const Coord& xx = mp.x ;
	    const Coord& yy = mp.y ;
	    const Coord& zz = mp.z ;

	    Tenseur r_bh(mp) ;
	    r_bh.set_etat_qcq() ;
	    r_bh.set() = pow( (xx+sepa)*(xx+sepa) + yy*yy + zz*zz, 0.5) ;
	    r_bh.set_std_base() ;

	    Tenseur msr(mp) ;
	    msr = ggrav * mass / r_bh ;
	    msr.set_std_base() ;

	    Tenseur lapse_bh2(mp) ;  // lapse_bh * lapse_bh
	    lapse_bh2 = 1. / (1.+2.*msr) ;
	    lapse_bh2.set_std_base() ;

	    Tenseur lapse_bh3(mp) ;  // lapse_bh * lapse_bh * lapse_bh 
	    lapse_bh3 = 1./pow(1.+2.*msr, 1.5) ;
	    lapse_bh3.set_std_base() ;

	    Tenseur xx_con(mp, 1, CON, mp.get_bvect_cart()) ;
	    xx_con.set_etat_qcq() ;
	    xx_con.set(0) = xx + sepa ;
	    xx_con.set(1) = yy ;
	    xx_con.set(2) = zz ;
	    xx_con.set_std_base() ;

	    Tenseur xsr_con(mp, 1, CON, mp.get_bvect_cart()) ;
	    xsr_con = xx_con / r_bh ;
	    xsr_con.set_std_base() ;

	    Tenseur xx_cov(mp, 1, COV, mp.get_bvect_cart()) ;
	    xx_cov.set_etat_qcq() ;
	    xx_cov.set(0) = xx + sepa ;
	    xx_cov.set(1) = yy ;
	    xx_cov.set(2) = zz ;
	    xx_cov.set_std_base() ;

	    Tenseur xsr_cov(mp, 1, COV, mp.get_bvect_cart()) ;
	    xsr_cov = xx_cov / r_bh ;
	    xsr_cov.set_std_base() ;

	    // X_i/r_bh in the spherical coordinate
	    const Coord& rr = mp.r ;
	    const Coord& st = mp.sint ;
	    const Coord& ct = mp.cost ;
	    const Coord& sp = mp.sinp ;
	    const Coord& cp = mp.cosp ;

	    Tenseur rr_spher_cov(mp, 1, COV, mp.get_bvect_spher()) ;
	    rr_spher_cov.set_etat_qcq() ;
	    rr_spher_cov.set(0) = rr + sepa*st*cp ;
	    rr_spher_cov.set(1) = sepa*ct*cp ;
	    rr_spher_cov.set(2) = - sepa*sp ;
	    rr_spher_cov.set_std_base() ;

	    Tenseur xsr_spher_cov(mp, 1, COV, mp.get_bvect_spher()) ;
	    xsr_spher_cov = rr_spher_cov / r_bh ;
	    xsr_spher_cov.set_std_base() ;

	    // l^j D_j H_ent
	    Tenseur ldent = flat_scalar_prod(xsr_con, ent.gradient()) ;

	    // l^j D_j beta_auto
	    Tenseur ldbeta = flat_scalar_prod(xsr_con, beta_auto.gradient()) ;

	    // l^j D_j psi0
	    Tenseur ldpsi0 = flat_scalar_prod(xsr_con, psi0.gradient()) ;

	    // l^i D_i (l^j D_j psi0)
	    Tenseur ldldpsi0 = flat_scalar_prod(xsr_con, ldpsi0.gradient()) ;

	    //-------------------------------------------------
	    // Source and coefficients a,b for poisson_compact
	    // (idenpendent from psi0)
	    //-------------------------------------------------

	    Cmp dndh_log = eos.der_nbar_ent(ent(), nzet) ;

	    // In order to avoid any division by zero in the computation of
	    //  zeta_h the value of dndh_log is set to 1
	    //  in the external domains:
	    for (int l=nzet; l <= nzm1; l++) {
	        dndh_log.set(l) = 1 ;
	    }

	    double erreur ;

	    Tenseur zeta_h( ent() / dndh_log ) ;
	    zeta_h.set_std_base() ;

	    Tenseur tmp_zeta = 1 - unsurc2 * zeta_h ;
	    tmp_zeta.set_std_base() ;

	    Tenseur bb = tmp_zeta * (ent.gradient_spher()
				     - 2.*lapse_bh2 * msr * ldent
				     * xsr_spher_cov)
	      + unsurc2 * zeta_h * (beta_auto.gradient_spher()
				    - 2.*lapse_bh2 * msr * ldbeta
				    * xsr_spher_cov)
	      - unsurc2 * 2. * zeta_h * lapse_bh2 * lapse_bh2 * msr / r_bh
	      * (1.+4.*msr) * xsr_spher_cov ;

	    Tenseur entmb = ent - beta_auto ;

	    Tenseur source = flat_scalar_prod(www - v_orb, ent.gradient())
	      + unsurc2*zeta_h*( flat_scalar_prod(v_orb, entmb.gradient())
				 + flat_scalar_prod(www, gam_euler.gradient())
				 / gam_euler )
	      + 2.*lapse_bh2 * msr * flat_scalar_prod(xsr_cov, v_orb)
	      * flat_scalar_prod(xsr_con, ent.gradient())
	      + unsurc2 * 2. * zeta_h
	      * (lapse_bh2*msr*(ldldpsi0
				- flat_scalar_prod(xsr_cov, v_orb)
				* (flat_scalar_prod(xsr_con, entmb.gradient())
				   - lapse_bh2 * (1.+4.*msr) / r_bh))
		 + a_car * hhh * gam_euler * lapse_bh3 * msr * (1.+3.*msr)
		 / r_bh) ;

	    source.annule(nzet, nzm1) ;

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
	        psi0.set_etat_qcq() ;
		psi0.set() = 0 ;
	    }

	    source.set().va.ylm() ;

	    mp.poisson_compact(source(), zeta_h(), bb, par, psi0.set() ) ;

	    //-----------------------
	    // Check of the solution
	    //-----------------------

	    Tenseur bb_dpsi0 = flat_scalar_prod( bb, psi0.gradient_spher() ) ;

	    Cmp oper = zeta_h() * psi0().laplacien() + bb_dpsi0() ;

	    source.set().va.ylm_i() ;

	    erreur = diffrel(oper, source())(0) ;

	    cout << "Check of the resolution of the continuity equation : "
		 << endl ;
	    cout << "norme(source) : " << norme(source())(0)
		 << "    diff oper/source : " << erreur << endl ;

	    //--------------------------
	    // Computation of grad(psi)
	    //--------------------------

	    // The computation is done component by component because
	    // psi0.gradient() is a covariant vector, whereas v_orb is a
	    // contravariant one.

	    d_psi.set_etat_qcq() ;

	    for (int i=0; i<3; i++) {
	        d_psi.set(i) = (psi0.gradient())(i) + v_orb(i) ;
	    }

	    d_psi.set_triad( *(v_orb.get_triad()) ) ;

	    // C^1 continuation of d_psi outside the star
	    //  (to ensure a smooth enthalpy field accross the stellar surface)
	    // ----------------------------------------------------------------

	    d_psi.annule(nzet, nzm1) ;
	    for (int i=0; i<3; i++) {
	        d_psi.set(i) = raccord_c1(d_psi(i), nzet) ;
	    }

	    assert( d_psi.get_triad() == &(mp.get_bvect_cart()) ) ;

	    d_psi.change_triad(ref_triad) ;

	    return erreur ;

	}

    }
    else {

        if (eos.identify() == 5 || eos.identify() == 4 ||
	    eos.identify() == 3) {

	    cout << "Etoile_bin::velocity_pot_extr" << endl ;
	    cout << "The code has not prepared for this kind of EOS!!!"
		 << endl;
	    abort() ;

	}  // End of strange stars case
	else {

	    int nzm1 = mp.get_mg()->get_nzone() - 1 ;

	    //--------------------------------
	    // Specific relativistic enthalpy		    ---> hhh
	    //--------------------------------

	    Tenseur hhh = exp(unsurc2 * ent) ;  // = 1 at the Newtonian limit
	    hhh.set_std_base() ;

	    //--------------------------------------------
	    // Computation of W^i = - A^2 h Gamma_n B^i/N
	    //--------------------------------------------

	    Tenseur www = - a_car * hhh * gam_euler * bsn ;

	    www.change_triad( mp.get_bvect_cart() ) ;
	                        // components on the mapping Cartesian basis

	    //-------------------------------------------------
	    // Constant value of W^i at the center of the star
	    //-------------------------------------------------

	    Tenseur v_orb(mp, 1, CON, mp.get_bvect_cart()) ;

	    v_orb.set_etat_qcq() ;
	    for (int i=0; i<3; i++) {
	        v_orb.set(i) = www(i)(0, 0, 0, 0) ;
	    }

	    v_orb.set_triad( *(www.get_triad()) ) ;
	    v_orb.set_std_base() ;

	    // Some auxiliary terms
	    // --------------------

	    const Coord& xx = mp.x ;
	    const Coord& yy = mp.y ;
	    const Coord& zz = mp.z ;

	    Tenseur r_bh(mp) ;
	    r_bh.set_etat_qcq() ;
	    r_bh.set() = pow( (xx+sepa)*(xx+sepa) + yy*yy + zz*zz, 0.5) ;
	    r_bh.set_std_base() ;

	    Tenseur msr(mp) ;
	    msr = ggrav * mass / r_bh ;
	    msr.set_std_base() ;

	    Tenseur xx_cov(mp, 1, COV, mp.get_bvect_cart()) ;
	    xx_cov.set_etat_qcq() ;
	    xx_cov.set(0) = xx + sepa ;
	    xx_cov.set(1) = yy ;
	    xx_cov.set(2) = zz ;
	    xx_cov.set_std_base() ;

	    // X_i in the spherical coordinate
	    const Coord& rr = mp.r ;
	    const Coord& st = mp.sint ;
	    const Coord& ct = mp.cost ;
	    const Coord& sp = mp.sinp ;
	    const Coord& cp = mp.cosp ;

	    Tenseur rr_spher_cov(mp, 1, COV, mp.get_bvect_spher()) ;
	    rr_spher_cov.set_etat_qcq() ;
	    rr_spher_cov.set(0) = rr + sepa*st*cp ;
	    rr_spher_cov.set(1) = sepa*ct*cp ;
	    rr_spher_cov.set(2) = - sepa*sp ;
	    rr_spher_cov.set_std_base() ;

	    Tenseur tmp_bh(mp) ;
	    tmp_bh = 0.5 * msr * msr / (1.-0.25*msr*msr) / r_bh / r_bh ;
	    tmp_bh.set_std_base() ;

	    //-------------------------------------------------
	    // Source and coefficients a,b for poisson_compact
	    // (idenpendent from psi0)
	    //-------------------------------------------------

	    Cmp dndh_log = eos.der_nbar_ent(ent(), nzet) ;

	    // In order to avoid any division by zero in the computation of
	    //  zeta_h the value of dndh_log is set to 1
	    //  in the external domains:
	    for (int l=nzet; l <= nzm1; l++) {
	        dndh_log.set(l) = 1 ;
	    }

	    double erreur ;

	    Tenseur zeta_h( ent() / dndh_log ) ;
	    zeta_h.set_std_base() ;

	    Tenseur tmp_zeta = 1 - unsurc2 * zeta_h ;
	    tmp_zeta.set_std_base() ;

	    Tenseur bb = tmp_zeta * ent.gradient_spher()
	      + unsurc2 * zeta_h * (beta_auto.gradient_spher()
				    + tmp_bh * rr_spher_cov) ;

	    Tenseur entmb = ent - beta_auto ;

	    Tenseur source = flat_scalar_prod(www - v_orb, ent.gradient())
	      + unsurc2*zeta_h*( flat_scalar_prod(v_orb, entmb.gradient())
				 - tmp_bh * flat_scalar_prod(v_orb, xx_cov)
				 + flat_scalar_prod(www, gam_euler.gradient())
				 / gam_euler ) ;

	    source.annule(nzet, nzm1) ;

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
	        psi0.set_etat_qcq() ;
		psi0.set() = 0 ;
	    }

	    source.set().va.ylm() ;

	    mp.poisson_compact(source(), zeta_h(), bb, par, psi0.set() ) ;

	    //-----------------------
	    // Check of the solution
	    //-----------------------

	    Tenseur bb_dpsi0 = flat_scalar_prod( bb, psi0.gradient_spher() ) ;

	    Cmp oper = zeta_h() * psi0().laplacien() + bb_dpsi0() ;

	    source.set().va.ylm_i() ;

	    erreur = diffrel(oper, source())(0) ;

	    cout << "Check of the resolution of the continuity equation : "
		 << endl ;
	    cout << "norme(source) : " << norme(source())(0)
		 << "    diff oper/source : " << erreur << endl ;

	    //--------------------------
	    // Computation of grad(psi)
	    //--------------------------

	    // The computation is done component by component because
	    // psi0.gradient() is a covariant vector, whereas v_orb is a
	    // contravariant one.

	    d_psi.set_etat_qcq() ;

	    for (int i=0; i<3; i++) {
	        d_psi.set(i) = (psi0.gradient())(i) + v_orb(i) ;
	    }

	    d_psi.set_triad( *(v_orb.get_triad()) ) ;

	    // C^1 continuation of d_psi outside the star
	    //  (to ensure a smooth enthalpy field accross the stellar surface)
	    // ----------------------------------------------------------------

	    d_psi.annule(nzet, nzm1) ;
	    for (int i=0; i<3; i++) {
	        d_psi.set(i) = raccord_c1(d_psi(i), nzet) ;
	    }

	    assert( d_psi.get_triad() == &(mp.get_bvect_cart()) ) ;

	    d_psi.change_triad(ref_triad) ;

	    return erreur ;

	}

    }

}
}
