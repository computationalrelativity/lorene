/*
 *  Methods of class Black_hole to compute global quantities
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

char blackhole_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2008/05/15 19:28:03  k_taniguchi
 * Change of some parameters.
 *
 * Revision 1.1  2007/06/22 01:19:51  k_taniguchi
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
#include "unites.h"
#include "utilitaires.h"

                    //-----------------------------------------//
                    //          Irreducible mass of BH         //
                    //-----------------------------------------//

double Black_hole::mass_irr() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_mass_irr == 0x0) {   // a new computation is required

        Scalar psi4(mp) ;
	psi4 = pow(confo, 4.) ;
	psi4.std_spectral_base() ;
	psi4.annule_domain(0) ;
	psi4.raccord(1) ;

	double radius_ah = mp.val_r(1,-1.,M_PI/2.,0.) ;

	Map_af& mp_aff= dynamic_cast<Map_af&>(mp) ;

	double a_ah = mp_aff.integrale_surface(psi4, radius_ah) ;
	double mirr = sqrt(a_ah/16./M_PI) / ggrav ;

	p_mass_irr = new double( mirr ) ;

    }

    return *p_mass_irr ;

}


                    //---------------------------//
                    //          ADM mass         //
                    //---------------------------//

double Black_hole::mass_adm() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_mass_adm == 0x0) {   // a new computation is required

        double madm ;
	double integ_s ;
	double integ_v ;

	double radius_ah = mp.val_r(1,-1.,M_PI/2.,0.) ;
	Map_af& mp_aff= dynamic_cast<Map_af&>(mp) ;

	Scalar source_surf(mp) ;
	source_surf.set_etat_qcq() ;
	Scalar source_volm(mp) ;
	source_volm.set_etat_qcq() ;

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

	Scalar lldconf = confo.dsdr() ;
	lldconf.std_spectral_base() ;

	if (kerrschild) {

	  cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	  abort() ;
	  /*
	  Scalar divshf(mp) ;
	  divshf = shift_rs(1).deriv(1) + shift_rs(2).deriv(2)
	    + shift_rs(3).deriv(3) ;
	  divshf.std_spectral_base() ;
	  divshf.dec_dzpuis(2) ;

	  Scalar llshift(mp) ;
	  llshift = ll(1)%shift_rs(1) + ll(2)%shift_rs(2)
	    + ll(3)%shift_rs(3) ;
	  llshift.std_spectral_base() ;

	  Scalar lldllsh = llshift.dsdr() ;
	  lldllsh.std_spectral_base() ;
	  lldllsh.dec_dzpuis(2) ;

	  double mass = ggrav * mass_bh ;

	  // Surface integral
	  // ----------------
	  source_surf = confo/rr - 2.*pow(confo,3.)*lapse_bh*mass/lapse/rr/rr
	    - pow(confo,3.)*(divshf - 3.*lldllsh
			     + 2.*mass*lapse_bh%lapse_bh%llshift/rr/rr
			     + 4.*mass*pow(lapse_bh,3.)*lapse_rs
			     *(1.+3.*mass/rr)/rr/rr)
	    /6./lapse/lapse_bh ;

	  source_surf.std_spectral_base() ;
	  source_surf.annule_domain(0) ;
	  source_surf.raccord(1) ;

	  integ_s = mp_aff.integrale_surface(source_surf, radius_ah) ;

	  // Volume integral
	  // ---------------
	  Scalar lldlldco = lldconf.dsdr() ;
	  lldlldco.std_spectral_base() ;

	  Scalar tmp1 = 2.*mass*mass*pow(lapse_bh,4.)*confo/pow(rr,4.) ;
	  tmp1.std_spectral_base() ;
	  tmp1.inc_dzpuis(4) ;

	  Scalar tmp2 = 2.*mass*mass*pow(lapse_bh,6.)
	    *(1.+3.*mass/rr)*(1.+3.*mass/rr)*pow(confo,5.)/3./pow(rr,4.) ;
	  tmp2.std_spectral_base() ;
	  tmp2.inc_dzpuis(4) ;

	  Scalar tmp3 = 4.*mass*lapse_bh*lapse_bh*lldlldco/rr ;
	  tmp3.std_spectral_base() ;
	  tmp3.inc_dzpuis(1) ;

	  Scalar tmp4 = 2.*mass*pow(lapse_bh,4.)*lldconf
	    *(3.+8.*mass/rr)/rr/rr ;
	  tmp4.std_spectral_base() ;
	  tmp4.inc_dzpuis(2) ;

	  source_volm = 0.25 * taij_quad / pow(confo,7.) - tmp1 - tmp2
	    - tmp3 - tmp4 ;

	  source_volm.annule_domain(0) ;
	  source_volm.std_spectral_base() ;

	  integ_v = source_volm.integrale() ;

	  // ADM mass
	  // --------
	  madm = mass_bh + integ_s / qpig + integ_v / qpig ;

	  // Another ADM mass
	  // ----------------
	  double mmm = mass_bh
	    - 2.*(mp_aff.integrale_surface_infini(lldconf))/qpig ;

	  cout << "Another ADM mass :   " << mmm / msol << endl ;
	  */
	}
	else {  // Isotropic coordinates with the maximal slicing

	  Scalar divshf(mp) ;
	  divshf = shift(1).deriv(1) + shift(2).deriv(2)
	    + shift(3).deriv(3) ;
	  divshf.std_spectral_base() ;
	  divshf.dec_dzpuis(2) ;

	  Scalar llshift(mp) ;
	  llshift = ll(1)%shift(1) + ll(2)%shift(2) + ll(3)%shift(3) ;
	  llshift.std_spectral_base() ;

	  Scalar lldllsh = llshift.dsdr() ;
	  lldllsh.std_spectral_base() ;
	  lldllsh.dec_dzpuis(2) ;

	  // Surface integral
	  // ----------------
	  source_surf = confo/rr
	    - pow(confo,4.) * (divshf - 3.*lldllsh) / lapconf / 6. ;

	  source_surf.std_spectral_base() ;
	  source_surf.annule_domain(0) ;
	  source_surf.raccord(1) ;

	  integ_s = mp_aff.integrale_surface(source_surf, radius_ah) ;

	  // Volume integral
	  // ---------------
	  source_volm = 0.25 * taij_quad / pow(confo,7.) ;

	  source_volm.std_spectral_base() ;
	  source_volm.annule_domain(0) ;

	  integ_v = source_volm.integrale() ;

	  // ADM mass
	  // --------
	  madm = integ_s / qpig + integ_v / qpig ;

	  // Another ADM mass
	  // ----------------
	  double mmm = - 2.*(mp_aff.integrale_surface_infini(lldconf))/qpig ;

	  cout << "Another ADM mass :   " << mmm / msol << endl ;

	}

	p_mass_adm = new double( madm ) ;

    }

    return *p_mass_adm ;

}

                    //-----------------------------//
                    //          Komar mass         //
                    //-----------------------------//

double Black_hole::mass_kom() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_mass_kom == 0x0) {   // a new computation is required

        double mkom ;
	double integ_s ;
	double integ_v ;

	double radius_ah = mp.val_r(1,-1.,M_PI/2.,0.) ;
	Map_af& mp_aff= dynamic_cast<Map_af&>(mp) ;

	Scalar source_surf(mp) ;
	source_surf.set_etat_qcq() ;
	Scalar source_volm(mp) ;
	source_volm.set_etat_qcq() ;

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

	Vector dlcnf(mp, CON, mp.get_bvect_cart()) ;
	dlcnf.set_etat_qcq() ;
	for (int i=1; i<=3; i++)
	  dlcnf.set(i) = confo.deriv(i) / confo ;

	dlcnf.std_spectral_base() ;

	if (kerrschild) {

	  cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	  abort() ;
	  /*
	  Scalar llshift(mp) ;
	  llshift = ll(1)%shift_rs(1) + ll(2)%shift_rs(2)
	    + ll(3)%shift_rs(3) ;
	  llshift.std_spectral_base() ;

	  Vector dlap(mp, CON, mp.get_bvect_cart()) ;
	  dlap.set_etat_qcq() ;

	  for (int i=1; i<=3; i++)
	    dlap.set(i) = lapse_rs.deriv(i) ;

	  dlap.std_spectral_base() ;

	  double mass = ggrav * mass_bh ;

	  // Surface integral
	  // ----------------
	  Scalar lldlap_bh(mp) ;
	  lldlap_bh = pow(lapse_bh,3.) * mass / rr / rr ;
	  lldlap_bh.std_spectral_base() ;
	  lldlap_bh.annule_domain(0) ;
	  lldlap_bh.inc_dzpuis(2) ;

	  Scalar lldlap_rs = lapse_rs.dsdr() ;
	  lldlap_rs.std_spectral_base() ;

	  source_surf = lldlap_rs + lldlap_bh ;
	  source_surf.std_spectral_base() ;
	  source_surf.dec_dzpuis(2) ;
	  source_surf.annule_domain(0) ;
	  source_surf.raccord(1) ;

	  integ_s = mp_aff.integrale_surface(source_surf, radius_ah) ;

	  // Volume integral
	  // ---------------
	  Scalar lldlldlap = lldlap_rs.dsdr() ;
	  lldlldlap.std_spectral_base() ;

	  Scalar lldlcnf = lconfo.dsdr() ;
	  lldlcnf.std_spectral_base() ;

	  Scalar tmp1(mp) ;
	  tmp1 = dlap(1)%dlcnf(1) + dlap(2)%dlcnf(2) + dlap(3)%dlcnf(3)
	    -2.*mass*lapse_bh%lapse_bh%lldlap_rs%lldlcnf/rr ;
	  tmp1.std_spectral_base() ;

	  Scalar tmp2(mp) ;
	  tmp2 = 4.*mass*mass*pow(lapse_bh,6.)*(1.+3.*mass/rr)*(1.+3.*mass/rr)
	    *lapse_rs*pow(confo,4.)/3./pow(rr,4.) ;
	  tmp2.std_spectral_base() ;
	  tmp2.inc_dzpuis(4) ;

	  Scalar tmp3(mp) ;
	  tmp3 = -2.*mass*pow(lapse_bh,5.)*llshift
	    *(2.+10.*mass/rr+9.*mass*mass/rr/rr)*pow(confo,4.)/pow(rr,3.) ;
	  tmp3.std_spectral_base() ;
	  tmp3.inc_dzpuis(4) ;

	  Scalar tmp4(mp) ;
	  tmp4 = 2.*mass*lapse_bh%lapse_bh%lldlldlap/rr ;
	  tmp4.std_spectral_base() ;
	  tmp4.inc_dzpuis(1) ;

	  Scalar tmp5(mp) ;
	  tmp5 = mass*pow(lapse_bh,4.)*lldlap_rs*(3.+8.*mass/rr)/rr/rr ;
	  tmp5.std_spectral_base() ;
	  tmp5.inc_dzpuis(2) ;

	  Scalar tmp6(mp) ;
	  tmp6 = -2.*pow(lapse_bh,5.)*mass*lldlcnf/rr/rr ;
	  tmp6.std_spectral_base() ;
	  tmp6.inc_dzpuis(2) ;

	  Scalar tmp_bh(mp) ;
	  tmp_bh = -pow(lapse_bh,7.)*mass*mass
	    *( 4.*(5.+24.*mass/rr+18.*mass*mass/rr/rr)*pow(confo,4.)/3.
	       + 1. - 6.*mass/rr) / pow(rr, 4.) ;
	  tmp_bh.std_spectral_base() ;
	  tmp_bh.inc_dzpuis(4) ;

	  source_volm = lapse % taij_quad / pow(confo,8.) - 2.*tmp1
	    + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp_bh ;

	  source_volm.annule_domain(0) ;
	  source_volm.std_spectral_base() ;

	  integ_v = source_volm.integrale() ;

	  // Komar mass
	  // ----------
	  mkom = integ_s / qpig + integ_v / qpig ;

	  // Another Komar mass
	  // ------------------
	  double mmm = (mp_aff.integrale_surface_infini(lldlap_rs+lldlap_bh))
	    / qpig ;

	  cout << "Another Komar mass : " << mmm / msol << endl ;
	  */
	}
	else {  // Isotropic coordinates with the maximal slicing

	  // Surface integral
	  // ----------------
	  Scalar lldlap = lapconf.dsdr() / confo
	    - lapconf * confo.dsdr() / confo / confo ;
	  lldlap.std_spectral_base() ;

	  source_surf = lldlap ;

	  source_surf.std_spectral_base() ;
	  source_surf.dec_dzpuis(2) ;
	  source_surf.annule_domain(0) ;
	  source_surf.raccord(1) ;

	  integ_s = mp_aff.integrale_surface(source_surf, radius_ah) ;

	  // Volume integral
	  // ---------------
	  Vector dlap(mp, CON, mp.get_bvect_cart()) ;
	  dlap.set_etat_qcq() ;

	  for (int i=1; i<=3; i++)
	    dlap.set(i) = lapconf.deriv(i) / confo
	      - lapconf * confo.deriv(i) / confo / confo ;

	  dlap.std_spectral_base() ;

	  source_volm = lapconf % taij_quad / pow(confo,9.)
	    - 2.*(dlap(1)%dlcnf(1) + dlap(2)%dlcnf(2) + dlap(3)%dlcnf(3)) ;

	  source_volm.std_spectral_base() ;
	  source_volm.annule_domain(0) ;

	  integ_v = source_volm.integrale() ;

	  // Komar mass
	  // ----------
	  mkom = integ_s / qpig + integ_v / qpig ;

	  // Another Komar mass
	  double mmm = (mp_aff.integrale_surface_infini(lldlap)) / qpig ;

	  cout << "Another Komar mass : " << mmm / msol << endl ;

	}

	p_mass_kom = new double( mkom ) ;

    }

    return *p_mass_kom ;

}

double Black_hole::rad_ah() const {

    if (p_rad_ah == 0x0) {   // a new computation is required

        Scalar rr(mp) ;
	rr = mp.r ;
	rr.std_spectral_base() ;

	double rad = rr.val_grid_point(1,0,0,0) ;

	p_rad_ah = new double( rad ) ;

    }

    return *p_rad_ah ;

}
