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
 * Revision 1.5  2014/10/13 08:52:46  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:02  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2008/07/02 20:45:58  k_taniguchi
 * Addition of routines to compute angular momentum.
 *
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
#include <cmath>

// Lorene headers
#include "blackhole.h"
#include "unites.h"
#include "utilitaires.h"

               //-----------------------------------------//
               //          Irreducible mass of BH         //
               //-----------------------------------------//

namespace Lorene {
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

	Map_af mp_aff(mp) ;

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
	Map_af mp_aff(mp) ;

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
	Map_af mp_aff(mp) ;

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


               //------------------------------------//
               //          Apparent horizon          //
               //------------------------------------//

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


          //-------------------------------------------//
          //     Quasi-local spin angular momentum     //
          //-------------------------------------------//

double Black_hole::spin_am_bh(bool bclapconf_nd, bool bclapconf_fs,
			      const Tbl& xi_i, const double& phi_i,
			      const double& theta_i, const int& nrk_phi,
			      const int& nrk_theta) const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_spin_am_bh == 0x0) {   // a new computation is required

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

	double radius_ah = mp.val_r(1,-1.,M_PI/2.,0.) ;

	if (kerrschild) {

	  cout << "Not yet prepared!!!" << endl ;
	  abort() ;

	}
	else {  // Isotropic coordinates

	  // Killing vector of the spherical components
	  Vector killing_spher(mp, COV, mp.get_bvect_spher()) ;
	  killing_spher.set_etat_qcq() ;
	  killing_spher = killing_vect_bh(xi_i, phi_i, theta_i,
					  nrk_phi, nrk_theta) ;
	  killing_spher.std_spectral_base() ;

	  killing_spher.set(2) = confo * confo * radius_ah * killing_spher(2) ;
	  killing_spher.set(3) = confo * confo * radius_ah * killing_spher(3) ;
	  // killing_spher(3) is already divided by sin(theta)
	  killing_spher.std_spectral_base() ;

	  // Killing vector of the Cartesian components
	  Vector killing(mp, COV, mp.get_bvect_cart()) ;
	  killing.set_etat_qcq() ;

	  killing.set(1) = (killing_spher(2) * ct * cp - killing_spher(3) * sp)
	    / radius_ah  ;
	  killing.set(2) = (killing_spher(2) * ct * sp + killing_spher(3) * cp)
	    / radius_ah ;
	  killing.set(3) = - killing_spher(2) * st / radius_ah ;
	  killing.std_spectral_base() ;

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  // Source terms in the surface integral
	  Scalar source_1(mp) ;
	  source_1 = (ll(1) * (taij(1,1) * killing(1)
			       + taij(1,2) * killing(2)
			       + taij(1,3) * killing(3))
		      + ll(2) * (taij(2,1) * killing(1)
				 + taij(2,2) * killing(2)
				 + taij(2,3) * killing(3))
		      + ll(3) * (taij(3,1) * killing(1)
				 + taij(3,2) * killing(2)
				 + taij(3,3) * killing(3)))
	    / pow(confo, 4.) ;
	  source_1.std_spectral_base() ;
	  source_1.dec_dzpuis(2) ;  // dzpuis : 2 -> 0

	  Scalar source_surf(mp) ;
	  source_surf = source_1 ;
	  source_surf.std_spectral_base() ;
	  source_surf.annule_domain(0) ;
	  source_surf.raccord(1) ;

	  Map_af mp_aff(mp) ;

	  double spin = mp_aff.integrale_surface(source_surf, radius_ah) ;
	  double spin_angmom = 0.5 * spin / qpig ;

	  p_spin_am_bh = new double( spin_angmom ) ;

	}

    }

    return *p_spin_am_bh ;

}

          //------------------------------------------//
          //          Total angular momentum          //
          //------------------------------------------//

const Tbl& Black_hole::angu_mom_bh() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_angu_mom_bh == 0x0) {   // a new computation is required

        p_angu_mom_bh = new Tbl(3) ;
	p_angu_mom_bh->annule_hard() ;  // fills the double array with zeros

	double integ_bh_s_x ;
	double integ_bh_s_y ;
	double integ_bh_s_z ;

	double radius_ah = mp.val_r(1,-1.,M_PI/2.,0.) ;
	Map_af mp_aff(mp) ;

	Scalar xx(mp) ;
	xx = mp.x ;
	xx.std_spectral_base() ;
	Scalar yy(mp) ;
	yy = mp.y ;
	yy.std_spectral_base() ;
	Scalar zz(mp) ;
	zz = mp.z ;
	zz.std_spectral_base() ;

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

	Vector jbh_x(mp, CON, mp.get_bvect_cart()) ;
	jbh_x.set_etat_qcq() ;
	for (int i=1; i<=3; i++)
	  jbh_x.set(i) = yy * taij(3,i) - zz * taij(2,i) ;

	jbh_x.std_spectral_base() ;

	Vector jbh_y(mp, CON, mp.get_bvect_cart()) ;
	jbh_y.set_etat_qcq() ;
	for (int i=1; i<=3; i++)
	  jbh_y.set(i) = zz * taij(1,i) - xx * taij(3,i) ;

	jbh_y.std_spectral_base() ;

	Vector jbh_z(mp, CON, mp.get_bvect_cart()) ;
	jbh_z.set_etat_qcq() ;
	for (int i=1; i<=3; i++)
	  jbh_z.set(i) = xx * taij(2,i) - yy * taij(1,i) ;

	jbh_z.std_spectral_base() ;

	/*
	if (kerrschild) {

	  cout << "Not yet prepared!!!" << endl ;
	  abort() ;

	}
	else {  // Isotropic coordinates

	  // Sets C/M^2 for each case of the lapse boundary condition
	  // --------------------------------------------------------
	  double cc ;

	  if (bclapconf_nd) {  // Neumann boundary condition
	    if (bclapconf_fs) {  // First condition
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
	    if (bclapconf_fs) {  // First condition
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
	    }
	  }

	  Scalar r_are(mp) ;
	  r_are = r_coord(bclapconf_nd, bclapconf_fs) ;
	  r_are.std_spectral_base() ;

	}
	*/

	//-------------//
	// X component //
	//-------------//

	//-------------------------------------
	//     Integration over the BH map
	//-------------------------------------

	// Surface integral <- dzpuis should be 0
	// --------------------------------------
	Scalar sou_bh_sx(mp) ;
	sou_bh_sx = jbh_x(1)%ll(1) + jbh_x(2)%ll(2) + jbh_x(3)%ll(3) ;
	sou_bh_sx.std_spectral_base() ;
	sou_bh_sx.dec_dzpuis(2) ;  // dzpuis : 2 -> 0
	sou_bh_sx.annule_domain(0) ;
	sou_bh_sx.raccord(1) ;

	integ_bh_s_x = mp_aff.integrale_surface(sou_bh_sx, radius_ah) ;

	p_angu_mom_bh->set(0) = 0.5 * integ_bh_s_x / qpig ;

	//-------------//
	// Y component //
	//-------------//

	//-------------------------------------
	//     Integration over the BH map
	//-------------------------------------

	// Surface integral <- dzpuis should be 0
	// --------------------------------------
	Scalar sou_bh_sy(mp) ;
	sou_bh_sy = jbh_y(1)%ll(1) + jbh_y(2)%ll(2) + jbh_y(3)%ll(3) ;
	sou_bh_sy.std_spectral_base() ;
	sou_bh_sy.dec_dzpuis(2) ;  // dzpuis : 2 -> 0
	sou_bh_sy.annule_domain(0) ;
	sou_bh_sy.raccord(1) ;

	integ_bh_s_y = mp_aff.integrale_surface(sou_bh_sy, radius_ah) ;

	p_angu_mom_bh->set(1) = 0.5 * integ_bh_s_y / qpig ;

	//-------------//
	// Z component //
	//-------------//

	//-------------------------------------
	//     Integration over the BH map
	//-------------------------------------

	// Surface integral <- dzpuis should be 0
	// --------------------------------------
	Scalar sou_bh_sz(mp) ;
	sou_bh_sz = jbh_z(1)%ll(1) + jbh_z(2)%ll(2) + jbh_z(3)%ll(3) ;
	sou_bh_sz.std_spectral_base() ;
	sou_bh_sz.dec_dzpuis(2) ;  // dzpuis : 2 -> 0
	sou_bh_sz.annule_domain(0) ;
	sou_bh_sz.raccord(1) ;

	integ_bh_s_z = mp_aff.integrale_surface(sou_bh_sz, radius_ah) ;

	p_angu_mom_bh->set(2) = 0.5 * integ_bh_s_z / qpig ;

    }

    return *p_angu_mom_bh ;

}
}
