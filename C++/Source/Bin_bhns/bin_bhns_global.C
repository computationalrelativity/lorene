/*
 *  Methods of class Bin_bhns to compute global quantities
 *
 *    (see file bin_bhns.h for documentation).
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

char bin_bhns_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:13:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/05/15 18:59:27  k_taniguchi
 * Introduction of new global quantities.
 *
 * Revision 1.1  2007/06/22 01:09:31  k_taniguchi
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
#include "bin_bhns.h"
#include "blackhole.h"
#include "unites.h"
#include "utilitaires.h"
#include "nbr_spx.h"

               //----------------------------//
               //          ADM mass          //
               //----------------------------//

double Bin_bhns::mass_adm_bhns_surf() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_mass_adm_bhns_surf == 0x0) {   // a new computation is required

        double madm ;

	const Map& mp_bh = hole.get_mp() ;
        const Map& mp_ns = star.get_mp() ;

	Map_af mp_aff(mp_bh) ;
	Map_af mp_ns_aff(mp_ns) ;

	Scalar rr(mp_bh) ;
	rr = mp_bh.r ;
	rr.std_spectral_base() ;

	double mass = ggrav * hole.get_mass_bh() ;

	if (hole.is_kerrschild()) {

	  cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	  abort() ;

	}
	else { // Isotropic coordinates with the maximal slicing

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Sets C/M^2 for each case of the lapse boundary condition
	  // --------------------------------------------------------
	  double cc ;

	  if (hole.has_bc_lapconf_nd()) {  // Neumann boundary condition
	      if (hole.has_bc_lapconf_fs()) {  // First condition
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
	      if (hole.has_bc_lapconf_fs()) {  // First condition
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
		  //	          cc = 2. * sqrt(2.) ;
	      }
	  }

	  Scalar r_are(mp_bh) ;
	  r_are = hole.r_coord(hole.has_bc_lapconf_nd(),
			       hole.has_bc_lapconf_fs()) ;
	  r_are.std_spectral_base() ;

	  // ADM mass by surface integral at infinity : dzpuis should be 2
	  // ----------------------------------------
	  const Scalar& confo_bh_auto_rs = hole.get_confo_auto_rs() ;

	  Scalar lldconf_iso = confo_bh_auto_rs.dsdr() ;  // dzpuis = 2
	  lldconf_iso.std_spectral_base() ;

	  Scalar anoth(mp_bh) ;
	  anoth = 0.5 * sqrt(r_are)
	    * (sqrt(1. -2.*mass/r_are/rr + cc*cc*pow(mass/r_are/rr,4.))
	       - 1.) / rr ;
	  anoth.std_spectral_base() ;
	  anoth.annule_domain(0) ;
	  anoth.raccord(1) ;
	  anoth.inc_dzpuis(2) ;

	  const Scalar& confo_ns_auto = star.get_confo_auto() ;

	  Scalar lldconf_ns = confo_ns_auto.dsdr() ;  // dzpuis = 2
	  lldconf_ns.std_spectral_base() ;

	  madm =
	    - 2.*(mp_aff.integrale_surface_infini(lldconf_iso+anoth))/qpig
	    - 2.*(mp_ns_aff.integrale_surface_infini(lldconf_ns))/qpig ;

	  cout << "ADM mass (surface) :   " << madm / msol << " [Mo]"
	       << endl ;

	}

	p_mass_adm_bhns_surf = new double( madm ) ;

    }

    return *p_mass_adm_bhns_surf ;

}


               //----------------------------//
               //          ADM mass          //
               //----------------------------//

double Bin_bhns::mass_adm_bhns_vol() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_mass_adm_bhns_vol == 0x0) {   // a new computation is required

        double madm ;
	double integ_bh_s ;
	double integ_bh_v ;
	double integ_ns_v ;

	const Map& mp_bh = hole.get_mp() ;
        const Map& mp_ns = star.get_mp() ;

	double radius_ah = mp_bh.val_r(1,-1.,M_PI/2.,0.) ;
	Map_af mp_aff(mp_bh) ;

	Map_af mp_ns_aff(mp_ns) ;

	Scalar source_bh_surf(mp_bh) ;
	source_bh_surf.set_etat_qcq() ;

	Scalar source_bh_volm(mp_bh) ;
	source_bh_volm.set_etat_qcq() ;

	Scalar source_ns_volm(mp_ns) ;
	source_ns_volm.set_etat_qcq() ;

	Scalar rr(mp_bh) ;
	rr = mp_bh.r ;
	rr.std_spectral_base() ;
	Scalar st(mp_bh) ;
	st = mp_bh.sint ;
	st.std_spectral_base() ;
	Scalar ct(mp_bh) ;
	ct = mp_bh.cost ;
	ct.std_spectral_base() ;
	Scalar sp(mp_bh) ;
	sp = mp_bh.sinp ;
	sp.std_spectral_base() ;
	Scalar cp(mp_bh) ;
	cp = mp_bh.cosp ;
	cp.std_spectral_base() ;

	Vector ll(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	ll.set_etat_qcq() ;
	ll.set(1) = st % cp ;
	ll.set(2) = st % sp ;
	ll.set(3) = ct ;
	ll.std_spectral_base() ;

	const Vector& shift_auto_rs = hole.get_shift_auto_rs() ;
	const Vector& shift_comp = hole.get_shift_comp() ;
	const Tensor& dshift_comp = hole.get_d_shift_comp() ;

	Scalar divshift(mp_bh) ;  // dzpuis = 2
	divshift = shift_auto_rs(1).deriv(1) + shift_auto_rs(2).deriv(2)
	  + shift_auto_rs(3).deriv(3) + dshift_comp(1,1)
	  + dshift_comp(2,2) + dshift_comp(3,3) ;
	divshift.std_spectral_base() ;

	Scalar llshift(mp_bh) ;   // dzpuis = 0
	llshift = ll(1) % (shift_auto_rs(1) + shift_comp(1))
	  + ll(2) % (shift_auto_rs(2) + shift_comp(2))
	  + ll(3) % (shift_auto_rs(3) + shift_comp(3)) ;
	llshift.std_spectral_base() ;

	Scalar llshift_auto(mp_bh) ;   // dzpuis = 0
	llshift_auto = ll(1)%shift_auto_rs(1) + ll(2)%shift_auto_rs(2)
	  + ll(3)%shift_auto_rs(3) ;
	llshift_auto.std_spectral_base() ;

	Scalar lldllsh = llshift_auto.dsdr()
	  + ll(1) * (ll(1) % dshift_comp(1,1) + ll(2) % dshift_comp(1,2)
		     + ll(3) % dshift_comp(1,3))
	  + ll(2) * (ll(1) % dshift_comp(2,1) + ll(2) % dshift_comp(2,2)
		     + ll(3) % dshift_comp(2,3))
	  + ll(3) * (ll(1) % dshift_comp(3,1) + ll(2) % dshift_comp(3,2)
		     + ll(3) % dshift_comp(3,3)) ;  // dzpuis = 2
	lldllsh.std_spectral_base() ;

	const Scalar& lapconf_bh = hole.get_lapconf_tot() ;
	const Scalar& lapconf_bh_auto_rs = hole.get_lapconf_auto_rs() ;
	const Scalar& lapconf_bh_comp = hole.get_lapconf_comp() ;
	const Scalar& confo_bh = hole.get_confo_tot() ;
	const Scalar& confo_bh_auto = hole.get_confo_auto() ;
	const Scalar& confo_bh_comp = hole.get_confo_comp() ;
	const Vector& dconfo_bh_comp = hole.get_d_confo_comp() ;
	const Scalar& taij_quad_tot_rs = hole.get_taij_quad_tot_rs() ;
	const Scalar& taij_quad_tot_rot = hole.get_taij_quad_tot_rot() ;

	const Scalar& taij_quad_auto_bh = hole.get_taij_quad_auto() ;
	const Scalar& taij_quad_comp_bh = hole.get_taij_quad_comp() ;

	const Scalar& confo_ns = star.get_confo_tot() ;
	const Scalar& confo_ns_auto = star.get_confo_auto() ;
	const Scalar& ener_euler = star.get_ener_euler() ;

	const Scalar& taij_quad_auto_ns = star.get_taij_quad_auto() ;

	Scalar lldconf = confo_bh_auto.dsdr() + ll(1)%dconfo_bh_comp(1)
	  + ll(2)%dconfo_bh_comp(2) + ll(3)%dconfo_bh_comp(3) ;  // dzpuis = 2
	lldconf.std_spectral_base() ;

	double mass = ggrav * hole.get_mass_bh() ;

	if (hole.is_kerrschild()) {

	  cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	  abort() ;

	  /*
	  Scalar lap_bh(mp_bh) ;
	  lap_bh = 1./sqrt(1.+2.*mass/rr) ;
	  lap_bh.std_spectral_base() ;

	  Scalar lap_bh2(mp_bh) ;
	  lap_bh2 = 1./(1.+2.*mass/rr) ;
	  lap_bh2.std_spectral_base() ;

	  Scalar omelld(mp_bh) ;
	  omelld = omega * (ll(2) * (mp_bh.get_ori_x() - x_rot)
			    - ll(1) * (mp_bh.get_ori_y() - y_rot)) ;
	  omelld.std_spectral_base() ;

	  Scalar lldlldconf(mp_bh) ;  // dzpuis = 3
	  lldlldconf = lldconf.dsdr() ;
	  lldlldconf.std_spectral_base() ;

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  Scalar divshift_zero(divshift) ;
	  divshift_zero.dec_dzpuis(2) ;

	  Scalar lldllsh_zero(lldllsh) ;
	  lldllsh_zero.dec_dzpuis(2) ;

	  source_bh_surf = confo_bh
	    * (1. - 2.*mass*lap_bh*confo_bh*confo_bh/lapse_bh/rr) / rr
	    - pow(confo_bh, 3.)
	    * ( divshift_zero - 3.*lldllsh_zero
		+ 2. * lap_bh2 * mass * (llshift + omelld) / rr / rr
		+ 4.*mass*lap_bh2*lap_bh*(1.+3.*mass/rr)
		*(lapse_bh_auto_rs + lapse_bh_comp)/rr/rr
		) / 6. / lap_bh / lapse_bh ;

	  source_bh_surf.std_spectral_base() ;
	  source_bh_surf.annule_domain(0) ;
	  source_bh_surf.raccord(1) ;

	  integ_bh_s = mp_aff.integrale_surface(source_bh_surf, radius_ah) ;

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  Scalar sou_bh1(mp_bh) ;
	  sou_bh1 = 2.*lap_bh2*mass*lldlldconf/rr ;
	  sou_bh1.std_spectral_base() ;
	  sou_bh1.inc_dzpuis(1) ;

	  Scalar sou_bh2(mp_bh) ;
	  sou_bh2 = lap_bh2*lap_bh2*mass*(3.+8.*mass/rr)*lldconf/rr/rr ;
	  sou_bh2.std_spectral_base() ;
	  sou_bh2.inc_dzpuis(2) ;

	  Scalar sou_bh3(mp_bh) ;
	  sou_bh3 = pow(lap_bh2,3.)*mass*mass*confo_bh
	    * ( (1.-lap_bh2/lapse_bh/lapse_bh)
		*(4.+12.*mass/rr+9.*mass*mass/rr/rr)*pow(confo_bh,4.)
		+3.*(1.+2.*mass/rr)*(1.-pow(confo_bh,4.)) )
	    /3./pow(rr,4.) ;
	  sou_bh3.std_spectral_base() ;
	  sou_bh3.inc_dzpuis(4) ;

	  source_bh_volm = 0.25 * (taij_quad_tot_rs + taij_quad_tot_rot)
	    / pow(confo_bh,7.)
	    - 2. * (sou_bh1 + sou_bh2 + sou_bh3) ;

	  source_bh_volm.std_spectral_base() ;
	  source_bh_volm.annule_domain(0) ;

	  integ_bh_v = source_bh_volm.integrale() ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  source_ns_volm = pow(confo_ns, 5.) * ener_euler ;
	  source_ns_volm.std_spectral_base() ;
	  source_ns_volm.inc_dzpuis(4) ;

	  integ_ns_v = source_ns_volm.integrale() ;

	  cout << "integ_bh_s : " << integ_bh_s/ qpig / msol
	       << "  integ_bh_v : "
	       << integ_bh_v/ qpig / msol
	       << "  integ_ns_v : " << integ_ns_v/ msol << endl ;

	  //------------------
	  //     ADM mass
	  //------------------
	  madm = hole.get_mass_bh()
	    + (integ_bh_s + integ_bh_v) / qpig + integ_ns_v ;

	  cout << "ADM mass :           " << madm / msol << " [Mo]"
	       << endl ;

	  // ADM mass by surface integral at infinity : dzpuis should be 2
	  // ----------------------------------------
	  double mmm = hole.get_mass_bh()
	    - 2.*(mp_aff.integrale_surface_infini(lldconf))/qpig ;

	  cout << "Another ADM mass :   " << mmm / msol << " [Mo]"
	       << endl ;
	  */
	}
	else { // Isotropic coordinates with the maximal slicing

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Sets C/M^2 for each case of the lapse boundary condition
	  // --------------------------------------------------------
	  double cc ;

	  if (hole.has_bc_lapconf_nd()) {  // Neumann boundary condition
	      if (hole.has_bc_lapconf_fs()) {  // First condition
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
	      if (hole.has_bc_lapconf_fs()) {  // First condition
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
		  //	          cc = 2. * sqrt(2.) ;
	      }
	  }

	  Scalar r_are(mp_bh) ;
	  r_are = hole.r_coord(hole.has_bc_lapconf_nd(),
			       hole.has_bc_lapconf_fs()) ;
	  r_are.std_spectral_base() ;

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  Scalar divshift_zero(divshift) ;
	  divshift_zero.dec_dzpuis(2) ;

	  Scalar lldllsh_zero(lldllsh) ;
	  lldllsh_zero.dec_dzpuis(2) ;

	  source_bh_surf = confo_bh / rr
	    - pow(confo_bh, 4.) * (divshift_zero - 3.*lldllsh_zero)
	    /6./lapconf_bh
	    - pow(confo_bh, 4.) * mass * mass * cc
	    * sqrt(1. -2.*mass/r_are/rr + cc*cc*pow(mass/r_are/rr,4.))
	    / lapconf_bh / pow(r_are*rr,3.) ;

	  source_bh_surf.std_spectral_base() ;
	  source_bh_surf.annule_domain(0) ;
	  source_bh_surf.raccord(1) ;

	  integ_bh_s = mp_aff.integrale_surface(source_bh_surf, radius_ah) ;

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  Scalar sou_bh1(mp_bh) ;
	  sou_bh1 = 1.5 * pow(confo_bh,7.) * pow(mass*mass*cc,2.)
	    * (1. - 2.*mass/r_are/rr + cc*cc*pow(mass/r_are/rr,4.))
	    / lapconf_bh / lapconf_bh / pow(r_are*rr,6.) ;
	  sou_bh1.std_spectral_base() ;
	  sou_bh1.inc_dzpuis(4) ;

	  source_bh_volm = 0.25 * taij_quad_auto_bh / pow(confo_bh,7.)
	    + 0.25 * taij_quad_comp_bh
	    * (pow(confo_bh/(confo_bh_comp+0.5),7.)
	       *pow((lapconf_bh_comp+0.5)/lapconf_bh,2.)
	       - 1.) / pow(confo_bh_comp+0.5,7.)
	    + sou_bh1 ;

	  source_bh_volm.std_spectral_base() ;
	  source_bh_volm.annule_domain(0) ;

	  integ_bh_v = source_bh_volm.integrale() ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  Scalar sou_ns1(mp_ns) ;
	  sou_ns1 = pow(confo_ns, 5.) * ener_euler ;
	  sou_ns1.std_spectral_base() ;
	  sou_ns1.inc_dzpuis(4) ;

	  source_ns_volm = 0.25 * taij_quad_auto_ns
	    / pow(confo_ns_auto+0.5, 7.) / qpig + sou_ns1 ;

	  source_ns_volm.std_spectral_base() ;

	  integ_ns_v = source_ns_volm.integrale() ;

	  cout << "integ_bh_s : " << integ_bh_s/ qpig / msol
	       << "  integ_bh_v : "
	       << integ_bh_v/ qpig / msol
	       << "  integ_ns_v : " << integ_ns_v/ msol << endl ;

	  //------------------
	  //     ADM mass
	  //------------------
	  madm = (integ_bh_s + integ_bh_v) / qpig + integ_ns_v ;

	  cout << "ADM mass (volume) :    " << madm / msol << " [Mo]"
	       << endl ;

	}

	p_mass_adm_bhns_vol = new double( madm ) ;

    }

    return *p_mass_adm_bhns_vol ;

}



               //------------------------------//
               //          Komar mass          //
               //------------------------------//

double Bin_bhns::mass_kom_bhns_surf() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_mass_kom_bhns_surf == 0x0) {   // a new computation is required

        double mkom ;

	const Map& mp_bh = hole.get_mp() ;
        const Map& mp_ns = star.get_mp() ;

	Map_af mp_aff(mp_bh) ;
	Map_af mp_ns_aff(mp_ns) ;

	Scalar rr(mp_bh) ;
	rr = mp_bh.r ;
	rr.std_spectral_base() ;

	double mass = ggrav * hole.get_mass_bh() ;

	if (hole.is_kerrschild()) {

	  cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	  abort() ;

	}
	else { // Isotropic coordinates with the maximal slicing

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Sets C/M^2 for each case of the lapse boundary condition
	  // --------------------------------------------------------
	  double cc ;

	  if (hole.has_bc_lapconf_nd()) {  // Neumann boundary condition
	      if (hole.has_bc_lapconf_fs()) {  // First condition
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
	      if (hole.has_bc_lapconf_fs()) {  // First condition
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
		  //	          cc = 2. * sqrt(2.) ;
	      }
	  }

	  Scalar r_are(mp_bh) ;
	  r_are = hole.r_coord(hole.has_bc_lapconf_nd(),
			       hole.has_bc_lapconf_fs()) ;
	  r_are.std_spectral_base() ;

	  // Komar mass by surface integral at infinity : dzpuis should be 2
	  // ------------------------------------------
	  const Scalar& lapconf_bh_auto_rs = hole.get_lapconf_auto_rs() ;
	  const Scalar& confo_bh_auto_rs = hole.get_confo_auto_rs() ;

	  Scalar lldlap_bh(mp_bh) ;
	  lldlap_bh = lapconf_bh_auto_rs.dsdr()
	    - confo_bh_auto_rs.dsdr() ; // dzpuis = 2
	  lldlap_bh.std_spectral_base() ;

	  Scalar anoth(mp_bh) ;
	  anoth = sqrt(r_are) * (1. - 1.5 *cc*cc*pow(mass/r_are/rr,4.)
				 - sqrt(1. - 2.*mass/r_are/rr
					+ cc*cc*pow(mass/r_are/rr,4.))) / rr ;
	  anoth.std_spectral_base() ;
	  anoth.annule_domain(0) ;
	  anoth.raccord(1) ;
	  anoth.inc_dzpuis(2) ;

	  const Scalar& lapconf_ns_auto = star.get_lapconf_auto() ;
	  const Scalar& confo_ns_auto = star.get_confo_auto() ;

	  Scalar lldlap_ns(mp_ns) ;
	  lldlap_ns = lapconf_ns_auto.dsdr() - confo_ns_auto.dsdr() ;
	  lldlap_ns.std_spectral_base() ; // dzpuis = 2

	  mkom =
	    (mp_aff.integrale_surface_infini(lldlap_bh+anoth))/qpig
	    + (mp_ns_aff.integrale_surface_infini(lldlap_ns))/qpig ;

	  cout << "Komar mass (surface) : " << mkom / msol << " [Mo]"
	       << endl ;

	}

	p_mass_kom_bhns_surf = new double( mkom ) ;

    }

    return *p_mass_kom_bhns_surf ;

}



               //------------------------------//
               //          Komar mass          //
               //------------------------------//

double Bin_bhns::mass_kom_bhns_vol() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_mass_kom_bhns_vol == 0x0) {   // a new computation is required

        double mkom ;
	double integ_bh_s ;
	double integ_bh_v ;
	double integ_ns_v ;

	const Map& mp_bh = hole.get_mp() ;
        const Map& mp_ns = star.get_mp() ;

	double radius_ah = mp_bh.val_r(1,-1.,M_PI/2.,0.) ;
	Map_af mp_aff(mp_bh) ;

	Map_af mp_ns_aff(mp_ns) ;

	Scalar source_bh_surf(mp_bh) ;
	source_bh_surf.set_etat_qcq() ;

	Scalar source_bh_volm(mp_bh) ;
	source_bh_volm.set_etat_qcq() ;

	Scalar source_ns_volm(mp_ns) ;
	source_ns_volm.set_etat_qcq() ;

	Scalar rr(mp_bh) ;
	rr = mp_bh.r ;
	rr.std_spectral_base() ;
	Scalar st(mp_bh) ;
	st = mp_bh.sint ;
	st.std_spectral_base() ;
	Scalar ct(mp_bh) ;
	ct = mp_bh.cost ;
	ct.std_spectral_base() ;
	Scalar sp(mp_bh) ;
	sp = mp_bh.sinp ;
	sp.std_spectral_base() ;
	Scalar cp(mp_bh) ;
	cp = mp_bh.cosp ;
	cp.std_spectral_base() ;

	Vector ll(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	ll.set_etat_qcq() ;
	ll.set(1) = st % cp ;
	ll.set(2) = st % sp ;
	ll.set(3) = ct ;
	ll.std_spectral_base() ;

	const Scalar& lapconf_bh = hole.get_lapconf_tot() ;
	const Scalar& lapconf_bh_auto_rs = hole.get_lapconf_auto_rs() ;
	const Scalar& lapconf_bh_comp = hole.get_lapconf_comp() ;
	const Vector& dlapconf_bh_comp = hole.get_d_lapconf_comp() ;
	const Scalar& confo_bh = hole.get_confo_tot() ;
	const Scalar& confo_bh_auto_rs = hole.get_confo_auto_rs() ;
	const Scalar& confo_bh_comp = hole.get_confo_comp() ;
	const Vector& dconfo_bh_comp = hole.get_d_confo_comp() ;
	const Scalar& taij_quad_tot_rs = hole.get_taij_quad_tot_rs() ;
	const Scalar& taij_quad_tot_rot = hole.get_taij_quad_tot_rot() ;

	const Scalar& taij_quad_auto_bh = hole.get_taij_quad_auto() ;
	const Scalar& taij_quad_comp_bh = hole.get_taij_quad_comp() ;

	const Vector& shift_auto_rs = hole.get_shift_auto_rs() ;
	const Vector& shift_comp = hole.get_shift_comp() ;
	const Tensor& dshift_comp = hole.get_d_shift_comp() ;

	Scalar divshift(mp_bh) ;  // dzpuis = 2
	divshift = shift_auto_rs(1).deriv(1) + shift_auto_rs(2).deriv(2)
	  + shift_auto_rs(3).deriv(3) + dshift_comp(1,1)
	  + dshift_comp(2,2) + dshift_comp(3,3) ;
	divshift.std_spectral_base() ;

	Scalar llshift_auto(mp_bh) ;   // dzpuis = 0
        llshift_auto = ll(1)%shift_auto_rs(1) + ll(2)%shift_auto_rs(2)
          + ll(3)%shift_auto_rs(3) ;
        llshift_auto.std_spectral_base() ;

	Scalar lldllsh = llshift_auto.dsdr()
          + ll(1) * (ll(1) % dshift_comp(1,1) + ll(2) % dshift_comp(1,2)
                     + ll(3) % dshift_comp(1,3))
          + ll(2) * (ll(1) % dshift_comp(2,1) + ll(2) % dshift_comp(2,2)
                     + ll(3) % dshift_comp(2,3))
          + ll(3) * (ll(1) % dshift_comp(3,1) + ll(2) % dshift_comp(3,2)
                     + ll(3) % dshift_comp(3,3)) ;  // dzpuis = 2
        lldllsh.std_spectral_base() ;

	const Scalar& lapconf_ns = star.get_lapconf_tot() ;
	const Scalar& ener_euler = star.get_ener_euler() ;
	const Scalar& s_euler = star.get_s_euler() ;

	const Scalar& confo_ns = star.get_confo_tot() ;
	const Scalar& lapconf_ns_auto = star.get_lapconf_auto() ;
	const Scalar& confo_ns_auto = star.get_confo_auto() ;
	const Vector& dconfo_ns_comp = star.get_d_confo_comp() ;
	const Scalar& taij_quad_auto_ns = star.get_taij_quad_auto() ;

	const Vector& dlapconf_bh_auto_rs = hole.get_d_lapconf_auto_rs() ;
	/*
	Vector dlc_bh_auto_rs(mp_bh, COV, mp_bh.get_bvect_cart()) ;
	dlc_bh_auto_rs.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	  dlc_bh_auto_rs.set(i) = lapconf_bh_auto_rs.deriv(i) ;
	}
	dlc_bh_auto_rs.std_spectral_base() ;
	*/

	const Vector& dconfo_bh_auto_rs = hole.get_d_confo_auto_rs() ;
	/*
	Vector dc_bh_auto_rs(mp_bh, COV, mp_bh.get_bvect_cart()) ;
	dc_bh_auto_rs.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	  dc_bh_auto_rs.set(i) = confo_bh_auto_rs.deriv(i) ;
	}
	dc_bh_auto_rs.std_spectral_base() ;
	*/

	const Vector& dlapconf_ns_auto = star.get_d_lapconf_auto() ;
	/*
	Vector dlc_ns_auto(mp_ns, COV, mp_ns.get_bvect_cart()) ;
	dlc_ns_auto.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	  dlc_ns_auto.set(i) = lapconf_ns_auto.deriv(i) ;
	}
	dlc_ns_auto.std_spectral_base() ;
	*/

	const Vector& dconfo_ns_auto = star.get_d_confo_auto() ;
	/*
	Vector dc_ns_auto(mp_ns, COV, mp_ns.get_bvect_cart()) ;
	dc_ns_auto.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	  dc_ns_auto.set(i) = confo_ns_auto.deriv(i) ;
	}
	dc_ns_auto.std_spectral_base() ;
	*/
	double mass = ggrav * hole.get_mass_bh() ;

	if (hole.is_kerrschild()) {

	  cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	  abort() ;
	  /*
	  const Vector& shift_auto_rs = hole.get_shift_auto_rs() ;
	  const Vector& shift_comp = hole.get_shift_comp() ;

	  Scalar llshift(mp_bh) ;  // dzpuis = 0
	  llshift = ll(1) % (shift_auto_rs(1) + shift_comp(1))
	    + ll(2) % (shift_auto_rs(2) + shift_comp(2))
	    + ll(3) % (shift_auto_rs(3) + shift_comp(3)) ;
	  llshift.std_spectral_base() ;

	  Scalar lldlldlap(mp_bh) ;  // dzpuis = 3
	  lldlldlap = lldlap.dsdr() ;
	  lldlldlap.std_spectral_base() ;

	  Scalar lap_bh2(mp_bh) ;
	  lap_bh2 = 1./(1.+2.*mass/rr) ;
	  lap_bh2.std_spectral_base() ;

	  Scalar omelld(mp_bh) ;
	  omelld = omega * (ll(2) * (mp_bh.get_ori_x() - x_rot)
			    - ll(1) * (mp_bh.get_ori_y() - y_rot)) ;
	  omelld.std_spectral_base() ;

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  source_bh_surf = lldlap ;

	  source_bh_surf.std_spectral_base() ;
	  source_bh_surf.annule_domain(0) ;
	  source_bh_surf.raccord(1) ;
	  source_bh_surf.dec_dzpuis(2) ;

	  integ_bh_s = mp_aff.integrale_surface(source_bh_surf, radius_ah) ;

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  Scalar sou_bh1(mp_bh) ;
	  sou_bh1 = -2. * pow(lap_bh2,2.5) * mass * lldlconf / rr / rr ;
	  sou_bh1.std_spectral_base() ;
	  sou_bh1.inc_dzpuis(2) ;

	  Scalar sou_bh2(mp_bh) ;
	  sou_bh2 = 4.*mass*mass*pow(lap_bh2,3.)*(1.+3.*mass/rr)
	    *(1.+3.*mass/rr)*(lapse_bh_auto_rs+lapse_bh_comp)
	    *pow(confo_bh,4.)/3./pow(rr,4.) ;
	  sou_bh2.std_spectral_base() ;
	  sou_bh2.inc_dzpuis(4) ;

	  Scalar sou_bh3(mp_bh) ;
	  sou_bh3 = - 2.*mass*pow(lap_bh2,2.5)
	    *(2.+10.*mass/rr+9.*mass*mass/rr/rr)
	    *pow(confo_bh,4.)*(llshift+omelld)/pow(rr,3.) ;
	  sou_bh3.std_spectral_base() ;
	  sou_bh3.inc_dzpuis(4) ;

	  Scalar sou_bh4(mp_bh) ;
	  sou_bh4 = 2. * lap_bh2 * mass * lldlldlap / rr ;
	  sou_bh4.std_spectral_base() ;
	  sou_bh4.inc_dzpuis(1) ;

	  Scalar sou_bh5(mp_bh) ;
	  sou_bh5 = lap_bh2*lap_bh2*mass*(3.+8.*mass/rr)*lldlap/rr/rr ;
	  sou_bh5.std_spectral_base() ;
	  sou_bh5.inc_dzpuis(2) ;

	  Scalar sou_bh6(mp_bh) ;
	  sou_bh6 = 4.*pow(lap_bh2,3.5)*mass*mass
	    *(2.*(sqrt(lap_bh2)/lapse_bh - 1.)*pow(confo_bh,4.)
	      *(4.+12.*mass/rr+9.*mass*mass/rr/rr) + 3.*(pow(confo_bh,4.)-1.))
	    /3./pow(rr,4.) ;
	  sou_bh6.std_spectral_base() ;
	  sou_bh6.inc_dzpuis(4) ;

	  source_bh_volm = lapse_bh * (taij_quad_tot_rs + taij_quad_tot_rot)
	    / pow(confo_bh,8.)
	    - 2. * dlapdlcf + 4. * lap_bh2 * mass * lldlap * lldlconf / rr
	    + sou_bh1 + sou_bh2 + sou_bh3 + sou_bh4 + sou_bh5 + sou_bh6 ;

	  source_bh_volm.std_spectral_base() ;
	  source_bh_volm.annule_domain(0) ;

	  integ_bh_v = source_bh_volm.integrale() ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  source_ns_volm = lapse_ns * psi4_ns * (ener_euler + s_euler) ;
	  source_ns_volm.std_spectral_base() ;
	  source_ns_volm.inc_dzpuis(4) ;

	  integ_ns_v = source_ns_volm.integrale() ;

	  cout << "integ_bh_s : " << integ_bh_s/ qpig / msol
	       << "  integ_bh_v : "
	       << integ_bh_v/ qpig / msol
	       << "  integ_ns_v : " << integ_ns_v/ msol << endl ;

	  //--------------------
	  //     Komar mass
	  //--------------------
	  mkom = hole.get_mass_bh()
	    + (integ_bh_s + integ_bh_v) / qpig + integ_ns_v ;

	  cout << "Komar mass :         " << mkom / msol << " [Mo]"
	       << endl ;

	  // Komar mass by surface integral at infinity : dzpuis should be 2
	  // ------------------------------------------
	  double mmm = hole.get_mass_bh()
	    + (mp_aff.integrale_surface_infini(lldlap))/qpig ;

	  cout << "Another Komar mass : " << mmm / msol << " [Mo]"  << endl ;
	  */
	}
	else { // Isotropic coordinates with the maximal slicing

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Sets C/M^2 for each case of the lapse boundary condition
	  // --------------------------------------------------------
	  double cc ;

	  if (hole.has_bc_lapconf_nd()) {  // Neumann boundary condition
	      if (hole.has_bc_lapconf_fs()) {  // First condition
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
	      if (hole.has_bc_lapconf_fs()) {  // First condition
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
		  //	          cc = 2. * sqrt(2.) ;
	      }
	  }

	  Scalar r_are(mp_bh) ;
	  r_are = hole.r_coord(hole.has_bc_lapconf_nd(),
			       hole.has_bc_lapconf_fs()) ;
	  r_are.std_spectral_base() ;

	  Scalar lldlapconf_is(mp_bh) ;
	  lldlapconf_is = ll(1)%dlapconf_bh_auto_rs(1)
	    + ll(2)%dlapconf_bh_auto_rs(2) + ll(3)%dlapconf_bh_auto_rs(3)
	    + ll(1)%dlapconf_bh_comp(1) + ll(2)%dlapconf_bh_comp(2)
	    + ll(3)%dlapconf_bh_comp(3) ;
	  // dzpuis = 2
	  lldlapconf_is.std_spectral_base() ;

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  Scalar divshift_zero(divshift) ;
	  divshift_zero.dec_dzpuis(2) ;

	  Scalar lldllsh_zero(lldllsh) ;
	  lldllsh_zero.dec_dzpuis(2) ;

	  Scalar sou_bh_s_psi(mp_bh) ;
	  sou_bh_s_psi = 0.5 * confo_bh / rr
	    - pow(confo_bh, 4.) * (divshift_zero - 3.*lldllsh_zero)
	    / 12. / lapconf_bh
	    - 0.5 * pow(confo_bh, 4.) * mass * mass * cc
	    * sqrt(1. -2.*mass/r_are/rr + cc*cc*pow(mass/r_are/rr,4.))
	    / lapconf_bh / pow(r_are*rr,3.) ;

	  sou_bh_s_psi.std_spectral_base() ;
	  sou_bh_s_psi.annule_domain(0) ;
	  sou_bh_s_psi.raccord(1) ;

	  if (hole.has_bc_lapconf_nd()) {  // Neumann boundary condition
	      if (hole.has_bc_lapconf_fs()) {  // First condition

		source_bh_surf = sou_bh_s_psi ;

		source_bh_surf.std_spectral_base() ;
		source_bh_surf.annule_domain(0) ;
		source_bh_surf.raccord(1) ;

	      }
	      else {

		Scalar sou_bh_s1(mp_bh) ;
		sou_bh_s1 = 0.5 * lapconf_bh / rr ;
		sou_bh_s1.std_spectral_base() ;
		sou_bh_s1.annule_domain(0) ;
		sou_bh_s1.raccord(1) ;

		source_bh_surf = sou_bh_s1 + sou_bh_s_psi ;

		source_bh_surf.std_spectral_base() ;
		source_bh_surf.annule_domain(0) ;
		source_bh_surf.raccord(1) ;

	      }
	  }
	  else {

	      Scalar sou_bh_s1(mp_bh) ;
	      sou_bh_s1 = lldlapconf_is ;
	      sou_bh_s1.std_spectral_base() ;
	      sou_bh_s1.dec_dzpuis(2) ;

	      Scalar sou_bh_s2(mp_bh) ;
	      sou_bh_s2 = 0.5 * sqrt(r_are)
		* (1. - 3.*cc*cc*pow(mass/r_are/rr,4.)
		   -sqrt(1. - 2.*mass/r_are/rr
			 + cc*cc*pow(mass/r_are/rr,4.))) / rr ;

	      sou_bh_s2.std_spectral_base() ;

	      source_bh_surf = sou_bh_s1 + sou_bh_s2 + sou_bh_s_psi ;

	      source_bh_surf.std_spectral_base() ;
	      source_bh_surf.annule_domain(0) ;
	      source_bh_surf.raccord(1) ;

	  }

	  integ_bh_s = mp_aff.integrale_surface(source_bh_surf, radius_ah) ;

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  Scalar sou_bh1(mp_bh) ;
	  sou_bh1 = 0.75 * pow(mass*mass*cc,2.)
	    * (1. - 2.*mass/r_are/rr + cc*cc*pow(mass/r_are/rr,4.))
	    * (7.*pow(confo_bh,6.)/lapconf_bh
	       + pow(confo_bh,7.)/lapconf_bh/lapconf_bh)
	    / pow(r_are*rr,6.) ;

	  sou_bh1.std_spectral_base() ;
	  sou_bh1.inc_dzpuis(4) ;

	  Scalar sou_bh2(mp_bh) ;
	  sou_bh2 = 0.125 * (7.*lapconf_bh/pow(confo_bh,8.)
			     + 1./pow(confo_bh,7.)) * taij_quad_auto_bh ;

	  sou_bh2.std_spectral_base() ;

	  Scalar sou_bh3(mp_bh) ;
	  sou_bh3 = 0.125 * (7.*((lapconf_bh_comp+0.5)/lapconf_bh
				 * pow(confo_bh/(confo_bh_comp+0.5),6.) - 1.)
			     * (lapconf_bh_comp+0.5)
			     / pow(confo_bh_comp+0.5,8.)
			     + (pow(confo_bh/(confo_bh_comp+0.5),7.)
				*pow((lapconf_bh_comp+0.5)/lapconf_bh,2.)
				- 1.) / pow(confo_bh_comp+0.5,7))
	    * taij_quad_comp_bh ;

	  sou_bh3.std_spectral_base() ;

	  source_bh_volm = sou_bh1 + sou_bh2 + sou_bh3 ;
	  source_bh_volm.std_spectral_base() ;
	  source_bh_volm.annule_domain(0) ;

	  integ_bh_v = source_bh_volm.integrale() ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  Scalar sou_ns1(mp_ns) ;
	  sou_ns1 = 0.5 * pow(confo_ns,4.) * (lapconf_ns + confo_ns)
	    * ener_euler + lapconf_ns * pow(confo_ns,4.) * s_euler ;
	  sou_ns1.std_spectral_base() ;
	  sou_ns1.inc_dzpuis(4) ;

	  Scalar sou_ns2(mp_ns) ;
	  sou_ns2 = 0.125 * (7.*(lapconf_ns_auto+0.5)/(confo_ns_auto+0.5)
			     + 1.) * taij_quad_auto_ns
	    / pow(confo_ns_auto+0.5, 7.) / qpig ;
	  sou_ns2.std_spectral_base() ;

	  source_ns_volm = sou_ns1 + sou_ns2 ;
	  source_ns_volm.std_spectral_base() ;

	  integ_ns_v = source_ns_volm.integrale() ;

	  cout << "integ_bh_s : " << integ_bh_s/ qpig / msol
	       << "  integ_bh_v : "
	       << integ_bh_v/ qpig / msol
	       << "  integ_ns_v : " << integ_ns_v/ msol << endl ;

	  //--------------------
	  //     Komar mass
	  //--------------------
	  mkom = (integ_bh_s + integ_bh_v) / qpig + integ_ns_v ;

	  cout << "Komar mass (volume) :  " << mkom / msol << " [Mo]"
	       << endl ;

	}

	p_mass_kom_bhns_vol = new double( mkom ) ;

    }

    return *p_mass_kom_bhns_vol ;

}


               //-----------------------------------------//
               //          Total linear momentum          //
               //-----------------------------------------//

const Tbl& Bin_bhns::line_mom_bhns() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_line_mom_bhns == 0x0) {   // a new computation is required

        p_line_mom_bhns = new Tbl(3) ;
	p_line_mom_bhns->annule_hard() ;  // fills the double array with zeros

	if (hole.is_kerrschild()) {

	  // Construction of a new grid and a new affine mapping
	  // ---------------------------------------------------
	  int nz = (hole.get_mp()).get_mg()->get_nzone() ;
	  double* bornes = new double[nz+1] ;
	  double radius = separ + 2. ;

	  for (int i=nz-1; i>0; i--) {
	    bornes[i] = radius ;
	    radius /= 2. ;
	  }
	  bornes[0] = 0. ;
	  bornes[nz] = __infinity ;

	  Map_af mp_aff(*((hole.get_mp()).get_mg()), bornes) ;
	  mp_aff.set_ori(0.,0.,0.) ;

	  delete [] bornes ;

	  // Construction of the extrinsic curvature
	  // ---------------------------------------
	  Vector shift_bh(mp_aff, CON, mp_aff.get_bvect_cart()) ;
	  shift_bh.set_etat_qcq() ;

	  Vector shift_ns(mp_aff, CON, mp_aff.get_bvect_cart()) ;
	  shift_ns.set_etat_qcq() ;

	  shift_bh.set(1).import(hole.get_shift_auto()(1)) ;
	  shift_bh.set(2).import(hole.get_shift_auto()(2)) ;
	  shift_bh.set(3).import(hole.get_shift_auto()(3)) ;

	  shift_ns.set(1).import(star.get_shift_auto()(1)) ;
	  shift_ns.set(2).import(star.get_shift_auto()(2)) ;
	  shift_ns.set(3).import(star.get_shift_auto()(3)) ;

	  Vector shift_tot = shift_bh + shift_ns ;
	  shift_tot.std_spectral_base() ;
	  shift_tot.annule(0, nz-2) ;

	  Scalar stc(mp_aff) ;
	  stc = mp_aff.sint ;
	  stc.std_spectral_base() ;
	  Scalar ctc(mp_aff) ;
	  ctc = mp_aff.cost ;
	  ctc.std_spectral_base() ;
	  Scalar spc(mp_aff) ;
	  spc = mp_aff.sinp ;
	  spc.std_spectral_base() ;
	  Scalar cpc(mp_aff) ;
	  cpc = mp_aff.cosp ;
	  cpc.std_spectral_base() ;

	  Vector lc(mp_aff, CON, mp_aff.get_bvect_cart()) ;
	  lc.set_etat_qcq() ;
	  lc.set(1) = stc * cpc ;
	  lc.set(2) = stc * spc ;
	  lc.set(3) = ctc ;
	  lc.std_spectral_base() ;

	  Vector kijlj(mp_aff, CON, mp_aff.get_bvect_cart()) ;
	  kijlj.set_etat_qcq() ;

	  Scalar rr(mp_aff) ;
	  rr = mp_aff.r ;
	  rr.std_spectral_base() ;

	  Scalar xbhsr(mp_aff) ;
	  xbhsr = (hole.get_mp()).get_ori_x() / rr ;
	  xbhsr.std_spectral_base() ;

	  Scalar ybhsr(mp_aff) ;
	  ybhsr = (hole.get_mp()).get_ori_y() / rr ;
	  ybhsr.std_spectral_base() ;

	  Scalar rl(mp_aff) ;
	  rl = sqrt(1. - 2.*xbhsr*lc(1) - 2.*ybhsr*lc(2)
		    + xbhsr*xbhsr + ybhsr*ybhsr) ;
	  rl.std_spectral_base() ;

	  Vector ll(mp_aff, CON, mp_aff.get_bvect_cart()) ;
	  ll.set_etat_qcq() ;
	  ll.set(1) = (lc(1) - xbhsr) / rl ;
	  ll.set(2) = (lc(2) - ybhsr)/ rl ;
	  ll.set(3) = lc(3) / rl ;
	  ll.std_spectral_base() ;

	  Scalar lcll(mp_aff) ;
	  lcll = lc(1)*ll(1) + lc(2)*ll(2) + lc(3)*ll(3) ;
	  lcll.std_spectral_base() ;

	  Scalar divshift(mp_aff) ;
	  divshift = shift_tot(1).deriv(1) + shift_tot(2).deriv(2)
	    + shift_tot(3).deriv(3) ;
	  divshift.std_spectral_base() ;

	  Scalar llshift(mp_aff) ;
	  llshift = ll(1)*shift_tot(1) + ll(2)*shift_tot(2)
	    + ll(3)*shift_tot(3) ;
	  llshift.std_spectral_base() ;

	  Scalar lcshift(mp_aff) ;
	  lcshift = lc(1)*shift_tot(1) + lc(2)*shift_tot(2)
	    + lc(3)*shift_tot(3) ;
	  lcshift.std_spectral_base() ;

	  Vector lcdshft(mp_aff, CON, mp_aff.get_bvect_cart()) ;
	  for (int i=1; i<=3; i++) {
	    lcdshft.set(i) = lc(1)*(shift_tot(1).deriv(i))
	      + lc(2)*(shift_tot(2).deriv(i))
	      + lc(3)*(shift_tot(3).deriv(i)) ;
	  }
	  lcdshft.std_spectral_base() ;

	  Vector dshift(mp_aff, CON, mp_aff.get_bvect_cart()) ;
	  for (int i=1; i<=3; i++) {
	      dshift.set(i) = shift_tot(i).dsdr() ;
	  }
	  dshift.std_spectral_base() ;

	  Vector lldshft(mp_aff, CON, mp_aff.get_bvect_cart()) ;
	  for (int i=1; i<=3; i++) {
	    lldshft.set(i) = ll(1)*(shift_tot(i).deriv(1))
	      + ll(2)*(shift_tot(i).deriv(2))
	      + ll(3)*(shift_tot(i).deriv(3)) ;
	  }
	  lldshft.std_spectral_base() ;

	  double mass = ggrav * hole.get_mass_bh() ;

	  Scalar lap_bh2(mp_aff) ;
	  lap_bh2 = 1./(1.+2.*mass/rl/rr) ;
	  lap_bh2.std_spectral_base() ;

	  Scalar lap2hbh(mp_aff) ;
	  lap2hbh = lap_bh2 * mass / rl / rr ;
	  lap2hbh.std_spectral_base() ;

	  Scalar omexsr(mp_aff) ;
	  omexsr = omega * ((hole.get_mp()).get_ori_x() - x_rot)
	    / rl / rr ;
	  omexsr.std_spectral_base() ;

	  Scalar omeysr(mp_aff) ;
	  omeysr = omega * ((hole.get_mp()).get_ori_y() - y_rot)
	    / rl / rr ;
	  omeysr.std_spectral_base() ;

	  Scalar kk(mp_aff) ;
	  kk = 4.*sqrt(lap_bh2)*lap2hbh*(1.+3.*mass/rl/rr)/3./rl/rr ;
	  kk.std_spectral_base() ;

	  //-----------------------------------------------------------
	  //     Surface integral at infinity : dzpuis should be 2
	  //-----------------------------------------------------------

	  // Source for X component
	  // ----------------------
	  Scalar kij_x1(mp_aff) ;
	  kij_x1 = omexsr*ll(1)*lc(2) - omeysr*(ll(1)*lc(1)+lcll)
	    + lcll*shift_tot(1)/rl/rr
	    + ll(1)*lcshift/rl/rr
	    + (lc(1)-lap_bh2*(9.+14.*mass/rl/rr)*ll(1)*lcll)
	    *(llshift/rl/rr + omexsr*ll(2) - omeysr*ll(1))/3. ;
	  kij_x1.std_spectral_base() ;
	  kij_x1.inc_dzpuis(2) ;

	  Scalar kij_x2(mp_aff) ;
	  kij_x2 = kk * (lc(1) - 2.*lap2hbh*ll(1)*lcll) ;
	  kij_x2.std_spectral_base() ;
	  kij_x2.inc_dzpuis(2) ;

	  kijlj.set(1) = lcdshft(1) + dshift(1) - 2.*lc(1)*divshift/3.
	    + 2.*lap2hbh * (-ll(1)*(ll(1)*lcdshft(1) + ll(2)*lcdshft(2)
				    + ll(3)*lcdshft(3))
			    - lcll*lldshft(1)
			    + 2.*ll(1)*lcll*divshift/3.
			    + kij_x1)
	    + kij_x2 ;

	  // Source for Y component
	  // ----------------------
	  Scalar kij_y1(mp_aff) ;
	  kij_y1 = omexsr*(ll(2)*lc(2)+lcll) - omeysr*ll(2)*lc(1)
	    + lcll*shift_tot(2)/rl/rr
	    + ll(2)*lcshift/rl/rr
	    + (lc(2)-lap_bh2*(9.+14.*mass/rl/rr)*ll(2)*lcll)
	    *(llshift/rl/rr + omexsr*ll(2) - omeysr*ll(1))/3. ;
	  kij_y1.std_spectral_base() ;
	  kij_y1.inc_dzpuis(2) ;

	  Scalar kij_y2(mp_aff) ;
	  kij_y2 = kk * (lc(2) - 2.*lap2hbh*ll(2)*lcll) ;
	  kij_y2.std_spectral_base() ;
	  kij_y2.inc_dzpuis(2) ;

	  kijlj.set(2) = lcdshft(2) + dshift(2) - 2.*lc(2)*divshift/3.
	    + 2.*lap2hbh * (-ll(2)*(ll(1)*lcdshft(1) + ll(2)*lcdshft(2)
				    + ll(3)*lcdshft(3))
			    - lcll*lldshft(2)
			    + 2.*ll(2)*lcll*divshift/3.
			    + kij_y1)
	    + kij_y2 ;

	  // Source for Z component
	  // ----------------------
	  Scalar kij_z1(mp_aff) ;
	  kij_z1 = omexsr*ll(3)*lc(2) - omeysr*ll(3)*lc(1)
	    + lcll*shift_tot(3)/rl/rr
	    + ll(3)*lcshift/rl/rr
	    + (lc(3)-lap_bh2*(9.+14.*mass/rl/rr)*ll(3)*lcll)
	    *(llshift/rl/rr + omexsr*ll(2) - omeysr*ll(1))/3. ;
	  kij_z1.std_spectral_base() ;
	  kij_z1.inc_dzpuis(2) ;

	  Scalar kij_z2(mp_aff) ;
	  kij_z2 = kk * (lc(3) - 2.*lap2hbh*ll(3)*lcll) ;
	  kij_z2.std_spectral_base() ;
	  kij_z2.inc_dzpuis(2) ;

	  kijlj.set(3) = lcdshft(3) + dshift(3) - 2.*lc(3)*divshift/3.
	    + 2.*lap2hbh * (-ll(3)*(ll(1)*lcdshft(1) + ll(2)*lcdshft(2)
				    + ll(3)*lcdshft(3))
			    - lcll*lldshft(3)
			    + 2.*ll(3)*lcll*divshift/3.
			    + kij_z1)
	    + kij_z2 ;

	  kijlj.std_spectral_base() ;

	  // X component  dzpuis should be 2
	  // -----------
	  double lm_x = mp_aff.integrale_surface_infini(kijlj(1)) ;
	  p_line_mom_bhns->set(0) = 0.25 * lm_x / qpig ;

	  // Y component
	  // -----------
	  double lm_y = mp_aff.integrale_surface_infini(kijlj(2)) ;
	  p_line_mom_bhns->set(1) = 0.25 * lm_y / qpig ;

	  // Z component
	  // -----------
	  double lm_z = mp_aff.integrale_surface_infini(kijlj(3)) ;
	  p_line_mom_bhns->set(2) = 0.25 * lm_z / qpig ;

	}
	else {  // Isotropic coordinates with the maximal slicing

	  /*
	  // Method by the sourface integral at infinity
          // -------------------------------------------

	  const Map& mp_bh = hole.get_mp() ;
	  Map_af mp_aff(mp_bh) ;

	  Scalar st(mp_bh) ;
	  st = mp_bh.sint ;
	  st.std_spectral_base() ;
	  Scalar ct(mp_bh) ;
	  ct = mp_bh.cost ;
	  ct.std_spectral_base() ;
	  Scalar sp(mp_bh) ;
	  sp = mp_bh.sinp ;
	  sp.std_spectral_base() ;
	  Scalar cp(mp_bh) ;
	  cp = mp_bh.cosp ;
	  cp.std_spectral_base() ;

	  Vector ll(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	  ll.set_etat_qcq() ;
	  ll.set(1) = st * cp ;
	  ll.set(2) = st * sp ;
	  ll.set(3) = ct ;
	  ll.std_spectral_base() ;

	  const Scalar& confo_bh = hole.get_confo_tot() ;
	  const Sym_tensor& taij_tot_rs = hole.get_taij_tot_rs() ;

	  Vector kijlj(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	  kijlj.set_etat_qcq() ;

	  for (int i=1; i<=3; i++) {
	    kijlj.set(i) = (taij_tot_rs(i,1)%ll(1)
			    + taij_tot_rs(i,2)%ll(2)
			    + taij_tot_rs(i,3)%ll(3)) / pow(confo_bh,10.) ;
	  }

	  kijlj.std_spectral_base() ;

	  // X component
	  // -----------
	  double lm_x = mp_aff.integrale_surface_infini(kijlj(1)) ;
	  p_line_mom_bhns->set(0) = 0.5 * lm_x / qpig ;

	  // Y component
	  // -----------
	  double lm_y = mp_aff.integrale_surface_infini(kijlj(2)) ;
	  p_line_mom_bhns->set(1) = 0.5 * lm_y / qpig ;

	  // Z component
	  // -----------
	  double lm_z = mp_aff.integrale_surface_infini(kijlj(3)) ;
	  p_line_mom_bhns->set(2) = 0.5 * lm_z / qpig ;
	  */

	  // Method by the volume integral and the surface integral
	  // at the BH horizon
	  // ------------------------------------------------------

	  const Map& mp_bh = hole.get_mp() ;
	  Map_af mp_aff(mp_bh) ;
	  const Map& mp_ns = star.get_mp() ;

	  const Sym_tensor& taij = hole.get_taij_tot() ;
	  const Scalar& confo_ns = star.get_confo_tot() ;
	  const Scalar& ee = star.get_ener_euler() ;
	  const Scalar& pp = star.get_press() ;
	  const Vector& uu = star.get_u_euler() ;

	  double radius_ah = mp_bh.val_r(1,-1.,M_PI/2.,0.) ;

	  Scalar st(mp_bh) ;
	  st = mp_bh.sint ;
	  st.std_spectral_base() ;
	  Scalar ct(mp_bh) ;
	  ct = mp_bh.cost ;
	  ct.std_spectral_base() ;
	  Scalar sp(mp_bh) ;
	  sp = mp_bh.sinp ;
	  sp.std_spectral_base() ;
	  Scalar cp(mp_bh) ;
	  cp = mp_bh.cosp ;
	  cp.std_spectral_base() ;

	  Vector ll(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	  ll.set_etat_qcq() ;
	  ll.set(1) = st % cp ;
	  ll.set(2) = st % sp ;
	  ll.set(3) = ct ;
	  ll.std_spectral_base() ;

	  // X component
	  // -----------

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  Scalar sou_bh_sx(mp_bh) ;
	  sou_bh_sx = taij(1,1) * ll(1) + taij(1,2) * ll(2)
	    + taij(1,3) * ll(3) ;
	  sou_bh_sx.std_spectral_base() ;
	  sou_bh_sx.dec_dzpuis(2) ;  // dzpuis : 2 -> 0

	  sou_bh_sx.annule_domain(0) ;
	  sou_bh_sx.raccord(1) ;

	  double integ_bh_s_x = mp_aff.integrale_surface(sou_bh_sx,
							 radius_ah) ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  Scalar sou_ns_vx(mp_ns) ;

	  sou_ns_vx = pow(confo_ns, 10.) * (ee + pp) * uu(1) ;

	  sou_ns_vx.std_spectral_base() ;
	  sou_ns_vx.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	  double integ_ns_v_x = sou_ns_vx.integrale() ;

	  p_line_mom_bhns->set(0) = integ_ns_v_x + 0.5*integ_bh_s_x/qpig ;

	  // Y component
	  // -----------

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  Scalar sou_bh_sy(mp_bh) ;
	  sou_bh_sy = taij(2,1) * ll(1) + taij(2,2) * ll(2)
	    + taij(2,3) * ll(3) ;
	  sou_bh_sy.std_spectral_base() ;
	  sou_bh_sy.dec_dzpuis(2) ;  // dzpuis : 2 -> 0

	  sou_bh_sy.annule_domain(0) ;
	  sou_bh_sy.raccord(1) ;

	  double integ_bh_s_y = mp_aff.integrale_surface(sou_bh_sy,
							 radius_ah) ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  Scalar sou_ns_vy(mp_ns) ;

	  sou_ns_vy = pow(confo_ns, 10.) * (ee + pp) * uu(2) ;

	  sou_ns_vy.std_spectral_base() ;
	  sou_ns_vy.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	  double integ_ns_v_y = sou_ns_vy.integrale() ;

	  p_line_mom_bhns->set(1) = integ_ns_v_y + 0.5*integ_bh_s_y/qpig ;


	  // Z component
	  // -----------

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  Scalar sou_bh_sz(mp_bh) ;
	  sou_bh_sz = taij(3,1) * ll(1) + taij(3,2) * ll(2)
	    + taij(3,3) * ll(3) ;
	  sou_bh_sz.std_spectral_base() ;
	  sou_bh_sz.dec_dzpuis(2) ;  // dzpuis : 2 -> 0

	  sou_bh_sz.annule_domain(0) ;
	  sou_bh_sz.raccord(1) ;

	  double integ_bh_s_z = mp_aff.integrale_surface(sou_bh_sz,
							 radius_ah) ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  Scalar sou_ns_vz(mp_ns) ;

	  sou_ns_vz = pow(confo_ns, 10.) * (ee + pp) * uu(3) ;

	  sou_ns_vz.std_spectral_base() ;
	  sou_ns_vz.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	  double integ_ns_v_z = sou_ns_vz.integrale() ;

	  p_line_mom_bhns->set(2) = integ_ns_v_z + 0.5*integ_bh_s_z/qpig ;

	}

    }

    return *p_line_mom_bhns ;

}


               //------------------------------------------//
               //          Total angular momentum          //
               //------------------------------------------//

const Tbl& Bin_bhns::angu_mom_bhns() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_angu_mom_bhns == 0x0) {   // a new computation is required

        p_angu_mom_bhns = new Tbl(3) ;
	p_angu_mom_bhns->annule_hard() ;  // fills the double array with zeros

	double integ_bh_s_x ;
	double integ_bh_s_y ;
	double integ_bh_s_z ;
	double integ_bh_v_x ;
	double integ_bh_v_y ;
	double integ_bh_v_z ;
	double integ_ns_v_x ;
	double integ_ns_v_y ;
	double integ_ns_v_z ;

	const Map& mp_bh = hole.get_mp() ;
	const Map& mp_ns = star.get_mp() ;

	double radius_ah = mp_bh.val_r(1,-1.,M_PI/2.,0.) ;
	Map_af mp_aff(mp_bh) ;

	Scalar source_bh_surf(mp_bh) ;
	source_bh_surf.set_etat_qcq() ;

	Scalar source_bh_volm(mp_bh) ;
	source_bh_volm.set_etat_qcq() ;

	Scalar source_ns_volm(mp_ns) ;
	source_ns_volm.set_etat_qcq() ;

	Scalar rr(mp_bh) ;
	rr = mp_bh.r ;
	rr.std_spectral_base() ;

	Scalar st(mp_bh) ;
	st = mp_bh.sint ;
	st.std_spectral_base() ;
	Scalar ct(mp_bh) ;
	ct = mp_bh.cost ;
	ct.std_spectral_base() ;
	Scalar sp(mp_bh) ;
	sp = mp_bh.sinp ;
	sp.std_spectral_base() ;
	Scalar cp(mp_bh) ;
	cp = mp_bh.cosp ;
	cp.std_spectral_base() ;

	Vector ll(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	ll.set_etat_qcq() ;
	ll.set(1) = st % cp ;
	ll.set(2) = st % sp ;
	ll.set(3) = ct ;
	ll.std_spectral_base() ;

	Scalar xx_bh(mp_bh) ;
	xx_bh = mp_bh.xa ;
	xx_bh.std_spectral_base() ;
	Scalar yy_bh(mp_bh) ;
	yy_bh = mp_bh.ya ;
	yy_bh.std_spectral_base() ;
	Scalar zz_bh(mp_bh) ;
	zz_bh = mp_bh.za ;
	zz_bh.std_spectral_base() ;

	Scalar xx_ns(mp_ns) ;
	xx_ns = mp_ns.xa ;
	xx_ns.std_spectral_base() ;
	Scalar yy_ns(mp_ns) ;
	yy_ns = mp_ns.ya ;
	yy_ns.std_spectral_base() ;
	Scalar zz_ns(mp_ns) ;
	zz_ns = mp_ns.za ;
	zz_ns.std_spectral_base() ;

	const Scalar& confo_bh = hole.get_confo_tot() ;
	const Sym_tensor& taij = hole.get_taij_tot() ;
	const Scalar& confo_ns = star.get_confo_tot() ;
	const Scalar& ee = star.get_ener_euler() ;
	const Scalar& pp = star.get_press() ;
	const Vector& uu = star.get_u_euler() ;

	Vector jbh_x(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	jbh_x.set_etat_qcq() ;
	for (int i=1; i<=3; i++)
	  jbh_x.set(i) = yy_bh * taij(3,i) - zz_bh * taij(2,i) ;

	jbh_x.std_spectral_base() ;

	Vector jbh_y(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	jbh_y.set_etat_qcq() ;
	for (int i=1; i<=3; i++)
	  jbh_y.set(i) = zz_bh * taij(1,i) - xx_bh * taij(3,i) ;

	jbh_y.std_spectral_base() ;

	Vector jbh_z(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	jbh_z.set_etat_qcq() ;
	for (int i=1; i<=3; i++)
	  jbh_z.set(i) = xx_bh * taij(2,i) - yy_bh * taij(1,i) ;

	jbh_z.std_spectral_base() ;

	double mass = ggrav * hole.get_mass_bh() ;

	if (hole.is_kerrschild()) {

	  double ori_bh = mp_bh.get_ori_x() ;

	  Scalar lap_bh2(mp_bh) ;
	  lap_bh2 = 1./(1.+2.*mass/rr) ;
	  lap_bh2.std_spectral_base() ;

	  Scalar lcnf(mp_bh) ;
	  lcnf = log(confo_bh) ;
	  lcnf.std_spectral_base() ;

	  Vector jbhsr_x(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	  jbhsr_x.set_etat_qcq() ;
	  for (int i=1; i<=3; i++)
	    jbhsr_x.set(i) = ll(2)*taij(3,i) - ll(3)*taij(2,i) ;

	  jbhsr_x.std_spectral_base() ;

	  Vector jbhsr_y(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	  jbhsr_y.set_etat_qcq() ;
	  for (int i=1; i<=3; i++)
	    jbhsr_y.set(i) = ll(3)*taij(1,i) - (ll(1)+ori_bh/rr)*taij(3,i) ;

	  jbhsr_y.std_spectral_base() ;

	  Vector jbhsr_z(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	  jbhsr_z.set_etat_qcq() ;
	  for (int i=1; i<=3; i++)
	    jbhsr_z.set(i) = (ll(1)+ori_bh/rr)*taij(2,i) - ll(2)*taij(1,i) ;

	  jbhsr_z.std_spectral_base() ;

	  Scalar tmp1(mp_bh) ;  // dzpuis = 0
	  tmp1 = 2. * pow(lap_bh2,2.5) * mass * (1.+3.*mass/rr)
	    * pow(confo_bh,6.) * ori_bh / 3. / rr / rr ;
	  tmp1.std_spectral_base() ;
	  tmp1.annule_domain(0) ;

	  Scalar tmp2(mp_bh) ;  // dzpuis = 0
	  tmp2 = 4. * sqrt(lap_bh2) * (1.+3.*mass/rr) * pow(confo_bh,6.) ;
	  tmp2.std_spectral_base() ;
	  tmp2.annule_domain(0) ;

	  Scalar lltaij(mp_bh) ;  // dzpuis = 2
	  lltaij = ll(1)*(ll(1)*taij(1,1)+ll(2)*taij(1,2)+ll(3)*taij(1,3))
	    + ll(2)*(ll(1)*taij(2,1)+ll(2)*taij(2,2)+ll(3)*taij(2,3))
	    + ll(3)*(ll(1)*taij(3,1)+ll(2)*taij(3,2)+ll(3)*taij(3,3)) ;
	  lltaij.std_spectral_base() ;
	  lltaij.dec_dzpuis(2) ; // dzpuis : 2 -> 0

	  Scalar dlcnf(mp_bh) ;   // dzpuis = 2
	  dlcnf = - 2. * lap_bh2 * tmp2 * mass * lcnf.dsdr() / rr ;
	  dlcnf.std_spectral_base() ;
	  dlcnf.dec_dzpuis(2) ;  // dzpuis : 2 -> 0
	  dlcnf.annule_domain(0) ;

	  Scalar tmp3(mp_bh) ;  // dzpuis = 0
	  tmp3 = -2.*pow(lap_bh2,2.5)
	    *(6.+32.*mass/rr+41.*mass*mass/rr/rr+24.*pow(mass,3.)/pow(rr,3.))
	    *pow(confo_bh,6.)/3./rr
	    + 3.* lltaij + dlcnf ;
	  tmp3.std_spectral_base() ;
	  tmp3.annule_domain(0) ;

	  Scalar tmp4(mp_bh) ;  // dzpuis = 0
	  tmp4 = lap_bh2 * mass / rr ;
	  tmp4.std_spectral_base() ;

	  //-------------//
	  // X component //
	  //-------------//

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  source_bh_surf = jbh_x(1)*ll(1) + jbh_x(2)*ll(2) + jbh_x(3)*ll(3) ;

	  source_bh_surf.std_spectral_base() ;
	  source_bh_surf.dec_dzpuis(2) ;  // dzpuis : 2 -> 0
	  source_bh_surf.annule_domain(0) ;
	  source_bh_surf.raccord(1) ;

	  integ_bh_s_x = mp_aff.integrale_surface(source_bh_surf, radius_ah) ;

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  source_bh_volm = tmp4
	    * ( jbhsr_x(1)*ll(1) + jbhsr_x(2)*ll(2) + jbhsr_x(3)*ll(3)
	       + tmp2 * ( ll(2)*lcnf.dsdz() - ll(3)*lcnf.dsdy() ) ) ;

	  source_bh_volm.std_spectral_base() ;
	  source_bh_volm.inc_dzpuis(2) ;  // dzpuis : 2 -> 4
	  source_bh_volm.annule_domain(0) ;

	  integ_bh_v_x = source_bh_volm.integrale() ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  source_ns_volm = pow(confo_ns, 10.) * (ee + pp)
	    * (yy_ns*uu(3) - zz_ns*uu(2)) ;

	  source_ns_volm.std_spectral_base() ;
	  source_ns_volm.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	  integ_ns_v_x = source_ns_volm.integrale() ;

	  p_angu_mom_bhns->set(0) = integ_ns_v_x
	    + 0.5 * (integ_bh_s_x + integ_bh_v_x) / qpig ;

	  //-------------//
	  // Y component //
	  //-------------//

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  Scalar jbh_y_ll(mp_bh) ;
	  jbh_y_ll = jbh_y(1)*ll(1) + jbh_y(2)*ll(2) + jbh_y(3)*ll(3) ;
	  jbh_y_ll.std_spectral_base() ;
	  jbh_y_ll.dec_dzpuis(2) ;  // dzpuis : 2 -> 0

	  source_bh_surf = jbh_y_ll - tmp1 * ll(3) ;
	  source_bh_surf.std_spectral_base() ;
	  source_bh_surf.annule_domain(0) ;
	  source_bh_surf.raccord(1) ;

	  integ_bh_s_y = mp_aff.integrale_surface(source_bh_surf, radius_ah) ;

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  Scalar tmp3_llz(mp_bh) ;
	  tmp3_llz = tmp3 * ll(3) * ori_bh / rr ;
	  tmp3_llz.std_spectral_base() ;
	  tmp3_llz.inc_dzpuis(2) ;  // dzpuis : 0 -> 2

	  source_bh_volm = tmp4
	    * ( jbhsr_y(1)*ll(1) + jbhsr_y(2)*ll(2) + jbhsr_y(3)*ll(3)
		+ tmp2 * ( ll(3)*lcnf.dsdx() - (ll(1)+ori_bh/rr)*lcnf.dsdz() )
		- tmp3_llz ) ;

	  source_bh_volm.std_spectral_base() ;
	  source_bh_volm.inc_dzpuis(2) ;  // dzpuis : 2 -> 4
	  source_bh_volm.annule_domain(0) ;

	  integ_bh_v_y = source_bh_volm.integrale() ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  source_ns_volm = pow(confo_ns, 10.) * (ee + pp)
	    * (zz_ns*uu(1) - xx_ns*uu(3)) ;

	  source_ns_volm.std_spectral_base() ;
	  source_ns_volm.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	  integ_ns_v_y = source_ns_volm.integrale() ;

	  p_angu_mom_bhns->set(1) = integ_ns_v_y
	    + 0.5 * (integ_bh_s_y + integ_bh_v_y) / qpig ;

	  //-------------//
	  // Z component //
	  //-------------//

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  Scalar jbh_z_ll(mp_bh) ;
	  jbh_z_ll = jbh_z(1)*ll(1) + jbh_z(2)*ll(2) + jbh_z(3)*ll(3) ;
	  jbh_z_ll.std_spectral_base() ;
	  jbh_z_ll.dec_dzpuis(2) ;  // dzpuis : 2 -> 0
	  source_bh_surf = jbh_z_ll + tmp1 * ll(2) ;

	  source_bh_surf.std_spectral_base() ;
	  source_bh_surf.annule_domain(0) ;
	  source_bh_surf.raccord(1) ;

	  integ_bh_s_z = mp_aff.integrale_surface(source_bh_surf, radius_ah) ;

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  Scalar tmp3_lly(mp_bh) ;
	  tmp3_lly = tmp3 * ll(2) * ori_bh / rr ;
	  tmp3_lly.std_spectral_base() ;
	  tmp3_lly.inc_dzpuis(2) ;  // dzpuis : 0 -> 2

	  source_bh_volm = tmp4
	    * ( jbhsr_z(1)*ll(1) + jbhsr_z(2)*ll(2) + jbhsr_z(3)*ll(3)
		+ tmp2 * ( (ll(1)+ori_bh/rr)*lcnf.dsdy() - ll(2)*lcnf.dsdx() )
		+ tmp3_lly ) ;

	  source_bh_volm.std_spectral_base() ;
	  source_bh_volm.inc_dzpuis(2) ;  // dzpuis : 2 -> 4
	  source_bh_volm.annule_domain(0) ;

	  integ_bh_v_z = source_bh_volm.integrale() ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  source_ns_volm = pow(confo_ns, 10.) * (ee + pp)
	    * (xx_ns*uu(2) - yy_ns*uu(1)) ;

	  source_ns_volm.std_spectral_base() ;
	  source_ns_volm.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	  integ_ns_v_z = source_ns_volm.integrale() ;

	  p_angu_mom_bhns->set(2) = integ_ns_v_z
	    + 0.5 * (integ_bh_s_z + integ_bh_v_z) / qpig ;

	}
	else {  // Isotropic coordinates with the maximal slicing

	  // Sets C/M^2 for each case of the lapse boundary condition
	  // --------------------------------------------------------
	  double cc ;

	  if (hole.has_bc_lapconf_nd()) {  // Neumann boundary condition
	      if (hole.has_bc_lapconf_fs()) {  // First condition
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
	      if (hole.has_bc_lapconf_fs()) {  // First condition
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
		  //	          cc = 2. * sqrt(2.) ;
	      }
	  }

	  Scalar r_are(mp_bh) ;
	  r_are = hole.r_coord(hole.has_bc_lapconf_nd(),
			       hole.has_bc_lapconf_fs()) ;
	  r_are.std_spectral_base() ;

	  const Scalar& lapconf_bh = hole.get_lapconf_tot() ;
	  const Scalar& confo_bh = hole.get_confo_tot() ;
	  const Sym_tensor& taij_rs = hole.get_taij_tot_rs() ;

	  Vector jbh_rs_x(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	  jbh_rs_x.set_etat_qcq() ;
	  for (int i=1; i<=3; i++)
	    jbh_rs_x.set(i) = yy_bh * taij_rs(3,i) - zz_bh * taij_rs(2,i) ;

	  jbh_rs_x.std_spectral_base() ;

	  Vector jbh_rs_y(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	  jbh_rs_y.set_etat_qcq() ;
	  for (int i=1; i<=3; i++)
	    jbh_rs_y.set(i) = zz_bh * taij_rs(1,i) - xx_bh * taij_rs(3,i) ;

	  jbh_rs_y.std_spectral_base() ;

	  Vector jbh_rs_z(mp_bh, CON, mp_bh.get_bvect_cart()) ;
	  jbh_rs_z.set_etat_qcq() ;
	  for (int i=1; i<=3; i++)
	    jbh_rs_z.set(i) = xx_bh * taij_rs(2,i) - yy_bh * taij_rs(1,i) ;

	  jbh_rs_z.std_spectral_base() ;

	  Scalar tmp(mp_bh) ;
	  tmp = - 2. * mass * mass * cc * pow(confo_bh,7.)
	    * sqrt(1. - 2.*mass/r_are/rr + cc*cc*pow(mass/r_are/rr,4.))
	    / lapconf_bh / pow(r_are*rr,3.) ;
	  tmp.std_spectral_base() ;

	  //-------------//
	  // X component //
	  //-------------//

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  Scalar sou_bh_sx1(mp_bh) ;
	  sou_bh_sx1 = jbh_rs_x(1)%ll(1) + jbh_rs_x(2)%ll(2)
	    + jbh_rs_x(3)%ll(3) ;
	  sou_bh_sx1.std_spectral_base() ;
	  sou_bh_sx1.dec_dzpuis(2) ;  // dzpuis : 2 -> 0

	  Scalar sou_bh_sx2(mp_bh) ;
	  sou_bh_sx2 = tmp * (yy_bh % ll(3) - zz_bh % ll(2)) ;
	  sou_bh_sx2.std_spectral_base() ;

	  source_bh_surf = sou_bh_sx1 + sou_bh_sx2 ;

	  source_bh_surf.std_spectral_base() ;
	  source_bh_surf.annule_domain(0) ;
	  source_bh_surf.raccord(1) ;

	  integ_bh_s_x = mp_aff.integrale_surface(source_bh_surf, radius_ah) ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  source_ns_volm = pow(confo_ns, 10.) * (ee + pp)
	    * (yy_ns*uu(3) - zz_ns*uu(2)) ;

	  source_ns_volm.std_spectral_base() ;
	  source_ns_volm.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	  integ_ns_v_x = source_ns_volm.integrale() ;

	  p_angu_mom_bhns->set(0) = integ_ns_v_x + 0.5*integ_bh_s_x/qpig ;


	  //-------------//
	  // Y component //
	  //-------------//

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  Scalar sou_bh_sy1(mp_bh) ;
	  sou_bh_sy1 = jbh_rs_y(1)%ll(1) + jbh_rs_y(2)%ll(2)
	    + jbh_rs_y(3)%ll(3) ;
	  sou_bh_sy1.std_spectral_base() ;
	  sou_bh_sy1.dec_dzpuis(2) ;  // dzpuis : 2 -> 0

	  Scalar sou_bh_sy2(mp_bh) ;
	  sou_bh_sy2 = tmp * (zz_bh % ll(1) - xx_bh % ll(3)) ;
	  sou_bh_sy2.std_spectral_base() ;

	  source_bh_surf = sou_bh_sy1 + sou_bh_sy2 ;

	  source_bh_surf.std_spectral_base() ;
	  source_bh_surf.annule_domain(0) ;
	  source_bh_surf.raccord(1) ;

	  integ_bh_s_y = mp_aff.integrale_surface(source_bh_surf, radius_ah) ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  source_ns_volm = pow(confo_ns, 10.) * (ee + pp)
	    * (zz_ns*uu(1) - xx_ns*uu(3)) ;

	  source_ns_volm.std_spectral_base() ;
	  source_ns_volm.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	  integ_ns_v_y = source_ns_volm.integrale() ;

	  p_angu_mom_bhns->set(1) = integ_ns_v_y + 0.5*integ_bh_s_y/qpig ;


	  //-------------//
	  // Z component //
	  //-------------//

	  //-------------------------------------
	  //     Integration over the BH map
	  //-------------------------------------

	  // Surface integral <- dzpuis should be 0
	  // --------------------------------------
	  Scalar sou_bh_sz1(mp_bh) ;
	  sou_bh_sz1 = jbh_rs_z(1)%ll(1) + jbh_rs_z(2)%ll(2)
	    + jbh_rs_z(3)%ll(3) ;
	  sou_bh_sz1.std_spectral_base() ;
	  sou_bh_sz1.dec_dzpuis(2) ;  // dzpuis : 2 -> 0

	  Scalar sou_bh_sz2(mp_bh) ;
	  sou_bh_sz2 = tmp * (xx_bh % ll(2) - yy_bh % ll(1)) ;
	  sou_bh_sz2.std_spectral_base() ;

	  source_bh_surf = sou_bh_sz1 + sou_bh_sz2 ;

	  source_bh_surf.std_spectral_base() ;
	  source_bh_surf.annule_domain(0) ;
	  source_bh_surf.raccord(1) ;

	  integ_bh_s_z = mp_aff.integrale_surface(source_bh_surf, radius_ah) ;

	  //-------------------------------------
	  //     Integration over the NS map
	  //-------------------------------------

	  // Volume integral <- dzpuis should be 4
	  // -------------------------------------
	  source_ns_volm = pow(confo_ns, 10.) * (ee + pp)
	    * (xx_ns*uu(2) - yy_ns*uu(1)) ;

	  source_ns_volm.std_spectral_base() ;
	  source_ns_volm.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	  integ_ns_v_z = source_ns_volm.integrale() ;

	  p_angu_mom_bhns->set(2) = integ_ns_v_z + 0.5*integ_bh_s_z/qpig ;

      }

    }

    return *p_angu_mom_bhns ;

}


               //-----------------------------------------------//
               //          Error on the virial theorem          //
               //-----------------------------------------------//

double Bin_bhns::virial_bhns_surf() const {

    if (p_virial_bhns_surf == 0x0) {   // a new computation is required

      double virial = 1. - mass_kom_bhns_surf() / mass_adm_bhns_surf() ;

      p_virial_bhns_surf = new double( virial ) ;

    }

    return *p_virial_bhns_surf ;

}


double Bin_bhns::virial_bhns_vol() const {

    if (p_virial_bhns_vol == 0x0) {   // a new computation is required

      double virial = 1. - mass_kom_bhns_vol() / mass_adm_bhns_vol() ;

      p_virial_bhns_vol = new double( virial ) ;

    }

    return *p_virial_bhns_vol ;

}

               //--------------------------------------------------//
               //          X coordinate of the barycenter          //
               //--------------------------------------------------//

double Bin_bhns::xa_barycenter() const {

    using namespace Unites ;

    if (p_xa_barycenter == 0x0) {    // a new computation is required

        double mass_b = star.mass_b_bhns(hole.is_kerrschild(),
					 hole.get_mass_bh(), separ) ;

	const Map& mp = star.get_mp() ;
	Scalar xxa(mp) ;
	xxa = mp.xa ;   // Absolute X coordinate
	xxa.std_spectral_base() ;

	Scalar tmp(mp) ;

	if (hole.is_kerrschild()) {

	    Scalar xx(mp) ;
	    xx = mp.x ;
	    xx.std_spectral_base() ;
	    Scalar yy(mp) ;
	    yy = mp.y ;
	    yy.std_spectral_base() ;
	    Scalar zz(mp) ;
	    zz = mp.z ;
	    zz.std_spectral_base() ;

	    double yns = (star.get_mp()).get_ori_y() ;

	    Scalar rbh(mp) ;
	    rbh = sqrt( (xx+separ)*(xx+separ) + (yy+yns)*(yy+yns) + zz*zz ) ;
	    rbh.std_spectral_base() ;

	    double mass = ggrav * hole.get_mass_bh() ;

	    tmp = sqrt( 1.+2.*mass/rbh ) ;
	    tmp.std_spectral_base() ;

	}
	else { // Isotropic coordinates with the maximal slicing

	    tmp = 1. ;
	    tmp.std_spectral_base() ;

	}

	Scalar confo = star.get_confo_tot() ;
	confo.std_spectral_base() ;

	Scalar g_euler = star.get_gam_euler() ;
	g_euler.std_spectral_base() ;

	Scalar nbary = star.get_nbar() ;
	nbary.std_spectral_base() ;

	Scalar dens = pow(confo, 6.) * g_euler * nbary * tmp * xxa ;
	dens.std_spectral_base() ;
	dens.inc_dzpuis(4) ;
	assert(dens.get_dzpuis() == 4) ;

	double xa_bary = dens.integrale() / mass_b ;

	p_xa_barycenter = new double( xa_bary ) ;

    }

    return *p_xa_barycenter ;

}


               //--------------------------------------------------//
               //          Y coordinate of the barycenter          //
               //--------------------------------------------------//

double Bin_bhns::ya_barycenter() const {

    using namespace Unites ;

    if (p_ya_barycenter == 0x0) {    // a new computation is required

        double mass_b = star.mass_b_bhns(hole.is_kerrschild(),
					 hole.get_mass_bh(), separ) ;

	const Map& mp = star.get_mp() ;
	Scalar yya(mp) ;
	yya = mp.ya ;   // Absolute Y coordinate
	yya.std_spectral_base() ;

	Scalar tmp(mp) ;

	if (hole.is_kerrschild()) {

	    Scalar xx(mp) ;
	    xx = mp.x ;
	    xx.std_spectral_base() ;
	    Scalar yy(mp) ;
	    yy = mp.y ;
	    yy.std_spectral_base() ;
	    Scalar zz(mp) ;
	    zz = mp.z ;
	    zz.std_spectral_base() ;

	    double yns = (star.get_mp()).get_ori_y() ;

	    Scalar rbh(mp) ;
	    rbh = sqrt( (xx+separ)*(xx+separ) + (yy+yns)*(yy+yns) + zz*zz ) ;
	    rbh.std_spectral_base() ;

	    double mass = ggrav * hole.get_mass_bh() ;

	    tmp = sqrt( 1.+2.*mass/rbh ) ;
	    tmp.std_spectral_base() ;

	}
	else { // Isotropic coordinates with the maximal slicing

	    tmp = 1. ;
	    tmp.std_spectral_base() ;

	}

	Scalar confo = star.get_confo_tot() ;
	confo.std_spectral_base() ;

	Scalar g_euler = star.get_gam_euler() ;
	g_euler.std_spectral_base() ;

	Scalar nbary = star.get_nbar() ;
	nbary.std_spectral_base() ;

	Scalar dens = pow(confo, 6.) * g_euler * nbary * tmp * yya ;
	dens.std_spectral_base() ;
	dens.inc_dzpuis(4) ;
	assert(dens.get_dzpuis() == 4) ;

	double ya_bary = dens.integrale() / mass_b ;

	p_ya_barycenter = new double( ya_bary ) ;

    }

    return *p_ya_barycenter ;

}
