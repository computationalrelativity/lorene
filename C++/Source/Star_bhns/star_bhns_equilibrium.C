/*
 *  Method of class Star_bhns to compute an equilibrium configuration
 *   of a neutron star in a black hole-neutron star binary
 *
 *    (see file star_bhns.h for documentation).
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

char star_bhns_equilibrium_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:13:16  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/05/15 19:13:45  k_taniguchi
 * Change of some parameters.
 *
 * Revision 1.1  2007/06/22 01:30:45  k_taniguchi
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
#include "star_bhns.h"
#include "param.h"
#include "cmp.h"
#include "tenseur.h"
#include "utilitaires.h"
#include "unites.h"
//#include "graphique.h"

void Star_bhns::equilibrium_bhns(double ent_c, const double& mass_bh,
				 const double& sepa, bool kerrschild,
				 int mer, int mermax_ns, int mermax_potvit,
				 int mermax_poisson, int filter_r,
				 int filter_r_s, int filter_p_s,
				 double relax_poisson, double relax_potvit,
				 double thres_adapt, double resize_ns,
				 const Tbl& fact_resize, Tbl& diff) {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    // Initializations
    // ---------------

    const Mg3d* mg = mp.get_mg() ;
    int nz = mg->get_nzone() ;    // Total number of domain
    //    int nt = mg->get_nt(0) ;

    // The following is required to initialize mp_prev as a Map_et:
    Map_et& mp_et = dynamic_cast<Map_et&>(mp) ;

    // Domain and radial indices of points at the surface of the star:
    int l_b = nzet - 1 ;
    int i_b = mg->get_nr(l_b) - 1 ;
    int k_b ;
    int j_b ;

    // Value of the enthalpy defining the surface of the star
    double ent_b = 0 ;

    // Error indicators
    // ----------------

    double& diff_ent = diff.set(0) ;
    double& diff_vel_pot = diff.set(1) ;
    double& diff_lapconf = diff.set(2) ;
    double& diff_confo = diff.set(3) ;
    double& diff_shift_x = diff.set(4) ;
    double& diff_shift_y = diff.set(5) ;
    double& diff_shift_z = diff.set(6) ;
    double& diff_dHdr = diff.set(7) ;
    double& diff_dHdr_min = diff.set(8) ;
    double& diff_phi_min = diff.set(9) ;
    double& diff_radius = diff.set(10) ;

    // Parameters for the function Map_et::adapt
    // -----------------------------------------

    Param par_adapt ;
    int nitermax = 100 ;
    int niter ;
    int adapt_flag = 1 ;   //  1 = performs the full computation,
                           //  0 = performs only the rescaling by
                           //      the factor alpha_r

    // Number of domains for searching the enthalpy isosurfaces
    int nz_search = nzet + 1 ;
    //    int nz_search = nzet ;

    double precis_secant = 1.e-14 ;
    double alpha_r ;
    double reg_map = 1. ; // 1 = regular mapping, 0 = contracting mapping

    Tbl ent_limit(nz) ;

    par_adapt.add_int(nitermax, 0) ;   // maximum number of iterations to
                                       // locate zeros by the secant method
    par_adapt.add_int(nzet, 1) ;       // number of domains where the
                                       // adjustment to the isosurfaces of
                                       // ent is to be performed
    par_adapt.add_int(nz_search, 2) ;  // number of domains to search for
                                       // the enthalpy isosurface
    par_adapt.add_int(adapt_flag, 3) ; // 1 = performs the full computation,
                                       // 0 = performs only the rescaling by
                                       //     the factor alpha_r
    par_adapt.add_int(j_b, 4) ;        // theta index of the collocation point
                                       // (theta_*, phi_*)
    par_adapt.add_int(k_b, 5) ;        // theta index of the collocation point
                                       // (theta_*, phi_*)
    par_adapt.add_int_mod(niter, 0) ;  // number of iterations actually used
                                       // in the secant method
    par_adapt.add_double(precis_secant, 0) ; // required absolute precision in
                                             // the determination of zeros by
                                             // the secant method
    par_adapt.add_double(reg_map, 1) ; // 1. = regular mapping,
                                       // 0 = contracting mapping
    par_adapt.add_double(alpha_r, 2) ; // factor by which all the radial
                                       // distances will be multiplied
    par_adapt.add_tbl(ent_limit, 0) ;  // array of values of the field ent
                                       // to define the isosurfaces

    Cmp ssjm1lapconf(ssjm1_lapconf) ;
    Cmp ssjm1confo(ssjm1_confo) ;

    ssjm1lapconf.set_etat_qcq() ;
    ssjm1confo.set_etat_qcq() ;

    double precis_poisson = 1.e-14 ;

    // Parameters for the function Scalar::poisson for lapconf_auto
    // ------------------------------------------------------------

    Param par_lapconf ;

    par_lapconf.add_int(mermax_poisson, 0) ;    // maximum number of iterations
    par_lapconf.add_double(relax_poisson, 0) ;  // relaxation parameter
    par_lapconf.add_double(precis_poisson, 1) ; // required precision
    par_lapconf.add_int_mod(niter, 0) ;  // number of iterations actually used
    par_lapconf.add_cmp_mod( ssjm1lapconf ) ;

    // Parameters for the function Scalar::poisson for confo_auto
    // ----------------------------------------------------------

    Param par_confo ;

    par_confo.add_int(mermax_poisson, 0) ;    // maximum number of iterations
    par_confo.add_double(relax_poisson, 0) ;  // relaxation parameter
    par_confo.add_double(precis_poisson, 1) ; // required precision
    par_confo.add_int_mod(niter, 0) ;  // number of iterations actually used
    par_confo.add_cmp_mod( ssjm1confo ) ;

    // Parameters for the function Vector::poisson for shift method 2
    // --------------------------------------------------------------

    Param par_shift2 ;

    par_shift2.add_int(mermax_poisson, 0) ;    // maximum number of iterations
    par_shift2.add_double(relax_poisson, 0) ;  // relaxation parameter
    par_shift2.add_double(precis_poisson, 1) ; // required precision
    par_shift2.add_int_mod(niter, 0) ; // number of iterations actually used

    Cmp ssjm1khi(ssjm1_khi) ;
    ssjm1khi.set_etat_qcq() ;

    Tenseur ssjm1wshift(mp, 1, CON, mp.get_bvect_cart()) ;
    ssjm1wshift.set_etat_qcq() ;
    for (int i=0; i<3; i++) {
        ssjm1wshift.set(i) = Cmp(ssjm1_wshift(i+1)) ;
    }

    par_shift2.add_cmp_mod(ssjm1khi) ;
    par_shift2.add_tenseur_mod(ssjm1wshift) ;

    Scalar ent_jm1 = ent ;  // Enthalpy at previous step

    Scalar lapconf_m1(mp) ;  // = lapconf_auto - 0.5
    Scalar confo_m1(mp) ;  // = confo_auto - 0.5

    Scalar source_lapconf(mp) ; // Source term in the equation for lapconf_auto
    source_lapconf.set_etat_qcq() ;
    Scalar source_confo(mp) ; // Source term in the equation for confo_auto
    source_confo.set_etat_qcq() ;
    Vector source_shift(mp, CON, mp.get_bvect_cart()) ; // Source term 
                            // in the equation for shift_auto
    source_shift.set_etat_qcq() ;

    //===========================================================
    //                    Start of iteration
    //===========================================================

    for (int mer_ns=0; mer_ns<mermax_ns; mer_ns++) {

        cout << "-----------------------------------------------" << endl ;
	cout << "step: " << mer_ns << endl ;
	cout << "diff_ent = " << diff_ent << endl ;

	//------------------------------------------------------
	// Resolution of the elliptic equation for the velocity
	// scalar potential
	//------------------------------------------------------

	if (irrotational) {
	    diff_vel_pot = velo_pot_bhns(mass_bh, sepa, kerrschild,
					 mermax_potvit, precis_poisson,
					 relax_potvit) ;
	}
	else {
	    diff_vel_pot = 0. ; // to avoid the warning
	}

	//-------------------------------------
	// Computation of the new radial scale
	//-------------------------------------

	// alpha_r (r = alpha_r r') is determined so that the enthalpy
	// takes the requested value ent_b at the stellar surface

	// Values at the center of the star:
	// lapconf_auto = lapconf_auto=0(r->infty) + 0.5
	// The term lapconf_auto=0(r->infty) should be rescaled.
	double lapconf_auto_c = lapconf_auto.val_grid_point(0,0,0,0) - 0.5 ;
	double lapconf_comp_c = lapconf_comp.val_grid_point(0,0,0,0) ;

	double confo_c = confo_tot.val_grid_point(0,0,0,0) ;

	double gam_c = gam.val_grid_point(0,0,0,0) ;
	double gam0_c = gam0.val_grid_point(0,0,0,0) ;

	double hhh_c = exp(ent_c) ;
	double hhh_b = exp(ent_b) ;

	// Search for the reference point (theta_*, phi_*) [notation of
	//  Bonazzola, Gourgoulhon & Marck PRD 58, 104020 (1998)]
	//  at the surface of the star
	// ------------------------------------------------------------
	double alpha_r2 = 0. ;

	for (int k=0; k<mg->get_np(l_b); k++) {
	    for (int j=0; j<mg->get_nt(l_b); j++) {

		double lapconf_auto_b =
		  lapconf_auto.val_grid_point(l_b,k,j,i_b) - 0.5 ;
		double lapconf_comp_b =
		  lapconf_comp.val_grid_point(l_b,k,j,i_b) ;

		double confo_b = confo_tot.val_grid_point(l_b,k,j,i_b) ;

		double gam_b = gam.val_grid_point(l_b,k,j,i_b) ;
		double gam0_b = gam0.val_grid_point(l_b,k,j,i_b) ;

		double aaa = (gam0_c*gam_b*hhh_b*confo_c)
		  / (gam0_b*gam_c*hhh_c*confo_b) ;

		// See Eq (100) from Gourgoulhon et al. (2001)
		double alpha_r2_jk = (aaa * lapconf_comp_b - lapconf_comp_c
				      + 0.5 * (aaa - 1.))
		  / (lapconf_auto_c - aaa * lapconf_auto_b ) ;

		if (alpha_r2_jk > alpha_r2) {
		    alpha_r2 = alpha_r2_jk ;
		    k_b = k ;
		    j_b = j ;
		}
	    }
	}

	alpha_r = sqrt(alpha_r2) ;

	cout << "k_b, j_b, alpha_r: " << k_b << "  " << j_b << "  "
	     << alpha_r << endl ;

	// New value of lapconf_auto
	// -------------------------

	lapconf_auto = alpha_r2 * (lapconf_auto - 0.5) + 0.5 ;
	Scalar lapconf_tot_tmp = lapconf_auto + lapconf_comp ;
	lapconf_tot_tmp.std_spectral_base() ;

	/*
	confo_auto = alpha_r2 * (confo_auto - 0.5) + 0.5 ;
	Scalar confo_tot_tmp = confo_auto + confo_comp ;
	confo_tot_tmp.std_spectral_base() ;
	*/
	//------------------------------------------------------------
	// Change the values of the inner points of the second domain
	// by those of the outer points of the first domain
	//------------------------------------------------------------

	lapconf_auto.set_spectral_va().smooth(nzet,lapconf_auto.set_spectral_va()) ;

	//--------------------------------------------
	// First integral  -->  enthalpy in all space
	// See Eq (98) from Gourgoulhon et al. (2001)
	//--------------------------------------------

	double log_lapconf_c = log(lapconf_tot_tmp.val_grid_point(0,0,0,0)) ;
	double log_confo_c = log(confo_tot.val_grid_point(0,0,0,0)) ;
	double loggam_c = loggam.val_grid_point(0,0,0,0) ;
	double pot_centri_c = pot_centri.val_grid_point(0,0,0,0) ;

	ent = (ent_c + log_lapconf_c - log_confo_c + loggam_c + pot_centri_c)
	  - log(lapconf_tot_tmp) + log(confo_tot) - loggam - pot_centri ;
	ent.std_spectral_base() ;


	//----------------------------------------------------------
	// Change the enthalpy field to be set its maximum position
	// at the coordinate center
	//----------------------------------------------------------

	double dentdx = ent.dsdx().val_grid_point(0,0,0,0) ;
	double dentdy = ent.dsdy().val_grid_point(0,0,0,0) ;

	cout << "dH/dx|_center = " << dentdx << endl ;
	cout << "dH/dy|_center = " << dentdy << endl ;

	double dec_fact_x = 1. ;
	double dec_fact_y = 1. ;

	Scalar func_in(mp) ;
	func_in = 1. - dec_fact_x * (dentdx/ent_c) * mp.x
	  - dec_fact_y * (dentdy/ent_c) * mp.y ;

	func_in.annule(nzet, nz-1) ;
	func_in.std_spectral_base() ;

	Scalar func_ex(mp) ;
	func_ex = 1. ;
	func_ex.annule(0, nzet-1) ;
	func_ex.std_spectral_base() ;

	// New enthalpy field
	// ------------------
	ent = ent * (func_in + func_ex) ;

	(ent.set_spectral_va()).smooth(nzet, ent.set_spectral_va()) ;

	double dentdx_new = ent.dsdx().val_grid_point(0,0,0,0) ;
	double dentdy_new = ent.dsdy().val_grid_point(0,0,0,0) ;
	cout << "dH/dx|_new    = " << dentdx_new << endl ;
	cout << "dH/dy|_new    = " << dentdy_new << endl ;

	//-----------------------------------------------------
	// Adaptation of the mapping to the new enthalpy field
	//-----------------------------------------------------

	double dent_eq = ent.dsdr().val_point(ray_eq_pi(),M_PI/2.,M_PI) ;
	double dent_pole = ent.dsdr().val_point(ray_pole(),0.,0.) ;
	double rap_dent = fabs( dent_eq / dent_pole ) ;
	cout << "| dH/dr_eq / dH/dr_pole | = " << rap_dent << endl ;
	diff_dHdr = rap_dent ;

	if ( rap_dent < thres_adapt ) {
	    adapt_flag = 0 ;	// No adaptation of the mapping
	    cout << "******* FROZEN MAPPING  *********" << endl ;
	}
	else {
	    adapt_flag = 1 ;	// The adaptation of the mapping is to be
	                        //  performed
	}

	ent_limit.set_etat_qcq() ;
	for (int l=0; l<nzet; l++) {	// loop on domains inside the star
	    ent_limit.set(l) = ent.val_grid_point(l, k_b, j_b, i_b) ;
	}
	ent_limit.set(nzet-1) = ent_b  ;

	Map_et mp_prev = mp_et ;

	Cmp ent_cmp(ent) ;
	mp.adapt(ent_cmp, par_adapt) ;
	ent = ent_cmp ;

	// Readjustment of the external boundary of domain l=nzet
	// to keep a fixed ratio with respect to star's surface

	double rr_in_1 = mp.val_r(1, -1., M_PI/2., 0.) ;

	// Resize of the outer boundary of the shell including the BH
	double rr_out_nm2 = mp.val_r(nz-2, 1., M_PI/2., 0.) ;
	mp.resize(nz-2, rr_in_1/rr_out_nm2 * fact_resize(1)) ;

	// Resize of the inner boundary of the shell including the BH
	double rr_out_nm3 = mp.val_r(nz-3, 1., M_PI/2., 0.) ;
	mp.resize(nz-3, rr_in_1/rr_out_nm3 * fact_resize(0)) ;

	//	if (mer % 2 == 0) {

	if (nz > 4) {

	  // Resize of the domain 1
	  double rr_out_1 = mp.val_r(1, 1., M_PI/2., 0.) ;
	  mp.resize(1, rr_in_1/rr_out_1 * resize_ns) ;

	  if (nz > 5) {

	    // Resize of the domain from 2 to N-4
	    double rr_out_nm3_new = mp.val_r(nz-3, 1., M_PI/2., 0.) ;

	    for (int i=1; i<nz-4; i++) {

	      double rr_out_i = mp.val_r(i, 1., M_PI/2., 0.) ;

	      double rr_mid = rr_out_i
		+ (rr_out_nm3_new - rr_out_i) / double(nz - 3 - i) ;

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

	  }  // End of (nz > 5) loop

	}  // End of (nz > 4) loop

	//	}  // End of (mer % 2) loop

	//----------------------------------------------------
	// Computation of the enthalpy at the new grid points
	//----------------------------------------------------

	mp_prev.homothetie(alpha_r) ;

	Cmp ent_cmp2 (ent) ;
	mp.reevaluate(&mp_prev, nzet+1, ent_cmp2) ;
	ent = ent_cmp2 ;

	double ent_s_max = -1 ;
	int k_s_max = -1 ;
	int j_s_max = -1 ;
	
	for (int k=0; k<mg->get_np(l_b); k++) {
	    for (int j=0; j<mg->get_nt(l_b); j++) {
	        double xx = fabs( ent.val_grid_point(l_b, k, j, i_b) ) ;
		if (xx > ent_s_max) {
		    ent_s_max = xx ;
		    k_s_max = k ;
		    j_s_max = j ;
		}
	    }
	}
	cout << "Max. abs(enthalpy) at the boundary between domains nzet-1"
	     << " and nzet : " << endl ;
	cout << "   " << ent_s_max << " reached for k = " << k_s_max
	     << " and j = " << j_s_max << endl ;

	//-------------------
	// Equation of state
	//-------------------

	equation_of_state() ; 	// computes new values for nbar (n), ener (e)
	                        // and press (p) from the new ent (H)

	//----------------------------------------------------------
	// Matter source terms in the gravitational field equations
	//----------------------------------------------------------

	hydro_euler_bhns(kerrschild, mass_bh, sepa) ;
	                      // computes new values for ener_euler (E),
	                      // s_euler (S) and u_euler (U^i)


	//-------------------------------------------------
	// Computation of the minimum of the indicator chi
	//-------------------------------------------------
	double azimu_min = phi_min() ;
	double rad_chi_min = radius_p(azimu_min) ;
	double chi_min = chi_rp(rad_chi_min, azimu_min) ;

	cout << "| dH/dr_eq / dH/dr_pole | (minimum) = " << chi_min << endl ;
	cout << "     phi =    " << azimu_min / M_PI << " [M_PI]" << endl ;
	cout << "     radius = " << rad_chi_min / km << " [km]" << endl ;

	diff_dHdr_min = chi_min ;
	diff_phi_min = azimu_min ;
	diff_radius = rad_chi_min ;

	//-----------------------------------
	// Poisson equation for lapconf_auto
	//-----------------------------------

	// Source
	//--------

	Scalar sou_lap1(mp) ;  // dzpuis = 0
	sou_lap1 = qpig * lapconf_tot_tmp * pow(confo_tot,4.)
	  * (0.5*ener_euler + s_euler) ;

	sou_lap1.std_spectral_base() ;
	sou_lap1.annule(nzet,nz-1) ;
	sou_lap1.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	Scalar sou_lap2(mp) ;  // dzpuis = 4
	sou_lap2 =  0.875 * (lapconf_auto+0.5) * taij_quad_auto
	  / pow(confo_auto+0.5,8.) ;
	sou_lap2.std_spectral_base() ;

	source_lapconf = sou_lap1 + sou_lap2 ;

	source_lapconf.std_spectral_base() ;
	//	source_lapse.annule(nzet,nz-1) ;

	if (filter_r != 0) {
	    if (source_lapconf.get_etat() != ETATZERO) {
	        source_lapconf.filtre(filter_r) ;
		//	        source_lapse.filtre_r(filter_r,0) ;
	    }
	}

	assert(source_lapconf.get_dzpuis() == 4) ;

	// Resolution of the Poisson equation (Outer BC : lapconf_m1 -> 0)
	// ----------------------------------

	lapconf_m1.set_etat_qcq() ;
	lapconf_m1 = lapconf_auto - 0.5 ;
	source_lapconf.poisson(par_lapconf, lapconf_m1) ;
	ssjm1_lapconf = ssjm1lapconf ;

	// Check: has the Poisson equation been correctly solved ?
	// -------------------------------------------------------

	Tbl tdiff_lapconf = diffrel(lapconf_m1.laplacian(), source_lapconf) ;
	cout <<
	  "Relative error in the resolution of the equation for lapconf_auto : "
	     << endl ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_lapconf(l) << "  " ;
	}
	cout << endl ;
	diff_lapconf = max(abs(tdiff_lapconf)) ;

	// Re-construction of the lapconf function
	// ---------------------------------------
	lapconf_auto = lapconf_m1 + 0.5 ;
	          // lapconf_tot = lapconf_auto + lapconf_comp
	          // lapconf_auto, _comp -> 0.5 (r -> inf)
	          // lapconf_tot -> 1 (r -> inf)

	//---------------------------------
	// Poisson equation for confo_auto
	//---------------------------------

	// Source
	//--------

	Scalar sou_con1(mp) ;  // dzpuis = 0
	sou_con1 = - 0.5 * qpig * pow(confo_tot,5.) * ener_euler ;
	sou_con1.std_spectral_base() ;
	sou_con1.annule(nzet,nz-1) ;
	sou_con1.inc_dzpuis(4) ;  // dzpuis : 0 -> 4

	Scalar sou_con2(mp) ;  // dzpuis = 4
	sou_con2 = - 0.125 * taij_quad_auto / pow(confo_auto+0.5,7.) ;
	sou_con2.std_spectral_base() ;

	source_confo = sou_con1 + sou_con2 ;

	source_confo.std_spectral_base() ;
	//	source_confo.annule(nzet,nz-1) ;

	if (filter_r != 0) {
	    if (source_confo.get_etat() != ETATZERO) {
	        source_confo.filtre(filter_r) ;
		//	        source_confo.filtre_r(filter_r,0) ;
	    }
	}

	assert(source_confo.get_dzpuis() == 4) ;

	// Resolution of the Poisson equation (Outer BC : confo_m1 -> 0)
	// ----------------------------------

	confo_m1.set_etat_qcq() ;
	confo_m1 = confo_auto - 0.5 ;
	source_confo.poisson(par_confo, confo_m1) ;
	ssjm1_confo = ssjm1confo ;

	// Check: has the Poisson equation been correctly solved ?
	// -------------------------------------------------------

	Tbl tdiff_confo = diffrel(confo_m1.laplacian(), source_confo) ;
	cout <<
	  "Relative error in the resolution of the equation for confo_auto : "
	     << endl ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_confo(l) << "  " ;
	}
	cout << endl ;
	diff_confo = max(abs(tdiff_confo)) ;

	// Re-construction of the conformal factor
	// ---------------------------------------
	confo_auto = confo_m1 + 0.5 ; // confo_tot = confo_auto + confo_comp
                                      // confo_auto, _comp -> 0.5 (r -> inf)
	                              // confo_tot -> 1 (r -> inf)

	//----------------------------------------
	// Vector Poisson equation for shift_auto
	//----------------------------------------

	// Source
	// ------

	Vector sou_shif1(mp, CON, mp.get_bvect_cart()) ;  // dzpuis = 0
	sou_shif1.set_etat_qcq() ;

	for (int i=1; i<=3; i++) {
	    sou_shif1.set(i) = 4.*qpig * lapconf_tot_tmp
	      * pow(confo_tot, 3.)
	      * (ener_euler + press) * u_euler(i) ;
	}

	sou_shif1.std_spectral_base() ;
	sou_shif1.annule(nzet, nz-1) ;

	for (int i=1; i<=3; i++) {
	    (sou_shif1.set(i)).inc_dzpuis(4) ;  // dzpuis: 0 -> 4
	}

	Vector sou_shif2(mp, CON, mp.get_bvect_cart()) ;  // dzpuis = 4
	sou_shif2.set_etat_qcq() ;
	for (int i=1; i<=3; i++) {
	    sou_shif2.set(i) = 2. *
	      (taij_auto(i,1)*(d_lapconf_auto(1)
			       -7.*(lapconf_auto+0.5)*d_confo_auto(1)
			       /(confo_auto+0.5))
	       +taij_auto(i,2)*(d_lapconf_auto(2)
				-7.*(lapconf_auto+0.5)*d_confo_auto(2)
				/(confo_auto+0.5))
	       +taij_auto(i,3)*(d_lapconf_auto(3)
				-7.*(lapconf_auto+0.5)*d_confo_auto(3)
				/(confo_auto+0.5))
	       ) / pow(confo_auto+0.5,7.) ;
	}
	sou_shif2.std_spectral_base() ;

	source_shift = sou_shif1 + sou_shif2 ;

	source_shift.std_spectral_base() ;
	//	source_shift.annule(nzet, nz-1) ;

	// Resolution of the Poisson equation
	// ----------------------------------

	// Filter for the source of shift vector

	if (filter_r_s != 0) {
	    for (int i=1; i<=3; i++) {
	        if (source_shift(i).get_etat() != ETATZERO) {
		    source_shift.set(i).filtre(filter_r_s) ;
		    //		    source_shift.set(i).filtre_r(filter_r_s, 0) ;
		}
	    }
	}

	if (filter_p_s != 0) {
	    for (int i=1; i<=3; i++) {
	        if (source_shift(i).get_etat() != ETATZERO) {
		    (source_shift.set(i)).filtre_phi(filter_p_s, nz-1) ;
		    //		    (source_shift.set(i)).filtre_phi(filter_p_s, 0) ;
		}
	    }
	}

	for (int i=1; i<=3; i++) {
	    if(source_shift(i).dz_nonzero()) {
	        assert( source_shift(i).get_dzpuis() == 4 ) ;
	    }
	    else {
	        (source_shift.set(i)).set_dzpuis(4) ;
	    }
	}

	double lambda = double(1) / double(3) ;

	Tenseur source_p(mp, 1, CON, mp.get_bvect_cart() ) ;
	source_p.set_etat_qcq() ;
	for (int i=0; i<3; i++) {
	    source_p.set(i) = Cmp(source_shift(i+1)) ;
	}

	Tenseur vect_auxi(mp, 1, CON, mp.get_bvect_cart()) ;
	vect_auxi.set_etat_qcq() ;
	for (int i=0; i<3 ;i++) {
	    vect_auxi.set(i) = 0. ;
	}
	Tenseur scal_auxi(mp) ;
	scal_auxi.set_etat_qcq() ;
	scal_auxi.set().annule_hard() ;
	scal_auxi.set_std_base() ;

	Tenseur resu_p(mp, 1, CON, mp.get_bvect_cart() ) ;
	resu_p.set_etat_qcq() ;

	source_p.poisson_vect(lambda, par_shift2, resu_p, vect_auxi,
			      scal_auxi) ;

	for (int i=1; i<=3; i++)
	    shift_auto.set(i) = resu_p(i-1) ;

	ssjm1_khi = ssjm1khi ;

	for (int i=0; i<3; i++) {
	    ssjm1_wshift.set(i+1) = ssjm1wshift(i) ;
	}

	// Check: has the equation for shift_auto been correctly solved ?
	// --------------------------------------------------------------

	Vector lap_shift = shift_auto.derive_con(flat).divergence(flat)
	  + lambda * shift_auto.divergence(flat).derive_con(flat) ;

	source_shift.dec_dzpuis() ;

	Tbl tdiff_shift_x = diffrel(lap_shift(1), source_shift(1)) ;
	Tbl tdiff_shift_y = diffrel(lap_shift(2), source_shift(2)) ;
	Tbl tdiff_shift_z = diffrel(lap_shift(3), source_shift(3)) ;

	cout <<
	  "Relative error in the resolution of the equation for shift_auto : "
	     << endl ;
	cout << "x component : " ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_shift_x(l) << "  " ;
	}
	cout << endl ;
	cout << "y component : " ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_shift_y(l) << "  " ;
	}
	cout << endl ;
	cout << "z component : " ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_shift_z(l) << "  " ;
	}
	cout << endl ;

	diff_shift_x = max(abs(tdiff_shift_x)) ;
	diff_shift_y = max(abs(tdiff_shift_y)) ;
	diff_shift_z = max(abs(tdiff_shift_z)) ;


	//-----------------------------
	// Relative change in enthalpy
	//-----------------------------

	Tbl diff_ent_tbl = diffrel( ent, ent_jm1 ) ;
	diff_ent = diff_ent_tbl(0) ;
	for (int l=0; l<nzet; l++) {
	    diff_ent += diff_ent_tbl(l) ;
	}
	diff_ent /= nzet ;

	ent_jm1 = ent ;

	/*
	des_profile( lapconf_auto, 0., 10.,
		     M_PI/2., M_PI, "Self lapconf function of NS",
		     "Lapconf (theta=pi/2, phi=0)" ) ;

	des_profile( lapconf_tot, 0., 10.,
		     M_PI/2., M_PI, "Total lapconf function seen by NS",
		     "Lapconf (theta=pi/2, phi=0)" ) ;

	des_profile( confo_auto, 0., 10.,
		     M_PI/2., M_PI, "Self conformal factor of NS",
		     "Confo (theta=pi/2, phi=0)" ) ;

	des_profile( confo_tot, 0., 10.,
		     M_PI/2., M_PI, "Total conformal factor seen by NS",
		     "Confo (theta=pi/2, phi=0)" ) ;

	des_coupe_vect_z( shift_auto, 0., -3., 0.5, 3,
			  "Self shift vector of NS") ;

	des_coupe_vect_z( shift_tot, 0., -3., 0.5, 3,
			  "Total shift vector seen by NS") ;
	*/
    } // End of main loop

    //=========================================================
    //                    End of iteration
    //=========================================================


}
