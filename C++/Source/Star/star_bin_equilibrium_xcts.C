/*
 *   Method of class Star_bin to compute an equilibrium configuration
 *  (see file star.h for documentation).
 */

/*
 *   Copyright (c) 2010 Michal Bejger
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.14  2016/12/05 16:18:14  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.13  2014/10/13 08:53:38  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.12  2014/10/06 15:13:16  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.11  2011/03/25 16:28:12  e_gourgoulhon
 * Still in progress
 *
 * Revision 1.10  2010/12/20 15:42:10  m_bejger
 * Various rearrangements of fields in Poissson equations
 *
 * Revision 1.9  2010/12/14 17:34:42  m_bejger
 * Improved iteration for beta_auto poisson_vect()
 *
 * Revision 1.8  2010/12/09 10:48:06  m_bejger
 * Testing the main equations
 *
 * Revision 1.7  2010/10/26 18:46:28  m_bejger
 * Added table fact_resize for domain resizing
 *
 * Revision 1.6  2010/10/18 19:08:14  m_bejger
 * Changed to allow for calculations with more than one domain in the star
 *
 * Revision 1.5  2010/06/23 20:40:56  m_bejger
 * Corrections in equations for Psi_auto, chi_auto and beta_auto
 *
 * Revision 1.4  2010/06/18 13:28:59  m_bejger
 * Adjusted the computation of the first integral, radial scale
 *
 * Revision 1.3  2010/06/17 17:05:06  m_bejger
 * Testing version
 *
 * Revision 1.2  2010/06/15 08:21:21  m_bejger
 * Minor changes; still not working properly
 *
 * Revision 1.1  2010/05/04 07:51:05  m_bejger
 * Initial version
 *
 * $Header$
 *
 */

// C headers
#include <cmath>

// Lorene headers
#include "cmp.h"
#include "tenseur.h"
#include "metrique.h"
#include "star.h"
#include "param.h"
#include "graphique.h"
#include "utilitaires.h"
#include "tensor.h"
#include "nbr_spx.h"
#include "unites.h"


namespace Lorene {
void Star_bin_xcts::equilibrium(double ent_c,
								int mermax,
								int mermax_potvit,
			   					int mermax_poisson,
			   					double relax_poisson,
			   					double relax_potvit,
			   					double thres_adapt,
			   					const Tbl& fact_resize,
			   					const Tbl* pent_limit,
			   					Tbl& diff) {

    // Fundamental constants and units
    // -------------------------------

  using namespace Unites ;

    // Initializations
    // ---------------

    const Mg3d* mg = mp.get_mg() ;
    int nz = mg->get_nzone() ;	    // total number of domains

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
    double& diff_psi = diff.set(2) ;
    double& diff_chi = diff.set(3) ;
    double& diff_beta_x = diff.set(4) ;
    double& diff_beta_y = diff.set(5) ;
    double& diff_beta_z = diff.set(6) ;

    // Parameters for the function Map_et::adapt
    // -----------------------------------------

    Param par_adapt ;
    int nitermax = 100 ;
    int niter ;
    int adapt_flag = 1 ;    //  1 = performs the full computation,
    //  0 = performs only the rescaling by
    //      the factor alpha_r
    //##    int nz_search = nzet + 1 ;  // Number of domains for searching the
    // enthalpy

	int nz_search = nzet ;					// Number of domains
											// for searching
											// the enthalpy isosurfaces

	double precis_secant = 1.e-14 ;
	double alpha_r ;
	double reg_map = 1. ;					// 1 = regular mapping,
    										// 0 = contracting mapping

	par_adapt.add_int(nitermax, 0) ; 		// maximum number of iterations to
    										// locate zeros by the secant method

	par_adapt.add_int(nzet, 1) ;			// number of domains where the adjustment
											// to the isosurfaces of ent is to be performed

	par_adapt.add_int(nz_search, 2) ;		// number of domains to search
											// the enthalpy isosurface

	par_adapt.add_int(adapt_flag, 3) ;		// 1 = performs the full computation,
											// 0 = performs only the rescaling by
											// the factor alpha_r

	par_adapt.add_int(j_b, 4) ;				// theta index of the collocation point
											//  (theta_*, phi_*)

	par_adapt.add_int(k_b, 5) ; 			// theta index of the collocation point
											//  (theta_*, phi_*)

	par_adapt.add_int_mod(niter, 0) ;		// number of iterations actually used in
											//  the secant method

	par_adapt.add_double(precis_secant, 0) ;// required absolute precision in
											// the determination of zeros by
											// the secant method
	par_adapt.add_double(reg_map, 1) ;		// 1. = regular mapping,  0 = contracting mapping

	par_adapt.add_double(alpha_r, 2) ;		// factor by which all the radial
											// distances will be multiplied

    // Enthalpy values for the adaptation
    Tbl ent_limit(nzet) ;
    if (pent_limit != 0x0) ent_limit = *pent_limit ;

	par_adapt.add_tbl(ent_limit, 0) ;		// array of values of the field ent
											// to define the isosurfaces.

    double precis_poisson = 1.e-16 ;

    Cmp ssjm1psi (ssjm1_psi) ;
    Cmp ssjm1chi (ssjm1_chi) ;

    // Parameters for the function Scalar::poisson for Psi_auto
    // ---------------------------------------------------------------

    Param par_psi ;

    par_psi.add_int(mermax_poisson,  0) ;  	// maximum number of iterations
    par_psi.add_double(relax_poisson,  0) ; // relaxation parameter
    par_psi.add_double(precis_poisson, 1) ; // required precision
    par_psi.add_int_mod(niter, 0) ; 		// number of iterations actually used
    par_psi.add_cmp_mod( ssjm1psi ) ;

    // Parameters for the function Scalar::poisson for chi_auto
    // ---------------------------------------------------------------

    Param par_chi ;

    par_chi.add_int(mermax_poisson,  0) ;   // maximum number of iterations
    par_chi.add_double(relax_poisson,  0) ; // relaxation parameter
    par_chi.add_double(precis_poisson, 1) ; // required precision
    par_chi.add_int_mod(niter, 0) ; 		// number of iterations actually used
    par_chi.add_cmp_mod( ssjm1chi ) ;

    // Parameters for the function Vector::poisson for beta
    // ----------------------------------------------------

    Param par_beta ;

    par_beta.add_int(mermax_poisson,  0) ;  	// maximum number of 
    											// iterations
    par_beta.add_double(relax_poisson,  0) ; 	// relaxation parameter
    par_beta.add_double(precis_poisson, 1) ; 	// required precision
    par_beta.add_int_mod(niter, 0) ; 			// number of iterations 
    											// actually used

    // Sources at the previous step, for a poisson_vect() solver 
    Cmp ssjm1khi (ssjm1_khi) ;

    Tenseur ssjm1wbeta(mp, 1, CON, mp.get_bvect_cart()) ;
    ssjm1wbeta.set_etat_qcq() ;
    for (int i=0; i<3; i++) ssjm1wbeta.set(i) = Cmp(ssjm1_wbeta(i+1)) ;

    par_beta.add_cmp_mod(ssjm1khi) ;
    par_beta.add_tenseur_mod(ssjm1wbeta) ;

    // Redefinition of external potential
    // See Eq (99) from Gourgoulhon et al. (2001)
    // logN = logN_auto + logn_ac_rest = log(chi_auto + 1.)
    // - log(Psi_auto + 1.) + logn_ac_rest
    //------------------------------------

    Scalar Psi_auto_p1 = Psi_auto + 1. ;
    Scalar chi_auto_p1 = chi_auto + 1. ;

    Scalar logn_auto = log(chi_auto_p1) - log(Psi_auto_p1) ;
    logn_auto.std_spectral_base() ;

    Scalar logn_ac_rest = log(1. + chi_comp/chi_auto_p1)
                        - log(1. + Psi_comp/Psi_auto_p1) ;

    logn_ac_rest.std_spectral_base() ;

    //cout << "logn_auto" << norme(logn_auto) << endl ;
    //cout << "logn_ac_rest" << norme(logn_ac_rest) << endl ;
    //cout << "pot_centri" << norme(pot_centri) << endl ;
    //cout << "loggam" << norme(loggam) << endl ;

    Scalar pot_ext = logn_ac_rest + pot_centri + loggam ;
    //cout << "pot_ext" << norme(pot_ext) << endl ;

    Scalar ent_jm1 = ent ;	// Enthalpy at previous step

    Scalar source_tot(mp) ; // source term in equations Psi_auto
							// and chi_auto

    Vector source_beta(mp, CON, mp.get_bvect_cart()) ; // source term
                                                       // in the equation
                                                       // for beta_auto

    //==================================================================
    // 			Start of iteration
    //==================================================================

    for(int mer=0 ; mer<mermax ; mer++ ) {

	cout << "-----------------------------------------------" << endl ;
	cout << "step: " << mer << endl ;
	cout << "diff_ent = " << diff_ent << endl ;

	//-----------------------------------------------------
	// Resolution of the elliptic equation for the velocity
	// scalar potential
	//-----------------------------------------------------

	if (irrotational) {

	    diff_vel_pot = velocity_potential(mermax_potvit, precis_poisson,
					      relax_potvit) ;
	}

	diff_vel_pot = 0. ; // to avoid the warning

	//-----------------------------------------------------
	// Computation of the new radial scale
	//--------------------------------------------------

	// alpha_r (r = alpha_r r') is determined so that the enthalpy
	// takes the requested value ent_b at the stellar surface

	// Values at the center of the star:
	double logn_auto_c  = logn_auto.val_grid_point(0, 0, 0, 0) ;
	double pot_ext_c  = pot_ext.val_grid_point(0, 0, 0, 0) ;

	// Search for the reference point (theta_*, phi_*) [notation of
	//  Bonazzola, Gourgoulhon & Marck PRD 58, 104020 (1998)]
	//  at the surface of the star
	// ------------------------------------------------------------
	double alpha_r2 = 0 ;
	for (int k=0; k<mg->get_np(l_b); k++) {
	    for (int j=0; j<mg->get_nt(l_b); j++) {

		double pot_ext_b  = pot_ext.val_grid_point(l_b, k, j, i_b) ;
		double logn_auto_b  = logn_auto.val_grid_point(l_b, k, j, i_b) ;
		// See Eq (100) from Gourgoulhon et al. (2001)
		double alpha_r2_jk = ( ent_c - ent_b + pot_ext_c - pot_ext_b) /
		    ( logn_auto_b - logn_auto_c ) ;

		// cout << "k, j, alpha_r2_jk : " << k << " " << j << " " << alpha_r2_jk << endl ; 
		    
		if (alpha_r2_jk > alpha_r2) {
		    alpha_r2 = alpha_r2_jk ;
		    k_b = k ;
		    j_b = j ;
		}

	    }
	}
  
	alpha_r = sqrt(alpha_r2) ;

	cout << "k_b, j_b, alpha_r: " << k_b << "  " << j_b << "  "
	     <<  alpha_r << endl ;

	// New value of logn_auto
	// ----------------------

    Psi_auto = pow(Psi_auto +1.,alpha_r2) - 1. ;
    chi_auto = pow(chi_auto +1.,alpha_r2) - 1. ;
    Psi_auto.std_spectral_base() ;
    chi_auto.std_spectral_base() ;

	//------------------------------------------------------------
	// Change the values of the inner points of the domain adjascent
	// to the surface of the star by those of the outer points of
	// the domain under the surface
	//------------------------------------------------------------

    Psi_auto.set_spectral_va().smooth(nzet, Psi_auto.set_spectral_va()) ;
    chi_auto.set_spectral_va().smooth(nzet, chi_auto.set_spectral_va()) ;

    logn_auto = log(chi_auto + 1.) - log(Psi_auto + 1.) ;
    logn_auto.std_spectral_base() ;

	logn_auto_c  = logn_auto.val_grid_point(0, 0, 0, 0) ;

	//------------------------------------------
	// First integral	-->  enthalpy in all space
	// See Eq (98) from Gourgoulhon et al. (2001)
	//-------------------------------------------

	ent = (ent_c + logn_auto_c + pot_ext_c) - logn_auto - pot_ext ;
	cout.precision(8) ;
	cout << "pot" << norme(pot_ext) << endl ;

	(ent.set_spectral_va()).smooth(nzet, ent.set_spectral_va()) ;

	//----------------------------------------------------
	// Adaptation of the mapping to the new enthalpy field
	//----------------------------------------------------

	// Shall the adaptation be performed (cusp) ?
	// ------------------------------------------

	double dent_eq = ent.dsdr().val_point(ray_eq(),M_PI/2.,0.) ;
	double dent_pole = ent.dsdr().val_point(ray_pole(),0.,0.) ;
	double rap_dent = fabs( dent_eq / dent_pole ) ;
	cout << "| dH/dr_eq / dH/dr_pole | = " << rap_dent << endl ;

	if ( rap_dent < thres_adapt ) {
	    adapt_flag = 0 ;	// No adaptation of the mapping
	    cout << "******* FROZEN MAPPING  *********" << endl ;
	}
	else{
	    adapt_flag = 1 ;	// The adaptation of the mapping is to be
	    //  performed
	}

	ent_limit.set_etat_qcq() ;
	for (int l=0; l<nzet; l++) {	// loop on domains inside the star
	    ent_limit.set(l) = ent.val_grid_point(l, k_b, j_b, i_b) ;

	} ent_limit.set(nzet-1) = ent_b  ;

	Map_et mp_prev = mp_et ;

	Cmp ent_cmp(ent) ;
	mp.adapt(ent_cmp, par_adapt) ;
	ent = ent_cmp ;

	// Readjustment of the external boundary of domain l=nzet
	// to keep a fixed ratio with respect to star's surface

	double rr_in_1 = mp.val_r(nzet, -1., M_PI/2., 0.) ;

	// Resizes the outer boundary of the shell including the comp. NS
	double rr_out_nm2 = mp.val_r(nz-2, 1., M_PI/2., 0.) ;

	mp.resize(nz-2, rr_in_1/rr_out_nm2 * fact_resize(1)) ;

	// Resizes the inner boundary of the shell including the comp. NS
	double rr_out_nm3 = mp.val_r(nz-3, 1., M_PI/2., 0.) ;

	mp.resize(nz-3, rr_in_1/rr_out_nm3 * fact_resize(0)) ;

	if (nz > nzet+3) {

	    // Resize of the domain from 1(nzet) to N-4
	    double rr_out_nm3_new = mp.val_r(nz-3, 1., M_PI/2., 0.) ;

	    for (int i=nzet-1; i<nz-4; i++) {

	        double rr_out_i = mp.val_r(i, 1., M_PI/2., 0.) ;

		double rr_mid = rr_out_i
		  + (rr_out_nm3_new - rr_out_i) / double(nz - 3 - i) ;

		double rr_2timesi = 2. * rr_out_i ;

		if (rr_2timesi < rr_mid) {

		    double rr_out_ip1 = mp.val_r(i+1, 1., M_PI/2., 0.) ;

		    mp.resize(i+1, rr_2timesi / rr_out_ip1) ;

		} else {

		    double rr_out_ip1 = mp.val_r(i+1, 1., M_PI/2., 0.) ;

		    mp.resize(i+1, rr_mid / rr_out_ip1) ;

		}  // End of else

	    }  // End of i loop

	}  // End of (nz > nzet+3) loop

//	if (nz>= 5) {
//
//	  double separation = 2. * fabs(mp.get_ori_x()) ;
//	  double ray_eqq = ray_eq() ;
//	  double ray_eqq_pi = ray_eq_pi() ;
//	  double new_rr_out_2 = (separation - ray_eqq) * 0.95 ;
//	  double new_rr_out_3 = (separation + ray_eqq_pi) * 1.05 ;
//
//	  double rr_in_1 = mp.val_r(1,-1., M_PI/2, 0.) ;
//	  double rr_out_1 = mp.val_r(1, 1., M_PI/2, 0.) ;
//	  double rr_out_2 = mp.val_r(2, 1., M_PI/2, 0.) ;
//	  double rr_out_3 = mp.val_r(3, 1., M_PI/2, 0.) ;
//
//	  mp.resize(1, 0.5*(new_rr_out_2 + rr_in_1) / rr_out_1) ;
//	  mp.resize(2, new_rr_out_2 / rr_out_2) ;
//	  mp.resize(3, new_rr_out_3 / rr_out_3) ;
//
//	  for (int dd=4; dd<=nz-2; dd++) {
//	    mp.resize(dd, new_rr_out_3 * pow(4., dd-3) /
//		      mp.val_r(dd, 1., M_PI/2, 0.)) ;
//	  }
//
//	}

	//else cout << "too small number of domains" << endl ;


	//------------------------------------------------------------------
	// Computation of the enthalpy at the new grid points
	//------------------------------------------------------------------

	mp_prev.homothetie(alpha_r) ;

	//Cmp ent_cmp2 (ent) ;
	//mp.reevaluate_symy(&mp_prev, nzet+1, ent_cmp2) ;
	//ent = ent_cmp2 ;

    mp.reevaluate_symy(&mp_prev, nzet+1, ent) ;

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
	cout << "   " << ent_s_max << " reached for k = " << k_s_max <<
	    " and j = " << j_s_max << endl ;

	//------------------------------------------------------------------
	// Equation of state
	//------------------------------------------------------------------

	equation_of_state() ; 	// computes new values for nbar (n), ener (e)
							// and press (p) from the new ent (H)

	//------------------------------------------------------------------
	// Matter source terms in the gravitational field equations
	//------------------------------------------------------------------

	hydro_euler() ;		// computes new values for ener_euler (E),
			  			// s_euler (S) and u_euler (U^i)

	// -------------------------------
	// AUXILIARY QUANTITIES
	// -------------------------------

	int nr = mp.get_mg()->get_nr(0) ;
	int nt = mp.get_mg()->get_nt(0) ;
	int np = mp.get_mg()->get_np(0) ;

    Scalar Psi3 = psi4 / Psi ;
    Psi3.std_spectral_base() ;

    Scalar sPsi7 = Psi * pow(psi4, -2.) ;
    sPsi7.std_spectral_base() ;

	//------------------------------------------------------------------
	// Poisson equation for Psi_auto (Eq. 8.127 of arXiv:gr-qc/0703035)
	//------------------------------------------------------------------

	// Source
	//--------

	source_tot = - 0.5 * qpig * psi4 % Psi % ener_euler
                 - 0.125 * sPsi7 * (hacar_auto + hacar_comp) ;
	source_tot.std_spectral_base() ;

	// Resolution of the Poisson equation
	// ----------------------------------

    cout << "Resolution of the Poisson equation for Psi_auto:" << endl ;
	source_tot.poisson(par_psi, Psi_auto) ;
	ssjm1_psi = ssjm1psi ;

	cout << "Psi_auto: " << endl << norme(Psi_auto/(nr*nt*np)) << endl ;

	// Check: has the Poisson equation been correctly solved ?
	// -----------------------------------------------------

	Tbl tdiff_psi = diffrel(Psi_auto.laplacian(), source_tot) ;
	cout <<
	    "Relative error in the resolution of the equation for Psi_auto : "
	     << endl ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_psi(l) << "  " ;
	}
	cout << endl
         << "==========================================================="
         << endl ;

	diff_psi = max(abs(tdiff_psi)) ;

    //------------------------------------------------------------------
	// Poisson equation for chi_auto (Eq. 8.129 of arXiv:gr-qc/0703035)
	//------------------------------------------------------------------

	// Source
	//--------

	source_tot = chi * 0.5 * qpig * psi4 % (ener_euler + 2.*s_euler) 
	           + 0.875 * nn % sPsi7 * (hacar_auto + hacar_comp)  ;
    source_tot.std_spectral_base() ;

	// Resolution of the Poisson equation
	// ----------------------------------

    cout << "Resolution of the Poisson equation for chi_auto:" << endl ;
	source_tot.poisson(par_chi, chi_auto) ;
	ssjm1_chi = ssjm1chi ;

	cout << "chi_auto: " << endl << norme(chi_auto/(nr*nt*np)) << endl ;

	// Check: has the Poisson equation been correctly solved ?
	// -----------------------------------------------------

	Tbl tdiff_chi = diffrel(chi_auto.laplacian(), source_tot) ;
	cout <<
	    "Relative error in the resolution of the equation for chi_auto : "
	     << endl ;

	for (int l=0; l<nz; l++) {
	    cout << tdiff_chi(l) << "  " ;
	}
    cout << endl
         << "==========================================================="
         << endl ;

	diff_chi = max(abs(tdiff_chi)) ;

	//------------------------------------------------------------------
	// Vector Poisson equation for beta_auto
	// (Eq. 8.128 of arXiv:gr-qc/0703035)
	//------------------------------------------------------------------

	// Source
	//--------

    source_beta = 4.* qpig * chi % Psi3 * (ener_euler + press) * u_euler 
	            + 2. * sPsi7 * contract(haij_auto, 1, dcov_chi, 0)
	            -14. * nn % sPsi7 * contract(haij_auto, 1, dcov_Psi, 0) ;
	source_beta.std_spectral_base() ;

    // Resolution of the Poisson equation
	// ----------------------------------

   	Tenseur source_p(mp, 1, CON, mp.get_bvect_cart() ) ;
	source_p.set_etat_qcq() ;
	for (int i=0; i<3; i++) {

	    source_p.set(i) = Cmp(source_beta(i+1)) ;
        //source_p.set(i).filtre(4) ;

	}

	Tenseur vect_auxi (mp, 1, CON, mp.get_bvect_cart()) ;
	vect_auxi.set_etat_qcq() ;
	for (int i=0; i<3 ;i++) vect_auxi.set(i) = Cmp(w_beta(i+1)) ;
	vect_auxi.set_std_base() ;


	Tenseur scal_auxi (mp) ;
	scal_auxi.set_etat_qcq() ;
	scal_auxi.set() = Cmp(khi) ; 
	scal_auxi.set_std_base() ;

	Tenseur resu_p(mp, 1, CON, mp.get_bvect_cart()) ;
	resu_p.set_etat_qcq() ;
	for (int i=0; i<3 ;i++) resu_p.set(i) = Cmp(beta_auto.set(i+1));
	resu_p.set_std_base() ;

   	// Resolution of the vector Poisson equation
	// -----------------------------------------

   	double lambda = double(1) / double(3) ;
    source_p.poisson_vect(lambda, par_beta, resu_p, vect_auxi, scal_auxi) ;
    //source_p.poisson_vect_oohara(lambda, par_beta, resu_p, scal_auxi) ;
  
	for (int i=1; i<=3; i++) { 	
		beta_auto.set(i) = resu_p(i-1) ;
		w_beta.set(i)    = vect_auxi(i-1) ; 
	}
	khi = scal_auxi() ;
	
    // Values of sources from the previous step 
    ssjm1_khi = ssjm1khi ;
	for (int i=0; i<3; i++) ssjm1_wbeta.set(i+1) = ssjm1wbeta(i) ;
	
	cout << "beta_auto_x" << endl << norme(beta_auto(1)/(nr*nt*np))
	     << endl ;
	cout << "beta_auto_y" << endl << norme(beta_auto(2)/(nr*nt*np))
	     << endl ;
	cout << "beta_auto_z" << endl << norme(beta_auto(3)/(nr*nt*np))
	     << endl ;

	// Check: has the equation for beta_auto been correctly solved ?
	// --------------------------------------------------------------

	Vector lap_beta = (beta_auto.derive_con(flat)).divergence(flat)
	    + lambda* beta_auto.divergence(flat).derive_con(flat) ;

	source_beta.dec_dzpuis() ;
	Tbl tdiff_beta_x = diffrel(lap_beta(1), source_beta(1)) ;
	Tbl tdiff_beta_y = diffrel(lap_beta(2), source_beta(2)) ;
	Tbl tdiff_beta_z = diffrel(lap_beta(3), source_beta(3)) ;

	cout <<
	    "Relative error in the resolution of the equation for beta_auto : "
	     << endl ;
	cout << "x component : " ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_beta_x(l) << "  " ;
	}
	cout << endl ;
	cout << "y component : " ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_beta_y(l) << "  " ;
	}
	cout << endl ;
	cout << "z component : " ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_beta_z(l) << "  " ;
	}
	cout << endl ;

	diff_beta_x = max(abs(tdiff_beta_x)) ;
	diff_beta_y = max(abs(tdiff_beta_y)) ;
	diff_beta_z = max(abs(tdiff_beta_z)) ;

	// End of relativistic equations

	//-------------------------------------------------
	//  Relative change in enthalpy
	//-------------------------------------------------

	Tbl diff_ent_tbl = diffrel( ent, ent_jm1 ) ;
	diff_ent = diff_ent_tbl(0) ;
	for (int l=1; l<nzet; l++) diff_ent += diff_ent_tbl(l) ;

	diff_ent /= nzet ;


	ent_jm1 = ent ;


  } // End of main loop

    //==================================================================
    // 			End of iteration
    //==================================================================

}







}
