
/*
 *   Method of class Star_bin to compute an equilibrium configuration
 *
 *  (see file star.h for documentation).
 */
/*
 *   Copyright (c) 2004 Francois Limousin
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

char star_bin_equilibrium_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.23  2005/09/14 12:48:02  f_limousin
 * Comment graphical outputs.
 *
 * Revision 1.22  2005/09/14 12:30:52  f_limousin
 * Saving of fields lnq and logn in class Star.
 *
 * Revision 1.21  2005/09/13 19:38:31  f_limousin
 * Reintroduction of the resolution of the equations in cartesian coordinates.
 *
 * Revision 1.20  2005/04/08 12:36:44  f_limousin
 * Just to avoid warnings...
 *
 * Revision 1.19  2005/02/24 16:27:21  f_limousin
 * Add mermax_poisson and relax_poisson in the parameters of the function.
 *
 * Revision 1.18  2005/02/24 16:04:13  f_limousin
 * Change the name of some variables (for instance dcov_logn --> dlogn).
 * Improve the resolution of the tensorial poisson equation for hh.
 *
 * Revision 1.17  2005/02/18 13:14:18  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
 * Revision 1.16  2005/02/17 17:32:53  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.15  2005/02/11 18:13:47  f_limousin
 * Important modification : all the poisson equations for the metric
 * quantities are now solved on an affine mapping.
 *
 * Revision 1.14  2004/12/17 16:23:19  f_limousin
 * Modif. comments.
 *
 * Revision 1.13  2004/06/22 12:49:12  f_limousin
 * Change qq, qq_auto and qq_comp to beta, beta_auto and beta_comp.
 *
 * Revision 1.12  2004/05/27 12:41:00  p_grandclement
 * correction of some shadowed variables
 *
 * Revision 1.11  2004/05/25 14:18:00  f_limousin
 * Include filters
 *
 * Revision 1.10  2004/05/10 10:26:22  f_limousin
 * Minor changes to avoid warnings in the compilation of Lorene
 *
 * Revision 1.9  2004/04/08 16:32:48  f_limousin
 * The new variable is ln(Q) instead of Q=psi^2*N. It improves the
 * convergence of the code.
 *
 * Revision 1.8  2004/03/25 10:29:26  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.7  2004/03/23 09:56:09  f_limousin
 * Many minor changes
 *
 * Revision 1.6  2004/02/27 21:16:32  e_gourgoulhon
 * Function contract_desal replaced by function contract with
 * argument desaliasing set to true.
 *
 * Revision 1.5  2004/02/27 09:51:51  f_limousin
 * Many minor changes.
 *
 * Revision 1.4  2004/02/21 17:05:13  e_gourgoulhon
 * Method Scalar::point renamed Scalar::val_grid_point.
 * Method Scalar::set_point renamed Scalar::set_grid_point.
 *
 * Revision 1.3  2004/01/20 15:17:48  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

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


void Star_bin::equilibrium(double ent_c, int mermax, int mermax_potvit, 
			   int mermax_poisson, double relax_poisson, 
			   double relax_potvit, double thres_adapt,
			   const Tbl& fact_resize, Tbl& diff, double om) {

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
    double& diff_logn = diff.set(2) ; 
    double& diff_lnq = diff.set(3) ; 
    double& diff_beta_x = diff.set(4) ; 
    double& diff_beta_y = diff.set(5) ; 
    double& diff_beta_z = diff.set(6) ; 
    double& diff_h11 = diff.set(7) ; 
    double& diff_h21 = diff.set(8) ; 
    double& diff_h31 = diff.set(9) ; 
    double& diff_h22 = diff.set(10) ; 
    double& diff_h32 = diff.set(11) ; 
    double& diff_h33 = diff.set(12) ; 



    // Parameters for te function Map_et::adapt
    // -----------------------------------------

    Param par_adapt ;
    int nitermax = 100 ;
    int niter ; 
    int adapt_flag = 1 ;    //  1 = performs the full computation, 
    //  0 = performs only the rescaling by 
    //      the factor alpha_r
    //##    int nz_search = nzet + 1 ;  // Number of domains for searching the 
    // enthalpy
    int nz_search = nzet ;	// Number of domains for searching the enthalpy
				//  isosurfaces

    double precis_secant = 1.e-14 ; 
    double alpha_r ; 
    double reg_map = 1. ; // 1 = regular mapping, 0 = contracting mapping

    Tbl ent_limit(nz) ; 


    par_adapt.add_int(nitermax, 0) ; // maximum number of iterations to 
    // locate zeros by the secant method
    par_adapt.add_int(nzet, 1) ;    // number of domains where the adjustment 
    // to the isosurfaces of ent is to be 
    // performed
    par_adapt.add_int(nz_search, 2) ;	// number of domains to search 	
                                        // the enthalpy isosurface
    par_adapt.add_int(adapt_flag, 3) ; //  1 = performs the full computation, 
    //  0 = performs only the rescaling by 
    //      the factor alpha_r
    par_adapt.add_int(j_b, 4) ; //  theta index of the collocation point 
    //  (theta_*, phi_*)
    par_adapt.add_int(k_b, 5) ; //  theta index of the collocation point 
    //  (theta_*, phi_*)

    par_adapt.add_int_mod(niter, 0) ;  // number of iterations actually used in
    //  the secant method
    
    par_adapt.add_double(precis_secant, 0) ; // required absolute precision in 
    // the determination of zeros by 
    // the secant method
    par_adapt.add_double(reg_map, 1)	;  // 1. = regular mapping, 
                                           // 0 = contracting mapping
    
    par_adapt.add_double(alpha_r, 2) ;	    // factor by which all the radial 
					    // distances will be multiplied 
    	   
    par_adapt.add_tbl(ent_limit, 0) ;	// array of values of the field ent 
				        // to define the isosurfaces. 
 

    Cmp ssjm1logn (ssjm1_logn) ;
    Cmp ssjm1lnq (ssjm1_lnq) ;
    Cmp ssjm1h11 (ssjm1_h11) ;
    Cmp ssjm1h21 (ssjm1_h21) ;
    Cmp ssjm1h31 (ssjm1_h31) ;
    Cmp ssjm1h22 (ssjm1_h22) ;
    Cmp ssjm1h32 (ssjm1_h32) ;
    Cmp ssjm1h33 (ssjm1_h33) ;
 

    ssjm1logn.set_etat_qcq() ;
    ssjm1lnq.set_etat_qcq() ;
    ssjm1h11.set_etat_qcq() ;
    ssjm1h21.set_etat_qcq() ;
    ssjm1h31.set_etat_qcq() ;
    ssjm1h22.set_etat_qcq() ;
    ssjm1h32.set_etat_qcq() ;
    ssjm1h33.set_etat_qcq() ;
   
    
    double precis_poisson = 1.e-16 ;     

    // Parameters for the function Scalar::poisson for logn_auto
    // ---------------------------------------------------------------
 
    Param par_logn ;    

    par_logn.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_logn.add_double(relax_poisson,  0) ; // relaxation parameter
    par_logn.add_double(precis_poisson, 1) ; // required precision
    par_logn.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_logn.add_cmp_mod( ssjm1logn ) ; 
    
    // Parameters for the function Scalar::poisson for lnq_auto
    // ---------------------------------------------------------------
    
    Param par_lnq ; 
    
    par_lnq.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_lnq.add_double(relax_poisson,  0) ; // relaxation parameter
    par_lnq.add_double(precis_poisson, 1) ; // required precision
    par_lnq.add_int_mod(niter, 0) ; // number of iterations actually used -
    par_lnq.add_cmp_mod( ssjm1lnq ) ; 
 
    // Parameters for the function Vector::poisson for beta method 2 
    // ---------------------------------------------------------------
    
    Param par_beta2 ; 
    
    par_beta2.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_beta2.add_double(relax_poisson,  0) ; // relaxation parameter
    par_beta2.add_double(precis_poisson, 1) ; // required precision
    par_beta2.add_int_mod(niter, 0) ; // number of iterations actually used 
  
    Cmp ssjm1khi (ssjm1_khi) ;
    Tenseur ssjm1wbeta(mp, 1, CON, mp.get_bvect_cart()) ;
    ssjm1wbeta.set_etat_qcq() ;
    for (int i=0; i<3; i++) {
      ssjm1wbeta.set(i) = Cmp(ssjm1_wbeta(i+1)) ;
    }
    
    par_beta2.add_cmp_mod(ssjm1khi) ; 
    par_beta2.add_tenseur_mod(ssjm1wbeta) ; 
    
    // Parameters for the function Scalar::poisson for h11_auto
    // -------------------------------------------------------------
    
    Param par_h11 ; 
    
    par_h11.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_h11.add_double(relax_poisson,  0) ; // relaxation parameter
    par_h11.add_double(precis_poisson, 1) ; // required precision
    par_h11.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_h11.add_cmp_mod( ssjm1h11 ) ; 
    
    // Parameters for the function Scalar::poisson for h21_auto
    // -------------------------------------------------------------
    
    Param par_h21 ; 
    
    par_h21.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_h21.add_double(relax_poisson,  0) ; // relaxation parameter
    par_h21.add_double(precis_poisson, 1) ; // required precision
    par_h21.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_h21.add_cmp_mod( ssjm1h21 ) ; 
    
    // Parameters for the function Scalar::poisson for h31_auto
    // -------------------------------------------------------------
    
    Param par_h31 ; 
    
    par_h31.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_h31.add_double(relax_poisson,  0) ; // relaxation parameter
    par_h31.add_double(precis_poisson, 1) ; // required precision
    par_h31.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_h31.add_cmp_mod( ssjm1h31 ) ; 
    
    // Parameters for the function Scalar::poisson for h22_auto
    // -------------------------------------------------------------
    
    Param par_h22 ; 
    
    par_h22.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_h22.add_double(relax_poisson,  0) ; // relaxation parameter
    par_h22.add_double(precis_poisson, 1) ; // required precision
    par_h22.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_h22.add_cmp_mod( ssjm1h22 ) ; 
    
    // Parameters for the function Scalar::poisson for h32_auto
    // -------------------------------------------------------------
    
    Param par_h32 ; 
    
    par_h32.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_h32.add_double(relax_poisson,  0) ; // relaxation parameter
    par_h32.add_double(precis_poisson, 1) ; // required precision
    par_h32.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_h32.add_cmp_mod( ssjm1h32 ) ; 
    
    // Parameters for the function Scalar::poisson for h33_auto
    // -------------------------------------------------------------
    
    Param par_h33 ; 
    
    par_h33.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_h33.add_double(relax_poisson,  0) ; // relaxation parameter
    par_h33.add_double(precis_poisson, 1) ; // required precision
    par_h33.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_h33.add_cmp_mod( ssjm1h33 ) ; 
    

    // External potential
    // See Eq (99) from Gourgoulhon et al. (2001)
    // ------------------
    
    cout << "logn_comp" << norme(logn_comp) << endl ;
    cout << "pot_centri" << norme(pot_centri) << endl ;
    cout << "loggam" << norme(loggam) << endl ;
    Scalar pot_ext = logn_comp + pot_centri + loggam ;
    
    Scalar ent_jm1 = ent ;	// Enthalpy at previous step
    
    Scalar source_tot(mp) ; // source term in the equation for hij_auto, 
                            // logn_auto and beta_auto
			    
    Vector source_beta(mp, CON, mp.get_bvect_cart()) ;  // source term 
    // in the equation for beta_auto



    //=========================================================================
    // 			Start of iteration
    //=========================================================================

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

	logn_auto = alpha_r2 * logn_auto ;
	logn_auto_c  = logn_auto.val_grid_point(0, 0, 0, 0) ;

	//------------------------------------------------------------
	// Change the values of the inner points of the second domain
	// by those of the outer points of the first domain
	//------------------------------------------------------------

	logn_auto.set_spectral_va().smooth(nzet, logn_auto.set_spectral_va()) ;

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
	}
	ent_limit.set(nzet-1) = ent_b  ; 

	Map_et mp_prev = mp_et ; 

	Cmp ent_cmp(ent) ;
	mp.adapt(ent_cmp, par_adapt) ; 
	ent = ent_cmp ;

	// Readjustment of the external boundary of domain l=nzet
	// to keep a fixed ratio with respect to star's surface
		
    	if (nz == 4 && nzet == 1) {
	double rr_in_1 = mp.val_r(1,-1., M_PI/2, 0.) ; 
	double rr_out_1 = mp.val_r(1, 1., M_PI/2, 0.) ; 
	double rr_out_2 = mp.val_r(2, 1., M_PI/2, 0.) ; 

	mp.resize(1, rr_in_1/rr_out_1 * fact_resize(0)) ; 
	mp.resize(2, rr_in_1/rr_out_2 * fact_resize(1)) ; 
	}
	else{

	    if (nz > 4 && nzet == 1) {
		double rr_in_1 = mp.val_r(1,-1., M_PI/2, 0.) ; 
		double rr_out_1 = mp.val_r(1, 1., M_PI/2, 0.) ; 
		double rr_out_2 = mp.val_r(2, 1., M_PI/2, 0.) ; 
		
		double fact_resize_0 ;
		if (fact_resize(0) > 2.4) fact_resize_0 = fact_resize(0)/2. ;
		else fact_resize_0 = fact_resize(0)/2. + 0.5 ;


		mp.resize(1, rr_in_1/rr_out_1 * fact_resize_0) ; 
		mp.resize(2, rr_in_1/rr_out_2 * fact_resize(0)) ; 
		for (int nnn=0; nnn<nz-4; nnn++){
		  double rr_out_nn = mp.val_r(3+nnn, 1., M_PI/2, 0.) ; 
		  mp.resize(3+nnn, (rr_in_1/rr_out_nn * fact_resize(1))
			    *pow(1.6,nnn)) ; 
		}
	    }
	    else{		
		int n_resize ;
		//      	if (nz > 4) {
		//       	  n_resize = nz - 4 ;
		if (nz > 4) {
		    n_resize = nz - 3 ;
		}
		else {
		    n_resize = nzet ;
		}
		
		double rr_in = mp.val_r(nzet,-1., M_PI/2, 0.) ; 
		double rr_out = mp.val_r(n_resize,1., M_PI/2, 0.) ; 
		
		mp.resize(n_resize, rr_in/rr_out * fact_resize(0)) ; 
	    }
	}

	//----------------------------------------------------
	// Computation of the enthalpy at the new grid points
	//----------------------------------------------------
	
	mp_prev.homothetie(alpha_r) ; 
	
	Cmp ent_cmp2 (ent) ;
	mp.reevaluate_symy(&mp_prev, nzet+1, ent_cmp2) ; 
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
	cout << "   " << ent_s_max << " reached for k = " << k_s_max <<
	    " and j = " << j_s_max << endl ; 

	//----------------------------------------------------
	// Equation of state  
	//----------------------------------------------------
	
	equation_of_state() ; 	// computes new values for nbar (n), ener (e) 
				// and press (p) from the new ent (H)
	
	//---------------------------------------------------------
	// Matter source terms in the gravitational field equations	
	//---------------------------------------------------------

	hydro_euler() ;		// computes new values for ener_euler (E), 
				// s_euler (S) and u_euler (U^i)

	
	// -------------------------------
	// AUXILIARY QUANTITIES
	// -------------------------------

	// Function exp(-(r-r_0)^2/sigma^2)
	// --------------------------------
	
	double r0 = mp.val_r(nz-2, 1, 0, 0) ;
	double sigma = r0 ;
	  
	Scalar rr (mp) ;
	rr = mp.r ;

	Scalar ff (mp) ;
	ff = exp( -(rr - r0)*(rr - r0)/sigma/sigma ) ;
	for (int ii=0; ii<nz-1; ii++)
	  ff.set_domain(ii) = 1. ;
	ff.set_outer_boundary(nz-1, 0) ;
	ff.std_spectral_base() ;
	
	//	ff.annule_domain(nz-1) ;
	//des_profile(ff, 0, 20, 0, 0) ;

	// Construction of Omega d/dphi
	// ----------------------------
	
	Vector omdsdp (mp, CON, mp.get_bvect_cart()) ;
	Scalar yya (mp) ;
	yya = mp.ya ;
	Scalar xxa (mp) ;
	xxa = mp.xa ;
	
	if (fabs(mp.get_rot_phi()) < 1e-10){ 
	  omdsdp.set(1) = - om * yya * ff ;
	  omdsdp.set(2) = om * xxa *ff ;
	  omdsdp.set(3).annule_hard() ;
	}
	else{
	  omdsdp.set(1) = om * yya * ff ;
	  omdsdp.set(2) = - om * xxa * ff ;
	  omdsdp.set(3).annule_hard() ;
	}

	omdsdp.set(1).set_outer_boundary(nz-1, 0) ;
	omdsdp.set(2).set_outer_boundary(nz-1, 0) ;

	omdsdp.set(1).set_spectral_va()
	  .set_base(*(mp.get_mg()->std_base_vect_cart()[0])) ;
	omdsdp.set(2).set_spectral_va()
	  .set_base(*(mp.get_mg()->std_base_vect_cart()[1])) ;
	omdsdp.set(3).set_spectral_va()
	  .set_base(*(mp.get_mg()->std_base_vect_cart()[2])) ;
	
	omdsdp.annule_domain(nz-1) ;

	//des_profile(omdsdp(1), 0, 20, 1, 1) ;
	//des_profile(omdsdp(2), 0, 20, 1, 1) ;
	

	// Derivatives of N and logN
	//--------------------------

	const Vector dcov_logn_auto = logn_auto.derive_cov(flat) ;
	const Vector dcon_logn_auto = logn_auto.derive_con(flat) ;
	
	Tensor dcovdcov_logn_auto = (logn_auto.derive_cov(flat))
	                                  .derive_cov(flat) ;
	dcovdcov_logn_auto.inc_dzpuis() ;

	// Derivatives of lnq, phi and Q
	//-------------------------------

	const Scalar phi (0.5 * (lnq - logn)) ;
	const Scalar phi_auto (0.5 * (lnq_auto - logn_auto)) ;

	const Vector dcov_phi_auto = phi_auto.derive_cov(flat) ;
	
	const Vector& dcov_lnq = 2*dcov_phi + dcov_logn ;
	const Vector& dcon_lnq = 2*dcon_phi + dcon_logn ;
	const Vector dcov_lnq_auto = lnq_auto.derive_cov(flat) ;
     	Tensor dcovdcov_lnq_auto = dcov_lnq_auto.derive_cov(flat) ;
	dcovdcov_lnq_auto.inc_dzpuis() ;

	Scalar qq = exp(lnq) ;
	qq.std_spectral_base() ;
	const Vector& dcov_qq = qq.derive_cov(flat) ;

	const Scalar divbeta_auto = beta_auto.divergence(flat) ;
	const Tensor dcov_beta_auto = beta_auto.derive_cov(flat) ;
	Tensor dcovdcov_beta_auto = beta_auto.derive_cov(flat)
	  .derive_cov(flat) ;
	dcovdcov_beta_auto.inc_dzpuis() ;


	// Derivatives of hij, gtilde... 
	//------------------------------

	Scalar psi2 (pow(psi4, 0.5)) ;
	psi2.std_spectral_base() ;

	const Tensor& dcov_hij = hij.derive_cov(flat) ;
	const Tensor& dcon_hij = hij.derive_con(flat) ;
	const Tensor dcov_hij_auto = hij_auto.derive_cov(flat) ;

	Tensor dcovdcov_hij_auto = dcov_hij_auto.derive_cov(flat) ;
	dcovdcov_hij_auto.inc_dzpuis() ;

	const Sym_tensor& gtilde_cov = gtilde.cov() ;
	const Sym_tensor& gtilde_con = gtilde.con() ;
	const Tensor& dcov_gtilde = gtilde_cov.derive_cov(flat) ;

	Connection gamijk (gtilde, flat) ;
	const Tensor& deltaijk = gamijk.get_delta() ;

	// H^i and its derivatives ( = O in Dirac gauge)
	// ---------------------------------------------

	double lambda_dirac = 0. ;

	const Vector hdirac = lambda_dirac * hij.divergence(flat) ;
	const Vector hdirac_auto = lambda_dirac * hij_auto.divergence(flat) ;

	Tensor dcov_hdirac = lambda_dirac * hdirac.derive_cov(flat) ;
	dcov_hdirac.inc_dzpuis() ;
	Tensor dcov_hdirac_auto = lambda_dirac * hdirac_auto.derive_cov(flat) ;
	dcov_hdirac_auto.inc_dzpuis() ;
	Tensor dcon_hdirac_auto = lambda_dirac * hdirac_auto.derive_con(flat) ;
	dcon_hdirac_auto.inc_dzpuis() ;


	//--------------------------------------------------------
	// Poisson equation for logn_auto (nu_auto)
	//--------------------------------------------------------
  
	// Source 
	//--------

	int nr = mp.get_mg()->get_nr(0) ;
	int nt = mp.get_mg()->get_nt(0) ;
	int np = mp.get_mg()->get_np(0) ;
    
	Scalar source1(mp) ;
	Scalar source2(mp) ;
	Scalar source3(mp) ;
	Scalar source4(mp) ;
	Scalar source5(mp) ;
	Scalar source6(mp) ;
	Scalar source7(mp) ;
	Scalar source8(mp) ;
	 
	source1 = qpig * psi4 % (ener_euler + s_euler) ; 

	source2 = psi4 % (kcar_auto + kcar_comp) ;

	source3 = - contract(dcov_logn_auto, 0, dcon_logn, 0, true) 
	    - 2. * contract(contract(gtilde_con, 0, dcov_phi, 0), 
			    0, dcov_logn_auto, 0, true) ;
		    
	source4 = - contract(hij, 0, 1, dcovdcov_logn_auto + 
			     dcov_logn_auto*dcov_logn, 0, 1) ;

	source5 = - contract(hdirac, 0, dcov_logn_auto, 0) ;

	source_tot = source1 + source2 + source3 + source4 + source5 ;
      
 
	cout << "moyenne de la source 1 pour logn_auto" << endl ;
	cout <<  norme(source1/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 2 pour logn_auto" << endl ;
	cout <<  norme(source2/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 3 pour logn_auto" << endl ;
	cout <<  norme(source3/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 4 pour logn_auto" << endl ;
	cout <<  norme(source4/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 5 pour logn_auto" << endl ;
	cout <<  norme(source5/(nr*nt*np)) << endl ;
	cout << "moyenne de la source pour logn_auto" << endl ;
	cout <<  norme(source_tot/(nr*nt*np)) << endl ;

	// Resolution of the Poisson equation 
	// ----------------------------------

	source_tot.poisson(par_logn, logn_auto) ; 
	ssjm1_logn = ssjm1logn ;

	cout << "logn_auto" << endl << norme(logn_auto/(nr*nt*np)) << endl ;
  
	// Check: has the Poisson equation been correctly solved ?
	// -----------------------------------------------------

	Tbl tdiff_logn = diffrel(logn_auto.laplacian(), source_tot) ;
	cout << 
	    "Relative error in the resolution of the equation for logn_auto : "
	     << endl ; 
	for (int l=0; l<nz; l++) {
	    cout << tdiff_logn(l) << "  " ; 
	}
	cout << endl ;
	diff_logn = max(abs(tdiff_logn)) ; 

	//--------------------------------------------------------
	// Poisson equation for lnq_auto
	//--------------------------------------------------------

	// Source
	//--------

	source1 = qpig * psi4 % s_euler ;

	source2 = 0.75 * psi4 % (kcar_auto + kcar_comp) ;

	source3 = - contract(dcon_lnq, 0, dcov_lnq_auto, 0, true) ;

	source4 = 2. * contract(contract(gtilde_con, 0, dcov_phi, 0), 0, 
				dcov_phi_auto + dcov_logn_auto, 0, true) ;

	source5 = 0.0625 * contract(gtilde_con, 0, 1, contract(
		 dcov_hij_auto, 0, 1, dcov_gtilde, 0, 1), 0, 1) ;
		
	source6 = - 0.125 * contract(gtilde_con, 0, 1, contract(
		       dcov_hij_auto, 0, 1, dcov_gtilde, 0, 2), 0, 1) ;

	source7 = - contract(hij, 0, 1, dcovdcov_lnq_auto + dcov_lnq_auto *
			     dcov_lnq, 0, 1) ;

	source8 = - 0.25 * contract(dcov_hdirac_auto, 0, 1) 
	  - contract(hdirac, 0, dcov_lnq_auto, 0) ;
	
	source_tot = source1 + source2 + source3 + source4 + source5 + 
	             source6 + source7 + source8 ;

  	
	cout << "moyenne de la source 1 pour lnq_auto" << endl ;
	cout <<  norme(source1/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 2 pour lnq_auto" << endl ;
	cout <<  norme(source2/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 3 pour lnq_auto" << endl ;
	cout <<  norme(source3/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 4 pour lnq_auto" << endl ;
	cout <<  norme(source4/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 5 pour lnq_auto" << endl ;
	cout <<  norme(source5/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 6 pour lnq_auto" << endl ;
	cout <<  norme(source6/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 7 pour lnq_auto" << endl ;
	cout <<  norme(source7/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 8 pour lnq_auto" << endl ;
	cout <<  norme(source8/(nr*nt*np)) << endl ;
	cout << "moyenne de la source pour lnq_auto" << endl ;
	cout <<  norme(source_tot/(nr*nt*np)) << endl ;
	

	// Resolution of the Poisson equation 
	// ----------------------------------

	source_tot.poisson(par_lnq, lnq_auto) ; 
	ssjm1_lnq = ssjm1lnq ;

	cout << "lnq_auto" << endl << norme(lnq_auto/(nr*nt*np)) << endl ;

	// Check: has the Poisson equation been correctly solved 
	// -----------------------------------------------------
    
	Tbl tdiff_lnq = diffrel(lnq_auto.laplacian(), source_tot) ;
	cout << 
	    "Relative error in the resolution of the equation for lnq : "
	     << endl ; 
	for (int l=0; l<nz; l++) {
	    cout << tdiff_lnq(l) << "  " ; 
	}
	cout << endl ;
	diff_lnq = max(abs(tdiff_lnq)) ; 

	//--------------------------------------------------------
	// Vector Poisson equation for beta_auto
	//--------------------------------------------------------

	// Source
	//--------

	Vector source1_beta(mp, CON, mp.get_bvect_cart()) ;
	Vector source2_beta(mp, CON, mp.get_bvect_cart()) ;
	Vector source3_beta(mp, CON, mp.get_bvect_cart()) ;
	Vector source4_beta(mp, CON, mp.get_bvect_cart()) ;
	Vector source5_beta(mp, CON, mp.get_bvect_cart()) ;
	Vector source6_beta(mp, CON, mp.get_bvect_cart()) ;
	Vector source7_beta(mp, CON, mp.get_bvect_cart()) ;

	source1_beta = (4.*qpig) * nn % psi4
	                   %(ener_euler + press) * u_euler ;

	source2_beta = 2. * nn * contract(tkij_auto, 1, 
					    dcov_logn - 6 * dcov_phi, 0)  ;

	source3_beta = - 2. * nn * contract(tkij_auto, 0, 1, deltaijk, 
    				      1, 2) ;

	source4_beta = - contract(hij, 0, 1, dcovdcov_beta_auto, 1, 2) ;
	
	source5_beta = - 0.3333333333333333 * contract(hij, 1, contract( 
					      dcovdcov_beta_auto, 0, 1), 0) ;
	
	source6_beta = hdirac_auto.derive_lie(omdsdp) ;
	source6_beta.inc_dzpuis() ;

	source7_beta = contract(beta, 0, dcov_hdirac_auto, 1) ;
	       + 2./3. * hdirac * divbeta_auto 
	       - contract(hdirac, 0, dcov_beta_auto, 1) ;

	source_beta = source1_beta + source2_beta + source3_beta 
	  + source4_beta + source5_beta + source6_beta + source7_beta ; 
	

	cout << "moyenne de la source 1 pour beta_auto" << endl ;
	cout <<  norme(source1_beta(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source1_beta(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source1_beta(3)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 2 pour beta_auto" << endl ;
	cout <<  norme(source2_beta(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source2_beta(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source2_beta(3)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 3 pour beta_auto" << endl ;
	cout <<  norme(source3_beta(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source3_beta(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source3_beta(3)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 4 pour beta_auto" << endl ;
	cout <<  norme(source4_beta(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source4_beta(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source4_beta(3)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 5 pour beta_auto" << endl ;
	cout <<  norme(source5_beta(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source5_beta(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source5_beta(3)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 6 pour beta_auto" << endl ;
	cout <<  norme(source6_beta(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source6_beta(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source6_beta(3)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 7 pour beta_auto" << endl ;
	cout <<  norme(source7_beta(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source7_beta(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source7_beta(3)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source pour beta_auto" << endl ;
	cout <<  norme(source_beta(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source_beta(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source_beta(3)/(nr*nt*np)) << endl ;

	// Resolution of the Poisson equation 
	// ----------------------------------

	// Filter for the source of beta vector
	
	for (int i=1; i<=3; i++) {
	  if (source_beta(i).get_etat() != ETATZERO)
	    source_beta.set(i).filtre(4) ;
	}
	
	for (int i=1; i<=3; i++) {
	  if(source_beta(i).dz_nonzero()) {
	    assert( source_beta(i).get_dzpuis() == 4 ) ; 
	  }
	  else{
	    (source_beta.set(i)).set_dzpuis(4) ; 
	  }
	}
	
	double lambda = double(1) / double(3) ; 

	Tenseur source_p(mp, 1, CON, mp.get_bvect_cart() ) ;
	source_p.set_etat_qcq() ;
	for (int i=0; i<3; i++) {
	  source_p.set(i) = Cmp(source_beta(i+1)) ;
	}
	
	Tenseur vect_auxi (mp, 1, CON, mp.get_bvect_cart()) ;
	vect_auxi.set_etat_qcq() ;
	for (int i=0; i<3 ;i++){
	  vect_auxi.set(i) = 0. ;
	}
	Tenseur scal_auxi (mp) ;
	scal_auxi.set_etat_qcq() ;
	scal_auxi.set().annule_hard() ;
	scal_auxi.set_std_base() ;
	
	Tenseur resu_p(mp, 1, CON, mp.get_bvect_cart() ) ;
	resu_p.set_etat_qcq() ;
	source_p.poisson_vect(lambda, par_beta2, resu_p, vect_auxi, 
			      scal_auxi) ;
	
	for (int i=1; i<=3; i++) 
	  beta_auto.set(i) = resu_p(i-1) ;

	ssjm1_khi = ssjm1khi ;
	for (int i=0; i<3; i++){
	    ssjm1_wbeta.set(i+1) = ssjm1wbeta(i) ;
	}
	
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


	if (!conf_flat) {
	   
	    //--------------------------------------------------------
	    // Poisson equation for hij
	    //--------------------------------------------------------
	 
	 
	    // Declaration of all sources 
	    //---------------------------

	    Scalar source_tot_hij(mp) ;
	    Scalar source_Sij(mp) ;
	    Scalar source_Rij(mp) ;

	    Tensor source_1 (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source_2 (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source_3 (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source_4 (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source_5 (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source_6 (mp, 2, CON, mp.get_bvect_cart()) ;

	    Tensor source1_Sij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source2_Sij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source3_Sij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source4_Sij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source5_Sij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source6_Sij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source7_Sij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source8_Sij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source9_Sij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source10_Sij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source11_Sij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source12_Sij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source13_Sij (mp, 2, CON, mp.get_bvect_cart()) ;

	    Tensor source1_Rij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source2_Rij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source3_Rij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source4_Rij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source5_Rij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source6_Rij (mp, 2, CON, mp.get_bvect_cart()) ;
	    Tensor source7_Rij (mp, 2, CON, mp.get_bvect_cart()) ;


	    // Source
	    //--------

	    source_1 = contract(dcon_hij, 1, dcov_lnq_auto, 0) ;

	    source_2 = - contract(dcon_hij, 2, dcov_lnq_auto, 0) 
	      - 2./3. * contract(hdirac, 0, dcov_lnq_auto, 0) * flat.con() ;
	    
	    source_3 = - tkij_auto.derive_lie(omdsdp) ;
	    source_3.inc_dzpuis() ;

	    source_4 = - tkij_auto.derive_lie(beta) ;
	    source_4.inc_dzpuis() ;
	    source_4 += - 2./3. * beta.divergence(flat) * tkij_auto ;

	    source_5 = dcon_hdirac_auto ;
 	    
	    source_6 = - 2./3. * hdirac_auto.divergence(flat) * flat.con() ;
	    source_6.inc_dzpuis() ;

	    // Source terms for Sij
	    //---------------------
	    
	    source1_Sij = 8. * nn / psi4 * phi_auto.derive_con(gtilde) * 
	      contract(gtilde_con, 0, dcov_phi, 0) ;

	    source2_Sij = 4. / psi4 * phi_auto.derive_con(gtilde) * 
	      nn * contract(gtilde_con, 0, dcov_logn, 0) ;

	    source3_Sij = - nn / (3.*psi4) * gtilde_con * 
	      ( 0.25 * contract(gtilde_con, 0, 1, contract(dcov_hij_auto, 0, 1,
					dcov_gtilde, 0, 1), 0, 1)
	       - 0.5 * contract(gtilde_con, 0, 1, contract(dcov_hij_auto, 0, 1,
					dcov_gtilde, 0, 2), 0, 1)) ;

	    source4_Sij = - 8.*nn / (3.*psi4) * gtilde_con * 
	 contract(dcov_phi_auto, 0, contract(gtilde_con, 0, dcov_phi, 0), 0) ;

	    source5_Sij = nn / (3.*psi4) * hdirac.divergence(flat)*hij_auto ;
	    source5_Sij.inc_dzpuis() ;

	    source6_Sij = - 8./(3.*psi4) * contract(dcov_phi_auto, 0,
		nn*contract(gtilde_con, 0, dcov_logn, 0), 0) * gtilde_con ;

	    source7_Sij = 2.*nn* contract(gtilde_cov, 0, 1, tkij_auto *
					   (tkij_auto+tkij_comp), 1, 3) ;

	    source8_Sij = - 2. * qpig * nn * ( psi4 * stress_euler 
			      - 0.33333333333333333 * s_euler * gtilde_con ) ; 
	    
	    source9_Sij = - 1./(psi4*psi2) * contract(gtilde_con, 1, 
			  contract(gtilde_con, 1, qq*dcovdcov_lnq_auto + 
				   qq*dcov_lnq_auto*dcov_lnq, 0), 1) ;

	    source10_Sij = - 0.5/(psi4*psi2) * contract(contract(hij, 1,
					dcov_hij_auto, 2), 1, dcov_qq, 0) ;
					
	    source11_Sij = 0.5/(psi4*psi2) * contract(contract(hij, 0,
					dcov_hij_auto, 2), 0, dcov_qq, 0) ;

	    source12_Sij = 1./(3.*psi4*psi2)*contract(gtilde_con, 0, 1,
                    qq*dcovdcov_lnq_auto + qq*dcov_lnq_auto*dcov_lnq, 0, 1)
	                          *gtilde_con ;

	    source13_Sij = 1./(3.*psi4*psi2) * contract(hdirac, 0, 
							dcov_qq, 0)*hij_auto ;

	    // Source terms for Rij
	    //---------------------

	    source1_Rij = contract(hij, 0, 1, dcovdcov_hij_auto, 2, 3) ;

	    source2_Rij = - contract(hij_auto, 1, dcov_hdirac, 1) ;
	    
	    source3_Rij = contract(hdirac, 0, dcov_hij_auto, 2) ;

	    source4_Rij = - contract(contract(dcov_hij_auto, 1, dcov_hij, 2),
				     1, 3) ;

	    source5_Rij = - contract(gtilde_cov, 0, 1, contract(contract(
		    gtilde_con, 0, dcov_hij_auto, 2), 0, dcov_hij, 2), 1, 3) ;

	    source6_Rij = contract(contract(contract(contract(gtilde_cov, 0, 
		 dcov_hij_auto, 1), 2, gtilde_con, 1), 0, dcov_hij, 1), 0, 3) ;

	    source7_Rij = 0.5 * contract(gtilde_con*gtilde_con, 1, 3, 
		   contract(dcov_hij_auto, 0, 1, dcov_gtilde, 0, 1), 0, 1) ;


	    for(int i=1; i<=3; i++) 
		for(int j=1; j<=i; j++) {

		    source_Sij = source1_Sij(i,j) + source2_Sij(i,j) +
		      source2_Sij(j,i) + source3_Sij(i,j) +
		      source4_Sij(i,j) + source5_Sij(i,j) +
		      source6_Sij(i,j) + source7_Sij(i,j) +
		      source8_Sij(i,j) + source9_Sij(i,j) +
		      source10_Sij(i,j) + source10_Sij(j,i) +
		      source11_Sij(i,j) + source12_Sij(i,j) +
		      source13_Sij(i,j) ;

		    source_Rij = 0.5*(source1_Rij(i,j) + source2_Rij(i,j) +
				      source2_Rij(j,i) + source3_Rij(i,j) +
				      source4_Rij(i,j) + source5_Rij(i,j) +
				      source6_Rij(i,j) + source6_Rij(j,i) +
				      source7_Rij(i,j)) ;

		    source_tot_hij = source_1(i,j) + source_1(j,i) 
			+ source_2(i,j) + 2.*psi4/nn * (source_3(i,j) + 
			    source_4(i,j) - source_Sij) - 2.* source_Rij +
		      source_5(i,j) + source_5(j,i) + source_6(i,j) ;
			
		    
		    cout << "source1 de Sij" << endl 
			 << norme(source1_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source2 de Sij" << endl 
			 << norme(source2_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source3 de Sij" << endl 
			 << norme(source3_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source4 de Sij" << endl 
			 << norme(source4_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source5 de Sij" << endl 
			 << norme(source5_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source6 de Sij" << endl 
			 << norme(source6_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source7 de Sij" << endl 
			 << norme(source7_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source8 de Sij" << endl 
			 << norme(source8_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source9 de Sij" << endl 
			 << norme(source9_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source10 de Sij" << endl 
			 << norme(source10_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source11 de Sij" << endl 
			 << norme(source11_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source12 de Sij" << endl 
			 << norme(source12_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source13 de Sij" << endl 
			 << norme(source13_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source1" << endl 
			 << norme(source_1(i,j)/(nr*nt*np)) << endl ;
		    cout << "source2" << endl 
			 << norme(source_2(i,j)/(nr*nt*np)) << endl ;
		    cout << "source3" << endl 
			 << norme(source_3(i,j)/(nr*nt*np)) << endl ;
		    cout << "source4" << endl 
			 << norme(source_4(i,j)/(nr*nt*np)) << endl ;
		    cout << "source5" << endl 
			 << norme(source_5(i,j)/(nr*nt*np)) << endl ;
		    cout << "source6" << endl 
			 << norme(source_6(i,j)/(nr*nt*np)) << endl ;
		    cout << "source_Rij" << endl 
			 << norme(source_Rij/(nr*nt*np)) << endl ;
		    cout << "source_Sij" << endl 
			 << norme(source_Sij/(nr*nt*np)) << endl ;
		    cout << "source_tot" << endl 
			 << norme(source_tot_hij/(nr*nt*np)) << endl ;
		    	
		    // Resolution of the Poisson equations and
		    // Check: has the Poisson equation been correctly solved ?
		    // -----------------------------------------------------

		    if(i==1 && j==1) {
		 
			source_tot_hij.poisson(par_h11, hij_auto.set(1,1)) ; 

			Tbl tdiff_h11 = diffrel(hij_auto(1,1).laplacian(), 
						source_tot_hij) ;  
			cout << "Relative error in the resolution of the equation for "
			     << "h11_auto : " << endl ; 
			for (int l=0; l<nz; l++) {
			    cout << tdiff_h11(l) << "  " ; 
			}
			cout << endl ;
			diff_h11 = max(abs(tdiff_h11)) ; 
		    }
	       	       
		    if(i==2 && j==1) {

			source_tot_hij.poisson(par_h21, hij_auto.set(2,1)) ; 
	    
			Tbl tdiff_h21 = diffrel(hij_auto(2,1).laplacian(), 
						source_tot_hij) ;
			cout << 
			    "Relative error in the resolution of the equation for " 
			     << "h21_auto : "  << endl ; 
			for (int l=0; l<nz; l++) {
			    cout << tdiff_h21(l) << "  " ; 
			}
			cout << endl ;
			diff_h21 = max(abs(tdiff_h21)) ; 
		    }
	       
		    if(i==3 && j==1) {
		 
			source_tot_hij.poisson(par_h31, hij_auto.set(3,1)) ; 

			Tbl tdiff_h31 = diffrel(hij_auto(3,1).laplacian(), 
						source_tot_hij) ;
			cout << 
			    "Relative error in the resolution of the equation for "
			     << "h31_auto : " << endl ; 
			for (int l=0; l<nz; l++) {
			    cout << tdiff_h31(l) << "  " ; 
			}
			cout << endl ;
			diff_h31 = max(abs(tdiff_h31)) ; 
		    }
	     
		    if(i==2 && j==2) {
		 
			source_tot_hij.poisson(par_h22, hij_auto.set(2,2)) ; 

			Tbl tdiff_h22 = diffrel(hij_auto(2,2).laplacian(), 
						source_tot_hij) ;
			cout << 
			    "Relative error in the resolution of the equation for "
			     << "h22_auto : " << endl ; 
			for (int l=0; l<nz; l++) {
			    cout << tdiff_h22(l) << "  " ; 
			}
			cout << endl ;
			diff_h22 = max(abs(tdiff_h22)) ; 
		    }
	       
		    if(i==3 && j==2) {
		 
			source_tot_hij.poisson(par_h32, hij_auto.set(3,2)) ; 
	
			Tbl tdiff_h32 = diffrel(hij_auto(3,2).laplacian(), 
						source_tot_hij) ;
			cout << 
			    "Relative error in the resolution of the equation for "
			     << "h32_auto : " << endl ; 
			for (int l=0; l<nz; l++) {
			    cout << tdiff_h32(l) << "  " ; 
			}
			cout << endl ;
			diff_h32 = max(abs(tdiff_h32)) ; 
		    }
	     
		    if(i==3 && j==3) {
		 
			source_tot_hij.poisson(par_h33, hij_auto.set(3,3)) ; 

			Tbl tdiff_h33 = diffrel(hij_auto(3,3).laplacian(), 
						source_tot_hij) ;
			cout << 
			    "Relative error in the resolution of the equation for "
			     << "h33_auto : " << endl ;
			for (int l=0; l<nz; l++) {
			    cout << tdiff_h33(l) << "  " ;
			}
			cout << endl ;
			diff_h33 = max(abs(tdiff_h33)) ;
		    }

		}
      
	    cout << "Tenseur hij auto in cartesian coordinates" << endl ;
	    for (int i=1; i<=3; i++)
		for (int j=1; j<=i; j++) {
		    cout << "  Comp. " << i << " " << j << " :  " ;
		    for (int l=0; l<nz; l++){
			cout << norme(hij_auto(i,j)/(nr*nt*np))(l) << " " ;
		    }
		    cout << endl ;
		}
	    cout << endl ;

	    ssjm1_h11 = ssjm1h11 ;
	    ssjm1_h21 = ssjm1h21 ;
	    ssjm1_h31 = ssjm1h31 ;
	    ssjm1_h22 = ssjm1h22 ;
	    ssjm1_h32 = ssjm1h32 ;
	    ssjm1_h33 = ssjm1h33 ;

	}

	// End of relativistic equations	
	   
      	   
	//-------------------------------------------------
	//  Relative change in enthalpy
	//-------------------------------------------------

	Tbl diff_ent_tbl = diffrel( ent, ent_jm1 ) ; 
	diff_ent = diff_ent_tbl(0) ; 
	for (int l=1; l<nzet; l++) {
	    diff_ent += diff_ent_tbl(l) ; 
	}
	diff_ent /= nzet ; 
	
	
	ent_jm1 = ent ; 


        } // End of main loop
    
    //=========================================================================
    // 			End of iteration
    //=========================================================================

}  





    
    
