
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
			   const Tbl& fact_resize, Tbl& diff, int step_coal,
			   double omega) {

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
    double& diff_beta = diff.set(3) ; 
    double& diff_shift_r = diff.set(4) ; 
    double& diff_shift_t = diff.set(5) ; 
    double& diff_shift_p = diff.set(6) ; 
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
 
        
    // Computation of all derived quantities. It must be performed before
    // the adaptation of the mapping to avoid discontinuities. 
    // Derivatives of N and logN
    //--------------------------
    
    Vector dcov_logn_auto = logn_auto.derive_cov(flat) ;
    Vector dcon_logn_auto = logn_auto.derive_con(flat) ;
    
    const Vector& dcov_nnn = nnn.derive_cov(flat) ;
    const Tensor& dcovdcov_nnn = dcov_nnn.derive_cov(flat) ;
    
    // Derivatives of beta, Q  and shift 
    //----------------------------------
    
    Vector dcov_lnpsi_auto = 0.5*(beta_auto - logn_auto).derive_cov(flat) ;
    Vector dcon_lnpsi_auto = 0.5*(beta_auto - logn_auto).derive_con(flat) ;
    
    const Vector& dcov_beta = beta.derive_cov(flat) ;
    const Vector& dcon_beta = beta.derive_con(flat) ;
    Tensor dcovdcov_beta = dcov_beta.derive_cov(flat) ;
    dcovdcov_beta.inc_dzpuis() ;
    Tensor dcondcov_beta = dcov_beta.derive_con(flat) ;
    dcondcov_beta.inc_dzpuis() ;
    
    const Vector& dcov_beta_auto = beta_auto.derive_cov(flat) ;
    
    const Tensor& dcov_shift = shift.derive_cov(flat) ;
    const Tensor& dcovdcov_shift = dcov_shift.derive_cov(flat) ;
    
    Scalar psi2 (pow(psi4, 0.5)) ;
    psi2.std_spectral_base() ;
    
    Scalar qq = exp(beta) ;
    qq.std_spectral_base() ;
    
    const Vector& dcov_qq = qq.derive_cov(flat) ;
    Tensor dcovdcov_qq = dcov_qq.derive_cov(flat) ;
    dcovdcov_qq.inc_dzpuis() ;
    Tensor dcondcov_qq = dcov_qq.derive_con(flat) ;
    dcondcov_qq.inc_dzpuis() ;
    
    
    // Derivatives of hij, gtilde... 
    //------------------------------
    
    const Tensor& dcov_hij = hij.derive_cov(flat) ;
    const Tensor& dcov_hij_auto = hij_auto.derive_cov(flat) ;
    const Tensor& dcon_hij_auto = hij_auto.derive_con(flat) ;
    
    Tensor dcovdcov_hij = dcov_hij.derive_cov(flat) ;
    dcovdcov_hij.inc_dzpuis() ;
    
    Sym_tensor gtilde_cov = gtilde.cov() ;
    Sym_tensor gtilde_con = gtilde.con() ;
    const Tensor& dcov_gtilde = gtilde_cov.derive_cov(flat) ;
    
    Connection gamijk (gtilde, flat) ;
    const Tensor& deltaijk = gamijk.get_delta() ;
    
    
    // External potential
    // See Eq (99) from Gourgoulhon et al. (2001)
    // ------------------
    
    Scalar pot_ext = logn_comp + pot_centri + loggam ;
    
    Scalar ent_jm1 = ent ;	// Enthalpy at previous step
    
    Scalar source_tot(mp) ; // source term in the equation for hij_auto, 
                            // logn_auto and beta_auto
			    
    Vector source_shift(mp, CON, mp.get_bvect_spher()) ;  // source term 
    // in the equation for shift_auto



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
	
	double precis_poisson = 1e-16 ;
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

	    if (nz == 5 && nzet == 1) {
		double rr_in_1 = mp.val_r(1,-1., M_PI/2, 0.) ; 
		double rr_out_1 = mp.val_r(1, 1., M_PI/2, 0.) ; 
		double rr_out_2 = mp.val_r(2, 1., M_PI/2, 0.) ; 
		double rr_out_3 = mp.val_r(3, 1., M_PI/2, 0.) ; 
		
		double fact_resize_0 ;
		if (fact_resize(0) > 2.4) fact_resize_0 = fact_resize(0)/2. ;
		else fact_resize_0 = fact_resize(0)/2. + 0.5 ;


		mp.resize(1, rr_in_1/rr_out_1 * fact_resize_0) ; 
		mp.resize(2, rr_in_1/rr_out_2 * fact_resize(0)) ; 
		mp.resize(3, rr_in_1/rr_out_3 * fact_resize(1)) ; 
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
	Scalar source9(mp) ;
	 
	source1 = qpig * psi4 % (ener_euler + s_euler) ; 

	source2 = psi4 % (kcar_auto + kcar_comp) ;

	source3 = - contract(dcov_logn_auto, 0, dcon_logn, 0, true) 
	-2. * contract(dcov_lnpsi, 0, dcon_logn_auto, 0, true) ;
		    
	source4 = - contract(hij_auto, 0, 1, dcovdcov_nnn/nnn, 0, 1) ;
	source4.inc_dzpuis(1) ;
	source4.annule_domain(nz-1) ;

	source5 = - 2. * contract(hij_auto, 0, 1, dcov_lnpsi 
				  * dcov_logn, 0, 1) ;

	source_tot = source1 + source2 + source3 + source4 + source5 ;

	Cmp source_cmp (source_tot) ;
	source_cmp.set_dzpuis(0) ;
	mp.reevaluate_symy(&mp_prev, nzet+1, source_cmp) ;
	source_cmp.set_dzpuis(4) ;
	source_tot.set_domain(0) = source_cmp(0) ;
	source_tot.set_domain(1) = source_cmp(1) ;

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

//	source_tot.filtre(4) ;
	
	// Construction of the affine mapping
	// ----------------------------------

	double* bornes = new double[nz+1];
	    bornes[0] = 0. ; 
	for (int l=1; l<nz; l++) 
	    bornes[l] = mp.val_r(l, -1., M_PI/2., 0.) ;
	    bornes[nz] = __infinity ; 
	Map_af mpaff (*mp.get_mg(), bornes) ;
	mpaff.set_ori(mp.get_ori_x(), mp.get_ori_y(), mp.get_ori_z()) ; 
	mpaff.set_rot_phi(mp.get_rot_phi()) ;

	// Source on the affine mapping
	Scalar source_aff(mpaff) ;
	Scalar temp (source_tot) ;
	temp.set_dzpuis(0) ;
	temp.annule(2, nz-1) ;
	source_aff.import(temp) ;
	source_aff.std_spectral_base() ;
	source_aff.set_dzpuis(4) ;

	for (int l=2; l<nz; l++)
	    for (int k=0; k<np; k++)
		for (int j=0; j<nt; j++)
		    for (int i=0; i<nr; i++)
			source_aff.set_grid_point(l,k,j,i) = 
			    source_tot.val_grid_point(l,k,j,i) ;

	Scalar logn_aff (mpaff) ;
	logn_aff = 0. ;

//	des_meridian(source_aff, 0., 10., "source_aff", 10) ; 

	// Resolution of the Poisson equation 
	// ----------------------------------

	logn_aff = source_aff.poisson() ; 
	logn_auto.import(logn_aff) ;

//	des_meridian(logn_aff, 0., 10., "logn_aff", 11) ; 

	cout << "logn_auto" << endl << norme(logn_auto/(nr*nt*np)) << endl ;
  
	// Check: has the Poisson equation been correctly solved ?
	// -----------------------------------------------------

	Tbl tdiff_logn = diffrel(logn_auto.laplacian(), source_tot) ;

	cout << "logn_auto.laplacian()" << endl 
	     << norme(logn_auto.laplacian())  << endl ;
	cout << "source" << endl << norme(source_tot)  << endl ;
//	des_profile(logn_auto.laplacian(), 0., 4., 1., 1., "logn.laplacian") ; 
//	des_profile(logn_auto.laplacian()- source_tot, 0., 4., 1., 1., "diff") ; 

	cout << 
	    "Relative error in the resolution of the equation for logn_auto : "
	     << endl ; 
	for (int l=0; l<nz; l++) {
	    cout << tdiff_logn(l) << "  " ; 
	}
	cout << endl ;
	diff_logn = max(abs(tdiff_logn)) ; 

	//--------------------------------------------------------
	// Poisson equation for beta_auto
	//--------------------------------------------------------

	// Source
	//--------

	source1 = qpig * psi4 % s_euler ;

	source2 = 0.75 * psi4 % (kcar_auto + kcar_comp) ;
	
	source3 = - 0.5 * contract(dcov_logn_auto, 0, dcon_logn, 0, true) ;

	source4 = - 0.5 * contract(dcov_beta_auto, 0, dcon_beta, 0, true) ;

	source5 = 0.0625 * contract(gtilde_con, 0, 1, contract(
			      dcov_hij_auto, 0, 1, dcov_gtilde, 0, 1), 0, 1) ;
			   
	source6 = - 0.125 * contract(gtilde_con, 0, 1, contract(dcov_hij_auto, 
					  0, 1, dcov_gtilde, 0, 2), 0, 1) ;

	source7 = 2. * contract(hij_auto, 0, 1, dcov_lnpsi * dcov_lnpsi, 
				0, 1) ;

	source8 = 2 * contract(hij_auto, 0, 1, dcov_lnpsi * dcov_logn, 0, 1) ;
	
	source9 = - contract(hij_auto, 0, 1, dcovdcov_qq, 0, 1) / qq ;

	source_tot = source1 + source2 + source3 + source4 + source5 + 
	             source6 + source7 + source8 + source9 ;

	source_cmp = source_tot ;
	source_cmp.set_dzpuis(0) ;
	mp.reevaluate_symy(&mp_prev, nzet+1, source_cmp) ;
	source_cmp.set_dzpuis(4) ;
	source_tot.set_domain(0) = source_cmp(0) ;
	source_tot.set_domain(1) = source_cmp(1) ;
	
	cout << "moyenne de la source 1 pour beta_auto" << endl ;
	cout <<  norme(source1/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 2 pour beta_auto" << endl ;
	cout <<  norme(source2/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 3 pour beta_auto" << endl ;
	cout <<  norme(source3/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 4 pour beta_auto" << endl ;
	cout <<  norme(source4/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 5 pour beta_auto" << endl ;
	cout <<  norme(source5/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 6 pour beta_auto" << endl ;
	cout <<  norme(source6/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 7 pour beta_auto" << endl ;
	cout <<  norme(source7/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 8 pour beta_auto" << endl ;
	cout <<  norme(source8/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 9 pour beta_auto" << endl ;
	cout <<  norme(source9/(nr*nt*np)) << endl ;
	cout << "moyenne de la source pour beta_auto" << endl ;
	cout <<  norme(source_tot/(nr*nt*np)) << endl ;
	
//	source_tot.filtre(4) ;

	// Source on the affine mapping
	temp = source_tot ;
	temp.set_dzpuis(0) ;
	temp.annule(2, nz-1) ;
	source_aff.import(temp) ;
	source_aff.std_spectral_base() ;
	source_aff.set_dzpuis(4) ;

	for (int l=2; l<nz; l++)
	    for (int k=0; k<np; k++)
		for (int j=0; j<nt; j++)
		    for (int i=0; i<nr; i++)
			source_aff.set_grid_point(l,k,j,i) = 
			    source_tot.val_grid_point(l,k,j,i) ;

	Scalar beta_aff (mpaff) ;
	beta_aff = 0. ;

//	des_meridian(source_aff, 0., 10., "source_aff", 14) ; 

	// Resolution of the Poisson equation 
	// ----------------------------------

	beta_aff = source_aff.poisson() ; 
	beta_auto.import(beta_aff) ;

//	des_meridian(beta_aff, 0., 10., "beta_aff", 11) ; 
	cout << "beta_auto" << endl << norme(beta_auto/(nr*nt*np)) << endl ;

	// Check: has the Poisson equation been correctly solved 
	// -----------------------------------------------------
    
	Tbl tdiff_beta = diffrel(beta_auto.laplacian(), source_tot) ;

	cout << "beta_auto.laplacian()" << endl 
	     << norme(beta_auto.laplacian())  << endl ;
	cout << "source" << endl << norme(source_tot)  << endl ;
//	des_profile(beta_auto.laplacian(), 0., 4., 1., 1., "beta.laplacian") ; 
//	des_profile(beta_auto.laplacian()- source_tot, 0., 4., 1., 1., "diff") ; 
	cout << 
	    "Relative error in the resolution of the equation for beta : "
	     << endl ; 
	for (int l=0; l<nz; l++) {
	    cout << tdiff_beta(l) << "  " ; 
	}
	cout << endl ;
	diff_beta = max(abs(tdiff_beta)) ; 

	//--------------------------------------------------------
	// Vector Poisson equation for shift_auto
	//--------------------------------------------------------

	// Source
	//-------

	Vector source1_shift(mp, CON, mp.get_bvect_spher()) ;
	Vector source2_shift(mp, CON, mp.get_bvect_spher()) ;
	Vector source3_shift(mp, CON, mp.get_bvect_spher()) ;
	Vector source4_shift(mp, CON, mp.get_bvect_spher()) ;
	Vector source5_shift(mp, CON, mp.get_bvect_spher()) ;

	u_euler.change_triad(mp.get_bvect_spher()) ;

	source1_shift = (4.*qpig) * nnn % psi4
	                   %(ener_euler + press) * u_euler ;
	
	source2_shift = 2. * nnn * contract(tkij_auto, 1, 
					dcov_logn - 6 * dcov_lnpsi, 0)  ;

	source3_shift = - 2. * nnn * contract(tkij_auto, 0, 1, deltaijk, 
					      1, 2) ;

	source4_shift = - contract(hij_auto, 0, 1, dcovdcov_shift, 1, 2) ;

	source4_shift.annule_domain(nz-1) ;
	source4_shift.inc_dzpuis(1) ;
	
	source5_shift = - 0.3333333333333333 * contract(contract(hij_auto, 1, 
					 dcovdcov_shift, 2), 1, 2) ;
	
	source5_shift.annule(nz-1, nz-1) ;
	source5_shift.inc_dzpuis(1) ;

	source_shift = source1_shift + source2_shift + source3_shift 
	    + source4_shift + source5_shift ;

/*	
	source_shift.change_triad(mp.get_bvect_cart()) ;
	Cmp source1_cmp (source_shift(1)) ;
	Cmp source2_cmp (source_shift(2)) ;
	Cmp source3_cmp (source_shift(3)) ;
	source1_cmp.set_dzpuis(0) ;
	source2_cmp.set_dzpuis(0) ;
	source3_cmp.set_dzpuis(0) ;
	mp.reevaluate(&mp_prev, nzet+1, source1_cmp) ;
	mp.reevaluate(&mp_prev, nzet+1, source2_cmp) ;
	mp.reevaluate(&mp_prev, nzet+1, source3_cmp) ;
	source1_cmp.set_dzpuis(4) ;
	source2_cmp.set_dzpuis(4) ;
	source3_cmp.set_dzpuis(4) ;
	source_shift.set(1) = source1_cmp ;
	source_shift.set(2) = source2_cmp ;
	source_shift.set(3) = source3_cmp ;
	source_shift.change_triad(mp.get_bvect_spher()) ;
*/

/*
	source1_shift.change_triad(mp.get_bvect_cart()) ;
	source2_shift.change_triad(mp.get_bvect_cart()) ;
	source_shift.change_triad(mp.get_bvect_cart()) ;
	
	cout << "moyenne de la source 1 pour shift_auto" << endl ;
	cout <<  norme(source1_shift(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source1_shift(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source1_shift(3)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 2 pour shift_auto" << endl ;
	cout <<  norme(source2_shift(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source2_shift(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source2_shift(3)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 3 pour shift_auto" << endl ;
	cout <<  norme(source3_shift(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source3_shift(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source3_shift(3)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 4 pour shift_auto" << endl ;
	cout <<  norme(source4_shift(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source4_shift(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source4_shift(3)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 5 pour shift_auto" << endl ;
	cout <<  norme(source5_shift(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source5_shift(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source5_shift(3)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source pour shift_auto" << endl ;
	cout <<  norme(source_shift(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source_shift(2)/(nr*nt*np)) << endl ;
	cout <<  norme(source_shift(3)/(nr*nt*np)) << endl ;

	source1_shift.change_triad(mp.get_bvect_spher()) ;
	source2_shift.change_triad(mp.get_bvect_spher()) ;
	source_shift.change_triad(mp.get_bvect_spher()) ;
*/	

	// Source on the affine mapping
	Vector temp_vect (source_shift) ;
	for(int i=1; i<=3; i++)
	    temp_vect.set(i).set_dzpuis(0) ;
	temp_vect.std_spectral_base() ;
	temp_vect.annule(2, nz-1) ;
	temp_vect.change_triad(mp.get_bvect_cart()) ;
	Vector source_aff_vect (mpaff, CON, mp.get_bvect_cart()) ;
	for(int i=1; i<=3; i++)
	    source_aff_vect.set(i).import(temp_vect(i)) ;
	source_aff_vect.std_spectral_base() ;
	source_aff_vect.change_triad(mpaff.get_bvect_spher()) ;
	for(int i=1; i<=3; i++)
	    source_aff_vect.set(i).set_dzpuis(4) ;

	for (int l=2; l<nz; l++)
	    for (int k=0; k<np; k++)
		for (int j=0; j<nt; j++)
		    for (int i=0; i<nr; i++)
			for (int lig=1; lig<=3; lig++)
			    source_aff_vect.set(lig).set_grid_point(l,k,j,i) = 
				source_shift(lig).val_grid_point(l,k,j,i) ;

        Vector shift_aff (mpaff, CON, mpaff.get_bvect_spher()) ;

//	des_meridian(source_aff_vect(1), 0., 10., "source_aff(1)", 22) ; 

	// Resolution of the Poisson equation 
	// ----------------------------------

	for (int i=1; i<=3; i++) {
	    if(source_shift(i).dz_nonzero()) {
		assert( source_shift(i).get_dzpuis() == 4 ) ; 
	    }
	    else{
		(source_shift.set(i)).set_dzpuis(4) ; 
	    }
	}

	double lambda = double(1) / double(3) ; 


	shift_aff = source_aff_vect.poisson(lambda, 2) ; 

	shift_auto.change_triad(mpaff.get_bvect_cart()) ;
	shift_aff.change_triad(mpaff.get_bvect_cart()) ;
	for(int i=1; i<=3; i++)
	    shift_auto.set(i).import(shift_aff(i)) ;
	shift_auto.std_spectral_base() ;
	shift_auto.change_triad(mp.get_bvect_spher()) ;

//	des_meridian(shift_aff(1), 0., 10., "beta_aff", 23) ; 

	cout << "shift_auto(r)" << endl << norme(shift_auto(1)/(nr*nt*np)) << endl ;
	cout << "shift_auto(t)" << endl << norme(shift_auto(2)/(nr*nt*np)) << endl ;
	cout << "shift_auto(p)" << endl << norme(shift_auto(3)/(nr*nt*np)) << endl ;
  

	// Check: has the equation for shift_auto been correctly solved ?
	// --------------------------------------------------------------
	
	Vector lap_shift = (shift_auto.derive_con(flat)).divergence(flat) 
	    + lambda* shift_auto.divergence(flat).derive_con(flat) ;
	
	source_shift.dec_dzpuis() ;
	Tbl tdiff_shift_r = diffrel(lap_shift(1), source_shift(1)) ; 
	Tbl tdiff_shift_t = diffrel(lap_shift(2), source_shift(2)) ; 
	Tbl tdiff_shift_p = diffrel(lap_shift(3), source_shift(3)) ; 


	cout << "norme de lap_shift" << endl << norme(lap_shift(1)) 
	     << norme(lap_shift(2)) << norme(lap_shift(3)) << endl ; 
	cout << "norme de source" << endl << norme(source_shift(1)) 
	     << norme(source_shift(2)) << norme(source_shift(3)) << endl ; 
//	des_profile(lap_shift(1), 0., 10., 1., 1., "lap_shift(1)") ;
//	des_profile(source_shift(1) - lap_shift(1), 0., 10., 1., 1., "diff_shift(1)") ;

	cout << 
	    "Relative error in the resolution of the equation for shift_auto : "
	     << endl ; 
	cout << "r component : " ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_shift_r(l) << "  " ; 
	}
	cout << endl ;
	cout << "t component : " ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_shift_t(l) << "  " ; 
	}
	cout << endl ;
	cout << "p component : " ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_shift_p(l) << "  " ; 
	}
	cout << endl ;
	
	diff_shift_r = max(abs(tdiff_shift_r)) ; 
	diff_shift_t = max(abs(tdiff_shift_t)) ; 
	diff_shift_p = max(abs(tdiff_shift_p)) ; 


	if (!conf_flat) {
	   
	    //--------------------------------------------------------
	    // Poisson equation for hij
	    //--------------------------------------------------------
	 
	 
	    // Declaration of all sources 
	    //---------------------------

	    Scalar source_Hij(mp) ;
	    Scalar source_Qij(mp) ;
	    Scalar source_Sij(mp) ;
	    Scalar source_tot_hij(mp) ;


	    Tensor source_1 (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source_2 (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source_3 (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source_4 (mp, 2, CON, mp.get_bvect_spher()) ;


	    Tensor source1_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source2_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source3_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source4_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source5_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source6_Hij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source7_Hij (mp, 2, CON, mp.get_bvect_spher()) ;

	    
	    Tensor source1_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source2_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source3_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source4_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source5_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source6_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source7_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source8_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source9_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source10_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source11_Qij (mp, 2, CON, mp.get_bvect_spher()) ;
	    

	    Tensor source1_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source2_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source3_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source4_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source5_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source6_Sij (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source7_Sij (mp, 2, CON, mp.get_bvect_spher()) ;


	    // Source
	    //--------

	    source_1 = contract(dcon_hij_auto, 1, dcov_logn + 2*dcov_lnpsi, 0) ;

	    source_2 = - contract(dcon_hij_auto, 2, dcov_logn + 2*dcov_lnpsi, 0) ;
	    
	    source_3 =  tkij_auto.derive_lie(shift) ;
	    source_3.inc_dzpuis() ;

	    source_4 = 0.666666666666667 * shift.divergence(flat) * tkij_auto ;

	    
	    // Source terms for Hij
	    //---------------------

	    source1_Hij = 8.*nnn*contract(hij_auto,1,dcov_lnpsi,0)*dcon_lnpsi ;

	    source2_Hij = - contract(hij_auto, 1 , dcondcov_qq, 0) / psi2 ;

	    source3_Hij = 4.*nnn*contract(hij_auto,1,dcov_logn,0)*dcon_lnpsi ; 

	    source4_Hij = 4.*nnn*contract(hij_auto,1,dcov_lnpsi,0)*dcon_logn ;

	    source5_Hij = ( contract(hij_auto, 0, 1, dcovdcov_qq, 0, 1)/psi2
			   - 8. * nnn * contract(contract(hij_auto,0
					  , dcov_lnpsi, 0), 0, dcov_logn, 0)
			   - 8. * nnn * contract(contract(hij_auto,0
					  , dcov_lnpsi, 0), 0, dcov_lnpsi, 0))
		* flat.con() * 0.3333333333333333 ; 

	    source6_Hij = ( qq.laplacian() / psi2 
			   - 8. * nnn * contract(dcov_lnpsi, 0, dcon_logn, 0)
			   - 8. * nnn * contract(dcov_lnpsi, 0, dcon_lnpsi, 0))
		* hij_auto * 0.3333333333333333 ;

	    source7_Hij = 0.6666666666666667 * qpig * nnn * s_euler * hij ;


	    // Source_terms for Qij
	    //---------------------

	    
	    source1_Qij = 0.5*nnn*contract(hij_auto, 0, 1,dcovdcov_hij, 2, 3) ;


	    source2_Qij = 0.5 * nnn * ( - contract(contract(dcov_hij_auto, 1,  
			  dcov_hij, 2), 1, 3)
	     - contract(contract(contract(gtilde_con, 0, 
		  dcov_hij_auto, 2), 0, dcov_hij, 2), 1, 3, gtilde_cov, 0, 1)
	     + contract(contract(contract(contract(gtilde_con, 1, 
		  dcov_hij_auto, 2), 2, dcov_hij, 2), 1, gtilde_cov, 0), 2, 3)
	     + contract(contract(dcov_gtilde, 0, 1, dcov_hij_auto, 0, 1), 
	     0, 1, gtilde_con * gtilde_con, 1, 3)*0.5 ) ;
	      
	    source3_Qij = 0.5 * nnn * contract(contract(contract(
                  contract(gtilde_con, 1, dcov_hij_auto, 2), 1, 
		  dcov_hij, 2), 1, gtilde_cov, 1), 2, 3) ;

	    source4_Qij = nnn * 8. * contract(contract(hij_auto, 1, 
				        dcov_lnpsi, 0)*hij, 2, dcov_lnpsi, 0) ;
	    
	    source5_Qij = - contract(contract(hij_auto, 1, 
			dcovdcov_qq, 1), 1, hij, 1) / psi2 ;

	    source6_Qij = - 0.5 * contract(contract(hij_auto, 1, dcov_hij, 2)
				     , 2, dcov_qq, 0) / psi2 ;


	    source7_Qij = 0.5 * contract(contract(hij_auto, 1, dcov_hij, 2), 
				     0, dcov_qq, 0) /psi2 ;


	    source8_Qij = 4. * nnn * contract(contract(hij_auto, 1, 
				       dcov_lnpsi, 0)*hij, 2, dcov_logn, 0)

		        + 4. * nnn * contract(contract(hij_auto, 1, 
				       dcov_logn, 0)*hij, 2, dcov_lnpsi, 0) ;

	    source9_Qij = nnn * (contract(contract(dcov_hij_auto, 0, 1, 
			 dcov_gtilde, 0, 2), 0, 1, gtilde_con, 0, 1)
				    -  contract(contract(dcov_hij_auto, 0, 1, 
		       	 dcov_gtilde, 0, 1), 0, 1, gtilde_con, 0, 1) * 0.5 )
		* 0.1666666666666667 * gtilde_con ;

	    
	    source10_Qij = - 2.666666666666667*nnn*( contract(contract(hij, 0, 
			        dcov_lnpsi, 0), 0, dcov_lnpsi, 0) * hij_auto 
					   + contract(contract(hij, 0, 
				dcov_lnpsi, 0), 0, dcov_logn, 0) * hij_auto ) ;


	    source11_Qij = contract(hij, 0, 1, dcovdcov_qq, 0, 1) / psi2
				    * hij_auto * 0.3333333333333333 ;


	    // Source terms for Sij
	    //---------------------


	    source1_Sij = 8. * nnn * dcon_lnpsi_auto * dcon_lnpsi  ;

	    source2_Sij = - qq.derive_con(flat).derive_con(flat) / psi2 * 0.5 ;
	    source2_Sij.inc_dzpuis(1) ;

	    source3_Sij = 4. * nnn * dcon_lnpsi * dcon_logn_auto ;

	    source4_Sij = 0.1666666666666667 * qq.laplacian()/psi2*flat.con() ;

	    source5_Sij = - ( nnn * contract(dcov_lnpsi, 0, 
						  dcon_logn_auto , 0)
	               + nnn * contract(dcov_lnpsi_auto, 0, dcon_lnpsi, 0))
		             * 2.666666666666667 * flat.con() ;

	    source6_Sij = 2. * nnn * contract(contract(gtilde_cov, 0, 
				 tkij_auto + tkij_comp, 1), 0, tkij_auto, 1) ;


	    source7_Sij = - 2. * qpig * nnn * ( psi4 * stress_euler 
			      - 0.33333333333333333 * s_euler * flat.con() ) ; 


	    source_1.change_triad(mp.get_bvect_cart()) ;
	    source_2.change_triad(mp.get_bvect_cart()) ;
	    source_3.change_triad(mp.get_bvect_cart()) ;
	    source_4.change_triad(mp.get_bvect_cart()) ;

	    source1_Hij.change_triad(mp.get_bvect_cart()) ;
	    source2_Hij.change_triad(mp.get_bvect_cart()) ;
	    source3_Hij.change_triad(mp.get_bvect_cart()) ;
	    source4_Hij.change_triad(mp.get_bvect_cart()) ;
	    source5_Hij.change_triad(mp.get_bvect_cart()) ;
	    source6_Hij.change_triad(mp.get_bvect_cart()) ;
	    source7_Hij.change_triad(mp.get_bvect_cart()) ;

	    source1_Qij.change_triad(mp.get_bvect_cart()) ;
	    source2_Qij.change_triad(mp.get_bvect_cart()) ;
	    source3_Qij.change_triad(mp.get_bvect_cart()) ;
	    source4_Qij.change_triad(mp.get_bvect_cart()) ;
	    source5_Qij.change_triad(mp.get_bvect_cart()) ;
	    source6_Qij.change_triad(mp.get_bvect_cart()) ;
	    source7_Qij.change_triad(mp.get_bvect_cart()) ;
	    source8_Qij.change_triad(mp.get_bvect_cart()) ;
	    source9_Qij.change_triad(mp.get_bvect_cart()) ;
	    source10_Qij.change_triad(mp.get_bvect_cart()) ;
	    source11_Qij.change_triad(mp.get_bvect_cart()) ;

	    source1_Sij.change_triad(mp.get_bvect_cart()) ;
	    source2_Sij.change_triad(mp.get_bvect_cart()) ;
	    source3_Sij.change_triad(mp.get_bvect_cart()) ;
	    source4_Sij.change_triad(mp.get_bvect_cart()) ;
	    source5_Sij.change_triad(mp.get_bvect_cart()) ;
	    source6_Sij.change_triad(mp.get_bvect_cart()) ;
	    source7_Sij.change_triad(mp.get_bvect_cart()) ;


	    hij_auto.change_triad(mp.get_bvect_cart()) ;

	    Coord& sin_theta = mp.sint ;
	    Scalar sint (mp) ;
	    sint = sin_theta ;
	    sint.std_spectral_base() ;

	    Sym_tensor source_dirac (mp, CON, mp.get_bvect_spher()) ;
	    // killing vector / rr
	    Vector killing (mp, CON, mp.get_bvect_spher()) ;
	    killing.set(1) = 0. ;
	    killing.set(2) = 0. ;
	    killing.set(3) = omega*sint ;
	    killing.std_spectral_base() ;
	    source_dirac =  - tkij_auto.derive_lie(killing) ;
	    for(int i=1; i<=3; i++) 
		for(int j=1; j<=i; j++) 
		    source_dirac.set(i,j).mult_r() ;
	    source_dirac.inc_dzpuis() ;

	    for(int i=1; i<=3; i++) 
		for(int j=1; j<=i; j++) {

		    source_Hij = ( source1_Hij(i,j) + source1_Hij(j,i) +
			           source2_Hij(j,i) + source2_Hij(i,j) +
				   source3_Hij(j,i) + source3_Hij(i,j) +
				   source4_Hij(j,i) + source4_Hij(i,j) +
				   source5_Hij(i,j) + source6_Hij(i,j)) / psi4
			        + source7_Hij(i,j) ;


		    source_Qij = ( source1_Qij(i,j) + source2_Qij(i,j) +
				   source3_Qij(j,i) + source4_Qij(i,j) +
				   source5_Qij(i,j) + source6_Qij(i,j) +
				   source6_Qij(j,i) + source7_Qij(i,j) +
				   source8_Qij(i,j) + source9_Qij(i,j) +
				   source10_Qij(i,j) + source11_Qij(i,j)) 
			      / psi4 ;
		   

		    source_Sij = ( source1_Sij(i,j) + source2_Sij(i,j) +
				   source3_Sij(i,j) + source3_Sij(j,i) +
				   source4_Sij(i,j) + source5_Sij(i,j)) / psi4
			+ source6_Sij(i,j) + source7_Sij(i,j) ;
			    

		    source_tot_hij = source_1(i,j) + source_1(j,i) 
			+ source_2(i,j) - 2 * psi4 / nnn * (
			    source_3(i,j) + source_4(i,j) +
			    source_Hij + source_Qij + source_Sij ) 
			+0* source_dirac(i,j) ;
			

//		    source_tot_hij.annule(nz-1, nz-1) ;

		    source_tot_hij.filtre(4) ;

		    

		    cout << "source1 de Sij" << endl << norme(source1_Sij(i,j)/ psi4/(nr*nt*np)) << endl ;
		    cout << "source2 de Sij" << endl << norme(source2_Sij(i,j)/ psi4/(nr*nt*np)) << endl ;
		    cout << "source3 de Sij" << endl << norme(2*source3_Sij(i,j)/ psi4/(nr*nt*np)) << endl ;
		    cout << "source4 de Sij" << endl << norme(source4_Sij(i,j)/ psi4/(nr*nt*np)) << endl ;
		    cout << "source5 de Sij" << endl << norme(source5_Sij(i,j)/ psi4/(nr*nt*np)) << endl ;
		    cout << "source6 de Sij" << endl << norme(source6_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source7 de Sij" << endl << norme(source7_Sij(i,j)/(nr*nt*np)) << endl ;
		    cout << "source1" << endl << norme(source_1(i,j)/(nr*nt*np)) << endl ;
		    cout << "source2" << endl << norme(source_2(i,j)/(nr*nt*np)) << endl ;
		    cout << "source3" << endl << norme(source_3(i,j)/(nr*nt*np)) << endl ;
		    cout << "source4" << endl << norme(source_4(i,j)/(nr*nt*np)) << endl ;
		    cout << "source_Hij" << endl << norme(source_Hij/(nr*nt*np)) << endl ;
		    cout << "source_Qij" << endl << norme(source_Qij/(nr*nt*np)) << endl ;
		    cout << "source_Sij" << endl << norme(source_Sij/(nr*nt*np)) << endl ;
		    cout << "source_tot" << endl << norme(source_tot_hij/(nr*nt*np)) << endl ;
		    

		    
		    // Resolution of the Poisson equations and
		    // Check: has the Poisson equation been correctly solved ?
		    // -----------------------------------------------------

		    if(i==1 && j==1) {
		 
			hij_auto.set(1,1) = source_tot_hij.poisson() ; 

			Tbl tdiff_h11 = diffrel(hij_auto(1,1).laplacian(), source_tot_hij) ;  
			cout << "Relative error in the resolution of the equation for "
			     << "h11_auto : " << endl ; 
			for (int l=0; l<nz; l++) {
			    cout << tdiff_h11(l) << "  " ; 
			}
			cout << endl ;
			diff_h11 = max(abs(tdiff_h11)) ; 
		    }
	       	       
		    if(i==2 && j==1) {

			hij_auto.set(2,1) = source_tot_hij.poisson() ; 
	    
			Tbl tdiff_h21 = diffrel(hij_auto(2,1).laplacian(), source_tot_hij) ;
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
		 
			hij_auto.set(3,1) = source_tot_hij.poisson() ; 

			Tbl tdiff_h31 = diffrel(hij_auto(3,1).laplacian(), source_tot_hij) ;
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
		 
			hij_auto.set(2,2) = source_tot_hij.poisson() ; 
			
			Tbl tdiff_h22 = diffrel(hij_auto(2,2).laplacian(), source_tot_hij) ;
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
		 
			hij_auto.set(3,2) = source_tot_hij.poisson() ; 
	
			Tbl tdiff_h32 = diffrel(hij_auto(3,2).laplacian(), source_tot_hij) ;
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
		 
			hij_auto.set(3,3) = source_tot_hij.poisson() ; 

			Tbl tdiff_h33 = diffrel(hij_auto(3,3).laplacian(), source_tot_hij) ;
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

	    hij_auto.change_triad(mp.get_bvect_spher()) ;

	    cout << "Tenseur hij auto in spherical coordinates" << endl ;
	    for (int i=1; i<=3; i++)
		for (int j=1; j<=i; j++) {
		    cout << "  Comp. " << i << " " << j << " :  " ;
		    for (int l=0; l<nz; l++){
			cout << norme(hij_auto(i,j)/(nr*nt*np))(l) << " " ;
		    }
		    cout << endl ;
		}
	    cout << endl ;

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





    
    
