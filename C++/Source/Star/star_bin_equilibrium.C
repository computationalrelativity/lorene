
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
#include "star.h"
#include "param.h"
#include "graphique.h"
#include "utilitaires.h"
#include "vector.h"
#include "nbr_spx.h"


void Star_bin::equilibrium(double ent_c, int mermax, int mermax_potvit, 
			   int mermax_poisson, double relax_poisson, 
			   double relax_potvit, double thres_adapt,
			   const Tbl& fact_resize, Tbl& diff) {


    // Fundamental constants and units
    // -------------------------------
#include "unites.h"	    
    // To avoid some compilation warnings
    if (ent_c < 0) {
	cout << f_unit << msol << km << mevpfm3 << endl ; 
    }    
    
    // Initializations
    // ---------------
    
    const Mg3d* mg = mp.get_mg() ; 
    int nz = mg->get_nzone() ;	    // total number of domains
    
    // The following is required to initialize mp_prev as a Map_et:
//    Map_et& mp_et = dynamic_cast<Map_et&>(mp) ; 
    
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
    double& diff_qq = diff.set(3) ; 
    double& diff_shift_x = diff.set(4) ; 
    double& diff_shift_y = diff.set(5) ; 
    double& diff_shift_z = diff.set(6) ; 
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
 

    // Parameters for the function Tensor::poisson for logn_auto
    // ---------------------------------------------------------
	     
    Cmp ssjm1logn (ssjm1_logn) ;
    Cmp ssjm1qq (ssjm1_qq) ;
    Cmp ssjm1h11 (ssjm1_h11) ;
    Cmp ssjm1h21 (ssjm1_h21) ;
    Cmp ssjm1h31 (ssjm1_h31) ;
    Cmp ssjm1h22 (ssjm1_h22) ;
    Cmp ssjm1h32 (ssjm1_h32) ;
    Cmp ssjm1h33 (ssjm1_h33) ;

    ssjm1logn.set_etat_qcq() ;
    ssjm1qq.set_etat_qcq() ;
    ssjm1h11.set_etat_qcq() ;
    ssjm1h21.set_etat_qcq() ;
    ssjm1h31.set_etat_qcq() ;
    ssjm1h22.set_etat_qcq() ;
    ssjm1h32.set_etat_qcq() ;
    ssjm1h33.set_etat_qcq() ;
 
    double precis_poisson = 1.e-16 ;     
    Param par_poisson1 ; 

    // Parameters for the function Scalar::poisson for logn
    // ---------------------------------------------------------------
 
    par_poisson1.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson1.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson1.add_double(precis_poisson, 1) ; // required precision
    par_poisson1.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_poisson1.add_cmp_mod( ssjm1logn ) ; 
    
    // Parameters for the function Scalar::poisson for qq_auto
    // ---------------------------------------------------------------
    
    Param par_poisson2 ; 
    
    par_poisson2.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson2.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson2.add_double(precis_poisson, 1) ; // required precision
    par_poisson2.add_int_mod(niter, 0) ; // number of iterations actually used -
    par_poisson2.add_cmp_mod( ssjm1qq ) ; 
 

    // Parameters for the function Vector::poisson for shift
    // ---------------------------------------------------------------
    
    Param par_shift ; 
    
    par_shift.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_shift.add_double(relax_poisson,  0) ; // relaxation parameter
    par_shift.add_double(precis_poisson, 1) ; // required precision
    par_shift.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_shift.add_tensor_mod(ssjm1_phi, 0) ; 
    par_shift.add_tensor_mod(ssjm1_khi, 1) ; 
    par_shift.add_tensor_mod(ssjm1_mu, 2) ; 
    
    // Parameters for the function Scalar::poisson for h11_auto
    // -------------------------------------------------------------
    
    Param par_poisson3 ; 
    
    par_poisson3.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson3.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson3.add_double(precis_poisson, 1) ; // required precision
    par_poisson3.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_poisson3.add_cmp_mod( ssjm1h11 ) ; 
    
    // Parameters for the function Scalar::poisson for h21_auto
    // -------------------------------------------------------------
    
    Param par_poisson4 ; 
    
    par_poisson4.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson4.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson4.add_double(precis_poisson, 1) ; // required precision
    par_poisson4.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_poisson4.add_cmp_mod( ssjm1h21 ) ; 
    
    // Parameters for the function Scalar::poisson for h31_auto
    // -------------------------------------------------------------
    
    Param par_poisson5 ; 
    
    par_poisson5.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson5.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson5.add_double(precis_poisson, 1) ; // required precision
    par_poisson5.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_poisson5.add_cmp_mod( ssjm1h31 ) ; 
    
    // Parameters for the function Scalar::poisson for h22_auto
    // -------------------------------------------------------------
    
    Param par_poisson6 ; 
    
    par_poisson6.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson6.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson6.add_double(precis_poisson, 1) ; // required precision
    par_poisson6.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_poisson6.add_cmp_mod( ssjm1h22 ) ; 
    
    // Parameters for the function Scalar::poisson for h32_auto
    // -------------------------------------------------------------
    
    Param par_poisson7 ; 
    
    par_poisson7.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson7.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson7.add_double(precis_poisson, 1) ; // required precision
    par_poisson7.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_poisson7.add_cmp_mod( ssjm1h32 ) ; 
    
    // Parameters for the function Scalar::poisson for h33_auto
    // -------------------------------------------------------------
    
    Param par_poisson8 ; 
    
    par_poisson8.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson8.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson8.add_double(precis_poisson, 1) ; // required precision
    par_poisson8.add_int_mod(niter, 0) ; // number of iterations actually used 
    par_poisson8.add_cmp_mod( ssjm1h33 ) ; 
    

    // External potential
    // See Eq (99) from Gourgoulhon et al. (2001)
    // ------------------
    

    Scalar pot_ext = logn_comp + pot_centri + loggam ;
    
    Scalar ent_jm1 = ent ;	// Enthalpy at previous step
    
    Scalar source_tot(mp) ; // source term in the equation for hij_auto, 
                            // logn_auto and qq_auto
			    
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
/*
	if (irrotational) {
	    diff_vel_pot = velocity_potential(mermax_potvit, precis_poisson, 
					      relax_potvit) ; 
	    
	}
*/
	diff_vel_pot = 0. ; // to avoid the warning 

	//-----------------------------------------------------
	// Computation of the new radial scale
	//--------------------------------------------------

	// alpha_r (r = alpha_r r') is determined so that the enthalpy
	// takes the requested value ent_b at the stellar surface
	
	// Values at the center of the star:
	double logn_auto_c  = logn_auto.val_grid_point(0, 0, 0, 0) ; 
	double pot_ext_c  = pot_ext.val_grid_point(0, 0, 0, 0) ; 


	cout << "ok2" << endl ;

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
	// First integral	--> enthalpy in all space
	// See Eq (98) from Gourgoulhon et al. (2001)
	//-------------------------------------------

	ent = (ent_c + logn_auto_c + pot_ext_c) - logn_auto - pot_ext ;

	(ent.set_spectral_va()).smooth(nzet, ent.set_spectral_va()) ;


/*
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
	
	int n_resize ;
	if (nz > 3) {
	    n_resize = nz - 3 ;
	}
	else {
	    n_resize = nzet ;
	}

	double rr_in = mp.val_r(nzet,-1., M_PI/2, 0.) ; 
	double rr_out = mp.val_r(n_resize,1., M_PI/2, 0.) ; 

	mp.resize(n_resize, rr_in/rr_out * fact_resize(0)) ; 

	//----------------------------------------------------
	// Computation of the enthalpy at the new grid points
	//----------------------------------------------------
	
	mp_prev.homothetie(alpha_r) ; 
	
	Cmp ent_cmp2 (ent) ;
	mp.reevaluate_symy(&mp_prev, nzet+1, ent_cmp2) ; 
	ent = ent_cmp ;

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
*/
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

	
	// Derivatives of N and logN
	//--------------------------

	const Vector& dcov_logn_auto = logn_auto.derive_cov(flat) ;
	const Vector& dcon_logn_auto = logn_auto.derive_con(flat) ;
	
	const Vector& dcov_nnn = nnn.derive_cov(flat) ;
	const Tensor& dcovdcov_nnn = dcov_nnn.derive_cov(flat) ;

	// Derivatives of Q and shift 
	//---------------------------

	const Vector& dcov_qq = qq.derive_cov(flat) ;
	
	const Tensor& dcovdcov_qq = dcov_qq.derive_cov(flat) ;
	const Tensor& dcondcov_qq = dcov_qq.derive_con(flat) ;

	const Vector& dcon_qq_auto = qq_auto.derive_con(flat) ;
	const Tensor& dcondcon_qq_auto = dcon_qq_auto.derive_con(flat) ;

	const Tensor& dcov_shift = shift.derive_cov(flat) ;
	const Tensor& dcovdcov_shift = dcov_shift.derive_cov(flat) ;

	Scalar psi2 (qq / nnn) ;

	Scalar laplacian_psi6 (mp) ;
	laplacian_psi6.set_etat_qcq() ;
	laplacian_psi6 = (psi4 * psi2).laplacian() ;
	Tensor dcondcov_psi6 = ((psi4 * psi2).derive_cov(flat))
	                                            .derive_con(flat) ;
	dcondcov_psi6.inc_dzpuis() ;


	// Derivatives of hij, gtilde... 
	//------------------------------

	const Tensor& dcov_hij = hij.derive_cov(flat) ;
	const Tensor& dcov_hij_auto = hij_auto.derive_cov(flat) ;
	const Tensor& dcon_hij_auto = hij_auto.derive_con(flat) ;

	const Tensor& dcovdcov_hij = dcov_hij.derive_cov(flat) ;

	Sym_tensor gtilde_cov = gtilde.cov() ;
	Sym_tensor gtilde_con = gtilde.con() ;
	gtilde_cov.std_spectral_base() ;
	gtilde_con.std_spectral_base() ;
	const Tensor& dcov_gtildeij = gtilde_cov.derive_cov(flat) ;

//	Connection gamijk (gtilde, flat) ;
//	const Tensor& deltaijk = gamijk.get_delta() ;

	const Tensor& dcov_tkij_auto = tkij_auto.derive_cov(flat) ;
   
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

	source3 = - contract_desal(dcov_logn_auto, 0, dcon_logn, 0) 
	-2. * contract_desal(dcov_lnpsi, 0, dcon_logn_auto, 0) ;
	
	source4 = - contract(contract(hij_auto, 0, dcovdcov_nnn/nnn, 0), 
			     0, 1) ;

	source4.inc_dzpuis(1) ;

	source5 = 2. * contract(contract(hij_auto, 0, dcov_lnpsi, 0), 
			       0, dcov_logn, 0) ;


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
    
	//   source.annule(nz-1) ;
 

//	logn_auto = source_tot.poisson() ; 
	source_tot.poisson(par_poisson1, logn_auto) ; 

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

	cout << "Rel error: " << norme(logn_auto.laplacian()- source_tot)
	    / max( max( source_tot )) << endl ; 

    
	//--------------------------------------------------------
	// Poisson equation for qq_auto
	//--------------------------------------------------------

	// Source
	//--------

	source1 = qpig * psi2 % psi4 % nnn % s_euler ;

	source2 = 3. * psi2 % psi4 % nnn % (kcar_auto + kcar_comp)/ 4. ;
	
	source3 = 2. * qq_auto % contract_desal(dcov_lnpsi, 0, dcon_lnpsi, 0) ;

	source4 = 2. * qq_auto % contract_desal(dcov_lnpsi, 0, dcon_logn, 0) ;

	source5 = qq_auto * contract(contract(contract(contract(gtilde_con, 0, 
		   dcov_hij, 2), 0, dcov_gtildeij, 2), 0, 2), 0, 1)/16. ;

	source6 = - qq_auto*contract(contract(contract(contract(gtilde_con, 0, 
		    dcov_hij, 2), 0, dcov_gtildeij, 1), 0, 2), 0, 1)/8. ;

	source7 = qq_auto * 2. * contract(contract(hij, 0, dcov_lnpsi, 0)
			       , 0, dcov_lnpsi, 0) ;

	source8 = 2 * qq_auto * contract(contract(hij, 0, dcov_lnpsi, 0)
				     , 0, dcov_logn, 0) ;
	
	source9 = - contract(contract(hij_auto, 0, dcovdcov_qq, 1), 0, 1) ;
	source9.inc_dzpuis(1) ;

	source_tot = source1 + source2 + source3 + source4 + source5 + 
	             source6 + source7 + source8 + source9 ;

  	
	cout << "moyenne de la source 1 pour qq_auto" << endl ;
	cout <<  norme(source1/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 2 pour qq_auto" << endl ;
	cout <<  norme(source2/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 3 pour qq_auto" << endl ;
	cout <<  norme(source3/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 4 pour qq_auto" << endl ;
	cout <<  norme(source4/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 5 pour qq_auto" << endl ;
	cout <<  norme(source5/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 6 pour qq_auto" << endl ;
	cout <<  norme(source6/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 7 pour qq_auto" << endl ;
	cout <<  norme(source7/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 8 pour qq_auto" << endl ;
	cout <<  norme(source8/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 9 pour qq_auto" << endl ;
	cout <<  norme(source9/(nr*nt*np)) << endl ;
	cout << "moyenne de la source pour qq_auto" << endl ;
	cout <<  norme(source_tot/(nr*nt*np)) << endl ;
	

	// Resolution of the Poisson equation 
	// ----------------------------------

	//source.set().filtre(4) ;
	//   source.annule(nz-1) ;
 
//	qq_auto = source_tot.poisson() ; 
	source_tot.poisson(par_poisson2, qq_auto) ; 

	qq_auto = qq_auto + decouple ;

	cout << "qq_auto" << endl << norme(qq_auto/(nr*nt*np)) << endl ;
	
         /*
	  des_profile(qq_auto(), 0, 20, 0, 0) ;
	  des_coef_xi(qq_auto().va, 0, 0, 0) ;
	  des_coef_xi(qq_auto().va, 1, 0, 0) ;
	  des_coef_xi(qq_auto().va, 2, 0, 0) ;
	*/
	// Check: has the Poisson equation been correctly solved 
	// -----------------------------------------------------
    
	Tbl tdiff_qq = diffrel(qq_auto.laplacian(), source_tot) ;
	cout << 
	    "Relative error in the resolution of the equation for qq : "
	     << endl ; 
	for (int l=0; l<nz; l++) {
	    cout << tdiff_qq(l) << "  " ; 
	}
	cout << endl ;
	diff_qq = max(abs(tdiff_qq)) ; 

	cout << "Rel error: " << norme(qq_auto.laplacian()- source_tot)
	    / max( max( source_tot ) ) << endl ; 

	//--------------------------------------------------------
	// Vector Poisson equation for shift_auto
	//--------------------------------------------------------

	// Source
	//--------


	Vector source1_shift(mp, CON, mp.get_bvect_spher()) ;
	Vector source2_shift(mp, CON, mp.get_bvect_spher()) ;
	Vector source3_shift(mp, CON, mp.get_bvect_spher()) ;
	Vector source4_shift(mp, CON, mp.get_bvect_spher()) ;
	Vector source5_shift(mp, CON, mp.get_bvect_spher()) ;

	cout << "-4*qpig*nnn" << endl << norme(-4*qpig*nnn) << endl ;
	cout << "psi4" << endl << norme(psi4) << endl ;
	cout << "ener_euler" << endl << norme(ener_euler) << endl ;
	cout << "press" << endl << norme(press) << endl ;
	u_euler.change_triad(mp.get_bvect_cart()) ;
	cout << "u_euler(1)" << endl << norme(u_euler(1)) << endl ;
	cout << "u_euler(2)" << endl << norme(u_euler(2)) << endl ;
	cout << "u_euler(3)" << endl << norme(u_euler(3)) << endl ;


	u_euler.change_triad(mp.get_bvect_spher()) ;

	// operator% only defined for scalar % scalar
	for (int i=1; i<=3; i++){
	    source1_shift.set(i) = (4.*qpig) * nnn % psi4
	                   %(ener_euler + press) * u_euler(i) ;
	}

	source2_shift = 2. * nnn * contract(tkij_auto, 1, 
					    dcov_logn - 6 * dcov_lnpsi, 0) ;

	source3_shift.set_etat_zero() ; // = - 2. * nnn * contract(contract(tkij_auto, 0, 
	// deltaijk, 1), 0, 2) ;

	source4_shift = - contract(contract(hij_auto, 0, dcovdcov_shift, 2)
				   , 0, 2) ;

	source4_shift.inc_dzpuis(1) ;

	source5_shift = - contract(contract(hij_auto, 1, dcovdcov_shift, 2)
				   , 1, 2) / 3. ;
	
	source5_shift.inc_dzpuis(1) ;

	source_shift = source1_shift + source2_shift + source3_shift 
	    + source4_shift + source5_shift ;

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
	
	// Resolution of the Poisson equation 
	// ----------------------------------

	// Filter for the source of shift vector
    
/*	for (int i=1; i<=3; i++) {
  	    if (source_shift(i).get_etat() != ETATZERO)
		source_shift.set(i).filtre(4) ;
	}
*/  
	for (int i=1; i<=3; i++) {
	    if(source_shift(i).dz_nonzero()) {
		assert( source_shift(i).get_dzpuis() == 4 ) ; 
	    }
	    else{
		(source_shift.set(i)).set_dzpuis(4) ; 
	    }
	}

	double lambda = 1./3. ;
	shift_auto = source_shift.poisson(lambda, flat, 2) ; 

	// Check: has the equation for shift_auto been correctly solved ?
	// --------------------------------------------------------------
	    
	// Divergence of shift_auto : 
	Scalar divn = shift_auto.divergence(flat) ;

	// Grad(div) : 
        Vector graddivn = divn.derive_cov(flat) ;
	graddivn.change_triad(mp.get_bvect_cart()) ;
	graddivn.inc_dzpuis() ;

	// Full operator :
	Vector lap_shift(mp, CON, mp.get_bvect_cart() ) ;
	shift_auto.change_triad(mp.get_bvect_cart()) ;
	for (int i=1; i<=3; i++) {
	    lap_shift.set(i) = shift_auto(i).laplacian() + lambda*graddivn(i);
	}

	cout << "shift_auto(1)" << endl << norme(shift_auto(1)/(nr*nt*np)) << endl ;
	cout << "shift_auto(2)" << endl << norme(shift_auto(2)/(nr*nt*np)) << endl ;
	cout << "shift_auto(3)" << endl << norme(shift_auto(3)/(nr*nt*np)) << endl ;


	shift_auto.change_triad(mp.get_bvect_spher()) ;
	lap_shift.change_triad(mp.get_bvect_spher()) ;
	graddivn.change_triad(mp.get_bvect_spher()) ;

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



	if (!conf_flat){
	   
	    //--------------------------------------------------------
	    // Poisson equation for hij
	    //--------------------------------------------------------
	 
	 
	    // Declaration of all sources 
	    //---------------------------

	    Scalar sol_poisson(mp) ;
	    Scalar source_Hij(mp) ;
	    Scalar source_Qij(mp) ;
	    Scalar source_Sij(mp) ;
	    Scalar source_tot_hij(mp) ;


	    Tensor source1 (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source2 (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source3 (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source4 (mp, 2, CON, mp.get_bvect_spher()) ;
	    Tensor source5 (mp, 2, CON, mp.get_bvect_spher()) ;


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

	    source1 = contract(dcon_hij_auto, 1, dcov_logn + 2*dcov_lnpsi, 0) ;

	    source2 = contract(dcon_hij_auto, 2, dcov_logn + 2*dcov_lnpsi, 0) ;

	    source3 = contract(shift, 0, dcov_tkij_auto, 2) ;
	    source3.inc_dzpuis(1) ;

	    source4 = - contract(tkij_auto, 1, dcov_shift, 1) ;

	    source5 = - 2. * contract(dcov_shift, 0, 1) * tkij_auto / 3. ;

	    
	    // Source terms for Hij
	    //---------------------

	    source1_Hij = 8.*nnn*contract(hij_auto,1,dcov_lnpsi,0)*dcon_lnpsi ;

//	    source2_Hij = - contract(hij_auto, 1 , dcondcov_qq, 0) / psi2 ;
//	    source2_Hij.inc_dzpuis(1) ;

	    source3_Hij = 4.*nnn*contract(hij_auto,1,dcov_logn,0)*dcon_lnpsi ; 

	    source4_Hij = 4.*nnn*contract(hij_auto,1,dcov_lnpsi,0)*dcon_logn ;

	    Tensor temp = contract(contract(hij_auto, 0
					    , dcovdcov_qq, 1), 0, 1) / psi2 ;
	    temp.inc_dzpuis(1) ;

	    source5_Hij = ( temp 
			   - 8. * nnn * contract(contract(hij_auto,0
					  , dcov_lnpsi, 0), 0, dcov_logn, 0)
			   - 8. * nnn * contract(contract(hij_auto,0
					  , dcov_lnpsi, 0), 0, dcov_lnpsi, 0))
		* flat.con() / 3. ; 

	    source6_Hij = ( qq.laplacian() / psi2
			   - 8. * nnn * contract(dcov_lnpsi, 0, dcon_logn, 0)
			   - 8. * nnn * contract(dcov_lnpsi, 0, dcon_lnpsi, 0))
		* hij_auto / 3. ;

	    source7_Hij = 2. * qpig * nnn * s_euler * hij_auto / 3. ;



	    // Source_terms for Qij
	    //---------------------

	    
	    source1_Qij = nnn * contract(contract(hij_auto, 0, 
						  dcovdcov_hij, 3), 0, 3) /2. ;
	    source1_Qij.inc_dzpuis(1) ;


	    source2_Qij = nnn / 2. * ( - contract(contract(dcov_hij_auto, 1, 
							   dcov_hij, 2), 1, 3)
	     - contract(contract(contract(contract(gtilde_cov, 0, 
		  dcov_hij_auto, 1), 0, dcov_hij, 1), 1, gtilde_con, 0), 2, 3)
	     + contract(contract(contract(contract(gtilde_con, 1, 
		  dcov_hij_auto, 2), 2, dcov_hij, 2), 1, gtilde_cov, 0), 2, 3)
	     + contract(contract(contract(contract(gtilde_con, 1, 
	        dcov_gtildeij, 2), 2, dcov_hij_auto, 1), 1, 2), 1
			, gtilde_con, 1) ) ;

	    source3_Qij = contract(contract(contract(contract(gtilde_con, 1, 
	         dcov_hij_auto, 2), 1, dcov_hij, 2), 1, gtilde_cov, 1), 2, 3) ;

	    source4_Qij = 16 * contract(contract(hij_auto, 1, 
				        dcov_lnpsi, 0)*hij, 2, dcov_lnpsi, 0) ;
	    
	    source5_Qij = - contract(contract(hij_auto, 1, 
					 dcovdcov_qq, 1), 1, hij, 1) / psi2 ;
	    source5_Qij.inc_dzpuis(1) ;


	    source6_Qij = - contract(contract(hij_auto, 1, dcov_hij, 2)
				     , 2, dcov_qq, 0) / (2.*psi2) ;


	    source7_Qij = - contract(contract(hij_auto, 1, dcov_hij, 2), 
				     0, dcov_qq, 0) / (2.*psi2) ;


	    source8_Qij = 4. * nnn * contract(contract(hij_auto, 1, 
				       dcov_lnpsi, 0)*hij, 2, dcov_logn, 0)

		        + 4. * nnn * contract(contract(hij_auto, 1, 
				       dcov_logn, 0)*hij, 2, dcov_lnpsi, 0) ;

	    source9_Qij = nnn / 6. * (contract(contract(contract(contract(
          gtilde_con, 0, dcov_hij_auto, 2), 0, dcov_gtildeij, 1), 0, 2), 0, 1)
		- contract(contract(contract(contract(gtilde_con, 0, 
		   dcov_hij_auto, 2), 0, dcov_gtildeij, 2), 0, 2), 0, 1) / 2.)
		* gtilde_con ;

	    
	    source10_Qij = - 8./3. * nnn * ( contract(contract(hij, 0, 
			        dcov_lnpsi, 0), 0, dcov_lnpsi, 0) * hij_auto 
					   + contract(contract(hij, 0, 
				dcov_lnpsi, 0), 0, dcov_logn, 0) * hij_auto ) ;


	    source11_Qij = contract(contract(hij, 0, dcovdcov_qq, 0), 0, 1) 
		* hij_auto / (3.*psi2) ;
	    
	    source11_Qij.inc_dzpuis(1) ;



	    // Source terms for Sij
	    //---------------------

	    
	    source1_Sij = 8. * qq_auto * dcon_lnpsi * dcon_lnpsi / psi2 ;

	    source2_Sij = - dcondcon_qq_auto / psi2 ;
	    source2_Sij.inc_dzpuis(1) ;

	    source3_Sij = 4. * nnn * dcon_lnpsi * dcon_logn_auto ;

	    source4_Sij = qq_auto.laplacian() / (3. * psi2) * flat.con() ;

	    source5_Sij = - 8. * ( nnn * contract(dcov_lnpsi, 0, dcon_logn, 0)
	                          - qq_auto / psi2
				  * contract(dcov_lnpsi, 0, dcon_lnpsi, 0))/3.
		           * flat.con() ;

	    source6_Sij = 2. * nnn * contract(contract(gtilde_cov, 0, 
				 tkij_auto + tkij_auto, 1), 0, tkij_auto, 1) ;


	    source7_Sij = - 2. * qpig * nnn * ( psi4 * stress_euler 
					       - s_euler * flat.con()/3. ) ; 



	    source1.change_triad(mp.get_bvect_cart()) ;
	    source2.change_triad(mp.get_bvect_cart()) ;
	    source3.change_triad(mp.get_bvect_cart()) ;
	    source4.change_triad(mp.get_bvect_cart()) ;
	    source5.change_triad(mp.get_bvect_cart()) ;

	    source1_Hij.change_triad(mp.get_bvect_cart()) ;
//	    source2_Hij.change_triad(mp.get_bvect_cart()) ;
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


	    for(int i=1; i<=3; i++) 
		for(int j=1; j<=i; j++) {

		    source_Hij = ( source1_Hij(i,j) + source1_Hij(j,i) +
				//     source2_Hij(j,i) + source2_Hij(i,j) +
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
				   source4_Sij(i,j) + source5_Sij(i,j) +
				   source6_Sij(i,j) + source7_Sij(i,j) )
			     / psi4 ;  


		    source_tot_hij = source1(i,j) + source1(j,i) - source2(i,j)
			+ source3(i,j) + source4(i,j) + source4(j,i) 
		    + source5(i,j) + source_Hij + source_Qij 
		 	+ source_Sij ;

		    
		    cout << "source1" << endl << norme(source1(i,j)/(nr*nt*np)) << endl ;
		    cout << "source2" << endl << norme(source2(i,j)/(nr*nt*np)) << endl ;
		    cout << "source3" << endl << norme(source3(i,j)/(nr*nt*np)) << endl ;
		    cout << "source4" << endl << norme(source4(i,j)/(nr*nt*np)) << endl ;
		    cout << "source5" << endl << norme(source5(i,j)/(nr*nt*np)) << endl ;
		    cout << "source_Hij" << endl << norme(source_Hij/(nr*nt*np)) << endl ;
		    cout << "source_Qij" << endl << norme(source_Qij/(nr*nt*np)) << endl ;
		    cout << "source_Sij" << endl << norme(source_Sij/(nr*nt*np)) << endl ;
		    cout << "source_tot" << endl << norme(source_tot_hij/(nr*nt*np)) << endl ;
		    



		    // Resolution of the Poisson equations and
		    // Check: has the Poisson equation been correctly solved ?
		    // -----------------------------------------------------

     
		    if(i==1 && j==1) {
		 
			sol_poisson = source_tot_hij.poisson() ; 
			hij_auto.set(1,1) = sol_poisson ;

			hij_auto.set(1,1).set_spectral_va().coef_i() ;

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

			sol_poisson = source_tot_hij.poisson() ; 
			hij_auto.set(2,1) = sol_poisson ;

			hij_auto.set(2,1).set_spectral_va().coef_i() ;
	    
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
		 
			sol_poisson = source_tot_hij.poisson() ; 
			hij_auto.set(3,1) = sol_poisson ;

			hij_auto.set(3,1).set_spectral_va().coef_i() ;

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
		 
			sol_poisson = source_tot_hij.poisson() ; 
			hij_auto.set(2,2) = sol_poisson ;

			hij_auto.set(2,2).set_spectral_va().coef_i() ;

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
		 
			sol_poisson = source_tot_hij.poisson() ; 
			hij_auto.set(3,2) = sol_poisson ;

			hij_auto.set(3,2).set_spectral_va().coef_i() ;

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
		 
			sol_poisson = source_tot_hij.poisson() ; 
			hij_auto.set(3,3) = sol_poisson ;

			hij_auto.set(3,3).set_spectral_va().coef_i() ;

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
      
	    hij_auto.std_spectral_base() ;
	    hij_auto.change_triad(mp.get_bvect_spher()) ;


	    ssjm1_h11 = ssjm1h11 ;
	    ssjm1_h21 = ssjm1h21 ;
	    ssjm1_h31 = ssjm1h31 ;
	    ssjm1_h22 = ssjm1h22 ;
	    ssjm1_h32 = ssjm1h32 ;
	    ssjm1_h33 = ssjm1h33 ;

      
	cout << "h11 auto"<<endl<< norme(hij_auto(1,1)/(nr*nt*np)) 
	     << endl << endl ;
	cout << "h21 auto"<<endl<< norme(hij_auto(2,1)/(nr*nt*np)) 
	     << endl << endl ;
	cout << "h31 auto"<<endl<< norme(hij_auto(3,1)/(nr*nt*np)) 
	     << endl << endl ;
	cout << "h22 auto"<<endl<< norme(hij_auto(2,2)/(nr*nt*np)) 
	     << endl << endl ;
	cout << "h32 auto"<<endl<< norme(hij_auto(3,2)/(nr*nt*np)) 
	     << endl << endl ;
	cout << "h33 auto"<<endl<< norme(hij_auto(3,3)/(nr*nt*np)) 
	     << endl << endl ;

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





    
    
