
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


void Star_bin::equilibrium(double ent_c, int mermax, int mermax_potvit, 
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
    double& diff_qq = diff.set(3) ; 
    double& diff_shift_x = diff.set(4) ; 
    double& diff_shift_y = diff.set(5) ; 
    double& diff_shift_z = diff.set(6) ; 
    double& diff_hij00 = diff.set(7) ; 
    double& diff_hij10 = diff.set(8) ; 
    double& diff_hij20 = diff.set(9) ; 
    double& diff_hij11 = diff.set(10) ; 
    double& diff_hij21 = diff.set(11) ; 
    double& diff_hij22 = diff.set(12) ; 



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
			   
    double precis_poisson = 1.e-16 ;     
   
    // External potential
    // See Eq (99) from Gourgoulhon et al. (2001)
    // ------------------
    

    Scalar pot_ext = logn_comp + pot_centri + loggam ;
    
    Scalar ent_jm1 = ent ;	// Enthalpy at previous step
    
    Scalar source(mp) ;    // source term in the equation for logn_auto,
    // qq_auto 

    Scalar source_tot(mp) ; // source term in the equation for gtildeij
			    
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
	//-----------------------------------------------------
	// Computation of the new radial scale
	//--------------------------------------------------

	// alpha_r (r = alpha_r r') is determined so that the enthalpy
	// takes the requested value ent_b at the stellar surface
	
	// Values at the center of the star:
	double logn_auto_c  = logn_auto.point(0, 0, 0, 0) ; 
	double pot_ext_c  = pot_ext.point(0, 0, 0, 0) ; 

	// Search for the reference point (theta_*, phi_*) [notation of
	//  Bonazzola, Gourgoulhon & Marck PRD 58, 104020 (1998)]
	//  at the surface of the star
	// ------------------------------------------------------------
	double alpha_r2 = 0 ; 
	for (int k=0; k<mg->get_np(l_b); k++) {
	    for (int j=0; j<mg->get_nt(l_b); j++) {
		
		double pot_ext_b  = pot_ext.point(l_b, k, j, i_b) ; 
		double logn_auto_b  = logn_auto.point(l_b, k, j, i_b) ; 
		

		// See Eq (100) from Gourgoulhon et al. (2001)
		double alpha_r2_jk = ( ent_c - ent_b + pot_ext_c - pot_ext_b) /
 
		    ( logn_auto_b - logn_auto_c ) ;

		//		cout << "k, j, alpha_r2_jk : " << k << "  " << j << "  " 
		//		     << alpha_r2_jk << endl ; 
		  
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
	logn_auto_c  = logn_auto.point(0, 0, 0, 0) ;

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
	    ent_limit.set(l) = ent.point(l, k_b, j_b, i_b) ; 
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
		double xx = fabs( ent.point(l_b, k, j, i_b) ) ;
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

	
	// Derivatives of N and logN
	//--------------------------

	const Vector& dcov_logn_auto = logn_auto.derive_cov(flat) ;
	const Vector& dcon_logn_auto = logn_auto.derive_con(flat) ;
	
	Vector dcov_nnn = nnn.derive_cov(flat) ;
	dcov_nnn.dec_dzpuis(2) ;
	Tensor dcovdcov_nnn = dcov_nnn.derive_cov(flat) ;
	dcov_nnn.inc_dzpuis(2) ;
	dcovdcov_nnn.inc_dzpuis(2) ;


	// Derivatives of Q and shift
	//---------------------------

	const Vector& dcov_qq = qq.derive_cov(flat) ;
	
	Vector dcov_qq2 = qq.derive_cov(flat) ;
	dcov_qq2.dec_dzpuis(2) ;
	Tensor dcovdcov_qq = dcov_qq2.derive_cov(flat) ;
	Tensor dcondcov_qq = dcov_qq2.derive_con(flat) ;
	dcov_qq2.inc_dzpuis(2) ;
	dcovdcov_qq.inc_dzpuis(2) ;
	dcondcov_qq.inc_dzpuis(2) ;
	
	Vector dcon_qq = qq.derive_con(flat) ;
	dcon_qq.dec_dzpuis(2) ;
        Tensor dcondcon_qq = dcon_qq.derive_con(flat) ;
	dcon_qq.inc_dzpuis(2) ;
	dcondcon_qq.inc_dzpuis(2) ;

	const Scalar divshift = shift.derive_cov(flat).scontract(0, 1) ;
	const Tensor dcov_divshift = divshift.derive_cov(flat) ;


	// Derivatives oh hij, gtilde... 
	//------------------------------

	const Tensor& dcov_hij = hij.derive_cov(flat) ;
	const Tensor& dcon_hij = hij.derive_con(flat) ;
	const Tensor& dcov_hij_auto = hij_auto.derive_cov(flat) ;

	Tensor dcov_hij2 =  hij.derive_cov(flat) ;
	dcov_hij2.dec_dzpuis(2) ;
	Tensor dcovdcov_hij = dcov_hij2.derive_cov(flat) ;
	dcov_hij2.inc_dzpuis(2) ;
	dcovdcov_hij.inc_dzpuis(2) ;

	const Tensor& dcov_gtildeij = (gtilde.cov()).derive_cov(flat) ;

   
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
	Scalar source10(mp) ;
 
 
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
	cout <<  norme(source/(nr*nt*np)) << endl ;
	
	
	source.std_spectral_base() ;
	/*
	  des_profile(source1(), 0, 20, 0, 0) ;
	  des_profile(source2(), 0, 20, 0, 0) ;
	  des_profile(source3(), 0, 20, 0, 0) ;
	  des_profile(source4(), 0, 20, 0, 0) ;
	  des_profile(source(), 0, 20, 0, 0) ;
	*/

	// Resolution of the Poisson equation 
	// ----------------------------------
    

	//source.filtre(4) ;
	//   source.annule(nz-1) ;
 

	logn_auto = source.poisson() ; 

	cout << "logn_auto" << endl << norme(logn_auto/(nr*nt*np)) << endl ;
  
	// Check: has the Poisson equation been correctly solved ?
	// -----------------------------------------------------

	Tbl tdiff_logn = diffrel(logn_auto.laplacian(), source) ;
	cout << 
	    "Relative error in the resolution of the equation for logn_auto : "
	     << endl ; 
	for (int l=0; l<nz; l++) {
	    cout << tdiff_logn(l) << "  " ; 
	}
	cout << endl ;
	diff_logn = max(abs(tdiff_logn)) ; 

	cout << "Rel error: " << norme(logn_auto.laplacian()- source)
	    / max( max( source )) << endl ; 

    
	//--------------------------------------------------------
	// Poisson equation for qq_auto
	//--------------------------------------------------------

	// Source
	//--------


  	
	cout << "moyenne de la source 1 pour qq_auto" << endl ;
	cout <<  norme(source1/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 2 pour qq_auto" << endl ;
	cout <<  norme(source2/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 3 pour qq_auto" << endl ;
	cout <<  norme(source3/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 4 pour qq_auto" << endl ;
	cout <<  norme(source4/(nr*nt*np)) << endl ;
	cout << "moyenne de la source pour qq_auto" << endl ;
	cout <<  norme(source/(nr*nt*np)) << endl ;
	

	source.std_spectral_base() ;

	// Resolution of the Poisson equation 
	// ----------------------------------

	//source.set().filtre(4) ;
	//   source.annule(nz-1) ;
 
	qq_auto = source.poisson() ; 
	/*
	  des_profile(qq_auto(), 0, 20, 0, 0) ;
	  des_coef_xi(qq_auto().va, 0, 0, 0) ;
	  des_coef_xi(qq_auto().va, 1, 0, 0) ;
	  des_coef_xi(qq_auto().va, 2, 0, 0) ;
	*/
	// Check: has the Poisson equation been correctly solved 
	// -----------------------------------------------------
    
	Tbl tdiff_qq = diffrel(qq_auto.laplacian(), source) ;
	cout << 
	    "Relative error in the resolution of the equation for qq : "
	     << endl ; 
	for (int l=0; l<nz; l++) {
	    cout << tdiff_qq(l) << "  " ; 
	}
	cout << endl ;
	diff_qq = max(abs(tdiff_qq)) ; 

	cout << "Rel error: " << norme(qq_auto.laplacian()- source)
	    / max( max( source ) ) << endl ; 

	qq_auto = qq_auto + decouple  ;

	cout << "qq_auto" << endl << norme(qq_auto/(nr*nt*np)) << endl ;

	//--------------------------------------------------------
	// Vector Poisson equation for shift_auto 
	//--------------------------------------------------------

	// Source
	//--------


	Vector source1_shift(mp, CON, mp.get_bvect_spher()) ;
	Vector source2_shift(mp, CON, mp.get_bvect_spher()) ;
	Vector source3_shift(mp, CON, mp.get_bvect_spher()) ;
	Vector source4_shift(mp, CON, mp.get_bvect_spher()) ;
 

       
	cout << "moyenne de la source 1 pour shift_auto" << endl ;
	cout <<  norme(source1_shift(0)/(nr*nt*np)) << endl ;
	cout <<  norme(source1_shift(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source1_shift(2)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 2 pour shift_auto" << endl ;
	cout <<  norme(source2_shift(0)/(nr*nt*np)) << endl ;
	cout <<  norme(source2_shift(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source2_shift(2)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 3 pour shift_auto" << endl ;
	cout <<  norme(source3_shift(0)/(nr*nt*np)) << endl ;
	cout <<  norme(source3_shift(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source3_shift(2)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source 4 pour shift_auto" << endl ;
	cout <<  norme(source4_shift(0)/(nr*nt*np)) << endl ;
	cout <<  norme(source4_shift(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source4_shift(2)/(nr*nt*np)) << endl ;
	cout << "moyenne de la source pour shift_auto" << endl ;
	cout <<  norme(source_shift(0)/(nr*nt*np)) << endl ;
	cout <<  norme(source_shift(1)/(nr*nt*np)) << endl ;
	cout <<  norme(source_shift(2)/(nr*nt*np)) << endl ;
	
	source_shift.std_spectral_base() ; 	
	
	// Resolution of the Poisson equation 
	// ----------------------------------

	// Filter for the source of shift vector
    
	for (int i=1; i<=3; i++) {
  	    if (source_shift(i).get_etat() != ETATZERO)
		source_shift.set(i).filtre(4) ;
	}
    
	for (int i=1; i<=3; i++) {
	    if(source_shift(i).dz_nonzero()) {
		assert( source_shift(i).get_dzpuis() == 4 ) ; 
	    }
	    else{
		(source_shift.set(i)).set_dzpuis(4) ; 
	    }
	}

	double lambda_shift = double(1)/double(3) ; 

	shift_auto = source_shift.poisson(lambda_shift) ;      


	// Check: has the equation for shift_auto been correctly solved ?
	// --------------------------------------------------------------
	// Divergence of shift_auto : 
	Scalar divn = shift_auto.derive_cov(flat).scontract(0, 1) ; 
	divn.dec_dzpuis(2) ;    // dzpuis 2 -> 0
	    
	// Grad(div) : 
        Vector graddivn = divn.derive_cov(flat) ;
	graddivn.inc_dzpuis(2) ;    // dzpuis 2 -> 4
	    
	// Full operator :
	Vector lap_shift(mp, CON, mp.get_bvect_spher() ) ;
	for (int i=1; i<=3; i++) {
	    lap_shift.set(i) = shift_auto(i).laplacian() 
		+ lambda_shift * graddivn(i) ;
	}

	Tbl tdiff_shift_x = diffrel(lap_shift(0), source_shift(0)) ; 
	Tbl tdiff_shift_y = diffrel(lap_shift(1), source_shift(1)) ; 
	Tbl tdiff_shift_z = diffrel(lap_shift(2), source_shift(2)) ; 

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
	   

	    // Source
	    //--------
	    
	    Scalar sol_poisson(mp) ;

	    for(int i=0; i<=2; i++) 
		for(int j=i; j<=2; j++) {
		    
		    cout << "moyenne de la source 1 pour hij_auto" << endl ;
		    cout <<  norme(source1/(nr*nt*np)) << endl ;
		    cout << "moyenne de la source 2 pour hij_auto" << endl ;
		    cout <<  norme(source2/(nr*nt*np)) << endl ;
		    cout << "moyenne de la source 3 pour hij_auto" << endl ;
		    cout <<  norme(source3/(nr*nt*np)) << endl ;
		    cout << "moyenne de la source 4 pour hij_auto" << endl ;
		    cout <<  norme(source4/(nr*nt*np)) << endl ;
		    cout << "moyenne de la source5 pour hij_auto" << endl ;
		    cout <<  norme(source5/(nr*nt*np)) << endl ;
		    cout << "moyenne de la source6 pour hij_auto" << endl ;
		    cout <<  norme(source6/(nr*nt*np)) << endl ;
		    cout << "moyenne de la source7 pour hij_auto" << endl ;
		    cout <<  norme(source7/(nr*nt*np)) << endl ;
		    cout << "moyenne de la source8 pour hij_auto" << endl ;
		    cout <<  norme(source8/(nr*nt*np)) << endl ;
		    cout << "moyenne de la source9 pour hij_auto" << endl ;
		    cout <<  norme(source9/(nr*nt*np)) << endl ;
		    cout << "moyenne de la source pour hij_auto" << endl ;
		    cout <<  norme(source_tot/(nr*nt*np)) << endl << endl ;
  

		    // Resolution of the Poisson equations and
		    // Check: has the Poisson equation been correctly solved ?
		    // -----------------------------------------------------

     
		    if(i==0 && j==0) {
		 
			source_tot.std_spectral_base() ;
			sol_poisson = source_tot.poisson() ; 
			hij_auto.set(0,0) = sol_poisson ;

			Tbl tdiff_hij00 = diffrel(hij_auto(0,0).laplacian(), source_tot) ;  
			cout << "Relative error in the resolution of the equation for "
			     << "hij00_auto : " << endl ; 
			for (int l=0; l<nz; l++) {
			    cout << tdiff_hij00(l) << "  " ; 
			}
			cout << endl ;
			diff_hij00 = max(abs(tdiff_hij00)) ; 
		    }
	       	       
		    if(i==0 && j==1) {

			source_tot.std_spectral_base() ;
			sol_poisson = source_tot.poisson() ; 
			hij_auto.set(0,1) = sol_poisson ;
	    
			Tbl tdiff_hij10 = diffrel(hij_auto(0,1).laplacian(), source_tot) ;
			cout << 
			    "Relative error in the resolution of the equation for " 
			     << "hij10_auto : "  << endl ; 
			for (int l=0; l<nz; l++) {
			    cout << tdiff_hij10(l) << "  " ; 
			}
			cout << endl ;
			diff_hij10 = max(abs(tdiff_hij10)) ; 
		    }
	       
		    if(i==0 && j==2) {
		 
	     
			source_tot.std_spectral_base() ;
			sol_poisson = source_tot.poisson() ; 
			hij_auto.set(0,2) = sol_poisson ;

			Tbl tdiff_hij20 = diffrel(hij_auto(0,2).laplacian(), source_tot) ;
			cout << 
			    "Relative error in the resolution of the equation for "
			     << "hij20_auto : " << endl ; 
			for (int l=0; l<nz; l++) {
			    cout << tdiff_hij20(l) << "  " ; 
			}
			cout << endl ;
			diff_hij20 = max(abs(tdiff_hij20)) ; 
		    }
	     
		    if(i==1 && j==1) {
		 
			source_tot.std_spectral_base() ;
			sol_poisson = source_tot.poisson() ; 
			hij_auto.set(1,1) = sol_poisson ;

			Tbl tdiff_hij11 = diffrel(hij_auto(1,1).laplacian(), source_tot) ;
			cout << 
			    "Relative error in the resolution of the equation for "
			     << "hij11_auto : " << endl ; 
			for (int l=0; l<nz; l++) {
			    cout << tdiff_hij11(l) << "  " ; 
			}
			cout << endl ;
			diff_hij11 = max(abs(tdiff_hij11)) ; 
		    }
	       
		    if(i==1 && j==2) {
		 
			source_tot.std_spectral_base() ;
			sol_poisson = source_tot.poisson() ; 
			hij_auto.set(1,2) = sol_poisson ;

			Tbl tdiff_hij21 = diffrel(hij_auto(1,2).laplacian(), source_tot) ;
			cout << 
			    "Relative error in the resolution of the equation for "
			     << "hij21_auto : " << endl ; 
			for (int l=0; l<nz; l++) {
			    cout << tdiff_hij21(l) << "  " ; 
			}
			cout << endl ;
			diff_hij21 = max(abs(tdiff_hij21)) ; 
		    }
	     
		    if(i==2 && j==2) {
		 
			source_tot.std_spectral_base() ;
			sol_poisson = source_tot.poisson() ; 
			hij_auto.set(2,2) = sol_poisson ;

			Tbl tdiff_hij22 = diffrel(hij_auto(2,2).laplacian(), source_tot) ;
			cout << 
			    "Relative error in the resolution of the equation for "
			     << "hij22_auto : " << endl ;
			for (int l=0; l<nz; l++) {
			    cout << tdiff_hij22(l) << "  " ;
			}
			cout << endl ;
			diff_hij22 = max(abs(tdiff_hij22)) ;
		    }

		}
	
      
	cout << "hij00 auto"<<endl<< norme(hij_auto(0,0)/(nr*nt*np)) 
	     << endl << endl ;
	cout << "hij10 auto"<<endl<< norme(hij_auto(0,1)/(nr*nt*np)) 
	     << endl << endl ;
	cout << "hij20 auto"<<endl<< norme(hij_auto(0,2)/(nr*nt*np)) 
	     << endl << endl ;
	cout << "hij11 auto"<<endl<< norme(hij_auto(1,1)/(nr*nt*np)) 
	     << endl << endl ;
	cout << "hij21 auto"<<endl<< norme(hij_auto(1,2)/(nr*nt*np)) 
	     << endl << endl ;
	cout << "hij22 auto"<<endl<< norme(hij_auto(2,2)/(nr*nt*np)) 
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





    
    
