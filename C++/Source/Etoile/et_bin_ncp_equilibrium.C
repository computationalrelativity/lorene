/*
 *   Method of class Et_bin_ncp to compute an equilibrium configuration
 *
 *  (see file et_bin_ncp.h for documentation).
 */

/*
 *   Copyright (c) 2002  Francois Limousin
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

char et_bin_ncp_equilibrium_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/02/06 17:26:17  f_limousin
 * *** empty log message ***
 *
 * Revision 1.2  2001/12/11 06:44:41  e_gourgoulhon
 * template files
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "et_bin_ncp.h"
#include "param.h"
#include "graphique.h"
#include "utilitaires.h"


void Et_bin_ncp::equilibrium(double ent_c, int mermax, int mermax_poisson, 
			 double relax_poisson, int mermax_potvit, 
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
    double& diff_loggamma = diff.set(3) ; 
    double& diff_shift_x = diff.set(4) ; 
    double& diff_shift_y = diff.set(5) ; 
    double& diff_shift_z = diff.set(6) ; 
    double& diff_gtilde00 = diff.set(7) ; 
    double& diff_gtilde10 = diff.set(8) ; 
    double& diff_gtilde20 = diff.set(9) ; 
    double& diff_gtilde11 = diff.set(10) ; 
    double& diff_gtilde21 = diff.set(11) ; 
    double& diff_gtilde22 = diff.set(12) ; 



    // Parameters for the function Map_et::adapt
    // -----------------------------------------
    
    Param par_adapt ; 
    int nitermax = 100 ;  
    int niter ; 
    int adapt_flag = 1 ;    //  1 = performs the full computation, 
			    //  0 = performs only the rescaling by 
			    //      the factor alpha_r
    //##    int nz_search = nzet + 1 ;  // Number of domains for searching the enthalpy
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
    par_adapt.add_int(nz_search, 2) ;	// number of domains to search 					// the enthalpy isosurface
    par_adapt.add_int(adapt_flag, 3) ; //  1 = performs the full computation, 
				       //  0 = performs only the rescaling by 
				       //      the factor alpha_r
    par_adapt.add_int(j_b, 4) ; //  theta index of the collocation point 
			        //  (theta_*, phi_*)
    par_adapt.add_int(k_b, 5) ; //  theta index of the collocation point 
			        //  (theta_*, phi_*)

    par_adapt.add_int_mod(niter, 0) ;  //  number of iterations actually used in 
				       //  the secant method
    
    par_adapt.add_double(precis_secant, 0) ; // required absolute precision in 
					     // the determination of zeros by 
					     // the secant method
    par_adapt.add_double(reg_map, 1)	;  // 1. = regular mapping, 0 = contracting mapping
    
    par_adapt.add_double(alpha_r, 2) ;	    // factor by which all the radial 
					    // distances will be multiplied 
    	   
    par_adapt.add_tbl(ent_limit, 0) ;	// array of values of the field ent 
				        // to define the isosurfaces. 			   

    // Parameters for the function Map_et::poisson for logn_auto
    // ---------------------------------------------------------

    double precis_poisson = 1.e-16 ;     

    Param par_poisson1 ; 

    par_poisson1.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson1.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson1.add_double(precis_poisson, 1) ; // required precision
    par_poisson1.add_int_mod(niter, 0) ;  //  number of iterations actually used 
    par_poisson1.add_cmp_mod( ssjm1_logn ) ; 
					   
    // Parameters for the function Map_et::poisson for log(gamma)
    // ----------------------------------------------------------

    Param par_poisson2 ; 

    par_poisson2.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson2.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson2.add_double(precis_poisson, 1) ; // required precision
    par_poisson2.add_int_mod(niter, 0) ;  //  number of iterations actually used 
    //    par_poisson1.add_cmp_mod( ssjm1_loggamma ) ; 






    // Parameters for the function Map_et::poisson for gtildeij
    // ---------------------------------------------------------

    Param par_poisson3 ; 

    par_poisson3.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson3.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson3.add_double(precis_poisson, 1) ; // required precision
    par_poisson3.add_int_mod(niter, 0) ;  //  number of iterations actually used 
 					   
    // Parameters for the function Tenseur::poisson_vect
    // -------------------------------------------------

    Param par_poisson_vect ; 

    par_poisson_vect.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson_vect.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson_vect.add_double(precis_poisson, 1) ; // required precision
    par_poisson_vect.add_cmp_mod( ssjm1_khi ) ; 
    par_poisson_vect.add_tenseur_mod( ssjm1_wshift ) ; 
    par_poisson_vect.add_int_mod(niter, 0) ;  

 					   
    // External potential
    // See Eq (99) from Gourgoulhon et al. (2001)
    // ------------------
    
    Tenseur pot_ext = logn_comp + pot_centri + loggam ;
//##
//	des_coupe_z(pot_ext(), 0., 1, "pot_ext", &(ent()) ) ; 
//##
    
    Tenseur ent_jm1 = ent ;	// Enthalpy at previous step
    
    Tenseur source(mp) ;    // source term in the equation for logn_auto
			    // and beta_auto
			    
    Tenseur source_shift(mp, 1, CON, ref_triad) ;  // source term in the equation 
						   //  for shift_auto



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

	//-----------------------------------------------------
	// Computation of the new radial scale
	//--------------------------------------------------

	// alpha_r (r = alpha_r r') is determined so that the enthalpy
	// takes the requested value ent_b at the stellar surface
	
	// Values at the center of the star:
	double logn_auto_c  = logn_auto()(0, 0, 0, 0) ; 
	double pot_ext_c  = pot_ext()(0, 0, 0, 0) ; 

	// Search for the reference point (theta_*, phi_*) [notation of
	//  Bonazzola, Gourgoulhon & Marck PRD 58, 104020 (1998)]
	//  at the surface of the star
	// ------------------------------------------------------------
	double alpha_r2 = 0 ; 
	for (int k=0; k<mg->get_np(l_b); k++) {
	    for (int j=0; j<mg->get_nt(l_b); j++) {
		
		double pot_ext_b  = pot_ext()(l_b, k, j, i_b) ; 
		double logn_auto_b  = logn_auto()(l_b, k, j, i_b) ; 


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
	logn_auto_regu = alpha_r2 * logn_auto_regu ;
	logn_auto_c  = logn_auto()(0, 0, 0, 0) ;


	//------------------------------------------------------------
	// Change the values of the inner points of the second domain
	// by those of the outer points of the first domain
	//------------------------------------------------------------

	(logn_auto().va).smooth(nzet, (logn_auto.set()).va) ;


	//------------------------------------------
	// First integral	--> enthalpy in all space
	// See Eq (98) from Gourgoulhon et al. (2001)
	//-------------------------------------------

	ent = (ent_c + logn_auto_c + pot_ext_c) - logn_auto - pot_ext ;

	(ent().va).smooth(nzet, (ent.set()).va) ;

	//----------------------------------------------------
	// Adaptation of the mapping to the new enthalpy field
	//----------------------------------------------------
    
	// Shall the adaptation be performed (cusp) ?
	// ------------------------------------------
	
	double dent_eq = ent().dsdr().val_point(ray_eq(),M_PI/2.,0.) ;
	double dent_pole = ent().dsdr().val_point(ray_pole(),0.,0.) ;
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
	    ent_limit.set(l) = ent()(l, k_b, j_b, i_b) ; 
	}
	ent_limit.set(nzet-1) = ent_b  ; 

	Map_et mp_prev = mp_et ; 


	mp.adapt(ent(), par_adapt) ; 

	// Readjustment of the external boundary of domain l=nzet
	// to keep a fixed ratio with respect to star's surface
	
	int n_resize ;
	//      	if (nz > 4) {
	//       	  n_resize = nz - 4 ;
    	if (nz > 3) {
   	  n_resize = nz - 3 ;
	}
	else {
	  n_resize = nzet ;
	}

	double rr_in = mp.val_r(nzet,-1., M_PI/2, 0.) ; 
	double rr_out = mp.val_r(n_resize,1., M_PI/2, 0.) ; 

	mp.resize(n_resize, rr_in/rr_out * fact_resize(0)) ; 

//##
//	des_coupe_z(ent(), 0., 1, "ent after adapt", &(ent()) ) ; 
//##
	//----------------------------------------------------
	// Computation of the enthalpy at the new grid points
	//----------------------------------------------------
	
	mp_prev.homothetie(alpha_r) ; 
	
	mp.reevaluate_symy(&mp_prev, nzet+1, ent.set()) ; 

//	des_coupe_z(ent(), 0., 1, "ent after reevaluate", &(ent()) ) ; 

	double ent_s_max = -1 ; 
	int k_s_max = -1 ; 
	int j_s_max = -1 ; 
	for (int k=0; k<mg->get_np(l_b); k++) {
	    for (int j=0; j<mg->get_nt(l_b); j++) {
		double xx = fabs( ent()(l_b, k, j, i_b) ) ;
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

	
	Tenseur dcon_loggamma_auto(loggamma_auto.derive_con(gtilde)) ;
	Tenseur dcov_loggamma_auto(loggamma_auto.derive_cov(gtilde)) ;
	Tenseur dcon_loggamma(loggamma.derive_con(gtilde)) ;
	Tenseur dcov_loggamma(loggamma.derive_cov(gtilde)) ;

	Tenseur dcov_lognauto(logn_auto.derive_cov(gtilde)) ;
	Tenseur dcon_lognauto(logn_auto.derive_con(gtilde)) ;
	Tenseur dcov_logn((logn_auto + logn_comp).derive_cov(gtilde)) ;
	Tenseur dcon_logn((logn_auto + logn_comp).derive_con(gtilde)) ;	


	source = - contract(dcon_loggamma, 0, dcov_lognauto, 0)/6. + pow(gamma,(1./3.)) * (qpig * (ener_euler + s_euler) + kcar_auto + kcar_comp) - contract(dcov_lognauto, 0, dcon_logn, 0) ;
	

	source.set_std_base() ;


	// Resolution of the Poisson equation 
	// ----------------------------------

	source().poisson(par_poisson1, logn_auto.set()) ; 

	// Construct logn_auto_regu for et_bin_ncp_upmetr.C
	// ------------------------------------------------

	logn_auto_regu = logn_auto ;

	// Check: has the Poisson equation been correctly solved ?
	// -----------------------------------------------------

	Tbl tdiff_logn = diffrel(logn_auto().laplacien(), source()) ;
	cout << 
	"Relative error in the resolution of the equation for logn_auto : "
	<< endl ; 
	for (int l=0; l<nz; l++) {
	    cout << tdiff_logn(l) << "  " ; 
	}
	cout << endl ;
	diff_logn = max(abs(tdiff_logn)) ; 


        //--------------------------------------------------------
        // Poisson equation for log(gamma)
        //--------------------------------------------------------

	// Source
	//--------


	Tenseur dcovdcon_loggamma_auto(loggamma_auto.derive_con(gtilde).derive_cov(gtilde)) ;


	//	source = 3./2. * gtilde.ricci_scal() - contract(dcon_loggamma_auto, 0, dcov_loggamma, 0)/12. - pow(gamma,(1./3.)) * (6.*qpig*ener_euler + 3./2. * (kcar_auto + kcar_comp)) + (loggamma_auto().laplacien() - dcovdcon_loggamma_auto) ;

	 // Resolution of the Poisson equation 
	 // ----------------------------------

         source().poisson(par_poisson2, loggamma.set()) ; 
	 gamma = exp(loggamma) ;
	 

         // Check: has the Poisson equation been correctly solved 
         // -----------------------------------------------------
    
         Tbl tdiff_loggamma = diffrel(loggamma().laplacien(), source()) ;
	 cout << 
         "Relative error in the resolution of the equation for loggamma : "
	      << endl ; 
         for (int l=0; l<nz; l++) {
	   cout << tdiff_loggamma(l) << "  " ; 
	    }
         cout << endl ;
         diff_loggamma = max(abs(tdiff_loggamma)) ; 



         //--------------------------------------------------------
         // Vector Poisson equation for shift_auto 
         //--------------------------------------------------------

         // Source
         //--------

	 Tenseur source1(mp) ;
	 Tenseur source2(mp) ;
	 Tenseur source3(mp) ;

	 Tenseur dcov_nnn(nnn.derive_cov(gtilde)) ;
	 Tenseur dconidcovj_shiftj_auto(contract(shift_auto.derive_cov(gtilde), 0, 1).derive_con(gtilde)) ;

	 source1 = (-4.*qpig) * nnn * pow(gamma,(1./3.))*(ener_euler + press) * u_euler ;

	 source2 = - dconidcovj_shiftj_auto/3. - contract(operator*(gtilde.con(), contract(operator*(gtilde.ricci(), shift_auto),1 ,2)), 1, 2) ;

	 source3 = - pow(gamma,(-1./3.))*contract(tkij_auto, 1, (nnn*dcov_loggamma - 2.* dcov_nnn), 0) ;

	 
	 for (int i=0; i<3; i++) {
	   source_shift.set(i) = source1(i) + source2(i) + source3(i)  
	   + shift_auto(i).laplacien() 
	      - (contract((shift_auto.derive_cov(gtilde)).derive_con(gtilde), 1, 2))(i) ; 
	 }

	 source_shift.set_std_base() ; 	

	 // Resolution of the Poisson equation 
	 // ----------------------------------

	 // Filter for the source of shift vector

	 for (int i=0; i<3; i++) {

           if (source_shift(i).get_etat() != ETATZERO)
	     source_shift.set(i).filtre(4) ;

	 }


	 // For Tenseur::poisson_vect, the triad must be the mapping triad,
	 // not the reference one:
	    
	 source_shift.change_triad( mp.get_bvect_cart() ) ; 

         for (int i=0; i<3; i++) {
	   if(source_shift(i).dz_nonzero()) {
	     assert( source_shift(i).get_dzpuis() == 4 ) ; 
	   }
	   else{
	     (source_shift.set(i)).set_dzpuis(4) ; 
	   }
	 }

         //##
         // source_shift.dec2_dzpuis() ;    // dzpuis 4 -> 
         double lambda_shift = double(0) ; 
	
         source_shift.poisson_vect(lambda_shift, par_poisson_vect, 
      			      shift_auto, w_shift, khi_shift) ;      


         // Check: has the equation for shift_auto been correctly solved ?
         // --------------------------------------------------------------
	             // Divergence of shift_auto : 
         Tenseur divn = contract(shift_auto.gradient(), 0, 1) ; 
         divn.dec2_dzpuis() ;    // dzpuis 2 -> 0
	    
         // Grad(div) : 
         Tenseur graddivn = divn.gradient() ; 
         graddivn.inc2_dzpuis() ;    // dzpuis 2 -> 4
	    
         // Full operator : 
         Tenseur lap_shift(mp, 1, CON, mp.get_bvect_cart() ) ;  
         lap_shift.set_etat_qcq() ; 
	 for (int i=0; i<3; i++) {
	   lap_shift.set(i) = shift_auto(i).laplacien() 
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

	 // Final result
         // ------------
         // The output of Tenseur::poisson_vect is on the mapping triad,
         // it should therefore be transformed to components on the
         // reference triad :
	    
         shift_auto.change_triad( ref_triad ) ; 


         //--------------------------------------------------------
	 // Poisson equation for gij
         //--------------------------------------------------------

         // Source
	 //--------
	 
	 
	 for(int i=0; i<=2; i++) {
	   for(int j=i; j<=2; j++) {

	     
	     Tenseur dcovidcovj_nnn(((nnn.derive_cov(gtilde)).derive_cov(gtilde))(i,j)) ;
	     Tenseur dcondcov_nnn(contract(nnn.derive_cov(gtilde).derive_con(gtilde),0,1)) ;
	     Tenseur dcovidcovj_loggamma_auto(((loggamma_auto.derive_cov(gtilde)).derive_cov(gtilde))(i,j)) ;
	     Tenseur dcondcov_loggamma_auto(contract((loggamma_auto.derive_cov(gtilde)).derive_con(gtilde), 0, 1)) ;

	     Cmp source4 (mp) ;
	     Cmp source5 (mp) ;
	     Cmp source6 (mp) ;
	     Cmp source_tot (mp) ;

	     source4 = - 1/nnn()*dcovidcovj_nnn() + (dcov_loggamma(i)*dcov_lognauto(j) + dcov_loggamma(j)*dcov_lognauto(i))/6. + (dcondcov_nnn()/(3*nnn()) - contract(dcon_loggamma, 0, dcov_lognauto, 0)()/9.)*(gtilde_auto.cov())(i,j) ;

	     source5 = -(dcovidcovj_loggamma_auto() - dcov_loggamma_auto(i)*dcov_loggamma(j)/6.)/6. - ((gtilde.ricci_scal())() - (dcondcov_loggamma_auto() - contract(dcon_loggamma_auto, 0, dcov_loggamma, 0)()/6.)/6.)/3.*(gtilde.cov())(i,j) ;

	     Tenseur atilde2( contract( contract( contract( contract( contract(  contract( operator*( operator*( operator*( operator*( operator*(tkij_auto, tkij_auto + tkij_comp), gtilde.cov() ), gtilde.cov() ), gtilde.cov() ), gtilde.cov() ), 9, 13), 7, 11), 5, 9), 3, 7), 0, 3), 0, 3) ) ;

	     source6 = -2.*pow(gamma(),(-2./3.))*atilde2() - 2.*qpig*(    - s_euler()*(gtilde.cov())(i,j)/3.) ;

	 
	     source_tot = 2.*source4 + 2.*source5 + 2.*pow(gamma(),(1./3.))*source6 + ((gtilde.ricci())() + ((gtilde_auto.cov())(i,j)).laplacien()/2.) ;

       	
	source.set_std_base() ; 

	// Resolution of the Poisson equations 
	// -----------------------------------
	
	Cmp gtildeij (mp) ;
	source().poisson(par_poisson3, gtildeij) ; 

	Cmp gtilde00 = gtilde.set_cov(0,0) ;
	Cmp gtilde10 = gtilde.set_cov(1,0) ;
	Cmp gtilde20 = gtilde.set_cov(2,0) ;
	Cmp gtilde11 = gtilde.set_cov(1,1) ;
	Cmp gtilde21 = gtilde.set_cov(2,1) ;
        Cmp gtilde22 = gtilde.set_cov(2,2) ;

	// Check: has the Poisson equation been correctly solved ?
	// -----------------------------------------------------

	if(i==0 && j==0) {
	  gtilde00 = gtildeij ; 
	  Tbl tdiff_gtilde00 = diffrel(gtilde00.laplacien(), source()) ;
	  cout << 
	    "Relative error in the resolution of the equation for gtilde00 : "
	       << endl ; 
	  for (int l=0; l<nz; l++) {
	    cout << tdiff_gtilde00(l) << "  " ; 
	  }
	  cout << endl ;
	  diff_gtilde00 = max(abs(tdiff_gtilde00)) ; 
	}


	if(i==1 && j==0) {
	  gtilde10 = gtildeij ; 
	  Tbl tdiff_gtilde10 = diffrel(gtilde10.laplacien(), source()) ;
	  cout << 
	    "Relative error in the resolution of the equation for gtilde10 : "
	       << endl ; 
	  for (int l=0; l<nz; l++) {
	    cout << tdiff_gtilde10(l) << "  " ; 
	  }
	  cout << endl ;
	  diff_gtilde10 = max(abs(tdiff_gtilde10)) ; 
	}

	if(i==2 && j==0) {
	  gtilde20 = gtildeij ; 
	  Tbl tdiff_gtilde20 = diffrel(gtilde20.laplacien(), source()) ;
	  cout << 
	    "Relative error in the resolution of the equation for gtilde20 : "
	       << endl ; 
	  for (int l=0; l<nz; l++) {
	    cout << tdiff_gtilde20(l) << "  " ; 
	  }
	  cout << endl ;
	  diff_gtilde20 = max(abs(tdiff_gtilde20)) ; 
	}

	if(i==1 && j==1) {
	  gtilde11 = gtildeij ; 
	  Tbl tdiff_gtilde11 = diffrel(gtilde11.laplacien(), source()) ;
	  cout << 
	    "Relative error in the resolution of the equation for gtilde11 : "
	       << endl ; 
	  for (int l=0; l<nz; l++) {
	    cout << tdiff_gtilde11(l) << "  " ; 
	  }
	  cout << endl ;
	  diff_gtilde11 = max(abs(tdiff_gtilde11)) ; 
	}

	if(i==2 && j==1) {
	  gtilde21 = gtildeij ; 
	  Tbl tdiff_gtilde21 = diffrel(gtilde21.laplacien(), source()) ;
	  cout << 
	    "Relative error in the resolution of the equation for gtilde21 : "
	       << endl ; 
	  for (int l=0; l<nz; l++) {
	    cout << tdiff_gtilde21(l) << "  " ; 
	  }
	  cout << endl ;
	  diff_gtilde21 = max(abs(tdiff_gtilde21)) ; 
	}

	if(i==2 && j==2) {
	  gtilde22 = gtildeij ; 
	  Tbl tdiff_gtilde22 = diffrel(gtilde22.laplacien(), source()) ;
	  cout << 
	    "Relative error in the resolution of the equation for gtilde22 : "
	       << endl ; 
	  for (int l=0; l<nz; l++) {
	    cout << tdiff_gtilde22(l) << "  " ; 
	  }
	  cout << endl ;
	  diff_gtilde22 = max(abs(tdiff_gtilde22)) ; 
	}
	
	   }
	 }            // End of relativistic equations	


	//-------------------------------------------------
	//  Relative change in enthalpy
	//-------------------------------------------------

	Tbl diff_ent_tbl = diffrel( ent(), ent_jm1() ) ; 
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



