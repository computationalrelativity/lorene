
/*
 *   Method of class Et_bin_ncp to compute an equilibrium configuration
 *
 *  (see file et_bin_ncp.h for documentation).
 */

/*
 *   Copyright (c) 2002-2003  Francois Limousin
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
 * Revision 1.8  2004/03/25 10:29:04  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.7  2003/12/05 14:50:26  j_novak
 * To suppress some warnings...
 *
 * Revision 1.6  2003/10/17 13:00:02  f_limousin
 * Changes from set_cov() to get_cov() for metrics.
 *
 * Revision 1.5  2003/10/13 10:29:56  f_limousin
 * *** empty log message ***
 *
 * Revision 1.4  2003/06/20 14:28:13  f_limousin
 * Many modif.
 *
 * Revision 1.3  2003/03/04 08:47:47  j_novak
 * Fixed a problem with display.
 *
 * Revision 1.2  2003/03/03 19:23:27  f_limousin
 * Many modifications. Add Cmp ssjm1_gtildeij, change of dzpuis...
 *
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
#include "unites.h"	    


void Et_bin_ncp::equilibrium(double ent_c, int mermax, int mermax_poisson, 
			     double relax_poisson, int mermax_potvit, 
			     double relax_potvit, double thres_adapt,
			     const Tbl& fact_resize, Tbl& diff) {


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
  double& diff_a_car = diff.set(3) ; 
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
			   

  // Parameters for the function Map_et::poisson for logn_auto
  // ---------------------------------------------------------

  double precis_poisson = 1.e-16 ;     

  Param par_poisson1 ; 

  par_poisson1.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson1.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson1.add_double(precis_poisson, 1) ; // required precision
  par_poisson1.add_int_mod(niter, 0) ; // number of iterations actually used 
  par_poisson1.add_cmp_mod( ssjm1_logn ) ; 
					   
  // Parameters for the function Map_et::poisson for a_car_auto
  // ---------------------------------------------------------------

  Param par_poisson2 ; 

  par_poisson2.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson2.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson2.add_double(precis_poisson, 1) ; // required precision
  par_poisson2.add_int_mod(niter, 0) ; // number of iterations actually used 
  par_poisson2.add_cmp_mod( ssjm1_a_car ) ; 

  // Parameters for the function Map_et::poisson for hij00_auto
  // -------------------------------------------------------------

  Param par_poisson3 ; 

  par_poisson3.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson3.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson3.add_double(precis_poisson, 1) ; // required precision
  par_poisson3.add_int_mod(niter, 0) ; // number of iterations actually used 
  par_poisson3.add_cmp_mod( ssjm1_hij00 ) ; 
 					   
  // Parameters for the function Map_et::poisson for hij10_auto
  // -------------------------------------------------------------

  Param par_poisson4 ; 

  par_poisson4.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson4.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson4.add_double(precis_poisson, 1) ; // required precision
  par_poisson4.add_int_mod(niter, 0) ; // number of iterations actually used 
  par_poisson4.add_cmp_mod( ssjm1_hij10 ) ; 

  // Parameters for the function Map_et::poisson for hij20_auto
  // -------------------------------------------------------------

  Param par_poisson5 ; 

  par_poisson5.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson5.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson5.add_double(precis_poisson, 1) ; // required precision
  par_poisson5.add_int_mod(niter, 0) ; // number of iterations actually used 
  par_poisson5.add_cmp_mod( ssjm1_hij20 ) ; 

  // Parameters for the function Map_et::poisson for hij11_auto
  // -------------------------------------------------------------

  Param par_poisson6 ; 

  par_poisson6.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson6.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson6.add_double(precis_poisson, 1) ; // required precision
  par_poisson6.add_int_mod(niter, 0) ; // number of iterations actually used 
  par_poisson6.add_cmp_mod( ssjm1_hij11 ) ; 

  // Parameters for the function Map_et::poisson for hij21_auto
  // -------------------------------------------------------------

  Param par_poisson7 ; 

  par_poisson7.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson7.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson7.add_double(precis_poisson, 1) ; // required precision
  par_poisson7.add_int_mod(niter, 0) ; // number of iterations actually used 
  par_poisson7.add_cmp_mod( ssjm1_hij21 ) ; 

  // Parameters for the function Map_et::poisson for hij22_auto
  // -------------------------------------------------------------

  Param par_poisson8 ; 

  par_poisson8.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson8.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson8.add_double(precis_poisson, 1) ; // required precision
  par_poisson8.add_int_mod(niter, 0) ; // number of iterations actually used 
  par_poisson8.add_cmp_mod( ssjm1_hij22 ) ; 

  // Parameters for the function Tenseur::poisson_vect
  // -------------------------------------------------

  Param par_poisson_vect ; 

  par_poisson_vect.add_int(mermax_poisson,  0) ;  // maximum number of 
  //  iterations
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
    
  Tenseur source(mp) ;    // source term in the equation for logn_auto,
  // loggamma_auto

  Tenseur source_tot(mp) ; // source term in the equation for gtildeij
			    
  Tenseur source_shift(mp, 1, CON, mp.get_bvect_cart()) ;  // source term 
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
  
  // Derivatives of shift
  // --------------------

  Tenseur dcov_shift_auto = shift_auto.derive_cov(flat) ;
  dcov_shift_auto.dec2_dzpuis() ;
  Tenseur dcovdcov_shift_auto = dcov_shift_auto.derive_cov(flat) ;
  dcovdcov_shift_auto.inc2_dzpuis() ;

  // Derivatives of gtilde
  // ---------------------

     const Tenseur& dcov_gtilde = (gtilde.get_cov()).derive_cov(flat) ;
  

  // Derivatives of hij
  // ------------------

  Tenseur dcov_hij = hij.derive_cov(flat) ;
  Tenseur dcov_hij_auto = hij_auto.derive_cov(flat) ;

  Tenseur dcov_hij_auto2 = hij.derive_cov(flat) ;
  dcov_hij_auto2.dec2_dzpuis() ; 
  Tenseur dcovdcov_hij_auto = dcov_hij_auto2.derive_cov(flat) ;
  dcovdcov_hij_auto.inc2_dzpuis() ;


    // Source 
    //--------

    int nr = mp.get_mg()->get_nr(0) ;
    int nt = mp.get_mg()->get_nt(0) ;
    int np = mp.get_mg()->get_np(0) ;
    
    Tenseur source1(mp) ;
    Tenseur source2(mp) ;
    Tenseur source3(mp) ;
    Tenseur source4(mp) ;
    Tenseur source5(mp) ;
    Tenseur source6(mp) ;
    Tenseur source7(mp) ;
    Tenseur source8(mp) ;
    Tenseur source9(mp) ;
    Tenseur source10(mp) ;

    source1 = qpig * pow(gamma,1./3.) % (ener_euler + s_euler) ;


    source2 = kcar_auto + kcar_comp ;
 
    source3 = - contract_desal( dcov_logn_auto, 0, dcon_logn + pow(gamma, -1./3.) % 
		   contract_desal(gtilde.con(), 0, dcov_acar, 0)*0.5, 0) ;
		
    source4 = - contract(contract(hij, 0, dcovdcov_logn_auto, 0), 0, 1) ;

    source5 = - contract(contract(hij, 0, dcov_logn, 0), 0,dcov_logn_auto, 0) ;

    source.set_etat_qcq() ;
    source.set() = source1() + source2() + source3() + source4() + source5();
 
    cout << "moyenne de la source 1 pour logn_auto" << endl ;
    cout <<  norme(source1()/(nr*nt*np)) << endl ;
    cout << "moyenne de la source 2 pour logn_auto" << endl ;
    cout <<  norme(source2()/(nr*nt*np)) << endl ;
    cout << "moyenne de la source 3 pour logn_auto" << endl ;
    cout <<  norme(source3()/(nr*nt*np)) << endl ;
    cout << "moyenne de la source 4 pour logn_auto" << endl ;
    cout <<  norme(source4()/(nr*nt*np)) << endl ;
    cout << "moyenne de la source 5 pour logn_auto" << endl ;
    cout <<  norme(source5()/(nr*nt*np)) << endl ;
    cout << "moyenne de la source pour logn_auto" << endl ;
    cout <<  norme(source()/(nr*nt*np)) << endl ;
	
	
    source.set_std_base() ;
    /*
    des_profile(source1(), 0, 20, 0, 0) ;
    des_profile(source2(), 0, 20, 0, 0) ;
    des_profile(source3(), 0, 20, 0, 0) ;
    des_profile(source4(), 0, 20, 0, 0) ;
    des_profile(source(), 0, 20, 0, 0) ;
    */

    // Resolution of the Poisson equation 
    // ----------------------------------
    

    //source.set().filtre(4) ;
    //   source.annule(nz-1) ;
 

    source().poisson(par_poisson1, logn_auto.set()) ; 

    cout << "logn_auto" << endl << norme(logn_auto()/(nr*nt*np)) << endl ;
    /*
    des_profile(logn_auto(), 0, 20, 0, 0) ;
    des_coef_xi(logn_auto().va, 0, 0, 0) ;
    des_coef_xi(logn_auto().va, 1, 0, 0) ;
    des_coef_xi(logn_auto().va, 2, 0, 0) ;
    */
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

    cout << "Rel error: " << norme(logn_auto().laplacien()- source())
      / max( max( source() ) ) << endl ; 

    
    //--------------------------------------------------------
    // Poisson equation for a_car
    //--------------------------------------------------------

    // Source
    //--------

	 
    source1 = a_car_auto * ricci_scal * 0.5 ;

    source2 = 3./4.*pow(gamma, -1./3.)%
	contract_desal(contract_desal(gtilde.con(), 0, dcov_acar, 0), 
		       0, dcov_acar_auto, 0) ;
		
    source3 =  -2*qpig*pow(gamma, 2./3.) % ener_euler 
	-2. * pow(gamma,1./3.) % (kcar_auto + kcar_comp) ;  //a voir 

    source4 = - contract(contract(hij, 0, dcovdcov_acar_auto, 0), 0, 1) ;

    source.set_etat_qcq() ;
    source.set() = source1() + source2() + source3() + source4() ; 
  	
    cout << "moyenne de la source 1 pour acar_auto" << endl ;
    cout <<  norme(source1()/(nr*nt*np)) << endl ;
    cout << "moyenne de la source 2 pour acar_auto" << endl ;
    cout <<  norme(source2()/(nr*nt*np)) << endl ;
    cout << "moyenne de la source 3 pour acar_auto" << endl ;
    cout <<  norme(source3()/(nr*nt*np)) << endl ;
    cout << "moyenne de la source 4 pour acar_auto" << endl ;
    cout <<  norme(source4()/(nr*nt*np)) << endl ;
    cout << "moyenne de la source pour acar_auto" << endl ;
    cout <<  norme(source()/(nr*nt*np)) << endl ;
	

    source.set_std_base() ;
    /*
    des_profile(source(), 0, 20, 0, 0) ;
    des_coef_xi(source().va, 0, 0, 0) ;
    des_coef_xi(source().va, 1, 0, 0) ;
    des_coef_xi(source().va, 2, 0, 0) ;
    */
    // des_profile(source(), 0, 20, 0, 0) ;

    // Resolution of the Poisson equation 
    // ----------------------------------

    //source.set().filtre(4) ;
    //   source.annule(nz-1) ;
 
    source().poisson(par_poisson2, a_car_auto.set()) ; 
    /*
    des_profile(a_car_auto(), 0, 20, 0, 0) ;
    des_coef_xi(a_car_auto().va, 0, 0, 0) ;
    des_coef_xi(a_car_auto().va, 1, 0, 0) ;
    des_coef_xi(a_car_auto().va, 2, 0, 0) ;
    */
    // Check: has the Poisson equation been correctly solved 
    // -----------------------------------------------------
    
    Tbl tdiff_a_car = diffrel(a_car_auto().laplacien(), source()) ;
    cout << 
      "Relative error in the resolution of the equation for a_car : "
	 << endl ; 
    for (int l=0; l<nz; l++) {
      cout << tdiff_a_car(l) << "  " ; 
    }
    cout << endl ;
    diff_a_car = max(abs(tdiff_a_car)) ; 

    cout << "Rel error: " << norme(a_car_auto().laplacien()- source())
      / max( max( source() ) ) << endl ; 

    a_car_auto.set() = a_car_auto() + decouple  ;

    cout << "acar_auto" << endl << norme(a_car_auto()/(nr*nt*np)) << endl ;

    //--------------------------------------------------------
    // Vector Poisson equation for shift_auto 
    //--------------------------------------------------------

    // Source
    //--------


    Tenseur source1_shift(mp, 1, CON, mp.get_bvect_cart()) ;
    Tenseur source2_shift(mp, 1, CON, mp.get_bvect_cart()) ;
    Tenseur source3_shift(mp, 1, CON, mp.get_bvect_cart()) ;
    Tenseur source4_shift(mp, 1, CON, mp.get_bvect_cart()) ;
    Tenseur source5_shift(mp, 1, CON, mp.get_bvect_cart()) ;

    Tenseur temp = contract(deltakij*shift, 2, 3) ;
    temp.dec2_dzpuis() ;
    Tenseur temp2 = temp.derive_cov(flat) ;
    temp2.inc2_dzpuis() ;

    source1_shift = (-4.*qpig) * nnn % pow(gamma,(1./3.))
      %(ener_euler + press) % u_euler ;

    source2_shift =  - contract(gtilde.con(), 1, contract(ricci_auto, 1,
					   	shift, 0), 0) ;
								
    source3_shift = + nnn % contract_desal(tkij_auto, 1
			        , (3.*dcov_acar%pow(a_car, -1.) - 2.* dcov_logn), 0) ;
		      
 	 
    source4_shift = 0;/*- contract(contract(gtilde.con()%temp2, 0, 4), 0, 1)
     - contract(contract_desal(  gtilde.con() % 
         (contract(deltakij%shift.derive_cov(flat), 2, 4) 
      - contract(deltakij%contract(deltakij%shift, 2, 3), 2, 3)), 0, 4), 0, 2)
     + contract(contract(  gtilde.con() %
         (contract(deltakij%shift.derive_cov(flat), 0, 3) 
     + contract(deltakij%contract(deltakij%shift, 2, 3), 0, 4)),0, 3), 0, 1) ;
		      */			          

    source5_shift = - contract(hij, 1, contract(dcovdcov_shift_auto,1,2),0)/3. 
	- contract( contract(hij, 0, dcovdcov_shift_auto, 0), 0, 1) ;

    

    source_shift.set_etat_qcq() ;
    for (int i=0; i<3; i++) {
      source_shift.set(i) = source1_shift(i) + source2_shift(i) 
	+ source3_shift(i) + source4_shift(i) + source5_shift(i) ; 

    }
      
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
    cout << "moyenne de la source 5 pour shift_auto" << endl ;
    cout <<  norme(source5_shift(0)/(nr*nt*np)) << endl ;
    cout <<  norme(source5_shift(1)/(nr*nt*np)) << endl ;
    cout <<  norme(source5_shift(2)/(nr*nt*np)) << endl ;
    cout << "moyenne de la source pour shift_auto" << endl ;
    cout <<  norme(source_shift(0)/(nr*nt*np)) << endl ;
    cout <<  norme(source_shift(1)/(nr*nt*np)) << endl ;
    cout <<  norme(source_shift(2)/(nr*nt*np)) << endl ;
	
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
    double lambda_shift = double(1)/double(3) ; 

    //    source_shift.annule(nz-1) ;
	
    source_shift.poisson_vect(lambda_shift, par_poisson_vect, 
      			      shift_auto, w_shift, khi_shift) ;      

    cout << "shift_auto" << endl << norme(shift_auto(0)/(nr*nt*np)) << endl << norme(shift_auto(1)/(nr*nt*np)) << endl << norme(shift_auto(2)/(nr*nt*np)) << endl ; 

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
	 
    /*
    des_profile(shift_auto(1), 0, 20, 0, 0) ;
    des_coef_xi(shift_auto(1).va, 0, 0, 0) ;
    des_coef_xi(shift_auto(1).va, 1, 0, 0) ;
    des_coef_xi(shift_auto(1).va, 2, 0, 0) ;
    */ 

    if (!conf_flat){
	   
      //--------------------------------------------------------
      // Poisson equation for hij
      //--------------------------------------------------------
	   
      // Source
      //--------
   
 
      const Tenseur& dcov_qq = qq.derive_cov(flat) ;
      
      Tenseur dcov_qq2 = qq.derive_cov(flat) ;
      dcov_qq2.dec2_dzpuis() ;
      Tenseur dcovdcov_qq = dcov_qq2.derive_cov(flat) ;
      dcovdcov_qq.inc2_dzpuis() ;

      Tenseur laplacien_qq(mp) ;
      laplacien_qq = qq().laplacien() ;
  
      tkij_auto.dec2_dzpuis() ;
      Tenseur dcov_tkij_auto = tkij_auto.derive_cov(flat) ;
      tkij_auto.inc2_dzpuis() ;
      dcov_tkij_auto.inc2_dzpuis() ;


      Cmp sol_poisson(mp) ;
      sol_poisson.set_etat_qcq() ;
      sol_poisson = 0 ;


      const Tenseur& gtilde_cov = gtilde.cov() ;
      const Tenseur& gtilde_con = gtilde.con() ;
 

    
      Tenseur source1_1 = contract(contract(hij,0, dcovdcov_hij_auto,0),0,1) ;
 
      Tenseur source2_1 = contract(contract(dcov_hij,0, dcov_hij_auto,2), 1, 2) ;
      Tenseur source2_2 = contract(contract(gtilde_cov, 0, 
					    contract(contract(gtilde_con,0,dcov_hij,0),0,dcov_hij_auto, 0), 1), 0, 3) ;
      Tenseur source2_3 = contract(gtilde_con, 1, contract(contract(
   contract(gtilde_cov,0,dcov_hij,1),0,dcov_hij_auto, 2),1, 2), 0) ;
    Tenseur source2_4 = contract(gtilde_con, 1, contract(contract(
   contract(gtilde_cov,1,dcov_hij,2),0,dcov_hij_auto, 2),1, 2), 0) ; 
     Tenseur source2_5 = contract(gtilde_con,1, contract(gtilde_con,1,
        contract(contract(dcov_gtilde,1,dcov_hij_auto,1), 1, 3)*0.5 
        + pow(gamma, -2./3.)*dcov_acar*dcov_acar_auto, 1), 1) ;


    Tenseur source3_1 = 2.*(pow(nnn, -1.)*pow(gamma,-1./6.)*contract(
       gtilde_con,1, contract(gtilde_auto.con(),1, dcovdcov_qq,1), 1)) ;


    Tenseur source4_1 = contract(contract(gtilde_con,1, dcov_hij_auto,0),2, 
				 dcov_qq,0) ;
    Tenseur source4_2 = contract(contract(gtilde_con,1, dcov_hij_auto,0),0, 
				 dcov_qq,0) ;
    Tenseur source4_3 = contract(gtilde_con,1,
	      	 contract(gtilde_con,1, dcov_acar*dcov_logn_auto,1),1) ;
    Tenseur source4_4 = contract(gtilde_con,1,
		 contract(gtilde_con,1, dcov_acar*dcov_logn_auto,0),1) ;


    Tenseur source5_1 = -2./3.*pow(gamma, -1./6.)*pow(nnn, -1.)
	*laplacien_qq*gtilde_auto.con() ;


    Tenseur source6_1 = -2./3.*pow(gamma, -1./6.)*pow(nnn, -1.)*contract
	(contract(hij,0, dcovdcov_qq,0),0,1)*gtilde_auto.con() ;


    Cmp source7_1 = - 2.*pow(a_car, -1.)()
	*contract(contract(gtilde_con,0, dcov_acar,0),0, dcov_logn,0)() ;
    Cmp source7_2 = - contract(contract(contract(contract(gtilde_con,0, 
      dcov_hij,0),0, dcov_gtilde,0), 0, 2), 0, 1)() *0.25 ;
    Cmp source7_3 = contract(contract(contract(contract(gtilde_con,0, 
        dcov_hij,0),0, dcov_gtilde,2), 0, 3), 0, 1)() *0.5 ;   
    Cmp source7_4 = - pow(gamma, -2./3.)()*contract(contract(gtilde_con,0,
        dcov_acar,0),0, dcov_acar,0)()*0.5 ;
    

    Tenseur source8_1 = -2.*pow(gamma,1./3.)*(2*kcar_con
		       - 2.*qpig*(pow(gamma, 1./3.)*stress
		       - s_euler*(gtilde.con())/3.)) ;


    Tenseur source9_1 = contract(shift,0, dcov_tkij_auto,0) ;
    Tenseur source9_2 = - contract(tkij_auto,0, shift.derive_cov(flat),0) ;
    Tenseur source9_3 = - contract(tkij_auto,1, shift.derive_cov(flat),0) ;
    Tenseur source9_4 = 2./3.*(contract(shift.derive_cov(flat),0, 1)
			       *tkij_auto) ;
			       

 
    hij_auto.set_etat_qcq() ;    

      for(int i=0; i<=2; i++) {
	for(int j=i; j<=2; j++) {

	    source1 = - source1_1(i,j);

	    source2 =  source2_1(i,j) + source2_2(i,j) - source2_3(i,j) 
	    - source2_4(j,i) ;//- source2_5(i,j) ;
	  
	    source3 = source3_1(i,j) ;
		      
	    source4 = 2.*pow(nnn, -1.)()*pow(gamma,-1./6.)()*
    (source4_1(i,j)*0.5 + source4_1(j,i)*0.5 - source4_2(i,j)*0.5)
	  - 2*pow(gamma,-1./3.)()*(source4_3(i,j) + source4_4(i,j)) ; 
		      
	  source5 = source5_1(i,j) ;
		                       
	  source6 = source6_1(i,j) ;
		                       
	  source7 =  -2./3.*(source7_1 )//+ source7_2 + source7_3 + source7_4)
	                   *gtilde_auto.con()(i,j) ;			   
		      
	  source8 =  source8_1(i,j) ;
		       
  	  source9 = - 2.*a_car()/nnn()*(
          source9_1(i,j)+ source9_2(j,i)+ source9_3(i,j) + source9_4(i,j)); 
		      

	  source1.annule(nz-1) ;
	  source2.annule(nz-1) ;
	  source3.annule(nz-1) ;
	  source4.annule(nz-1) ;
	  source5.annule(nz-1) ;
	  source6.annule(nz-1) ;
	  source7.annule(nz-1) ;
	  source8.annule(nz-1) ;
	  source9.annule(nz-1) ;


	  source_tot.set_etat_qcq() ;
	  source_tot.set() = source1() + source2() + source3() + source4() 
	      + source5() + source6() + source7() + source8() + source9() ;
 

	  //  source_tot.set_std_base() ;
	  
	  //source_tot.set().filtre(4) ;
//	  	  source_tot.annule(nz-1) ;
/*
	  source3.set_std_base() ; 
	  source4.set_std_base() ; 
	  des_profile(source_tot(), 0, 20, 0, 0) ;
	  des_coef_xi(source_tot().va, 0, 0, 0) ;
	  des_coef_xi(source_tot().va, 1, 0, 0) ;
	  des_coef_xi(source_tot().va, 2, 0, 0) ;
*/	  
  
	  cout << "moyenne de la source 1 pour hij_auto" << endl ;
	  cout <<  norme(source1()/(nr*nt*np)) << endl ;
	  cout << "moyenne de la source 2 pour hij_auto" << endl ;
	  cout <<  norme(source2()/(nr*nt*np)) << endl ;
	  cout << "moyenne de la source 3 pour hij_auto" << endl ;
	  cout <<  norme(source3()/(nr*nt*np)) << endl ;
	  cout << "moyenne de la source 4 pour hij_auto" << endl ;
	  cout <<  norme(source4()/(nr*nt*np)) << endl ;
	  cout << "moyenne de la source5 pour hij_auto" << endl ;
	  cout <<  norme(source5()/(nr*nt*np)) << endl ;
	  cout << "moyenne de la source6 pour hij_auto" << endl ;
	  cout <<  norme(source6()/(nr*nt*np)) << endl ;
	  cout << "moyenne de la source7 pour hij_auto" << endl ;
	  cout <<  norme(source7()/(nr*nt*np)) << endl ;
	  cout << "moyenne de la source8 pour hij_auto" << endl ;
	  cout <<  norme(source8()/(nr*nt*np)) << endl ;
	  cout << "moyenne de la source9 pour hij_auto" << endl ;
	  cout <<  norme(source9()/(nr*nt*np)) << endl ;
	  cout << "moyenne de la source pour hij_auto" << endl ;
	  cout <<  norme(source_tot()/(nr*nt*np)) << endl << endl ;
  

	  // Resolution of the Poisson equations and
	  // Check: has the Poisson equation been correctly solved ?
	  // -----------------------------------------------------

     
	  if(i==0 && j==0) {
		 
	    source_tot.set_std_base() ;
	    source_tot().poisson(par_poisson3, sol_poisson) ; 
	    hij_auto.set(0,0) = sol_poisson ;

	    gtilde_auto.set_con(0,0) = hij_auto(0,0) + decouple ;
	    gtilde_auto_con.set(0,0) = hij_auto(0,0) + decouple ;
	    
	    Tbl tdiff_hij00 = diffrel(hij_auto(0,0).laplacien(), source_tot()) ;  
	    cout << "Relative error in the resolution of the equation for "
		 << "hij00_auto : " << endl ; 
	    for (int l=0; l<nz; l++) {
	      cout << tdiff_hij00(l) << "  " ; 
	    }
	    cout << endl ;
	    diff_hij00 = max(abs(tdiff_hij00)) ; 
	  }
	       	       
	  if(i==0 && j==1) {

	    source_tot.set_std_base() ;
	    source_tot().poisson(par_poisson4, sol_poisson) ; 
	    hij_auto.set(0,1) = sol_poisson ;
	    
	    gtilde_auto.set_con(0,1) = hij_auto(0,1) ;
	    gtilde_auto_con.set(0,1) = hij_auto(0,1) ;

	    Tbl tdiff_hij10 = diffrel(hij_auto(0,1).laplacien(), source_tot()) ;
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
		 
	    source_tot.set_std_base() ;
	    source_tot().poisson(par_poisson5, sol_poisson) ; 
	    hij_auto.set(0,2) = sol_poisson ;

	    gtilde_auto.set_con(0,2) = hij_auto(0,2) ;
	    gtilde_auto_con.set(0,2) = hij_auto(0,2) ;
		 
	    Tbl tdiff_hij20 = diffrel(hij_auto(0,2).laplacien(), source_tot()) ;
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
		 
	    source_tot.set_std_base() ;
	    source_tot().poisson(par_poisson6, sol_poisson) ; 
	    hij_auto.set(1,1) = sol_poisson ;

	    gtilde_auto.set_con(1,1) = hij_auto(1,1) + decouple ;
	    gtilde_auto_con.set(1,1) = hij_auto(1,1) + decouple ;
		 
	    Tbl tdiff_hij11 = diffrel(hij_auto(1,1).laplacien(), source_tot()) ;
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
		 
	    source_tot.set_std_base() ;
	    source_tot().poisson(par_poisson7, sol_poisson) ; 
	    hij_auto.set(1,2) = sol_poisson ;

	    gtilde_auto.set_con(1,2) = hij_auto(1,2) ;
	    gtilde_auto_con.set(1,2) = hij_auto(1,2) ;

	    Tbl tdiff_hij21 = diffrel(hij_auto(1,2).laplacien(), source_tot()) ;
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
		 
	    source_tot.set_std_base() ;
	    source_tot().poisson(par_poisson8, sol_poisson) ; 
	    hij_auto.set(2,2) = sol_poisson ;

	    gtilde_auto.set_con(2,2) = hij_auto(2,2) + decouple ;
	    gtilde_auto_con.set(2,2) = hij_auto(2,2) + decouple ;
	 
	    Tbl tdiff_hij22 = diffrel(hij_auto(2,2).laplacien(), source_tot()) ;
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
      }
      
      cout << "gtilde00 auto"<<endl<< norme(gtilde_auto.con()(0,0)/(nr*nt*np)) 
	   << endl << endl ;
      cout << "gtilde10 auto"<<endl<< norme(gtilde_auto.con()(0,1)/(nr*nt*np)) 
	   << endl << endl ;
      cout << "gtilde20 auto"<<endl<< norme(gtilde_auto.con()(0,2)/(nr*nt*np)) 
	   << endl << endl ;
      cout << "gtilde11 auto"<<endl<< norme(gtilde_auto.con()(1,1)/(nr*nt*np)) 
	   << endl << endl ;
      cout << "gtilde21 auto"<<endl<< norme(gtilde_auto.con()(1,2)/(nr*nt*np)) 
	   << endl << endl ;
      cout << "gtilde22 auto"<<endl<< norme(gtilde_auto.con()(2,2)/(nr*nt*np)) 
	   << endl << endl ;

      //gtilde_auto.set_std_base() ;
     
//##       double nr = (*(mp.get_mg())).get_nr(0) ;
//##       double nt = (*(mp.get_mg())).get_nt(0) ;
//##       double np = (*(mp.get_mg())).get_np(0) ;
	   
      // Only valuable if nr is the same in all domains

    }
    else {
      for(int i=0; i<=2; i++) 
	for(int j=i; j<=2; j++) {
	  if (i == j) {
	    gtilde_auto_con.set(i,i) = decouple   ; 
	  }
	  else {
	 gtilde_auto_con.set(i,j) = 0. ;
	  }
	}
    }

    /*
    des_profile(gtilde_auto.cov()(0,0), 0, 20, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(0,0).va, 0, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(0,0).va, 1, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(0,0).va, 2, 0, 0) ;

    des_profile(gtilde_auto.cov()(1,0), 0, 20, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(1,0).va, 0, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(1,0).va, 1, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(1,0).va, 2, 0, 0) ;

    des_profile(gtilde_auto.cov()(2,0), 0, 20, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(2,0).va, 0, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(2,0).va, 1, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(2,0).va, 2, 0, 0) ;

    des_profile(gtilde_auto.cov()(1,1), 0, 20, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(1,1).va, 0, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(1,1).va, 1, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(1,1).va, 2, 0, 0) ;

    des_profile(gtilde_auto.cov()(2,1), 0, 20, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(2,1).va, 0, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(2,1).va, 1, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(2,1).va, 2, 0, 0) ;

    des_profile(gtilde_auto.cov()(2,2), 0, 20, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(2,2).va, 0, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(2,2).va, 1, 0, 0) ;
    des_coef_xi(gtilde_auto.cov()(2,2).va, 2, 0, 0) ;
    */

    // End of relativistic equations	
	   
	   
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



