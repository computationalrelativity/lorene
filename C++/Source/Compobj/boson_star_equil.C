/*
 * Method Boson_star::equilibrium
 *
 *    (see file boson_star.h for documentation).
 */

/*
 *   Copyright (c) 2012 Claire Some, Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


char boson_star_equil_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/06 15:13:04  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2013/04/03 12:10:13  e_gourgoulhon
 * Added member kk to Compobj; suppressed tkij
 *
 * Revision 1.3  2012/12/03 15:27:30  c_some
 * Small changes
 *
 * Revision 1.2  2012/11/23 15:43:05  c_some
 * Small changes
 *
 * Revision 1.1  2012/11/22 16:04:12  c_some
 * New class Boson_star
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene
#include "boson_star.h"
#include "param.h"
#include "tenseur.h"

#include "graphique.h"
#include "utilitaires.h"
#include "unites.h"

void Boson_star::equilibrium(double, double,
			     int nzadapt, const Tbl& phi_limit, const Itbl& icontrol,
			     const Tbl& control, Tbl& diff, Param*) {
			     
    // Fundamental constants and units
    // -------------------------------

    using namespace Unites ;	    
    
    // For the display 
    // ---------------
    char display_bold[]="x[1m" ; display_bold[0] = 27 ;
    char display_normal[] = "x[0m" ; display_normal[0] = 27 ;

    // Grid parameters
    // ---------------
    
    const Mg3d* mg = mp.get_mg() ; 
    int nz = mg->get_nzone() ;	    // total number of domains
    
    // The following is required to initialize mp_prev as a Map_et:
    Map_et& mp_et = dynamic_cast<Map_et&>(mp) ; 
        
    // Index of the point at phi=0, theta=pi/2 at the surface of the star:

    int nzet = nzadapt ; //## to be checked 
    
    assert(mg->get_type_t() == SYM) ; 
    int l_b = nzet - 1 ; 
    int j_b = mg->get_nt(l_b) - 1 ; 
    int k_b = 0 ; 
    
    // Value of the enthalpy defining the surface of the star
    // double ent_b = phi_limit(nzet-1) ;
    
    // Parameters to control the iteration
    // -----------------------------------
    
    int mer_max = icontrol(0) ; 
//    int mer_rot = icontrol(1) ;
//    int mer_change_omega = icontrol(2) ; 
//    int mer_fix_omega = icontrol(3) ; 
//##    int mer_mass = icontrol(4) ; 
    int mermax_poisson = icontrol(5) ; 
    int mer_triax = icontrol(6) ; 
//##    int delta_mer_kep = icontrol(7) ; 


    double precis = control(0) ; 
//    double omega_ini = control(1) ; 
    double relax = control(2) ;
    double relax_prev = double(1) - relax ;  
    double relax_poisson = control(3) ; 
//    double thres_adapt = control(4) ; 
    double ampli_triax = control(5) ; 
    double precis_adapt = control(6) ; 


    // Error indicators
    // ----------------
    
    diff.set_etat_qcq() ; 
    double& diff_phi = diff.set(0) ; 
    double& diff_nuf = diff.set(1) ; 
    double& diff_nuq = diff.set(2) ; 
//    double& diff_dzeta = diff.set(3) ; 
//    double& diff_ggg = diff.set(4) ; 
    double& diff_shift_x = diff.set(5) ; 
    double& diff_shift_y = diff.set(6) ; 
    double& vit_triax = diff.set(7) ; 
    
    // Parameters for the function Map_et::adapt
    // -----------------------------------------
    
    Param par_adapt ; 
    int nitermax = 100 ;  
    int niter ; 
    int adapt_flag = 1 ;    //  1 = performs the full computation, 
			    //  0 = performs only the rescaling by 
			    //      the factor alpha_r
    int nz_search = nzet + 1 ;  // Number of domains for searching the enthalpy
				//  isosurfaces
    double alpha_r ; 
    double reg_map = 1. ; // 1 = regular mapping, 0 = contracting mapping

    par_adapt.add_int(nitermax, 0) ; // maximum number of iterations to
				     // locate zeros by the secant method
    par_adapt.add_int(nzadapt, 1) ; // number of domains where the adjustment 
				    // to the isosurfaces of ent is to be 
				    // performed
    par_adapt.add_int(nz_search, 2) ;	// number of domains to search for
					// the enthalpy isosurface
    par_adapt.add_int(adapt_flag, 3) ; //  1 = performs the full computation, 
				       //  0 = performs only the rescaling by 
				       //      the factor alpha_r
    par_adapt.add_int(j_b, 4) ; //  theta index of the collocation point 
			        //  (theta_*, phi_*)
    par_adapt.add_int(k_b, 5) ; //  theta index of the collocation point 
			        //  (theta_*, phi_*)

    par_adapt.add_int_mod(niter, 0) ;  //  number of iterations actually used in 
				       //  the secant method
    
    par_adapt.add_double(precis_adapt, 0) ; // required absolute precision in 
					     // the determination of zeros by 
					     // the secant method
    par_adapt.add_double(reg_map, 1)	;  // 1. = regular mapping, 0 = contracting mapping
    
    par_adapt.add_double(alpha_r, 2) ;	// factor by which all the radial 
					   // distances will be multiplied 
    	   
    par_adapt.add_tbl(phi_limit, 0) ;	// array of values of the field Phi
				        // to define the isosurfaces. 			   

    // Parameters for the function Map_et::poisson for nuf
    // ----------------------------------------------------

    double precis_poisson = 1.e-16 ;     

	// Preparations 
	// ------------
	
	// Cartesian components of the shift vector are required
	beta.change_triad( mp.get_bvect_cart() ) ; 


    Cmp cssjm1_nuf(ssjm1_nuf) ; 
    Cmp cssjm1_nuq(ssjm1_nuq) ; 
    Cmp cssjm1_tggg(ssjm1_tggg) ; 
    Cmp cssjm1_khi(ssjm1_khi) ; 
    Tenseur cssjm1_wshift(mp, 1, CON, mp.get_bvect_cart() ) ;
    cssjm1_wshift.set_etat_qcq() ;
    for (int i=1; i<=3; i++) {
	cssjm1_wshift.set(i-1) = ssjm1_wshift(i) ; 
    }

    Tenseur cshift(mp, 1, CON, mp.get_bvect_cart() ) ;
    cshift.set_etat_qcq() ;
    for (int i=1; i<=3; i++) {
      cshift.set(i-1) = -beta(i) ; 
    }

    Tenseur cw_shift(mp, 1, CON, mp.get_bvect_cart() ) ;
    cw_shift.set_etat_qcq() ;
    for (int i=1; i<=3; i++) {
      cw_shift.set(i-1) = w_shift(i) ; 
    }

    Tenseur ckhi_shift(mp) ;
    ckhi_shift.set_etat_qcq() ;
    ckhi_shift.set() = khi_shift ; 

    Param par_poisson_nuf ; 
    par_poisson_nuf.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson_nuf.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson_nuf.add_double(precis_poisson, 1) ; // required precision
    par_poisson_nuf.add_int_mod(niter, 0) ;  //  number of iterations actually used 
    par_poisson_nuf.add_cmp_mod( cssjm1_nuf ) ; 
					   
    Param par_poisson_nuq ; 
    par_poisson_nuq.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson_nuq.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson_nuq.add_double(precis_poisson, 1) ; // required precision
    par_poisson_nuq.add_int_mod(niter, 0) ;  //  number of iterations actually used 
    par_poisson_nuq.add_cmp_mod( cssjm1_nuq ) ; 
					   
    Param par_poisson_tggg ; 
    par_poisson_tggg.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson_tggg.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson_tggg.add_double(precis_poisson, 1) ; // required precision
    par_poisson_tggg.add_int_mod(niter, 0) ;  //  number of iterations actually used 
    par_poisson_tggg.add_cmp_mod( cssjm1_tggg ) ; 
    double lambda_tggg ;
    par_poisson_tggg.add_double_mod( lambda_tggg ) ; 
    
    Param par_poisson_dzeta ; 
    double lbda_grv2 ;
    par_poisson_dzeta.add_double_mod( lbda_grv2 ) ; 
 					   
    // Parameters for the function Scalar::poisson_vect
    // -------------------------------------------------

    Param par_poisson_vect ; 

    par_poisson_vect.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson_vect.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson_vect.add_double(precis_poisson, 1) ; // required precision
    par_poisson_vect.add_cmp_mod( cssjm1_khi ) ; 
    par_poisson_vect.add_tenseur_mod( cssjm1_wshift ) ; 
    par_poisson_vect.add_int_mod(niter, 0) ;   

 					   
    // Initializations
    // ---------------

	//## Spherical components of the shift vector are restored
	beta.change_triad( mp.get_bvect_spher() ) ; 
	update_metric() ; 
	//## Back to Cartesian components 
	beta.change_triad( mp.get_bvect_cart() ) ; 
	

    // Quantities at the previous step : 	
    Map_et mp_prev = mp_et ; 
    Scalar rphi_prev = rphi ;	    
    Scalar logn_prev = logn ;	    
    Scalar dzeta_prev = dzeta ;	    
    
    // Creation of uninitialized tensors:
    Scalar source_nuf(mp) ;    // source term in the equation for nuf
    Scalar source_nuq(mp) ;    // source term in the equation for nuq
    Scalar source_dzf(mp) ;	// matter source term in the eq. for dzeta
    Scalar source_dzq(mp) ;	// quadratic source term in the eq. for dzeta
    Scalar source_tggg(mp) ;	// source term in the eq. for tggg
    Vector source_shift(mp, CON, mp.get_bvect_cart()) ;  
						    // source term for shift
    
		    
    ofstream fichconv("convergence.d") ;    // Output file for diff_phi
    fichconv << "#     diff_phi     GRV2       max_triax      vit_triax" << endl ; 
    
    
    ofstream fichevol("evolution.d") ;    // Output file for various quantities
    fichevol << 
    "#       |dH/dr_eq/dH/dr_pole|      r_pole/r_eq	rphi_c" 
    << endl ; 
    
    diff_phi = 1 ; 
    double err_grv2 = 1 ; 
    double max_triax_prev = 0 ;	    // Triaxial amplitude at previous step

    //=========================================================================
    // 			Start of iteration
    //=========================================================================

    for(int mer=0 ; (diff_phi > precis) && (mer<mer_max) ; mer++ ) {

 	cout << "-----------------------------------------------" << endl ;
	cout << "step: " << mer << endl ;
	cout << "diff_phi = " << display_bold << diff_phi << display_normal
	     << endl ;    
	cout << "err_grv2 = " << err_grv2 << endl ;    
	fichconv << mer ;
	fichevol << mer ;


	//-----------------------------------------------
	//  Sources of the Poisson equations
	//-----------------------------------------------
	
	// Source for nu
	// -------------
	Scalar bet = log(bbb) ; 
	bet.std_spectral_base() ; 

	Vector d_logn = logn.derive_cov( mp.flat_met_spher() ) ; 
	Vector d_bet = bet.derive_cov( mp.flat_met_spher() ) ; 

	Scalar s_euler = stress_euler.trace(gamma) ; 
	source_nuf =  qpig * a_car *( ener_euler + s_euler ) ; 

	source_nuq = ak_car - d_logn(1)*(d_logn(1)+d_bet(1))
		- d_logn(2)*(d_logn(2)+d_bet(2))
		- d_logn(3)*(d_logn(3)+d_bet(3)) ; 

	source_nuf.std_spectral_base() ; 	
	source_nuq.std_spectral_base() ; 	

	// Source for dzeta
	// ----------------
	source_dzf = 2 * qpig * a_car * b_car * stress_euler(3,3) ;
	source_dzf.std_spectral_base() ; 
  
	source_dzq = 1.5 * ak_car 
	    - d_logn(1)*d_logn(1) - d_logn(2)*d_logn(2) - d_logn(3)*d_logn(3) ;	    
	source_dzq.std_spectral_base() ; 	
	
	// Source for tggg
	// ---------------
	
	source_tggg = 2 * qpig * nn * a_car * bbb * (s_euler - b_car * stress_euler(3,3)) ;
	source_tggg.std_spectral_base() ; 
	
	source_tggg.mult_rsint() ; 
	

	// Source for shift
	// ----------------
	    
	// Matter term: 
	Vector mom_euler_cart = mom_euler ;
	mom_euler_cart.change_triad(mp.get_bvect_cart()) ;
	source_shift = (-4*qpig) * nn * a_car  * mom_euler_cart ;

	// Quadratic terms:
	Vector vtmp =  6 * bet.derive_con( mp.flat_met_spher() ) 
	    - 2 * logn.derive_con( mp.flat_met_spher() ) ;

	Vector squad  = nn * contract(kk, 1, vtmp, 0) / b_car ; 
    squad.change_triad(mp.get_bvect_cart()) ;
    
	source_shift = source_shift + squad.up(0, mp.flat_met_cart() ) ; 

	//----------------------------------------------
	// Resolution of the Poisson equation for nuf 
	//----------------------------------------------

	source_nuf.poisson(par_poisson_nuf, nuf) ; 
	
	cout << "Test of the Poisson equation for nuf :" << endl ; 
	Tbl err = source_nuf.test_poisson(nuf, cout, true) ; 
	diff_nuf = err(0, 0) ; 

	//---------------------------------------
	// Triaxial perturbation of nuf
	//---------------------------------------
	
	if (mer == mer_triax) {
	
	    if ( mg->get_np(0) == 1 ) {
		cout << 
		"Boson_star::equilibrium: np must be stricly greater than 1"
		<< endl << " to set a triaxial perturbation !" << endl ; 
		abort() ; 
	    }
	
	    const Coord& phi = mp.phi ; 
	    const Coord& sint = mp.sint ; 
	    Scalar perturb(mp) ; 
	    perturb = 1 + ampli_triax * sint*sint * cos(2*phi) ;
	    nuf = nuf * perturb ;  
	    
	    nuf.std_spectral_base() ;    // set the bases for spectral expansions
				    //  to be the standard ones for a 
				    //  scalar field

	}
	
	// Monitoring of the triaxial perturbation
	// ---------------------------------------
	
	const Valeur& va_nuf = nuf.get_spectral_va() ; 
	va_nuf.coef() ;		// Computes the spectral coefficients
	double max_triax = 0 ; 

	if ( mg->get_np(0) > 1 ) {

	    for (int l=0; l<nz; l++) {	// loop on the domains
		for (int j=0; j<mg->get_nt(l); j++) {
		    for (int i=0; i<mg->get_nr(l); i++) {
		
			// Coefficient of cos(2 phi) : 
			double xcos2p = (*(va_nuf.c_cf))(l, 2, j, i) ; 

			// Coefficient of sin(2 phi) : 
			double xsin2p = (*(va_nuf.c_cf))(l, 3, j, i) ; 
		    
			double xx = sqrt( xcos2p*xcos2p + xsin2p*xsin2p ) ; 
		    
			max_triax = ( xx > max_triax ) ? xx : max_triax ; 		    
		    }
		} 
	    }
	    
	}
	
	cout << "Triaxial part of nuf : " << max_triax << endl ; 

	    //----------------------------------------------
	    // Resolution of the Poisson equation for nuq 
	    //----------------------------------------------

	    source_nuq.poisson(par_poisson_nuq, nuq) ; 
	    
	    cout << "Test of the Poisson equation for nuq :" << endl ; 
	    err = source_nuq.test_poisson(nuq, cout, true) ;
	    diff_nuq = err(0, 0) ; 
	
	    //---------------------------------------------------------
	    // Resolution of the vector Poisson equation for the shift
	    //---------------------------------------------------------


	    for (int i=1; i<=3; i++) {
	      if(source_shift(i).get_etat() != ETATZERO) {
		    if(source_shift(i).dz_nonzero()) {
			assert( source_shift(i).get_dzpuis() == 4 ) ; 
		    }
		    else{
			(source_shift.set(i)).set_dzpuis(4) ; 
		    }
		}
	    }

	    double lambda_shift = double(1) / double(3) ; 

	    if ( mg->get_np(0) == 1 ) {
		lambda_shift = 0 ; 
	    }
	
	    Tenseur csource_shift(mp, 1, CON, mp.get_bvect_cart() ) ;
	    csource_shift.set_etat_qcq() ;
	    for (int i=1; i<=3; i++) {
		csource_shift.set(i-1) = source_shift(i) ; 
	    }
	    csource_shift.set(2).set_etat_zero() ;  //## bizarre...

	    csource_shift.poisson_vect(lambda_shift, par_poisson_vect, 
	    			      cshift, cw_shift, ckhi_shift) ;      	    

	    for (int i=1; i<=3; i++) {
		beta.set(i) = - cshift(i-1) ;
		beta.set(i).set_dzpuis(0) ;     //## bizarre...
		w_shift.set(i) = cw_shift(i-1) ; 
	    }
	    khi_shift = ckhi_shift() ; 

	    cout << "Test of the Poisson equation for shift_x :" << endl ; 
	    err = source_shift(1).test_poisson(-beta(1), cout, true) ;
	    diff_shift_x = err(0, 0) ; 
	
	    cout << "Test of the Poisson equation for shift_y :" << endl ; 
	    err = source_shift(2).test_poisson(-beta(2), cout, true) ;
	    diff_shift_y = err(0, 0) ; 
	    
	    // Computation of tnphi and nphi from the Cartesian components
	    //  of the shift
	    // -----------------------------------------------------------
	    
	    fait_nphi() ; 
	


	//----------------------------------------------------
	// Adaptation of the mapping to the new enthalpy field
	//----------------------------------------------------
    
	// Shall the adaptation be performed ?
	// ---------------------------------
	
	adapt_flag = 0 ;	// No adaptation of the mapping 

	mp_prev = mp_et ; 

	
	//---------------------------------------------------------
	// Matter source terms in the gravitational field equations	
	//---------------------------------------------------------

	//## Computation of tnphi and nphi from the Cartesian components
	//  of the shift for the test in hydro_euler():
	    
    fait_nphi() ; 

	// Update of the energy-momentum tensor
	update_ener_mom() ; 
	

	    //-------------------------------------------------------
	    //	2-D Poisson equation for tggg
	    //-------------------------------------------------------

	    Cmp csource_tggg(source_tggg) ; 
	    Cmp ctggg(tggg) ; 
	    mp.poisson2d(csource_tggg, mp.cmp_zero(), par_poisson_tggg,
			 ctggg) ; 
	    tggg = ctggg ; 

	    
	    //-------------------------------------------------------
	    //	2-D Poisson equation for dzeta
	    //-------------------------------------------------------

	    Cmp csource_dzf(source_dzf) ; 
	    Cmp csource_dzq(source_dzq) ; 
	    Cmp cdzeta(dzeta) ; 
	    mp.poisson2d(csource_dzf, csource_dzq, par_poisson_dzeta,
			 cdzeta) ;
	    dzeta = cdzeta ; 
	    
	    err_grv2 = lbda_grv2 - 1; 
	    cout << "GRV2: " << err_grv2 << endl ; 
	    

	//---------------------------------------
	// Computation of the metric coefficients (except for N^phi)
	//---------------------------------------

	// Relaxations on nu and dzeta :  

	if (mer >= 10) {
		logn = relax * logn + relax_prev * logn_prev ;

		dzeta = relax * dzeta + relax_prev * dzeta_prev ; 
	}

	// Update of the metric coefficients N, A, B and computation of K_ij :

	//## Spherical components of the shift vector are restored
	beta.change_triad( mp.get_bvect_spher() ) ; 
	update_metric() ; 
	//## Back to Cartesian components 
	beta.change_triad( mp.get_bvect_cart() ) ; 
	
	//-----------------------
	//  Informations display
	//-----------------------

	cout << *this << endl ; 


	//------------------------------------------------------------
	//  Relative change in Phi with respect to previous step 
	//------------------------------------------------------------

	Tbl diff_phi_tbl = diffrel( rphi, rphi_prev ) ; 
	diff_phi = diff_phi_tbl(0) ; 
	for (int l=1; l<nzet; l++) {
	    diff_phi += diff_phi_tbl(l) ; 
	}
	diff_phi /= nzet ; 
	
	fichconv << "  " << log10( fabs(diff_phi) + 1.e-16 ) ;
	fichconv << "  " << log10( fabs(err_grv2) + 1.e-16 ) ;
	fichconv << "  " << log10( fabs(max_triax) + 1.e-16 ) ;

	vit_triax = 0 ;
	if ( (mer > mer_triax+1) && (max_triax_prev > 1e-13) ) {
		vit_triax = (max_triax - max_triax_prev) / max_triax_prev ; 
	}

	fichconv << "  " << vit_triax ;
	
	//------------------------------
	//  Recycling for the next step
	//------------------------------
	
	rphi_prev = rphi ; 
	logn_prev = logn ; 
	dzeta_prev = dzeta ; 
	max_triax_prev = max_triax ; 
	
	fichconv << endl ;
	fichevol << endl ;
	fichconv.flush() ; 
	fichevol.flush() ; 
	
    } // End of main loop
    
    //=========================================================================
    // 			End of iteration
    //=========================================================================

    ssjm1_nuf = cssjm1_nuf ; 
    ssjm1_nuq = cssjm1_nuq ; 
    ssjm1_tggg = cssjm1_tggg ; 
    ssjm1_khi = cssjm1_khi ; 
    for (int i=1; i<=3; i++) {
		ssjm1_wshift.set(i) = cssjm1_wshift(i-1) ; 
    }

	// Spherical components of the shift vector are restored
	beta.change_triad( mp.get_bvect_spher() ) ; 

    fichconv.close() ; 
    fichevol.close() ; 
    
}
