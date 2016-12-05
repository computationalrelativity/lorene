/*
 * Method Star_rot::equilibrium
 *
 * (see file star_h.h for documentation)
 *
 */

/*
 *   Copyright (c) 2010 Eric Gourgoulhon
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2016/12/05 16:18:15  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:53:39  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:13:17  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2010/01/26 16:49:53  e_gourgoulhon
 * Reformated some outputs to the screen.
 *
 * Revision 1.3  2010/01/26 10:49:43  e_gourgoulhon
 * First working version.
 *
 * Revision 1.2  2010/01/25 22:32:38  e_gourgoulhon
 * Start of real implementation
 *
 * Revision 1.1  2010/01/25 18:15:52  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene
#include "star_rot.h"
#include "param.h"
#include "tenseur.h"

#include "graphique.h"
#include "utilitaires.h"
#include "unites.h"

namespace Lorene {
void Star_rot::equilibrium(double ent_c, double omega0, double fact_omega, 
			     int nzadapt, const Tbl& ent_limit, const Itbl& icontrol,
			     const Tbl& control, double mbar_wanted, 
			     double aexp_mass, Tbl& diff, Param*) {
			     
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
    int nzm1 = nz - 1 ; 
    
    // The following is required to initialize mp_prev as a Map_et:
    Map_et& mp_et = dynamic_cast<Map_et&>(mp) ; 
        
    // Index of the point at phi=0, theta=pi/2 at the surface of the star:
    assert(mg->get_type_t() == SYM) ; 
    int l_b = nzet - 1 ; 
    int i_b = mg->get_nr(l_b) - 1 ; 
    int j_b = mg->get_nt(l_b) - 1 ; 
    int k_b = 0 ; 
    
    // Value of the enthalpy defining the surface of the star
    double ent_b = ent_limit(nzet-1) ;
    
    // Parameters to control the iteration
    // -----------------------------------
    
    int mer_max = icontrol(0) ; 
    int mer_rot = icontrol(1) ;
    int mer_change_omega = icontrol(2) ; 
    int mer_fix_omega = icontrol(3) ; 
    int mer_mass = icontrol(4) ; 
    int mermax_poisson = icontrol(5) ; 
    int mer_triax = icontrol(6) ; 
    int delta_mer_kep = icontrol(7) ; 

    // Protections:
    if (mer_change_omega < mer_rot) {
	cout << "Star_rot::equilibrium: mer_change_omega < mer_rot !" << endl ;
	cout << " mer_change_omega = " << mer_change_omega << endl ; 
	cout << " mer_rot = " << mer_rot << endl ; 
	abort() ; 
    }
    if (mer_fix_omega < mer_change_omega) {
	cout << "Star_rot::equilibrium: mer_fix_omega < mer_change_omega !" 
	     << endl ;
	cout << " mer_fix_omega = " << mer_fix_omega << endl ; 
	cout << " mer_change_omega = " << mer_change_omega << endl ; 
	abort() ; 
    }

    // In order to converge to a given baryon mass, shall the central
    // enthalpy be varied or Omega ?
    bool change_ent = true ; 
    if (mer_mass < 0) {
	change_ent = false ; 
	mer_mass = abs(mer_mass) ;
    }

    double precis = control(0) ; 
    double omega_ini = control(1) ; 
    double relax = control(2) ;
    double relax_prev = double(1) - relax ;  
    double relax_poisson = control(3) ; 
    double thres_adapt = control(4) ; 
    double ampli_triax = control(5) ; 
    double precis_adapt = control(6) ; 


    // Error indicators
    // ----------------
    
    diff.set_etat_qcq() ; 
    double& diff_ent = diff.set(0) ; 
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
    	   
    par_adapt.add_tbl(ent_limit, 0) ;	// array of values of the field ent 
				        // to define the isosurfaces. 			   

    // Parameters for the function Map_et::poisson for nuf
    // ----------------------------------------------------

    double precis_poisson = 1.e-16 ;     

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

    // Initial angular velocity
    omega = 0 ; 
  
    double accrois_omega = (omega0 - omega_ini) / 
			    double(mer_fix_omega - mer_change_omega) ; 


    update_metric() ;	// update of the metric coefficients

    equation_of_state() ;	// update of the density, pressure, etc...
    
    hydro_euler() ;	// update of the hydro quantities relative to the 
			//  Eulerian observer

    // Quantities at the previous step : 	
    Map_et mp_prev = mp_et ; 
    Scalar ent_prev = ent ;	    
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
    Scalar mlngamma(mp) ;	// centrifugal potential
    
		    
    ofstream fichconv("convergence.d") ;    // Output file for diff_ent
    fichconv << "#     diff_ent     GRV2       max_triax      vit_triax" << endl ; 
    
    ofstream fichfreq("frequency.d") ;    // Output file for  omega
    fichfreq << "#       f [Hz]" << endl ; 
    
    ofstream fichevol("evolution.d") ;    // Output file for various quantities
    fichevol << 
    "#       |dH/dr_eq/dH/dr_pole|      r_pole/r_eq	ent_c" 
    << endl ; 
    
    diff_ent = 1 ; 
    double err_grv2 = 1 ; 
    double max_triax_prev = 0 ;	    // Triaxial amplitude at previous step

    //=========================================================================
    // 			Start of iteration
    //=========================================================================

    for(int mer=0 ; (diff_ent > precis) && (mer<mer_max) ; mer++ ) {

 	cout << "-----------------------------------------------" << endl ;
	cout << "step: " << mer << endl ;
	cout << "diff_ent = " << display_bold << diff_ent << display_normal
	     << endl ;    
	cout << "err_grv2 = " << err_grv2 << endl ;    
	fichconv << mer ;
	fichfreq << mer ;
	fichevol << mer ;

	if (mer >= mer_rot) {
	
	    if (mer < mer_change_omega) {
		omega = omega_ini ; 
	    }
	    else {
		if (mer <= mer_fix_omega) {
		    omega = omega_ini + accrois_omega * 
			      (mer - mer_change_omega) ;
		}
	    }

	}

	//-----------------------------------------------
	//  Sources of the Poisson equations
	//-----------------------------------------------
	
	// Source for nu
	// -------------
	Scalar bet = log(bbb) ; 
	bet.std_spectral_base() ; 

	Vector d_logn = logn.derive_cov( mp.flat_met_spher() ) ; 
	Vector d_bet = bet.derive_cov( mp.flat_met_spher() ) ; 

	if (relativistic) {
	    source_nuf =  qpig * a_car *( ener_euler + s_euler ) ; 

	    source_nuq = ak_car - d_logn(1)*(d_logn(1)+d_bet(1))
		- d_logn(2)*(d_logn(2)+d_bet(2))
		- d_logn(3)*(d_logn(3)+d_bet(3)) ; 
	}
	else {
	    source_nuf = qpig * nbar ; 

	    source_nuq = 0 ; 
	}
	source_nuf.std_spectral_base() ; 	
	source_nuq.std_spectral_base() ; 	

	// Source for dzeta
	// ----------------
	source_dzf = 2 * qpig * a_car * (press + (ener_euler+press) * uuu*uuu ) ;
	source_dzf.std_spectral_base() ; 
  
	source_dzq = 1.5 * ak_car 
	    - d_logn(1)*d_logn(1) - d_logn(2)*d_logn(2) - d_logn(3)*d_logn(3) ;	    
	source_dzq.std_spectral_base() ; 	
	
	// Source for tggg
	// ---------------
	
	source_tggg = 4 * qpig * nn * a_car * bbb * press ;
	source_tggg.std_spectral_base() ; 
	
	source_tggg.mult_rsint() ; 
	

	// Source for shift
	// ----------------
	    
	// Matter term: 
	source_shift = (-4*qpig) * nn * a_car  * (ener_euler + press)
				* u_euler ;

	// Quadratic terms:
	Vector vtmp =  6 * bet.derive_con( mp.flat_met_spher() ) 
	    - 2 * logn.derive_con( mp.flat_met_spher() ) ;
	vtmp.change_triad(mp.get_bvect_cart()) ; 

	Vector squad  = nn * contract(tkij, 1, vtmp, 0) ; 

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
		"Star_rot::equilibrium: np must be stricly greater than 1"
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

	if (relativistic) {
	    
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
	
	}

	//-----------------------------------------
	// Determination of the fluid velociy U
	//-----------------------------------------
	
	if (mer > mer_fix_omega + delta_mer_kep) {
		
    	    omega *= fact_omega ;  // Increase of the angular velocity if 
	}			   //  fact_omega != 1
	
	bool omega_trop_grand = false ; 
	bool kepler = true ; 

	while ( kepler ) {
	
	    // Possible decrease of Omega to ensure a velocity < c 
	
	    bool superlum = true ; 
	
	    while ( superlum ) {
	
		// New fluid velocity U :
	
		Scalar tmp = omega - nphi ; 
		tmp.annule_domain(nzm1) ; 
		tmp.std_spectral_base() ;
    
		tmp.mult_rsint() ;	    //  Multiplication by r sin(theta)

		uuu = bbb / nn * tmp ; 
    
		if (uuu.get_etat() == ETATQCQ) {
		    // Same basis as (Omega -N^phi) r sin(theta) :
		    (uuu.set_spectral_va()).set_base( tmp.get_spectral_va().get_base() ) ;   
		}

		// Is the new velocity larger than c in the equatorial plane ?
	
		superlum = false ; 
	
		for (int l=0; l<nzet; l++) {
		    for (int i=0; i<mg->get_nr(l); i++) {
		
			double u1 = uuu.val_grid_point(l, 0, j_b, i) ; 
			if (u1 >= 1.) {	    // superluminal velocity
			    superlum = true ; 
			    cout << "U > c  for l, i : " << l << "  " << i 
				 << "   U = " << u1 << endl ;  
			}
		    }
		}
		if ( superlum ) {
		    cout << "**** VELOCITY OF LIGHT REACHED ****" << endl ; 
		    omega /= fact_omega ;    // Decrease of Omega
		    cout << "New rotation frequency : " 
			 << omega/(2.*M_PI) * f_unit <<  " Hz" << endl ; 
		    omega_trop_grand = true ;  
		}
	    }	// end of while ( superlum )

	
	    // New computation of U (which this time is not superluminal)
	    //  as well as of gam_euler, ener_euler, etc...
	    // -----------------------------------

	    hydro_euler() ; 
	

	    //------------------------------------------------------
	    //	First integral of motion 
	    //------------------------------------------------------
	
	    // Centrifugal potential : 
	    if (relativistic) {
		mlngamma = - log( gam_euler ) ;
	    }
	    else {
		mlngamma = - 0.5 * uuu*uuu ; 
	    }
	
	    // Equatorial values of various potentials :
	    double nuf_b  = nuf.val_grid_point(l_b, k_b, j_b, i_b) ; 
	    double nuq_b  = nuq.val_grid_point(l_b, k_b, j_b, i_b) ; 
	    double mlngamma_b  = mlngamma.val_grid_point(l_b, k_b, j_b, i_b) ; 

	    // Central values of various potentials :
	    double nuf_c = nuf.val_grid_point(0,0,0,0) ; 
	    double nuq_c = nuq.val_grid_point(0,0,0,0) ; 
	    double mlngamma_c = 0 ;
	
	    // Scale factor to ensure that the enthalpy is equal to ent_b at 
	    //  the equator
	    double alpha_r2 = ( ent_c - ent_b + mlngamma_c - mlngamma_b
				+ nuq_c - nuq_b) / ( nuf_b - nuf_c  ) ;
	    alpha_r = sqrt(alpha_r2) ;
	    cout << "alpha_r = " << alpha_r << endl ; 

	    // Readjustment of nu :
	    // -------------------

	    logn = alpha_r2 * nuf + nuq ;
	    double nu_c =  logn.val_grid_point(0,0,0,0) ;

	    // First integral	--> enthalpy in all space
	    //-----------------

	    ent = (ent_c + nu_c + mlngamma_c) - logn - mlngamma ;

	    // Test: is the enthalpy negative somewhere in the equatorial plane
	    //  inside the star ? If yes, this means that the Keplerian velocity
	    //  has been overstep.

	    kepler = false ; 
	    for (int l=0; l<nzet; l++) {
		int imax = mg->get_nr(l) - 1 ;
		if (l == l_b) imax-- ;	// The surface point is skipped
		for (int i=0; i<imax; i++) { 
		    if ( ent.val_grid_point(l, 0, j_b, i) < 0. ) {
			kepler = true ;
			cout << "ent < 0 for l, i : " << l << "  " << i 
			     << "   ent = " << ent.val_grid_point(l, 0, j_b, i) << endl ;  
		    } 
		}
	    }
	
	    if ( kepler ) {
		cout << "**** KEPLERIAN VELOCITY REACHED ****" << endl ; 
		omega /= fact_omega ;    // Omega is decreased
		cout << "New rotation frequency : " 
			 << omega/(2.*M_PI) * f_unit << " Hz" << endl ; 
		omega_trop_grand = true ;  
	    }

	}   // End of while ( kepler )
	
	if ( omega_trop_grand ) {	// fact_omega is decreased for the
					//  next step 
	    fact_omega = sqrt( fact_omega ) ; 
	    cout << "**** New fact_omega : " << fact_omega << endl ; 
	}

	//----------------------------------------------------
	// Adaptation of the mapping to the new enthalpy field
	//----------------------------------------------------
    
	// Shall the adaptation be performed (cusp) ?
	// ------------------------------------------
	
	double dent_eq = ent.dsdr().val_grid_point(l_b, k_b, j_b, i_b) ; 
	double dent_pole = ent.dsdr().val_grid_point(l_b, k_b, 0, i_b) ;
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

	mp_prev = mp_et ; 

	Cmp cent(ent) ; 
	mp.adapt(cent, par_adapt) ; 

 	//----------------------------------------------------
	// Computation of the enthalpy at the new grid points
	//----------------------------------------------------
	
	mp_prev.homothetie(alpha_r) ; 
	
	mp.reevaluate(&mp_prev, nzet+1, cent) ;
	ent = cent ; 

	//----------------------------------------------------
	// Equation of state  
	//----------------------------------------------------
	
	equation_of_state() ; 	// computes new values for nbar (n), ener (e) 
				// and press (p) from the new ent (H)
	
	//---------------------------------------------------------
	// Matter source terms in the gravitational field equations	
	//---------------------------------------------------------

	//## Computation of tnphi and nphi from the Cartesian components
	//  of the shift for the test in hydro_euler():
	    
        fait_nphi() ; 

	hydro_euler() ;		// computes new values for ener_euler (E), 
				// s_euler (S) and u_euler (U^i)

	if (relativistic) {

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
	    
	}
	else {
		err_grv2 = grv2() ; 
	}


	//---------------------------------------
	// Computation of the metric coefficients (except for N^phi)
	//---------------------------------------

	// Relaxations on nu and dzeta :  

	if (mer >= 10) {
		logn = relax * logn + relax_prev * logn_prev ;

		dzeta = relax * dzeta + relax_prev * dzeta_prev ; 
	}

	// Update of the metric coefficients N, A, B and computation of K_ij :

	update_metric() ; 
	
	//-----------------------
	//  Informations display
	//-----------------------

	partial_display(cout) ; 
	fichfreq << "  " << omega / (2*M_PI) * f_unit ; 
	fichevol << "  " << rap_dent ; 
	fichevol << "  " << ray_pole() / ray_eq() ; 
	fichevol << "  " << ent_c ; 

	//-----------------------------------------
	// Convergence towards a given baryon mass 
	//-----------------------------------------

	if (mer > mer_mass) {
	    
	    double xx ; 
	    if (mbar_wanted > 0.) {
		xx = mass_b() / mbar_wanted - 1. ;
		cout << "Discrep. baryon mass <-> wanted bar. mass : " << xx 
		     << endl ; 
	    }
	    else{
		xx = mass_g() / fabs(mbar_wanted) - 1. ;
		cout << "Discrep. grav. mass <-> wanted grav. mass : " << xx 
		     << endl ; 
	    }
	    double xprog = ( mer > 2*mer_mass) ? 1. : 
				 double(mer-mer_mass)/double(mer_mass) ; 
	    xx *= xprog ; 
	    double ax = .5 * ( 2. + xx ) / (1. + xx ) ; 
	    double fact = pow(ax, aexp_mass) ; 
	    cout << "  xprog, xx, ax, fact : " << xprog << "  " <<
			xx << "  " << ax << "  " << fact << endl ; 

	    if ( change_ent ) {
		ent_c *= fact ; 
	    }
	    else {
		if (mer%4 == 0) omega *= fact ; 
	    }
	}
	

	//------------------------------------------------------------
	//  Relative change in enthalpy with respect to previous step 
	//------------------------------------------------------------

	Tbl diff_ent_tbl = diffrel( ent, ent_prev ) ; 
	diff_ent = diff_ent_tbl(0) ; 
	for (int l=1; l<nzet; l++) {
	    diff_ent += diff_ent_tbl(l) ; 
	}
	diff_ent /= nzet ; 
	
	fichconv << "  " << log10( fabs(diff_ent) + 1.e-16 ) ;
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
	
	ent_prev = ent ; 
	logn_prev = logn ; 
	dzeta_prev = dzeta ; 
	max_triax_prev = max_triax ; 
	
	fichconv << endl ;
	fichfreq << endl ;
	fichevol << endl ;
	fichconv.flush() ; 
	fichfreq.flush() ; 
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

    fichconv.close() ; 
    fichfreq.close() ; 
    fichevol.close() ; 
    
}
}
