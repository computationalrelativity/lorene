/*
 * Function Et_rot_diff::equilibrium
 *
 * (see file etoile.h for documentation)
 *
 */

/*
 *   Copyright (c) 2001-2003 Eric Gourgoulhon
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


char et_rot_diff_equil_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/10/27 10:54:43  e_gourgoulhon
 * Changed local variable name lambda_grv2 to lbda_grv2 in order not
 * to shadow method name.
 *
 * Revision 1.3  2003/05/25 19:59:02  e_gourgoulhon
 * Added the possibility to choose the factor a = R_eq / R0, instead of R0
 * in the differential rotation law.
 *
 * Revision 1.2  2002/10/16 14:36:36  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2001/10/19  08:18:16  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */



// Headers C
#include <math.h>

// Headers Lorene
#include "et_rot_diff.h"
#include "param.h"

#include "graphique.h"
#include "utilitaires.h"

void Et_rot_diff::equilibrium(double ent_c, double omega_c0, double fact_omega, 
			     int nzadapt, const Tbl& ent_limit, 
			     const Itbl& icontrol,
			     const Tbl& control, double mbar_wanted, 
			     double aexp_mass, Tbl& diff) {
			     
    // Fundamental constants and units
    // -------------------------------
    #include "unites.h"	    
    // To avoid some compilation warnings
    if (ent_c < 0) {
	cout << f_unit << msol << km << mevpfm3 << endl ; 
    }    
    
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
	cout << "Et_rot_diff::equilibrium: mer_change_omega < mer_rot !" << endl ;
	cout << " mer_change_omega = " << mer_change_omega << endl ; 
	cout << " mer_rot = " << mer_rot << endl ; 
	abort() ; 
    }
    if (mer_fix_omega < mer_change_omega) {
	cout << "Et_rot_diff::equilibrium: mer_fix_omega < mer_change_omega !" 
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

    Param par_poisson_nuf ; 
    par_poisson_nuf.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson_nuf.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson_nuf.add_double(precis_poisson, 1) ; // required precision
    par_poisson_nuf.add_int_mod(niter, 0) ;  //  number of iterations actually used 
    par_poisson_nuf.add_cmp_mod( ssjm1_nuf ) ; 
					   
    Param par_poisson_nuq ; 
    par_poisson_nuq.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson_nuq.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson_nuq.add_double(precis_poisson, 1) ; // required precision
    par_poisson_nuq.add_int_mod(niter, 0) ;  //  number of iterations actually used 
    par_poisson_nuq.add_cmp_mod( ssjm1_nuq ) ; 
					   
    Param par_poisson_tggg ; 
    par_poisson_tggg.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson_tggg.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson_tggg.add_double(precis_poisson, 1) ; // required precision
    par_poisson_tggg.add_int_mod(niter, 0) ;  //  number of iterations actually used 
    par_poisson_tggg.add_cmp_mod( ssjm1_tggg ) ; 
    double lambda_tggg ;
    par_poisson_tggg.add_double_mod( lambda_tggg ) ; 
    
    Param par_poisson_dzeta ; 
    double lbda_grv2 ;
    par_poisson_dzeta.add_double_mod( lbda_grv2 ) ; 
 					   
    // Parameters for the function Tenseur::poisson_vect
    // -------------------------------------------------

    Param par_poisson_vect ; 

    par_poisson_vect.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson_vect.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson_vect.add_double(precis_poisson, 1) ; // required precision
    par_poisson_vect.add_cmp_mod( ssjm1_khi ) ; 
    par_poisson_vect.add_tenseur_mod( ssjm1_wshift ) ; 
    par_poisson_vect.add_int_mod(niter, 0) ;   

 					   
    // Initializations
    // ---------------

    // Initial central angular velocity
    double omega_c = 0 ; 
  
    double accrois_omega = (omega_c0 - omega_ini) / 
			    double(mer_fix_omega - mer_change_omega) ; 


    update_metric() ;	// update of the metric coefficients

    equation_of_state() ;	// update of the density, pressure, etc...
    
    hydro_euler() ;	// update of the hydro quantities relative to the 
			//  Eulerian observer

    // Quantities at the previous step : 	
    Map_et mp_prev = mp_et ; 
    Tenseur ent_prev = ent ;	    
    Tenseur logn_prev = logn ;	    
    Tenseur dzeta_prev = dzeta ;	    
    
    // Creation of uninitialized tensors:
    Tenseur source_nuf(mp) ;    // source term in the equation for nuf
    Tenseur source_nuq(mp) ;    // source term in the equation for nuq
    Tenseur source_dzf(mp) ;	// matter source term in the eq. for dzeta
    Tenseur source_dzq(mp) ;	// quadratic source term in the eq. for dzeta
    Tenseur source_tggg(mp) ;	// source term in the eq. for tggg
    Tenseur source_shift(mp, 1, CON, mp.get_bvect_cart()) ;  
						    // source term for shift
    Tenseur mlngamma(mp) ;	// centrifugal potential
    
    // Preparations for the Poisson equations:
    // --------------------------------------
    if (nuf.get_etat() == ETATZERO) {
	nuf.set_etat_qcq() ; 
	nuf.set() = 0 ; 
    }
    
    if (relativistic) {
	if (nuq.get_etat() == ETATZERO) {
	    nuq.set_etat_qcq() ; 
	    nuq.set() = 0 ; 
	}

	if (tggg.get_etat() == ETATZERO) {
	    tggg.set_etat_qcq() ; 
	    tggg.set() = 0 ; 
	}
	
	if (dzeta.get_etat() == ETATZERO) {
	    dzeta.set_etat_qcq() ; 
	    dzeta.set() = 0 ; 
	}
    }
		    
    ofstream fichconv("convergence.d") ;    // Output file for diff_ent
    fichconv << "#     diff_ent     GRV2       max_triax      vit_triax" << endl ; 
    
    ofstream fichfreq("frequency.d") ;    // Output file for  omega_c
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
		omega_c = omega_ini ; 
	    }
	    else {
		if (mer <= mer_fix_omega) {
		    omega_c = omega_ini + accrois_omega * 
			      (mer - mer_change_omega) ;
		}
	    }

	}

	//-----------------------------------------------
	//  Sources of the Poisson equations
	//-----------------------------------------------
	
	// Source for nu
	// -------------
	Tenseur beta = log(bbb) ; 
	beta.set_std_base() ; 

	if (relativistic) {
	    source_nuf =  qpig * a_car *( ener_euler + s_euler ) ; 

	    source_nuq = ak_car - flat_scalar_prod(logn.gradient_spher(), 
			logn.gradient_spher() + beta.gradient_spher()) ; 
	}
	else {
	    source_nuf = qpig * nbar ; 

	    source_nuq = 0 ; 
	}
	source_nuf.set_std_base() ; 	
	source_nuq.set_std_base() ; 	

	// Source for dzeta
	// ----------------
	source_dzf = 2 * qpig * a_car * (press + (ener_euler+press) * uuu*uuu ) ;
	source_dzf.set_std_base() ; 
  
	source_dzq = 1.5 * ak_car - flat_scalar_prod(logn.gradient_spher(),
						     logn.gradient_spher() ) ;	    
	source_dzq.set_std_base() ; 	
	
	// Source for tggg
	// ---------------
	
	source_tggg = 4 * qpig * nnn * a_car * bbb * press ;
	source_tggg.set_std_base() ; 
	
	(source_tggg.set()).mult_rsint() ; 
	

	// Source for shift
	// ----------------
	    
	// Matter term: 
	source_shift = (-4*qpig) * nnn * a_car  * (ener_euler + press)
				* u_euler ;

	// Quadratic terms:
	Tenseur vtmp =  6 * beta.gradient_spher() - 2 * logn.gradient_spher() ;
	vtmp.change_triad(mp.get_bvect_cart()) ; 

	Tenseur squad  = nnn * flat_scalar_prod(tkij, vtmp) ;     

	// The addition of matter terms and quadratic terms is performed
	//  component by component because u_euler is contravariant,
	//  while squad is covariant. 
		
	if (squad.get_etat() == ETATQCQ) {
	    for (int i=0; i<3; i++) {
		source_shift.set(i) += squad(i) ; 
	    }
	}

	source_shift.set_std_base() ; 	

	//----------------------------------------------
	// Resolution of the Poisson equation for nuf 
	//----------------------------------------------

	source_nuf().poisson(par_poisson_nuf, nuf.set()) ; 
	
	cout << "Test of the Poisson equation for nuf :" << endl ; 
	Tbl err = source_nuf().test_poisson(nuf(), cout, true) ; 
	diff_nuf = err(0, 0) ; 
	
	//---------------------------------------
	// Triaxial perturbation of nuf
	//---------------------------------------
	
	if (mer == mer_triax) {
	
	    if ( mg->get_np(0) == 1 ) {
		cout << 
		"Et_rot_diff::equilibrium: np must be stricly greater than 1"
		<< endl << " to set a triaxial perturbation !" << endl ; 
		abort() ; 
	    }
	
	    const Coord& phi = mp.phi ; 
	    const Coord& sint = mp.sint ; 
	    Cmp perturb(mp) ; 
	    perturb = 1 + ampli_triax * sint*sint * cos(2*phi) ;
	    nuf.set() = nuf() * perturb ;  
	    
	    nuf.set_std_base() ;    // set the bases for spectral expansions
				    //  to be the standard ones for a 
				    //  scalar field

	}
	
	// Monitoring of the triaxial perturbation
	// ---------------------------------------
	
	Valeur& va_nuf = nuf.set().va ; 
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

	    source_nuq().poisson(par_poisson_nuq, nuq.set()) ; 
	    
	    cout << "Test of the Poisson equation for nuq :" << endl ; 
	    err = source_nuq().test_poisson(nuq(), cout, true) ;
	    diff_nuq = err(0, 0) ; 
	
	    //---------------------------------------------------------
	    // Resolution of the vector Poisson equation for the shift
	    //---------------------------------------------------------


	    if (source_shift.get_etat() != ETATZERO) {

		for (int i=0; i<3; i++) {
		    if(source_shift(i).dz_nonzero()) {
			assert( source_shift(i).get_dzpuis() == 4 ) ; 
		    }
		    else{
			(source_shift.set(i)).set_dzpuis(4) ; 
		    }
		}

	    }
	    //##
	    // source_shift.dec2_dzpuis() ;    // dzpuis 4 -> 2

	    double lambda_shift = double(1) / double(3) ; 

	    if ( mg->get_np(0) == 1 ) {
		lambda_shift = 0 ; 
	    }
	
	    source_shift.poisson_vect(lambda_shift, par_poisson_vect, 
				      shift, w_shift, khi_shift) ;      
	    
	    cout << "Test of the Poisson equation for shift_x :" << endl ; 
	    err = source_shift(0).test_poisson(shift(0), cout, true) ;
	    diff_shift_x = err(0, 0) ; 
	
	    cout << "Test of the Poisson equation for shift_y :" << endl ; 
	    err = source_shift(1).test_poisson(shift(1), cout, true) ;
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
		
    	    omega_c *= fact_omega ;  // Increase of the angular velocity if 
	}			   //  fact_omega != 1
	
	
	bool omega_trop_grand = false ; 
	bool kepler = true ; 

	while ( kepler ) {
	
	    // Possible decrease of Omega to ensure a velocity < c 
	
	    bool superlum = true ; 
	
	    while ( superlum ) {
	
		// Computation of Omega(r,theta)

		if (omega_c == 0.) {
		    omega_field = 0 ; 
		}
		else {
		    par_frot.set(0) = omega_c ; 
		    if (par_frot(2) != double(0)) {  // fixed a = R_eq / R_0
			par_frot.set(1) = ray_eq() / par_frot(2) ; 
		    }
		    double omeg_min = 0 ; 
		    double omeg_max = omega_c ; 
		    double precis1 = 1.e-14 ;
		    int nitermax1 = 100 ; 

		    fait_omega_field(omeg_min, omeg_max, precis1, nitermax1) ;
		}

		// New fluid velocity U :
	
		Cmp tmp = omega_field() - nphi() ; 
		tmp.annule(nzm1) ; 
		tmp.std_base_scal() ;
    
		tmp.mult_rsint() ;	    //  Multiplication by r sin(theta)

		uuu = bbb() / nnn() * tmp ; 
    
		if (uuu.get_etat() == ETATQCQ) {
		    // Same basis as (Omega -N^phi) r sin(theta) :
		    ((uuu.set()).va).set_base( (tmp.va).base ) ;   
		}


		// Is the new velocity larger than c in the equatorial plane ?
	
		superlum = false ; 
	
		for (int l=0; l<nzet; l++) {
		    for (int i=0; i<mg->get_nr(l); i++) {
		
			double u1 = uuu()(l, 0, j_b, i) ; 
			if (u1 >= 1.) {	    // superluminal velocity
			    superlum = true ; 
			    cout << "U > c  for l, i : " << l << "  " << i 
				 << "   U = " << u1 << endl ;  
			}
		    }
		}
		if ( superlum ) {
		    cout << "**** VELOCITY OF LIGHT REACHED ****" << endl ; 
		    omega_c /= fact_omega ;    // Decrease of Omega_c
		    cout << "New central rotation frequency : " 
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
	    double nuf_b  = nuf()(l_b, k_b, j_b, i_b) ; 
	    double nuq_b  = nuq()(l_b, k_b, j_b, i_b) ; 
	    double mlngamma_b  = mlngamma()(l_b, k_b, j_b, i_b) ; 
	    double primf_b = prim_field()(l_b, k_b, j_b, i_b) ;
	    
	    
	    // Central values of various potentials :
	    double nuf_c = nuf()(0,0,0,0) ; 
	    double nuq_c = nuq()(0,0,0,0) ; 
	    double mlngamma_c = 0 ;
	    double primf_c = prim_field()(0,0,0,0) ; 
	
	    // Scale factor to ensure that the enthalpy is equal to ent_b at 
	    //  the equator
	    double alpha_r2 = ( ent_c - ent_b + mlngamma_c - mlngamma_b
				+ nuq_c - nuq_b + primf_c - primf_b) 
				/ ( nuf_b - nuf_c  ) ;
	    alpha_r = sqrt(alpha_r2) ;
	    cout << "alpha_r = " << alpha_r << endl ; 

	    // Readjustment of nu :
	    // -------------------

	    logn = alpha_r2 * nuf + nuq ;
	    double nu_c =  logn()(0,0,0,0) ;

	    // First integral	--> enthalpy in all space
	    //-----------------

	    ent = (ent_c + nu_c + mlngamma_c) - logn - mlngamma - prim_field ;

	    // Test: is the enthalpy negative somewhere in the equatorial plane
	    //  inside the star ? If yes, this means that the Keplerian velocity
	    //  has been overstep.

	    kepler = false ; 
	    for (int l=0; l<nzet; l++) {
		int imax = mg->get_nr(l) - 1 ;
		if (l == l_b) imax-- ;	// The surface point is skipped
		for (int i=0; i<imax; i++) { 
		    if ( ent()(l, 0, j_b, i) < 0. ) {
			kepler = true ;
			cout << "ent < 0 for l, i : " << l << "  " << i 
			     << "   ent = " << ent()(l, 0, j_b, i) << endl ;  
		    } 
		}
	    }
	
	    if ( kepler ) {
		cout << "**** KEPLERIAN VELOCITY REACHED ****" << endl ; 
		omega_c /= fact_omega ;    // Omega is decreased
		cout << "New central rotation frequency : " 
			 << omega_c/(2.*M_PI) * f_unit << " Hz" << endl ; 
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
	
	double dent_eq = ent().dsdr()(l_b, k_b, j_b, i_b) ; 
	double dent_pole = ent().dsdr()(l_b, k_b, 0, i_b) ;
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

	mp.adapt(ent(), par_adapt) ; 

 	//----------------------------------------------------
	// Computation of the enthalpy at the new grid points
	//----------------------------------------------------
	
	mp_prev.homothetie(alpha_r) ; 
	
	mp.reevaluate(&mp_prev, nzet+1, ent.set()) ; 


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

	    mp.poisson2d(source_tggg(), mp.cmp_zero(), par_poisson_tggg,
			 tggg.set()) ; 
	    
	    //-------------------------------------------------------
	    //	2-D Poisson equation for dzeta
	    //-------------------------------------------------------

	    mp.poisson2d(source_dzf(), source_dzq(), par_poisson_dzeta,
			 dzeta.set()) ; 
	    
	    err_grv2 = lbda_grv2 - 1; 
	    cout << "GRV2: " << err_grv2 << endl ; 
	    
	}


	//---------------------------------------
	// Computation of the metric coefficients (except for N^phi)
	//---------------------------------------

	// Relaxations on nu and dzeta :  

	logn = relax * logn + relax_prev * logn_prev ;

	dzeta = relax * dzeta + relax_prev * dzeta_prev ; 

	// Update of the metric coefficients N, A, B and computation of K_ij :

	update_metric() ; 
	
	//-----------------------
	//  Informations display
	//-----------------------

	partial_display(cout) ; 
	fichfreq << "  " << omega_c / (2*M_PI) * f_unit ; 
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
		if (mer%4 == 0) omega_c *= fact ; 
	    }
	}
	

	//------------------------------------------------------------
	//  Relative change in enthalpy with respect to previous step 
	//------------------------------------------------------------

	Tbl diff_ent_tbl = diffrel( ent(), ent_prev() ) ; 
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

    fichconv.close() ; 
    fichfreq.close() ; 
    fichevol.close() ; 
    

}
