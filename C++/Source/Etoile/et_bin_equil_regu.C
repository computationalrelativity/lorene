/*
 *  Method of class Etoile_bin to compute an equilibrium configuration
 *  by regularizing source.
 *
 *  (see file etoile.h for documentation).
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Keisuke Taniguchi
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


char et_bin_equil_regu_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/12/11 12:51:26  k_taniguchi
 * Change the multiplication "*" to "%"
 *   and flat_scalar_prod to flat_scalar_prod_desal.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.17  2001/08/07  09:49:00  keisuke
 * Change of the method to set the longest radius of a star
 *  on the first domain.
 * Addition of a new argument in Etoile_bin::equil_regular : Tbl fact.
 *
 * Revision 2.16  2001/06/22  08:54:53  keisuke
 * Set the inner values of the second domain of ent
 *   by using the outer ones of the first domain.
 *
 * Revision 2.15  2001/05/17  12:22:26  keisuke
 * Change of the method to calculate chi from setting position in map
 *  to val_point.
 *
 * Revision 2.14  2001/02/07  09:47:28  eric
 * unsgam1 est desormais donne par Eos::der_nbar_ent (cas newtonien)
 *   ou Eos::der_ener_ent (cas relativiste).
 *
 * Revision 2.13  2001/01/16  17:02:32  keisuke
 * *** empty log message ***
 *
 * Revision 2.12  2001/01/16  16:58:08  keisuke
 * Change the method to set the values on the surface.
 *
 * Revision 2.11  2001/01/10  16:45:34  keisuke
 * Set the inner values of the second domain of logn_auto
 *   by using the outer ones of the first domain.
 *
 * Revision 2.10  2000/12/20  10:33:14  eric
 * Changement important : nz_search = nzet ---> nz_search = nzet + 1
 *
 * Revision 2.9  2000/10/25  14:01:03  keisuke
 * Modif de Map_et::adapt: on y rentre desormais avec nz_search
 *  (dans le cas present nz_search = nzet).
 *
 * Revision 2.8  2000/10/06  15:29:01  keisuke
 * Change poisson_vect into poisson_vect_regu.
 *
 * Revision 2.7  2000/09/25  15:01:10  keisuke
 * Suppress "int" from the declaration of k_div.
 *
 * Revision 2.6  2000/09/22  15:51:39  keisuke
 * d_logn_auto est desormais calcule en dehors (dans update_metric).
 *
 * Revision 2.5  2000/09/13  09:50:33  keisuke
 * Minor change on change_triad.
 *
 * Revision 2.4  2000/09/08  15:57:31  keisuke
 * Change the basis of d_logn_auto_div from the spherical coordinate
 *  to the Cartesian one with respect to ref_triad.
 *
 * Revision 2.3  2000/09/07  15:47:19  keisuke
 * Minor change.
 *
 * Revision 2.2  2000/09/07  15:43:41  keisuke
 * Add a new argument in poisson_regular and suppress logn_auto_total.
 *
 * Revision 2.1  2000/08/29  14:01:43  keisuke
 * Modify the arguments of poisson_regular.
 *
 * Revision 2.0  2000/08/29  11:39:02  eric
 * Version provisoire.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <math.h>

// Headers Lorene
#include "etoile.h"
#include "param.h"
#include "eos.h"
//#include "graphique.h"
#include "utilitaires.h"

//********************************************************************

void Etoile_bin::equil_regular(double ent_c, int mermax, int mermax_poisson, 
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

    k_div = 2 ;  // Regularity parameter for poisson_regular

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
    double& diff_shift_x = diff.set(4) ;
    double& diff_shift_y = diff.set(5) ;
    double& diff_shift_z = diff.set(6) ;

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
    par_adapt.add_int(nz_search, 2) ;	// number of domains to search for
					// the enthalpy isosurface
    par_adapt.add_int(adapt_flag, 3) ; //  1 = performs the full computation, 
				       //  0 = performs only the rescaling by 
				       //      the factor alpha_r
    par_adapt.add_int(j_b, 4) ; //  theta index of the collocation point 
			        //  (theta_*, phi_*)
    par_adapt.add_int(k_b, 5) ; //  theta index of the collocation point 
			        //  (theta_*, phi_*)
    par_adapt.add_int_mod(niter, 0) ;  //  number of iterations actually
                                       //  used in the secant method
    par_adapt.add_double(precis_secant, 0) ; // required absolute precision in
					     // the determination of zeros by
					     // the secant method
    par_adapt.add_double(reg_map, 1)	;  // 1. = regular mapping,
                                           // 0 = contracting mapping
    par_adapt.add_double(alpha_r, 2) ;	    // factor by which all the radial 
					    // distances will be multiplied 
    par_adapt.add_tbl(ent_limit, 0) ;	// array of values of the field ent 
				        // to define the isosurfaces.

    // Parameters for the function Map_et::poisson_regular for logn_auto
    // -----------------------------------------------------------------

    double precis_poisson = 1.e-16 ;     

    Param par_poisson1 ; 

    par_poisson1.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson1.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson1.add_double(precis_poisson, 1) ; // required precision
    par_poisson1.add_int_mod(niter, 0) ;  //  number of iterations actually
                                          //  used 
    par_poisson1.add_cmp_mod( ssjm1_logn ) ;

    // Parameters for the function Map_et::poisson for beta_auto
    // ---------------------------------------------------------

    Param par_poisson2 ;

    par_poisson2.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson2.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson2.add_double(precis_poisson, 1) ; // required precision
    par_poisson2.add_int_mod(niter, 0) ;  //  number of iterations actually
                                          //  used 
    par_poisson2.add_cmp_mod( ssjm1_beta ) ;

    // Parameters for the function Tenseur::poisson_vect_regu
    // ------------------------------------------------------

    Param par_poisson_vect ;

    par_poisson_vect.add_int(mermax_poisson,  0) ;  // maximum number of
                                                    // iterations
    par_poisson_vect.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson_vect.add_double(precis_poisson, 1) ; // required precision
    par_poisson_vect.add_cmp_mod( ssjm1_khi ) ;
    par_poisson_vect.add_tenseur_mod( ssjm1_wshift ) ;
    par_poisson_vect.add_int_mod(niter, 0) ;

	   
    // External potential
    // ------------------

    Tenseur pot_ext = logn_comp + pot_centri + loggam ;
//##
//	des_coupe_z(pot_ext(), 0., 1, "pot_ext", &(ent()) ) ;
//##

    Tenseur ent_jm1 = ent ;	// Enthalpy at previous step

    Tenseur source(mp) ;    // source term in the equation for logn_auto
			    // and beta_auto

    Tenseur source_shift(mp, 1, CON, ref_triad) ;  // source term in the
                                                   // equation for shift_auto

    Cmp source_regu(mp) ;
    Cmp source_div(mp) ;

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
	//-----------------------------------------------------

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


	//--------------------
	// First integral	--> enthalpy in all space
	//--------------------

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

//##
//	des_coupe_z(ent(), 0., 1, "ent before adapt", &(ent()) ) ;
//##

	mp.adapt(ent(), par_adapt) ;

	// Readjustment of the external boundary of domain l=nzet
	// to keep a fixed ratio with respect to star's surface
	
	double rr_in = mp.val_r(nzet,-1., M_PI/2, 0.) ; 
	double rr_out = mp.val_r(nzet,1., M_PI/2, 0.) ; 

	mp.resize(nzet, rr_in/rr_out * fact_resize(0)) ; 

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
	// ------

	double unsgam1 ;    // effective power of H in the source 
			    // close to the surface

	if (relativistic) {
	    source = qpig * a_car % (ener_euler + s_euler)
		    + akcar_auto + akcar_comp
		    - flat_scalar_prod_desal(d_logn_auto,
				       d_beta_auto + d_beta_comp) ;
				       
	    // 1/(gam-1) = dln(e)/dln(H) close to the surface : 
	    unsgam1 = eos.der_ener_ent_p(ent_b + 1e-10*(ent_c-ent_b)) ; 
	
	}
	else {
	    source = qpig * nbar ;
	    
	    // 1/(gam-1) = dln(n)/dln(H) close to the surface : 
	    unsgam1 = eos.der_nbar_ent_p(ent_b + 1e-10*(ent_c-ent_b)) ; 
	}

	source.set_std_base() ;

	// Resolution of the Poisson equation
	// ----------------------------------

	logn_auto_regu.set_etat_qcq() ;
	logn_auto_div.set_etat_qcq() ;
	d_logn_auto_div.set_etat_qcq() ;

	source_regu.std_base_scal() ;
	source_div.std_base_scal() ;

	source().poisson_regular(k_div, nzet, unsgam1, par_poisson1,
				 logn_auto.set(), logn_auto_regu.set(),
				 logn_auto_div.set(),
				 d_logn_auto_div,
				 source_regu, source_div) ;

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


	if (relativistic) {

	    //--------------------------------------------------------
	    // Poisson equation for beta_auto
	    //--------------------------------------------------------

	    // Source
	    // ------

	    source = qpig * a_car % s_euler
		    + .75 * ( akcar_auto + akcar_comp )
		    - .5 * flat_scalar_prod_desal(d_logn_auto,
					    d_logn_auto + d_logn_comp)
		    - .5 * flat_scalar_prod_desal(d_beta_auto,
					    d_beta_auto + d_beta_comp) ;

	    source.set_std_base() ;

	    // Resolution of the Poisson equation
	    // ----------------------------------

	    source().poisson(par_poisson2, beta_auto.set()) ;


	    // Check: has the Poisson equation been correctly solved ?
	    // -----------------------------------------------------

	    Tbl tdiff_beta = diffrel(beta_auto().laplacien(), source()) ;
	    cout <<
	    "Relative error in the resolution of the equation for beta_auto : "
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
	    // ------

	    Tenseur vtmp =  6. * ( d_beta_auto + d_beta_comp )
			   -8. * ( d_logn_auto + d_logn_comp ) ;

	    source_shift = (-4.*qpig) * nnn % a_car % (ener_euler + press)
	                        % u_euler
			   + nnn % flat_scalar_prod_desal(tkij_auto, vtmp) ;

	    source_shift.set_std_base() ;

	    // Resolution of the Poisson equation
	    // ----------------------------------

	    // Filter for the source of shift vector

	    for (int i=0; i<3; i++) {

	      if (source_shift(i).get_etat() != ETATZERO)
		source_shift.set(i).filtre(4) ;

	    }

	    // For Tenseur::poisson_vect_regu,
	    // the triad must be the mapping triad,
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
	    // source_shift.dec2_dzpuis() ;    // dzpuis 4 -> 2

	    double lambda_shift = double(1) / double(3) ;

	    source_shift.poisson_vect_regu(k_div, nzet, unsgam1,
					   lambda_shift, par_poisson_vect,
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
	      "Relative error in the resolution of the equation "
	      "for shift_auto : "
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


	}   // End of relativistic equations


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
