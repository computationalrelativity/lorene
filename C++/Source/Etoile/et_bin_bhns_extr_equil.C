/*
 *  Method of class Et_bin_bhns_extr to compute an equilibrium configuration
 *  of a BH-NS binary system with an extreme mass ratio
 *
 *    (see file et_bin_bhns_extr.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Keisuke Taniguchi
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

char et_bin_bhns_extr_equil_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/11/30 20:48:45  k_taniguchi
 * *** empty log message ***
 *
 * Revision 1.4  2003/10/19 20:01:10  e_gourgoulhon
 * Template file
 *
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "et_bin_bhns_extr.h"
#include "etoile.h"
#include "map.h"
#include "coord.h"
#include "param.h"
#include "eos.h"
#include "graphique.h"
#include "utilitaires.h"
#include "unites.h"

void Et_bin_bhns_extr::equil_bhns_extr(double ent_c, const double& mass,
				       const double& sepa, int mermax,
				       int mermax_poisson,
				       double relax_poisson, int mermax_potvit,
				       double relax_potvit, double thres_adapt,
				       Tbl& diff) {

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
    int nz_search = nzet ;	// Number of domains for searching
				//  the enthalpy isosurfaces
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
    par_adapt.add_double(reg_map, 1)	;  // 1 = regular mapping,
                                           // 0 = contracting mapping
    par_adapt.add_double(alpha_r, 2) ;	    // factor by which all the radial
					    // distances will be multiplied
    par_adapt.add_tbl(ent_limit, 0) ;	// array of values of the field ent
				        // to define the isosurfaces

    // Parameters for the function Map_et::poisson for logn_auto
    // ---------------------------------------------------------

    double precis_poisson = 1.e-16 ;

    Param par_poisson1 ;

    par_poisson1.add_int(mermax_poisson,  0) ; // maximum number of iterations
    par_poisson1.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson1.add_double(precis_poisson, 1) ; // required precision
    par_poisson1.add_int_mod(niter, 0) ;  // number of iterations actually
                                          // used 
    par_poisson1.add_cmp_mod( ssjm1_logn ) ;

    // Parameters for the function Map_et::poisson for beta_auto
    // ---------------------------------------------------------

    Param par_poisson2 ;

    par_poisson2.add_int(mermax_poisson,  0) ; // maximum number of iterations
    par_poisson2.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson2.add_double(precis_poisson, 1) ; // required precision
    par_poisson2.add_int_mod(niter, 0) ;  // number of iterations actually
                                          // used 
    par_poisson2.add_cmp_mod( ssjm1_beta ) ;

    // Parameters for the function Tenseur::poisson_vect
    // -------------------------------------------------

    Param par_poisson_vect ;

    par_poisson_vect.add_int(mermax_poisson, 0) ; // maximum number of
                                                  // iterations
    par_poisson_vect.add_double(relax_poisson, 0) ; // relaxation parameter
    par_poisson_vect.add_double(precis_poisson, 1) ; // required precision
    par_poisson_vect.add_cmp_mod( ssjm1_khi ) ;
    par_poisson_vect.add_tenseur_mod( ssjm1_wshift ) ;
    par_poisson_vect.add_int_mod(niter, 0) ;

    // External potential
    // See Eq (99) from Gourgoulhon et al. (2001)
    // -----------------------------------------

    Tenseur pot_ext = logn_comp + pot_centri + loggam ;

    Tenseur ent_jm1 = ent ;  // Enthalpy at previous step

    Tenseur source(mp) ;    // source term in the equation for logn_auto
			    // and beta_auto

    Tenseur source_shift(mp, 1, CON, ref_triad) ;  // source term in the
						   // equation for shift_auto

    //==========================================================//
    //                    Start of iteration                    //
    //==========================================================//

    for(int mer=0 ; mer<mermax ; mer++ ) {

        cout << "-----------------------------------------------" << endl ;
	cout << "step: " << mer << endl ;
	cout << "diff_ent = " << diff_ent << endl ;

	//------------------------------------------------------
	// Resolution of the elliptic equation for the velocity
	// scalar potential
	//------------------------------------------------------

	if (irrotational) {
	    diff_vel_pot = velocity_pot_extr(mass, sepa, mermax_potvit,
					     precis_poisson, relax_potvit) ;
	}

	//-------------------------------------
	// Computation of the new radial scale
	//-------------------------------------

	// alpha_r (r = alpha_r r') is determined so that the enthalpy
	// takes the requested value ent_b at the stellar surface

	// Values at the center of the star:
	double logn_auto_c  = logn_auto()(0, 0, 0, 0) ;
	double pot_ext_c  = pot_ext()(0, 0, 0, 0) ;

	// Search for the reference point (theta_*, phi_*)
	// [notation of Bonazzola, Gourgoulhon & Marck PRD 58, 104020 (1998)]
	// at the surface of the star
	// ------------------------------------------------------------------
	double alpha_r2 = 0 ;
	for (int k=0; k<mg->get_np(l_b); k++) {
	    for (int j=0; j<mg->get_nt(l_b); j++) {

	        double pot_ext_b  = pot_ext()(l_b, k, j, i_b) ;
		double logn_auto_b  = logn_auto()(l_b, k, j, i_b) ;

		// See Eq (100) from Gourgoulhon et al. (2001)
		double alpha_r2_jk = ( ent_c - ent_b + pot_ext_c - pot_ext_b)
		  / ( logn_auto_b - logn_auto_c ) ;

		if (alpha_r2_jk > alpha_r2) {
		    alpha_r2 = alpha_r2_jk ;
		    k_b = k ;
		    j_b = j ;
		}

	    }
	}

	alpha_r = sqrt(alpha_r2) ;

	cout << "k_b, j_b, alpha_r: " << k_b << "  " << j_b << "  " 
	     << alpha_r << endl ;

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

	//--------------------------------------------
	// First integral --> enthalpy in all space
	// See Eq (98) from Gourgoulhon et al. (2001)
	//--------------------------------------------

	ent = (ent_c + logn_auto_c + pot_ext_c) - logn_auto - pot_ext ;

	(ent().va).smooth(nzet, (ent.set()).va) ;

	//-----------------------------------------------------
	// Adaptation of the mapping to the new enthalpy field
	//----------------------------------------------------

	// Shall the adaptation be performed (cusp) ?
	// ------------------------------------------

	double dent_eq = ent().dsdr().val_point(ray_eq_pi(),M_PI/2.,M_PI) ;
	double dent_pole = ent().dsdr().val_point(ray_pole(),0.,0.) ;
	double rap_dent = fabs( dent_eq / dent_pole ) ;
	cout << "| dH/dr_eq / dH/dr_pole | = " << rap_dent << endl ;

	if ( rap_dent < thres_adapt ) {
	    adapt_flag = 0 ;	// No adaptation of the mapping
	    cout << "******* FROZEN MAPPING  *********" << endl ;
	}
	else{
	    adapt_flag = 1 ;	// The adaptation of the mapping is to be
	                        // performed
	}

	ent_limit.set_etat_qcq() ;
	for (int l=0; l<nzet; l++) {	// loop on domains inside the star
	    ent_limit.set(l) = ent()(l, k_b, j_b, i_b) ;
	}

	ent_limit.set(nzet-1) = ent_b ;
	Map_et mp_prev = mp_et ;

	mp.adapt(ent(), par_adapt) ;

	double fact_resize = 1. / alpha_r ;
	mp.resize(nz-2, fact_resize) ;
	mp_et.resize_extr(fact_resize) ;

	//----------------------------------------------------
	// Computation of the enthalpy at the new grid points
	//----------------------------------------------------

	mp_prev.homothetie(alpha_r) ;

	mp.reevaluate(&mp_prev, nzet+1, ent.set()) ;

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
	cout << "   " << ent_s_max << " reached for k = " << k_s_max
	     << " and j = " << j_s_max << endl ;

	//-------------------
	// Equation of state
	//-------------------

	equation_of_state() ; 	// computes new values for nbar (n), ener (e)
	                        // and press (p) from the new ent (H)

	//----------------------------------------------------------
	// Matter source terms in the gravitational field equations
	//---------------------------------------------------------

	hydro_euler_extr(mass, sepa) ;	// computes new values for
	                                // ener_euler (E), s_euler (S)
	                                // and u_euler (U^i)


	//----------------------------------------------------
	// Auxiliary terms for the source of Poisson equation
	//----------------------------------------------------

	const Coord& xx = mp.x ;
	const Coord& yy = mp.y ;
	const Coord& zz = mp.z ;

	Tenseur r_bh(mp) ;
	r_bh.set_etat_qcq() ;
	r_bh.set() = pow( (xx+sepa)*(xx+sepa) + yy*yy + zz*zz, 0.5) ;
	r_bh.set_std_base() ;

	Tenseur xx_cov(mp, 1, COV, ref_triad) ;
	xx_cov.set_etat_qcq() ;
	xx_cov.set(0) = xx + sepa ;
	xx_cov.set(1) = yy ;
	xx_cov.set(2) = zz ;
	xx_cov.set_std_base() ;

	Tenseur xsr_cov(mp, 1, COV, ref_triad) ;
	xsr_cov = xx_cov / r_bh ;
	xsr_cov.set_std_base() ;

	Tenseur msr(mp) ;
	msr = ggrav * mass / r_bh ;
	msr.set_std_base() ;

	Tenseur lapse_bh(mp) ;
	lapse_bh = 1. / sqrt( 1.+2.*msr ) ;
	lapse_bh.set_std_base() ;

	Tenseur lapse_bh2(mp) ;  // lapse_bh * lapse_bh
	lapse_bh2 = 1. / (1.+2.*msr) ;
	lapse_bh2.set_std_base() ;

	Tenseur ldnu(mp) ;
	ldnu = flat_scalar_prod_desal(xsr_cov, d_logn_auto) ;
	ldnu.set_std_base() ;

	Tenseur ldbeta(mp) ;
	ldbeta = flat_scalar_prod_desal(xsr_cov, d_beta_auto) ;
	ldbeta.set_std_base() ;

	Tenseur lltkij(mp) ;
	lltkij.set_etat_qcq() ;
	lltkij.set() = 0 ;
	lltkij.set_std_base() ;

	for (int i=0; i<3; i++)
	    for (int j=0; j<3; j++)
	        lltkij.set() += xsr_cov(i) * xsr_cov(j) * tkij_auto(i, j) ;

	Tenseur lshift(mp) ;
	lshift = flat_scalar_prod_desal(xsr_cov, shift_auto) ;
	lshift.set_std_base() ;

	Tenseur d_ldnu(mp, 1, COV, ref_triad) ;
	d_ldnu = ldnu.gradient() ;        // (d/dx, d/dy, d/dz)
	d_ldnu.change_triad(ref_triad) ;  // --> (d/dX, d/dY, d/dZ)

	Tenseur ldldnu(mp) ;
	ldldnu = flat_scalar_prod_desal(xsr_cov, d_ldnu) ;
	ldldnu.set_std_base() ;

	Tenseur d_ldbeta(mp, 1, COV, ref_triad) ;
	d_ldbeta = ldbeta.gradient() ;      // (d/dx, d/dy, d/dz)
	d_ldbeta.change_triad(ref_triad) ;  // --> (d/dX, d/dY, d/dZ)

	Tenseur ldldbeta(mp) ;
	ldldbeta = flat_scalar_prod_desal(xsr_cov, d_ldbeta) ;
	ldldbeta.set_std_base() ;

	//------------------------------------------
	// Poisson equation for logn_auto (nu_auto)
	//------------------------------------------

	// Source
	// ------

	if (relativistic) {
	    source = qpig * a_car % (ener_euler + s_euler) + akcar_auto
	      - flat_scalar_prod_desal(d_logn_auto, d_beta_auto)
	      + 2.*lapse_bh2 % msr % (ldnu % ldbeta + ldldnu)
	      + lapse_bh2 % lapse_bh2 % msr % (2.*(ldnu + 4.*msr % ldnu)
					       - ldbeta) / r_bh
	      - (4.*a_car % lapse_bh2 % lapse_bh2 % msr / 3. / nnn / r_bh)
	      * (2.+3.*msr) * (3.+4.*msr) * lltkij
	      + (2.*a_car % lapse_bh2 % lapse_bh2 % lapse_bh % msr
		 / nnn / r_bh / r_bh) * (2.+10.*msr+9.*msr%msr) * lshift
	      + (4.*pow(lapse_bh2, 3.) % msr % msr / 3. / r_bh / r_bh)
	      % (2.*(a_car%lapse_bh2/nnn/nnn - 1.) * pow(2.+3.*msr, 2.)
		 + (a_car - 1.) % pow(1.+3.*msr, 2.)
		 - 3.*(a_car%lapse_bh/nnn - 1.)*(2.+10.*msr+9.*msr%msr)) ;
	}
	else {
	    cout << "The computation of BH-NS binary systems"
		 << " should be relativistic !!!" << endl ;
	    abort() ;
	}

	source.set_std_base() ;

	// Resolution of the Poisson equation
	// ----------------------------------

	int k_falloff = 1 ;
	source().poisson_falloff(par_poisson1, logn_auto.set(), k_falloff) ;

	// Construct logn_auto_regu for et_bin_upmetr_extr.C
	// -------------------------------------------------

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

	if (relativistic) {

	    //--------------------------------
	    // Poisson equation for beta_auto
	    //--------------------------------

	    // Source
	    // ------

	    source = qpig * a_car % s_euler + 0.75 * akcar_auto
	      - 0.5 * flat_scalar_prod_desal(d_logn_auto, d_logn_auto)
	      - 0.5 * flat_scalar_prod_desal(d_beta_auto, d_beta_auto)
	      + lapse_bh2 % msr % (ldnu%ldnu + ldbeta%ldbeta + 2.*ldldbeta)
	      + lapse_bh2 % lapse_bh2 % msr % (2.*(1.+4.*msr) * ldbeta
					       - ldnu) / r_bh
	      - (a_car % lapse_bh2 % lapse_bh2 % msr / nnn / r_bh)
	      * (2.+3.*msr) * (3.+4.*msr) * lltkij
	      + (2.*a_car % lapse_bh2 % lapse_bh2 % lapse_bh % msr
		 / nnn / r_bh / r_bh) * (2.+10.*msr+9.*msr%msr) * lshift
	      + (2.*pow(lapse_bh2, 3.) % msr % msr / r_bh / r_bh)
	      % ((a_car%lapse_bh2/nnn/nnn - 1.) * pow(2.+3.*msr, 2.)
		 + (a_car - 1.) * pow(1.+3.*msr, 2.)
		 - 2.*(a_car%lapse_bh/nnn - 1.)*(2.+10.*msr+9.*msr%msr)) ;

	    source.set_std_base() ;

	    // Resolution of the Poisson equation
	    // ----------------------------------

	    source().poisson_falloff(par_poisson2, beta_auto.set(),
				     k_falloff) ;

	    // Check: has the Poisson equation been correctly solved ?
	    // -----------------------------------------------------

	    Tbl tdiff_beta = diffrel(beta_auto().laplacien(), source()) ;

	    cout << "Relative error in the resolution of the equation for "
		 << "beta_auto : " << endl ;
	    for (int l=0; l<nz; l++) {
	        cout << tdiff_beta(l) << "  " ;
	    }
	    cout << endl ;
	    diff_beta = max(abs(tdiff_beta)) ;

	    //----------------------------------------
	    // Vector Poisson equation for shift_auto
	    //----------------------------------------

	    // Some auxiliary terms for the source
	    // -----------------------------------

	    Tenseur xx_con(mp, 1, CON, ref_triad) ;
	    xx_con.set_etat_qcq() ;
	    xx_con.set(0) = xx + sepa ;
	    xx_con.set(1) = yy ;
	    xx_con.set(2) = zz ;
	    xx_con.set_std_base() ;

	    Tenseur xsr_con(mp, 1, CON, ref_triad) ;
	    xsr_con = xx_con / r_bh ;
	    xsr_con.set_std_base() ;

	    // Components of shift_auto with respect to the Cartesian triad
	    //  (d/dx, d/dy, d/dz) of the mapping :

	    Tenseur shift_auto_local = shift_auto ;
	    shift_auto_local.change_triad( mp.get_bvect_cart() ) ;

	    // Gradient (partial derivatives with respect to the Cartesian
	    //           coordinates of the mapping)
	    // dn(i, j) = D_i N^j

	    Tenseur dn = shift_auto_local.gradient() ;

	    // Return to the absolute reference frame
	    dn.change_triad(ref_triad) ;

	    // Trace of D_i N^j = divergence of N^j :
	    Tenseur divn = contract(dn, 0, 1) ;

	    // l^j D_j N^i
	    Tenseur ldn_con = contract(xsr_con, 0, dn, 0) ;

	    // D_j (l^k D_k N^i): dldn(j, i)
	    Tenseur ldn_local = ldn_con ;
	    ldn_local.change_triad( mp.get_bvect_cart() ) ;
	    Tenseur dldn = ldn_local.gradient() ;
	    dldn.change_triad(ref_triad) ;

	    // l^j D_j (l^k D_k N^i)
	    Tenseur ldldn = contract(xsr_con, 0, dldn, 0) ;

	    // l_k D_j N^k
	    Tenseur ldn_cov = contract(xsr_cov, 0, dn, 1) ;

	    // l^j l_k D_j N^k
	    Tenseur lldn_cov = contract(xsr_con, 0, ldn_cov, 0) ;

	    // eta^{ij} l_k D_j N^k
	    Tenseur eldn(mp, 1, CON, ref_triad) ;
	    eldn.set_etat_qcq() ;
	    eldn.set(0) = ldn_cov(0) ;
	    eldn.set(1) = ldn_cov(1) ;
	    eldn.set(2) = ldn_cov(2) ;
	    eldn.set_std_base() ;

	    // l^i D_j N^j
	    Tenseur ldivn = xsr_con % divn ;

	    // D_j (l^i D_k N^k): dldivn(j, i)
	    Tenseur ldivn_local = ldivn ;
	    ldivn_local.change_triad( mp.get_bvect_cart() ) ;
	    Tenseur dldivn = ldivn_local.gradient() ;
	    dldivn.change_triad(ref_triad) ;

	    // l^j D_j (l^i D_k N^k)
	    Tenseur ldldivn = contract(xsr_con, 0, dldivn, 0) ;

	    // l_j N^j
	    Tenseur ln = contract(xsr_cov, 0, shift_auto, 0) ;

	    Tenseur vtmp =  6. * d_beta_auto - 8. * d_logn_auto ;

	    Tenseur lvtmp = contract(xsr_con, 0, vtmp, 0) ;

	    // eta^{ij} vtmp_j
	    Tenseur evtmp(mp, 1, CON, ref_triad) ;
	    evtmp.set_etat_qcq() ;
	    evtmp.set(0) = vtmp(0) ;
	    evtmp.set(1) = vtmp(1) ;
	    evtmp.set(2) = vtmp(2) ;
	    evtmp.set_std_base() ;

	    // lapse_ns
	    Tenseur lapse_ns(mp) ;
	    lapse_ns = exp(logn_auto) ;
	    lapse_ns.set_std_base() ;

	    // Source
	    // ------

	    source_shift = (-4.*qpig) * nnn % a_car % (ener_euler + press)
	      % u_euler
	      + nnn % flat_scalar_prod_desal(tkij_auto, vtmp)
	      - 2.*nnn % lapse_bh2 * msr / r_bh
	      % flat_scalar_prod_desal(tkij_auto, xsr_cov)
	      + 2.*lapse_bh2 * msr * (3.*ldldn + ldldivn) / 3.
	      - lapse_bh2 * msr / r_bh
	      * (4.*ldivn - lapse_bh2 % (3.*ldn_con + 8.*msr * ldn_con)
		 - (eldn + 2.*lapse_bh2*(9.+11.*msr)*lldn_cov%xsr_con) / 3.)
	      - 2.*lapse_bh2 % lapse_bh2 * msr / r_bh / r_bh
	      * ( (4.+11.*msr) * shift_auto
		  - lapse_bh2 * (12.+51.*msr+46.*msr*msr) * ln % xsr_con )
	      / 3.
	      + 8.*pow(lapse_bh2, 4.) * msr / r_bh / r_bh
	      % (lapse_ns - 1.) * (2.+10.*msr+9.*msr*msr) * xsr_con / 3.
	      + 2.*pow(lapse_bh2, 3.) * msr / r_bh * (2.+3.*msr)
	      * ( (1.+2.*msr) * evtmp - (3.+2.*msr) * lvtmp * xsr_con) / 3. ;

	    source_shift.set_std_base() ;

	    // Resolution of the Poisson equation
	    // ----------------------------------

	    // Filter for the source of shift vector : No filter (tempo)
	    /*
	    for (int i=0; i<3; i++) {

	        if (source_shift(i).get_etat() != ETATZERO)
		    source_shift.set(i).filtre(4) ;

	    }
	    */

	    // For Tenseur::poisson_vect, the triad must be the mapping
	    // triad, not the reference one:

	    source_shift.change_triad( mp.get_bvect_cart() ) ;

	    for (int i=0; i<3; i++) {
	        if(source_shift(i).dz_nonzero()) {
		    assert( source_shift(i).get_dzpuis() == 4 ) ;
		}
		else {
		    (source_shift.set(i)).set_dzpuis(4) ;
		}
	    }

	    source_shift.dec2_dzpuis() ;    // dzpuis 4 -> 2

	    double lambda_shift = double(1) / double(3) ;

	    source_shift.poisson_vect_falloff(lambda_shift, par_poisson_vect,
					      shift_auto, w_shift,
					      khi_shift, k_falloff) ;

	    // Check: has the equation for shift_auto been correctly solved ?
	    // --------------------------------------------------------------

	    // Divergence of shift_auto :
	    Tenseur divna = contract(shift_auto.gradient(), 0, 1) ;
	    divna.dec2_dzpuis() ;    // dzpuis 2 -> 0

	    // Grad(div) :
	    Tenseur graddivn = divna.gradient() ;
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

	    cout << "Relative error in the resolution of the equation for "
		 << "shift_auto : " << endl ;
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
	    // The output of Tenseur::poisson_vect_falloff is on the mapping
	    // triad, it should therefore be transformed to components on the
	    // reference triad :

	    shift_auto.change_triad( ref_triad ) ;

	}   // End of relativistic equations

	//------------------------------
	//  Relative change in enthalpy
	//------------------------------

	Tbl diff_ent_tbl = diffrel( ent(), ent_jm1() ) ;
	diff_ent = diff_ent_tbl(0) ;
	for (int l=1; l<nzet; l++) {
	    diff_ent += diff_ent_tbl(l) ;
	}
	diff_ent /= nzet ;

	ent_jm1 = ent ;

    } // End of main loop

    //========================================================//
    //                    End of iteration                    //
    //========================================================//

}

