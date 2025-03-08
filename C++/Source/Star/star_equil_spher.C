/*
 * Method of class Star to compute a static spherical configuration.
 *
 * (see file star.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Francois Limousin
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
 * Revision 1.17  2017/10/20 13:54:54  j_novak
 * Adaptation of the method to be used within nrotstar.
 *
 * Revision 1.16  2016/12/05 16:18:15  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.15  2014/10/13 08:53:39  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.14  2010/10/18 20:16:10  m_bejger
 * Commented-out the reevaluation of the mapping for the case of many domains inside the star
 *
 * Revision 1.13  2010/10/18 18:56:33  m_bejger
 * Changed to allow initial data with more than one domain in the star
 *
 * Revision 1.12  2005/09/14 12:30:52  f_limousin
 * Saving of fields lnq and logn in class Star.
 *
 * Revision 1.11  2005/09/13 19:38:31  f_limousin
 * Reintroduction of the resolution of the equations in cartesian coordinates.
 *
 * Revision 1.10  2005/02/17 17:32:35  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.9  2004/06/22 12:45:31  f_limousin
 * Improve the convergence
 *
 * Revision 1.8  2004/06/07 16:27:47  f_limousin
 * Add the computation of virial error.
 *
 * Revision 1.6  2004/04/08 16:33:42  f_limousin
 * The new variable is ln(Q) instead of Q=psi^2*N. It improves the
 * convergence of the code.
 *
 * Revision 1.5  2004/03/25 10:29:26  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.4  2004/02/27 09:57:55  f_limousin
 * We now update the metric gamma at the end of this routine for
 * the calculus of mass_b and mass_g.
 *
 * Revision 1.3  2004/02/21 17:05:13  e_gourgoulhon
 * Method Scalar::point renamed Scalar::val_grid_point.
 * Method Scalar::set_point renamed Scalar::set_grid_point.
 *
 * Revision 1.2  2004/01/20 15:20:35  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "tenseur.h"
#include "star_rot.h"
#include "param.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"	    

namespace Lorene {
void Star::equilibrium_spher(double ent_c, double precis, const Tbl* pent_limit){
    
    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;
    
    // Initializations
    // ---------------
    
    const Mg3d* mg = mp.get_mg() ; 
    int nz = mg->get_nzone() ;	    // total number of domains
    
    // Index of the point at phi=0, theta=pi/2 at the surface of the star:
    int l_b = nzet - 1 ; 
    int i_b = mg->get_nr(l_b) - 1 ; 
    int j_b = mg->get_nt(l_b) - 1 ; 
    int k_b = 0 ; 
    
    // Value of the enthalpy defining the surface of the star
    double ent_b = 0 ; 
    
    // Initialization of the enthalpy field to the constant value ent_c :
    
    ent = ent_c ; 
    ent.annule(nzet, nz-1) ; 
    
    
    // Corresponding profiles of baryon density, energy density and pressure
    
    equation_of_state() ; 

    Scalar a_car(mp) ;
    a_car = 1 ;	    
    lnq = 0 ;
    lnq.std_spectral_base() ;

    // Auxiliary quantities
    // --------------------
    
    // Affine mapping for solving the Poisson equations
    Map_af mpaff(mp);	
        
    Param par_nul   ;	 // Param (null) for Map_af::poisson.

    Scalar ent_jm1(mp) ;	// Enthalpy at previous step
    ent_jm1 = 0 ; 
    
    Scalar source(mp) ; 
    Scalar logn_quad(mp) ; 
    Scalar logn_mat(mp) ;
    logn_mat = 0 ;
    logn = 0 ; 
    logn_quad = 0 ; 

    Scalar dlogn(mp) ; 
    Scalar dlnq(mp) ; 
    
    double diff_ent = 1 ;     
    int mermax = 200 ;	    // Max number of iterations
    
    //=========================================================================
    // 			Start of iteration
    //=========================================================================

    for(int mer=0 ; (diff_ent > precis) && (mer<mermax) ; mer++ ) {

      double alpha_r = 1 ; 
	
	cout << "-----------------------------------------------" << endl ;
	cout << "step: " << mer << endl ;
	cout << "alpha_r: " << alpha_r << endl ;
	cout << "diff_ent = " << diff_ent << endl ;    

	//-----------------------------------------------------
	// Resolution of Poisson equation for ln(N)
	//-----------------------------------------------------

	// Matter part of ln(N)
	// --------------------
	
	source = a_car * (ener + 3.*press) ;
    	
	source.set_dzpuis(4) ; 
	
	Cmp source_logn_mat (source) ;
	Cmp logn_mat_cmp (logn_mat) ;
	logn_mat_cmp.set_etat_qcq() ;
	
	mpaff.poisson(source_logn_mat, par_nul, logn_mat_cmp) ; 

	logn_mat = logn_mat_cmp ;

	// NB: at this stage logn is in machine units, not in c^2

	// Quadratic part of ln(N)
	// -----------------------

	mpaff.dsdr(logn, dlogn) ;
	mpaff.dsdr(lnq, dlnq) ;
		
	source = - dlogn * dlnq ; 

	Cmp source_logn_quad (source) ;
	Cmp logn_quad_cmp (logn_quad) ;
	logn_quad_cmp.set_etat_qcq() ;
	    
	mpaff.poisson(source_logn_quad, par_nul, logn_quad_cmp) ; 

	logn_quad = logn_quad_cmp ;

	    	    
	//-----------------------------------------------------
	// Computation of the new radial scale
	//-----------------------------------------------------

	// alpha_r (r = alpha_r r') is determined so that the enthalpy
	// takes the requested value ent_b at the stellar surface
	
	double nu_mat0_b  = logn_mat.val_grid_point(l_b, k_b, j_b, i_b) ; 
	double nu_mat0_c  = logn_mat.val_grid_point(0, 0, 0, 0) ; 

	double nu_quad0_b  = logn_quad.val_grid_point(l_b, k_b, j_b, i_b) ; 
	double nu_quad0_c  = logn_quad.val_grid_point(0, 0, 0, 0) ; 

	double alpha_r2 = ( ent_c - ent_b - nu_quad0_b + nu_quad0_c )
				/ ( qpig*(nu_mat0_b - nu_mat0_c) ) ;

	alpha_r = sqrt(alpha_r2) ;
	
    // One domain inside the star:
    // --------------------------- 
	if(nzet==1)  mpaff.homothetie( alpha_r ) ; 

	//--------------------
	// First integral
	//--------------------

	// Gravitation potential in units c^2 :
	logn_mat = alpha_r2*qpig * logn_mat ;
	logn = logn_mat + logn_quad ;

	// Enthalpy in all space
	double logn_c = logn.val_grid_point(0, 0, 0, 0) ;
	ent = ent_c - logn + logn_c ;

    // Two or more domains inside the star:
    // ------------------------------------ 

	if(nzet>1) { 

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

    int nzadapt = nzet ; 

    //cout << "no. of domains where the ent adjustment will be done: " << nzet << endl ;
    //cout << "ent limits: " << ent_limit << endl ; 

    double precis_adapt = 1.e-14 ; 
 
    double reg_map = 1. ; // 1 = regular mapping, 0 = contracting mapping

    par_adapt.add_int(nitermax, 0) ;   // maximum number of iterations to
				       // locate zeros by the secant method
    par_adapt.add_int(nzadapt, 1) ;    // number of domains where the adjustment 
				       // to the isosurfaces of ent is to be 
				       // performed
    par_adapt.add_int(nz_search, 2) ;  // number of domains to search for
				       // the enthalpy isosurface
    par_adapt.add_int(adapt_flag, 3) ; // 1 = performs the full computation, 
				       // 0 = performs only the rescaling by 
				       // the factor alpha_r
    par_adapt.add_int(j_b, 4) ;        // theta index of the collocation point 
			               // (theta_*, phi_*)
    par_adapt.add_int(k_b, 5) ;        // theta index of the collocation point 
			               // (theta_*, phi_*)

    par_adapt.add_int_mod(niter, 0) ;  // number of iterations actually used in 
				       // the secant method
    
    par_adapt.add_double(precis_adapt, 0) ; // required absolute precision in 
					// the determination of zeros by 
					// the secant method
    par_adapt.add_double(reg_map, 1) ;  // 1. = regular mapping, 0 = contracting mapping
    
    par_adapt.add_double(alpha_r, 2) ;	// factor by which all the radial 
					// distances will be multiplied 
    	   
    // Enthalpy values for the adaptation
    Tbl ent_limit(nzet) ; 
    if (pent_limit != 0x0) ent_limit = *pent_limit ; 
    	   
    par_adapt.add_tbl(ent_limit, 0) ;	// array of values of the field ent 
				        // to define the isosurfaces. 			   

      double* bornes = new double[nz+1] ;
      bornes[0] = 0. ; 

      for(int l=0; l<nz; l++) {
 
	bornes[l+1] = mpaff.get_alpha()[l]  + mpaff.get_beta()[l] ;     

      } 
      bornes[nz] = __infinity ; 
 
      Map_et mp0(*mg, bornes) ;
 
        mp0 = mpaff; 
        mp0.adapt(ent, par_adapt) ; 

        //Map_af mpaff_prev (mpaff) ; 

        double alphal, betal ; 

        for(int l=0; l<nz; l++) { 

			alphal = mp0.get_alpha()[l] ; 
			betal  = mp0.get_beta()[l] ;
 
			mpaff.set_alpha(alphal, l) ; 
			mpaff.set_beta(betal, l) ;

        } 
 

		// testing printout    
		//int num_r1 = mg->get_nr(0) - 1;        
		//cout << "Pressure difference:" << press.val_grid_point(0,0,0,num_r1) 
		//- press.val_grid_point(1,0,0,0) << endl ; 
		//cout << "Difference in enthalpies at the domain boundary:" << endl ;
		//cout << ent.val_grid_point(0,0,0,num_r1) << endl ; 
		//cout << ent.val_grid_point(1,0,0,0) << endl ;

		//cout << "Enthalpy difference: " << ent.val_grid_point(0,0,0,num_r1)
		//- ent.val_grid_point(1,0,0,0) << endl ;

		// Computation of the enthalpy at the new grid points
		//----------------------------------------------------                     
		//mpaff.reevaluate(&mpaff_prev, nzet+1, ent) ;

    } 


	//---------------------
	// Equation of state
	//---------------------
	
	equation_of_state() ; 
	
	//---------------------
	// Equation for lnq_auto
	//---------------------
	
	mpaff.dsdr(logn, dlogn) ;
	mpaff.dsdr(lnq, dlnq) ;
	
	source = 3 * qpig * a_car * press ;
	
	source = source - 0.5 * ( dlnq * dlnq + dlogn * dlogn ) ;

	source.std_spectral_base() ;
	Cmp source_lnq (source) ;
	Cmp lnq_cmp (logn_quad) ;
	lnq_cmp.set_etat_qcq() ;

	mpaff.poisson(source_lnq, par_nul, lnq_cmp) ; 

	lnq = lnq_cmp ;

	// Metric coefficient psi4 update
	
	nn = exp( logn ) ;
	nn.std_spectral_base() ;

	Scalar qq = exp( lnq ) ;
	qq.std_spectral_base() ;

	a_car = qq * qq / ( nn * nn ) ;
	a_car.std_spectral_base() ;

    // Relative difference with enthalpy at the previous step
    // ------------------------------------------------------
    
	diff_ent = norme( diffrel(ent, ent_jm1) ) / nzet ; 

    // Next step
    // ---------
    
    ent_jm1 = ent ; 


    }  // End of iteration loop 
    
    //=========================================================================
    // 			End of iteration
    //=========================================================================
 
    // The mapping is transfered to that of the star:
    // ----------------------------------------------
    mp = mpaff ; 

    // Sets value to all the Tenseur's of the star
    // -------------------------------------------
    
    nn = exp( logn ) ;
    nn.std_spectral_base() ;
   
    Scalar qq = exp( lnq ) ;
    qq.std_spectral_base() ;
    
    a_car = qq * qq / ( nn * nn ) ;
   
    Sym_tensor gamma_cov(mp, COV, mp.get_bvect_cart()) ;
    gamma_cov.set_etat_zero() ;

    for (int i=1; i<=3; i++){
	gamma_cov.set(i,i) = a_car ;
    }
    gamma = gamma_cov ;

    // ... hydro 
    ent.annule(nzet, nz-1) ;	// enthalpy set to zero at the exterior of 
				// the star
    ener_euler = ener ; 
    s_euler = 3 * press ;
    gam_euler = 1 ;  
    for(int i=1; i<=3; i++) u_euler.set(i) = 0 ;

    Star_rot* p_star_rot = dynamic_cast<Star_rot*>(this) ;
    if ( p_star_rot != 0x0) {
      p_star_rot->a_car = a_car ;
      p_star_rot->bbb = qq / nn ;
      p_star_rot->b_car = p_star_rot->bbb * p_star_rot->bbb ;
      p_star_rot->dzeta = lnq ;
    }
    
    // Info printing
    // -------------
    
    cout << endl 
     << "Characteristics of the star obtained by Etoile::equilibrium_spher : " 
     << endl 
     << "-----------------------------------------------------------------" 
     << endl ; 

    double ray = mp.val_r(l_b, 1., M_PI/2., 0) ; 
    cout << "Coordinate radius  :       " << ray / km << " km" << endl ; 

    double rcirc = ray * sqrt(a_car.val_grid_point(l_b, k_b, j_b, i_b) ) ; 
	
    double compact = qpig/(4.*M_PI) * mass_g() / rcirc ; 

    cout << "Circumferential radius R : " << rcirc/km  << " km" << endl ;
    cout << "Baryon mass :	     " << mass_b()/msol << " Mo" << endl ;
    cout << "Gravitational mass M :  " << mass_g()/msol << " Mo" << endl ;
    cout << "Compacity parameter GM/(c^2 R) : " << compact << endl ;
     
    //-----------------
    // Virial theorem
    //-----------------
    
    //... Pressure term

    source = qpig * a_car * sqrt(a_car) * s_euler ; 
    source.std_spectral_base() ;	    
    double vir_mat = source.integrale() ; 
    
    //... Gravitational term

    Scalar tmp = lnq - logn ; 

    source =  - ( logn.dsdr() * logn.dsdr() 
		      - 0.5 * tmp.dsdr() * tmp.dsdr() ) 
		    * sqrt(a_car)  ; 

    source.std_spectral_base() ;	    
    double vir_grav = source.integrale() ; 

    //... Relative error on the virial identity GRV3
    
    double grv3 = ( vir_mat + vir_grav ) / vir_mat ;

    cout << "Virial theorem GRV3 : " << endl ; 
    cout << "     3P term    : " << vir_mat << endl ; 
    cout << "     grav. term : " << vir_grav << endl ; 
    cout << "     relative error : " << grv3 << endl ; 
    
    // To be compatible with class Etoile :
    //logn = logn - logn_quad ;

 
}
}
