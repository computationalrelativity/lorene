/*
 * Method of class Etoile to compute a static spherical configuration
 *  with the outer boundary condition at a finite radius
 *
 * (see file etoile.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Keisuke Taniguchi
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


char etoile_eqsph_falloff_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/13 08:52:58  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2004/11/30 20:52:24  k_taniguchi
 * *** empty log message ***
 *
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "etoile.h"
#include "param.h"
#include "unites.h"	    

namespace Lorene {
void Etoile::equil_spher_falloff(double ent_c, double precis) {
    
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
    
    // Initial metric 
    a_car = 1 ;	     // this value will remain unchanged in the Newtonian case
    beta_auto = 0 ;  // this value will remain unchanged in the Newtonian case
    

    // Auxiliary quantities
    // --------------------
    
    // Affine mapping for solving the Poisson equations
    Map_af mpaff(mp);	
        
    Param par_nul   ;	 // Param (null) for Map_af::poisson.

    Tenseur ent_jm1(mp) ;	// Enthalpy at previous step
    ent_jm1 = 0 ; 
    
    Tenseur source(mp) ; 
    Tenseur logn(mp) ; 
    Tenseur logn_quad(mp) ; 
    logn = 0 ; 
    logn_quad = 0 ; 

    Cmp dlogn(mp) ; 
    Cmp dbeta(mp) ; 
    
    double diff_ent = 1 ;     
    int mermax = 200 ;	    // Max number of iterations
    
    double alpha_r = 1 ; 
    int k_falloff = 1 ;

    //=========================================================================
    // 			Start of iteration
    //=========================================================================

    for(int mer=0 ; (diff_ent > precis) && (mer<mermax) ; mer++ ) {

	cout << "-----------------------------------------------" << endl ;
	cout << "step: " << mer << endl ;
	cout << "alpha_r: " << alpha_r << endl ;
	cout << "diff_ent = " << diff_ent << endl ;    

	//-----------------------------------------------------
	// Resolution of Poisson equation for ln(N)
	//-----------------------------------------------------

	// Matter part of ln(N)
	// --------------------
	if (relativistic) {
	    source = a_car * (ener + 3*press) ;
	}
	else {
	    source = nbar ; 
	}
	
	(source.set()).set_dzpuis(4) ; 
	
	source.set_std_base() ;	    // Sets the standard spectral bases. 
	
	logn_auto.set_etat_qcq() ; 

	mpaff.poisson_falloff(source(), par_nul, logn_auto.set(), k_falloff) ; 

	// NB: at this stage logn_auto is in machine units, not in c^2

	// Quadratic part of ln(N)
	// -----------------------

	if (relativistic) {
	    
	    mpaff.dsdr(logn(), dlogn) ; 
	    mpaff.dsdr(beta_auto(), dbeta) ; 
	    
	    source = - dlogn * dbeta ; 
		      
	    logn_quad.set_etat_qcq() ; 
	    
	    mpaff.poisson_falloff(source(), par_nul, logn_quad.set(),
				  k_falloff) ; 
	    	    
	}

	//-----------------------------------------------------
	// Computation of the new radial scale
	//-----------------------------------------------------

	// alpha_r (r = alpha_r r') is determined so that the enthalpy
	// takes the requested value ent_b at the stellar surface
	
	double nu_mat0_b  = logn_auto()(l_b, k_b, j_b, i_b) ; 
	double nu_mat0_c  = logn_auto()(0, 0, 0, 0) ; 

	double nu_quad0_b  = logn_quad()(l_b, k_b, j_b, i_b) ; 
	double nu_quad0_c  = logn_quad()(0, 0, 0, 0) ; 

	double alpha_r2 = ( ent_c - ent_b - nu_quad0_b + nu_quad0_c )
				/ ( qpig*(nu_mat0_b - nu_mat0_c) ) ;

	alpha_r = sqrt(alpha_r2) ;
	
	// New radial scale
	mpaff.homothetie( alpha_r ) ; 

	//--------------------
	// First integral
	//--------------------

	// Gravitation potential in units c^2 :
	logn_auto = alpha_r2*qpig * logn_auto ;
	logn = logn_auto + logn_quad ;

	// Enthalpy in all space
	double logn_c = logn()(0, 0, 0, 0) ;
	ent = ent_c - logn() + logn_c ;

	//---------------------
	// Equation of state
	//---------------------
	
	equation_of_state() ; 
	
	if (relativistic) {
	    
	    //----------------------------
	    // Equation for beta = ln(AN)
	    //----------------------------
	    
	    mpaff.dsdr(logn(), dlogn) ; 
	    mpaff.dsdr(beta_auto(), dbeta) ; 
	    
	    source = 3 * qpig * a_car * press ;
	    
	    source = source() 
		     - 0.5 * (  dlogn * dlogn + dbeta * dbeta ) ;
	
	    source.set_std_base() ;   // Sets the standard spectral bases. 

	    beta_auto.set_etat_qcq() ; 
	         
	    mpaff.poisson_falloff(source(), par_nul, beta_auto.set(),
				  k_falloff) ; 
	    

	    // Metric coefficient A^2 update
	    
	    a_car = exp(2*(beta_auto - logn)) ;
   	    
	}
	
    // Relative difference with enthalpy at the previous step
    // ------------------------------------------------------
    
    diff_ent = norme( diffrel(ent(), ent_jm1()) ) / nzet ; 

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
    
    // ... hydro 
    ent.annule(nzet, nz-1) ;	// enthalpy set to zero at the exterior of 
				// the star
    ener_euler = ener ; 
    s_euler = 3 * press ;
    gam_euler = 1 ;  
    u_euler = 0 ; 
    
    // ... metric
    nnn = exp( unsurc2 * logn ) ; 
    nnn.set_std_base() ;
    shift = 0 ; 
    a_car.set_std_base() ;

    // Info printing
    // -------------
    
    cout << endl 
     << "Characteristics of the star obtained by Etoile::equil_spher_falloff : " 
     << endl 
     << "-------------------------------------------------------------------" 
     << endl ; 

    double ray = mp.val_r(l_b, 1., M_PI/2., 0) ; 
    cout << "Coordinate radius  :       " << ray / km << " km" << endl ; 

    double rcirc = ray * sqrt( a_car()(l_b, k_b, j_b, i_b) ) ; 
	
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
    source.set_std_base() ;	    
    double vir_mat = source().integrale() ; 
    
    //... Gravitational term

    Cmp tmp = beta_auto() - logn() ; 

    source =  - ( logn().dsdr() * logn().dsdr() 
		      - 0.5 * tmp.dsdr() * tmp.dsdr() ) 
		    * sqrt(a_car())  ; 

    source.set_std_base() ;

    double vir_grav = source().integrale() ;

    //... Relative error on the virial identity GRV3
    
    double grv3 = ( vir_mat + vir_grav ) / vir_mat ;

    cout << "Virial theorem GRV3 : " << endl ; 
    cout << "     3P term    : " << vir_mat << endl ; 
    cout << "     grav. term : " << vir_grav << endl ; 
    cout << "     relative error : " << grv3 << endl ; 
    
}
}
