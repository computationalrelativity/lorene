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


char star_equil_spher_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
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
#include "star.h"
#include "param.h"
#include "graphique.h"

void Star::equilibrium_spher(double ent_c, double precis){
    
    // Fundamental constants and units
    // -------------------------------
    #include "unites.h"	    
    // To avoid some compilation warnings
    if (ent_c < 0) {
	cout << f_unit << msol << mevpfm3 << endl ; 
    }    
    
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
    qq = 1 ;
    qq.std_spectral_base() ;

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
    Scalar dqq(mp) ; 
    
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
	
	source.std_spectral_base() ;	    // Sets the standard spectral bases
 	
	Cmp source_logn_mat (source) ;
	Cmp logn_mat_cmp (logn_mat) ;
	logn_mat_cmp.set_etat_qcq() ;
	
	mpaff.poisson(source_logn_mat, par_nul, logn_mat_cmp) ; 

	logn_mat = logn_mat_cmp ;

	// NB: at this stage logn is in machine units, not in c^2

	// Quadratic part of ln(N)
	// -----------------------

	dlogn = logn.dsdr() ; 
	dqq = qq.dsdr() ; 
		
	source = - dlogn * dqq / qq ; 

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
	
	double nu_mat0_b  = logn_mat.point(l_b, k_b, j_b, i_b) ; 
	double nu_mat0_c  = logn_mat.point(0, 0, 0, 0) ; 

	double nu_quad0_b  = logn_quad.point(l_b, k_b, j_b, i_b) ; 
	double nu_quad0_c  = logn_quad.point(0, 0, 0, 0) ; 

	double alpha_r2 = ( ent_c - ent_b - nu_quad0_b + nu_quad0_c )
				/ ( qpig*(nu_mat0_b - nu_mat0_c) ) ;

	alpha_r = sqrt(alpha_r2) ;
	
	// New radial scale
	mpaff.homothetie( alpha_r ) ; 

	//--------------------
	// First integral
	//--------------------

	// Gravitation potential in units c^2 :
	logn_mat = alpha_r2*qpig * logn_mat ;
	logn = logn_mat + logn_quad ;

	// Enthalpy in all space
	double logn_c = logn.point(0, 0, 0, 0) ;
	ent = ent_c - logn + logn_c ;

	//---------------------
	// Equation of state
	//---------------------
	
	equation_of_state() ; 
	
	//---------------------
	// Equation for qq_auto
	//---------------------
	
	dlogn = logn.dsdr() ; 
	dqq = qq.dsdr() ; 
	
	source = 3 * qpig * a_car * qq * press ;
	
	source = source + ( dqq * dqq / qq - qq * dlogn * dlogn ) / 2. ;

	source.std_spectral_base() ;	  // Sets the standard spectral bases. 

	Cmp source_qq (source) ;
	Cmp qq_cmp (logn_quad) ;
	qq_cmp.set_etat_qcq() ;
	    
	mpaff.poisson(source_qq, par_nul, qq_cmp) ; 

	qq = qq_cmp ;

	qq = qq + 1 ;
    
      	qq.std_spectral_base() ;

	// Metric coefficient psi4 update
	
	nnn = exp( logn ) ; 
	a_car = qq * qq / ( nnn * nnn ) ;


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
    
    // ... hydro 
    ent.annule(nzet, nz-1) ;	// enthalpy set to zero at the exterior of 
				// the star
    ener_euler = ener ; 
    s_euler = 3 * press ;
    gam_euler = 1 ;  
    for(int i=1; i<=3; i++) u_euler.set(i) = 0 ; 
    
    // ... metric
    nnn = exp( logn ) ;
    nnn.std_spectral_base() ;
    for(int i=1; i<=3; i++) shift.set(i) = 0 ; 
 
    Sym_tensor tens_gamma(mp, COV, mp.get_bvect_spher()) ;
    for (int i=1; i<=3; i++)
	for (int j=1; j<=3; j++){
	    if (i == j) tens_gamma.set(i,i) = a_car ;
	    else tens_gamma.set(i,j) = 0. ;
	}
    gamma = tens_gamma ;

 
    // Info printing
    // -------------
    
    cout << endl 
     << "Characteristics of the star obtained by Etoile::equilibrium_spher : " 
     << endl 
     << "-----------------------------------------------------------------" 
     << endl ; 

    double ray = mp.val_r(l_b, 1., M_PI/2., 0) ; 
    cout << "Coordinate radius  :       " << ray / km << " km" << endl ; 

    double rcirc = ray * sqrt(a_car.point(l_b, k_b, j_b, i_b) ) ; 
	
    double compact = qpig/(4.*M_PI) * mass_g() / rcirc ; 

    cout << "Circumferential radius R : " << rcirc/km  << " km" << endl ;
    cout << "Baryon mass :	     " << mass_b()/msol << " Mo" << endl ;
    cout << "Gravitational mass M :  " << mass_g()/msol << " Mo" << endl ;
    cout << "Compacity parameter GM/(c^2 R) : " << compact << endl ;
     
}
