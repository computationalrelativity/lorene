/*
 * Method of class Etoile to compute a static spherical configuration.
 *
 * (see file etoile.h for documentation).
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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
 * Revision 1.7  2016/12/05 16:17:54  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:52:59  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2008/11/14 13:48:06  e_gourgoulhon
 * Added parameter pent_limit to force the enthalpy values at the
 * boundaries between the domains, in case of more than one domain inside
 * the star.
 *
 * Revision 1.4  2004/05/07 12:13:15  k_taniguchi
 * Change the position of the initialization of alpha_r.
 *
 * Revision 1.3  2004/03/25 10:29:06  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.2  2003/04/23 15:09:38  j_novak
 * Standard basis is set to a_car and nnn before exiting.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/02/02  09:23:51  eric
 * Ajout du theoreme du viriel.
 * Affichage quantites globales a la fin.
 *
 * Revision 2.3  2000/01/28  17:18:36  eric
 * Modifs mineures.
 *
 * Revision 2.2  2000/01/27  16:47:16  eric
 * Premiere version qui tourne !
 *
 * Revision 2.1  2000/01/26  13:18:19  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/01/24  17:13:56  eric
 * *** empty log message ***
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
#include "nbr_spx.h"
#include "graphique.h"

namespace Lorene {
void Etoile::equilibrium_spher(double ent_c, double precis, const Tbl* pent_limit){
    
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

	mpaff.poisson(source(), par_nul, logn_auto.set()) ; 

	// NB: at this stage logn_auto is in machine units, not in c^2

	// Quadratic part of ln(N)
	// -----------------------

	if (relativistic) {
	    
	    mpaff.dsdr(logn(), dlogn) ; 
	    mpaff.dsdr(beta_auto(), dbeta) ; 
	    
	    source = - dlogn * dbeta ; 
		      
	    logn_quad.set_etat_qcq() ; 
	    
	    mpaff.poisson(source(), par_nul, logn_quad.set()) ; 
	    	    
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


        // One domain inside the star:
        // --------------------------- 
	if(nzet==1) { 

        mpaff.homothetie( alpha_r ) ; 

        } 
				      
	//--------------------
	// First integral
	//--------------------

	// Gravitation potential in units c^2 :
	logn_auto = alpha_r2*qpig * logn_auto ;
	logn = logn_auto + logn_quad ;

	// Enthalpy in all space
	double logn_c = logn()(0, 0, 0, 0) ;
	ent = ent_c - logn() + logn_c ;
	
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
        mp0.adapt(ent(), par_adapt) ; 

        //Map_af mpaff_prev (mpaff) ; 

        double alphal, betal ; 

        for(int l=0; l<nz; l++) { 

          alphal = mp0.get_alpha()[l] ; 
          betal  = mp0.get_beta()[l] ;
 
	  mpaff.set_alpha(alphal, l) ; 
          mpaff.set_beta(betal, l) ;

        } 
 

        //mbtest 
        int num_r1 = mg->get_nr(0) - 1;        
  
        cout << "Pressure difference:" << get_press()()(0,0,0,num_r1) - get_press()()(1,0,0,0) << endl ; 
	cout << "Difference in enthalpies at the domain boundary:" << endl ;
        cout << get_ent()()(0,0,0,num_r1) << endl ; 
        cout << get_ent()()(1,0,0,0) << endl ;

        cout << "Enthalpy difference: " << get_ent()()(0,0,0,num_r1) - get_ent()()(1,0,0,0) << endl ;

	// Computation of the enthalpy at the new grid points
        //----------------------------------------------------                     

        //mpaff.reevaluate(&mpaff_prev, nzet+1, ent.set()) ;

    } 

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
	
	    source.set_std_base() ;	    // Sets the standard spectral bases. 

	    beta_auto.set_etat_qcq() ; 
	         
	    mpaff.poisson(source(), par_nul, beta_auto.set()) ; 
	    

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
     << "Characteristics of the star obtained by Etoile::equilibrium_spher : " 
     << endl 
     << "-----------------------------------------------------------------" 
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
    
    if (nzet > 1) {
      cout.precision(10) ;

      for (int ltrans = 0; ltrans < nzet-1; ltrans++) {
	cout << endl << "Values at boundary between domains no. " << ltrans << " and " << 	ltrans+1 << " for theta = pi/2 and phi = 0 :" << endl ;

	double rt1 = mp.val_r(ltrans, 1., M_PI/2, 0.) ; 
	double rt2 = mp.val_r(ltrans+1, -1., M_PI/2, 0.) ; 
	cout << "   Coord. r [km] (left, right, rel. diff) : " 
	  << rt1 / km << "  " << rt2 / km << "  " << (rt2 - rt1)/rt1  << endl ; 
  
	int ntm1 = mg->get_nt(ltrans) - 1; 
	int nrm1 = mg->get_nr(ltrans) - 1 ; 
	double ent1 = ent()(ltrans, 0, ntm1, nrm1) ; 
	double ent2 = ent()(ltrans+1, 0, ntm1, 0) ; 
	  cout << "   Enthalpy (left, right, rel. diff) : " 
	    << ent1 << "  " << ent2 << "  " << (ent2-ent1)/ent1  << endl ; 

	double press1 = press()(ltrans, 0, ntm1, nrm1) ; 
	double press2 = press()(ltrans+1, 0, ntm1, 0) ; 
	  cout << "   Pressure (left, right, rel. diff) : " 
	  << press1 << "  " << press2 << "  " << (press2-press1)/press1 << endl ; 

	double nb1 = nbar()(ltrans, 0, ntm1, nrm1) ; 
	double nb2 = nbar()(ltrans+1, 0, ntm1, 0) ; 
	  cout << "   Baryon density (left, right, rel. diff) : " 
	  << nb1 << "  " << nb2 << "  " << (nb2-nb1)/nb1  << endl ; 
      }
    }


/*    double r_max = 1.2 * ray_eq() ; 
    des_profile(nbar(), 0., r_max, M_PI/2, 0., "n", "Baryon density") ; 
    des_profile(ener(), 0., r_max, M_PI/2, 0., "e", "Energy density") ; 
    des_profile(press(), 0., r_max, M_PI/2, 0., "p", "Pressure") ; 
    des_profile(ent(), 0., r_max, M_PI/2, 0., "H", "Enthalpy") ; 
*/
   
}
}
