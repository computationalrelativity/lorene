/*
 * Method of class Etoile to compute a static spherical configuration
 *  by regularizing source.
 *
 * (see file etoile.h for documentation).
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


char et_equil_spher_regu_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.14  2001/03/06  16:34:04  keisuke
 * Change the regularization degree k_div=1 to the arbitrary one.
 *
 * Revision 2.13  2001/02/07  09:46:11  eric
 * unsgam1 est desormais donne par Eos::der_nbar_ent (cas newtonien)
 *   ou Eos::der_ener_ent (cas relativiste).
 *
 * Revision 2.12  2000/09/26  15:42:35  keisuke
 * Correction erreur: la triade de duu_div doit etre celle de mp et non celle
 *  de l'objet temporaire mpaff.
 *
 * Revision 2.11  2000/09/26  06:56:04  keisuke
 * Suppress "int" from the declaration of k_div.
 *
 * Revision 2.10  2000/09/22  15:53:36  keisuke
 * d_logn_auto_div est desormais un membre de la classe Etoile.
 *
 * Revision 2.9  2000/09/08  12:23:17  keisuke
 * Correct a typological error in the equation.
 *
 * Revision 2.8  2000/09/07  15:37:23  keisuke
 * Add a new argument in poisson_regular
 *  and suppress logn_auto_total.
 *
 * Revision 2.7  2000/09/06  12:47:35  keisuke
 * Suppress #include "graphique.h".
 *
 * Revision 2.6  2000/09/06  12:39:59  keisuke
 * Save the map in the every iterative step after operating "homothetie".
 *
 * Revision 2.5  2000/09/04  16:15:05  keisuke
 * Change the argument of Map_af::poisson_regular.
 *
 * Revision 2.4  2000/08/31  16:05:26  keisuke
 * Modify the arguments of Map_af::poisson_regular.
 *
 * Revision 2.3  2000/08/29  11:38:25  eric
 * Ajout des membres k_div et logn_auto_div a la classe Etoile.
 *
 * Revision 2.2  2000/08/28  16:10:39  keisuke
 * Add "nzet" in the argumant of poisson_regular.
 *
 * Revision 2.1  2000/08/25  14:58:15  keisuke
 * Modif (Virial theorem).
 *
 * Revision 2.0  2000/08/25  09:01:33  keisuke
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <math.h>

// Headers Lorene
#include "etoile.h"
#include "eos.h"
#include "param.h"

//********************************************************************

void Etoile::equil_spher_regular(double ent_c, double precis){

    // Fundamental constants and units
    // -------------------------------
    #include "unites.h"
    // To avoid some compilation warnings
    if (ent_c < 0) {
      cout << f_unit << msol << mevpfm3 <<endl ;
    }

    // Initializations
    // ---------------

    //    k_div = 1 ;	    // Regularization index

    cout << "Input the regularization degree (k_div) : " ;
    cin >> k_div ;          // Regularization index

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
    Map_af mpaff(mp) ;	

    Param par_nul   ;	 // Param (null) for Map_af::poisson.

    Tenseur ent_jm1(mp) ;	// Enthalpy at previous step
    ent_jm1 = 0 ;

    Tenseur source(mp) ;
    Tenseur logn(mp) ;
    Tenseur logn_quad(mp) ;
    logn = 0 ;
    logn_quad = 0 ;

    Tenseur dlogn_auto_regu(mp, 1, COV, mp.get_bvect_spher()) ;

    Cmp source_regu(mp) ;
    Cmp source_div(mp) ;

    Cmp dlogn(mp) ;
    dlogn = 0 ;
    Cmp dbeta(mp) ;

    Tenseur dlogn_auto(mp, 1, COV, mp.get_bvect_spher()) ;
    Tenseur dlogn_quad(mp) ;
    dlogn_quad = 0 ;

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

      double unsgam1 ;    // effective power of H in the source 
			  // close to the surface

      if (relativistic) {
	
	source = a_car * (ener + 3*press) ;

	// 1/(gam-1) = dln(e)/dln(H) close to the surface : 
	unsgam1 = eos.der_ener_ent_p(ent_b + 1e-10*(ent_c-ent_b)) ; 
      }
      else {
	
	source = nbar ;
	
	// 1/(gam-1) = dln(n)/dln(H) close to the surface : 
	unsgam1 = eos.der_nbar_ent_p(ent_b + 1e-10*(ent_c-ent_b)) ; 
      }

      (source.set()).set_dzpuis(4) ; 

      source.set_std_base() ;	    // Sets the standard spectral bases. 

      logn_auto.set_etat_qcq() ;
      logn_auto_regu.set_etat_qcq() ;
      logn_auto_div.set_etat_qcq() ;

      source_regu.std_base_scal() ;
      source_div.std_base_scal() ;

      mpaff.poisson_regular(source(), k_div, nzet, unsgam1, par_nul,
			    logn_auto.set(), logn_auto_regu.set(),
			    logn_auto_div.set(),
			    d_logn_auto_div, source_regu, source_div) ;

      dlogn_auto_regu = logn_auto_regu.gradient_spher() ;

      dlogn_auto = dlogn_auto_regu + d_logn_auto_div ;

      // NB: at this stage logn_auto is in machine units, not in c^2

      // Quadratic part of ln(N)
      // -----------------------

      if (relativistic) {

	mpaff.dsdr(beta_auto(), dbeta) ;

 	source = - dlogn * dbeta ;

	logn_quad.set_etat_qcq() ;

	mpaff.poisson(source(), par_nul, logn_quad.set()) ; 

	dlogn_quad.set_etat_qcq() ;

	mpaff.dsdr(logn_quad(), dlogn_quad.set()) ;

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
      // -----------------
      mpaff.homothetie( alpha_r ) ;

      // The mapping is transfered to that of the star:
      // ----------------------------------------------
      mp = mpaff ;

      d_logn_auto_div.set_triad( mp.get_bvect_spher() ) ;  // Absolutely necessary !!!

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

      // derivative of gravitation potential in units c^2 :
      dlogn_auto = alpha_r*qpig * dlogn_auto ;
      dlogn = dlogn_auto(0) + dlogn_quad() ;

      if (relativistic) {

	//----------------------------
	// Equation for beta = ln(AN)
	//----------------------------

	mpaff.dsdr(beta_auto(), dbeta) ;

	source = 3 * qpig * a_car * press ;

	source = source() 
	  - 0.5 * (  dlogn * dlogn + dbeta * dbeta ) ;

	source.set_std_base() ;    // Sets the standard spectral bases.

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
    shift = 0 ;

    // Info printing
    // -------------

    cout << endl
     << "Characteristics of the star obtained by Etoile::equil_spher_regular : "
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

    Cmp tmp = beta_auto().dsdr() - dlogn ;

    source =  - ( dlogn * dlogn - 0.5 * tmp * tmp ) * sqrt(a_car())  ;

    source.set_std_base() ;
    double vir_grav = source().integrale() ;

    //... Relative error on the virial identity GRV3

    double grv3 = ( vir_mat + vir_grav ) / vir_mat ;

    cout << "Virial theorem GRV3 : " << endl ;
    cout << "     3P term    : " << vir_mat << endl ;
    cout << "     grav. term : " << vir_grav << endl ;
    cout << "     relative error : " << grv3 << endl ;


}
