/*
 *  Method of class Etoile to compute a static spherical configuration
 *   of a neutron star in a NS-BH binary system.
 *
 *  (see file etoile.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003 Keisuke Taniguchi
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

char name_of_this_file_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/10/20 12:20:26  k_taniguchi
 * Computation of a static spherical configuration of a neutron star
 * in a NS-BH binary system.
 *
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "etoile.h"
#include "param.h"

void Etoile::nsbh_equilspher(double ent_c, double precis){

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

    // Initial metric
    n_auto = 1. ;   // Initialization of the laspe function
    confpsi_auto = 1. ;   // Initialization of the conformal factor

    // Auxiliary quantities
    // --------------------

    // Affine mapping for solving the Poisson equations
    Map_af mpaff(mp) ;

    Param par_nul ;   // Param (null) for Map_af::poisson.

    Tenseur ent_jm1(mp) ;   // Enthalpy at previous step
    ent_jm1 = 0 ;

    Tenseur source(mp) ;
    Tenseur lapse(mp) ;
    Tenseur n_quad(mp) ;
    Tenseur logn(mp) ;

    lapse = 1. ;
    n_quad = 1. ;

    Cmp dlapse(mp) ;
    Cmp dconfpsi(mp) ;
    Tenseur confpsi_c(mp) ;

    double diff_ent = 1 ;
    int mermax = 200 ;   // Max number of iterations

    //==============================================================
    // 			Start of iteration
    //==============================================================

    for (int mer=0 ; (diff_ent > precis) && (mer<mermax) ; mer++ ) {

	double alpha_r = 1. ;

	cout << "-----------------------------------------------" << endl ;
	cout << "step: " << mer << endl ;
	cout << "alpha_r: " << alpha_r << endl ;
	cout << "diff_ent = " << diff_ent << endl ;

	//--------------------------------------------------
	// Resolution of Poisson equation for the lapse (N)
	//--------------------------------------------------

	// Matter part of N
	// ----------------
	if (relativistic) {
	    source = nnn * confpsi_q * (ener + 3.*press) ;
	}
	else {
	    source = nbar ;
	}

	(source.set()).set_dzpuis(4) ;

	source.set_std_base() ;   // Sets the standard spectral bases.

	n_auto.set_etat_qcq() ;

	mpaff.poisson(source(), par_nul, n_auto.set()) ;
	// NB: at this stage n_auto is in machine units, not in c^2

	// Quadratic part of N
	// -------------------
	if (relativistic) {

	    mpaff.dsdr(lapse(), dlapse) ;
	    mpaff.dsdr(confpsi(), dconfpsi) ;

	    source = -2. * dlapse * dconfpsi / confpsi ;

	    n_quad.set_etat_qcq() ;

	    mpaff.poisson(source(), par_nul, n_quad.set()) ;

	}

	//-------------------------------------
	// Computation of the new radial scale
	//-------------------------------------

	// alpha_r (r = alpha_r r') is determined so that the enthalpy
	// takes the requested value ent_b at the stellar surface

	double hhh_b = exp( ent_b ) ;
	double hhh_c = exp( ent_c ) ;

	double n_mat0_b  = n_auto()(l_b, k_b, j_b, i_b) ; 
	double n_mat0_c  = n_auto()(0, 0, 0, 0) ; 

	double n_quad0_b  = n_quad()(l_b, k_b, j_b, i_b) ; 
	double n_quad0_c  = n_quad()(0, 0, 0, 0) ; 

	double alpha_r2 = ( hhh_c * n_quad0_c - hhh_b * n_quad0_b )
	    / ( qpig*(hhh_b * n_mat0_b - hhh_c * n_mat0_c) ) ;

	alpha_r = sqrt(alpha_r2) ;
	
	// New radial scale
	mpaff.homothetie( alpha_r ) ;

	//--------------------------------------
	// First integral of the Euler equation
	//--------------------------------------

	// Gravitation potential in units c^2 :
	n_auto = alpha_r2 * qpig * n_auto ;
	lapse = n_auto + n_quad ;

	logn = log( lapse ) ;

	// Enthalpy in all space
	double logn_c = logn()(0, 0, 0, 0) ;
	ent = ent_c - logn() + logn_c ;

	//-------------------
	// Equation of state
	//-------------------

	equation_of_state() ;

	if (relativistic) {

	    //---------------------------------------------
	    // Equation for the conformal factor (confpsi)
	    //---------------------------------------------

	    confpsi_c = pow(confpsi, 5.) ;
	    source = - 0.5 * qpig * confpsi_c * ener ;

	    source.set_std_base() ;   // Sets the standard spectral bases.

	    confpsi_auto.set_etat_qcq() ;

	    mpaff.poisson(source(), par_nul, confpsi_auto.set()) ;

	}

	// Relative difference with enthalpy at the previous step
	// ------------------------------------------------------

	diff_ent = norme( diffrel(ent(), ent_jm1()) ) / nzet ;

	// Next step
	// ---------

	ent_jm1 = ent ;

    }    // End of iteration loop 

    //==========================================================
    // 			End of iteration
    //==========================================================

    // The mapping is transfered to that of the star:
    // ----------------------------------------------
    mp = mpaff ;

    // Sets value to all the Tenseur's of the star
    // -------------------------------------------

    // ... hydro
    ent.annule(nzet, nz-1) ;   // enthalpy set to zero at the exterior of
                               // the star
    ener_euler = ener ;
    s_euler = 3 * press ;
    gam_euler = 1 ;
    u_euler = 0 ;

    // ... metric
    nnn = lapse ;
    nnn.set_std_base() ;
    shift = 0 ;
    confpsi.set_std_base() ;

    // Info printing
    // -------------

    cout << endl
	 << "Characteristics of the star obtained by Etoile::nsbh_equilspher : "
	 << endl
	 << "-----------------------------------------------------------------"
	 << endl ;

    double ray = mp.val_r(l_b, 1., M_PI/2., 0) ;
    cout << "Coordinate radius  :       " << ray / km << " km" << endl ;

    double rcirc = ray * sqrt( confpsi_q()(l_b, k_b, j_b, i_b) ) ;

    double compact = qpig/(4.*M_PI) * mass_g() / rcirc ;

    cout << "Circumferential radius R : " << rcirc/km  << " km" << endl ;
    cout << "Baryon mass :	     " << mass_b()/msol << " Mo" << endl ;
    cout << "Gravitational mass M :  " << mass_g()/msol << " Mo" << endl ;
    cout << "Compacity parameter GM/(c^2 R) : " << compact << endl ;

    //-----------------
    // Virial theorem
    //-----------------

    //... Pressure term

    source = qpig * pow(confpsi, 6.) * s_euler ;
    source.set_std_base() ;
    double vir_mat = source().integrale() ;

    //... Gravitational term

    Cmp tmp = beta_auto() - logn() ;

    source =  - ( logn().dsdr() * logn().dsdr()
		  - 0.5 * tmp.dsdr() * tmp.dsdr() )
	* sqrt(a_car()) ;

    source.set_std_base() ;
    double vir_grav = source().integrale() ;

    //... Relative error on the virial identity GRV3

    double grv3 = ( vir_mat + vir_grav ) / vir_mat ;

    cout << "Virial theorem GRV3 : " << endl ;
    cout << "     3P term    : " << vir_mat << endl ;
    cout << "     grav. term : " << vir_grav << endl ;
    cout << "     relative error : " << grv3 << endl ;

}
