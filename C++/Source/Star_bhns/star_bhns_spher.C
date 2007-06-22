/*
 *  Method of class Star_bhns to compute a spherical star configuration
 *
 *    (see file star_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Keisuke Taniguchi
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

char star_bhns_spher_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:32:19  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
#include <math.h>

// Lorene headers
#include "star_bhns.h"
#include "param.h"
#include "cmp.h"
#include "tenseur.h"
#include "unites.h"

void Star_bhns::equil_spher_bhns(double ent_c, double precis) {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    // Initializations
    // ---------------

    const Mg3d* mg = mp.get_mg() ;
    int nz = mg->get_nzone() ;

    // Index of the point at phi=0, theta=pi/2 at the surface of the star:
    int l_b = nzet - 1 ;
    int i_b = mg->get_nr(l_b) - 1 ;
    int j_b = mg->get_nt(l_b) - 1 ;
    int k_b = 0 ;

    // Value of the enthalpy defining the surface of the star
    double ent_b = 0. ;

    // Initialization of the enthalpy field to the constant value ent_c :
    ent = ent_c ;
    ent.annule(nzet, nz-1) ;

    // Corresponding profiles of baryon density, energy density and pressure
    equation_of_state() ;

    // Initial metric
    lapse_auto = 1. ;
    lapse_auto.std_spectral_base() ;
    confo_auto = 1. ;
    confo_auto.std_spectral_base() ;

    // Auxiliary quantities
    // --------------------

    // Affine mapping for solving the Poisson equations
    Map_af mpaff(mp) ;

    Param par_nul ;  // Param (null) for Map_af::poisson

    Scalar ent_jm1(mp) ;  // Enthalpy at previous step
    ent_jm1 = ent_c ;
    ent_jm1.std_spectral_base() ;
    ent_jm1.annule(nzet, nz-1) ;

    Scalar source_lapse(mp) ;
    Scalar source_confo(mp) ;
    Scalar lapse_mat_m1(mp) ;
    Scalar lapse_quad_m1(mp) ;
    Scalar confo_auto_m1(mp) ;
    lapse_mat_m1 = 0. ;
    lapse_quad_m1 = 0. ;
    confo_auto_m1 = 0. ;
    Scalar lconfo(mp) ;
    Scalar dlapse(mp) ;
    Scalar dlconf(mp) ;

    double diff_ent = 1. ;
    int mermax = 200 ;  // Maximum number of iterations

    double alpha_r = 1. ;

    //==========================================================//
    //                    Start of iteration                    //
    //==========================================================//

    for (int mer=0 ; (diff_ent > precis) && (mer<mermax) ; mer++ ) {

      cout << "-----------------------------------------------" << endl ;
      cout << "step: " << mer << endl ;
      cout << "alpha_r: " << alpha_r << endl ;
      cout << "diff_ent = " << diff_ent << endl ;

      //---------------------------------------------------
      // Resolution of Poisson equation for lapse function
      //---------------------------------------------------

      // Matter part of lapse
      // --------------------
      source_lapse = lapse_auto * pow(confo_auto,4.) * (ener + 3.*press) ;

      source_lapse.inc_dzpuis(4-source_lapse.get_dzpuis()) ;
      source_lapse.std_spectral_base() ;

      Cmp sou_lap_mat(source_lapse) ;
      Cmp lap_mat_cmp(lapse_mat_m1) ;
      lap_mat_cmp.set_etat_qcq() ;

      mpaff.poisson(sou_lap_mat, par_nul, lap_mat_cmp) ;

      // Re-construction of a scalar
      lapse_mat_m1 = lap_mat_cmp ;

      // Quadratic part of lapse
      // -----------------------
      confo_auto.set_etat_qcq() ;
      lconfo = log(confo_auto) ;
      lconfo.std_spectral_base() ;

      lapse_auto.set_etat_qcq() ;
      lconfo.set_etat_qcq() ;

      mpaff.dsdr(lapse_auto, dlapse) ;
      mpaff.dsdr(lconfo, dlconf) ;

      source_lapse = - 2. * dlapse * dlconf ;
      source_lapse.inc_dzpuis(4-source_lapse.get_dzpuis()) ;
      source_lapse.std_spectral_base() ;

      Cmp sou_lap_quad(source_lapse) ;
      Cmp lap_quad_cmp(lapse_quad_m1) ;
      lap_quad_cmp.set_etat_qcq() ;

      mpaff.poisson(sou_lap_quad, par_nul, lap_quad_cmp) ;

      // Re-construction of a scalar
      lapse_quad_m1 = lap_quad_cmp ;

      //-------------------------------------
      // Computation of the new radial scale
      //-------------------------------------

      double exp_ent_c = exp(ent_c) ;
      double exp_ent_b = exp(ent_b) ;

      double lap_mat_c = lapse_mat_m1.val_grid_point(0,0,0,0) ;
      double lap_mat_b = lapse_mat_m1.val_grid_point(l_b,k_b,j_b,i_b) ;

      // +1 comes from the total lapse
      double lap_quad_c = lapse_quad_m1.val_grid_point(0,0,0,0) + 1. ;
      double lap_quad_b = lapse_quad_m1.val_grid_point(l_b,k_b,j_b,i_b) + 1. ;

      double alpha_r2 = (exp_ent_b*lap_quad_b - exp_ent_c*lap_quad_c)
	/ ( qpig*(exp_ent_c*lap_mat_c - exp_ent_b*lap_mat_b) ) ;

      alpha_r = sqrt(alpha_r2) ;

      // New radial scale
      mpaff.homothetie( alpha_r ) ;

      //----------------
      // First integral
      //----------------

      // Lapse function
      lapse_mat_m1 = alpha_r2 * qpig * lapse_mat_m1 ;
      lapse_auto = lapse_mat_m1 + lapse_quad_m1 + 1. ;

      // Enthalpy in all space
      double lap_c = lapse_auto.val_grid_point(0,0,0,0) ;
      ent = ent_c + log(lap_c) - log(lapse_auto) ;
      ent.std_spectral_base() ;

      //-------------------
      // Equation of state
      //-------------------

      equation_of_state() ;

      //-----------------------------------------------------
      // Resolution of Poisson equation for conformal factor
      //-----------------------------------------------------

      source_confo = - 0.5 * qpig * pow(confo_auto,5.) * ener ;
      source_confo.inc_dzpuis(4-source_confo.get_dzpuis()) ;
      source_confo.std_spectral_base() ;

      Cmp sou_confo(source_confo) ;
      Cmp cnf_auto_cmp(confo_auto_m1) ;
      cnf_auto_cmp.set_etat_qcq() ;

      mpaff.poisson(sou_confo, par_nul, cnf_auto_cmp) ;

      // Re-construction of a scalr
      confo_auto_m1 = cnf_auto_cmp ;

      confo_auto = confo_auto_m1 + 1. ;

      // Relative difference with enthalphy at the previous step
      // -------------------------------------------------------

      diff_ent = norme( diffrel(ent, ent_jm1) ) / nzet ;

      // Next step
      // ---------

      ent_jm1 = ent ;

    } // End of iteration loop

    //========================================================//
    //                    End of iteration                    //
    //========================================================//

    // The mapping is transfered to that of the star
    // ---------------------------------------------
    mp = mpaff ;

    // Sets values
    // -----------

    // ... hydro
    ent.annule(nzet, nz-1) ;

    ener_euler = ener ;
    s_euler = 3. * press ;
    gam_euler = 1. ;
    for(int i=1; i<=3; i++)
      u_euler.set(i) = 0 ;

    // ... metric
    lapse_tot = lapse_auto ;
    confo_tot = confo_auto ;
    psi4 = pow(confo_auto, 4.) ;
    for (int i=1; i<=3; i++)
      shift_auto.set(i) = 0. ;

    // Info printing
    // -------------

    cout << endl
	 << "Characteristics of the star obtained by Star_bhns::equil_spher_bhns : "
	 << endl
	 << "-------------------------------------------------------------------   "
	 << endl ;

    cout.precision(16) ;
    double ray = mp.val_r(l_b, 1., M_PI/2., 0) ;
    cout << "Coordinate radius :               "
	 << ray / km << " [km]" << endl ;

    double rcirc = ray * sqrt(psi4.val_grid_point(l_b, k_b, j_b, i_b) ) ;
    double compact = qpig/(4.*M_PI) * mass_g_bhns() / rcirc ;

    cout << "Circumferential radius R :        "
	 << rcirc/km  << " [km]" << endl ;
    cout << "Baryon mass :                     "
	 << mass_b_bhns(0,0.,1.)/msol << " [Mo]" << endl ;
    cout << "Gravitational mass M :            "
	 << mass_g_bhns()/msol << " [Mo]" << endl ;
    cout << "Compaction parameter GM/(c^2 R) : " << compact << endl ;

    //----------------
    // Virial theorem
    //----------------

    //... Pressure term
    Scalar source(mp) ;
    source = qpig * pow(confo_auto,6.) * s_euler ;
    source.std_spectral_base() ;
    double vir_mat = source.integrale() ;

    //... Gravitational term
    Scalar tmp1(mp) ;
    tmp1 = log(lapse_auto) ;
    tmp1.std_spectral_base() ;

    Scalar tmp2(mp) ;
    tmp2 = log(confo_auto) ;
    tmp2.std_spectral_base() ;

    source = confo_auto * confo_auto
      * ( 2. * tmp2.dsdr() * tmp2.dsdr() - tmp1.dsdr() * tmp1.dsdr() ) ;
    source.std_spectral_base() ;	    
    double vir_grav = source.integrale() ;

    //... Relative error on the virial identity GRV3
    double grv3 = ( vir_mat + vir_grav ) / vir_mat ;

    cout << "Virial theorem GRV3 : " << endl ;
    cout << "     3P term :        " << vir_mat << endl ;
    cout << "     grav. term :     " << vir_grav << endl ;
    cout << "     relative error : " << grv3 << endl ;

}
