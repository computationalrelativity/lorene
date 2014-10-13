/*
 *  Method of class Star_bhns to compute a spherical star configuration
 *
 *    (see file star_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005,2007 Keisuke Taniguchi
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
 * Revision 1.4  2014/10/13 08:53:41  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:16  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/05/15 19:16:54  k_taniguchi
 * Change of some parameters.
 *
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
#include <cmath>

// Lorene headers
#include "star_bhns.h"
#include "param.h"
#include "cmp.h"
#include "tenseur.h"
#include "unites.h"

namespace Lorene {
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
    lapconf_auto = 1. ;
    lapconf_auto.std_spectral_base() ;
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

    Scalar source_lapconf(mp) ;
    Scalar source_confo(mp) ;
    Scalar lapconf_auto_m1(mp) ;
    Scalar confo_auto_m1(mp) ;
    lapconf_auto_m1 = 0. ;
    confo_auto_m1 = 0. ;

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

      //----------------------------------------------------
      // Resolution of Poisson equation for lapconf function
      //----------------------------------------------------

      // Matter part of lapconf
      // ----------------------
      source_lapconf = qpig * lapconf_auto * pow(confo_auto,4.)
	* (0.5*ener + 3.*press) ;

      source_lapconf.inc_dzpuis(4-source_lapconf.get_dzpuis()) ;
      source_lapconf.std_spectral_base() ;

      Cmp sou_lapconf(source_lapconf) ;
      Cmp lapconf_cmp(lapconf_auto_m1) ;
      lapconf_cmp.set_etat_qcq() ;

      mpaff.poisson(sou_lapconf, par_nul, lapconf_cmp) ;

      // Re-construction of a scalar
      lapconf_auto_m1 = lapconf_cmp ;

      //-------------------------------------
      // Computation of the new radial scale
      //-------------------------------------

      double exp_ent_c = exp(ent_c) ;
      double exp_ent_b = exp(ent_b) ;

      double lap_auto_c = lapconf_auto_m1.val_grid_point(0,0,0,0) ;
      double lap_auto_b = lapconf_auto_m1.val_grid_point(l_b,k_b,j_b,i_b) ;

      double confo_c = confo_auto.val_grid_point(0,0,0,0) ;
      double confo_b = confo_auto.val_grid_point(l_b,k_b,j_b,i_b) ;

      double alpha_r2 = (exp_ent_b*confo_c - exp_ent_c*confo_b)
	/ ( exp_ent_c*confo_b*lap_auto_c
	    - exp_ent_b*confo_c*lap_auto_b )  ;

      alpha_r = sqrt(alpha_r2) ;

      // New radial scale
      mpaff.homothetie( alpha_r ) ;

      //----------------
      // First integral
      //----------------

      // Lapconf function
      lapconf_auto_m1 = alpha_r2 * lapconf_auto_m1 ;
      lapconf_auto = lapconf_auto_m1 + 1. ;

      // Enthalpy in all space
      double lapconfo_c = lapconf_auto.val_grid_point(0,0,0,0) ;
      confo_c = confo_auto.val_grid_point(0,0,0,0) ;
      ent = ent_c + log(lapconfo_c) - log(confo_c)
	- log(lapconf_auto) + log(confo_auto) ;
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
    lapconf_tot = lapconf_auto ;
    lapse_auto = lapconf_auto / confo_auto ;
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
    tmp1 = confo_auto.dsdr() ;
    tmp1.std_spectral_base() ;

    Scalar tmp2(mp) ;
    tmp2 = confo_auto * lapconf_auto.dsdr() / lapconf_auto  - tmp1 ;
    tmp2.std_spectral_base() ;

    source = 2. * tmp1 * tmp1 - tmp2 * tmp2 ;

    /*
    Scalar tmp1(mp) ;
    tmp1 = log(lapse_auto) ;
    tmp1.std_spectral_base() ;

    Scalar tmp2(mp) ;
    tmp2 = log(confo_auto) ;
    tmp2.std_spectral_base() ;

    source = confo_auto * confo_auto
      * ( 2. * tmp2.dsdr() * tmp2.dsdr() - tmp1.dsdr() * tmp1.dsdr() ) ;
    */
    source.std_spectral_base() ;	    
    double vir_grav = source.integrale() ;

    //... Relative error on the virial identity GRV3
    double grv3 = ( vir_mat + vir_grav ) / vir_mat ;

    cout << "Virial theorem GRV3 : " << endl ;
    cout << "     3P term :        " << vir_mat << endl ;
    cout << "     grav. term :     " << vir_grav << endl ;
    cout << "     relative error : " << grv3 << endl ;

}
}
