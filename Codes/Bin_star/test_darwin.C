/*
 * Test program for solving Darwin-like binaries
 * Fixing the orbital separation and the central value of the star
 * The case of the identical star binary
 */
 
/*
 *   Copyright (c) 1999-2001 Keisuke Taniguchi
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

char test_darwin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:31  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.3  1999/12/24  12:15:51  keisuke
 * Minor changes.
 *
 * Revision 1.2  1999/12/23  16:41:11  keisuke
 * Change of the criterion for the convergence.
 *
 * Revision 1.1  1999/12/23  13:44:55  keisuke
 * Initial revision
 *
 *
 * $Header$
 *
 */

// version of 09.12.1999
// version of 16.12.1999
// version of 22.12.1999
// version of 23.12.1999

// headers C++
#include <iostream.h>
#include <fstream.h>

// headers C
#include <stdlib.h>
#include <math.h>

// headers Lorene
#include "type_parite.h"
#include "cmp.h"
#include "utilitaires.h"
#include "graphique.h"
#include "param.h"

// Local prototypes:
Cmp eos1_local(const Cmp& ent, int nzet, double n_index) ;
Cmp eos2_local(const Cmp& ent, int nzet, double n_index) ;
double funct_zero_ent(double r, const Param& par) ;

//**********************************************************************

void main(){
    
    // Identification of all the subroutines called by the code : 
    
    system("ident test_darwin > identif_darwin.d") ; 

    //-----------------------------------------------------------------------
    //			Input from file "part.d"
    //-----------------------------------------------------------------------
    
    int type_t, type_p, nt, np, nz, l;
    char blabla[80];
    double n_index, R_sep ;

    //----------------------------------------------------
    //             Declear the polytropic index
    //----------------------------------------------------
    cout << "Input the polytropic index: " ;
    cin >> n_index ;
    cout << "Input the orbital separation: " ;
    cin >> R_sep ;

    ifstream fich("part.d") ;
    fich >> nt; fich.getline(blabla, 80);
    fich >> np; fich.getline(blabla, 80);
    fich >> nz; fich.getline(blabla, 80);

    cout << "nb de points en phi : np = " << np << endl;
    cout << "nb de points en theta : nt = " << nt << endl;
    cout << "nb de zones : nz = " << nz << endl;

    // initialisation des tableaux decrivant chaque zone 

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
    int* type_r = new int[nz];
     
    fich >> bornes[nz] ; fich.getline(blabla, 80) ;
    for (l=0; l<nz; l++) {
	fich >> nr[l]; 
	fich >> bornes[l]; 
	fich >> type_r[l]; fich.getline(blabla, 80);
	np_tab[l] = np ; 
	nt_tab[l] = nt ; 
    }
    if (type_r[nz-1]==UNSURR) bornes[nz] = 1./bornes[nz] ;

    fich.close();
   

    //-----------------------------------------------------------------------
    //		Construction of 2 multi-grids
    //-----------------------------------------------------------------------
    
    type_t = SYM ; 
    type_p = NONSYM ; 

    const Mg3d mg1(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;
    const Mg3d mg2(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

     //-----------------------------------------------------------------------
    //		Construction of 2 mappings
    //-----------------------------------------------------------------------
    
    Map_af mp1(mg1, bornes) ;
    double xo1 = - R_sep/2 ; 
    mp1.set_ori(xo1, 0, 0) ; 

    Map_af mp2(mg2, bornes) ;
    double xo2 = R_sep/2 ; 
    mp2.set_ori(xo2, 0, 0) ; 
    mp2.set_rot_phi(M_PI) ; 
    
    
    cout << "Mapping mp1 : " << mp1 << endl ; 
    cout << "Mapping mp2 : " << mp2 << endl ; 

    //-----------------------------------------------------------------------
    //		Construction of Cmp's
    //-----------------------------------------------------------------------

    //-------------------------------------
    // Construction of two spherical stars
    //-------------------------------------

    Cmp rho1(mp1) ; 
    const Coord& x1 = mp1.x ; 
    const Coord& y1 = mp1.y ;
    
    rho1 = 1  ;
    rho1.annule(nz-1) ;
    rho1.set_dzpuis(4) ; 
    rho1.std_base_scal() ;

    Cmp rho2(mp2) ;
    const Coord& x2 = mp2.x ;
    const Coord& y2 = mp2.y ;

    rho2 = 1 ;
    rho2.annule(nz-1) ;
    rho2.set_dzpuis(4) ;
    rho2.std_base_scal() ;

    Cmp pot_star1(mp1) ;
    Cmp rho_prev_star1(mp1) ;
    Cmp ent_star1(mp1) ;
    double ent_star1_c = 1 ;
    double ent_star1_s = 0 ;

    Cmp pot_star2(mp2) ;
    Cmp rho_prev_star2(mp2) ;
    Cmp ent_star2(mp2) ;
    double ent_star2_c = 1 ;
    double ent_star2_s = 0 ;

    int nr_star1 = mg1.get_nr(nz-2) ;
    int nr_star2 = mg2.get_nr(nz-2) ;

    // Start the do loop for constructing spherical stars
    //----------------------------------------------------

    // For star 1
    //------------
    do {

      // Solving the Poisson equation
      //------------------------------
      pot_star1 = - rho1.poisson() ;

      // Obtaining "lambda" in order to rescale the mapping
      //----------------------------------------------------
      double pot_star1_c = pot_star1(0,0,0,0) ;
      double pot_star1_s = pot_star1(nz-2,0,0,nr_star1-1) ;

      double lambda2_star1 = (ent_star1_s - ent_star1_c) /
	(pot_star1_s - pot_star1_c) ; 
      cout << "lambda2_star1 : " << lambda2_star1 << endl ; 

      // New potential after rescaling
      //-------------------------------
      pot_star1 = lambda2_star1 * pot_star1 ; 
      pot_star1_c *= lambda2_star1 ; 

      // First integral of motion
      // -------------------------
      ent_star1 = ent_star1_c + pot_star1 - pot_star1_c ;
      // enthalpy is just the same as the Lane-Emden function Theta

      // Rescale of the mapping
      //------------------------
      mp1.homothetie( sqrt(lambda2_star1) ) ; 

      // EOS : rho = rho(H) 
      //--------------------
      rho_prev_star1 = rho1 ;
      rho1 = pow(abs(ent_star1), n_index) ; 

      rho1.annule(nz-1) ;
      rho1.set_dzpuis(4) ;
      rho1.std_base_scal() ;

      cout << max(abs(rho1 - rho_prev_star1)) << endl;
      cout << "Maximum difference between rho1 and rho1_prev : "
	   << max(diffrel(rho1, rho_prev_star1)) << endl ;

    } while(max(diffrel(rho1, rho_prev_star1)) > 1.e-6) ;

    // For star 2
    //------------

    do {

      // Solving the Poisson equation
      //------------------------------
      pot_star2 = - rho2.poisson() ;

      // Obtaining "lambda" in order to rescale the mapping
      //----------------------------------------------------
      double pot_star2_c = pot_star2(0,0,0,0) ;
      double pot_star2_s = pot_star2(nz-2,0,0,nr_star2-1) ;

      double lambda2_star2 = (ent_star2_s - ent_star2_c) /
	(pot_star2_s - pot_star2_c) ; 
      cout << "lambda2_star2 : " << lambda2_star2 << endl ; 

      // New potential after rescaling
      //-------------------------------
      pot_star2 = lambda2_star2 * pot_star2 ; 
      pot_star2_c *= lambda2_star2 ; 

      // First integral of motion
      // -------------------------
      ent_star2 = ent_star2_c + pot_star2 - pot_star2_c ;
      // enthalpy is just the same as the Lane-Emden function Theta

      // Rescale of the mapping
      //------------------------
      mp2.homothetie( sqrt(lambda2_star2) ) ; 

      // EOS : rho = rho(H) 
      //--------------------
      rho_prev_star2 = rho2 ;
      rho2 = pow(abs(ent_star2), n_index) ; 

      rho2.annule(nz-1) ;
      rho2.set_dzpuis(4) ;
      rho2.std_base_scal() ;

      cout << max(abs(rho2 - rho_prev_star2)) << endl;
      cout << "Maximum difference between rho2 and rho2_prev : "
	   << max(diffrel(rho2, rho_prev_star2)) << endl ;

    } while(max(diffrel(rho2, rho_prev_star2)) > 1.e-6) ;

    arrete() ;

    // Masses of stars
    //-----------------
    double mass1 = rho1.integrale() ;
    double mass2 = rho2.integrale() ;


    //-----------------------------------------------
    // Construction of the Darwin-like binary system
    //-----------------------------------------------

    int nr1_s = mg1.get_nr(nz-2) ;
    int nr2_s = mg2.get_nr(nz-2) ;

    Cmp pot11(mp1) ;
    Cmp pot12(mp2) ;
    Cmp pot22(mp2) ;
    Cmp pot21(mp1) ;
    pot11.std_base_scal() ;
    pot22.std_base_scal() ;
    pot12.std_base_scal() ;
    pot21.std_base_scal() ;

    Cmp rho1_prev(mp1) ;
    Cmp rho2_prev(mp2) ;
    Cmp rot1(mp1) ;
    Cmp rot2(mp2) ;

    Cmp ent1(mp1) ;
    Cmp ent2(mp2) ;
    double ent1_c = 1 ;
    double ent1_s = 0 ;
    double ent2_c = 1 ;
    double ent2_s = 0 ;

    Cmp dpot21(mp1) ;       // x derivative of pot21
    Cmp disp1(mp1) ;        // Cmp description of x 
    Cmp rho1_dpot21(mp1) ;  // multiply rho1 by dpot21
    Cmp rho1_disp1(mp1) ;   // multiply rho1 by disp1

    Cmp dpot12(mp2) ;
    Cmp disp2(mp2) ;
    Cmp rho2_dpot12(mp2) ;
    Cmp rho2_disp2(mp2) ;

    cout << "Mass of star 1 = " << mass1 << endl ;
    cout << "Mass of star 2 = " << mass2 << endl ;

    double mass1_prev ;
    double rel_mass1 ;
    double mass2_prev ;
    double rel_mass2 ;

    double omega2_1 ;
    double omega2_prev1 ;
    double rel_omega1 ;
    omega2_1 = (mass1 + mass2) / (M_PI * R_sep * R_sep * R_sep) ;

    double omega2_2 ;
    double omega2_prev2 ;
    double rel_omega2 ;
    omega2_2 = (mass1 + mass2) / (M_PI * R_sep * R_sep * R_sep) ;

    cout << "omega^2 = " << omega2_1 << endl ;

    double cc1 ;
    double cc1_prev ;
    double rel_cc1 ;
    cc1 = 0 ;

    double cc2 ;
    double cc2_prev ;
    double rel_cc2 ;
    cc2 = 0 ;

    double r1_x ;
    double r2_x ;

    arrete() ;

    do {

      // Solving the Poisson equation
      //------------------------------
      pot11 = - rho1.poisson() ;
      pot22 = - rho2.poisson() ;

      pot12.import(nz-1, pot11) ;
      pot21.import(nz-1, pot22) ;

      // Angular velocity
      //------------------

      // For star 1
      omega2_prev1 = omega2_1 ;
      dpot21 = pot21.dsdx() ;
      disp1 = x1 ;

      rho1_dpot21 = rho1 * dpot21 ;
      rho1_dpot21.annule(nz-1) ;
      rho1_dpot21.set_dzpuis(4) ;
      rho1_dpot21.std_base_scal() ;

      rho1_disp1 = rho1 * disp1 ;
      rho1_disp1.annule(nz-1) ;
      rho1_disp1.set_dzpuis(4) ;
      rho1_disp1.std_base_scal() ;

      double integ_r1_dp21 = rho1_dpot21.integrale() ;
      double integ_r1_dis1 = rho1_disp1.integrale() ;

      omega2_1 = 8. * integ_r1_dp21 / (mass1 * R_sep - 2. * integ_r1_dis1) ;

      cout << "Omega^2/(pi G rho_c) for star 1 : " << omega2_1 << endl ;

      rel_omega1 = fabs(1. - omega2_prev1 / omega2_1) ;
      cout << "Relative error in Omega^2 for star 1 : " << rel_omega1
	   << endl ;


      // For star 2
      omega2_prev2 = omega2_2 ;
      dpot12 = pot12.dsdx() ;
      disp2 = x2 ;

      rho2_dpot12 = rho2 * dpot12 ;
      rho2_dpot12.annule(nz-1) ;
      rho2_dpot12.set_dzpuis(4) ;
      rho2_dpot12.std_base_scal() ;

      rho2_disp2 = rho2 * disp2 ;
      rho2_disp2.annule(nz-1) ;
      rho2_disp2.set_dzpuis(4) ;
      rho2_disp2.std_base_scal() ;

      double integ_r2_dp12 = rho2_dpot12.integrale() ;
      double integ_r2_dis2 = rho2_disp2.integrale() ;

      omega2_2 = 8. * integ_r2_dp12 / (mass2 * R_sep - 2. * integ_r2_dis2) ;

      cout << "Omega^2/(pi G rho_c) for star 2 : " << omega2_2 << endl ;

      rel_omega2 = fabs(1. - omega2_prev2 / omega2_2) ;
      cout << "Relative error in Omega^2 for star 2 : " << rel_omega2
	   << endl ;

      //      arrete() ;

      // Obtaining "lambda" for each star in order to rescale mappings
      //---------------------------------------------------------------
      double pot11_c = pot11(0,0,0,0) ;
      double pot22_c = pot22(0,0,0,0) ;
      double pot21_c = pot21(0,0,0,0) ;
      double pot12_c = pot12(0,0,0,0) ;

      r1_x = mp1.val_r(nz-2,1.,M_PI/2,0) ;
      r2_x = mp2.val_r(nz-2,1.,M_PI/2,0) ;

      cout << "Radius x of star 1 : " << r1_x << endl
	   << "Radius x of star 2 : " << r2_x << endl ;

      double pot11_sx = pot11.val_point(r1_x,M_PI/2,0) ;
      double pot21_sx = pot21.val_point(r1_x,M_PI/2,0) ;
      double pot22_sx = pot22.val_point(r2_x,M_PI/2,0) ;
      double pot12_sx = pot12.val_point(r2_x,M_PI/2,0) ;

      // lambda2
      //---------
      double lambda2_1 = (ent1_c - ent1_s - pot21_c + pot21_sx +
			  0.125 * omega2_1 * r1_x * (r1_x - R_sep)) /
	(pot11_c - pot11_sx) ;

      double lambda2_2 = (ent2_c - ent2_s - pot12_c + pot12_sx +
			  0.125 * omega2_2 * r2_x * (r2_x - R_sep)) /
	(pot22_c - pot22_sx) ;

      cout << "lambda2 for star 1 : " << lambda2_1 << endl
	   << "lambda2 for star 2 : " << lambda2_2 << endl ;

      // New potentials after rescaling
      //--------------------------------
      pot11 = lambda2_1 * pot11 ;
      pot11_c = lambda2_1 * pot11_c ;
      pot11_sx = lambda2_1 * pot11_sx ;

      pot22 = lambda2_2 * pot22 ;
      pot22_c = lambda2_2 * pot22_c ;
      pot22_sx = lambda2_2 * pot22_sx ;

      // Integration constant
      //----------------------
      cc1_prev = cc1 ;
      cc1 = ent1_c - pot11_c - pot21_c - omega2_1 * R_sep * R_sep /32. ;
      cc2_prev = cc2 ;
      cc2 = ent2_c - pot22_c - pot12_c - omega2_2 * R_sep * R_sep /32. ;

      cout << "Integration constant for star 1 : " << cc1 << endl
	   << "Integration constant for star 2 : " << cc2 << endl ;

      rel_cc1 = fabs(1. - cc1_prev / cc1) ;
      rel_cc2 = fabs(1. - cc2_prev / cc2) ;
      cout << "Relative error in cc1 : " << rel_cc1 << endl
	   << "Relative error in cc2 : " << rel_cc2 << endl ;

      //      arrete() ;


      // First integration of motion
      //-----------------------------
      rot1 = 0.125 * omega2_1 * (0.25 * R_sep * R_sep - R_sep * x1 +
				 x1 * x1 + y1 * y1) ;
      rot1.annule(nz-1) ;

      rot2 = 0.125 * omega2_2 * (0.25 * R_sep * R_sep - R_sep * x2 +
				 x2 * x2 + y2 * y2) ;
      rot2.annule(nz-1) ;

      ent1 = pot11 + pot21 + rot1 + cc1 ;
      ent1.annule(nz-1) ;

      ent2 = pot22 + pot12 + rot2 + cc2 ;
      ent2.annule(nz-1) ;

      // Rescaling of mappings
      //-----------------------
      mp1.homothetie( sqrt(lambda2_1) ) ;
      mp2.homothetie( sqrt(lambda2_2) ) ;

      rho1_prev = rho1 ;
      rho2_prev = rho2 ;

      // EOS : rho = rho(H)
      //--------------------
      rho1 = eos1_local(ent1, nz-1, n_index) ;
      rho2 = eos2_local(ent2, nz-1, n_index) ;

      rho1.set_dzpuis(4) ;
      rho2.set_dzpuis(4) ;
      rho1.std_base_scal() ;
      rho2.std_base_scal() ;

      cout << max(abs(rho1 - rho1_prev)) << endl
	   << max(abs(rho2 - rho2_prev)) << endl ;

      // Masses of stars
      //-----------------

      // For star 1
      mass1_prev = mass1 ;
      mass1 = rho1.integrale() ;
      cout << "Mass/rho_c for star 1 : " << mass1 << endl ;

      rel_mass1 = fabs(1. - mass1_prev / mass1) ;
      cout << "Relative error in Mass for star 1 : " << rel_mass1 << endl ;

      // For star 2
      mass2_prev = mass2 ;
      mass2 = rho2.integrale() ;
      cout << "Mass/rho_c for star 2 : " << mass2 << endl ;

      rel_mass2 = fabs(1. - mass2_prev / mass2) ;
      cout << "Relative error in Mass for star 2 : " << rel_mass2 << endl ;

      cout << "Maximum difference between rho1 and rho1_prev : "
	   << max(diffrel(rho1, rho1_prev)) << endl
	   << "Maximum difference between rho2 and rho2_prev : "
	   << max(diffrel(rho2, rho2_prev)) << endl ;

      //      arrete() ;

    } while(max(diffrel(rho1, rho1_prev)) > 1.e-4) ;

    cout << "Coef of rho1 : " << endl ; 
    rho1.affiche_seuil(cout) ; 
    
    cout << "Coef of pot11 : " << endl ; 
    pot11.affiche_seuil(cout) ; 

    cout << "Value of pot11 at the origin : " << pot11(0,0,0,0) << endl ;
    cout << "Value of pot11 at the surface (r=ray) : "
         << pot11(nz-1,0,0,0) << endl ;
    cout << "Value of pot11 at the surface (r=ray) : "
         << pot11(nz-2,0,0,nr1_s-1) << endl ; 
    arrete() ; 

    cout << "Coef of rho2 : " << endl ; 
    rho2.affiche_seuil(cout) ; 
    
    cout << "Coef of pot22 : " << endl ; 
    pot22.affiche_seuil(cout) ; 

    cout << "Value of pot22 at the origin : " << pot22(0,0,0,0) << endl ;
    cout << "Value of pot22 at the surface (r=ray) : "
         << pot22(nz-1,0,0,0) << endl ;
    cout << "Value of pot22 at the surface (r=ray) : "
         << pot22(nz-2,0,0,nr2_s-1) << endl ; 
    arrete() ; 

    // Checking whether or not the potentials satisfy their Poisson equations
    //------------------------------------------------------------------------
    Cmp lap1 = - pot11.laplacien() ; 
    cout << "max( |lap(pot11) - rho1| ) " << max(abs(lap1 - rho1)) << endl ; 
    cout << "relative error for star 1 : " << diffrel(lap1, rho1) << endl ;

    Cmp lap2 = - pot22.laplacien() ;
    cout << "max( |lap(pot22) - rho2| ) " << max(abs(lap2 - rho2)) << endl ; 
    cout << "relative error for star 2: " << diffrel(lap2, rho2) << endl ;

    // Searching the radius to the x-axis direction for star 1
    //---------------------------------------------------------
    Param par_funct_zero_x ;
    par_funct_zero_x.add_cmp(ent1) ; 
    double phi_search_x = 0 ;
    double theta_search_x = M_PI/2 ;
    par_funct_zero_x.add_double(theta_search_x, 0) ;
    par_funct_zero_x.add_double(phi_search_x, 1) ;

    double r_min_search_x = 0 ;
    double r_max_search_x = r1_x ;
    double precis_x = 1e-10 ;
    int nitermax_x = 100 ;
    int niter_x ;

    double radius_x = zerosec(funct_zero_ent, par_funct_zero_x, 
			 r_min_search_x, r_max_search_x, precis_x,
			 nitermax_x, niter_x) ;

    cout << "Radius of star 1 in the x-axis direction : " << radius_x
	 << endl ;

    // Searching the radius opposite to the x-axis direction for star 1
    //------------------------------------------------------------------
    Param par_funct_zero_xopp ;
    par_funct_zero_xopp.add_cmp(ent1) ; 
    double phi_search_xopp = - M_PI ;
    double theta_search_xopp = M_PI/2 ;
    par_funct_zero_xopp.add_double(theta_search_xopp, 0) ;
    par_funct_zero_xopp.add_double(phi_search_xopp, 1) ;

    double r_min_search_xopp = 0 ;
    double r_max_search_xopp = r1_x ;
    double precis_xopp = 1e-10 ;
    int nitermax_xopp = 100 ;
    int niter_xopp ;

    double radius_xopp = zerosec(funct_zero_ent, par_funct_zero_xopp, 
			 r_min_search_xopp, r_max_search_xopp, precis_xopp,
			 nitermax_xopp, niter_xopp) ;

    cout << "Radius of star 1 opposite to the x-axis direction : "
	 << radius_xopp << endl ;

    // Searching the radius to the y-axis direction for star 1
    //---------------------------------------------------------
    Param par_funct_zero_y ;
    par_funct_zero_y.add_cmp(ent1) ; 
    double phi_search_y = M_PI/2 ;
    double theta_search_y = M_PI/2 ;
    par_funct_zero_y.add_double(theta_search_y, 0) ;
    par_funct_zero_y.add_double(phi_search_y, 1) ;

    double r_min_search_y = 0 ;
    double r_max_search_y = r1_x ;
    double precis_y = 1e-10 ;
    int nitermax_y = 100 ;
    int niter_y ;

    double radius_y = zerosec(funct_zero_ent, par_funct_zero_y, 
			 r_min_search_y, r_max_search_y, precis_y,
			 nitermax_y, niter_y) ;

    cout << "Radius of star 1 in the y-axis direction : " << radius_y
	 << endl ;

    // Searching the radius to the north pole for star 1
    //---------------------------------------------------
    Param par_funct_zero_z ;
    par_funct_zero_z.add_cmp(ent1) ; 
    double phi_search_z = 0 ;
    double theta_search_z = 0 ;
    par_funct_zero_z.add_double(theta_search_z, 0) ;
    par_funct_zero_z.add_double(phi_search_z, 1) ;

    double r_min_search_z = 0 ;
    double r_max_search_z = r1_x ;
    double precis_z = 1e-10 ;
    int nitermax_z = 100 ;
    int niter_z ;

    double radius_z = zerosec(funct_zero_ent, par_funct_zero_z, 
			 r_min_search_z, r_max_search_z, precis_z,
			 nitermax_z, niter_z) ;

    cout << "Radius of star 1 in the z-axis direction : " << radius_z
	 << endl ;

    // Radius of star 1
    //------------------
    double ray1 = mp1.val_r(nz-2,1.,M_PI/2,0) ;
    cout << "Radius of star 1 in the x-axis direction : " << ray1 << endl ;

    cout << "The ratios of the axes : " << endl
	 << "a_2/a_1 : " << radius_y/radius_x << endl
	 << "a_3/a_1 : " << radius_z/radius_x << endl ;

    // Radius of star 1
    //------------------
    double ray2 = mp2.val_r(nz-2,1.,M_PI/2,0) ;
    cout << "Radius of star 2 in the x-axis direction : " << ray2 << endl ;

    // Angular velocity
    //------------------
    double omega1 = sqrt(omega2_1) ;
    double omega2 = sqrt(omega2_2) ;

    cout << "Angular velocity of the star Omega/(pi G rho_c)^{1/2} : "
	 << omega1 << endl ;
    cout << "Angular velocity Omega^2/(pi G rho_c) : " << omega2_1 << endl ;

    arrete() ;

    // Total energy of the binary system
    //-----------------------------------

    // Internal energy
    //-----------------

    // For star 1
    Cmp rho1_n(mp1) ;

    rho1_n = pow(abs(rho1), 1 + 1/n_index) ;
    rho1_n.annule(nz-1) ;
    rho1_n.set_dzpuis(4) ;
    rho1_n.std_base_scal() ;

    double energy1_internal = 4. * M_PI *
      (1 / (1 + n_index)) * rho1_n.integrale() * ray1 / (mass1 * mass1) ;

    cout << "Internal energy of star 1 : "
	 << n_index * energy1_internal << endl ;

    // For star 2
    Cmp rho2_n(mp2) ;

    rho2_n = pow(abs(rho2), 1 + 1/n_index) ;
    rho2_n.annule(nz-1) ;
    rho2_n.set_dzpuis(4) ;
    rho2_n.std_base_scal() ;

    double energy2_internal = 4. * M_PI *
      (1 / (1 + n_index)) * rho2_n.integrale() * ray2 / (mass2 * mass2) ;

    cout << "Internal energy of star 2 : "
	 << n_index * energy2_internal << endl ;

    cout << "Internal energy of the Darwin binary system : E_int/(GM^2/ray) : "
	 << n_index * energy1_internal + n_index * energy2_internal << endl ;

    // Self-gravity energy
    //---------------------

    // For star 1
    Cmp rho1_pot11(mp1) ;
    rho1_pot11 = rho1 * pot11 ;
    rho1_pot11.annule(nz-1) ;
    rho1_pot11.set_dzpuis(4) ;
    rho1_pot11.std_base_scal() ;
    double energy1_selfgrav = -2. * M_PI * rho1_pot11.integrale() * ray1 /
      (mass1 * mass1) ;

    cout << "Self-gravity energy of star 1 : " << energy1_selfgrav << endl ;

    // For star 2
    Cmp rho2_pot22(mp2) ;
    rho2_pot22 = rho2 * pot22 ;
    rho2_pot22.annule(nz-1) ;
    rho2_pot22.set_dzpuis(4) ;
    rho2_pot22.std_base_scal() ;
    double energy2_selfgrav = -2. * M_PI * rho2_pot22.integrale() * ray2 /
      (mass2 * mass2) ;

    cout << "Self-gravity energy of star 2 : " << energy2_selfgrav << endl ;

    cout <<
      "Self-gravity energy of the Darwin binary system : E_self/(GM^2/ray) : "
	 << energy1_selfgrav + energy2_selfgrav << endl ;

    // Interaction energy
    //--------------------

    // For star 1
    Cmp rho1_pot21(mp1) ;
    rho1_pot21 = rho1 * pot21 ;
    rho1_pot21.annule(nz-1) ;
    rho1_pot21.set_dzpuis(4) ;
    rho1_pot21.std_base_scal() ;
    double energy1_interact = -2. * M_PI * rho1_pot21.integrale() * ray1 /
      (mass1 * mass1) ;

    cout << "Interaction energy for star 1 : " << energy1_interact << endl ;

    // For star 2
    Cmp rho2_pot12(mp2) ;
    rho2_pot12 = rho2 * pot12 ;
    rho2_pot12.annule(nz-1) ;
    rho2_pot12.set_dzpuis(4) ;
    rho2_pot12.std_base_scal() ;
    double energy2_interact = -2. * M_PI * rho2_pot12.integrale() * ray2 /
      (mass2 * mass2) ;

    cout << "Interaction energy for star 2 : " << energy2_interact << endl ;

    cout <<
      "Interaction energy of the Darwin binary system : E_ext/(GM^2/ray) : "
	 << energy1_interact + energy2_interact << endl ;

    // Kinetic energy
    //----------------

    // For star 1
    Cmp rr1(mp1) ;
    Cmp rho1_quad1(mp1) ;
    rr1 = x1 * x1 + y1 * y1 - R_sep * x1 ;
    rho1_quad1 = rho1 * rr1 ;
    rho1_quad1.annule(nz-1) ;
    rho1_quad1.set_dzpuis(4) ;
    rho1_quad1.std_base_scal() ;
    double energy1_kinetic = 0.5 * M_PI * omega2_1 *
      (0.25 * R_sep * R_sep + rho1_quad1.integrale() / mass1 ) *
      ray1 / mass1 ;

    cout << "Kinetic energy for star 1 : " << energy1_kinetic << endl ;

    // For star 2
    Cmp rr2(mp2) ;
    Cmp rho2_quad2(mp2) ;
    rr2 = x2 * x2 + y2 * y2 - R_sep * x2 ;
    rho2_quad2 = rho2 * rr2 ;
    rho2_quad2.annule(nz-1) ;
    rho2_quad2.set_dzpuis(4) ;
    rho2_quad2.std_base_scal() ;
    double energy2_kinetic = 0.5 * M_PI * omega2_2 *
      (0.25 * R_sep * R_sep + rho2_quad2.integrale() / mass2 ) *
      ray2 / mass2 ;

    cout << "Kinetic energy for star 2 : " << energy2_kinetic << endl ;

    cout << "Kinetic energy of the Darwin binary system : "
	 << "E_kinet/(GM^2/ray) : "
	 << energy1_kinetic + energy2_kinetic << endl ;

    // Total energy
    //--------------

    // For star 1
    double energy1_total = n_index * energy1_internal + energy1_selfgrav
      + energy1_kinetic + energy1_interact ;

    cout << "Total energy for star 1 : E_1/(GM_1^2/ray1) : "
	 << energy1_total << endl ;

    // For star 2
    double energy2_total = n_index * energy2_internal + energy2_selfgrav
      + energy2_kinetic + energy2_interact ;

    cout << "Total energy for star 2 : E_2/(GM_2^2/ray) : "
	 << energy2_total << endl ;

    // For binary system
    cout << "Total energy of the Darwin binary system : E/(GM^2/ray) : "
	 << energy1_total + energy2_total << endl ;

    // Angular momentum
    //------------------

    // For star 1
    double angmom1 = sqrt(M_PI) * omega1 *
      ( 0.25 * mass1 * R_sep * R_sep + rho1_quad1.integrale() ) /
      sqrt( mass1 * mass1 * mass1 * ray1 ) ;

    cout << "Angular momentum for star 1 : " << angmom1 << endl ;

    // For star 2
    double angmom2 = sqrt(M_PI) * omega2 *
      ( 0.25 * mass2 * R_sep * R_sep + rho2_quad2.integrale() ) /
      sqrt( mass2 * mass2 * mass2 * ray2 ) ;

    cout << "Angular momentum for star 2 : " << angmom2 << endl ;

    cout <<
      "Angular momentum of the Darwin binary system : J/(GM^3 ray)^{1/2} : "
	 << angmom1 + angmom2 << endl ;

    // Virial relation
    //-----------------

    // For star 1
    double virial1 = 3. * energy1_internal + energy1_selfgrav
      + 2. * energy1_kinetic + energy1_interact ;

    cout << "Virial relation for star 1 : " << virial1 << endl ;

    // For star 2
    double virial2 = 3. * energy2_internal + energy2_selfgrav
      + 2. * energy2_kinetic + energy2_interact ;

    cout << "Virial relation for star 2 : " << virial2 << endl ;

    cout << "Virial relation of the Darwin binary system : "
	 << virial1 + virial2 << endl ;

    arrete() ;

    //------------------------------------------
    //          Plot the figures
    //------------------------------------------
    double rmax = R_sep/2 ;
    //    cout << "r_max ?" << endl ;
    //    cin >> rmax ;

    des_profile(rho1,0.,rmax,M_PI/2,0.,"rho1/rho_c","rho1 (x direction)") ;
    des_profile(rho1,0.,rmax,M_PI/2,M_PI/2,"rho1/rho_c","rho1 (y direction)") ;
    des_profile(rho1,0.,rmax,0.,0.,"rho1/rho_c","rho1 (z direction)") ;

    //    des_coupe_x(rho1, 0., -rmax, rmax, -rmax, rmax, "rho1 (y-z plane)") ;
    des_coupe_x(rho1, 0., -rmax, rmax, -rmax, rmax, "rho1 (y-z plane)",&ent1) ;
    des_coupe_y(rho1, 0., -rmax, rmax, -rmax, rmax, "rho1 (z-x plane)",&ent1) ;
    des_coupe_z(rho1, 0., -rmax, rmax, -rmax, rmax, "rho1 (x-y plane)",&ent1) ;

    des_coef_xi(rho1.va, 0, 0, 0, 1.e-14, "log|c_i|", "domain no. 0") ;
    des_coef_theta(rho1.va, 0, 0, 0) ;
    des_coef_phi(rho1.va, 0, 0, 0) ;

    des_coef_xi(rho1.va, 1, 0, 0, 1.e-14, "log|c_i|", "domain no. 1") ;
    des_coef_theta(rho1.va, 0, 0, 0) ;
    des_coef_phi(rho1.va, 0, 0, 0) ;

    exit(-1) ; 
        
}

