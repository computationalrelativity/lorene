/*
 * Test program for solving Roche-like binaries
 * Fixing the orbital separation and the central value of the star
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

char test_roche_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:31  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.1  1999/12/24  11:08:16  keisuke
 * Initial revision
 *
 *
 * $Header$
 *
 */

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
Cmp eos_local(const Cmp& ent, int nzet, double n_index) ;
double funct_zero_ent(double r, const Param& par) ;

//**********************************************************************

void main(){
    
    // Identification of all the subroutines called by the code : 
    
    system("ident test_roche > identif_roche.d") ; 

    //-----------------------------------------------------------------------
    //			Input from file "part.d"
    //-----------------------------------------------------------------------
    
    int type_t, type_p, nt, np, nz, l;
    char blabla[80];
    double n_index, R_sep, p_ratio ;

    //----------------------------------------------------
    //             Declear the polytropic index
    //----------------------------------------------------
    cout << "Input the polytropic index : " ;
    cin >> n_index ;
    cout << "Input the orbital separation : " ;
    cin >> R_sep ;
    cout << "Input the ratio M_1/M_2 : " ;
    cin >> p_ratio ;

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
    //		Construction of a multi-grid
    //-----------------------------------------------------------------------
    
    type_t = SYM ; 
    type_p = NONSYM ; 
    
    Mg3d mg(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    cout << endl << "Grid mg : " << mg << endl ; 
    
    //-----------------------------------------------------------------------
    //		Construction of a mapping
    //-----------------------------------------------------------------------
    
    Map_af mp(mg, bornes) ;
    cout << "Mapping mp : " << mp << endl ; 

    //-----------------------------------------------------------------------
    //		Construction of a Cmp
    //-----------------------------------------------------------------------

    //----------------------
    // Spherical Star
    //----------------------
    Cmp rho(mp) ; 

    // Initial density profile
    //-------------------------
    rho = 1  ;		        // rho = 1 in internal domains
    rho.annule(nz-1) ;		// rho = 0 in external domain

    rho.set_dzpuis(4) ; 
    rho.std_base_scal() ; // Sets the standard basis for spectral expansions

    Cmp pot_star(mp) ;
    Cmp rho_prev_star(mp) ;

    Cmp ent_star(mp) ; 
    double ent_star_c = 1 ; 	// central value of enthalpy
    double ent_star_s = 0 ; 	// surface value of enthalpy

    int nr_star = mg.get_nr(nz-2) ; 

    do {

      // Solving the Poisson equation
      //------------------------------
      pot_star = - rho.poisson() ;

      // Obtaining "lambda" in order to rescale the mapping
      //----------------------------------------------------
      double pot_star_c = pot_star(0,0,0,0) ;
      double pot_star_s = pot_star(nz-2,0,0,nr_star-1) ;

      double lambda2_star = (ent_star_s - ent_star_c) /
	(pot_star_s - pot_star_c) ; 
      cout << "lambda2_star : " << lambda2_star << endl ; 

      // New potential after rescaling
      //-------------------------------
      pot_star = lambda2_star * pot_star ; 
      pot_star_c *= lambda2_star ; 

      // First integral of motion
      // -------------------------
      ent_star = ent_star_c + pot_star - pot_star_c ;
      // enthalpy is just the same as the Lane-Emden function Theta

      // Rescale of the mapping
      //------------------------
      mp.homothetie( sqrt(lambda2_star) ) ; 

      // EOS : rho = rho(H) 
      //--------------------
      rho_prev_star = rho ;
      rho = pow(abs(ent_star), n_index) ; 

      rho.annule(nz-1) ;
      rho.set_dzpuis(4) ;
      rho.std_base_scal() ;

      cout << max(abs(rho - rho_prev_star)) << endl;
      cout << "Maximum difference between rho and rho_prev : "
	   << max(diffrel(rho, rho_prev_star)) << endl ;

    } while(max(diffrel(rho, rho_prev_star)) > 1.e-6) ;

    arrete() ;

    // Mass of the star
    //------------------
    double mass = rho.integrale() ; 

    //---------------
    // Binary
    //---------------
    const Coord& r = mp.r ; 
    const Coord& x = mp.x ; 
    const Coord& y = mp.y ;
    const Coord& z = mp.z ;
    const Coord& sint = mp.sint ; 

    Cmp pot(mp) ;
    Cmp pot2(mp) ;
    Cmp pot3(mp) ;
    Cmp disp(mp) ;
    Cmp rho_pot3(mp) ;
    Cmp rho_disp(mp) ;
    Cmp rho_prev(mp) ;
    Cmp rot(mp) ;

    Cmp ent(mp) ;
    double ent_c = 1 ;
    double ent_s = 0 ;

    double cc ;
    double cc_prev ;
    double rel_cc ;
    cc = 0 ;

    cout << "mass = " << mass << endl ;
    double mass_prev ;
    double rel_mass ;

    double omega2 ;
    double omega2_prev ;
    double rel_omega ;
    omega2 = mass * (1 + p_ratio) / (M_PI * p_ratio * R_sep * R_sep * R_sep) ;
    cout << "omega^2 = " << omega2 << endl ;

    double r_x ;
    int nr_s = mg.get_nr(nz-2) ; 

      arrete() ;

    do {

      // Solving the Poisson equation
      //------------------------------
      pot = - rho.poisson() ;

      // Angular velocity
      //------------------
      omega2_prev = omega2 ;
      pot3 = (R_sep - x) / pow((R_sep-x) * (R_sep-x) + y*y + z*z, 1.5) ;
      disp = x ;

      rho_pot3 = rho * pot3 ;
      rho_pot3.annule(nz-1) ;
      rho_pot3.set_dzpuis(4) ;
      rho_pot3.std_base_scal() ;

      rho_disp = rho * disp ;
      rho_disp.annule(nz-1) ;
      rho_disp.set_dzpuis(4) ;
      rho_disp.std_base_scal() ;

      double integ1 = rho_pot3.integrale() ;
      double integ2 = rho_disp.integrale() ;

      //      cout << "Integrales : " << endl
      //	   << integ1 << endl << integ2 << endl ;

      omega2 = (1/(p_ratio * M_PI)) * integ1 /
	(R_sep/(1+p_ratio) - integ2 / mass) ;

      cout << "Omega^2/(pi G rho_c) : " << omega2 << endl ;

      rel_omega = fabs(1. - omega2_prev / omega2) ;
      cout << "Relative error in Omega^2 : " << rel_omega << endl ;

      // Obtaining "lambda" in order to rescale the mapping
      //----------------------------------------------------
      double pot_c = pot(0,0,0,0) ;
      r_x = mp.val_r(nz-2,1.,M_PI/2,0) ;
      cout << "Radius x : " << r_x << endl ;
      double pot_sx = pot.val_point(r_x,M_PI/2,0) ;
      cout << pot_c << endl << pot_sx << endl ;

      double lambda2 = (ent_c - ent_s +
			0.25 * mass *
			(1/(R_sep - r_x) - 1/R_sep) / (p_ratio * M_PI) +
			0.125 * omega2 *
			(-2. * R_sep * r_x / (1+p_ratio) + r_x * r_x)) /
	(pot_c - pot_sx) ;

      cout << "lambda2 : " << lambda2 << endl ;

      // Integration constant
      //----------------------
      cc_prev = cc ;
      cc = ent_c - pot_c - 0.25 * mass / (p_ratio * M_PI * R_sep)
	-0.125 * omega2 * R_sep * R_sep / ((1+p_ratio) * (1+p_ratio)) ;
      cout << "Integration constant : " << cc << endl ;

      rel_cc = fabs(1. - cc_prev / cc) ;
      cout << "Relative error in cc : " << rel_cc << endl ;

      // New potential after rescaling
      //-------------------------------
      pot = lambda2 * pot ;
      pot_c = lambda2 * pot_c ;
      pot_sx = lambda2 * pot_sx ;

      // First integration of motion
      //-----------------------------
      rot = 0.125 * omega2 * (- 2*R_sep*x/(1+p_ratio) + x*x + y*y) ;
      rot.annule(nz-1) ;
      pot2 = 1./pow( (R_sep - x) * (R_sep - x) + y*y + z*z, 0.5) - 1/R_sep ;
      //      pot2.set_dzpuis(4) ;
      //      pot2.std_base_scal() ;

      ent = ent_c - pot_c + pot +
	0.25 * mass * pot2 / (p_ratio * M_PI) + rot ;
      ent.annule(nz-1) ;
      //      ent.set_dzpuis(4) ;
      //      ent.std_base_scal() ;


      // Rescaling of the map
      //----------------------
      mp.homothetie( sqrt(lambda2) ) ;

      rho_prev = rho ;

      // EOS : rho = rho(H)
      //--------------------
      rho = eos_local(ent, nz-1, n_index) ;

      rho.set_dzpuis(4) ;
      rho.std_base_scal() ;

      cout << max(abs(rho - rho_prev)) << endl ;

      // Mass of the star
      //------------------
      mass_prev = mass ;
      mass = rho.integrale() ;
      cout << "Mass/rho_c : " << mass << endl ;

      rel_mass = fabs(1. - mass_prev/mass) ;
      cout << "Relative error in Mass : " << rel_mass << endl ;

      cout << "Maximum difference between rho and rho_prev : "
	   << max(diffrel(rho, rho_prev)) << endl ;

    } while(max(diffrel(rho, rho_prev)) > 1.e-4) ;

    cout << "Coef of rho : " << endl ; 
    rho.affiche_seuil(cout) ; 
    
    cout << "Coef of pot : " << endl ; 
    pot.affiche_seuil(cout) ; 

    arrete() ;
    cout << "Value of pot at the origin : " << pot(0,0,0,0) << endl ;
    cout << "Value of pot at the surface (r=ray) : "
         << pot(nz-1,0,0,0) << endl ;
    cout << "Value of pot at the surface (r=ray) : "
         << pot(nz-2,0,0,nr_s-1) << endl ; 
    //    arrete() ; 

    // Checking whether or not the potential satisfies the Poisson equation
    //----------------------------------------------------------------------
    Cmp lap = - pot.laplacien() ; 

    cout << "max( |lap(pot) - rho| ) " << max(abs(lap - rho)) << endl ; 
    cout << "relative error : " << diffrel(lap, rho) << endl ;

    // Searching the radius to the x-axis direction
    //----------------------------------------------
    Param par_funct_zero_x ;
    par_funct_zero_x.add_cmp(ent) ; 
    double phi_search_x = 0 ;
    double theta_search_x = M_PI/2 ;
    par_funct_zero_x.add_double(theta_search_x, 0) ;
    par_funct_zero_x.add_double(phi_search_x, 1) ;

    double r_min_search_x = 0 ;
    double r_max_search_x = r_x ;
    double precis_x = 1e-10 ;
    int nitermax_x = 100 ;
    int niter_x ;

    double radius_x = zerosec(funct_zero_ent, par_funct_zero_x, 
			 r_min_search_x, r_max_search_x, precis_x,
			 nitermax_x, niter_x) ;

    cout << "Radius of the star in the x-axis direction : " << radius_x
	 << endl ;

    // Searching the radius opposite to the x-axis direction
    //-------------------------------------------------------
    Param par_funct_zero_xopp ;
    par_funct_zero_xopp.add_cmp(ent) ; 
    double phi_search_xopp = M_PI ;
    double theta_search_xopp = M_PI/2 ;
    par_funct_zero_xopp.add_double(theta_search_xopp, 0) ;
    par_funct_zero_xopp.add_double(phi_search_xopp, 1) ;

    double r_min_search_xopp = 0 ;
    double r_max_search_xopp = r_x ;
    double precis_xopp = 1e-10 ;
    int nitermax_xopp = 100 ;
    int niter_xopp ;

    double radius_xopp = zerosec(funct_zero_ent, par_funct_zero_xopp, 
			 r_min_search_xopp, r_max_search_xopp, precis_xopp,
			 nitermax_xopp, niter_xopp) ;

    cout << "Radius of the star opposite to the x-axis direction : "
	 << radius_xopp << endl ;

    // Searching the radius to the y-axis direction
    //----------------------------------------------
    Param par_funct_zero_y ;
    par_funct_zero_y.add_cmp(ent) ; 
    double phi_search_y = M_PI/2 ;
    double theta_search_y = M_PI/2 ;
    par_funct_zero_y.add_double(theta_search_y, 0) ;
    par_funct_zero_y.add_double(phi_search_y, 1) ;

    double r_min_search_y = 0 ;
    double r_max_search_y = r_x ;
    double precis_y = 1e-10 ;
    int nitermax_y = 100 ;
    int niter_y ;

    double radius_y = zerosec(funct_zero_ent, par_funct_zero_y, 
			 r_min_search_y, r_max_search_y, precis_y,
			 nitermax_y, niter_y) ;

    cout << "Radius of the star in the y-axis direction : " << radius_y
	 << endl ;

    // Searching the radius to the north pole
    //----------------------------------------
    Param par_funct_zero_z ;
    par_funct_zero_z.add_cmp(ent) ; 
    double phi_search_z = 0 ;
    double theta_search_z = 0 ;
    par_funct_zero_z.add_double(theta_search_z, 0) ;
    par_funct_zero_z.add_double(phi_search_z, 1) ;

    double r_min_search_z = 0 ;
    double r_max_search_z = r_x ;
    double precis_z = 1e-10 ;
    int nitermax_z = 100 ;
    int niter_z ;

    double radius_z = zerosec(funct_zero_ent, par_funct_zero_z, 
			 r_min_search_z, r_max_search_z, precis_z,
			 nitermax_z, niter_z) ;

    cout << "Radius of the star in the z-axis direction : " << radius_z
	 << endl ;

    // Radius of the star
    //--------------------
    double ray = mp.val_r(nz-2,1.,M_PI/2,0) ;
    cout << "Radius of the star in the x-axis direction : " << ray << endl ;

    cout << "The ratios of the axes : " << endl
	 << "a_2/a_1 : " << radius_y/radius_x << endl
	 << "a_3/a_1 : " << radius_z/radius_x << endl ;

    // Angular velocity
    //------------------
    double omega = sqrt(omega2) ;

    cout << "Angular velocity of the star Omega/(pi G rho_c)^{1/2} : "
	 << omega << endl ;
    cout << "Angular velocity Omega^2/(pi G rho_c) : " << omega2 << endl ;

    arrete() ;

    // Total energy of the binary system
    //-----------------------------------

    // Internal energy
    Cmp rho_n(mp) ;

    rho_n = pow(abs(rho), 1 + 1/n_index) ;
    rho_n.annule(nz-1) ;
    rho_n.set_dzpuis(4) ;
    rho_n.std_base_scal() ;

    double energy_internal = 4. * M_PI *
      (1 / (1 + n_index)) * rho_n.integrale() * ray / (mass * mass) ;
    cout << "Internal energy of the Roche binary system : E_int/(GM^2/ray) : "
	 << n_index * energy_internal << endl ;

    // Self-gravity energy
    Cmp rho_pot(mp) ;
    rho_pot = rho * pot ;
    rho_pot.annule(nz-1) ;
    rho_pot.set_dzpuis(4) ;
    rho_pot.std_base_scal() ;
    double energy_selfgrav = -2. * M_PI * rho_pot.integrale() * ray /
      (mass * mass) ;
    cout <<
      "Self-gravity energy of the Roche binary system : E_self/(GM^2/ray) : "
	 << energy_selfgrav << endl ;

    // Interaction energy
    Cmp pot21(mp) ;
    Cmp rho_pot21(mp) ;
    pot21 = 1/pow((R_sep-x) * (R_sep-x) + y*y + z*z, 0.5) ;
    rho_pot21 = rho * pot21 ;
    rho_pot21.annule(nz-1) ;
    rho_pot21.set_dzpuis(4) ;
    rho_pot21.std_base_scal() ;
    double energy_interact = -0.5 * rho_pot21.integrale() * ray /
      (mass * p_ratio) ;
    cout <<
      "Interaction energy of the Roche binary system : E_ext/(GM^2/ray) : "
	 << energy_interact << endl ;

    // Kinetic energy
    Cmp rr(mp) ;
    Cmp rho_quad(mp) ;
    rr = r * r * sint * sint - 2.*R_sep*x/(1+p_ratio) ;
    rho_quad = rho * rr ;
    rho_quad.annule(nz-1) ;
    rho_quad.set_dzpuis(4) ;
    rho_quad.std_base_scal() ;
    double energy_kinetic = 0.5 * M_PI * omega2 *
      ( pow(R_sep/(1+p_ratio), 2) + rho_quad.integrale() / mass ) * ray/mass ;
    cout << "Kinetic energy of the Roche binary system : E_kinet/(GM^2/ray) : "
	 << energy_kinetic << endl ;

    double energy_total = n_index * energy_internal + energy_selfgrav
      + energy_kinetic + energy_interact ;

    cout << "Total energy of the Roche binary system : E/(GM^2/ray) : "
	 << energy_total << endl ;

    // Angular momentum
    //------------------
    double angmom = sqrt(M_PI) * omega *
      ( mass * pow(R_sep/(1+p_ratio), 2) + rho_quad.integrale() ) /
      sqrt( mass * mass * mass * ray ) ;

    cout <<
      "Angular momentum of the Roche binary system : J/(GM^3 ray)^{1/2} : "
	 << angmom << endl ;

    // Virial relation
    //-----------------
    double virial = 3. * energy_internal + energy_selfgrav
      + 2. * energy_kinetic + energy_interact ;
    cout << "Virial relation : " << virial << endl ;

    //------------------------------------------
    //          Plot the figures
    //------------------------------------------
    double rmax ;
    cout << "r_max ?" << endl ;
    cin >> rmax ;

    des_profile(rho, 0., rmax, M_PI/2, 0., "rho/rho_c","rho (x direction)") ;
    des_profile(rho, 0.,rmax,M_PI/2,M_PI/2,"rho/rho_c","rho (y direction)") ;
    des_profile(rho, 0., rmax, 0., 0., "rho/rho_c", "rho (z direction)") ;

    des_profile(ent,0.,rmax,M_PI/2,0.,"enthalpy","enthalpy (x direction)") ;
    des_profile(ent,0.,rmax,M_PI/2,M_PI/2,"enthalpy",
		"enthalpy (y direction)") ;
    des_profile(ent,0.,rmax,0.,0.,"enthalpy","enthalpy (z direction)") ;

    /*
    des_profile(rho_n,0.,rmax,M_PI/2,0.,"rho_n","rho_n (x direction)") ;
    des_profile(rho_n,0.,rmax,M_PI/2,M_PI/2,"rho_n","rho_n (y direction)") ;
    des_profile(rho_n,0.,rmax,0.,0.,"rho_n","rho_n (z direction)") ;

    des_profile(rho_pot,0.,rmax,M_PI/2,0.,"rho_pot","rho_pot (x direction)") ;
    des_profile(rho_pot,0.,rmax,M_PI/2,M_PI/2,"rho_pot",
		"rho_pot (y direction)") ;
    des_profile(rho_pot,0.,rmax,0.,0.,"rho_pot","rho_pot (z direction)") ;

    des_profile(rho_pot21,0.,rmax,M_PI/2,0.,"rho_pot21",
		"rho_pot21 (x direction)") ;
    des_profile(rho_pot21,0.,rmax,M_PI/2,M_PI/2,"rho_pot21",
		"rho_pot21 (y direction)") ;
    des_profile(rho_pot21,0.,rmax,0.,0.,"rho_pot21",
		"rho_pot21 (z direction)") ;

    des_profile(rho_quad,0.,rmax,M_PI/2,0.,"rho_quad",
		"rho_quad (x direction)") ;
    des_profile(rho_quad,0.,rmax,M_PI/2,M_PI/2,"rho_quad",
		"rho_quad (y direction)") ;
    des_profile(rho_quad,0.,rmax,0.,0.,"rho_quad","rho_quad (z direction)") ;
    */

    des_coupe_x(rho, 0., -rmax, rmax, -rmax, rmax, "rho (y-z plane)", &ent) ;
    des_coupe_y(rho, 0., -rmax, rmax, -rmax, rmax, "rho (z-x plane)", &ent) ;
    des_coupe_z(rho, 0., -rmax, rmax, -rmax, rmax, "rho (x-y plane)", &ent) ;

    des_coef_xi(rho.va, 0, 0, 0, 1.e-14, "log|c_i|", "domain no. 0") ;
    des_coef_theta(rho.va, 0, 0, 0) ;
    des_coef_phi(rho.va, 0, 0, 0) ;

    des_coef_xi(rho.va, 1, 0, 0, 1.e-14, "log|c_i|", "domain no. 1") ;
    des_coef_theta(rho.va, 0, 0, 0) ;
    des_coef_phi(rho.va, 0, 0, 0) ;

    exit(-1) ; 
        
}

