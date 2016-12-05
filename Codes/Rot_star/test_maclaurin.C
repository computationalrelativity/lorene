/*
 * Test program for solving Maclaurin-like figures
 * One zone for the internal of the star
 *
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:18:26  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:59  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:09:46  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/01/09 11:07:51  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  1999/12/23  13:39:23  keisuke
 * Non significant change (for test only).
 *
 * Revision 1.1  1999/12/23  13:35:01  keisuke
 * Initial revision
 *
 *
 * $Header$
 *
 */

// headers C
#include <cstdlib>
#include <cmath>

// headers Lorene
#include "type_parite.h"
#include "cmp.h"
#include "graphique.h"
#include "utilitaires.h"

namespace Lorene {
// Local prototypes:
Cmp eos_local(const Cmp& ent, int nzet, double n_index) ;
}
//**********************************************************************

int main(){
    
    // Identification of all the subroutines called by the code : 
    
    system("ident test_maclaurin > identif_maclaurin.d") ; 

    //-----------------------------------------------------------------------
    //			Input from file "part.d"
    //-----------------------------------------------------------------------
    
    int type_t, type_p, nt, np, nz, l;
    char blabla[80];
    //    const double criteria = 1.e-8 ;
    double n_index, radius ;

    //----------------------------------------------------
    //             Declear the polytropic index
    //----------------------------------------------------
    cout << "Input a polytropic index: " ;
    cin >> n_index ;
    cout << "Input an equatorial radius of the star : " ;
    cin >> radius ;

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

    Cmp rho(mp) ; 

    const Coord& r = mp.r ; 
    //    const Coord& x = mp.x ; 
    //    const Coord& y = mp.y ; 
    const Coord& sint = mp.sint ; 
    
    rho = 1  ;
    rho.annule(nz-1) ;		// rho = 1 in internal domains
				// rho = 0 in external domain

    rho.set_dzpuis(4) ; 
    rho.std_base_scal() ; // Sets the standard basis for spectral expansions

    int nr_s = mg.get_nr(nz-2) ;

    Cmp pot(mp) ;
    Cmp rho_prev(mp) ;
    Cmp rot(mp) ;

    Cmp ent(mp) ;
    double ent_c = 1 ;  // central value of enthalpy
    double ent_s = 0 ;  // surface value of enthalpy

    double omega2 ;
    double omega2_prev ;
    double rel_omega ;
    omega2 = 0 ;

    do {

      // Solving the Poisson equation
      //------------------------------
      pot = - rho.poisson() ;

      // Obtaining "lambda" in order to rescale the mapping
      //----------------------------------------------------
      double pot_c = pot(0,0,0,0) ; // central value of potential

      // surface value of potential in the equatorial plane
      double pot_seq = pot(nz-2,0,nt-1,nr_s-1) ;

      double radi = mp.val_r(nz-2,1.,M_PI/2,0) ;

      double lambda = radius / radi ;
      cout << "lambda : " << lambda << endl ;

      // Angular velocity
      //------------------
      omega2_prev = omega2 ;
      omega2 = 8.*(ent_s - ent_c - pot_seq + pot_c) / (radius * radius) ;
      cout << "Omega^2/(pi G rho_c) : " << omega2 << endl ;

      rel_omega = fabs(1. - omega2_prev / omega2) ;
      cout << "Relative error in Omega^2 : " << rel_omega << endl ;

      // New potential after rescaling
      //-------------------------------
      pot = lambda * lambda * pot ;
      pot_c = lambda * lambda * pot_c ;
      pot_seq = lambda * lambda * pot_seq ;

      // First integral of motion
      //--------------------------
      rot = 0.125 * omega2 * r * r * sint * sint ;
      //      rot.annule(nz-1) ;
      ent = ent_c + pot - pot_c + rot ; // enthalpy
      ent.annule(nz-1) ;

      // Rescale of the map
      //--------------------
      mp.homothetie( lambda ) ; 

      rho_prev = rho ;

      // EOS : rho = rho(H)
      //--------------------
      rho = eos_local(ent, nz-1, n_index) ; 

      rho.set_dzpuis(4) ;
      rho.std_base_scal() ;

      cout << max(abs(rho - rho_prev)) << endl ;
      cout << "Maximum error in rho : "
	   << max(diffrel(rho, rho_prev)) << endl ;

    } while(max(diffrel(rho, rho_prev)) > 1.e-10) ;

    arrete() ;

    cout << "Coef of rho : " << endl ; 
    rho.affiche_seuil(cout) ;
    
    cout << "Coef of pot : " << endl ; 
    pot.affiche_seuil(cout) ; 

    arrete() ;
    cout << "Value of pot at the origin : " << pot(0,0,0,0) << endl ;
    cout << "Value of pot at the surface (r=radius_eq) : "
	 << pot(nz-1,0,0,0) << endl ;
    cout << "Value of pot at the surface (r=radius_eq) : "
	 << pot(nz-2,0,0,nr_s-1) << endl ;


    // Checking whether or not the potential satisfies the Poisson equation
    //----------------------------------------------------------------------
    Cmp lap = - pot.laplacien() ; 
    
    cout << "max( |lap(pot) - rho| ) " << max(abs(lap - rho)) << endl ; 
    cout << "relative error : " << diffrel(lap, rho) << endl ;


    // Serching the radius to the north pole
    //---------------------------------------
    double zaxis = 0 ;
    do {
      zaxis += 1.e-1 ;
    } while(ent.val_point(zaxis,0,0) >= 0) ;
    zaxis = zaxis - 1.e-1 ;
    do {
      zaxis += 1.e-2 ;
    } while(ent.val_point(zaxis,0,0) >= 0) ;
    zaxis = zaxis - 1.e-2 ;
    do {
      zaxis += 1.e-3 ;
    } while(ent.val_point(zaxis,0,0) >= 0) ;
    zaxis = zaxis - 1.e-3 ;
    do {
      zaxis += 1.e-4 ;
    } while(ent.val_point(zaxis,0,0) >= 0) ;
    zaxis = zaxis - 1.e-4 ;
    do {
      zaxis += 1.e-5 ;
    } while(ent.val_point(zaxis,0,0) >= 0) ;
    zaxis = zaxis - 1.e-5 ;
    do {
      zaxis += 1.e-6 ;
    } while(ent.val_point(zaxis,0,0) >= 0) ;
    zaxis = zaxis - 1.e-6 ;

    arrete() ;

    double radius_z = zaxis ;
    cout << "Radius of the star at the north pole : " << radius_z
	 << endl ; 

    // Radius of the star
    //--------------------
    double radius_eq = mp.val_r(nz-2,1.,M_PI/2,0) ;
    cout << "Radius of the star in the equatorial plane : "
	 << radius_eq << endl ;

    cout << "Ratio of the axes a_3/a_1 : " << radius_z/radius_eq << endl ;
    cout << "Eccentricity : " << sqrt(1. - pow(radius_z / radius_eq, 2))
	 << endl ;

    // Angular velocity
    //------------------
    double omega = sqrt(omega2) ;
    cout << "Angular velocity of the star Omega/(pi G rho_c)^{1/2} : "
         << omega << endl ;
    cout << "Omega^2/(pi G rho_c) : " << omega2 << endl ;

    arrete() ;

    // Mass of the star
    //------------------
    double mass = rho.integrale() ;
    cout << "Mass/rho_c : " << mass << endl ;

    // Total energy of the star
    //--------------------------
    Cmp rho_n(mp) ;

    rho_n = pow(abs(rho), 1 + 1/n_index) ;
    rho_n.annule(nz-1) ;
    rho_n.set_dzpuis(4) ;
    rho_n.std_base_scal() ;

    double energy_internal = 4. * M_PI *
      (1 / (1 + n_index)) * rho_n.integrale() * radius_eq / (mass * mass) ;
    cout << "Internal energy of the axisymmetric star : "
	 << "E_int/(GM^2/radius_eq) : "
         << n_index * energy_internal << endl ;

    Cmp rho_pot(mp) ;
    rho_pot = rho * pot ;
    rho_pot.annule(nz-1) ;
    rho_pot.set_dzpuis(4) ;
    rho_pot.std_base_scal() ;
    double energy_selfgrav = -2. * M_PI * rho_pot.integrale() * radius_eq /
      (mass * mass) ;
    cout << "Self-gravity energy of the axisymmetric star : "
	 << "E_self/(GM^2/radius_eq) : "
         << energy_selfgrav << endl ;

    Cmp rr(mp) ;
    Cmp rho_quad(mp) ;
    rr = r * r * sint * sint ;
    rr.annule(nz-1) ;
    rr.set_dzpuis(4) ;
    rr.std_base_scal() ;
    rho_quad = rho * rr ;
    rho_quad.annule(nz-1) ;
    rho_quad.set_dzpuis(4) ;
    rho_quad.std_base_scal() ;
    double energy_kinetic = 0.5 * M_PI * omega2 * rho_quad.integrale() *
      radius_eq / (mass * mass) ;
    cout << "Kinetic energy of the axisymmetric star : "
	 << "E_kinet/(GM^2/radius_eq) : "
         << energy_kinetic << endl ;

    double energy_total = n_index * energy_internal + energy_selfgrav
      + energy_kinetic ;

    cout << "Total energy of the axisymmetric star : "
	 << "E/(GM^2/radius_eq) : "
         << energy_total << endl ;

    // Angular momentum
    //------------------
    double angmom = sqrt(M_PI) * omega * rho_quad.integrale() /
      sqrt( mass * mass * mass * radius_eq ) ;

    cout << "Angular momentum of the axisymmetric star : "
	 << "J/(GM^3 radius_eq)^{1/2} : "
	 << angmom << endl ;

    // Ratio T/|W|
    //-------------
    cout << "Ratio T/|W| : "
	 << energy_kinetic / fabs(energy_selfgrav) << endl ;

    // Virial relation
    //-----------------
    double virial = 3. * energy_internal + energy_selfgrav
      + 2. * energy_kinetic ;
    cout << "Virial relation : " << virial << endl ;

    //------------------------------------------
    //          Plot the figures
    //------------------------------------------
    double rmax ;
    cout << "r_max ?" << endl ;
    cin >> rmax ;

    des_profile(rho, 0., rmax, 0., 0., "rho/rho_c", "rho (z direction)") ;
    des_profile(rho, 0., rmax, M_PI/2, 0., "rho/rho_c", "rho (x-y plane)") ;
    des_profile(ent,0.,rmax,0.,0., "enthalpy", "enthalpy (z direction)") ;
    des_profile(ent,0.,rmax,M_PI/2,0., "enthalpy", "enthalpy (x-y plane)") ;
    des_coupe_x(rho, 0., -rmax, rmax, -rmax, rmax, "rho (y-z plane)") ;
    des_coupe_y(rho, 0., -rmax, rmax, -rmax, rmax, "rho (z-x plane)") ;
    des_coupe_z(rho, 0., -rmax, rmax, -rmax, rmax, "rho (x-y plane)") ;

    des_profile(rho_n,0.,rmax,0.,0.,"internal","internal (z direction)") ;
    des_profile(rho_n,0.,rmax,M_PI/2,0.,"internal","internal (x-y plane)") ;
    des_profile(rho_pot,0.,rmax,0.,0.,"selfgrav","selfgrav (z direction)") ;
    des_profile(rho_pot,0.,rmax,M_PI/2,0.,"selfgrav","selfgrav (x-y plane)") ;
    des_profile(rho_quad,0.,rmax,0.,0.,"kinetic","kinetic (z direction)") ;
    des_profile(rho_quad,0.,rmax,M_PI/2,0.,"kinetic","kinetic (x-y plane)") ;

    des_coef_xi(rho.va, 0, 0, 0, 1.e-14, "log|c_i|", "domain no. 0") ;
    des_coef_theta(rho.va, 0, 0, 0) ;
    des_coef_phi(rho.va, 0, 0, 0) ;

    des_coef_xi(rho.va, 1, 0, 0, 1.e-14, "log|c_i|", "domain no. 1") ;
    des_coef_theta(rho.va, 0, 0, 0) ;
    des_coef_phi(rho.va, 0, 0, 0) ;

    //    des_coef_xi(rho.va, 2, 0, 0, 1.e-14, "log|c_i|", "domain no. 2") ;
    //    des_coef_theta(rho.va, 0, 0, 0) ;
    //    des_coef_phi(rho.va, 0, 0, 0) ;

    exit(-1) ; 
        
}
