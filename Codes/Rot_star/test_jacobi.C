/*
 * Test program for solving Jacobi-like figures
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
 * Revision 1.1  1999/12/24  13:32:14  keisuke
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
#include "param.h"

namespace Lorene {
// Local prototypes:
Cmp eos_local(const Cmp& ent, int nzet, double n_index) ;
double funct_zero_ent(double r, const Param& par) ;
}
//**********************************************************************

int main(){

    // Identification of all the subroutines called by the code : 

    system("ident test_jacobi > identif_jacobi.d") ; 

    //-----------------------------------------------------------------------
    //			Input from file "part.d"
    //-----------------------------------------------------------------------
    
    int type_t, type_p, nt, np, nz, l;
    char blabla[80];
    double n_index, radius_x, radius_y ;

    //----------------------------------------------------
    //             Declear the polytropic index
    //----------------------------------------------------
    cout << "Input a polytropic index: ";
    cin >> n_index;

    cout << "Input the radius_x of the star : " ;
    cin >> radius_x ;
    cout << "Input the radius_y of the star (< radius_x ) : " ;
    cin >> radius_y ;

    ifstream fich("part.d") ;
    fich >> nt; fich.getline(blabla, 80);
    fich >> np; fich.getline(blabla, 80);
    fich >> nz; fich.getline(blabla, 80);

    cout << "nb of points in phi : np = " << np << endl;
    cout << "nb of points in theta : nt = " << nt << endl;
    cout << "nb of zones : nz = " << nz << endl;

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

    // Obtaining "lambda" in order to rescale the mapping
    //----------------------------------------------------
    double lambda = radius_x ;

    // Rescale of the map
    //--------------------
    mp.homothetie( lambda ) ; 

    Cmp rho(mp) ; 
    Cmp rho1(mp) ;
    Cmp rho2(mp) ;
    Cmp circle1(mp) ;
    Cmp circle2(mp) ;

    const Coord& r = mp.r ; 
    const Coord& x = mp.x ; 
    //    const Coord& y = mp.y ; 
    const Coord& sint = mp.sint ; 

    circle1 = pow(abs(r*r - x + 0.25), 0.5) ;
    circle1.annule(nz-1) ;
    circle1.set_dzpuis(4) ;
    circle1.std_base_scal() ;

    rho1.set_etat_qcq() ;
    (rho1.va).set_etat_c_qcq() ;
    (rho1.va).c->set_etat_qcq() ;

    circle2 = pow(abs(r*r + x + 0.25), 0.5) ;
    circle2.annule(nz-1) ;
    circle2.set_dzpuis(4) ;
    circle2.std_base_scal() ;

    rho2.set_etat_qcq() ;
    (rho2.va).set_etat_c_qcq() ;
    (rho2.va).c->set_etat_qcq() ;

    int nr_s = mg.get_nr(nz-2) ;

    //-------------------------------------------
    //         Initial densty profile
    //-------------------------------------------

    for (int l=0; l<nz; l++) {

      (rho1.va).c->t[l]->set_etat_qcq() ;
      (rho2.va).c->t[l]->set_etat_qcq() ;

      for (int k=0; k<np; k++) {
	for (int j=0; j<nt; j++) {
	  for (int i=0; i<nr_s; i++) {

	    double circ1 = circle1(l, k, j, i) ;

	    if (circ1 <= 0.5) {
	      rho1.set(l, k, j, i) = 1 ;
	    }
	    else {
	      rho1.set(l, k, j, i) = 0 ;
	    }

	    double circ2 = circle2(l, k, j, i) ;
	    if (circ2 <= 0.5) {
	      rho2.set(l, k, j, i) = 1 ;
	    }
	    else {
	      rho2.set(l, k, j, i) = 0 ;
	    }

	  }
	}
      }
    }

    rho = rho1 + rho2  ;
    rho.annule(nz-1) ;

    rho.set_dzpuis(4) ; 
    rho.std_base_scal() ; // Sets the standard basis for spectral expansions

    Cmp pot(mp) ;
    Cmp rho_prev(mp) ;
    Cmp rot(mp) ;

    Cmp ent(mp) ;
    ent = pow(rho, 1/n_index) ;
    ent.annule(nz-1) ;
    ent.set_dzpuis(4) ;
    ent.std_base_scal() ;

    double omega2 ;
    double omega2_prev ;
    double rel_omega ;
    omega2 = 0 ;

    double cc ;
    double cc_prev ;
    double rel_cc ;
    cc = 0 ;

    double rmax ;
    cout << "r_max ?" << endl ;
    cin >> rmax ;

    do {

      /*
      des_coupe_x(rho,0.,-rmax,rmax,-rmax,rmax,"rho (y-z plane)",&ent) ;
      des_coupe_y(rho,0.,-rmax,rmax,-rmax,rmax,"rho (z-x plane)",&ent) ;
      des_coupe_z(rho,0.,-rmax,rmax,-rmax,rmax,"rho (x-y plane)",&ent) ;
      */

      // Solving the Poisson equation
      //------------------------------
      pot = - rho.poisson() ;


      // Surface potentials
      //--------------------
      // surface value of potential in the x-axis direction
      //      double pot_sx_2 = pot(nz-2,0,nt-1,nr_s-1) ;
      double pot_sx = pot.val_point(radius_x,M_PI/2,0) ;

      // Computation of the radius to the y-axis
      //-----------------------------------------

      Param par_funct_zero_y ;
      par_funct_zero_y.add_cmp(ent) ; 
      double phi_search_y = M_PI/2 ;
      double theta_search_y = M_PI/2 ;
      par_funct_zero_y.add_double(theta_search_y, 0) ;
      par_funct_zero_y.add_double(phi_search_y, 1) ;

      double r_min_search_y = 0 ;
      double r_max_search_y = radius_x ;
      double precis_y = 1e-10 ;
      int nitermax_y = 100 ;
      int niter_y ;

      double r_y = zerosec(funct_zero_ent, par_funct_zero_y, 
			   r_min_search_y, r_max_search_y, precis_y,
			   nitermax_y, niter_y) ;

      cout << "Radius to the y-axis direction : " << r_y << endl ;

      // surface value of potential in the y-axis direction
      double pot_sy = pot.val_point(radius_y,M_PI/2,M_PI/2) ;

      /*
      double rel_rad_y = fabs(1. - radius_y / r_y) ;
      cout << "Relative error in radius_y : " << rel_rad_y << endl ;

      cout << "The surface values of the potential : " << endl
	   << pot_sx << endl << pot_sx_2 << endl << pot_sy << endl ;

      double radi_x = mp.val_r(nz-2,1.,M_PI/2,0) ;
      double rel_rad_x = fabs(1. - radius_x / radi_x) ;
      cout << "Relative error in radius_x : " << rel_rad_x << endl ;
      */


      // Angular velocity
      //------------------
      omega2_prev = omega2 ;

      omega2 = 8. * (- pot_sx + pot_sy) /
	(radius_x * radius_x - radius_y * radius_y) ;
      cout << "Omega^2/(pi G rho_c) : " << omega2 << endl ;

      rel_omega = fabs(1. - omega2_prev / omega2) ;
      cout << "Relative error in Omega^2 : " << rel_omega << endl ;

      // Integration constant
      //----------------------
      cc_prev = cc ;

      /*
      cc = (radius_y * radius_y * pot_sx
	    - radius_x * radius_x * pot_sy) /
	(radius_x * radius_x - radius_y * radius_y) ;
	*/

      cc = - pot_sx - 0.125 * omega2 * radius_x * radius_x ;

      cout << "Integration constant : " << cc << endl ;

      rel_cc = fabs(1. - cc_prev / cc) ;
      cout << "Relative error in cc : " << rel_cc << endl ;

      // First integral of motion
      //--------------------------
      rot = 0.125 * omega2 * r * r * sint * sint ;
      rot.annule(nz-1) ; 
      ent = pot + rot + cc ; // enthalpy
      ent.annule(nz-1) ;
      ent.set_dzpuis(4) ;
      ent.std_base_scal() ;

      cout << "The central value of the enthalpy : " << ent(0,0,0,0)
	   << endl ;
      /*
      cout << "The surface values of the enthalpy : " << endl
	   << ent.val_point(radius_x,M_PI/2,0) << endl
	   << ent(nz-2,0,nt-1,nr_s-1) << endl
	   << ent(nz-1,0,nt-1,0) << endl
	   << ent.val_point(radius_y,M_PI/2,M_PI/2) << endl
	   << ent.val_point(radius_y,M_PI/2,0) << endl ;
	   */

      rho_prev = rho ;

      // EOS : rho = rho(H)
      //--------------------
      rho = eos_local(ent, nz-1, n_index) ;
      rho.set_dzpuis(4) ;
      rho.std_base_scal() ;

      cout << max(abs(rho - rho_prev)) << endl ;

      cout << "Maximum difference between rho and rho_prev : "
	   << max(diffrel(rho, rho_prev)) << endl ;

      arrete() ;

    } while(max(diffrel(rho, rho_prev)) > 1.e-4) ;

    arrete() ;

    cout << "Coef of rho : " << endl ; 
    rho.affiche_seuil(cout) ; 

    cout << "Coef of pot : " << endl ; 
    pot.affiche_seuil(cout) ; 

    arrete() ;
    cout << "Value of pot at the origin : " << pot(0,0,0,0) << endl ;
    cout << "Value of pot at the surface (r=ray) : " << pot(nz-1,0,0,0)
	 << endl ;
    cout << "Value of pot at the surface (r=ray) : " << pot(nz-2,0,0,nr_s-1) 
	 << endl ; 

    // Checking whether or not the potential satisfies the Poisson equation
    //----------------------------------------------------------------------
    Cmp lap = - pot.laplacien() ; 
    
    cout << "max( |lap(pot) - rho| ) " << max(abs(lap - rho)) << endl ; 
    cout << "relative error : " << diffrel(lap, rho) << endl ;

    // Searching the radius to the z-axis direction
    //----------------------------------------------
    Param par_funct_zero_z ;
    par_funct_zero_z.add_cmp(ent) ; 
    double phi_search_z = 0 ;
    double theta_search_z = 0 ;
    par_funct_zero_z.add_double(theta_search_z, 0) ;
    par_funct_zero_z.add_double(phi_search_z, 1) ;

    double r_min_search_z = 0 ;
    double r_max_search_z = radius_x ;
    double precis_z = 1e-10 ;
    int nitermax_z = 100 ;
    int niter_z ;

    double r_z = zerosec(funct_zero_ent, par_funct_zero_z, 
			   r_min_search_z, r_max_search_z, precis_z,
			   nitermax_z, niter_z) ;

    // Radius of the star
    //--------------------
    double ray = mp.val_r(nz-2,1.,M_PI/2,0) ;
    cout << "Radius of the star in the x-axis direction : " << ray << endl ;

    cout << "Radius of the star in the y-axis direction : " << radius_y
	 << endl ;

    cout << "Radius of the star in the pole : " << r_z
	 << endl ;

    cout << "The ratio of the axis : a_2/a_1 : " << radius_y/ray << endl
	 << "The ratio of the axis : a_3/a_1 : " << r_z/ray << endl ;

    // Angular velocity
    //------------------
    double omega = sqrt(omega2) ;

    cout << "Angular velocity of the star Omega/(pi G rho_c)^{1/2} : "
	 << omega << endl ;
    cout << "Angular velocity Omega^2/(pi G rho_c) : "
	 << omega2 << endl ;

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

    double energy_internal = 4 * M_PI *
      (1 / (1 + n_index)) * rho_n.integrale() * ray / (mass * mass) ;
    cout << "Internal energy of the triaxial star : E_int/(GM^2/ray) : "
	 << n_index * energy_internal << endl ;

    Cmp rho_pot(mp) ;
    rho_pot = rho * pot ;
    rho_pot.annule(nz-1) ;
    rho_pot.set_dzpuis(4) ;
    rho_pot.std_base_scal() ;
    double energy_selfgrav = -2. * M_PI * rho_pot.integrale() * ray /
      (mass * mass) ;
    cout <<
      "Self-gravity energy of the triaxial star : E_self/(GM^2/ray) : "
	 << energy_selfgrav << endl ;

    Cmp rr(mp) ;
    Cmp rho_quad(mp) ;
    rr = r * r * sint * sint ;
    rho_quad = rho * rr ;
    rho_quad.annule(nz-1) ;
    rho_quad.set_dzpuis(4) ;
    rho_quad.std_base_scal() ;
    double energy_kinetic = 0.5 * M_PI * omega2 * rho_quad.integrale() *
      ray / (mass * mass) ;
    cout << "Kinetic energy of the triaxial star : E_kinet/(GM^2/ray) : "
	 << energy_kinetic << endl ;

    double energy_total = n_index * energy_internal + energy_selfgrav
      + energy_kinetic ;

    cout << "Total energy of the triaxial star : E/(GM^2/ray) : "
	 << energy_total << endl ;

    // Angular momentum
    //------------------
    double angmom = sqrt(M_PI) * omega * rho_quad.integrale() /
      sqrt( mass * mass * mass * ray) ;

    cout <<
      "Angular momentum of the triaxial star : J/(GM^3 ray)^{1/2} : "
	 << angmom << endl ;

    // Virial relation
    //-----------------
    double virial = 3. * energy_internal + energy_selfgrav
      + 2. * energy_kinetic ;
    cout << "Virial relation : " << virial << endl ;

    //------------------------------------------
    //          Plot the figures
    //------------------------------------------
    //    double rmax ;
    //    cout << "r_max ?" << endl ;
    //    cin >> rmax ;

    des_profile(rho,0.,rmax,M_PI/2, 0.,"rho/rho_c","rho (x direction)") ;
    des_profile(rho,0.,rmax,M_PI/2,M_PI/2,"rho/rho_c","rho (y direction)") ;
    des_profile(rho,0.,rmax,0.,0.,"rho/rho_c","rho (z direction)") ;
    des_profile(ent,0.,rmax,M_PI/2,0.,"ent","enthalpy (x direction)") ;
    des_profile(ent,0.,rmax,M_PI/2,M_PI/2,"ent","enthalpy (y direction)") ;
    des_profile(ent,0.,rmax,0.,0.,"ent","enthalpy (z direction)") ;

    des_coupe_x(rho,0.,-rmax,rmax,-rmax,rmax,"rho (x-z plane)",&ent) ;
    des_coupe_y(rho,0.,-rmax,rmax,-rmax,rmax,"rho (y-z plane)",&ent) ;
    des_coupe_z(rho,0.,-rmax,rmax,-rmax,rmax,"rho (x-y plane)",&ent) ;

    des_coef_xi(rho.va, 0, 0, 0, 1.e-14, "log|c_i|", "domain no. 0") ;
    des_coef_theta(rho.va, 0, 0, 0) ;
    des_coef_phi(rho.va, 0, 0, 0) ;

    /*
    des_coef_xi(rho.va, 1, 0, 0, 1.e-14, "log|c_i|", "domain no. 1") ;
    des_coef_theta(rho.va, 1, 0, 0) ;
    des_coef_phi(rho.va, 1, 0, 0) ;
    */
    exit(-1) ; 

}

