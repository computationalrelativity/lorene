/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *                 2012      Jerome Novak
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

char examrot_dirac_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/13 08:53:58  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2012/05/10 09:05:30  j_novak
 * New code examrot_dirac for reading the results of
 * rotstar_dirac. Simplification of the parrot.d parameter file for
 * rotstar_dirac.
 *
 *
 * $Header$
 *
 */

// headers C
#include <cmath>

// headers Lorene
#include "star_rot_dirac.h"
#include "cmp.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "unites.h"

namespace Lorene {
// Local prototype (for drawings only)
Cmp raccord_c1(const Cmp& uu, int l1) ;
}

//******************************************************************************

using namespace Lorene ;

int main(int argc, char** argv){

  using namespace Unites ;

    if (argc < 2) {
	cout <<
	    "examrot_dirac : the name of a file containing a stellar configuration"
	    << endl << " must be given in argument !" << endl ;
	abort() ;
    }

    char* nomresu = argv[1] ;
    cout << "Name of the file to be read : " << nomresu << endl ;

    cout << endl <<
    "Do you want to draw the boundaries of the various domains (y/n) ? [y]"
	 << endl ;
    char rep ;
    cin.get(rep) ;
    bool draw_bound = !(rep == 'n') ;


    FILE* fich = fopen(nomresu, "r") ;   // open binary file in readonly

    Mg3d mg(fich) ;
    Map_af mp(mg, fich) ;

    Eos* peos = Eos::eos_from_file(fich) ;

    Star_rot_Dirac star(mp, *peos, fich) ;

    fclose(fich) ;

    star.update_metric() ;
    star.equation_of_state() ;
    star.hydro_euler() ;


    cout.precision(10) ;
    cout << star << endl ;


    // Formula (150) of Shapiro & Zane (1998) --> tsw1
    double tcin = 0.5 * star.get_omega() * star.angu_mom() ;
    double tsw1 = tcin / ( tcin + star.mass_b() - star.mass_g() ) ;

    cout << "tsw, tsw1 : " << star.tsw() << "  " << tsw1 << endl ;


    // Print of the enthalpy field
    cout << "Enthalpy : " << endl ;
    cout << star.get_ent() << endl ;
    arrete() ;

    // Print of the energy density field
    cout << "Proper energy density : " << endl ;
    cout << star.get_ener() << endl ;
    arrete() ;


    int nzet = star.get_nzet() ;

    // 2-D drawings
    // ------------

    // Cmp defining the surface of the star (via the enthalpy field)
    Cmp surf = star.get_ent() ;
    Cmp surf_ext(mp) ;
    surf_ext = - 0.2 * surf(0, 0, 0, 0) ;
    surf_ext.annule(0, star.get_nzet()-1) ;
    surf.annule(star.get_nzet(), mg.get_nzone()-1) ;
    surf = surf + surf_ext ;
    surf = raccord_c1(surf, star.get_nzet()) ;

    double zoom = 1.2 ;
	
    // Various fields
    des_coupe_y(star.get_ent(), 0., nzet, "Enthalpy", &surf, zoom,
		    draw_bound) ;

    if (mg.get_np(0) > 1) {
      des_coupe_z(star.get_ent(), 0., nzet, "Enthalpy (equatorial plane)",
		  &surf, zoom, draw_bound) ;
    }
	
    des_coupe_y(star.get_logn(), 0., nzet, "Gravitational potential \\gn", &surf, zoom,
		draw_bound) ;
	
	
    des_coupe_y(star.get_beta()(3), 0., nzet, "Azimuthal shift \\gb\\u\\gf", &surf, zoom,
		draw_bound) ;
	
    des_coupe_y(star.get_lnq(), 0., nzet, "Metric potential ln(Q)", &surf, zoom,
		    draw_bound) ;
	

    des_coupe_y(star.get_aa_quad(), 0., nzet, "A\\uij\\d A\\dij", &surf, zoom,
		draw_bound) ;
    
    des_coupe_y(star.get_hh()(1,1), 0., nzet, "Metric potential h\\urr", &surf) ; 
    des_coupe_y(star.get_hh()(1,2), 0., nzet, "Metric potential h\\ur\\gh", &surf) ; 
    des_coupe_y(star.get_hh()(2,2), 0., nzet, "Metric potential h\\u\\gh\\gh", &surf) ; 
    des_coupe_y(star.get_hh()(3,3), 0., nzet, "Metric potential h\\u\\gf\\gf", &surf) ; 


    // 1-D plots
    // ---------
    
    // Energy density
    double r_min = 0 ; 	// plot starts at center
    double r_max = star.ray_eq() ;
    double theta = M_PI / 2 ; 	
    double phi = 0 ;
    des_profile(star.get_ener(), r_min, r_max, theta, phi,
		"Proper energy density [rho_nuc]") ;
    
    // Outputs in file for Xmgr, Gnuplot, ...
    
    ofstream fdata("profiles.d") ;
    fdata << "#   r [km]      e(theta=pi/2) 	e(theta=0)" << endl ;
    
    r_min = 0 ;
    r_max = star.ray_eq() ;
    int npt = 200 ;
    
    for (int i=0; i<npt; i++) {
      
      double r = r_min + double(i) * (r_max -r_min) / double(npt-1) ;
      
      theta = M_PI / 2. ;
      double uu1 = star.get_ener().val_point(r,theta,phi) ;
      
      theta = 0 ;
      double uu2 = star.get_ener().val_point(r,theta,phi) ;
      
      fdata << r / km << "  " << uu1 << "   " << uu2 ;
      
      
      fdata << endl ;
    }
    
    fdata.close() ;
    
    
    
    // Cleaning
    // --------
    
    delete peos ;

    exit(EXIT_SUCCESS) ;

    return EXIT_SUCCESS ;

}
