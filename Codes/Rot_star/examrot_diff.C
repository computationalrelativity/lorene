/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 2021      Jerome Novak
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
 * Revision 1.1  2021/04/14 16:39:54  j_novak
 * Added a code to read rotdiff output
 *
 *
 * $Header$
 *
 */

// headers C
#include <cstdlib>
#include <cstring>
#include <cmath>

// headers Lorene
#include "et_rot_diff.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "unites.h"

namespace Lorene {

// Function defining the rotation profile
double frotlin(double omega, const Lorene::Tbl& par) ; 
double primfrotlin(double omega, const Lorene::Tbl& par) ; 

// Local prototype (for drawings only)
Cmp raccord_c1(const Cmp& uu, int l1) ;

}
//******************************************************************************

using namespace Lorene ;

int main(int argc, char** argv){
  using namespace Unites ;

  if (argc < 2) {
       cout <<
     	"examrot_diff : the name of a file containing a stellar configuration"
     	   << endl << " must be given in argument !" << endl ;
       abort() ;
     }

     char* nomresu = argv[1] ;
     //  const char* nomresu = "resu.d" ;
    cout << "Name of the file to be read : " << nomresu << endl ;

    cout << endl <<
    "Do you want to draw the boundaries of the various domains (y/n) ? [y]"
	 << endl ;
    char rep ;
    cin.get(rep) ;
    bool draw_bound = !(rep == 'n') ;


    FILE* fich = fopen(nomresu, "r") ;   // open binary file in readonly

    Mg3d mg(fich) ;
    cout << mg ;
    Map_et mp(mg, fich) ;

    Eos* peos = Eos::eos_from_file(fich) ;

    Et_rot_diff star(mp, *peos, fich, &frotlin, &primfrotlin) ;

    fclose(fich) ;

    star.update_metric() ;
    star.equation_of_state() ;
    star.hydro_euler() ;

    cout.precision(10) ;
    cout << star ;

        // Formula (150) of Shapiro & Zane (1998) --> tsw1
    double tcin = 0.5 * star.get_omega_c() * star.angu_mom() ;
    double tsw1 = tcin / ( tcin + star.mass_b() - star.mass_g() ) ;

    cout << "tsw, tsw1 : " << star.tsw() << "  " << tsw1 << endl ;

    int nzet = star.get_nzet() ;

    des_profile(star.get_ak_car()(), 0., 3., 1, 1) ;

    // Drawing of the mapping
    // ----------------------

    for (int l=0; l<nzet; l++) {
    	des_map_et(mp, l) ;
    }	

    // 2-D drawings
    // ------------

	// Cmp defining the surface of the star (via the enthalpy field)
	Cmp surf = star.get_ent()() ;
	Cmp surf_ext(mp) ;
	surf_ext = - 0.2 * surf(0, 0, 0, 0) ;
	surf_ext.annule(0, star.get_nzet()-1) ;
	surf.annule(star.get_nzet(), mg.get_nzone()-1) ;
	surf = surf + surf_ext ;
	surf = raccord_c1(surf, star.get_nzet()) ;

	double zoom = 1.2 ;
	
	// Only the surface and possible transitions]
	cout << "Value of the enthalpy defining the transition ?" << endl ;
	double ent_trans ;
	cin >> ent_trans ;
	Cmp transit = star.get_ent()() - ent_trans ;

	// ... surface of the star
	int newgraph = 1 ; 	// opens the device
	double ray = zoom * star.ray_eq() ;
	des_surface_y(surf, 0., "?", newgraph, -ray, ray, -ray, ray,
		      "x [km]", "z [km]") ;

	// ... transition surface
	newgraph = 2 ; 	// closes the device at the end of the plot
	des_surface_y(transit, 0., "?", newgraph, -ray, ray, -ray, ray,
		      "x [km]", "z [km]") ;
		

	// Various fields
	des_coupe_y(star.get_ent()(), 0., nzet, "Enthalpy", &surf, zoom,
		    draw_bound) ;

	if (mg.get_np(0) > 1) {
	    des_coupe_z(star.get_ent()(), 0., nzet, "Enthalpy (equatorial plane)",
			&surf, zoom, draw_bound) ;
	}
	
	des_coupe_y(star.get_logn()(), 0., nzet, "\\gn", &surf, zoom,
		    draw_bound) ;	

	des_coupe_y(star.get_nphi()(), 0., nzet, "Azimuthal shift N\\u\\gf",
		    &surf, zoom, draw_bound) ;
	
	des_coupe_y(star.get_dzeta()(), 0., nzet, "Metric potential \\gz",
		    &surf, zoom, draw_bound) ;
	
	des_coupe_y(star.get_tggg()(), 0., nzet,
		    "Metric potential (NB-1) r sin\\gh", &surf, zoom,
		    draw_bound) ;
	
	des_coupe_y(star.get_ak_car()(), 0., nzet,
		    "A\\u2\\d K\\uij\\d K\\dij\\u", &surf, zoom,
		    draw_bound) ;

	// 1-D plots
	// ---------
	
	// Energy density
	double r_min = 0 ; 	// plot starts at center
	double r_max = star.ray_eq() ;
	double theta = M_PI / 2 ; 	
	double phi = 0 ;
	des_profile(star.get_ener()(), r_min, r_max, theta, phi,
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
		double uu1 = star.get_ener()().val_point(r,theta,phi) ;
		
		theta = 0 ;
		double uu2 = star.get_ener()().val_point(r,theta,phi) ;
		
		fdata << r / km << "  " << uu1 << "   " << uu2 ;
		
	
		fdata << endl ;
	 }

	 fdata.close() ;
	
	
    // Cleaning
    // --------

    delete peos ;


    return EXIT_SUCCESS ;

}

namespace Lorene {

                //------------------------------//
		//	Function F(Omega)	//
		//------------------------------//

  double frotlin(double omega, const Tbl& par){
    
    double omega_c = par(0) ; 
    double rrot = par(1) ; 
    
    return rrot*rrot* (omega_c - omega) ;  
  }

		//------------------------------//
		//	Primitive of  F(Omega)	//
		//------------------------------//

  double primfrotlin(double omega, const Tbl& par){
    
    double omega_c = par(0) ; 
    double rrot = par(1) ; 
    
    return - 0.5 * rrot*rrot* (omega_c - omega)*(omega_c - omega) ;  
  }
}
