/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char des_coupe_z_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.8  2000/02/12  11:19:20  eric
 * Ajout de la version avec determination automatique des bornes de la fenetre
 * graphique.
 *
 * Revision 1.7  2000/02/11  18:44:13  eric
 * Ajout de l'argument draw_bound.
 *
 * Revision 1.6  2000/02/11  16:53:36  eric
 * Utilisation des coordonnees cartesiennes abolues (X,Y,Z) et non plus
 * des coordonnees relatives (x,y,z).
 *
 * Revision 1.5  1999/12/28  09:02:38  eric
 * *** empty log message ***
 *
 * Revision 1.4  1999/12/27  12:18:51  eric
 * Ajout du dessin des frontieres de domains.
 *
 * Revision 1.3  1999/12/24  13:01:02  eric
 * On appelle desormais la routine des_surface_y pour le dessin de la
 *   surface.
 *
 * Revision 1.2  1999/12/23  16:16:16  eric
 * Ajout du dessin de la surface (argument defsurf).
 *
 * Revision 1.1  1999/12/09  16:38:18  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Header C
#include <math.h>

// Header Lorene
#include "cmp.h"
#include "graphique.h"
#include "param.h"
#include "utilitaires.h"

//******************************************************************************

void des_coupe_z(const Cmp& uu, double z0, int nzdes, char* title, 
		 const Cmp* defsurf, double zoom, bool draw_bound, 
		 int ncour, int nx, int ny) {
		     
    const Map& mp = *(uu.get_mp()) ; 

    double a1 = mp.val_r(nzdes-1, 1., M_PI/2., 0.) ; 		 
    double a2 = mp.val_r(nzdes-1, 1., M_PI/2., M_PI/2.) ; 		 
    double a3 = mp.val_r(nzdes-1, 1., M_PI/2., M_PI) ; 		 
    double ray = mp.val_r(nzdes-1, 1., 0., 0.) ; 
    
    ray = ( a1 > ray ) ? a1 : ray ; 
    ray = ( a2 > ray ) ? a2 : ray ; 
    ray = ( a3 > ray ) ? a3 : ray ; 
    		 
    ray *= zoom ; 
    
    double x_min = mp.get_ori_x() - ray ; 
    double x_max = mp.get_ori_x() + ray ; 
    double y_min = mp.get_ori_y() - ray ; 
    double y_max = mp.get_ori_y() + ray ; 

    des_coupe_z(uu, z0, x_min, x_max, y_min, y_max, title, defsurf, draw_bound, 
		ncour, nx, ny) ;

}
//******************************************************************************


void des_coupe_z(const Cmp& uu, double z0, double x_min, double x_max, 
		 double y_min, double y_max, char* title, const Cmp* defsurf, 
		 bool draw_bound, int ncour, int nx, int ny) {
		
#include "unites.h"
    // To avoid some compiler warnings :
    if (uu.get_etat() == ETATNONDEF) {
	cout << f_unit << qpig << msol << mevpfm3 << endl ;
    }
  
    const Map& mp = *(uu.get_mp()) ; 

    // Plot of isocontours
    // -------------------
       
    float* uutab = new float[ny*nx] ; 
    
    double hy = (y_max - y_min) / double(ny-1) ; 
    double hx = (x_max - x_min) / double(nx-1) ; 

    for (int j=0; j<ny; j++) {
	
	double y = y_min + hy * j ; 
	
	for (int i=0; i<nx; i++) {
    
	    double x = x_min + hx * i ; 
	    
	    // Computation of (r,theta,phi) : 	    
	    double r, theta, phi ; 
	    mp.convert_absolute(x, y, z0, r, theta, phi) ; 
		
	    uutab[nx*j+i] = uu.val_point(r, theta, phi) ; 
	}
    }
    
    float ymin1 = y_min / km ;
    float ymax1 = y_max / km ;
    float xmin1 = x_min / km ;
    float xmax1 = x_max / km ;
    
    char* nomy = "y [km]" ; 
    char* nomx = "x [km]" ; 
    
    if (title == 0x0) {
	title = "" ;
    }
    
    char* device = 0x0 ; 
    int newgraph = ( (defsurf != 0x0) || draw_bound ) ? 1 : 3 ; 
    
    des_equipot(uutab, nx, ny, xmin1, xmax1, ymin1, ymax1, ncour, nomx, nomy,
		title, device, newgraph) ;    

    delete [] uutab ; 
    
        
    // Plot of the surface
    // -------------------
    
    if (defsurf != 0x0) {

	assert(defsurf->get_mp() == uu.get_mp()) ; 

	newgraph = draw_bound ? 0 : 2 ;  
	
	des_surface_z(*defsurf, z0, device, newgraph) ; 
	
    }  // End of the surface drawing

    // Plot of the domains outer boundaries
    // ------------------------------------
    
    if (draw_bound) {

	int ndom = mp.get_mg()->get_nzone() ;  // total number of domains
    
	for (int l=0; l<ndom-1; l++) {	// loop on the domains (except the
					//  last one)

	    newgraph = (l == ndom-2) ? 2 : 0 ; 
	
	    des_domaine_z(mp, l, z0, device, newgraph) ; 
	}
    }


} 
