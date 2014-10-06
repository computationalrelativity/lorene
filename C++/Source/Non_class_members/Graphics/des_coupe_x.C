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


char des_coupe_x_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/06 15:16:04  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.3  2004/03/25 10:29:24  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.2  2002/09/06 15:18:52  e_gourgoulhon
 * Changement du nom de la variable "hz" en "hza"
 * pour assurer la compatibilite avec le compilateur xlC_r
 * sur IBM Regatta (le preprocesseur de ce compilateur remplace
 * "hz" par "100").
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.8  2000/02/12  11:18:27  eric
 * Ajout de la version avec determination automatique des bornes de la fenetrraphique
 *
 * Revision 1.7  2000/02/11  18:43:53  eric
 * Ajout de l'argument draw_bound.
 *
 * Revision 1.6  2000/02/11  16:52:42  eric
 * Utilisation des coordonnees cartesiennes abolues (X,Y,Z) et non plus
 * des coordonnees relatives (x,y,z).
 *
 * Revision 1.5  1999/12/28  09:02:24  eric
 * *** empty log message ***
 *
 * Revision 1.4  1999/12/27  12:18:18  eric
 * Ajout du dessin des frontieres de domains.
 *
 * Revision 1.3  1999/12/24  13:00:20  eric
 * On appelle desormais la routine des_surface_x pour le dessin de la
 *  surface.
 *
 * Revision 1.2  1999/12/23  16:15:52  eric
 * Ajout du dessin de la surface (argument defsurf).
 *
 * Revision 1.1  1999/12/09  16:38:03  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Header C
#include <cmath>

// Header Lorene
#include "cmp.h"
#include "graphique.h"
#include "param.h"
#include "utilitaires.h"
#include "unites.h"

//******************************************************************************

void des_coupe_x(const Cmp& uu, double x0, int nzdes, const char* title, 
		 const Cmp* defsurf, double zoom, bool draw_bound, 
		 int ncour, int ny, int nz) {
		     
    const Map& mp = *(uu.get_mp()) ; 

    double a1 = mp.val_r(nzdes-1, 1., M_PI/2., 0.) ; 		 
    double a2 = mp.val_r(nzdes-1, 1., M_PI/2., M_PI/2.) ; 		 
    double a3 = mp.val_r(nzdes-1, 1., M_PI/2., M_PI) ; 		 
    double ray = mp.val_r(nzdes-1, 1., 0., 0.) ; 
    
    ray = ( a1 > ray ) ? a1 : ray ; 
    ray = ( a2 > ray ) ? a2 : ray ; 
    ray = ( a3 > ray ) ? a3 : ray ; 
    		 
    ray *= zoom ; 
    
    double y_min = mp.get_ori_y() - ray ; 
    double y_max = mp.get_ori_y() + ray ; 
    double z_min = mp.get_ori_z() - ray ; 
    double z_max = mp.get_ori_z() + ray ; 

    des_coupe_x(uu, x0, y_min, y_max, z_min, z_max, title, defsurf, draw_bound, 
		ncour, ny, nz) ;

}



//******************************************************************************

void des_coupe_x(const Cmp& uu, double x0, double y_min, double y_max, 
		 double z_min, double z_max, const char* title, const Cmp* defsurf, 
		 bool draw_bound, int ncour, int ny, int nz) {
		
using namespace Unites ;

    const Map& mp = *(uu.get_mp()) ; 

    // Plot of isocontours
    // -------------------
       
    float* uutab = new float[ny*nz] ; 
    
    double hy = (y_max - y_min) / double(ny-1) ; 
    double hza = (z_max - z_min) / double(nz-1) ; 
    
    for (int j=0; j<nz; j++) {
	
	double z = z_min + hza * j ; 
	
	for (int i=0; i<ny; i++) {
    
	    double y = y_min + hy * i ; 

	    // Computation of (r,theta,phi) : 	    
	    double r, theta, phi ; 
	    mp.convert_absolute(x0, y, z, r, theta, phi) ; 
	
	    uutab[ny*j+i] = float(uu.val_point(r, theta, phi)) ; 
	}
    }
    
    float ymin1 = float(y_min / km) ;
    float ymax1 = float(y_max / km) ;
    float zmin1 = float(z_min / km) ;
    float zmax1 = float(z_max / km) ;
    
    const char* nomy = "y [km]" ; 
    const char* nomz = "z [km]" ; 
    
    if (title == 0x0) {
	title = "" ;
    }
    
    const char* device = 0x0 ; 
    int newgraph = ( (defsurf != 0x0) || draw_bound ) ? 1 : 3 ; 
    
    des_equipot(uutab, ny, nz, ymin1, ymax1, zmin1, zmax1, ncour, nomy, nomz,
		title, device, newgraph) ;    

    delete [] uutab ;     
    
    // Plot of the surface
    // -------------------
    
    if (defsurf != 0x0) {

	assert(defsurf->get_mp() == uu.get_mp()) ; 

	newgraph = draw_bound ? 0 : 2 ;  
	
	des_surface_x(*defsurf, x0, device, newgraph) ; 
	
    }  // End of the surface drawing
    

    // Plot of the domains outer boundaries
    // ------------------------------------
    
    if (draw_bound) {

	int ndom = mp.get_mg()->get_nzone() ;  // total number of domains
    
	for (int l=0; l<ndom-1; l++) {	// loop on the domains (except the
					//  last one)

	    newgraph = (l == ndom-2) ? 2 : 0 ; 
	
	    des_domaine_x(mp, l, x0, device, newgraph) ; 
	}
    }

        
} 
