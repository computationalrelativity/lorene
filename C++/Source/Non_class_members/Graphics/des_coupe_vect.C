/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char des_coupe_vect_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2014/10/13 08:53:22  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
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
 * Revision 2.0  2000/03/01  16:11:40  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Header C
#include <cmath>

// Header Lorene
#include "tenseur.h"
#include "graphique.h"
#include "param.h"
#include "utilitaires.h"
#include "unites.h"

namespace Lorene {
//******************************************************************************

void des_coupe_vect_x(const Tenseur& vv, double x0, double scale, double sizefl,
		      int nzdes, const char* title, const Cmp* defsurf, double zoom, 
		      bool draw_bound, int ny, int nz) {
		     
    const Map& mp = *(vv.get_mp()) ; 

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

    des_coupe_vect_x(vv, x0, scale, sizefl, y_min, y_max, z_min, z_max, title, 
		    defsurf, draw_bound, ny, nz) ;

}



//******************************************************************************

void des_coupe_vect_x(const Tenseur& vv, double x0, double scale, double
		      sizefl, double y_min, double y_max, double z_min, 
		      double z_max, const char* title, const Cmp* defsurf, 
		      bool draw_bound, int ny, int nz) {

  using namespace Unites ;
    
    const Map& mp = *(vv.get_mp()) ; 

    if (vv.get_valence() != 1) {
	cout << 
    "des_coupe_vect_x: the Tenseur must be of valence 1 (vector) !" << endl ;
	abort() ; 
    }

    if ( vv.get_triad()->identify() != mp.get_bvect_cart().identify() ) {
	cout << 
    "des_coupe_vect_x: the vector must be given in Cartesian components !" 
	<< endl ;
	abort() ; 
    }

    
    // Plot of the vector field
    // ------------------------
       
    float* vvy = new float[ny*nz] ; 
    float* vvz = new float[ny*nz] ; 
    
    double hy = (y_max - y_min) / double(ny-1) ; 
    double hza = (z_max - z_min) / double(nz-1) ; 
    
    for (int j=0; j<nz; j++) {
	
	double z = z_min + hza * j ; 
	
	for (int i=0; i<ny; i++) {
    
	    double y = y_min + hy * i ; 

	    // Computation of (r,theta,phi) : 	    
	    double r, theta, phi ; 
	    mp.convert_absolute(x0, y, z, r, theta, phi) ; 
	
	    vvy[ny*j+i] = float(vv(1).val_point(r, theta, phi)) ; 
	    vvz[ny*j+i] = float(vv(2).val_point(r, theta, phi)) ; 
	    
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
    
    des_vect(vvy, vvz, ny, nz, ymin1, ymax1, zmin1, zmax1,
	     scale,  sizefl, nomy, nomz, title, device, newgraph) ; 
		 
		 
    delete [] vvy ;     
    delete [] vvz ;     
    
    // Plot of the surface
    // -------------------
    
    if (defsurf != 0x0) {

	assert(defsurf->get_mp() == vv.get_mp()) ; 

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



//******************************************************************************

void des_coupe_vect_y(const Tenseur& vv, double y0, double scale, double sizefl,
		      int nzdes, const char* title, const Cmp* defsurf, double zoom, 
		      bool draw_bound, int nx, int nz) {
		     
    const Map& mp = *(vv.get_mp()) ; 

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
    double z_min = mp.get_ori_z() - ray ; 
    double z_max = mp.get_ori_z() + ray ; 


    des_coupe_vect_y(vv, y0, scale, sizefl, x_min, x_max, z_min, z_max, title, 
		    defsurf, draw_bound, nx, nz) ;

}



//******************************************************************************

void des_coupe_vect_y(const Tenseur& vv, double y0, double scale, double
		      sizefl, double x_min, double x_max, double z_min, 
		      double z_max, const char* title, const Cmp* defsurf, 
		      bool draw_bound, int nx, int nz) {
		
  using namespace Unites ;
    
    const Map& mp = *(vv.get_mp()) ; 

    if (vv.get_valence() != 1) {
	cout << 
    "des_coupe_vect_y: the Tenseur must be of valence 1 (vector) !" << endl ;
	abort() ; 
    }

    if ( vv.get_triad()->identify() != mp.get_bvect_cart().identify() ) {
	cout << 
    "des_coupe_vect_y: the vector must be given in Cartesian components !" 
	<< endl ;
	abort() ; 
    }

    
    // Plot of the vector field
    // ------------------------
       
    float* vvx = new float[nx*nz] ; 
    float* vvz = new float[nx*nz] ; 
    
    double hx = (x_max - x_min) / double(nx-1) ; 
    double hza = (z_max - z_min) / double(nz-1) ; 
    
    for (int j=0; j<nz; j++) {
	
	double z = z_min + hza * j ; 
	
	for (int i=0; i<nx; i++) {
    
	    double x = x_min + hx * i ; 
	    
	    // Computation of (r,theta,phi) : 	    
	    double r, theta, phi ; 
	    mp.convert_absolute(x, y0, z, r, theta, phi) ; 
	
	    vvx[nx*j+i] = float(vv(0).val_point(r, theta, phi)) ; 
	    vvz[nx*j+i] = float(vv(2).val_point(r, theta, phi)) ; 
	    
	}
    }
    
    float xmin1 = float(x_min / km) ;
    float xmax1 = float(x_max / km) ;
    float zmin1 = float(z_min / km) ;
    float zmax1 = float(z_max / km) ;
    
    const char* nomx = "x [km]" ; 
    const char* nomz = "z [km]" ; 
    
    if (title == 0x0) {
	title = "" ;
    }
    

    const char* device = 0x0 ; 
    int newgraph = ( (defsurf != 0x0) || draw_bound ) ? 1 : 3 ; 

    des_vect(vvx, vvz, nx, nz, xmin1, xmax1, zmin1, zmax1,
	     scale,  sizefl, nomx, nomz, title, device, newgraph) ; 
		 
		 
    delete [] vvx ;     
    delete [] vvz ;     
    
    // Plot of the surface
    // -------------------
    
    if (defsurf != 0x0) {

	assert(defsurf->get_mp() == vv.get_mp()) ; 

	newgraph = draw_bound ? 0 : 2 ;  
	
	des_surface_y(*defsurf, y0, device, newgraph) ; 
	
    }  // End of the surface drawing


    // Plot of the domains outer boundaries
    // ------------------------------------
    
    if (draw_bound) {

	int ndom = mp.get_mg()->get_nzone() ;  // total number of domains
    
	for (int l=0; l<ndom-1; l++) {	// loop on the domains (except the
					//  last one)

	    newgraph = (l == ndom-2) ? 2 : 0 ; 
	
	    des_domaine_y(mp, l, y0, device, newgraph) ; 
	}
    }

    
} 


//******************************************************************************

void des_coupe_vect_z(const Tenseur& vv, double z0, double scale, double sizefl,
		      int nzdes, const char* title, const Cmp* defsurf, double zoom, 
		      bool draw_bound, int nx, int ny) {
		     
    const Map& mp = *(vv.get_mp()) ; 

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

    des_coupe_vect_z(vv, z0, scale, sizefl, x_min, x_max, y_min, y_max, title, 
		    defsurf, draw_bound, nx, ny) ;

}



//******************************************************************************

void des_coupe_vect_z(const Tenseur& vv, double z0, double scale, double
		      sizefl, double x_min, double x_max, double y_min, 
		      double y_max, const char* title, const Cmp* defsurf, 
		      bool draw_bound, int nx, int ny) {
		
  using namespace Unites ;
    
    const Map& mp = *(vv.get_mp()) ; 

    if (vv.get_valence() != 1) {
	cout << 
    "des_coupe_vect_y: the Tenseur must be of valence 1 (vector) !" << endl ;
	abort() ; 
    }

    if ( vv.get_triad()->identify() != mp.get_bvect_cart().identify() ) {
	cout << 
    "des_coupe_vect_y: the vector must be given in Cartesian components !" 
	<< endl ;
	abort() ; 
    }

    
    // Plot of the vector field
    // ------------------------
       
    float* vvx = new float[nx*ny] ; 
    float* vvy = new float[nx*ny] ; 
    
    double hy = (y_max - y_min) / double(ny-1) ; 
    double hx = (x_max - x_min) / double(nx-1) ; 

    for (int j=0; j<ny; j++) {
	
	double y = y_min + hy * j ; 
	
	for (int i=0; i<nx; i++) {
    
	    double x = x_min + hx * i ; 
	    
	    // Computation of (r,theta,phi) : 	    
	    double r, theta, phi ; 
	    mp.convert_absolute(x, y, z0, r, theta, phi) ; 
	
	    vvx[nx*j+i] = float(vv(0).val_point(r, theta, phi)) ; 
	    vvy[nx*j+i] = float(vv(1).val_point(r, theta, phi)) ; 
	    
	}
    }
    
    float ymin1 = float(y_min / km) ;
    float ymax1 = float(y_max / km) ;
    float xmin1 = float(x_min / km) ;
    float xmax1 = float(x_max / km) ;
    
    const char* nomy = "y [km]" ; 
    const char* nomx = "x [km]" ; 
    
    if (title == 0x0) {
	title = "" ;
    }
    
    const char* device = 0x0 ; 
    int newgraph = ( (defsurf != 0x0) || draw_bound ) ? 1 : 3 ; 
    
    des_vect(vvx, vvy, nx, ny, xmin1, xmax1, ymin1, ymax1,
	     scale,  sizefl, nomx, nomy, title, device, newgraph) ; 
		 
		 
    delete [] vvx ;     
    delete [] vvy ;     
    
    // Plot of the surface
    // -------------------
    
    if (defsurf != 0x0) {

	assert(defsurf->get_mp() == vv.get_mp()) ; 

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
}
