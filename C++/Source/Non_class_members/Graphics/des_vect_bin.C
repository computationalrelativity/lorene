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


 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2016/12/05 16:18:07  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:53:23  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:16:05  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.3  2004/03/25 10:29:25  j_novak
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
 * Revision 2.1  2000/10/09  12:27:20  eric
 * Ajout du test sur l'identite des triades de vv1 et vv2.
 *
 * Revision 2.0  2000/03/02  10:34:24  eric
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

void des_vect_bin_x(const Tenseur& vv1, const Tenseur& vv2, double x0, 
		    double scale, double sizefl, double y_min, double y_max, 
		    double z_min, double z_max, const char* title, 
		    const Cmp* defsurf1, const Cmp* defsurf2, 
		    bool draw_bound, int ny, int nz) {
		
  using namespace Unites ;

    const Map& mp1 = *(vv1.get_mp()) ; 
    const Map& mp2 = *(vv2.get_mp()) ; 

    // Protections
    // -----------

    if (vv1.get_valence() != 1) {
	cout << 
    "des_vect_bin_x: the Tenseur vv1 must be of valence 1 (vector) !" << endl ;
	abort() ; 
    }

    if ( vv1.get_triad()->identify() != mp1.get_bvect_cart().identify() ) {
	cout << 
    "des_vect_bin_x: the vector vv1 must be given in Cartesian components !" 
	<< endl ;
	abort() ; 
    }

    if (vv2.get_valence() != 1) {
	cout << 
    "des_vect_bin_x: the Tenseur vv2 must be of valence 1 (vector) !" << endl ;
	abort() ; 
    }

    if ( vv2.get_triad()->identify() != mp2.get_bvect_cart().identify() ) {
	cout << 
    "des_vect_bin_x: the vector vv2 must be given in Cartesian components !" 
	<< endl ;
	abort() ; 
    }

    if ( *(vv1.get_triad()) != *(vv2.get_triad()) ) {
	cout << 
    "des_vect_bin_x: the components of the two vectors are not defined"
	    << endl << "on the same triad !" << endl ; 
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

	    double r, theta, phi ; 

	    mp1.convert_absolute(x0, y, z, r, theta, phi) ; 
	    double vv1y = vv1(1).val_point(r, theta, phi) ; 
	    double vv1z = vv1(2).val_point(r, theta, phi) ; 
	    
	    mp2.convert_absolute(x0, y, z, r, theta, phi) ; 
	    double vv2y = vv2(1).val_point(r, theta, phi) ; 
	    double vv2z = vv2(2).val_point(r, theta, phi) ; 
	     	
	    vvy[ny*j+i] = float(vv1y + vv2y) ; 
	    vvz[ny*j+i] = float(vv1z + vv2z) ; 
	    
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
    int newgraph = ( (defsurf1 != 0x0) || (defsurf2 != 0x0) || draw_bound ) ?
		    1 : 3 ; 
    
    des_vect(vvy, vvz, ny, nz, ymin1, ymax1, zmin1, zmax1,
	     scale,  sizefl, nomy, nomz, title, device, newgraph) ; 
		 
    delete [] vvy ;     
    delete [] vvz ;     
    
    // Plot of the surfaces
    // --------------------
    
    if (defsurf1 != 0x0) {

	assert(defsurf1->get_mp() == vv1.get_mp()) ; 
	newgraph = ( (defsurf2 != 0x0) || draw_bound ) ? 0 : 2 ;  
	des_surface_x(*defsurf1, x0, device, newgraph) ; 
    }
    
    if (defsurf2 != 0x0) {

	assert(defsurf2->get_mp() == vv2.get_mp()) ; 
	newgraph = draw_bound ? 0 : 2 ;  
	des_surface_x(*defsurf2, x0, device, newgraph) ; 
    }
    

    // Plot of the domains outer boundaries
    // ------------------------------------
    
    if (draw_bound) {

	int ndom1 = mp1.get_mg()->get_nzone() ;  
	int ndom2 = mp2.get_mg()->get_nzone() ;  
    
	for (int l=0; l<ndom1-1; l++) {	// loop on the domains (except the
					//  last one)
	    newgraph = 0 ;  
	    des_domaine_x(mp1, l, x0, device, newgraph) ; 
	}

	for (int l=0; l<ndom2-1; l++) {	// loop on the domains (except the
					//  last one)

	    newgraph = (l == ndom2-2) ? 2 : 0 ; 
	
	    des_domaine_x(mp2, l, x0, device, newgraph) ; 
	}

    }
} 


//******************************************************************************

void des_vect_bin_y(const Tenseur& vv1, const Tenseur& vv2, double y0, 
		    double scale, double sizefl, double x_min, double x_max, 
		    double z_min, double z_max, const char* title, 
		    const Cmp* defsurf1, const Cmp* defsurf2, 
		    bool draw_bound, int nx, int nz) {
		
  using namespace Unites ;

    const Map& mp1 = *(vv1.get_mp()) ; 
    const Map& mp2 = *(vv2.get_mp()) ; 

    // Protections
    // -----------

    if (vv1.get_valence() != 1) {
	cout << 
    "des_vect_bin_y: the Tenseur vv1 must be of valence 1 (vector) !" << endl ;
	abort() ; 
    }

    if ( vv1.get_triad()->identify() != mp1.get_bvect_cart().identify() ) {
	cout << 
    "des_vect_bin_y: the vector vv1 must be given in Cartesian components !" 
	<< endl ;
	abort() ; 
    }

    if (vv2.get_valence() != 1) {
	cout << 
    "des_vect_bin_y: the Tenseur vv2 must be of valence 1 (vector) !" << endl ;
	abort() ; 
    }

    if ( vv2.get_triad()->identify() != mp2.get_bvect_cart().identify() ) {
	cout << 
    "des_vect_bin_y: the vector vv2 must be given in Cartesian components !" 
	<< endl ;
	abort() ; 
    }

    if ( *(vv1.get_triad()) != *(vv2.get_triad()) ) {
	cout << 
    "des_vect_bin_y: the components of the two vectors are not defined"
	    << endl << "on the same triad !" << endl ; 
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
	    
	    double r, theta, phi ; 

	    mp1.convert_absolute(x, y0, z, r, theta, phi) ; 
	    double vv1x = vv1(0).val_point(r, theta, phi) ; 
	    double vv1z = vv1(2).val_point(r, theta, phi) ; 
	    
	    mp2.convert_absolute(x, y0, z, r, theta, phi) ; 
	    double vv2x = vv2(0).val_point(r, theta, phi) ; 
	    double vv2z = vv2(2).val_point(r, theta, phi) ; 
	     		
	    vvx[nx*j+i] = float(vv1x + vv2x) ; 
	    vvz[nx*j+i] = float(vv1z + vv2z) ; 
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
    int newgraph = ( (defsurf1 != 0x0) || (defsurf2 != 0x0) || draw_bound ) ?
		    1 : 3 ; 

    des_vect(vvx, vvz, nx, nz, xmin1, xmax1, zmin1, zmax1,
	     scale,  sizefl, nomx, nomz, title, device, newgraph) ; 
		 
    delete [] vvx ;     
    delete [] vvz ;     
    
    // Plot of the surfaces
    // --------------------
    
    if (defsurf1 != 0x0) {

	assert(defsurf1->get_mp() == vv1.get_mp()) ; 
	newgraph = ( (defsurf2 != 0x0) || draw_bound ) ? 0 : 2 ;  
	des_surface_y(*defsurf1, y0, device, newgraph) ; 
    }
    
    if (defsurf2 != 0x0) {

	assert(defsurf2->get_mp() == vv2.get_mp()) ; 
	newgraph = draw_bound ? 0 : 2 ;  
	des_surface_y(*defsurf2, y0, device, newgraph) ; 
    }
    

    // Plot of the domains outer boundaries
    // ------------------------------------
    
    if (draw_bound) {

	int ndom1 = mp1.get_mg()->get_nzone() ;  
	int ndom2 = mp2.get_mg()->get_nzone() ;  
    
	for (int l=0; l<ndom1-1; l++) {	// loop on the domains (except the
					//  last one)
	    newgraph = 0 ; 
    	    des_domaine_y(mp1, l, y0, device, newgraph) ; 
	}

	for (int l=0; l<ndom2-1; l++) {	// loop on the domains (except the
					//  last one)

	    newgraph = (l == ndom2-2) ? 2 : 0 ; 

	    des_domaine_y(mp2, l, y0, device, newgraph) ; 
	}

    }
} 


//******************************************************************************

void des_vect_bin_z(const Tenseur& vv1, const Tenseur& vv2, double z0, 
		    double scale, double sizefl, double x_min, double x_max, 
		    double y_min, double y_max, const char* title, 
		    const Cmp* defsurf1, const Cmp* defsurf2, 
		    bool draw_bound, int nx, int ny) {
		
  using namespace Unites ;

    const Map& mp1 = *(vv1.get_mp()) ; 
    const Map& mp2 = *(vv2.get_mp()) ; 

    // Protections
    // -----------

    if (vv1.get_valence() != 1) {
	cout << 
    "des_vect_bin_z: the Tenseur vv1 must be of valence 1 (vector) !" << endl ;
	abort() ; 
    }

    if ( vv1.get_triad()->identify() != mp1.get_bvect_cart().identify() ) {
	cout << 
    "des_vect_bin_z: the vector vv1 must be given in Cartesian components !" 
	<< endl ;
	abort() ; 
    }

    if (vv2.get_valence() != 1) {
	cout << 
    "des_vect_bin_z: the Tenseur vv2 must be of valence 1 (vector) !" << endl ;
	abort() ; 
    }

    if ( vv2.get_triad()->identify() != mp2.get_bvect_cart().identify() ) {
	cout << 
    "des_vect_bin_z: the vector vv2 must be given in Cartesian components !" 
	<< endl ;
	abort() ; 
    }

    if ( *(vv1.get_triad()) != *(vv2.get_triad()) ) {
	cout << 
    "des_vect_bin_z: the components of the two vectors are not defined"
	    << endl << "on the same triad !" << endl ; 
	abort() ; 
    }

    // Plot of the vector field
    // ------------------------
       
    float* vvx = new float[nx*ny] ; 
    float* vvy = new float[nx*ny] ; 
    
    double hx = (x_max - x_min) / double(nx-1) ; 
    double hy = (y_max - y_min) / double(ny-1) ; 

    for (int j=0; j<ny; j++) {
	
	double y = y_min + hy * j ; 
	
	for (int i=0; i<nx; i++) {
    
	    double x = x_min + hx * i ; 
	    
	    double r, theta, phi ; 

	    mp1.convert_absolute(x, y, z0, r, theta, phi) ; 
	    double vv1x = vv1(0).val_point(r, theta, phi) ; 
	    double vv1y = vv1(1).val_point(r, theta, phi) ; 
	    
	    mp2.convert_absolute(x, y, z0, r, theta, phi) ; 
	    double vv2x = vv2(0).val_point(r, theta, phi) ; 
	    double vv2y = vv2(1).val_point(r, theta, phi) ; 
		
	    vvx[nx*j+i] = float(vv1x + vv2x) ; 
	    vvy[nx*j+i] = float(vv1y + vv2y) ; 
	}
    }

    float xmin1 = float(x_min / km) ;
    float xmax1 = float(x_max / km) ;
    float ymin1 = float(y_min / km) ;
    float ymax1 = float(y_max / km) ;
    
    const char* nomx = "x [km]" ; 
    const char* nomy = "y [km]" ; 
    
    if (title == 0x0) {
	title = "" ;
    }
    
    const char* device = 0x0 ; 
    int newgraph = ( (defsurf1 != 0x0) || (defsurf2 != 0x0) || draw_bound ) ?
		    1 : 3 ; 

    des_vect(vvx, vvy, nx, ny, xmin1, xmax1, ymin1, ymax1,
	     scale,  sizefl, nomx, nomy, title, device, newgraph) ; 
		 
    delete [] vvx ;     
    delete [] vvy ;     
    
    // Plot of the surfaces
    // --------------------
    
    if (defsurf1 != 0x0) {

	assert(defsurf1->get_mp() == vv1.get_mp()) ; 
	newgraph = ( (defsurf2 != 0x0) || draw_bound ) ? 0 : 2 ;  
	des_surface_z(*defsurf1, z0, device, newgraph) ; 
    }
    
    if (defsurf2 != 0x0) {

	assert(defsurf2->get_mp() == vv2.get_mp()) ; 
	newgraph = draw_bound ? 0 : 2 ;  
	des_surface_z(*defsurf2, z0, device, newgraph) ; 
    }
    

    // Plot of the domains outer boundaries
    // ------------------------------------
    
    if (draw_bound) {

	int ndom1 = mp1.get_mg()->get_nzone() ;  
	int ndom2 = mp2.get_mg()->get_nzone() ;  
    
	for (int l=0; l<ndom1-1; l++) {	// loop on the domains (except the
					//  last one)
	    newgraph = 0 ; 
	    des_domaine_z(mp1, l, z0, device, newgraph) ; 
	}

	for (int l=0; l<ndom2-1; l++) {	// loop on the domains (except the
					//  last one)

	    newgraph = (l == ndom2-2) ? 2 : 0 ; 
	
	    des_domaine_z(mp2, l, z0, device, newgraph) ; 
	}

    }
} 
}
