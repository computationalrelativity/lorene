/*
 *  Basic routine for drawing the surface of a star
 */

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


char des_surface_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.3  2000/03/21  11:14:45  eric
 * Le parametre precis est mis a 1e-8.
 *
 * Revision 1.2  2000/02/11  16:54:00  eric
 * Utilisation des coordonnees cartesiennes abolues (X,Y,Z) et non plus
 * des coordonnees relatives (x,y,z).
 *
 * Revision 1.1  1999/12/24  13:01:21  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// C headers:
#include <math.h>

// PGPLOT headers:
#include <cpgplot.h>

// Lorene headers
#include "cmp.h"
#include "param.h"
#include "utilitaires.h"

// Local prototypes
double fonc_des_surface_x(double, const Param&) ; 
double fonc_des_surface_y(double, const Param&) ; 
double fonc_des_surface_z(double, const Param&) ; 

//******************************************************************************
 
void des_surface_x(const Cmp& defsurf, double x0, char* device, int newgraph, 
		   double y_min, double y_max, double z_min, double z_max, 
		   char* nomy, char* nomz, char* title, int nxpage, int nypage)
{
#include "unites.h"
    // To avoid some compiler warnings :
    if (defsurf.get_etat() == ETATNONDEF) {
	cout << f_unit << qpig << msol << mevpfm3 << endl ;
    }

    assert(defsurf.get_etat() == ETATQCQ) ; 


    const Map* mp =  defsurf.get_mp(); 	
	
    double khi ;
	
    Param parzerosec ;
    parzerosec.add_double_mod(x0, 0) ; 	
    parzerosec.add_double_mod(khi, 1) ; 	
    parzerosec.add_cmp(defsurf) ;

    double rhomin = 0 ; 
    double rhomax = 2 * 
		    mp->val_r(mp->get_mg()->get_nzone() - 1, -1., 0., 0.) ;  
    double precis = 1.e-8 ; 
    int nitermax = 100 ; 
    int niter ; 
	
    const int np = 101 ; 
    float yg[np] ; 
    float zg[np] ; 

    double hkhi = 2 * M_PI / (np-1) ; 
    
    bool coupe_surface = true ; 
	
    for (int i=0; i< np; i++) {

	khi = hkhi * i ; 
	
	// Search for the interval [rhomin0, rhomax0] which contains 
	//  the first zero of des_surf:
	
	double rhomin0 ;
	double rhomax0 ;

	if ( zero_premier(fonc_des_surface_x, parzerosec, rhomin, rhomax, 100, 
		     rhomin0, rhomax0) == false ) {
	    cout << 
   "des_surface_x : WARNING : no interval containing a zero of defsurf"
		<< endl ; 
	    cout << "  has been found for khi = " << khi << " !" << endl ; 

	    coupe_surface = false ; 
	    break ; 

	}
		
		     
	// Search for the zero in the interval [rhomin0, rhomax0] :
	
	double rho = zerosec(fonc_des_surface_x, parzerosec, rhomin0, rhomax0, 
			     precis, nitermax, niter) ;
			       
	yg[i] = ( rho * cos(khi) + mp->get_ori_y() ) / km ; 	    
	zg[i] = ( rho * sin(khi) + mp->get_ori_z() ) / km ; 	    
    }
	
    // Graphics display
    // ----------------

    if ( (newgraph == 1) || (newgraph == 3) ) {

	if (device == 0x0) device = "?" ; 
   
	int ier = cpgbeg(0, device, nxpage, nypage) ;
	if (ier != 1) {
	cout << "des_surface_x: problem in opening PGPLOT display !" << endl ;
	} 
	
	// Taille des caracteres:
	float size = 1.3 ;
	cpgsch(size) ;
    
	// Epaisseur des traits:
	int lepais = 1 ; 
	cpgslw(lepais) ;
    
	cpgscf(2) ; // Fonte axes: caracteres romains
	
	float ymin1 = y_min / km ;
	float ymax1 = y_max / km ;
	float zmin1 = z_min / km ;
	float zmax1 = z_max / km ;
	
	cpgenv(ymin1, ymax1, zmin1, zmax1, 1, 0 ) ; 

	if (nomy == 0x0) nomy = "y [km]" ;
	if (nomz == 0x0) nomz = "z [km]" ; 
	if (title == 0x0) title = " " ; 
	cpglab(nomy,nomz,title) ;

    }

    if (coupe_surface) {
	cpgsls(1) ;		// lignes en trait plein
	cpgslw(6) ;		// traits gras
	cpgline(np, yg, zg) ;
	cpgslw(1) ;		// traits normaux
    }
    
    
    // Closing graphic display
    // -----------------------

    if ( (newgraph == 2) || (newgraph == 3) ) {    
	cpgend() ; 
    }

}

//******************************************************************************
 
void des_surface_y(const Cmp& defsurf, double y0, char* device, int newgraph, 
		   double x_min, double x_max, double z_min, double z_max, 
		   char* nomx, char* nomz, char* title, int nxpage, int nypage)
{
#include "unites.h"
    // To avoid some compiler warnings :
    if (defsurf.get_etat() == ETATNONDEF) {
	cout << f_unit << qpig << msol << mevpfm3 << endl ;
    }

    assert(defsurf.get_etat() == ETATQCQ) ; 


    const Map* mp =  defsurf.get_mp(); 	
	
    double khi ;
	
    Param parzerosec ;
    parzerosec.add_double_mod(y0, 0) ; 	
    parzerosec.add_double_mod(khi, 1) ; 	
    parzerosec.add_cmp(defsurf) ;

    double rhomin = 0 ; 
    double rhomax = 2 * 
		    mp->val_r(mp->get_mg()->get_nzone() - 1, -1., 0., 0.) ;  
    double precis = 1.e-8 ; 
    int nitermax = 100 ; 
    int niter ; 
	
    const int np = 101 ; 
    float xg[np] ; 
    float zg[np] ; 

    double hkhi = 2 * M_PI / (np-1) ; 
    
    bool coupe_surface = true ; 
	
    for (int i=0; i< np; i++) {

	khi = hkhi * i ; 
	
	// Search for the interval [rhomin0, rhomax0] which contains 
	//  the first zero of des_surf:
	
	double rhomin0 ;
	double rhomax0 ;

	if ( zero_premier(fonc_des_surface_y, parzerosec, rhomin, rhomax, 100, 
		     rhomin0, rhomax0) == false ) {
	    cout << 
   "des_surface_y : WARNING : no interval containing a zero of defsurf"
		<< endl ; 
	    cout << "  has been found for khi = " << khi << " !" << endl ; 

	    coupe_surface = false ; 
	    break ; 

	}
		
		     
	// Search for the zero in the interval [rhomin0, rhomax0] :
	
	double rho = zerosec(fonc_des_surface_y, parzerosec, rhomin0, rhomax0, 
			     precis, nitermax, niter) ;
			       
	xg[i] = ( rho * cos(khi) + mp->get_ori_x() ) / km ; 	    
	zg[i] = ( rho * sin(khi) + mp->get_ori_z() ) / km ; 	    
    }
	
    // Graphics display
    // ----------------

    if ( (newgraph == 1) || (newgraph == 3) ) {

	if (device == 0x0) device = "?" ; 
   
	int ier = cpgbeg(0, device, nxpage, nypage) ;
	if (ier != 1) {
	cout << "des_surface_y: problem in opening PGPLOT display !" << endl ;
	} 
	
	// Taille des caracteres:
	float size = 1.3 ;
	cpgsch(size) ;
    
	// Epaisseur des traits:
	int lepais = 1 ; 
	cpgslw(lepais) ;
    
	cpgscf(2) ; // Fonte axes: caracteres romains
	
	float xmin1 = x_min / km ;
	float xmax1 = x_max / km ;
	float zmin1 = z_min / km ;
	float zmax1 = z_max / km ;
	
	cpgenv(xmin1, xmax1, zmin1, zmax1, 1, 0 ) ; 

	if (nomx == 0x0) nomx = "x [km]" ;
	if (nomz == 0x0) nomz = "z [km]" ; 
	if (title == 0x0) title = " " ; 
	cpglab(nomx,nomz,title) ;

    }

    if (coupe_surface) {
	cpgsls(1) ;		// lignes en trait plein
	cpgslw(6) ;		// traits gras
	cpgline(np, xg, zg) ;
	cpgslw(1) ;		// traits normaux
    }
    
    
    // Closing graphic display
    // -----------------------

    if ( (newgraph == 2) || (newgraph == 3) ) {    
	cpgend() ; 
    }

}

//******************************************************************************
 
void des_surface_z(const Cmp& defsurf, double z0, char* device, int newgraph, 
		   double x_min, double x_max, double y_min, double y_max, 
		   char* nomx, char* nomy, char* title, int nxpage, int nypage)
{
#include "unites.h"
    // To avoid some compiler warnings :
    if (defsurf.get_etat() == ETATNONDEF) {
	cout << f_unit << qpig << msol << mevpfm3 << endl ;
    }

    assert(defsurf.get_etat() == ETATQCQ) ; 


    const Map* mp =  defsurf.get_mp(); 	
	
    double khi ;
	
    Param parzerosec ;
    parzerosec.add_double_mod(z0, 0) ; 	
    parzerosec.add_double_mod(khi, 1) ; 	
    parzerosec.add_cmp(defsurf) ;

    double rhomin = 0 ; 
    double rhomax = 2 * 
		    mp->val_r(mp->get_mg()->get_nzone() - 1, -1., 0., 0.) ;  
    double precis = 1.e-8 ; 
    int nitermax = 100 ; 
    int niter ; 
	
    const int np = 101 ; 
    float xg[np] ; 
    float yg[np] ; 

    double hkhi = 2 * M_PI / (np-1) ; 
    
    bool coupe_surface = true ; 
	
    for (int i=0; i< np; i++) {

	khi = hkhi * i ; 
	
	// Search for the interval [rhomin0, rhomax0] which contains 
	//  the first zero of des_surf:
	
	double rhomin0 ;
	double rhomax0 ;

	if ( zero_premier(fonc_des_surface_z, parzerosec, rhomin, rhomax, 100, 
		     rhomin0, rhomax0) == false ) {
	    cout << 
   "des_surface_z : WARNING : no interval containing a zero of defsurf"
		<< endl ; 
	    cout << "  has been found for khi = " << khi << " !" << endl ; 

	    coupe_surface = false ; 
	    break ; 

	}
		
		     
	// Search for the zero in the interval [rhomin0, rhomax0] :
	
	double rho = zerosec(fonc_des_surface_z, parzerosec, rhomin0, rhomax0, 
			     precis, nitermax, niter) ;
			       
	xg[i] = ( rho * cos(khi) + mp->get_ori_x() ) / km ; 	    
	yg[i] = ( rho * sin(khi) + mp->get_ori_y() ) / km ; 	    
    }
	
    // Graphics display
    // ----------------

    if ( (newgraph == 1) || (newgraph == 3) ) {

	if (device == 0x0) device = "?" ; 
   
	int ier = cpgbeg(0, device, nxpage, nypage) ;
	if (ier != 1) {
	cout << "des_surface_z: problem in opening PGPLOT display !" << endl ;
	} 
	
	// Taille des caracteres:
	float size = 1.3 ;
	cpgsch(size) ;
    
	// Epaisseur des traits:
	int lepais = 1 ; 
	cpgslw(lepais) ;
    
	cpgscf(2) ; // Fonte axes: caracteres romains
	
	float xmin1 = x_min / km ;
	float xmax1 = x_max / km ;
	float ymin1 = y_min / km ;
	float ymax1 = y_max / km ;
	
	cpgenv(xmin1, xmax1, ymin1, ymax1, 1, 0 ) ; 

	if (nomx == 0x0) nomx = "x [km]" ;
	if (nomy == 0x0) nomy = "y [km]" ; 
	if (title == 0x0) title = " " ; 
	cpglab(nomx,nomy,title) ;

    }

    if (coupe_surface) {
	cpgsls(1) ;		// lignes en trait plein
	cpgslw(6) ;		// traits gras
	cpgline(np, xg, yg) ;
	cpgslw(1) ;		// traits normaux
    }
    
    
    // Closing graphic display
    // -----------------------

    if ( (newgraph == 2) || (newgraph == 3) ) {    
	cpgend() ; 
    }

}










//*****************************************************************************

double fonc_des_surface_x(double vrho, const Param& par) {
    
    double x = par.get_double_mod(0) ; 
    double khi = par.get_double_mod(1) ; 
    const Cmp& defsurf = par.get_cmp() ; 
    const Map& mp = *(defsurf.get_mp()) ; 
    
    // Absolute Cartesian coordinates: 
    double y = vrho * cos(khi) + mp.get_ori_y() ;     
    double z = vrho * sin(khi) + mp.get_ori_z() ;     

    // Spherical coordinates of the mapping:
    double r, theta, phi ; 
    mp.convert_absolute(x, y, z, r, theta, phi) ;
    
    return defsurf.val_point(r, theta, phi) ; 
    
}

//*****************************************************************************

double fonc_des_surface_y(double vrho, const Param& par) {
    
    double y = par.get_double_mod(0) ; 
    double khi = par.get_double_mod(1) ; 
    const Cmp& defsurf = par.get_cmp() ; 
    const Map& mp = *(defsurf.get_mp()) ; 
    
    // Absolute Cartesian coordinates: 
    double x = vrho * cos(khi) + mp.get_ori_x() ; 
    double z = vrho * sin(khi) + mp.get_ori_z() ; 

    // Spherical coordinates of the mapping:
    double r, theta, phi ; 
    mp.convert_absolute(x, y, z, r, theta, phi) ;
    
    return defsurf.val_point(r, theta, phi) ; 
    
}

//*****************************************************************************

double fonc_des_surface_z(double vrho, const Param& par) {
    
    double z = par.get_double_mod(0) ; 
    double khi = par.get_double_mod(1) ; 
    const Cmp& defsurf = par.get_cmp() ; 
    const Map& mp = *(defsurf.get_mp()) ; 
    
    // Absolute Cartesian coordinates: 
    double x = vrho * cos(khi) + mp.get_ori_x() ; 
    double y = vrho * sin(khi) + mp.get_ori_y() ; 

    // Spherical coordinates of the mapping:
    double r, theta, phi ; 
    mp.convert_absolute(x, y, z, r, theta, phi) ;
    
    return defsurf.val_point(r, theta, phi) ; 
    
}
