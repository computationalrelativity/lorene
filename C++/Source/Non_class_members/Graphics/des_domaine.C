/*
 *  Basic routines for drawing the boundaries between domains.
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


char des_domaine_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.4  2001/02/28  09:45:02  eric
 * Correction erreur affichage "des_domaine_y" dans des_domaine_z.
 *
 * Revision 1.3  2000/02/11  16:53:50  eric
 * Utilisation des coordonnees cartesiennes abolues (X,Y,Z) et non plus
 * des coordonnees relatives (x,y,z).
 *
 * Revision 1.2  1999/12/27  12:20:59  eric
 * *** empty log message ***
 *
 * Revision 1.1  1999/12/27  12:19:03  eric
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
#include "map.h"
#include "param.h"
#include "utilitaires.h"

// Local prototypes
double fonc_des_domaine_x(double, const Param&) ; 
double fonc_des_domaine_y(double, const Param&) ; 
double fonc_des_domaine_z(double, const Param&) ; 

//******************************************************************************
 
void des_domaine_x(const Map& mp, int l0, double x0, char* device, int newgraph, 
		   double y_min, double y_max, double z_min, double z_max, 
		   char* nomy, char* nomz, char* title, int nxpage, int nypage)
{
#include "unites.h"
    // To avoid some compiler warnings :
    if (&mp == 0x0) {
	cout << f_unit << qpig << msol << mevpfm3 << endl ;
    }

    double khi ;
	
    Param parzerosec ;
    parzerosec.add_int(l0, 0) ; 	
    parzerosec.add_double_mod(x0, 0) ; 	
    parzerosec.add_double_mod(khi, 1) ; 	
    parzerosec.add_map(mp, 0) ; 	

    double rhomin = 0 ; 
    double rhomax = 2 * 
		    mp.val_r(mp.get_mg()->get_nzone() - 1, -1., 0., 0.) ;  
    double precis = 1.e-14 ; 
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

	if ( zero_premier(fonc_des_domaine_x, parzerosec, rhomin, rhomax, 100, 
		     rhomin0, rhomax0) == false ) {
	    cout << 
   "des_domaine_x : WARNING : no crossing with the domain boundary"
		<< endl ; 
	    cout << "  has been found for khi = " << khi << " !" << endl ; 

	    coupe_surface = false ; 
	    break ; 

	}
		
		     
	// Search for the zero in the interval [rhomin0, rhomax0] :
	
	double rho = zerosec(fonc_des_domaine_x, parzerosec, rhomin0, rhomax0, 
			     precis, nitermax, niter) ;
			       
	yg[i] = ( rho * cos(khi) + mp.get_ori_y() ) / km ; 	    
	zg[i] = ( rho * sin(khi) + mp.get_ori_z() ) / km ; 	    

    }
	
    // Graphics display
    // ----------------

    if ( (newgraph == 1) || (newgraph == 3) ) {

	if (device == 0x0) device = "?" ; 
   
	int ier = cpgbeg(0, device, nxpage, nypage) ;
	if (ier != 1) {
	cout << "des_domaine_x: problem in opening PGPLOT display !" << endl ;
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
	cpgsls(3) ;		// lignes en trait mixte
	cpgsci(3) ;		// couleur verte
	cpgline(np, yg, zg) ;
	cpgsls(1) ;		// retour aux lignes en trait plein
	cpgsci(1) ;		// couleur noire
    }
    
    
    // Closing graphic display
    // -----------------------

    if ( (newgraph == 2) || (newgraph == 3) ) {    
	cpgend() ; 
    }

}


//******************************************************************************

void des_domaine_y(const Map& mp, int l0, double y0, char* device, int newgraph, 
		   double x_min, double x_max, double z_min, double z_max, 
		   char* nomx, char* nomz, char* title, int nxpage, int nypage)
{
#include "unites.h"
    // To avoid some compiler warnings :
    if (&mp == 0x0) {
	cout << f_unit << qpig << msol << mevpfm3 << endl ;
    }

    double khi ;
	
    Param parzerosec ;
    parzerosec.add_int(l0, 0) ; 	
    parzerosec.add_double_mod(y0, 0) ; 	
    parzerosec.add_double_mod(khi, 1) ; 	
    parzerosec.add_map(mp, 0) ; 	

    double rhomin = 0 ; 
    double rhomax = 2 * 
		    mp.val_r(mp.get_mg()->get_nzone() - 1, -1., 0., 0.) ;  
    double precis = 1.e-14 ; 
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

	if ( zero_premier(fonc_des_domaine_y, parzerosec, rhomin, rhomax, 100, 
		     rhomin0, rhomax0) == false ) {
	    cout << 
   "des_domaine_y : WARNING : no crossing with the domain boundary"
		<< endl ; 
	    cout << "  has been found for khi = " << khi << " !" << endl ; 

	    coupe_surface = false ; 
	    break ; 

	}
		
		     
	// Search for the zero in the interval [rhomin0, rhomax0] :
	
	double rho = zerosec(fonc_des_domaine_y, parzerosec, rhomin0, rhomax0, 
			     precis, nitermax, niter) ;
			       
	xg[i] = ( rho * cos(khi) + mp.get_ori_x() ) / km ; 	    
	zg[i] = ( rho * sin(khi) + mp.get_ori_z() ) / km ; 	    

    }
	
    // Graphics display
    // ----------------

    if ( (newgraph == 1) || (newgraph == 3) ) {

	if (device == 0x0) device = "?" ; 
   
	int ier = cpgbeg(0, device, nxpage, nypage) ;
	if (ier != 1) {
	cout << "des_domaine_y: problem in opening PGPLOT display !" << endl ;
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
	cpgsls(3) ;		// lignes en trait mixte
	cpgsci(3) ;		// couleur verte
	cpgline(np, xg, zg) ;
	cpgsls(1) ;		// retour aux lignes en trait plein
	cpgsci(1) ;		// couleur noire
    }
    
    
    // Closing graphic display
    // -----------------------

    if ( (newgraph == 2) || (newgraph == 3) ) {    
	cpgend() ; 
    }

}

//******************************************************************************

void des_domaine_z(const Map& mp, int l0, double z0, char* device, int newgraph, 
		   double x_min, double x_max, double y_min, double y_max, 
		   char* nomx, char* nomy, char* title, int nxpage, int nypage)
{
#include "unites.h"
    // To avoid some compiler warnings :
    if (&mp == 0x0) {
	cout << f_unit << qpig << msol << mevpfm3 << endl ;
    }

    double khi ;
	
    Param parzerosec ;
    parzerosec.add_int(l0, 0) ; 	
    parzerosec.add_double_mod(z0, 0) ; 	
    parzerosec.add_double_mod(khi, 1) ; 	
    parzerosec.add_map(mp, 0) ; 	

    double rhomin = 0 ; 
    double rhomax = 2 * 
		    mp.val_r(mp.get_mg()->get_nzone() - 1, -1., 0., 0.) ;  
    double precis = 1.e-14 ; 
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

	if ( zero_premier(fonc_des_domaine_z, parzerosec, rhomin, rhomax, 100, 
		     rhomin0, rhomax0) == false ) {
	    cout << 
   "des_domaine_z : WARNING : no crossing with the domain boundary"
		<< endl ; 
	    cout << "  has been found for khi = " << khi << " !" << endl ; 

	    coupe_surface = false ; 
	    break ; 

	}
		
		     
	// Search for the zero in the interval [rhomin0, rhomax0] :
	
	double rho = zerosec(fonc_des_domaine_z, parzerosec, rhomin0, rhomax0, 
			     precis, nitermax, niter) ;
			       
	xg[i] = ( rho * cos(khi) + mp.get_ori_x() ) / km ; 	    
	yg[i] = ( rho * sin(khi) + mp.get_ori_y() ) / km ; 	    

    }
	
    // Graphics display
    // ----------------

    if ( (newgraph == 1) || (newgraph == 3) ) {

	if (device == 0x0) device = "?" ; 
   
	int ier = cpgbeg(0, device, nxpage, nypage) ;
	if (ier != 1) {
	cout << "des_domaine_z: problem in opening PGPLOT display !" << endl ;
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
	cpgsls(3) ;		// lignes en trait mixte
	cpgsci(3) ;		// couleur verte
	cpgline(np, xg, yg) ;
	cpgsls(1) ;		// retour aux lignes en trait plein
	cpgsci(1) ;		// couleur noire
    }
    
    
    // Closing graphic display
    // -----------------------

    if ( (newgraph == 2) || (newgraph == 3) ) {    
	cpgend() ; 
    }

}




//*****************************************************************************

double fonc_des_domaine_x(double vrho, const Param& par) {
    
    int l = par.get_int(0) ; 
    double x = par.get_double_mod(0) ; 
    double khi = par.get_double_mod(1) ; 
    const Map& mp = par.get_map(0) ; 
    
    // Absolute Cartesian coordinates: 
    double y = vrho * cos(khi) + mp.get_ori_y() ; 
    double z = vrho * sin(khi) + mp.get_ori_z() ; 

    // Spherical coordinates of the mapping:
    double r, theta, phi ; 
    mp.convert_absolute(x, y, z, r, theta, phi) ;
    
    return  r - mp.val_r(l, 1., theta, phi) ; 
    
}

//*****************************************************************************

double fonc_des_domaine_y(double vrho, const Param& par) {
    
    int l = par.get_int(0) ; 
    double y = par.get_double_mod(0) ; 
    double khi = par.get_double_mod(1) ; 
    const Map& mp = par.get_map(0) ; 
    
    // Absolute Cartesian coordinates: 
    double x = vrho * cos(khi) + mp.get_ori_x() ; 
    double z = vrho * sin(khi) + mp.get_ori_z() ; 

    // Spherical coordinates of the mapping:
    double r, theta, phi ; 
    mp.convert_absolute(x, y, z, r, theta, phi) ;
    
    return  r - mp.val_r(l, 1., theta, phi) ; 
    
}

//*****************************************************************************

double fonc_des_domaine_z(double vrho, const Param& par) {
    
    int l = par.get_int(0) ; 
    double z = par.get_double_mod(0) ; 
    double khi = par.get_double_mod(1) ; 
    const Map& mp = par.get_map(0) ; 
    
    // Absolute Cartesian coordinates: 
    double x = vrho * cos(khi) + mp.get_ori_x() ; 
    double y = vrho * sin(khi) + mp.get_ori_y() ; 

    // Spherical coordinates of the mapping:
    double r, theta, phi ; 
    mp.convert_absolute(x, y, z, r, theta, phi) ;
    
    return  r - mp.val_r(l, 1., theta, phi) ; 
    
}

