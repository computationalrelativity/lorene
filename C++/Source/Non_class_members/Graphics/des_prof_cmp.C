/*
 * Draws the profile of a {\tt Cmp} along some radial axis determined by
 *  a fixed value of $(\theta, \phi)$.
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


char des_prof_cmp_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/06/03 09:59:35  e_gourgoulhon
 * Added a new des_profile with scale and nomx in arguments
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  1999/12/09  16:38:31  eric
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

void des_profile(const Cmp& uu, double r_min, double r_max, 
		     double theta, double phi, char* nomy, char* title) {
		
#include "unites.h"
    // To avoid some compiler warnings :
    if (r_min < 0) {
	cout << f_unit << qpig << msol << mevpfm3 << endl ;
    }
  

    const int npt = 400 ;   // Number of points along the axis
    
    float uutab[npt] ;	    // Value of uu at the npt points
    
    double hr = (r_max - r_min) / double(npt-1) ; 
    
    for (int i=0; i<npt; i++) {
    
	double r = hr * i + r_min ; 
	
	uutab[i] = uu.val_point(r, theta, phi) ; 
	
    }
    
    float xmin = r_min / km ;
    float xmax = r_max / km ;
    
    char* nomx = "r [km]" ; 
    
    if (title == 0x0) {
	title = "" ;
    }

    if (nomy == 0x0) {
	nomy = "" ;
    }
    
    
    des_profile(uutab, npt, xmin, xmax, nomx, nomy, title) ; 
    
} 

void des_profile(const Cmp& uu, double r_min, double r_max, double scale,
		     double theta, double phi, char* nomx, char* nomy, char* title) {
		

    const int npt = 400 ;   // Number of points along the axis
    
    float uutab[npt] ;	    // Value of uu at the npt points
    
    double hr = (r_max - r_min) / double(npt-1) ; 
    
    for (int i=0; i<npt; i++) {
    
	double r = hr * i + r_min ; 
	
	uutab[i] = uu.val_point(r, theta, phi) ; 
	
    }
    
    float xmin = r_min * scale ;
    float xmax = r_max * scale ;
    
    
    if (title == 0x0) {
	title = "" ;
    }

    if (nomx == 0x0) {
	nomx = "" ;
    }
    
    if (nomy == 0x0) {
	nomy = "" ;
    }
    
    
    des_profile(uutab, npt, xmin, xmax, nomx, nomy, title) ; 
    
} 
