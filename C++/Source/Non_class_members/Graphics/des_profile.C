/*
 *  Basic routine for drawing profiles.
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


char des_profile_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:57  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  1999/12/09  16:38:41  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */


// C++ headers:
#include"headcpp.h"

// C headers:
#include <math.h>
#include <cpgplot.h>

//******************************************************************************

void des_profile(float* uutab, int nx, float xmin, float xmax, 
		 char* nomx, char* nomy, char* title, char* device) {
		 
    // Search for the extremal values of the field : 
    // -------------------------------------------

    float uumin = uutab[0] ;
    float uumax = uutab[0] ;
    for (int i=1; i<nx; i++) {
	uumin = (uutab[i] < uumin) ? uutab[i] : uumin ;
	uumax = (uutab[i] > uumax) ? uutab[i] : uumax ;	
    }

    cout << "  " << nomy << " : min, max : " << uumin << "   " << uumax 
         << endl ; 

    // Points abscisses : 
    // ----------------
    
    float* xx = new float[nx] ; 
    float hx = (xmax-xmin)/float(nx-1) ;
    for(int i=0; i<nx; i++) {
	xx[i] = xmin + i * hx ; 
    }
         
    // Graphics display
    // ----------------
    
    if (device == 0x0) {
	device = "?" ; 
    }
    
    int ier = cpgbeg(0, device, 1, 1) ;
    if (ier != 1) {
	cout << "des_profile: problem in opening PGPLOT display !" << endl ;
    }
    
    // Taille des caracteres:
    float size = 1.3 ;
    cpgsch(size) ;
    
    // Epaisseur des traits:
    int lepais = 1 ; 
    cpgslw(lepais) ;
    
    // Fonte axes: caracteres romains:
    cpgscf(2) ;

    // Cadre de la figure
    float uuamp = uumax - uumin ; 
    float uumin1 = uumin - 0.05 * uuamp ; 
    float uumax1 = uumax + 0.05 * uuamp ; 
    cpgenv(xmin, xmax, uumin1, uumax1, 0, 0 ) ; 
    cpglab(nomx,nomy,title) ;
     
    cpgline(nx, xx, uutab) ; 
    
    cpgend() ; 
    
    delete [] xx ; 

}
