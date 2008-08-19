/*
 * Basic routine for drawing isocontours.
 *
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


char des_equipot_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.2  2002/10/16 14:36:57  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  1999/12/23  16:15:19  eric
 * Ajout des arguments newgraph, nxpage, nypage et device.
 *
 * Revision 1.1  1999/12/09  16:38:24  eric
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

// PGPLOT headers:
#include <cpgplot.h>

//******************************************************************************

void des_equipot(float* uutab, int nx, int ny, float xmin, float xmax, 
		 float ymin, float ymax, int ncour, const char* nomx, const char* nomy, 
		 const char* title, const char* device, int newgraph, int nxpage, 
		 int nypage) {
		 
    // Search for the extremal values of the field : 
    // -------------------------------------------

    float uumin = uutab[0] ;
    float uumax = uutab[0] ;
    for (int i=1; i<nx*ny; i++) {
	uumin = (uutab[i] < uumin) ? uutab[i] : uumin ;
	uumax = (uutab[i] > uumax) ? uutab[i] : uumax ;	
    }

    cout << "  " << title << " : min, max : " << uumin << "   " << uumax 
         << endl ; 

    // Values of equipotentials
    // -------------------------
 
    float* isopot = new float [ncour] ;
    float hh = (uumax-uumin) / float(ncour) ; 
    for (int i=0; i<ncour; i++) {
	isopot[i] = uumin + hh * float(i) ;
    }
    
    // Array defining the grid for pgcont_
    // -----------------------------------
    float hx = (xmax - xmin)/float(nx-1) ; 
    float hy = (ymax - ymin)/float(ny-1) ; 

    float tr[6] ;
    tr[0] = xmin - hx ;
    tr[1] = hx ;
    tr[2] = 0 ;
    tr[3] = ymin - hy ; 
    tr[4] = 0 ;
    tr[5] = hy ;
     
    // Graphics display
    // ----------------

    if ( (newgraph == 1) || (newgraph == 3) ) {

	if (device == 0x0) device = "?" ; 
   
	int ier = cpgbeg(0, device, nxpage, nypage) ;
	if (ier != 1) {
	cout << "des_equipot: problem in opening PGPLOT display !" << endl ;
	}

    }

    // Taille des caracteres:
    float size = float(1.3) ;
    cpgsch(size) ;
    
    // Epaisseur des traits:
    int lepais = 1 ; 
    cpgslw(lepais) ;
    
    // Fonte axes: caracteres romains:
    cpgscf(2) ;
    
    // Cadre de la figure
    cpgenv(xmin, xmax, ymin, ymax, 1, 0 ) ; 
    cpglab(nomx,nomy,title) ;

    // On n'effectue le dessin que si la dynamique est suffisante
    
    float dynamique = float(fabs(uumax - uumin)) ; 

    if (dynamique > 1.e-14) {
    
	cpgcont(uutab, nx, ny, 1, nx, 1, ny, isopot, ncour, tr) ;
	
    }
    
    // Closing the graphical output
    // ----------------------------

    if ( (newgraph == 2) || (newgraph == 3) ) {    
	cpgend() ; 
    }
    
    
    delete [] isopot ; 

}
