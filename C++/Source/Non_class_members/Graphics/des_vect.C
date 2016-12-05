/*
 * Basic routine for plotting a vector field.
 *
 */

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
 * Revision 1.6  2016/12/05 16:18:07  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:53:23  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:16:05  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.2  2002/10/16 14:36:58  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/03/01  16:12:02  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */



// C++ headers:
#include"headcpp.h"

// C headers:
#include <cmath>

// PGPLOT headers:
#include <cpgplot.h>

namespace Lorene {
//******************************************************************************

void des_vect(float* vvx, float* vvy, int nx, int ny, float xmin, float xmax, 
		 float ymin, float ymax, double scale,  double sizefl, 
		 const char* nomx, const char* nomy, const char* title, const char* device, 
		 int newgraph, int nxpage, int nypage) {
		 

    
    // Array defining the grid for pgvect_
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
     

    // Determination de la taille des fleches representant les vecteurs:
    // -----------------------------------------------------------------

    if (scale < 0) {		    // si scale >= 0, on garde la valeur de 
				    //  scale donnee en argument
	double norme_max = 0 ; 
	for(int ix = 0; ix < nx; ix++) {
	    for(int iy = 0; iy < ny; iy++) {
		double xxx = sqrt(  vvx[iy*nx+ix]*vvx[iy*nx+ix] +
				    vvy[iy*nx+ix]*vvy[iy*nx+ix]  ) ; 
		norme_max = (xxx > norme_max) ? xxx : norme_max ; 
	    }
	}
	
	if (norme_max < 1.e-14) {
	    scale = 1 ; 
	}
	else{
	    double pas_max = (hx > hy) ? hx : hy ;     
	    scale = fabs(scale) * pas_max / norme_max ; 
	}
	cout << "des_vect: norme_max = " << norme_max << endl ; 
	cout << "des_vect: scale = " << scale << endl ; 
    }


    // Graphics display
    // ----------------

    if ( (newgraph == 1) || (newgraph == 3) ) {

	if (device == 0x0) device = "?" ; 
   
	int ier = cpgbeg(0, device, nxpage, nypage) ;
	if (ier != 1) {
	cout << "des_vect: problem in opening PGPLOT display !" << endl ;
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



    float sizefl1 = float(sizefl) ; 
    cpgsch(sizefl1) ;    //	controle la taille des extremites des fleches 

    float blank = 0 ; 
    float scale1 = float(scale) ; 
    int nc = 1 ; 
    
    cpgvect(vvx, vvy, nx, ny, 1, nx, 1, ny, scale1, nc, tr, blank) ;   

    cpgsch(size) ;	// restauration 


    
    // Closing the graphical output
    // ----------------------------

    if ( (newgraph == 2) || (newgraph == 3) ) {    
	cpgend() ; 
    }
    
    
}
}
