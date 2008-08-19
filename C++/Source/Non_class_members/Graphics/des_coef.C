/* Basic routine for plot of spectral coefficients.
 *
 * (see file graphique.h for the documentation).
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


char des_coef_C[] = "$Header$" ;


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
 * Revision 1.2  1999/12/20  10:57:17  eric
 * Modif commentaires.
 *
 * Revision 1.1  1999/12/10  12:14:37  eric
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
 
void des_coef(const double* cf, int n, double pzero,
	      const char* nomx, const char* nomy, const char* title, const char* device, 
	      int newgraph, int nxpage, int nypage) {

    float xdes[2], ydes[2] ;

    double pzerol = log10(pzero) ; 
    ydes[0] = float(pzerol) ; 

    // Figure frame
    // ------------

    double xmin = - 1 ; 
    double xmax = n ; 
    double ymin = pzerol ; 
    double ymax = pzerol ; 
    
    for (int i=0; i<n; i++) {
	double yl = log10(fabs(cf[i])) ;
	if ( yl > ymax ) ymax = yl ; 
    }
    double yamp = ymax - ymin ;
    ymax = ymax + 0.05 * yamp ; 
    ymin = ymin - 0.05 * yamp ; 
    
    if (ymax <= pzerol) ymax = 2*pzerol ;   // assure que le cadre existe
    
    // Graphics display
    // ----------------

    if ( (newgraph == 1) || (newgraph == 3) ) {

	if (device == 0x0) device = "?" ; 
   
	int ier = cpgbeg(0, device, nxpage, nypage) ;
	if (ier != 1) {
	cout << "des_coef: problem in opening PGPLOT display !" << endl ;
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

    // Figure frame
    float xmin1 = float(xmin) ; 
    float xmax1 = float(xmax) ; 
    float ymin1 = float(ymin) ; 
    float ymax1 = float(ymax) ; 
    cpgenv(xmin1, xmax1, ymin1, ymax1, 0, 0 ) ; 
    cpglab(nomx, nomy, title) ;
    
    // Drawing of vertical lines to represent the coefficients
    for (int i=0; i<n; i++) {
	int lstyle ; 
	xdes[0] = float( i ) ;
	xdes[1] = xdes[0] ;
	if ( fabs(cf[i]) < pzero ) {
	    ydes[1] = float(pzerol) ; 
	    lstyle = 1 ; 
	}
	else {
	    ydes[1] = float( log10(fabs(cf[i])) ) ;
	    if (cf[i] < 0) lstyle = 2 ; 
	    else lstyle = 1 ; 
	}
	cpgsls(lstyle) ;	// line en trait plein (lstyle=1) 
			        //   ou pointilles (lstyle=2)
	cpgline(2, xdes, ydes) ; 
    }
       
    cpgsls(1) ;	    // restitution de la valeur par defaut 
    
    // Closing the graphical output
    // ----------------------------

    if ( (newgraph == 2) || (newgraph == 3) ) {    
	cpgend() ; 
    }
    
}


