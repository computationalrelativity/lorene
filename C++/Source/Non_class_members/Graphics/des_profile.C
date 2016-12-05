/*
 *  Basic routine for drawing profiles.
 */

/*
 *   Copyright (c) 1999-2004 Eric Gourgoulhon
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
 * Revision 1.12  2016/12/05 16:18:06  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.11  2014/10/13 08:53:22  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.10  2014/10/06 15:16:05  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.9  2012/01/17 10:35:46  j_penner
 * added point plot
 *
 * Revision 1.8  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.7  2005/03/25 19:56:28  e_gourgoulhon
 * Added plot of domain boundaries (new arguments nbound and xbound).
 *
 * Revision 1.6  2004/02/17 22:19:22  e_gourgoulhon
 * Changed prototype of des_profile_mult.
 * Added version of des_profile_mult with arbitrary x sampling.
 *
 * Revision 1.5  2004/02/16 10:54:08  e_gourgoulhon
 * Added #include <stdlib.h>.
 *
 * Revision 1.4  2004/02/15 21:56:49  e_gourgoulhon
 * des_profile_mult: added call to cpgask(0).
 *
 * Revision 1.3  2004/02/12 16:21:57  e_gourgoulhon
 * Added new function des_profile_mult.
 *
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
#include <cstdlib>
#include <cmath>
#include <cpgplot.h>


namespace Lorene {
//******************************************************************************
//      Single profile, single device, uniform sampling 
//******************************************************************************

void des_profile(const float* uutab, int nx, float xmin, float xmax, 
		 const char* nomx, const char* nomy, 
                 const char* title, const char* device,
                 int nbound, float* xbound) {
		 
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
	xx[i] = xmin + float(i) * hx ; 
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
    float size = float(1.3) ;
    cpgsch(size) ;
    
    // Epaisseur des traits:
    int lepais = 1 ; 
    cpgslw(lepais) ;
    
    // Fonte axes: caracteres romains:
    cpgscf(2) ;

    // Cadre de la figure
    float uuamp = uumax - uumin ; 
    float uumin1 = uumin - float(0.05) * uuamp ; 
    float uumax1 = uumax + float(0.05) * uuamp ; 
    cpgenv(xmin, xmax, uumin1, uumax1, 0, 0 ) ; 
    cpglab(nomx,nomy,title) ;
    
    // Drawing of curve 
    cpgline(nx, xx, uutab) ; 
    
    
    // Plot of domain boundaries
    // -------------------------
    
    if (nbound > 0) {
        float xb[2] ; 
        float yb[2] ; 
        yb[0] = uumin1 ; 
        yb[1] = uumax1 ; 
        cpgsls(3) ;		// lignes en trait mixte
        cpgsci(3) ;		// couleur verte
        for (int i=0; i<nbound; i++) {
            xb[0] = xbound[i] ; 
            xb[1] = xbound[i] ; 
            cpgline(2, xb, yb) ; 
        }
        cpgsls(1) ;		// retour aux lignes en trait plein
        cpgsci(1) ;		// couleur noire
    }

    cpgend() ; 
    
    delete [] xx ; 

}



//******************************************************************************
//      Multiple profiles, multiple device, uniform sampling 
//******************************************************************************

void des_profile_mult(const float* uutab, int nprof, int nx,
            float xmin, float xmax, const char* nomx, const char* nomy, 
	    const char* title, const int* line_style, 
            int ngraph, bool closeit, const char* device,
            int nbound, float* xbound) {

    const int ngraph_max = 100 ; 
    static int graph_list[ngraph_max] ; 
    static bool first_call = true ; 
    
    // First call operations
    // ---------------------
        
    if (first_call) {       // initialization of all the graphic devices to 0 :
        for (int i=0; i<ngraph_max; i++) {
            graph_list[i] = 0 ; 
        } 
        first_call = false ; 
    }
		          
         
    // Search for the extremal values of the field : 
    // -------------------------------------------

    int ntot = nprof * nx ; 
    float uumin = uutab[0] ;
    float uumax = uutab[0] ;
    for (int i=1; i<ntot; i++) {
	    if (uutab[i] < uumin) uumin = uutab[i] ;
	    if (uutab[i] > uumax) uumax = uutab[i] ;
    }

    cout << "  " << nomy << " : min, max : " << uumin << "   " << uumax 
         << endl ; 

    // Points abscisses : 
    // ----------------
    
    float* xx = new float[nx] ; 
    float hx = (xmax-xmin)/float(nx-1) ;
    for(int i=0; i<nx; i++) {
	xx[i] = xmin + float(i) * hx ; 
    }
         
    // Graphics display
    // ----------------
    
    // Opening of the device
    
    if ( (ngraph < 0) || (ngraph >= ngraph_max) ) {
        cerr << "des_profile_mult : graph index out of range ! \n" ;
        cerr << " ngraph = " << ngraph << "  while range = 0, " 
            << ngraph_max-1 << endl ; 
        abort() ;
    }
    
    if (graph_list[ngraph] == 0) { // opening is required
                                   // -------------------
    
        if (device == 0x0) device = "?" ; 
   
        graph_list[ngraph] = cpgopen(device) ; 
        
        if ( graph_list[ngraph] <= 0 ) {
	        cerr << "des_profile_mult: problem in opening PGPLOT display !\n" ;
            abort() ; 
        }
        
        cpgask(0) ;  // Disables the ``Type RETURN for next page:'' prompt
        
    }
    else {   // the graphic device has been opened previously   

        cpgslct( graph_list[ngraph] ) ; // selects the appropriate device
    }
    
    // Drawing
    // -------
     
    // Taille des caracteres:
    float size = float(1.3) ;
    cpgsch(size) ;
    
    // Epaisseur des traits:
    int lepais = 1 ; 
    cpgslw(lepais) ;
    
    // Fonte axes: caracteres romains:
    cpgscf(2) ;

    // Cadre de la figure
    float uuamp = uumax - uumin ; 
    float uumin1 = uumin - float(0.05) * uuamp ; 
    float uumax1 = uumax + float(0.05) * uuamp ; 
    cpgsls(1) ;  
    cpgenv(xmin, xmax, uumin1, uumax1, 0, 0 ) ; 
    cpglab(nomx,nomy,title) ;
     
    
    for (int i=0; i<nprof; i++) {
	const float* uudes = uutab + i*nx ;
        
        if (line_style == 0x0) cpgsls(i%5 + 1) ; 
        else cpgsls(line_style[i]) ;  
    	
        cpgline(nx, xx, uudes) ; 
    }
	
    // Plot of domain boundaries
    // -------------------------
    
    if (nbound > 0) {
        float xb[2] ; 
        float yb[2] ; 
        yb[0] = uumin1 ; 
        yb[1] = uumax1 ; 
        cpgsls(3) ;		// lignes en trait mixte
        cpgsci(3) ;		// couleur verte
        for (int i=0; i<nbound; i++) {
            xb[0] = xbound[i] ; 
            xb[1] = xbound[i] ; 
            cpgline(2, xb, yb) ; 
        }
        cpgsls(1) ;		// retour aux lignes en trait plein
        cpgsci(1) ;		// couleur noire
    }


    if (closeit) {
        cpgclos() ; 
        graph_list[ngraph] = 0 ; 
    }
    
    delete [] xx ; 

}


//******************************************************************************
//      Single profile, single device, arbitrary sampling 
//******************************************************************************

void des_profile(const float* uutab, int nx, const float *xtab, 
		 const char* nomx, const char* nomy, 
                 const char* title, const char* device,
                 int nbound, float* xbound) {
		 
    // Search for the extremal values of the field : 
    // -------------------------------------------

    float uumin = uutab[0] ;
    float uumax = uutab[0] ;
    float xmin = xtab[0] ;
    float xmax = xtab[0] ;
    for (int i=1; i<nx; i++) {
	uumin = (uutab[i] < uumin) ? uutab[i] : uumin ;
	uumax = (uutab[i] > uumax) ? uutab[i] : uumax ;	
	xmin = (xtab[i] < xmin) ? xtab[i] : xmin ;
	xmax = (xtab[i] > xmax) ? xtab[i] : xmax ;	
    }

    cout << "  " << nomy << " : min, max : " << uumin << "   " << uumax << endl;
    cout << "  " << "domain: " << "min, max : " << xmin << "   " << xmax << endl ; 

    // Points abscisses : 
    // ----------------
/*    
    float* xx = new float[nx] ; 
    float hx = (xmax-xmin)/float(nx-1) ;
    for(int i=0; i<nx; i++) {
	xx[i] = xmin + float(i) * hx ; 
    }
*/         
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
    float size = float(1.3) ;
    cpgsch(size) ;
    
    // Epaisseur des traits:
    int lepais = 1 ; 
    cpgslw(lepais) ;
    
    // Fonte axes: caracteres romains:
    cpgscf(2) ;

    // Cadre de la figure
    float uuamp = uumax - uumin ; 
    float uumin1 = uumin - float(0.05) * uuamp ; 
    float uumax1 = uumax + float(0.05) * uuamp ; 
    cpgenv(xmin, xmax, uumin1, uumax1, 0, 0 ) ; 
    cpglab(nomx,nomy,title) ;
    
    // Drawing of curve 
    cpgline(nx, xtab, uutab) ; 
    
    
    // Plot of domain boundaries
    // -------------------------
    
    if (nbound > 0) {
        float xb[2] ; 
        float yb[2] ; 
        yb[0] = uumin1 ; 
        yb[1] = uumax1 ; 
        cpgsls(3) ;		// lignes en trait mixte
        cpgsci(3) ;		// couleur verte
        for (int i=0; i<nbound; i++) {
            xb[0] = xbound[i] ; 
            xb[1] = xbound[i] ; 
            cpgline(2, xb, yb) ; 
        }
        cpgsls(1) ;		// retour aux lignes en trait plein
        cpgsci(1) ;		// couleur noire
    }

    cpgend() ; 

}



//******************************************************************************
//      Multiple profiles, multiple device, arbitrary sampling 
//******************************************************************************

void des_profile_mult(const float* uutab, int nprof, int nx, const float* xtab, 
            const char* nomx, const char* nomy, const char* title, 
            const int* line_style, int ngraph, bool closeit,
            const char* device, int nbound, float* xbound) {

    const int ngraph_max = 100 ; 
    static int graph_list[ngraph_max] ; 
    static bool first_call = true ; 
    
    // First call operations
    // ---------------------
        
    if (first_call) {       // initialization of all the graphic devices to 0 :
        for (int i=0; i<ngraph_max; i++) {
            graph_list[i] = 0 ; 
        } 
        first_call = false ; 
    }
		          
         
    // Search for the extremal values of x and of the field : 
    // ----------------------------------------------------

    int ntot = nprof * nx ; 
    float uumin = uutab[0] ;
    float uumax = uutab[0] ;
    for (int i=1; i<ntot; i++) {
	    if (uutab[i] < uumin) uumin = uutab[i] ;
	    if (uutab[i] > uumax) uumax = uutab[i] ;
    }

    float xmin = xtab[0] ;
    float xmax = xtab[0] ;
    for (int i=1; i<ntot; i++) {
	    if (xtab[i] < xmin) xmin = xtab[i] ;
	    if (xtab[i] > xmax) xmax = xtab[i] ;
    }

    cout << "  " << nomy << " : min, max : " << uumin << "   " << uumax 
         << endl ; 

         
    // Graphics display
    // ----------------
    
    // Opening of the device
    
    if ( (ngraph < 0) || (ngraph >= ngraph_max) ) {
        cerr << "des_profile_mult : graph index out of range ! \n" ;
        cerr << " ngraph = " << ngraph << "  while range = 0, " 
            << ngraph_max-1 << endl ; 
        abort() ;
    }
    
    if (graph_list[ngraph] == 0) { // opening is required
                                   // -------------------
    
        if (device == 0x0) device = "?" ; 
   
        graph_list[ngraph] = cpgopen(device) ; 
        
        if ( graph_list[ngraph] <= 0 ) {
	        cerr << "des_profile_mult: problem in opening PGPLOT display !\n" ;
            abort() ; 
        }
        
        cpgask(0) ;  // Disables the ``Type RETURN for next page:'' prompt
        
    }
    else {   // the graphic device has been opened previously   

        cpgslct( graph_list[ngraph] ) ; // selects the appropriate device
    }
    
    // Drawing
    // -------
     
    // Taille des caracteres:
    float size = float(1.3) ;
    cpgsch(size) ;
    
    // Epaisseur des traits:
    int lepais = 1 ; 
    cpgslw(lepais) ;
    
    // Fonte axes: caracteres romains:
    cpgscf(2) ;

    // Draw the figure
    float uuamp = uumax - uumin ; 
    float uumin1 = uumin - float(0.05) * uuamp ; 
    float uumax1 = uumax + float(0.05) * uuamp ; 
    cpgsls(1) ;  
    cpgenv(xmin, xmax, uumin1, uumax1, 0, 0 ) ; 
    cpglab(nomx,nomy,title) ;
     
    
    for (int i=0; i<nprof; i++) {
	const float* uudes = uutab + i*nx ;
	const float* xdes = xtab + i*nx ;
        
        if (line_style == 0x0) cpgsls(i%5 + 1) ; 
        else cpgsls(line_style[i]) ;  
    	
        cpgline(nx, xdes, uudes) ; 
    }
	
    // Plot of domain boundaries
    // -------------------------
    
    if (nbound > 0) {
        float xb[2] ; 
        float yb[2] ; 
        yb[0] = uumin1 ; 
        yb[1] = uumax1 ; 
        cpgsls(3) ;		// lignes en trait mixte
        cpgsci(3) ;		// couleur verte
        for (int i=0; i<nbound; i++) {
            xb[0] = xbound[i] ; 
            xb[1] = xbound[i] ; 
            cpgline(2, xb, yb) ; 
        }
        cpgsls(1) ;		// retour aux lignes en trait plein
        cpgsci(1) ;		// couleur noire
    }

    if (closeit) {
        cpgclos() ; 
        graph_list[ngraph] = 0 ; 
    }
    
}

}

