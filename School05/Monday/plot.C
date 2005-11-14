/*
 * Definition of elementary graphical functions
 */

/*
 *   Copyright (c) 2005 Eric Gourgoulhon
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
 * Revision 1.1  2005/11/14 01:57:00  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

#include <iostream>

using namespace std ;

#include <stdlib.h>

#include <cpgplot.h>

#include "plot.h"


namespace Plotid {
    const int nfig_max = 100 ; 
    int fig_id[nfig_max] ;     
}

//===========================================================================

void plot_init() {
    using namespace Plotid ; 

    static bool never_called = true ; 
    
    if (never_called) {
        for (int i=0; i<nfig_max; i++) {
            fig_id[i] = 0 ; 
        }     
        never_called = false ; 
    }
    
}

//===========================================================================

void plot_open(int nfig, double ymin, double ymax, const char* title,    
    const char* label_y, const char* device) {
    using namespace Plotid ; 

    assert( fig_id[nfig] == 0 ) ; 
    
    if (device == 0x0) device = "?" ; 
   
    fig_id[nfig] = cpgopen(device) ; 
        
    if ( fig_id[nfig] <= 0 ) {
        cerr << "plot_profile: problem in opening PGPLOT display !\n" ;
        abort() ; 
    }
        
    cpgask(0) ;  // Disables the ``Type RETURN for next page:'' prompt
            
    cpgsch(1.3) ;  // Character size
    
    cpgslw(4) ;  // Line width
    
    cpgscf(2) ;  // Axis fonts = Roman

    cpgsls(1) ;  // Continous lines

    // Sets the window :         
    cpgenv(-1., 1., float(ymin), float(ymax), 0, 0 ) ; 
    
    // Title and labels
    if (label_y == 0x0) label_y = "y" ; 
    if (title == 0x0) title = " " ; 
    cpglab("x", label_y, title) ;
    
}



//===========================================================================

void plot_close(int nfig) {

    using namespace Plotid ;     
    plot_init() ; 

    if ( (nfig < 0) || (nfig >= nfig_max) ) {
        cerr << "plot_close : figure index out of range ! \n" ;
        cerr << " nfig = " << nfig << "  while range = 0, " 
            << nfig_max-1 << endl ; 
        abort() ;
    }
    
    if (fig_id[nfig] != 0) { 
        cpgslct( fig_id[nfig] ) ; // selects the appropriate device        
        cpgclos() ; 
        fig_id[nfig] = 0 ; 
    }

}

//===========================================================================

void plot_close_all() {

    using namespace Plotid ;     
    plot_init() ; 

    for (int nfig = 0; nfig < nfig_max; nfig++) {
        if (fig_id[nfig] != 0) { 
            cpgslct( fig_id[nfig] ) ; // selects the appropriate device        
            cpgclos() ; 
            fig_id[nfig] = 0 ; 
        }
    }
}

//===========================================================================

void plot_profile(const double* yy, int nx, int color, int style,    
    int nfig, double ymin, double ymax, const char* title,    
    const char* label_y, const char* device) {

    using namespace Plotid ; 
    plot_init() ; 
		          
    // Conversion double -> float
    // --------------------------
    float* yyf = new float[nx] ; 
    for (int i=0; i<nx; i++) {
        yyf[i] = yy[i] ; 
    }
         
    // Points abscisses : 
    // ----------------
    float xmin = -1. ; 
    float xmax = 1. ; 
    float* xx = new float[nx] ; 
    float hx = (xmax-xmin)/float(nx-1) ;
    for(int i=0; i<nx; i++) {
	    xx[i] = xmin + i * hx ; 
    }
         
    // Graphics display
    // ----------------
    
    // Opening of the device
    
    if ( (nfig < 0) || (nfig >= nfig_max) ) {
        cerr << "plot_profile : graph index out of range ! \n" ;
        cerr << " nfig = " << nfig << "  while range = 0, " 
            << nfig_max-1 << endl ; 
        abort() ;
    }
    
    if (fig_id[nfig] == 0) { // opening is required
                               
        plot_open(nfig, ymin, ymax, title, label_y, device ) ; 

    }
    else {   // the graphic device has been opened previously   

        cpgslct( fig_id[nfig] ) ; // selects the appropriate device
    }  
    
    // Drawing of the curve
    cpgsci(color) ; 
    cpgsls(style) ;  
    cpgline(nx, xx, yyf) ; 	

    delete [] xx ; 
    delete [] yyf ; 

}

//===========================================================================

void plot_point(double x, double y, int color, int nfig, double ymin, 
    double ymax, const char* title, const char* label_y, const char* device) {

    using namespace Plotid ; 
    plot_init() ; 

    if ( (nfig < 0) || (nfig >= nfig_max) ) {
        cerr << "plot_point : figure index out of range ! \n" ;
        cerr << " nfig = " << nfig << "  while range = 0, " 
            << nfig_max-1 << endl ; 
        abort() ;
    }
    
    if (fig_id[nfig] == 0) { // opening is required
                               
        plot_open(nfig, ymin, ymax, title, label_y, device ) ; 

    }
    else {   // the graphic device has been opened previously   

        cpgslct( fig_id[nfig] ) ; // selects the appropriate device
    }

    cpgsci(color) ; 

    cpgpt1(float(x), float(y), 24) ; 
    
}

//===========================================================================

void plot_point_set(int np, const double* xx, const double* yy, int color, 
                    int nfig, double ymin, double ymax, const char* title, 
                    const char* label_y, const char* device) {

    using namespace Plotid ; 
    plot_init() ; 

    if ( (nfig < 0) || (nfig >= nfig_max) ) {
        cerr << "plot_point : figure index out of range ! \n" ;
        cerr << " nfig = " << nfig << "  while range = 0, " 
            << nfig_max-1 << endl ; 
        abort() ;
    }
    
    if (fig_id[nfig] == 0) { // opening is required
                               
        plot_open(nfig, ymin, ymax, title, label_y, device ) ; 

    }
    else {   // the graphic device has been opened previously   

        cpgslct( fig_id[nfig] ) ; // selects the appropriate device
    }

    cpgsci(color) ; 
    
    float* x1 = new float[np] ; 
    float* y1 = new float[np] ; 
    
    for (int i=0; i<np; i++) {
        x1[i] = xx[i] ;
        y1[i] = yy[i] ;
    }
    cpgpt(np, x1, y1, 24) ; 

    delete [] x1 ; 
    delete [] y1 ; 
}
