/*
 * Code to exhibit Runge's phenomenon 
 *
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
 * Revision 1.2  2005/11/14 14:17:02  e_gourgoulhon
 * Added #include <stdio.h>
 *
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
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "grid.h"
#include "plot.h"

void petite_pause() ; 
double ff(double x) ; 

int main() {

    // For the plots
    int ng = 1000 ;      // Number of points to draw the curves
    double* xg = new double[ng] ; 
    double* yg = new double[ng] ; 
    for (int j=0; j<ng; j++) xg[j] = -1. + 2. * double(j) / double(ng-1) ; 
    
    
    for (int j=0; j<ng; j++) yg[j] = ff(xg[j]) ; 
    
    int nfig0 = 0 ;  // Figure index
    plot_profile(yg, ng, 1, 1, nfig0, -0.5, 1.5, "Interpolation of 1/(1+25x\\u2\\d)") ;     

    // Loop on the number of nodes
    
    for (int k = 0; k<6; k++) {
    
        int nb_nodes = 4*(k+1) + 1 ; 
        int nn = nb_nodes - 1 ; 

        int nfig1 = k+1 ;  // Figure index
               
        Grid_uniform xcoloc(nb_nodes) ;    // Construction of the grid 
        
        // Grid_Chebyshev_GL xcoloc(nb_nodes) ;  // other example (free of Runge phenomenon) 
    
        cout << "Grid: " << xcoloc << endl ; 
        
        // Interpolating polynomial 

        for (int j=0; j<ng; j++) yg[j] = ff(xg[j]) ; 
    
        char title[180] ; 
        strcpy(title, "Interpolation of 1/(1+25x\\u2\\d) , N = ") ; 
        char string_nn[4] ;
        sprintf(string_nn, "%d", nn) ; 
        strcat(title, string_nn) ; 
        plot_profile(yg, ng, 1, 1, nfig1, -0.5, 1.5, title) ;     

        for (int j=0; j<ng; j++) yg[j] = xcoloc.interpole(ff, xg[j]) ; 

        plot_profile(yg, ng, k%15+2, 2, nfig0) ; 
        plot_profile(yg, ng, 2, 1, nfig1) ; 

        xcoloc.plot(2, nfig1) ;   

    }  // end of loop on the number of nodes
    
    petite_pause() ; 

    plot_close_all() ; 
    
    delete [] xg ; 
    delete [] yg ; 

    return EXIT_SUCCESS ; 
}

void petite_pause() {
    cout.flush() ;
    cout << "Continue = 'return'" << endl ;
    char cret ; 
    cin.get(cret) ;
}


double ff(double x) {
//    return 2*x*x - 1. ; 
//     return 1./ (cos(6*x) + 2.) ;  // other example
    return 1. / (1. + 25.*x*x) ; 
//    if (x == double(0)) return 0.;
//    else return exp(-1. / (25*x*x)) ; 
}

