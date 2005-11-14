/*
 * Test code for class Grid
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
 * Revision 1.1  2005/11/14 01:56:58  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */


#include <iostream>

using namespace std ;

#include <stdlib.h>
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
    
    // Nodes
    int nb_nodes = 1 ; 
    while (nb_nodes > 0) {
        cout << "Number of nodes ? (negative value = exit) " << endl ; 
        cin >> nb_nodes ; 
        while ( cin.get()!='\n' ) ;
        if (nb_nodes <= 0) continue ; 
        
        int nn = nb_nodes - 1 ; 

    double* xx = new double[nb_nodes] ;
    for (int i=0; i<=nn; i++) xx[i] = -1 + 2 * double(i) / double(nn) ; 
       
    Grid xcoloc(nb_nodes,xx) ;    // Construction of the grid 
    
    cout << "Grid: " << xcoloc << endl ; 
    cout << "Lebesgue constant : " << xcoloc.lebesgue_constant() << endl ; 
    int nfig = 0 ;  // index of next figure
    // xcoloc.plot(3, nfig) ;  // draws the point in 3=green
    
    // Lagrange polynomials 
    
    for (int i=0; i<=nn; i++) {
        for (int j=0; j<ng; j++) yg[j] = xcoloc.lagrange(i, xg[j]) ; 
        // nfig++ ;  // new figure
        plot_profile(yg, ng, i%15+2, 1, nfig, -3., 3., "Lagrange polynomials") ; 
        xcoloc.plot(1, nfig) ;
        plot_point(xcoloc(i), 0., i%15+2, nfig) ; 
        double xp = xcoloc(i) ;
        double yp = xcoloc.lagrange(i, xp) ;
        plot_point_set(1, &xp, &yp, i%15+2, nfig) ; 
        // plot_close(nfig) ;   // Figure needs to be closed before the
                                // next one if it's an EPS file. 
    }
         
    // Nodal polynomial
    
    for (int j=0; j<ng; j++) yg[j] = xcoloc.nodal_polynomial(xg[j]) ; 

    nfig++ ;   // new figure
    plot_profile(yg, ng, 2, 1, nfig, -0.1, 0.1, "Nodal polynomial") ; 
    xcoloc.plot(2, nfig) ; 

    // Interpolating polynomial 

    for (int j=0; j<ng; j++) {
        yg[j] = ff(xg[j]) ; 
    }
    nfig++ ;   // new figure
 //   plot_profile(yg, ng, 3, 1, nfig, -0.2, 1.1, "Interpolation of 1/(1+25x\\u2\\d)") ;     
    plot_profile(yg, ng, 3, 1, nfig, -1.2, 1.2, "Interpolation of cos(2 exp(x))") ;     

    for (int j=0; j<ng; j++) {
        yg[j] = xcoloc.interpole(ff, xg[j]) ; 
    }
    plot_profile(yg, ng, 2, 1, nfig) ; 
    xcoloc.plot(2, nfig) ; 
    
    petite_pause() ; 

    plot_close_all() ; 
    
    delete [] xx ; 

    }  // end of while

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
//    return 1. / (1. + 25.*x*x) ; 
    return cos(exp(2*x)) ; 
}
