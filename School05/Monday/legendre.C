/*
 * Code for Legendre expansion
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
 * Revision 1.2  2014/10/06 15:09:48  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2005/11/14 01:56:59  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

#include <iostream>

using namespace std ;

#include <cstdlib>
#include <cstring>
#include <cmath>

#include "ortho_poly.h"
#include "grid.h"
#include "plot.h"

void petite_pause() ; 
double ff(double x) ; 

int main() {

    int nn = 8 ; 
    
    Legendre_poly leg(nn) ;
    
    cout << leg << endl ; 
    // cout << leg.gauss_nodes() << endl ; 
    cout << leg.gauss_lobatto_nodes() << endl ; 

    // For the plots
    int ng = 1000 ;      // Number of points to draw the curves
    double* xg = new double[ng] ; 
    double* yg = new double[ng] ; 
    for (int j=0; j<ng; j++) xg[j] = -1. + 2. * double(j) / double(ng-1) ; 
    
    double* xc = new double[nn+1] ;     // to draw points
    double* yc = new double[nn+1] ;     
    
    int nfig = 0 ;  // Figure index
    
    for (int i=0; i<=nn; i++) {
        
        for (int j=0; j<ng; j++) yg[j] = leg(i, xg[j]) ; 
        plot_profile(yg, ng, i%15 + 1, 1, nfig, -1.1, 1.1, 
            "Legendre polynomials up to N=8") ;     

        // petite_pause() ; 

    }

    
    double* cf = new double[nn+1] ; 
    
    leg.coef_interpolant_GL(ff, cf) ; 

    cout << endl <<
     "Coefficients of f interpolant (Gauss-Lobatto)" 
     << endl ; 
    for (int i=0; i<=nn; i++) {
        cout << i << " :  " << cf[i] << endl ;  
    }

    nfig++ ; 
    for (int j=0; j<ng; j++) yg[j] = ff(xg[j]) ; 
    plot_profile(yg, ng, 2, 1, nfig, -0.1, 1.1, "Approximation of a function") ;     

    for (int j=0; j<ng; j++) yg[j] = leg.series(cf, xg[j]) ; 
    plot_profile(yg, ng, 3, 1, nfig) ;     
    leg.gauss_lobatto_nodes().plot(3, nfig) ; 

    for (int i=0; i<=nn; i++) {
        xc[i] = leg.gauss_lobatto_nodes()(i) ; 
        yc[i] = leg.series(cf, xc[i]) ; 
    }
    plot_point_set(nn+1, xc, yc, 3, nfig) ; 
    
    // for (int j=0; j<ng; j++) yg[j] = leg.series(cf2, xg[j]) ; 
    // plot_profile(yg, ng, 7, 1, nfig) ;     
    // leg.gauss_nodes().plot(7, nfig) ; 

    double* cf_proj = new double[nn+1] ; 
    leg.coef_projection(ff, cf_proj) ; 
    for (int j=0; j<ng; j++) yg[j] = leg.series(cf_proj, xg[j]) ; 
    plot_profile(yg, ng, 5, 1, nfig) ;     
   

    // Check: comparison between Grid::interpole and Ortho_poly::series
    cout << endl <<  
    "Maximum difference between Grid::interpole and Ortho_poly::series :"
    << endl ; 
    double diff = 0 ;
    for (int j=0; j<ng; j++) {
        double diff0 = fabs( leg.series(cf, xg[j]) 
            - leg.gauss_lobatto_nodes().interpole(ff, xg[j]) ) ; 
        if (diff0 > diff) diff = diff0 ; 
    }    
    cout << "  Gauss-Lobatto case : " << diff << endl ; 

    // diff = 0 ;
    // for (int j=0; j<ng; j++) {
    //    double diff0 = fabs( leg.series(cf2, xg[j]) 
    //        - leg.gauss_nodes().interpole(ff, xg[j]) ) ; 
    //    if (diff0 > diff) diff = diff0 ; 
    // }    
    // cout << "  Gauss case : " << diff << endl ; 

    petite_pause() ; 

    plot_close_all() ; 
    
    delete [] cf ; 
   // delete [] cf2 ; 
    delete [] cf_proj ; 
    
    delete [] xc ; 
    delete [] yc ; 
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
