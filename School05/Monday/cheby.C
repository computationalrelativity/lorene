/*
 * Code for Chebyshev expansion
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
#include <fstream>

using namespace std ;

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ortho_poly.h"
#include "grid.h"
#include "plot.h"

void petite_pause() ; 
double ff(double x) ; 

int main() {

    int nn = 4 ; 
    
    Chebyshev_poly cheb(nn) ;
    
    cout << cheb << endl ; 
    cout << cheb.gauss_nodes() << endl ; 
    cout << cheb.gauss_lobatto_nodes() << endl ; 

    // For the plots
    int ng = 1000 ;      // Number of points to draw the curves
    double* xg = new double[ng] ; 
    double* yg = new double[ng] ; 
    for (int j=0; j<ng; j++) xg[j] = -1. + 2. * double(j) / double(ng-1) ; 
    
    double* xc = new double[nn+1] ;     // to draw points
    double* yc = new double[nn+1] ;     
    
    int nfig = 0 ;  // Figure index
    
    for (int i=0; i<=nn; i++) {
        
        for (int j=0; j<ng; j++) yg[j] = cheb(i, xg[j]) ; 
            plot_profile(yg, ng, i%15 + 1, 1, nfig, -1.1, 1.1, 
                     "Chebyshev polynomials up to N=8") ;     

        // petite_pause() ; 

    }

    Chebyshev_poly cheb1(nn+1) ;
    nfig++ ; 
    for (int j=0; j<ng; j++) yg[j] = cheb1(nn+1, xg[j]) ; 
    plot_profile(yg, ng, 2, 1, nfig, -1.1, 1.1, "Chebyshev Gauss nodes") ;     

    cheb.gauss_nodes().plot(3, nfig) ; 
    // cheb.gauss_lobatto_nodes().plot(4, nfig) ; 
    
    double* cf = new double[nn+1] ; 
    double* cf1 = new double[nn+1] ; 
    
    cheb.coef_interpolant_GL(ff, cf) ; 
    cheb.coef_interpolant_GL_FFT(ff, cf1) ; 
    
    cout << "Coefficients of f interpolant (Gauss-Lobatto nodes): " << endl ; 
    for (int i=0; i<=nn; i++) {
        cout << i << " :  " << cf[i] << "  " << cf1[i] 
        << "  " << cf[i] - cf1[i] << endl ;  
    }

    double* cf2 = new double[nn+1] ; 
    cheb.coef_interpolant_Gauss(ff, cf2) ; 
    cout << endl <<
     "Coefficients of f interpolant (Gauss-Lobatto versus Gauss nodes)" 
     << endl ; 
    for (int i=0; i<=nn; i++) {
        cout << i << " :  " << cf[i] << "  " << cf2[i] 
        << "  " << cf[i] - cf2[i] << endl ;  
    }

    nfig++ ; 
    for (int j=0; j<ng; j++) yg[j] = ff(xg[j]) ; 
    plot_profile(yg, ng, 2, 1, nfig, -1.1, 1.1) ;     

    for (int j=0; j<ng; j++) yg[j] = cheb.series(cf, xg[j]) ; 
    plot_profile(yg, ng, 3, 1, nfig) ;     
    cheb.gauss_lobatto_nodes().plot(3, nfig) ; 

    for (int i=0; i<=nn; i++) {
        xc[i] = cheb.gauss_lobatto_nodes()(i) ; 
        yc[i] = cheb.series(cf, xc[i]) ; 
    }
    plot_point_set(nn+1, xc, yc, 3, nfig) ; 
    
    for (int j=0; j<ng; j++) yg[j] = cheb.series(cf2, xg[j]) ; 
    //plot_profile(yg, ng, 7, 1, nfig) ;     
    //cheb.gauss_nodes().plot(7, nfig) ; 

    double* cf_proj = new double[nn+1] ; 
    cheb.coef_projection(ff, cf_proj) ; 
    for (int j=0; j<ng; j++) yg[j] = cheb.series(cf_proj, xg[j]) ; 
    plot_profile(yg, ng, 5, 1, nfig) ;     
    
    ofstream fcoef("coef_proj.d") ;
    for (int i=0; i<=nn; i++) 
        fcoef << i << "  " << fabs(cf_proj[i]) << endl ; 
    fcoef.close() ;    
   
    
    // Check: comparison between Grid::interpole and Ortho_poly::series
    cout << endl <<  
    "Maximum difference between Grid::interpole and Ortho_poly::series :"
    << endl ; 
    double diff = 0 ;
    for (int j=0; j<ng; j++) {
        double diff0 = fabs( cheb.series(cf, xg[j]) 
            - cheb.gauss_lobatto_nodes().interpole(ff, xg[j]) ) ; 
        if (diff0 > diff) diff = diff0 ; 
    }    
    cout << "  Gauss-Lobatto case : " << diff << endl ; 

    diff = 0 ;
    for (int j=0; j<ng; j++) {
        double diff0 = fabs( cheb.series(cf2, xg[j]) 
            - cheb.gauss_nodes().interpole(ff, xg[j]) ) ; 
        if (diff0 > diff) diff = diff0 ; 
    }    
    cout << "  Gauss case : " << diff << endl ; 

    petite_pause() ; 

    plot_close_all() ; 
    
    delete [] cf ; 
    delete [] cf1 ; 
    delete [] cf2 ; 
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

//    return 1. / (1. + 25.*x*x) ; 

    return cos(exp(2*x)) ; 
    
//    if (x == double(0)) return 0.;
//    else return exp(-1. / (25*x*x)) ; 
}
