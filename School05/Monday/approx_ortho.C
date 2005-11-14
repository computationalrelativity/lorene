/*
 * Approximation of a function by orthogonal polynomials
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
#include <math.h>

#include "grid.h"
#include "ortho_poly.h"
#include "plot.h"

double ff(double x) ; 

int main() {

    // For the plots
    int ng = 1000 ;      // Number of points to draw the curves
    double* xg = new double[ng] ; 
    double* yg = new double[ng] ; 
    double* yg0 = new double[ng] ; 
    double* yg1 = new double[ng] ; 
    
    for (int j=0; j<ng; j++) xg[j] = -1. + 2. * double(j) / double(ng-1) ; 
    
    int nfig = -1 ; // Figure index

    ofstream fproj("err_proj.d") ;
    ofstream finter("err_inter.d") ;
    
    int nn = 0 ; 
    while (nn >= 0) {
        cout << "Maximum degree N of polynomials ? (negative value = exit) " 
            << endl ; 
        cin >> nn ; 
        while ( cin.get()!='\n' ) ;
        if (nn < 0) continue ; 
        
        double* xc = new double[nn+1] ;     // to draw points
        double* yc = new double[nn+1] ;     

        // Plot of f
        nfig++ ; 
        for (int j=0; j<ng; j++) yg0[j] = ff(xg[j]) ; 
        plot_profile(yg0, ng, 2, 1, nfig, -1.2, 1.2) ;     

        // Basis of orthogonal polynomials:
        
        Chebyshev_poly poly(nn) ;

        // Coefficients of the projection of f onto P_N:

        double* cf_proj = new double[nn+1] ; 
        poly.coef_projection(ff, cf_proj) ; 

        double err = 0 ; 
        for (int j=0; j<ng; j++) {
            double y = poly.series(cf_proj, xg[j]) ; 
            double diff = fabs(y - yg0[j]) ; 
            if (diff > err) err = diff ; 
            yg1[j] = y ; 
        }
        cout << "Error projection : " << err << endl ; 
        fproj << nn << "  " << err << endl ; 
        plot_profile(yg1, ng, 5, 1, nfig) ;     
        
        // Coefficients of the Gauss-Lobatto interpolant of f
        double* cf_inter = new double[nn+1] ; 
        poly.coef_interpolant_GL(ff, cf_inter) ; 

        err = 0 ; 
        for (int j=0; j<ng; j++) {
            double y = poly.series(cf_inter, xg[j]) ; 
            double diff = fabs(y - yg0[j]) ; 
            if (diff > err) err = diff ; 
            yg[j] = y ; 
        }
        cout << "Error GL interpolant : " << err << endl ; 
        finter << nn << "  " << err << endl ; 
        plot_profile(yg, ng, 3, 1, nfig) ;     
        
        for (int i=0; i<=nn; i++) {
            xc[i] = poly.gauss_lobatto_nodes()(i) ; 
            yc[i] = poly.gauss_lobatto_nodes().interpole(ff, xc[i]) ; 
        }
        plot_point_set(nn+1, xc, yc, 3, nfig) ; 
        
        plot_close(nfig) ; // closing required if EPS figure 

        // Aliasing error

        // for (int j=0; j<ng; j++) yg[j] -= yg1[j] ; 
        // nfig++ ; 
        // plot_profile(yg, ng, 4, 1, nfig, -0.2, 0.2, "Aliasing error") ;     
                
        
        // plot_close(nfig) ; // closing required if EPS figure 

        delete [] xc ; 
        delete [] yc ; 
        delete [] cf_proj ; 
        delete [] cf_inter ; 
       
    }
    
    plot_close_all() ; 
    
    fproj.close() ; 
    finter.close() ; 

    delete [] xg ; 
    delete [] yg ; 
    delete [] yg0 ; 
    delete [] yg1 ; 

    return EXIT_SUCCESS ; 
}



double ff(double x) {
//    return 1. / (1. + 25.*x*x) ; 
//    return 1. / (1. + 16.*x*x) ; 
    return cos(exp(2*x)) ; 
}

