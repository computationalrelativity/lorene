/*
 * Comparison of various interpolation of a function
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
 * Revision 1.2  2005/11/14 14:11:22  e_gourgoulhon
 * log(int ) -> log(double)
 * + suppressed plot_close.
 *
 * Revision 1.1  2005/11/14 01:56:59  e_gourgoulhon
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
#include "plot.h"

double ff(double x) ; 

int main() {

    // For the plots
    int ng = 1000 ;      // Number of points to draw the curves
    double* xg = new double[ng] ; 
    double* yg = new double[ng] ; 
    for (int j=0; j<ng; j++) xg[j] = -1. + 2. * double(j) / double(ng-1) ; 
    
    int nfig = -1 ; // Figure index

    ofstream funif("err_unif.d") ;
    ofstream fleg("err_leg.d") ;
    ofstream fcheb("err_cheb.d") ;
    
    int nn = 0 ; 
    while (nn >= 0) {
        cout << "Maximum degree N of polynomials ? (negative value = exit) " 
            << endl ; 
        cin >> nn ; 
        while ( cin.get()!='\n' ) ;
        if (nn < 0) continue ; 
        
        double* xc = new double[nn+1] ;     // to draw points
        double* yc = new double[nn+1] ;     

        // The three types of grid:
        
        Grid_uniform x_unif(nn+1) ; 
        Grid_Legendre_GL x_leg(nn+1) ; 
        Grid_Chebyshev_Gauss x_cheb(nn+1) ;
        
        cout << "Lebesgue constants: " << endl ; 
        cout << "   uniform grid : " << x_unif.lebesgue_constant() << endl ; 
        cout << "   Legendre GL grid : " << x_leg.lebesgue_constant() << endl ; 
        cout << "   Chebyshev Gauss grid : " << x_cheb.lebesgue_constant() 
            << " <-> theor. lower / upper bound: " 
            << 2./M_PI * log(double(nn+1)) + 0.9625 << " / " 
            << 2./M_PI * log(double(nn+1)) + 1 << endl ; 
        
        // Interpolation through the uniform grid :
        
        nfig++ ; 
        for (int j=0; j<ng; j++) yg[j] = ff(xg[j]) ; 
        plot_profile(yg, ng, 2, 1, nfig, -1.2, 1.2) ;     
        
        double err = 0 ; 
        for (int j=0; j<ng; j++) {
            double y = x_unif.interpole(ff, xg[j]) ; 
            double diff = fabs(y - yg[j]) ; 
            if (diff > err) err = diff ; 
            yg[j] = y ; 
        }
        cout << "Error uniform interpolation: " << err << endl ; 
        funif << nn << "  " << err << endl ; 
        plot_profile(yg, ng, 3, 1, nfig) ;     
        
        for (int i=0; i<=nn; i++) {
            xc[i] = x_unif(i) ; 
            yc[i] = x_unif.interpole(ff, xc[i]) ; 
        }
        plot_point_set(nn+1, xc, yc, 3, nfig) ; 
        
        // plot_close(nfig) ; // closing required if EPS figure 
        
        // Interpolation through the Legendre GL grid :
        
        nfig++ ; 
        for (int j=0; j<ng; j++) yg[j] = ff(xg[j]) ; 
        plot_profile(yg, ng, 2, 1, nfig, -0.2, 1.2) ;     
        
        err = 0 ; 
        for (int j=0; j<ng; j++) {
            double y = x_leg.interpole(ff, xg[j]) ; 
            double diff = fabs(y - yg[j]) ; 
            if (diff > err) err = diff ; 
            yg[j] = y ; 
        }
        cout << "Error Legendre GL interpolation: " << err << endl ; 
        fleg << nn << "  " << err << endl ; 
        plot_profile(yg, ng, 3, 1, nfig) ;     
        
        for (int i=0; i<=nn; i++) {
            xc[i] = x_leg(i) ; 
            yc[i] = x_leg.interpole(ff, xc[i]) ; 
        }
        plot_point_set(nn+1, xc, yc, 3, nfig) ; 
        
        // plot_close(nfig) ; // closing required if EPS figure 

        // Interpolation through the Chebyshev Gauss grid :
        
        nfig++ ; 
        for (int j=0; j<ng; j++) yg[j] = ff(xg[j]) ; 
        plot_profile(yg, ng, 2, 1, nfig, -0.2, 1.2) ;     
        
        err = 0 ; 
        for (int j=0; j<ng; j++) {
            double y = x_cheb.interpole(ff, xg[j]) ; 
            double diff = fabs(y - yg[j]) ; 
            if (diff > err) err = diff ; 
            yg[j] = y ; 
        }
        cout << "Error Chebyshev Gauss interpolation: " << err << endl ; 
        fcheb << nn << "  " << err << endl ; 
        plot_profile(yg, ng, 3, 1, nfig) ;     
        
        for (int i=0; i<=nn; i++) {
            xc[i] = x_cheb(i) ; 
            yc[i] = x_cheb.interpole(ff, xc[i]) ; 
        }
        plot_point_set(nn+1, xc, yc, 3, nfig) ; 
        
        // plot_close(nfig) ; // closing required if EPS figure 

        delete [] xc ; 
        delete [] yc ; 
       
    }
    
    plot_close_all() ; 
    
    funif.close() ; 
    fleg.close() ; 
    fcheb.close() ; 

    delete [] xg ; 
    delete [] yg ; 

    return EXIT_SUCCESS ; 
}

double ff(double x) {
//    return 1. / (1. + 25.*x*x) ; 
    return 1. / (1. + 16.*x*x) ; 
//    return cos(exp(2*x)) ; 
}
