/*
 * Definition of class Chebyshev_poly
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
 * Revision 1.3  2014/10/06 15:09:47  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2005/11/14 14:12:09  e_gourgoulhon
 * Added include <assert.h>
 *
 * Revision 1.1  2005/11/14 01:56:58  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */


#include <iostream>

using namespace std ;

#include <cstdlib>
#include <cmath>
#include <cassert>

#include "ortho_poly.h"
#include "grid.h"

void coef_cheb_fft(int np, const double* ff, double* cf) ;

//----------------------//
// Standard constructor //
//----------------------//

Chebyshev_poly::Chebyshev_poly(int ni) : Ortho_poly(ni) {}

//--------------------//
//  Copy constructor  //
//--------------------//

Chebyshev_poly::Chebyshev_poly(const Chebyshev_poly& pp) : Ortho_poly(pp) {}


//--------------//
//  Destructor  //
//--------------//

Chebyshev_poly::~Chebyshev_poly() {}


//--------------//
//  Assignment  //
//--------------//

void Chebyshev_poly::operator=(const Chebyshev_poly& ) {

    cerr << "Chebyshev_poly::operator= not implemented !" << endl ; 
    abort() ; 
    
}

//-----------//
//  Display  //
//-----------//

ostream& operator<<(ostream& ost, const Chebyshev_poly& pp) {

    ost << "Basis of Chebyshev polynomials up to degree " << pp.n() << endl ; 
    
    return ost ; 
}

//-------------------//
//  Weight function  //
//-------------------//

double Chebyshev_poly::weight(double x) const {
    
    return 1. / sqrt( 1. - x*x) ; 
    
} 


//-------------------------------------------//
//  Evaluation of the Chebyshev polynomials  //
//-------------------------------------------//

double Chebyshev_poly::operator()(int i, double x) const {

    assert( i >= 0 ) ; 
    assert( i <= nn ) ; 
    
    if (i==0) return 1. ; 

    if (i==1) return x ; 
    
    double tjm2 = 1. ;  // value at step j - 2
    double tjm1 = x  ;  // value at step j - 1
    
    for (int j=2; j<=i; j++) {
	    double tj = 2*x* tjm1 - tjm2 ;
        tjm2 = tjm1 ; 
        tjm1 = tj ;       	    
    }
    
    return tjm1 ;     
    
}

//------------------------------------//
//  Computation of nodes and weights  //
//------------------------------------//


const Grid& Chebyshev_poly::gauss_nodes() const {
    
    if (p_gauss_nodes == 0x0) {  // the nodes must initialized

        p_gauss_nodes = new Grid_Chebyshev_Gauss(nn+1) ;
    }
    
    return *p_gauss_nodes ;
}


double Chebyshev_poly::gauss_weight(int i) const {

    if (p_gauss_weights == 0x0) {  // the weights must be computed

        p_gauss_weights = new double[nn+1] ;
        
        for (int j=0; j<=nn; j++) p_gauss_weights[j] = M_PI / double(nn+1) ; 
    }
    
    return p_gauss_weights[i] ;

}


double Chebyshev_poly::gauss_gamma(int i) const {

    if (p_gauss_gamma == 0x0) {  // the weights must be computed

        p_gauss_gamma = new double[nn+1] ;
        
        p_gauss_gamma[0] = M_PI ; 
        
        for (int j=1; j<=nn; j++) p_gauss_gamma[j] = 0.5 * M_PI  ; 
    }
    
    return p_gauss_gamma[i] ;

}


const Grid& Chebyshev_poly::gauss_lobatto_nodes() const {
    
    if (p_gauss_lobatto_nodes == 0x0) {  // the nodes must initialized

        p_gauss_lobatto_nodes = new Grid_Chebyshev_GL(nn+1) ; 
    }
    
    return *p_gauss_lobatto_nodes ;
}


double Chebyshev_poly::gauss_lobatto_weight(int i) const {

    if (p_gauss_lobatto_weights == 0x0) {  // the weights must be computed

        p_gauss_lobatto_weights = new double[nn+1] ;
        
        p_gauss_lobatto_weights[0] = M_PI / double(2*nn) ;  

        for (int j=1; j<nn; j++) 
            p_gauss_lobatto_weights[j] = M_PI / double(nn) ; 

        p_gauss_lobatto_weights[nn] = M_PI / double(2*nn) ;  
    }
    
    return p_gauss_lobatto_weights[i] ;

}


double Chebyshev_poly::gauss_lobatto_gamma(int i) const {

    if (p_gauss_lobatto_gamma == 0x0) {  // the weights must be computed

        p_gauss_lobatto_gamma = new double[nn+1] ;
        
        p_gauss_lobatto_gamma[0] = M_PI ; 
        
        for (int j=1; j<nn; j++) p_gauss_lobatto_gamma[j] = 0.5 * M_PI  ; 

        p_gauss_lobatto_gamma[nn] = M_PI ; 
        
    }
    
    return p_gauss_lobatto_gamma[i] ;

}


//---------------------------------------------------//
//  Gauss-Lobatto interpolant polynomial via a FFT   //
//---------------------------------------------------//

void Chebyshev_poly::coef_interpolant_GL_FFT(double (*f)(double), 
                                         double* cf) const {

    // Values of the function at the Gauss-Lobatto nodes
    
    double* ff = new double[nn+1] ; 
    
    for (int i=0; i<=nn; i++) ff[i] = f( gauss_lobatto_nodes()(i) ) ; 
    
    // Chebyshev transform via a FFT
    
    coef_cheb_fft(nn+1, ff, cf) ; 
        
    delete [] ff ; 

}


//---------------------------------------------//
//  Coefficient of the orthogonal projection   //
//---------------------------------------------//

void Chebyshev_poly::coef_projection(double (*f)(double), double* cf) const {

    // The computation is an approximate one:
    // it returns the coefficients of an interpolating polynomial with$
    // a large number of Gauss-Lobatto nodes
    
    int n_large = 128 ; 
    
    if (nn > n_large) {
        cerr << "Chebyshev_poly::coef_projection : nn > n_large !"
            << endl << "  nn = " << nn << "  n_large = " << n_large << endl ; 
        abort() ; 
    }
    
    Chebyshev_poly cheb_large(n_large) ;
    
    double* cf_large = new double[n_large + 1] ;
    
    cheb_large.coef_interpolant_GL(f, cf_large) ; 
    
    for (int i=0; i<=nn; i++) cf[i] = cf_large[i] ;  
    
    delete [] cf_large ; 
    
}












