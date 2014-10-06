/*
 * Definition of class Legendre_poly
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
 * Revision 1.3  2014/10/06 15:09:48  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2005/11/14 14:12:10  e_gourgoulhon
 * Added include <assert.h>
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

#include <cstdlib>
#include <cmath>
#include <cassert>

#include "ortho_poly.h"
#include "grid.h"

void legendre_nodes_weight_GL(int n, double* coloc, double* weight, 
                              double prec = 1.e-12, int itemax = 100) ; 
     

//----------------------//
// Standard constructor //
//----------------------//

Legendre_poly::Legendre_poly(int ni) : Ortho_poly(ni) {}

//--------------------//
//  Copy constructor  //
//--------------------//

Legendre_poly::Legendre_poly(const Legendre_poly& pp) : Ortho_poly(pp) {}


//--------------//
//  Destructor  //
//--------------//

Legendre_poly::~Legendre_poly() {}


//--------------//
//  Assignment  //
//--------------//

void Legendre_poly::operator=(const Legendre_poly& ) {

    cerr << "Legendre_poly::operator= not implemented !" << endl ; 
    abort() ; 
    
}

//-----------//
//  Display  //
//-----------//

ostream& operator<<(ostream& ost, const Legendre_poly& pp) {

    ost << "Basis of Legendre polynomials up to degree " << pp.n() << endl ; 
    
    return ost ; 
}

//-------------------//
//  Weight function  //
//-------------------//

double Legendre_poly::weight(double ) const {
    
    return 1. ; 
    
} 


//-------------------------------------------//
//  Evaluation of the Legendre polynomials   //
//-------------------------------------------//

double Legendre_poly::operator()(int i, double x) const {

    assert( i >= 0 ) ; 
    assert( i <= nn ) ; 
    
    if (i==0) return 1. ; 

    if (i==1) return x ; 
    
    double tjm2 = 1. ;  // value at step j - 2
    double tjm1 = x  ;  // value at step j - 1
    
    for (int j=2; j<=i; j++) {
	    double tj = ( (2*j-1)*x* tjm1 - (j-1)*tjm2 ) / double(j) ;
        tjm2 = tjm1 ; 
        tjm1 = tj ;       	    
    }
    
    return tjm1 ;     
    
}

//------------------------------------//
//  Computation of nodes and weights  //
//------------------------------------//


const Grid& Legendre_poly::gauss_nodes() const {
    
    if (p_gauss_nodes == 0x0) {  // the nodes must initialized
        
        p_gauss_nodes = new Grid_Legendre_Gauss(nn+1) ;
    }
    
    return *p_gauss_nodes ;
}


double Legendre_poly::gauss_weight(int i) const {

    if (p_gauss_weights == 0x0) {  // the weights must be computed

        cerr << "Legendre_poly::gauss_weight : not implemented yet !"
            << endl ; 
        abort() ; 
    }
    
    return p_gauss_weights[i] ;

}


double Legendre_poly::gauss_gamma(int i) const {

    if (p_gauss_gamma == 0x0) {  // the weights must be computed

        p_gauss_gamma = new double[nn+1] ;
        
        for (int j=0; j<=nn; j++) 
            p_gauss_gamma[j] = 1. / (double(j) + 0.5) ; 
    }
    
    return p_gauss_gamma[i] ;

}


const Grid& Legendre_poly::gauss_lobatto_nodes() const {
    
    if (p_gauss_lobatto_nodes == 0x0) {  // the nodes must initialized

        p_gauss_lobatto_nodes = new Grid_Legendre_GL(nn+1) ; 
    }
    
    return *p_gauss_lobatto_nodes ;
}


double Legendre_poly::gauss_lobatto_weight(int i) const {

    if (p_gauss_lobatto_weights == 0x0) {  // the weights must be computed

        p_gauss_lobatto_weights = new double[nn+1] ;
        
        double* xx_work = new double[nn+1] ;
        
        double precis = 1.e-12 ;    // required precision
        int nitermax = 200 ;     // Maximum number of iterations
    
        legendre_nodes_weight_GL(nn+1, xx_work, p_gauss_lobatto_weights, 
                                 precis, nitermax) ; 
        
        delete [] xx_work ; 
    }
    
    return p_gauss_lobatto_weights[i] ;

}


double Legendre_poly::gauss_lobatto_gamma(int i) const {

    if (p_gauss_lobatto_gamma == 0x0) {  // the weights must be computed

        p_gauss_lobatto_gamma = new double[nn+1] ;
        
        for (int j=0; j<nn; j++) 
            p_gauss_lobatto_gamma[j] = 1. / (double(j) + 0.5) ; 

        p_gauss_lobatto_gamma[nn] = 2. / double(nn) ; 
        
    }
    
    return p_gauss_lobatto_gamma[i] ;

}


//---------------------------------------------//
//  Coefficient of the orthogonal projection   //
//---------------------------------------------//

void Legendre_poly::coef_projection(double (*f)(double), double* cf) const {

    // The computation is an approximate one:
    // it returns the coefficients of an interpolating polynomial with$
    // a large number of Gauss-Lobatto nodes
    
    int n_large = 128 ; 
    
    if (nn > n_large) {
        cerr << "Legendre_poly::coef_projection : nn > n_large !"
            << endl << "  nn = " << nn << "  n_large = " << n_large << endl ; 
        abort() ; 
    }
    
    Legendre_poly leg_large(n_large) ;
    
    double* cf_large = new double[n_large + 1] ;
    
    leg_large.coef_interpolant_GL(f, cf_large) ; 
    
    for (int i=0; i<=nn; i++) cf[i] = cf_large[i] ;  
    
    delete [] cf_large ; 
    
}












