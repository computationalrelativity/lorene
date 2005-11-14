/*
 * Definition of class Ortho_poly
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
#include <math.h>

#include "ortho_poly.h"
#include "grid.h"


//--------------------//
//  Copy constructor  //
//--------------------//

Ortho_poly::Ortho_poly(const Ortho_poly& pp) : nn(pp.nn) {

    // Initialization of the pointers of the derived quantities to 0x0
    p_gauss_nodes = 0x0 ;             
    p_gauss_weights = 0x0 ;            
    p_gauss_gamma = 0x0 ;            
    p_gauss_lobatto_nodes = 0x0 ;    
    p_gauss_lobatto_weights = 0x0 ;    
    p_gauss_lobatto_gamma = 0x0 ;    
    
}

//---------------------------------//
// Constructor for derived classes //
//---------------------------------//

Ortho_poly::Ortho_poly(int ni) : nn(ni) {
    
    assert(ni >= 0) ; 
    
    // Initialization of the pointers of the derived quantities to 0x0
    p_gauss_nodes = 0x0 ;             
    p_gauss_weights = 0x0 ;            
    p_gauss_gamma = 0x0 ;            
    p_gauss_lobatto_nodes = 0x0 ;    
    p_gauss_lobatto_weights = 0x0 ;    
    p_gauss_lobatto_gamma = 0x0 ;    
    
}


//--------------//
//  Destructor  //
//--------------//

Ortho_poly::~Ortho_poly() {

    if (p_gauss_nodes != 0x0) delete p_gauss_nodes ; 
    if (p_gauss_weights != 0x0) delete [] p_gauss_weights ; 
    if (p_gauss_gamma != 0x0) delete [] p_gauss_gamma ; 
    if (p_gauss_lobatto_nodes != 0x0) delete p_gauss_lobatto_nodes ; 
    if (p_gauss_lobatto_weights != 0x0) delete [] p_gauss_lobatto_weights ; 
    if (p_gauss_lobatto_gamma != 0x0) delete [] p_gauss_lobatto_gamma ; 

}

//--------------//
//  Assignment  //
//--------------//

void Ortho_poly::operator=(const Ortho_poly& ) {

    cerr << "Ortho_poly::operator= not implemented !" << endl ; 
    abort() ; 
    
}

//-------------//
//  Accessors  //
//-------------//

int Ortho_poly::n() const {

    return nn ; 
    
}


//-------------------------------------------------//
//  Coefficients of Gauss interpolant polynomial   //
//-------------------------------------------------//

void Ortho_poly::coef_interpolant_Gauss(double (*f)(double), 
                                         double* cf) const {

    // Values of the function at the Gauss nodes
    
    double* ff = new double[nn+1] ; 
    
    for (int i=0; i<=nn; i++) ff[i] = f( gauss_nodes()(i) ) ; 
    
    // Computation of the coefficients
    
    for (int i=0; i<=nn; i++) {
        double sum = 0 ; 
        for (int j=0; j<=nn; j++) {
            sum += ff[j] * operator()(i, gauss_nodes()(j) )
                    * gauss_weight(j) ;
        }
        cf[i] = sum / gauss_gamma(i) ; 
    }
        
    delete [] ff ; 

}

//--------------------------------------------------------//
//  Coefficients of Gauss-Lobatto interpolant polynomial  //
//--------------------------------------------------------//

void Ortho_poly::coef_interpolant_GL(double (*f)(double), 
                                         double* cf) const {

    // Values of the function at the Gauss-Lobatto nodes
    
    double* ff = new double[nn+1] ; 
    
    for (int i=0; i<=nn; i++) ff[i] = f( gauss_lobatto_nodes()(i) ) ; 
    
    // Computation of the coefficients
    
    for (int i=0; i<=nn; i++) {
        double sum = 0 ; 
        for (int j=0; j<=nn; j++) {
            sum += ff[j] * operator()(i, gauss_lobatto_nodes()(j) )
                    * gauss_lobatto_weight(j) ;
        }
        cf[i] = sum / gauss_lobatto_gamma(i) ; 
    }
        
    delete [] ff ; 

}

//----------------------------------------//
//  Values of the interpolant polynomial  //
//----------------------------------------//

double Ortho_poly::series(const double* a, double x) const {

    double sum = 0 ;     
    for (int i=0; i<=nn; i++) {
        sum += a[i] * operator()(i, x) ;
    }
        
    return sum ; 

}


