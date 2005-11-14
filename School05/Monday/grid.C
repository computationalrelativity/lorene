/*
 * Definition of class Grid
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

//----------------------//
// Standard constructor //
//----------------------//

Grid::Grid(int nb_nodes, double* xi) : nn(nb_nodes-1) {
    
    assert(nb_nodes > 1) ; 

    xx = new double[nb_nodes] ; 
    
    for (int i=0; i<=nn; i++) {
        double xi1 = xi[i] ; 
        if ((xi1 < -1) || (xi1 > 1)) {
            cerr << "Grid::Grid: node no. " << i << " is not in [-1,1] !" << endl ;
            abort() ;  
        }
        for (int j=0; j<i; j++) {
            if (xi1 == xx[j]) {
                cerr << "Grid::Grid: the nodes must be different !" << endl ; 
                abort() ; 
            }
        }
        
        xx[i] = xi1 ; 
    }
}


//--------------------//
//  Copy constructor  //
//--------------------//

Grid::Grid(const Grid& nod) : nn(nod.nn) {

    xx = new double[nn+1] ; 
    for (int i=0; i<=nn; i++) {
        xx[i] = nod.xx[i] ; 
    }
}

//---------------------------------//
// Constructor for derived classes //
//---------------------------------//

Grid::Grid(int nb_nodes) : nn(nb_nodes-1) {
    
    assert(nb_nodes > 1) ; 
    
    xx = new double[nb_nodes] ; 
    
    // Only the memory allocation for the array xx is performed here,
    // setting the values in the array is left to the derived classes
}


//--------------//
//  Destructor  //
//--------------//

Grid::~Grid() {

    delete [] xx ; 
    
}


//--------------//
//  Assignment  //
//--------------//

void Grid::operator=(const Grid& nod) {

    assert( nod.nn == nn ) ;
    for (int i=0; i<=nn; i++) {
        xx[i] = nod.xx[i] ; 
    }
    
}


//-------------//
//  Accessors  //
//-------------//

int Grid::n() const {

    return nn ; 
    
}

double Grid::operator()(int i) const {

    assert( (i >= 0) && (i<=nn) ) ; 

    return xx[i] ; 
}


//-----------//
//  Display  //
//-----------//

ostream& operator<<(ostream& ost, const Grid& xn) {

    int nn = xn.n() ; 
    ost << "Number of nodes : " << nn + 1 << endl ; 
    for (int i=0; i<=nn; i++) {
        ost << "   x(" << i << ") = " << xn(i) << endl ; 
    }
    
    return ost ; 
}


void Grid::plot(int color, int nfig, double ymin, double ymax, 
            const char* title, const char* label_y, 
            const char* device) const {
    
    for (int i=0; i<=nn; i++) {

        plot_point(xx[i], 0., color, nfig, ymin, ymax, title, label_y, device) ;
        
    }    

}


//-----------------------//
// Lagrange polynomials  //
//-----------------------//

double Grid::lagrange(int i, double x) const {

    double y = 1 ; 
    
    for (int j = 0; j<=nn ; j++) {
        if (j==i) continue ; 
        y *= (x - xx[j]) / (xx[i] - xx[j]) ;         
    }
    
    return y ; 
}


//-------------------//
// Nodal polynomial  //
//-------------------//

double Grid::nodal_polynomial(double x) const {
    
    double y = 1 ; 
    
    for (int i = 0; i<=nn ; i++) {
        y *= x - xx[i] ;         
    }
    
    return y ; 
    
}

//---------------------------//
// Interpolating polynomial  //
//---------------------------//
        
double Grid::interpole(double (*f)(double), double x) const {

    double y = 0 ; 
    
    for (int i = 0; i<=nn ; i++) {
        y += f(xx[i]) * lagrange(i, x) ;         
    }
    
    return y ; 
        
} 

//--------------------//
// Lebesgue constant  //
//--------------------//
        
double Grid::lebesgue_constant() const {

    int np = 10000 ; 
    double lbc = 0 ; 
    
    for (int j=0; j<np; j++) {  // scan of [-1,1]
        double x = -1. + 2. * double(j) / double(np-1) ; 
        double som = fabs( lagrange(0, x) ) ; 
        for (int i=1; i<=nn; i++) {
            som += fabs( lagrange(i, x) ) ; 
        }
        if (som > lbc) lbc = som ; 
    }
    
    return lbc ; 
    
} 
