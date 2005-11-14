/*
 * Definition of class Grid_Legendre_GL
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

void legendre_nodes_weight_GL(int n, double* coloc, double* weight, 
                              double prec = 1.e-12, int itemax = 100) ; 
     

//----------------------//
// Standard constructor //
//----------------------//

Grid_Legendre_GL::Grid_Legendre_GL(int nb_nodes) : Grid(nb_nodes) {
    
    double* weights = new double[nn+1] ;
    
    double precis = 1.e-12 ;    // required precision
    int nitermax = 200 ;     // Maximum number of iterations
    
    legendre_nodes_weight_GL(nn+1, xx, weights, precis, nitermax) ; 
    
    delete [] weights ; 
}


//--------------------//
//  Copy constructor  //
//--------------------//

Grid_Legendre_GL::Grid_Legendre_GL(const Grid_Legendre_GL& gi) : Grid(gi) {}


//--------------//
//  Destructor  //
//--------------//

Grid_Legendre_GL::~Grid_Legendre_GL() {}


//--------------//
//  Assignment  //
//--------------//

void Grid_Legendre_GL::operator=(const Grid_Legendre_GL& gi) {

    Grid::operator=(gi) ;     
}



//-----------//
//  Display  //
//-----------//

ostream& operator<<(ostream& ost, const Grid_Legendre_GL& xn) {

    ost << "Legendre Gauss-Lobatto nodes : " ;
    int nn = xn.n() ; 
    ost << "number of nodes : " << nn + 1 << endl ; 
    for (int i=0; i<=nn; i++) {
        ost << "   x(" << i << ") = " << xn(i) << endl ; 
    }
    
    return ost ; 
}

