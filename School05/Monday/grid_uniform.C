/*
 * Definition of class Grid_uniform
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

#include "grid.h"

//----------------------//
// Standard constructor //
//----------------------//

Grid_uniform::Grid_uniform(int nb_nodes) : Grid(nb_nodes) {
    
    for (int i=0; i<=nn; i++) xx[i] = -1. + 2 * double(i) / double(nn) ; 

}


//--------------------//
//  Copy constructor  //
//--------------------//

Grid_uniform::Grid_uniform(const Grid_uniform& gi) : Grid(gi) {}


//--------------//
//  Destructor  //
//--------------//

Grid_uniform::~Grid_uniform() {}


//--------------//
//  Assignment  //
//--------------//

void Grid_uniform::operator=(const Grid_uniform& gi) {

    Grid::operator=(gi) ;     
}



//-----------//
//  Display  //
//-----------//

ostream& operator<<(ostream& ost, const Grid_uniform& xn) {

    ost << "Uniform grid : " ;
    int nn = xn.n() ; 
    ost << "number of nodes : " << nn + 1 << endl ; 
    for (int i=0; i<=nn; i++) {
        ost << "   x(" << i << ") = " << xn(i) << endl ; 
    }
    
    return ost ; 
}

