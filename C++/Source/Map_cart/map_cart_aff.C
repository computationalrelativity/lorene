/*
 *  Methods of class Map_cart_aff
 *
 */

/*
 *   Copyright (c) 2002 Nicolas Chamel
 *   Copyright (c) 2002 Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char map_cart_aff_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/03/15 13:16:23  n_chamel
 * Mapping between grid coordinates and physical coordinates
 *
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include <iostream.h>

// C headers

// Lorene headers
#include "map_cart.h"
#include "utilitaires.h"

                                                                //-------------------------//
                                                                //              Constructors                      //
                                                               //-------------------------//


// Standard constructor
// ---------------

Map_cart_aff::Map_cart_aff(const Grille_cart&  grid_i,   double alpha1_i, double beta1_i,
                                                                                double alpha2_i, double beta2_i,
                                                                                double alpha3_i, double beta3_i) :
        Map_cart(grid_i),
        alpha1(alpha1_i),
        beta1(beta1_i),
        alpha2(alpha2_i),
        beta2(beta2_i),
        alpha3(alpha3_i),
        beta3(beta3_i)
        {}



// Copy constructor
// ------------

Map_cart_aff:: Map_cart_aff(const Map_cart_aff& map_i) :
        Map_cart(map_i),
        alpha1(map_i.alpha1),
        beta1(map_i.beta1),
        alpha2(map_i.alpha2),
        beta2(map_i.beta2),
        alpha3(map_i.alpha3),
        beta3(map_i.beta3)
        {}


// Constructor  from file
// ---------------

Map_cart_aff:: Map_cart_aff(const Grille_cart&  grid_i,   FILE* fich)  :  Map_cart(grid_i)
{

                fread_be(&alpha1, sizeof(double), 1, fich) ;
                fread_be(&beta1, sizeof(double), 1, fich) ;
                fread_be(&alpha2, sizeof(double), 1, fich) ;
                fread_be(&beta2, sizeof(double), 1, fich) ;
                fread_be(&alpha3, sizeof(double), 1, fich) ;
                fread_be(&beta3, sizeof(double), 1, fich) ;

}





                                                              //-------------------------//
                                                               //              Destructors                        //
                                                               //-------------------------//

Map_cart_aff::~Map_cart_aff(){}


                                                       //--------------------//
                                                        //              Outputs                  //
                                                       //--------------------//



// Operateurs >>
ostream & Map_cart_aff::operator>>(ostream & ost) const {

    ost << "Affine Cartesian mapping (class Map_cart_aff)" << endl ;

    ost << "  alpha1 = " << alpha1 << "  beta1 = " << beta1 << endl ;
    ost << "  alpha2 = " << alpha2 << "  beta2 = " << beta2 << endl ;
    ost << "  alpha3 = " << alpha3 << "  beta3 = " << beta3 << endl ;

    return ost ;
}

// Save in a file
void Map_cart_aff::sauve(FILE* fich) const {

                fwrite_be(&alpha1, sizeof(double), 1, fich) ;
                fwrite_be(&beta1, sizeof(double), 1, fich) ;
                fwrite_be(&alpha2, sizeof(double), 1, fich) ;
                fwrite_be(&beta2, sizeof(double), 1, fich) ;
                fwrite_be(&alpha3, sizeof(double), 1, fich) ;
                fwrite_be(&beta3, sizeof(double), 1, fich) ;

}
