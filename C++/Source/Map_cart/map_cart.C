/*
 *  Methods of class Map_cart
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

char map_cart_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:42  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1  2002/03/15 13:16:23  n_chamel
 * Mapping between grid coordinates and physical coordinates
 *
 *
 *
 *
 * $Header$
 *
 */

// C headers

// Lorene headers
#include "headcpp.h"
#include "map_cart.h"

                                                                //-------------------------//
                                                                //              Constructors                      //
                                                               //-------------------------//


// Standard constructor
// ---------------

Map_cart:: Map_cart(const Grille_cart& grid_i) : grid(&grid_i) {}

// Copy constructor
// ------------

Map_cart:: Map_cart(const Map_cart& map_i) : grid(map_i.grid) {}

                                                               //-------------------------//
                                                               //              Destructors                        //
                                                               //-------------------------//

Map_cart::~Map_cart(){}


                                                       //--------------------//
                                                        //              Outputs                  //
                                                       //--------------------//



  // Operateurs <<
ostream& operator<<(ostream& ost, const Map_cart & map)  {
    map >> ost ;
    return ost ;
}

