/*
 *  Method of class Champ_cart   to compute the spectral coefficients
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

char champ_cart_coef_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/03/28 10:47:29  n_chamel
 * New class Champ_cart for fields on Cartesian grids.
 *
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include <iostream.h>

// Lorene headers
#include "champ_cart.h"
#include "type_parite.h"
#include "utilitaires.h"

void Champ_cart::coef(int base_x, int base_y, int base_z) const {

        base[0] = base_x ;
        base[1] = base_y ;
        base[2] = base_z ;

}
