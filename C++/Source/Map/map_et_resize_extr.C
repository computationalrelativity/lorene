/*
 *  Method of the class Map_et to compute the rescale of the outermost domain
 *  in the case of non-compactified external domain.
 *
 *    (see file map.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Keisuke Taniguchi
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

char map_et_resize_extr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/11/30 20:54:24  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include <assert.h>

// Lorene headers
#include "map.h"

void Map_et::resize_extr(double lambda) {

    // Protections
    // -----------
    int l = mg->get_nzone() - 1 ;

    if (mg->get_type_r(l) != FIN) {
        cout << "Map_et::resize_extr can be applied only to a shell !"
	     << endl ;
	abort() ;
    }

    // Assertion
    // ---------
    assert(mg->get_nzone() >= 3) ;  // The outermost domain should be
                                    //  a spherical shell in this method.

    // New values of alpha and beta in the outermost domain :
    // ----------------------------------------------------
    double n_alpha = 0.5 * ( (lambda + 1.) * alpha[l]
			     + (lambda - 1.) * beta[l] ) ;

    double n_beta = 0.5 * ( (lambda - 1.) * alpha[l]
			    + (lambda + 1.) * beta[l] ) ;

    alpha[l] = n_alpha ;
    beta[l] = n_beta ;

    // The coords are no longer up to date :
    reset_coord() ;

}
