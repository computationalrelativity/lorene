/*
 *  Method for vector Poisson equation inverting eqs. for V^r and eta as a block.
 *
 *    (see file vector.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005  Jerome Novak
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

char vector_poisson_block_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2005/02/14 13:01:50  j_novak
 * p_eta and p_mu are members of the class Vector. Most of associated functions
 * have been moved from the class Vector_divfree to the class Vector.
 *
 *
 * $Header$
 *
 */

// C headers
#include <assert.h>
#include <stdlib.h>
#include <math.h>

// Lorene headers
#include "tensor.h"
#include "diff.h"
#include "proto.h"

void Vector::poisson_block(double, Vector& ) const {

    const Map_af* mpaff = dynamic_cast<const Map_af*>(mp) ;
#ifndef NDEBUG 
    for (int i=0; i<3; i++)
	assert(cmp[i]->check_dzpuis(4)) ;
    // All this has a meaning only for spherical components:
    const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
    assert(bvs != 0x0) ; 
    //## ... and affine mapping, for the moment!
    assert( mpaff != 0x0) ;
#endif

    cout << *mpaff << endl ;
    cout << "Not ready yet ..." << endl ;
    abort() ;

  return ;

}
