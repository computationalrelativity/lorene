/*
 *  Methods for the Diff class.
 *
 *    (see file diff.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Jerome Novak
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

char diff_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2005/01/10 16:34:52  j_novak
 * New class for 1D mono-domain differential operators.
 *
 *
 * $Header$
 *
 */

// C headers
#include <assert.h>

// Lorene headers
#include "diff.h"


Diff::Diff(int base_r, int nr) : base(base_r >> TRA_R), npoints(nr) {

    assert (base < MAX_BASE) ;
    assert (npoints < max_points) ;

}

Diff::Diff(const Diff& diff_in) : base(diff_in.base), 
				  npoints(diff_in.npoints) {
    assert (base < MAX_BASE) ;
    assert (npoints < max_points) ;

}    

Diff::~Diff() {}

void Diff::operator=(const Diff& diff_in) {

    base = diff_in.base ;
    npoints = diff_in.npoints ;
    assert (base < MAX_BASE) ;
    assert (npoints < max_points) ;

}

ostream& operator<<(ostream& ost, const Diff& ope) {

    ost << "Differential operator : " ;
    
    ope >> ost ;

    ost << "Radial base: " ;

    switch (ope.base) {

	case R_CHEB >> TRA_R :
	    ost << "Chebyshev polynomials (R_CHEB)"  ;
	    break ;

	case R_CHEBP >> TRA_R :
	    ost << "Even Chebyshev polynomials (R_CHEBP)" ;
	    break ;

	case R_CHEBI >> TRA_R :
	    ost << "Odd Chebyshev polynomials (R_CHEBI)"  ;
	    break ;

	case R_CHEBU >> TRA_R :
	    ost << "Chebyshev polynomials / compactified domain (R_CHEBU)" ;
	    break ;

	default:
	    ost << "unknown!" << endl ;
    }

    ost << " with " << ope.npoints << " coefficients." << endl ;
    ost << endl ;

    return ost ;
}
