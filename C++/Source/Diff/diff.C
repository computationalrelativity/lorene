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

 

/*
 * $Id$
 * $Log$
 * Revision 1.6  2016/12/05 16:17:50  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:52:50  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:04  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2007/12/11 15:28:11  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.2  2005/02/09 09:53:24  j_novak
 * Removed irrelevant asserts on number of points.
 *
 * Revision 1.1  2005/01/10 16:34:52  j_novak
 * New class for 1D mono-domain differential operators.
 *
 *
 * $Header$
 *
 */

// C headers
#include <cassert>

// Lorene headers
#include "diff.h"


namespace Lorene {
Diff::Diff(int base_r, int nr) : base(base_r >> TRA_R), npoints(nr) {

    assert (base < MAX_BASE) ;

}

Diff::Diff(const Diff& diff_in) : base(diff_in.base), 
				  npoints(diff_in.npoints) {
    assert (base < MAX_BASE) ;

}    

Diff::~Diff() {}

void Diff::operator=(const Diff& diff_in) {

    base = diff_in.base ;
    npoints = diff_in.npoints ;
    assert (base < MAX_BASE) ;

}

ostream& operator<<(ostream& ost, const Diff& ope) {

    ost << "Differential operator : " ;
    
    ope >> ost ;

    ost << "Radial base: " ;

    switch (ope.base) {

	case R_CHEB >> TRA_R :
	    ost << "Chebyshev polynomials (R_CHEB)"  ;
	    break ;

	case R_JACO02 >> TRA_R :
	    ost << "Jacobi(0,2) polynomials (R_JACO02)" ;
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
}
