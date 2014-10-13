/*
 *  Methods of class Bin_bhns_extr to set analytical value to omega
 *
 *    (see file bin_bhns_extr.h for documentation).
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

char bin_bhns_extr_omegaana_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:52:42  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/11/30 20:46:36  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include <cmath>

// Lorene headers
#include "bin_bhns_extr.h"
#include "unites.h"

namespace Lorene {
void Bin_bhns_extr::analytical_omega() {

  using namespace Unites ;

    // BH-NS binary systems should be relativistic
    // -------------------------------------------
    if ( !star.is_relativistic() ) {

        cout << "BH-NS binary systems should be relativistic !!!" << endl ;
        abort() ;
    }

    double rr = separ ;
    double mtot = mass_bh ; // Approximates the extreme mass ratio

    // Compaction factor
    double compact = ggrav * mtot / rr ;

    double omega2 ;

    if ( star.is_irrotational() ) {

        // Irrotational case
        // -----------------

        omega2 = ggrav * mtot / pow(rr, 3.)
	  * (1. - 2.75 * compact + 8.625 * compact*compact ) ;

    }
    else {
        // Corotating case
        // ---------------

        // a0/R
        double a0sr = star.ray_eq() / rr ;

	// Rescaled moment of inertia 5 I / (2 M a0^2)
	double ired = double(5)/double(3) * ( 1. - double(6)/M_PI/M_PI ) ;
	omega2 = ggrav * mtot / pow(rr, 3.)
	  * (1. - compact * ( 2.75 + 2.*a0sr*a0sr * ired
			      - 0.48*pow(a0sr, 4) * ired*ired )
	     + compact*compact * ( 8.625 + 2.75*a0sr*a0sr * ired
				   + 2.*pow(a0sr, 4) * ired*ired ) ) ;

    }

    omega = sqrt( omega2 ) ;

    // The derived quantities are obsolete:
    del_deriv() ;

}
}
