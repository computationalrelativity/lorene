/*
 *  Methods of class Bin_ns_bh to set analytical value to omega
 *
 *    (see file bin_ns_bh.h for documentation).
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

char bin_ns_bh_omegaana_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2005/11/30 11:09:06  p_grandclement
 * Changes for the Bin_ns_bh project
 *
 * Revision 1.1  2004/06/09 06:27:40  k_taniguchi
 * First revision.
 *
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "bin_ns_bh.h"
#include "unites.h"

void Bin_ns_bh::analytical_omega() {

    // NS-BH binary systems should be relativistic
    // -------------------------------------------
    if ( !star.is_relativistic() ) {
        abort() ;
    }

    using namespace Unites ;

    double rr = separation() ;
    double mtot = star.mass_g() + hole.masse_adm_seul() / ggrav ;

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

    set_omega (sqrt( omega2 )) ;
    
    // The derived quantities are obsolete:
    del_deriv() ;

}
