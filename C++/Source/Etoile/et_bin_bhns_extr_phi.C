/*
 *  Method of class Et_bin_bhns_extr to search the position of the longest
 *   radius from the position of the maximum enthalpy
 *  The code returns the position of "phi" only because "xi=1" and
 *   "theta=pi/2".
 *
 *    (see file et_bin_bhns_extr.h for documentation).
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

char et_bin_bhns_extr_phi_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/11/30 20:50:48  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "et_bin_bhns_extr.h"
#include "utilitaires.h"

double Et_bin_bhns_extr::phi_longest_rad(double x_max, double y_max) const {

    //----------------------------------------------------------------//
    //          Construct the surface function S(theta, phi)          //
    //----------------------------------------------------------------//

    // The following is required to access functions of Map_et
    Map_et& mp_et = dynamic_cast<Map_et&>(mp) ;

    const Valeur& ff0 = mp_et.get_ff() ;
    const Valeur& gg0 = mp_et.get_gg() ;

    Valeur fff = ff0 ;
    Valeur ggg = gg0 ;
    Valeur dff = fff.dsdp() ;
    Valeur dgg = ggg.dsdp() ;

    double ppp = M_PI/2. ; // Initial position of the phi-coordinate
    double ptmp ;
    int mm ;          // Number of steps to the phi-direction
    double dp = 1. ;  // Step interval to the phi-direction, initialized to 1
    double diff ;
    double diff_prev ;
    double ss ;

    while ( dp > 1.e-15 ) {

        diff = 1. ;
	mm = 0 ;
	dp = 0.1 * dp ;

	diff_prev = ( dff.val_point(0,1.,M_PI/2.,ppp)
		      + dgg.val_point(0,1.,M_PI/2.,ppp) )
	  * ( 1. + ff0.val_point(0,1.,M_PI/2.,ppp)
	      + gg0.val_point(0,1.,M_PI/2.,ppp)
	      - x_max * cos(ppp) - y_max * sin(ppp) )
	  - ( 1. + ff0.val_point(0,1.,M_PI/2.,ppp)
	      + gg0.val_point(0,1.,M_PI/2.,ppp) )
	  * ( - x_max * sin(ppp) + y_max * cos(ppp) ) ;

	if ( diff_prev > 0. ) {
	    ss = 1. ;
	}
	else {
	    ss = -1. ;
	}

	while ( diff > 1.e-15 ) {

	    mm++ ;
	    ptmp = ppp + mm * dp ;

	    diff = ss * ( ( dff.val_point(0,1.,M_PI/2.,ptmp)
			    + dgg.val_point(0,1.,M_PI/2.,ptmp) )
			  * ( 1. + ff0.val_point(0,1.,M_PI/2.,ptmp)
			      + gg0.val_point(0,1.,M_PI/2.,ptmp)
			      - x_max * cos(ptmp) - y_max * sin(ptmp) )
			  - ( 1. + ff0.val_point(0,1.,M_PI/2.,ptmp)
			      + gg0.val_point(0,1.,M_PI/2.,ptmp) )
			  * ( - x_max * sin(ptmp) + y_max * cos(ptmp) ) ) ;

	}
	ppp += ss * (mm - 1) * dp ;

    }

    return ppp ;

}
