/*
 *  Method of class Et_bin_bhns_extr to search the position of the maximum
 *   enthalpy
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

char et_bin_bhns_extr_max_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/11/30 20:50:24  k_taniguchi
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
//#include "utilitaires.h"

void Et_bin_bhns_extr::ent_max_search(double& xx, double& yy) const {

    //------------------------------------------------------------------//
    //          Calculate the derivative of the enthalpy field          //
    //------------------------------------------------------------------//

    const Tenseur& dent = ent.gradient() ;

    double xxp = 0. ;  // Position of the x-coordinate, initialized to zero
    double yyp = 0. ;  // Position of the y-coordinate, initialized to zero
    double xtmp, ytmp ;
    int mm, nn ;       // Number of steps to the x and y directions
    double rr = 0. ;   // r coordinate, initialized to zero
    double pp = M_PI/2. ;    // phi coordinate, initialized to M_PI/2.
    double dval_x ;    // Direction of dent(0) (1 or -1)
    double dval_y ;    // Direction of dent(1) (1 or -1)
    double ss ;

    while ( fabs(dent(0).val_point(rr, M_PI/2., pp)) > 1.e-15 ||
	    fabs(dent(1).val_point(rr, M_PI/2., pp)) > 5.e-15) {

        double dx = 1. ; // Step interval to the x-direction, initialized to 1
	double dy = 1. ; // Step interval to the y-direction, initialized to 1
	double diff_dent_x ;
	double diff_dent_y ;

	while ( dy > 1.e-15 ) {

	    diff_dent_y = 1. ;
	    nn = 0 ;
	    dy = 0.1 * dy ;

	    rr = sqrt( xxp*xxp + yyp*yyp ) ;

	    if ( xxp == 0. ) {
	        pp = M_PI/2. ;  // There is a possibility of (pp = 1.5*M_PI)
	    }
	    else {
	        pp = acos( xxp / rr ) ;
	    }

	    dval_y = dent(1).val_point(rr, M_PI/2., pp) ;

	    if ( dval_y > 0. ) {
	        ss = 1. ;
	    }
	    else {
	        ss = -1. ;
	    }

	    while( diff_dent_y > 1.e-15 ) {

	        nn++ ;
		ytmp = yyp + ss * nn * dy ;

		rr = sqrt( xxp*xxp + ytmp*ytmp ) ;

		if ( xxp == 0. ) {
		    if ( ss > 0. ) {
		        pp = M_PI/2. ;
		    }
		    else {
		        pp = 1.5*M_PI ;
		    }
		}
		else {
		    pp = acos( xxp / rr ) ;
		}

		diff_dent_y = ss * dent(1).val_point(rr, M_PI/2., pp) ;

	    }
	    yyp += ss * (nn - 1) * dy ;

	}

	while ( dx > 1.e-15 ) {

	    diff_dent_x = 1. ;
	    mm = 0 ;
	    dx = 0.1 * dx ;

	    rr = sqrt( xxp*xxp + yyp*yyp ) ;

	    if ( xxp == 0. ) {
	        pp = M_PI/2. ;  // There is a possibility of (pp = 1.5*M_PI)
	    }
	    else {
	        pp = acos( xxp / rr ) ;
	    }

	    dval_x = dent(0).val_point(rr, M_PI/2., pp) ;

	    if ( dval_x > 0. ) {
	        ss = 1. ;
	    }
	    else {
	        ss = -1. ;
	    }

	    while( diff_dent_x > 1.e-15 ) {

	        mm++ ;
		xtmp = xxp + ss * mm * dx ;

		rr = sqrt( xtmp*xtmp + yyp*yyp ) ;

		if ( xtmp == 0. ) {
		    if ( ss > 0. ) {
		        pp = M_PI/2. ;
		    }
		    else {
		        pp = 1.5*M_PI ;
		    }
		}
		else {
		    pp = acos( xtmp / rr ) ;
		}

		diff_dent_x = ss * dent(0).val_point(rr, M_PI/2., pp) ;

	    }
	    xxp += ss * (mm - 1) * dx ;

	}
    }

    xx = xxp ;
    yy = yyp ;

}
