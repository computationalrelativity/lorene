/*
 *  Methods of class Black_hole to compute the radius of the apparent
 *  horizon in isotropic coordinates
 *
 *    (see file blackhole.h for documentation).
 *
 */

/*
 *   Copyright (c) 2006-2007 Keisuke Taniguchi
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

char blackhole_rah_iso_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2008/05/15 19:30:35  k_taniguchi
 * Change of some parameters.
 *
 * Revision 1.1  2007/06/22 01:20:13  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
#include <math.h>

// Lorene headers
#include "blackhole.h"
#include "utilitaires.h"

// Local function
double ff(double, const double) ;

          //----------------------------------------------------------//
          //     Radius of the apparent horizon (excised surface)     //
          //----------------------------------------------------------//

double Black_hole::rah_iso(bool neumann, bool first) const {

    // Sets C/M^2 for each case of the lapse boundary condition
    // --------------------------------------------------------
    double cc ;

    if (neumann) {  // Neumann boundary condition
        if (first) {  // First condition
	  // d(\alpha \psi)/dr = 0
	  // ---------------------
	  cc = 2. * (sqrt(13.) - 1.) / 3. ;
	}
	else {  // Second condition
	  // d(\alpha \psi)/dr = (\alpha \psi)/(2 rah)
	  // -----------------------------------------
	  cc = 4. / 3. ;
	}
    }
    else {  // Dirichlet boundary condition
       if (first) {  // First condition
	 // (\alpha \psi) = 1/2
	 // -------------------
	 cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
	 abort() ;
       }
       else {  // Second condition
	 // (\alpha \psi) = 1/sqrt(2.) \psi_KS
	 // ----------------------------------
	 cout << "!!!!! WARNING: Not yet prepared !!!!!" << endl ;
	 abort() ;
       }
    }

    int nn = 1000 ;
    double hh = 1./double(nn) ;
    double integ = 0. ;
    double rah ;  // rah [M]

    int mm ;
    double x1, x2, x3, x4, x5 ;

    // Boole's Rule (Newton-Cotes Integral) for integration
    // ----------------------------------------------------

    assert(nn%4 == 0) ;
    mm = nn/4 ;

    for (int i=0; i<mm; i++) {

        x1 = hh * double(4*i) ;
	x2 = hh * double(4*i+1) ;
	x3 = hh * double(4*i+2) ;
	x4 = hh * double(4*i+3) ;
	x5 = hh * double(4*i+4) ;

	integ += (hh/45.) * (14.*ff(x1,cc) + 64.*ff(x2,cc)
			     + 24.*ff(x3,cc) + 64.*ff(x4,cc)
			     + 14.*ff(x5,cc)) ;

    }

    rah = 2. * exp(integ) ;  // rah : normalized by M

    return rah ;

}

//*****************************************************************

double ff(double xx, const double cc) {

    double tcc2 = cc*cc/16. ;
    double tmp = sqrt(1. - xx + tcc2*pow(xx, 4.)) ;

    double resu = (-1. + tcc2 * pow(xx, 3.)) / tmp / (1. + tmp) ;

    return resu ;

}
