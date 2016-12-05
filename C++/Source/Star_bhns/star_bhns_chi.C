/*
 *  Method of class Star_bhns to compute a sensitve indicator of
 *   the mass-shedding and quantities related to the indicator
 *
 *    (see file star_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2006 Keisuke Taniguchi
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
 * Revision 1.4  2016/12/05 16:18:16  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:40  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:16  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2007/06/22 01:30:27  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
#include <cmath>

// Lorene headers
#include "star_bhns.h"
#include "param.h"
#include "tbl.h"
#include "utilitaires.h"

namespace Lorene {
double Star_bhns::chi_rp(double radius, double phi) {

    const Scalar& dent = ent.dsdr() ;

    double dent_pole = dent.val_point(ray_pole(), 0., 0.) ;
    double dent_eq = dent.val_point(radius, M_PI/2., phi) ;

    double chi = fabs( dent_eq / dent_pole ) ;

    return chi ;

}

double Star_bhns::radius_p(double phi) {

    double rad = mp.val_r(nzet-1, 1., M_PI/2., phi) ;
    // We assume that the stellar surface is fitted to the domain (nzet-1)

    return rad ;

}

double Star_bhns::phi_min() {

    const Mg3d* mg = mp.get_mg() ;
    int np = mg->get_np(0) ;
    int nps2 = np/2 ;

    Tbl phi(nps2+1) ;
    phi.set_etat_qcq() ;

    for (int i=0; i<=nps2; i++) {
        phi.set(i) = 2.*M_PI*i/np + 0.5*M_PI ;
    }

    Tbl chi(nps2+1) ;
    chi.set_etat_qcq() ;

    for (int i=0; i<=nps2; i++) {

        double phi_i = phi_local_min( phi(i) ) ;
	double rad_i = radius_p( phi_i ) ;

	chi.set(i) = chi_rp(rad_i, phi_i) ;

    }

    for (int i=0; i<=nps2; i++) {

        cout.precision(16) ;
	cout << "chi(" << i << ") = " << chi(i)
	     << "  phi = " << phi_local_min( phi(i) ) / M_PI
	     << " [M_PI]" << endl ;
	cout.precision(4) ;

    }

    double chi_ini = chi(0) ; // Initialization
    double delta_chi ;
    int jj = 0 ; // Initialization

    for (int i=0; i<nps2; i++) {

        if ( chi(i+1) < 1.e-12 )
	    chi.set(i+1) = 1. ;

	delta_chi = chi_ini - chi(i+1) ;

	if ( delta_chi > 0. ) {

	    chi_ini = chi(i+1) ;
	    jj = i+1 ;

	}

    }

    double phi_glob_min = phi_local_min( phi(jj) ) ;

    return phi_glob_min ;

}


double Star_bhns::phi_local_min(double phi_ini) {

    int mm ;                // Number of steps to the phi direction
    double ppp = phi_ini ;  // Initial position of the phi coordinate
    double diff ;           // Difference between two succesive chi
    double dp = M_PI/2. ;   // Step interval to the phi direction
    double ptmp ;

    double rad1, rad2 ;

    double init_check = chi_rp(radius_p(phi_ini), phi_ini)
      - chi_rp(radius_p(phi_ini+1.e-10*dp), phi_ini+1.e-10*dp) ;

    if ( init_check >= 0. ) {

        while ( dp > 1.e-15 ) {

	    diff = 1. ;
	    mm = 0 ;
	    dp = 0.1 * dp ;

	    while ( diff > 0. && (ppp+mm*dp) < 2.*M_PI ) {

	        mm++ ;
		ptmp = ppp + mm * dp ;

		rad1 = radius_p(ptmp-dp) ;
		rad2 = radius_p(ptmp) ;

		diff = chi_rp(rad1, ptmp-dp) - chi_rp(rad2, ptmp) ;

	    }

	    ppp += (mm - 2) * dp ;

	}

	if ( (ppp+2.*dp) >= 2.*M_PI ) {

	    cout << "No minimum for phi > " << phi_ini / M_PI
		 << " [M_PI]" << endl ;

	}

    }
    else {

        while ( dp > 1.e-15 ) {

	    diff = 1. ;
	    mm = 0 ;
	    dp = 0.1 * dp ;

	    while ( diff > 0. && (ppp-mm*dp) > 0. ) {

	        mm++ ;
		ptmp = ppp - mm * dp ;

		rad1 = radius_p(ptmp+dp) ;
		rad2 = radius_p(ptmp) ;

		diff = chi_rp(rad1, ptmp+dp) - chi_rp(rad2, ptmp) ;

	    }

	    ppp -= (mm - 2) * dp ;

	}

	if ( (ppp-2.*dp) < 0. ) {

	    cout << "No minimum for phi < " << phi_ini / M_PI
		 << " [M_PI]" << endl ;

	}

    }

    return ppp ;

}
}
