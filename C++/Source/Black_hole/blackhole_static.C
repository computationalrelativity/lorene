/*
 *  Method of class Black_hole to set metric quantities to a spherical,
 *   static, analytic  solution
 *
 *    (see file blackhole.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005-2007 Keisuke Taniguchi
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

char blackhole_static_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2008/05/15 19:31:17  k_taniguchi
 * Change of some parameters.
 *
 * Revision 1.1  2007/06/22 01:20:50  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
//#include <math.h>

// Lorene headers
#include "blackhole.h"
#include "unites.h"
#include "utilitaires.h"

void Black_hole::static_bh(bool neumann, bool first) {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    double mass = ggrav * mass_bh ;

    Scalar rr(mp) ;
    rr = mp.r ;
    rr.std_spectral_base() ;

    //-------------------------------------//
    //          Metric quantities          //
    //-------------------------------------//

    Scalar st(mp) ;
    st = mp.sint ;
    st.std_spectral_base() ;
    Scalar ct(mp) ;
    ct = mp.cost ;
    ct.std_spectral_base() ;
    Scalar sp(mp) ;
    sp = mp.sinp ;
    sp.std_spectral_base() ;
    Scalar cp(mp) ;
    cp = mp.cosp ;
    cp.std_spectral_base() ;

    Vector ll(mp, CON, mp.get_bvect_cart()) ;
    ll.set_etat_qcq() ;
    ll.set(1) = st * cp ;
    ll.set(2) = st * sp ;
    ll.set(3) = ct ;
    ll.std_spectral_base() ;

    if (kerrschild) {

        // Lapconf function
        // ----------------
        lapconf = 1./sqrt(1.+2.*mass/rr) ;
	lapconf.annule_domain(0) ;
	lapconf.std_spectral_base() ;
	lapconf.raccord(1) ;

	// Conformal factor
	// ----------------
	confo = 1. ;
	confo.std_spectral_base() ;

	// Shift vector
	// ------------
	for (int i=1; i<=3; i++)
            shift.set(i) = 2. * mass * lapconf % lapconf % ll(i) / rr ;

	shift.annule_domain(0) ;
	shift.std_spectral_base() ;

    }
    else {  // Isotropic coordinates with the maximal slicing

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

        Scalar r_are(mp) ;
	r_are = r_coord(neumann, first) ;
	r_are.std_spectral_base() ;

        // Lapconf function
        // ----------------
	lapconf = sqrt(1. - 2.*mass/r_are/rr
		       + cc*cc*pow(mass/r_are/rr,4.)) * sqrt(r_are) ;

	lapconf.std_spectral_base() ;
	lapconf.annule_domain(0) ;
	lapconf.raccord(1) ;

        // Conformal factor
	// ----------------
	confo = sqrt(r_are) ;
	confo.std_spectral_base() ;
	confo.annule_domain(0) ;
	confo.raccord(1) ;

	// Lapse function
	// --------------
	lapse = lapconf / confo ;
	lapse.std_spectral_base() ;
	lapse.annule_domain(0) ;
	lapse.raccord(1) ;

        // Shift vector
	// ------------
	for (int i=1; i<=3; i++) {
	    shift.set(i) = mass * mass * cc * ll(i) / rr / rr
	      / pow(r_are,3.) ;
	}

	shift.std_spectral_base() ;

	for (int i=1; i<=3; i++) {
	    shift.set(i).annule_domain(0) ;
	    shift.set(i).raccord(1) ;
	}

    }

}
