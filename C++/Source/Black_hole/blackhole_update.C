/*
 *  Methods of class Black_hole to compute the background metric quantities
 *
 *    (see file blackhole.h for documentation).
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

char blackhole_update_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:21:07  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
//#include <>

// Lorene headers
#include "blackhole.h"
#include "utilitaires.h"
#include "param.h"
#include "unites.h"

void Black_hole::update_metric() {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (kerrschild) {

        double mass = ggrav * mass_bh ;

        Scalar rr(mp) ;
	rr = mp.r ;
	rr.std_spectral_base() ;
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

	Vector ll(mp, COV, mp.get_bvect_cart()) ;
	ll.set_etat_qcq() ;
	ll.set(1) = st * cp ;
	ll.set(2) = st * sp ;
	ll.set(3) = ct ;
	ll.std_spectral_base() ;

	lapse_bh = 1. / sqrt(1.+2.*mass/rr) ;
	lapse_bh.std_spectral_base() ;
	lapse_bh.annule_domain(0) ;
	lapse_bh.raccord(1) ;

	for (int i=1; i<=3; i++) {
	    shift_bh.set(i) = 2. * lapse_bh * lapse_bh * mass * ll(i) / rr ;
	}
	shift_bh.std_spectral_base() ;

    }
    else {

        cout << "!!! WARNING : Black_hole::update_metric() is"
	     << "              for the Kerr-Schild case !!!" << endl ;
	abort() ;

    }

    // Deletes the derived quantities
    // ------------------------------

    del_deriv() ;

}
