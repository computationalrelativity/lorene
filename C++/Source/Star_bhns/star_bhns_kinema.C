/*
 *  Method of class Star_bhns to compute kinematic quantities
 *
 *    (see file star_bhns.h for documentation).
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

char star_bhns_kinema_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2008/05/15 19:16:06  k_taniguchi
 * Change of a parameter.
 *
 * Revision 1.1  2007/06/22 01:32:00  k_taniguchi
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
#include "star_bhns.h"
#include "unites.h"

void Star_bhns::kinema_bhns(bool kerrschild, const double& mass_bh,
			    const double& sepa, double omega,
			    double x_rot, double y_rot) {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    int nz = mp.get_mg()->get_nzone() ;
    int nzm1 = nz - 1 ;

    //----------------------
    // Computation of B^i/N
    //----------------------

    // 1/ Computation of omega m^i

    const Coord& xa = mp.xa ;
    const Coord& ya = mp.ya ;

    //    bsn.change_triad(mp.get_bvect_cart()) ;

    if (fabs(mp.get_rot_phi()) < 1.e-10) {

        bsn.set(1) = - omega * (ya - y_rot) ;
	bsn.set(2) = omega * (xa - x_rot) ;
	bsn.set(3) = 0. ;

    }
    else {

        bsn.set(1) = omega * (ya - y_rot) ;
	bsn.set(2) = - omega * (xa - x_rot) ;
	bsn.set(3) = 0. ;

    }

    bsn.std_spectral_base() ;
    bsn.annule_domain(nzm1) ;

    // 2/ Addition of shift_tot and division by lapse

    for (int i=1; i<=3; i++) {
        bsn.set(i) = confo_tot * ( bsn(i) + shift_tot(i) ) / lapconf_tot ;
    }
    bsn.std_spectral_base() ;
    bsn.annule_domain(nzm1) ;

    //-----------------------------------------------------
    // Lorentz factor between the co-orbiting               ---> gam0
    // observer and the Eulerian one
    // See Eq (23) and (24) from Gourgoulhon et al. (2001)
    //-----------------------------------------------------

    Scalar bsn2(mp) ;
    bsn2 = bsn(1)%bsn(1) + bsn(2)%bsn(2) + bsn(3)%bsn(3) ;
    bsn2.std_spectral_base() ;

    if (kerrschild) {

        double mass = ggrav * mass_bh ;

	Scalar xx(mp) ;
	xx = mp.x ;
	xx.std_spectral_base() ;
	Scalar yy(mp) ;
	yy = mp.y ;
	yy.std_spectral_base() ;
	Scalar zz(mp) ;
	zz = mp.z ;
	zz.std_spectral_base() ;

	double yns = mp.get_ori_y() ;

	Scalar rbh(mp) ;
	rbh = sqrt( (xx+sepa)*(xx+sepa) + (yy+yns)*(yy+yns) + zz*zz ) ;
	rbh.std_spectral_base() ;

	Vector ll(mp, CON, mp.get_bvect_cart()) ;
	ll.set_etat_qcq() ;
	ll.set(1) = (xx+sepa) / rbh ;
	ll.set(2) = (yy+yns) / rbh ;
	ll.set(3) = zz / rbh ;
	ll.std_spectral_base() ;

	Scalar msr(mp) ;
	msr = mass / rbh ;
	msr.std_spectral_base() ;

	Scalar llbsn(mp) ;
	llbsn = ll(1)%bsn(1) + ll(2)%bsn(2) + ll(3)%bsn(3) ;
	llbsn.std_spectral_base() ;

	Scalar tmp1(mp) ;
	tmp1 = 2. * msr % llbsn % llbsn ;
	tmp1.std_spectral_base() ;

	gam0 = 1. / sqrt(1. - psi4*(bsn2+tmp1)) ;
	gam0.std_spectral_base() ;

    }
    else { // Isotropic coordinates with the maximal slicing

        gam0 = 1. / sqrt(1. - psi4%bsn2) ;
	gam0.std_spectral_base() ;

    }

    //-----------------------
    // Centrifugal potential
    //-----------------------

    pot_centri = - log( gam0 ) ;
    pot_centri.annule_domain(nzm1) ;
    pot_centri.std_spectral_base() ;

    // The derived quantities are obsolete
    // -----------------------------------

    del_deriv() ;

}
