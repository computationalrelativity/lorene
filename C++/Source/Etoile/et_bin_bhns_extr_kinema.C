/*
 *  Method Et_bin_bhns_extr::kinematics_extr_ks
 *  and Et_bin_bhns_extr::kinematics_extr_cf
 *
 *    (see file et_bin_bhns_extr.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004-2005 Keisuke Taniguchi
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

char et_bin_bhns_extr_kinema_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2005/02/28 23:15:09  k_taniguchi
 * Modification to include the case of the conformally flat background metric
 *
 * Revision 1.1  2004/11/30 20:49:58  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "et_bin_bhns_extr.h"
#include "etoile.h"
#include "coord.h"
#include "unites.h"

void Et_bin_bhns_extr::kinematics_extr(double omega, const double& mass,
				       const double& sepa) {

  using namespace Unites ;

    if (kerrschild) {

        int nz = mp.get_mg()->get_nzone() ;
	int nzm1 = nz - 1 ;

	// --------------------
	// Computation of B^i/N
	// --------------------

	//  1/ Computation of  - omega m^i

	const Coord& xa = mp.xa ;
	const Coord& ya = mp.ya ;

	bsn.set_etat_qcq() ;

	bsn.set(0) = omega * ya ;
	bsn.set(1) = - omega * xa ;
	bsn.set(2) = 0 ;

	bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC

	//	2/ Addition of shift and division by lapse

	bsn = ( bsn + shift ) / nnn ;

	bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC
	bsn.set_std_base() ;   // set the bases for spectral expansions

	//-------------------------
	// Centrifugal potentatial
	//-------------------------

	const Coord& xx = mp.x ;
	const Coord& yy = mp.y ;
	const Coord& zz = mp.z ;

	Tenseur r_bh(mp) ;
	r_bh.set_etat_qcq() ;
	r_bh.set() = pow( (xx+sepa)*(xx+sepa) + yy*yy + zz*zz, 0.5) ;
	r_bh.set_std_base() ;

	Tenseur xx_cov(mp, 1, COV, ref_triad) ;
	xx_cov.set_etat_qcq() ;
	xx_cov.set(0) = xx + sepa ;
	xx_cov.set(1) = yy ;
	xx_cov.set(2) = zz ;
	xx_cov.set_std_base() ;

	Tenseur xsr_cov(mp, 1, COV, ref_triad) ;
	xsr_cov = xx_cov / r_bh ;
	xsr_cov.set_std_base() ;

	Tenseur msr(mp) ;
	msr = ggrav * mass / r_bh ;
	msr.set_std_base() ;

	if (relativistic) {

	    // Lorentz factor between the co-orbiting observer and
	    // the Eulerian one

	    Tenseur tmp1(mp) ;
	    tmp1.set_etat_qcq() ;
	    tmp1.set() = 0 ;
	    tmp1.set_std_base() ;

	    for (int i=0; i<3; i++)
	        tmp1.set() += xsr_cov(i) % bsn(i) ;

	    tmp1.set_std_base() ;

	    Tenseur tmp2 = 2.*msr % tmp1 % tmp1 ;
	    tmp2.set_std_base() ;

	    for (int i=0; i<3; i++)
	        tmp2.set() += bsn(i) % bsn(i) ;

	    tmp2 = a_car % tmp2 ;

	    Tenseur gam0 = 1 / sqrt( 1 - tmp2 ) ;

	    pot_centri = - log( gam0 ) ;

	}
	else {
	    cout << "BH-NS binary system should be relativistic !!!" << endl ;
	    abort() ;
	}

	pot_centri.annule(nzm1, nzm1) ;	// set to zero in the external domain
	pot_centri.set_std_base() ;   // set the bases for spectral expansions

	// The derived quantities are obsolete
	// -----------------------------------

	Etoile_bin::del_deriv() ;

    }
    else {

        int nz = mp.get_mg()->get_nzone() ;
	int nzm1 = nz - 1 ;

	// --------------------
	// Computation of B^i/N
	// --------------------

	//  1/ Computation of  - omega m^i

	const Coord& xa = mp.xa ;
	const Coord& ya = mp.ya ;

	bsn.set_etat_qcq() ;

	bsn.set(0) = omega * ya ;
	bsn.set(1) = - omega * xa ;
	bsn.set(2) = 0 ;

	bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC

	//	2/ Addition of shift and division by lapse

	bsn = ( bsn + shift ) / nnn ;

	bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC
	bsn.set_std_base() ;   // set the bases for spectral expansions

	//-------------------------
	// Centrifugal potentatial
	//-------------------------

	if (relativistic) {

	    // Lorentz factor between the co-orbiting observer and
	    // the Eulerian one

	    Tenseur tmp(mp) ;
	    tmp.set_etat_qcq() ;
	    tmp.set() = 0. ;
	    tmp.set_std_base() ;

	    for (int i=0; i<3; i++)
	        tmp.set() += bsn(i) % bsn(i) ;

	    tmp = a_car % tmp ;

	    Tenseur gam0 = 1 / sqrt( 1 - tmp ) ;

	    pot_centri = - log( gam0 ) ;

	}
	else {
	    cout << "BH-NS binary system should be relativistic !!!" << endl ;
	    abort() ;
	}

	pot_centri.annule(nzm1, nzm1) ;	// set to zero in the external domain
	pot_centri.set_std_base() ;   // set the bases for spectral expansions

	// The derived quantities are obsolete
	// -----------------------------------

	Etoile_bin::del_deriv() ;

    }

}
