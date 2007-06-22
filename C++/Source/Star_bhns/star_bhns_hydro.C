/*
 *  Method of class Star_bhns to compute hydro quantities
 *
 *    (see file star_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Keisuke Taniguchi
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

char star_bhns_hydro_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:31:42  k_taniguchi
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
#include "utilitaires.h"
#include "unites.h"

void Star_bhns::hydro_euler_bhns(bool kerrschild, const double& mass_bh,
				 const double& sepa) {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    int nz = mp.get_mg()->get_nzone() ;
    int nzm1 = nz - 1 ;

    //--------------------------------
    // Specific relativistic enthalpy          ---> hhh
    //--------------------------------

    Scalar hhh = exp(ent) ;  // = 1 at the Newtonian limit
    hhh.std_spectral_base() ;

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

	//--------------------------------------------
	// Lorentz factor and 3-velocity of the fluid
	//  with respect to the Eulerian observer
	//--------------------------------------------

	if (irrotational) {

	    d_psi.std_spectral_base() ;

	    Scalar dpsi2(mp) ;
	    dpsi2 = d_psi(1)%d_psi(1) + d_psi(2)%d_psi(2)
	      + d_psi(3)%d_psi(3) ;
	    dpsi2.std_spectral_base() ;

	    Scalar lldpsi(mp) ;
	    lldpsi = ll(1)%d_psi(1) + ll(2)%d_psi(2) + ll(3)%d_psi(3) ;
	    lldpsi.std_spectral_base() ;

	    Scalar lap_bh2(mp) ;
	    lap_bh2 = 1. / (1.+2.*msr) ;
	    lap_bh2.std_spectral_base() ;

	    Scalar tmp1(mp) ;
	    tmp1 = 2. * msr % lap_bh2 % lldpsi % lldpsi ;
	    tmp1.std_spectral_base() ;

	    gam_euler = sqrt( 1.+(dpsi2-tmp1)/(hhh%hhh)/psi4 ) ;
	    gam_euler.std_spectral_base() ;

	    u_euler.set_etat_qcq() ;
	    assert( u_euler.get_triad() == bsn.get_triad() ) ;

	    for (int i=1; i<=3; i++)
	      u_euler.set(i) = (d_psi(i)-2.*msr%lap_bh2%lldpsi%ll(i))
		/ (hhh % gam_euler) / psi4 ;

	    u_euler.std_spectral_base() ;

	}
	else {

	    // Rigid rotation
	    // --------------

	    gam_euler = gam0 ;
	    gam_euler.std_spectral_base() ;
	    u_euler = bsn ;

	}

	//--------------------------------------------------------
	// Energy density E with respect to the Eulerian observer
	// See Eq (53) from Gourgoulhon et al. (2001)
	//--------------------------------------------------------

	ener_euler = gam_euler % gam_euler % (ener + press) - press ;
	ener_euler.std_spectral_base() ;

	//---------------------------------------------------
	// Trace of the stress-energy tensor with respect to
	// the Eulerian observer
	// See Eq (54) from Gourgoulhon et al. (2001)
	//---------------------------------------------------

	Scalar ueuler2(mp) ;
	ueuler2 = u_euler(1)%u_euler(1) + u_euler(2)%u_euler(2)
	  + u_euler(3)%u_euler(3) ;
	ueuler2.std_spectral_base() ;

	Scalar llueuler(mp) ;
	llueuler = ll(1)%u_euler(1) + ll(2)%u_euler(2) + ll(3)%u_euler(3) ;
	llueuler.std_spectral_base() ;

	s_euler = 3. * press + (ener_euler + press) * psi4
	  * (ueuler2 + 2.*msr%llueuler%llueuler) ;
	s_euler.std_spectral_base() ;

	//--------------------------------------------
	// Lorentz factor between the fluid and        ---> gam
	//  co-orbiting observers
	// See Eq (58) from Gourgoulhon et al. (2001)
	//--------------------------------------------

	if (irrotational) {

	    Scalar bsnueuler(mp) ;
	    bsnueuler = bsn(1)%u_euler(1) + bsn(2)%u_euler(2)
	      + bsn(3)%u_euler(3) ;
	    bsnueuler.std_spectral_base() ;

	    Scalar llbsn(mp) ;
	    llbsn = ll(1)%bsn(1) + ll(2)%bsn(2) + ll(3)%bsn(3) ;
	    llbsn.std_spectral_base() ;

	    Scalar tmp2(mp) ;
	    tmp2 = 1. - psi4 * (bsnueuler + 2.*msr%llueuler%llbsn) ;
	    tmp2.std_spectral_base() ;

	    gam = gam0 % gam_euler % tmp2 ;
	    gam.std_spectral_base() ;

	    //--------------------------------------------
	    // Spetial projection of the fluid 3-velocity
	    // with respect to the co-orbiting observer
	    //--------------------------------------------

	    wit_w = gam_euler / gam * u_euler - gam0 * bsn ;
	    wit_w.std_spectral_base() ;
	    wit_w.annule_domain(nzm1) ;

	    //-----------------------------------------
	    // Logarithm of the Lorentz factor between
	    // the fluid and co-orbiting observer
	    //-----------------------------------------

	    loggam = log( gam ) ;
	    loggam.std_spectral_base() ;

	    //-------------------------------------------------
	    // Velocity fields set to zero in external domains
	    //-------------------------------------------------

	    loggam.annule_domain(nzm1) ;
	    wit_w.annule_domain(nzm1) ;
	    u_euler.annule_domain(nzm1) ;

	    loggam.set_dzpuis(0) ;

	}
	else {  // Rigid rotation

	    gam = 1. ;
	    gam.std_spectral_base() ;
	    loggam = 0. ;
	    wit_w.set_etat_zero() ;

	}

    }  // End of Kerr-Schild
    else {  // Isotropic coordinates with the maximal slicing

	//--------------------------------------------
	// Lorentz factor and 3-velocity of the fluid
	//  with respect to the Eulerian observer
	//--------------------------------------------

	if (irrotational) {

	    d_psi.std_spectral_base() ;

	    Scalar dpsi2(mp) ;
	    dpsi2 = d_psi(1)%d_psi(1) + d_psi(2)%d_psi(2)
	      + d_psi(3)%d_psi(3) ;
	    dpsi2.std_spectral_base() ;

	    gam_euler = sqrt( 1. + dpsi2/(hhh%hhh)/psi4 ) ;
	    gam_euler.std_spectral_base() ;

	    u_euler.set_etat_qcq() ;
	    for (int i=1; i<=3; i++)
	      u_euler.set(i) = d_psi(i)/(hhh%gam_euler)/psi4 ;

	    u_euler.std_spectral_base() ;

	}
	else {

	    // Rigid rotation
	    // --------------

	    gam_euler = gam0 ;
	    gam_euler.std_spectral_base() ;
	    u_euler = bsn ;

	}

	//--------------------------------------------------------
	// Energy density E with respect to the Eulerian observer
	// See Eq (53) from Gourgoulhon et al. (2001)
	//--------------------------------------------------------

	ener_euler = gam_euler % gam_euler % (ener + press) - press ;
	ener_euler.std_spectral_base() ;

	//---------------------------------------------------
	// Trace of the stress-energy tensor with respect to
	// the Eulerian observer
	// See Eq (54) from Gourgoulhon et al. (2001)
	//---------------------------------------------------

	Scalar ueuler2(mp) ;
	ueuler2 = u_euler(1)%u_euler(1) + u_euler(2)%u_euler(2)
	  + u_euler(3)%u_euler(3) ;
	ueuler2.std_spectral_base() ;

	s_euler = 3.*press + (ener_euler+press)*psi4*ueuler2 ;
	s_euler.std_spectral_base() ;


	//--------------------------------------------
	// Lorentz factor between the fluid and        ---> gam
	//  co-orbiting observers
	// See Eq (58) from Gourgoulhon et al. (2001)
	//--------------------------------------------

	if (irrotational) {

	    Scalar bsnueuler(mp) ;
	    bsnueuler = bsn(1)%u_euler(1) + bsn(2)%u_euler(2)
	      + bsn(3)%u_euler(3) ;
	    bsnueuler.std_spectral_base() ;

	    Scalar tmp3(mp) ;
	    tmp3 = 1. - psi4 % bsnueuler ;
	    tmp3.std_spectral_base() ;

	    gam = gam0 % gam_euler % tmp3 ;
	    gam.std_spectral_base() ;

	    //--------------------------------------------
	    // Spetial projection of the fluid 3-velocity
	    // with respect to the co-orbiting observer
	    //--------------------------------------------

	    wit_w = gam_euler / gam * u_euler - gam0 * bsn ;
	    wit_w.std_spectral_base() ;
	    wit_w.annule_domain(nzm1) ;

	    //-----------------------------------------
	    // Logarithm of the Lorentz factor between
	    // the fluid and co-orbiting observer
	    //-----------------------------------------

	    loggam = log( gam ) ;
	    loggam.std_spectral_base() ;

	    //-------------------------------------------------
	    // Velocity fields set to zero in external domains
	    //-------------------------------------------------

	    loggam.annule_domain(nzm1) ;
	    wit_w.annule_domain(nzm1) ;
	    u_euler.annule_domain(nzm1) ;

	    loggam.set_dzpuis(0) ;

	}
	else {  // Rigid rotation

	    gam = 1. ;
	    gam.std_spectral_base() ;
	    loggam = 0. ;
	    wit_w.set_etat_zero() ;

	}

    }   // End of Isotropic coordinates

    // The derived quantities are obsolete
    // -----------------------------------

    del_deriv() ;

}
