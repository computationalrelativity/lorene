/*
 *  Methods of the class Et_bin_bhns_extr for computing hydro quantities
 *  in the Kerr-Schild background metric or in the conformally flat one
 *  with extreme mass ratio
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

char et_bin_bhns_extr_hydro_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:52:55  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2005/02/28 23:14:16  k_taniguchi
 * Modification to include the case of the conformally flat background metric
 *
 * Revision 1.1  2004/11/30 20:49:34  k_taniguchi
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

namespace Lorene {
void Et_bin_bhns_extr::hydro_euler_extr(const double& mass,
					const double& sepa) {

  using namespace Unites ;

    if (kerrschild) {

        int nz = mp.get_mg()->get_nzone() ;
	int nzm1 = nz - 1 ;

	//--------------------------------
	// Specific relativistic enthalpy		    ---> hhh
	//--------------------------------

	Tenseur hhh = exp(unsurc2 * ent) ;  // = 1 at the Newtonian limit
	hhh.set_std_base() ;

	//----------------------------------------
	// Lorentz factor between the co-orbiting	          ---> gam0
	// observer and the Eulerian one
	//----------------------------------------

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

	Tenseur gam0 = 1 / sqrt( 1 - unsurc2*tmp2 ) ;
	gam0.set_std_base() ;

	//--------------------------------------------
	// Lorentz factor and 3-velocity of the fluid
	//  with respect to the Eulerian observer
	//--------------------------------------------

	if (irrotational) {

	    d_psi.set_std_base() ;

	    Tenseur xx_con(mp, 1, CON, ref_triad) ;
	    xx_con.set_etat_qcq() ;
	    xx_con.set(0) = xx + sepa ;
	    xx_con.set(1) = yy ;
	    xx_con.set(2) = zz ;
	    xx_con.set_std_base() ;

	    Tenseur xsr_con(mp, 1, CON, ref_triad) ;
	    xsr_con = xx_con / r_bh ;
	    xsr_con.set_std_base() ;

	    Tenseur tmp3(mp) ;
	    tmp3.set_etat_qcq() ;
	    tmp3.set() = 0 ;
	    tmp3.set_std_base() ;

	    for (int i=0; i<3; i++)
	        tmp3.set() += xsr_con(i) % d_psi(i) ;

	    tmp3.set_std_base() ;

	    Tenseur tmp4 = -2.*msr % tmp3 % tmp3 / (1.+2.*msr) ;
	    tmp4.set_std_base() ;

	    for (int i=0; i<3; i++)
	        tmp4.set() += d_psi(i) % d_psi(i) ;

	    tmp4 = tmp4 / a_car ;

	    gam_euler = sqrt( 1 + unsurc2 * tmp4 / (hhh%hhh) ) ;

	    gam_euler.set_std_base() ;  // sets the standard spectral bases
	                                // for a scalar field

	    Tenseur vtmp1 = d_psi / ( hhh % gam_euler % a_car ) ;
	                // COV (a_car correction) -> CON
	    Tenseur vtmp2 = -2.* msr % tmp3 % xsr_con / (1.+2.*msr)
	      / ( hhh % gam_euler % a_car ) ;
	                // CON

	    // The assignment of u_euler is performed component by component
	    //  because u_euler is contravariant and d_psi is covariant
	    if (vtmp1.get_etat() == ETATZERO) {
	        u_euler.set_etat_zero() ;
	    }
	    else {
	        assert(vtmp1.get_etat() == ETATQCQ) ;
		u_euler.set_etat_qcq() ;
		for (int i=0; i<3; i++) {
		    u_euler.set(i) = vtmp1(i) + vtmp2(i) ;
		}
		u_euler.set_triad( *(vtmp1.get_triad()) ) ;
	    }

	    u_euler.set_std_base() ;

	}
	else {          // Rigid rotation
	                // --------------

	    gam_euler = gam0 ;
	    gam_euler.set_std_base() ;  // sets the standard spectral bases
	                                // for a scalar field

	    u_euler = - bsn ;

	}

	//--------------------------------------------------------
	// Energy density E with respect to the Eulerian observer
	//--------------------------------------------------------

	ener_euler = gam_euler % gam_euler % ( ener + press ) - press ;

	//------------------------------------------------------------------
	// Trace of the stress tensor with respect to the Eulerian observer
	//------------------------------------------------------------------

	Tenseur stmp1(mp) ;
	stmp1.set_etat_qcq() ;
	stmp1.set() = 0 ;
	stmp1.set_std_base() ;

	for (int i=0; i<3; i++)
	    stmp1.set() += xsr_cov(i) % u_euler(i) ;

	stmp1.set_std_base() ;

	Tenseur stmp2 = 2.*msr % stmp1 % stmp1 ;
	stmp2.set_std_base() ;

	for (int i=0; i<3; i++)
	    stmp2.set() += u_euler(i) % u_euler(i) ;

	stmp2 = a_car % stmp2 ;

	s_euler = 3 * press + ( ener_euler + press ) * stmp2 ;
	s_euler.set_std_base() ;

	//--------------------------------------
	// Lorentz factor between the fluid and		---> gam
	//	co-orbiting observers
	//--------------------------------------

	Tenseur gtmp = 2.*msr % tmp1 % stmp1 ;  //<- bsn^i = - U_0^i
	gtmp.set_std_base() ;

	for (int i=0; i<3; i++)
	    gtmp.set() += bsn(i) % u_euler(i) ; //<- bsn^i = - U_0^i

	gtmp = a_car % gtmp ;

	Tenseur tmp = ( 1+unsurc2*gtmp ) ; //<- (minus->plus) because of U_0^i
	tmp.set_std_base() ;
	Tenseur gam = gam0 % gam_euler % tmp ;

	//--------------------------------------------
	// Spatial projection of the fluid 3-velocity
	//  with respect to the co-orbiting observer
	//--------------------------------------------

	wit_w = gam_euler / gam % u_euler + gam0 % bsn ;

	wit_w.set_std_base() ;  // set the bases for spectral expansions

	wit_w.annule(nzm1) ;	// zero in the ZEC

	//-----------------------------------------
	// Logarithm of the Lorentz factor between
	//	the fluid and co-orbiting observers
	//-----------------------------------------

	if (relativistic) {

	    loggam = log( gam ) ;
	    loggam.set_std_base() ;   // set the bases for spectral expansions

	}
	else {
	    cout << "BH-NS binary systems should be relativistic !!!" << endl ;
	    abort() ;
	}

	//-------------------------------------------------
	// Velocity fields set to zero in external domains
	//-------------------------------------------------

	loggam.annule(nzm1) ;	    // zero in the ZEC only

	wit_w.annule(nzm1) ;	    // zero outside the star

	u_euler.annule(nzm1) ;	    // zero outside the star

	if (loggam.get_etat() != ETATZERO) {
	    (loggam.set()).set_dzpuis(0) ;
	}

	if (!irrotational) {

	    gam = 1 ;
	    loggam = 0 ;
	    wit_w = 0 ;

	}

	// The derived quantities are obsolete
	// -----------------------------------

	Etoile_bin::del_deriv() ;

    }
    else {

        int nz = mp.get_mg()->get_nzone() ;
	int nzm1 = nz - 1 ;

	//--------------------------------
	// Specific relativistic enthalpy		    ---> hhh
	//--------------------------------

	Tenseur hhh = exp(unsurc2 * ent) ;  // = 1 at the Newtonian limit
	hhh.set_std_base() ;

	//----------------------------------------
	// Lorentz factor between the co-orbiting	          ---> gam0
	// observer and the Eulerian one
	//----------------------------------------

	Tenseur gam0 = 1 / sqrt( 1 - unsurc2*sprod(bsn, bsn) ) ;
	gam0.set_std_base() ;

	//--------------------------------------------
	// Lorentz factor and 3-velocity of the fluid
	//  with respect to the Eulerian observer
	//--------------------------------------------

	if (irrotational) {

	    d_psi.set_std_base() ;

	    gam_euler = sqrt( 1 + unsurc2 * sprod(d_psi, d_psi)
			      / (hhh%hhh) ) ;

	    gam_euler.set_std_base() ;  // sets the standard spectral bases
	                                // for a scalar field

	    Tenseur vtmp = d_psi / ( hhh % gam_euler % a_car ) ;
	                // COV (a_car correction) -> CON

	    // The assignment of u_euler is performed component by component
	    //  because u_euler is contravariant and d_psi is covariant
	    if (vtmp.get_etat() == ETATZERO) {
	        u_euler.set_etat_zero() ;
	    }
	    else {
	        assert(vtmp.get_etat() == ETATQCQ) ;
		u_euler.set_etat_qcq() ;
		for (int i=0; i<3; i++) {
		    u_euler.set(i) = vtmp(i) ;
		}
		u_euler.set_triad( *(vtmp.get_triad()) ) ;
	    }

	    u_euler.set_std_base() ;

	}
	else {          // Rigid rotation
	                // --------------

	    gam_euler = gam0 ;
	    gam_euler.set_std_base() ;  // sets the standard spectral bases
	                                // for a scalar field

	    u_euler = - bsn ;

	}

	//--------------------------------------------------------
	// Energy density E with respect to the Eulerian observer
	//--------------------------------------------------------

	ener_euler = gam_euler % gam_euler % ( ener + press ) - press ;

	//------------------------------------------------------------------
	// Trace of the stress tensor with respect to the Eulerian observer
	//------------------------------------------------------------------

	s_euler = 3 * press + ( ener_euler + press )
	  * sprod(u_euler, u_euler) ;
	s_euler.set_std_base() ;

	//--------------------------------------
	// Lorentz factor between the fluid and		---> gam
	//	co-orbiting observers
	//--------------------------------------

	Tenseur tmp = ( 1+unsurc2*sprod(bsn, u_euler) ) ;
	                                //<- (minus->plus) because of U_0^i
	tmp.set_std_base() ;
	Tenseur gam = gam0 % gam_euler % tmp ;

	//--------------------------------------------
	// Spatial projection of the fluid 3-velocity
	//  with respect to the co-orbiting observer
	//--------------------------------------------

	wit_w = gam_euler / gam % u_euler + gam0 % bsn ;

	wit_w.set_std_base() ;  // set the bases for spectral expansions

	wit_w.annule(nzm1) ;	// zero in the ZEC

	//-----------------------------------------
	// Logarithm of the Lorentz factor between
	//	the fluid and co-orbiting observers
	//-----------------------------------------

	if (relativistic) {

	    loggam = log( gam ) ;
	    loggam.set_std_base() ;  // set the bases for spectral expansions

	}
	else {
	    cout << "BH-NS binary systems should be relativistic !!!" << endl ;
	    abort() ;
	}

	//-------------------------------------------------
	// Velocity fields set to zero in external domains
	//-------------------------------------------------

	loggam.annule(nzm1) ;	    // zero in the ZEC only

	wit_w.annule(nzm1) ;	    // zero outside the star

	u_euler.annule(nzm1) ;	    // zero outside the star

	if (loggam.get_etat() != ETATZERO) {
	    (loggam.set()).set_dzpuis(0) ;
	}

	if (!irrotational) {

	    gam = 1 ;
	    loggam = 0 ;
	    wit_w = 0 ;

	}

	// The derived quantities are obsolete
	// -----------------------------------

	Etoile_bin::del_deriv() ;

    }

}
}
