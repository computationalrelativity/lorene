/*
 *  Methods Et_bin_bhns_extr::update_metric_extr_ks
 *  and Et_bin_bhns_extr::update_metric_extr_cf
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:17:52  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:52:55  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:08  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2005/02/28 23:16:37  k_taniguchi
 * Modification to include the case of the conformally flat background metric
 *
 * Revision 1.1  2004/11/30 20:51:32  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include <cmath>

// Lorene headers
#include "et_bin_bhns_extr.h"
#include "etoile.h"
#include "coord.h"
#include "unites.h"

          //-----------------------------------------------------------//
          //          No relaxation for a fixed BH background          //
          //-----------------------------------------------------------//

namespace Lorene {
void Et_bin_bhns_extr::update_metric_extr(const double& mass,
					  const double& sepa)
{

  using namespace Unites ;

    if (kerrschild) {

        // Computation of quantities coming from the companion (K-S BH)
        // ------------------------------------------------------------

        const Coord& xx = mp.x ;
	const Coord& yy = mp.y ;
	const Coord& zz = mp.z ;

	Tenseur r_bh(mp) ;
	r_bh.set_etat_qcq() ;
	r_bh.set() = pow( (xx+sepa)*(xx+sepa) + yy*yy + zz*zz, 0.5) ;
	r_bh.set_std_base() ;

	Tenseur xx_con(mp, 1, CON, ref_triad) ;
	xx_con.set_etat_qcq() ;
	xx_con.set(0) = xx + sepa ;
	xx_con.set(1) = yy ;
	xx_con.set(2) = zz ;
	xx_con.set_std_base() ;

	Tenseur xsr_con(mp, 1, CON, ref_triad) ;
	xsr_con = xx_con / r_bh ;
	xsr_con.set_std_base() ;

	Tenseur msr(mp) ;
	msr = ggrav * mass / r_bh ;
	msr.set_std_base() ;

	Tenseur lapse_bh(mp) ;
	lapse_bh = 1. / sqrt( 1.+2.*msr ) ;
	lapse_bh.set_std_base() ;

	logn_comp.set_etat_qcq() ;
	logn_comp.set() = log( lapse_bh() ) ;
	logn_comp.set_std_base() ;

	beta_comp.set_etat_qcq() ;
	beta_comp.set() = log( lapse_bh() ) ; 
	                             // conformal factor of KS-BH is unity
	beta_comp.set_std_base() ;

	shift_comp.set_etat_qcq() ;

	shift_comp.set(0) = -2.*lapse_bh()*lapse_bh()*msr()*xsr_con(0) ;
	shift_comp.set(1) = -2.*lapse_bh()*lapse_bh()*msr()*xsr_con(1) ;
	shift_comp.set(2) = -2.*lapse_bh()*lapse_bh()*msr()*xsr_con(2) ;

	shift_comp.set_std_base() ;
	shift_comp.set_triad( ref_triad ) ;

	// Lapse function N
	// ----------------

	nnn = exp( unsurc2 * logn_auto ) * lapse_bh ;

	nnn.set_std_base() ;

	// Conformal factor A^2
	// --------------------

	a_car = exp ( 2.*unsurc2*(beta_auto - logn_auto) ) ;

	a_car.set_std_base() ;

	// Shift vector N^i
	// ----------------

	shift = shift_auto + shift_comp ;

	// Derivative of metric coefficients
	// ----------------------------------

	// ... (d/dX,d/dY,d/dZ)(logn_auto) :
	d_logn_auto_regu = logn_auto_regu.gradient() ;    // (d/dx, d/dy, d/dz)
	d_logn_auto_regu.change_triad(ref_triad) ;   // -->  (d/dX, d/dY, d/dZ)

	if ( *(d_logn_auto_div.get_triad()) != ref_triad ) {

	    // Change the basis from spherical coordinate to Cartesian one
	    d_logn_auto_div.change_triad( mp.get_bvect_cart() ) ;

	    // Change the basis from mapping coordinate to absolute one
	    d_logn_auto_div.change_triad( ref_triad ) ;

	}

	d_logn_auto = d_logn_auto_regu + d_logn_auto_div ;

	// ... (d/dX,d/dY,d/dZ)(beta_auto) :
	d_beta_auto = beta_auto.gradient() ;    // (d/dx, d/dy, d/dz)
	d_beta_auto.change_triad(ref_triad) ;   // -->  (d/dX, d/dY, d/dZ)

	if (relativistic) {
	    // ... extrinsic curvature (tkij_auto and akcar_auto)
	    extrinsic_curv_extr(mass, sepa) ;
	}

	// The derived quantities are obsolete
	// -----------------------------------

	Etoile_bin::del_deriv() ;

    }
    else {

        // Computation of quantities coming from the companion (CF Sch. BH)
        // ----------------------------------------------------------------

        const Coord& xx = mp.x ;
	const Coord& yy = mp.y ;
	const Coord& zz = mp.z ;

	Tenseur r_bh(mp) ;
	r_bh.set_etat_qcq() ;
	r_bh.set() = pow( (xx+sepa)*(xx+sepa) + yy*yy + zz*zz, 0.5) ;
	r_bh.set_std_base() ;

	Tenseur msr(mp) ;
	msr = ggrav * mass / r_bh ;
	msr.set_std_base() ;

	Tenseur lapse_bh(mp) ;
	lapse_bh = (1.-0.5*msr) / (1.+0.5*msr) ;
	lapse_bh.set_std_base() ;

	logn_comp.set_etat_qcq() ;
	logn_comp.set() = log( lapse_bh() ) ;
	logn_comp.set_std_base() ;

	Tenseur lappsi(mp) ;
	lappsi = 1. - 0.25*msr*msr ;
	lappsi.set_std_base() ;

	beta_comp.set_etat_qcq() ;
	beta_comp.set() = log( lappsi() ) ;
	beta_comp.set_std_base() ;

	shift_comp.set_etat_qcq() ;

	shift_comp.set(0) = 0. ;
	shift_comp.set(1) = 0. ;
	shift_comp.set(2) = 0. ;

	shift_comp.set_std_base() ;
	shift_comp.set_triad( ref_triad ) ;

	// Lapse function N
	// ----------------

	nnn = exp( unsurc2 * logn_auto ) * lapse_bh ;

	nnn.set_std_base() ;

	// Conformal factor A^2
	// --------------------

	a_car = exp ( 2.*unsurc2*(beta_auto + beta_comp
				  - logn_auto - logn_comp) ) ;

	a_car.set_std_base() ;

	// Shift vector N^i
	// ----------------

	shift = shift_auto + shift_comp ;

	// Derivative of metric coefficients
	// ----------------------------------

	// ... (d/dX,d/dY,d/dZ)(logn_auto) :
	d_logn_auto_regu = logn_auto_regu.gradient() ;  // (d/dx, d/dy, d/dz)
	d_logn_auto_regu.change_triad(ref_triad) ; // -->  (d/dX, d/dY, d/dZ)

	if ( *(d_logn_auto_div.get_triad()) != ref_triad ) {

	    // Change the basis from spherical coordinate to Cartesian one
	    d_logn_auto_div.change_triad( mp.get_bvect_cart() ) ;

	    // Change the basis from mapping coordinate to absolute one
	    d_logn_auto_div.change_triad( ref_triad ) ;

	}

	d_logn_auto = d_logn_auto_regu + d_logn_auto_div ;

	// ... (d/dX,d/dY,d/dZ)(beta_auto) :
	d_beta_auto = beta_auto.gradient() ;    // (d/dx, d/dy, d/dz)
	d_beta_auto.change_triad(ref_triad) ;   // -->  (d/dX, d/dY, d/dZ)

	if (relativistic) {
	    // ... extrinsic curvature (tkij_auto and akcar_auto)
	    extrinsic_curv_extr(mass, sepa) ;
	}

	// The derived quantities are obsolete
	// -----------------------------------

	Etoile_bin::del_deriv() ;

    }

}
}
