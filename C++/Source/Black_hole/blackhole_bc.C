/*
 *  Methods of class Black_hole to compute the inner boundary condition
 *  at the excised surface
 *
 *    (see file blackhole.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005-2006 Keisuk Taniguchi
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

char blackhole_bc_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:18:23  k_taniguchi
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
#include "valeur.h"
#include "grilles.h"
#include "unites.h"

                    //----------------------------------//
                    //     Inner boundary condition     //
                    //----------------------------------//

const Valeur Black_hole::bc_lapse(bool neumann, bool first) const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    const Mg3d* mg_angu = mg->get_angu() ;
    Valeur bc(mg_angu) ;

    if (kerrschild) {

        if (neumann) {

	    if (first) {

	        Scalar rr(mp) ;
		rr = mp.r ;
		rr.std_spectral_base() ;

		int nt = mg->get_nt(0) ;
		int np = mg->get_np(0) ;

		Scalar tmp(mp) ;

	        tmp = - pow(lapse_bh,3.) * ggrav * mass_bh / rr / rr ;
		// dlapse/dr = 0

		bc = 1. ;
		for (int j=0; j<nt; j++) {
		    for (int k=0; k<np; k++) {
		        bc.set(0,k,j,0) = tmp.val_grid_point(1,k,j,0) ;
		    }
		}
	    }
	    else {

	        bc = 0. ;  // dlapse/dr = 0.25*lapse/rr

	    }
	}
	else {
	    if (first) {  // The poisson solver in LORENE assumes the
	                  // asymptotic behavior of the function -> 0
	        bc = 0.5 - 1./sqrt(2.) ;  // <- bc of the real function = 0.5

	    }
	    else {

	        bc = 0. ; // <- bc of the real function = 1./sqrt(2.)

	    }
	}
    }
    else {  // Isotropic coordinates with the maximal slicing

        if (neumann) {

	    if (first) {

	        bc = 0. ;  // dlapse/dr = 0

	    }
	    else {

	        Scalar rr(mp) ;
		rr = mp.r ;
		rr.std_spectral_base() ;

	        int nt = mg->get_nt(0) ;
		int np = mg->get_np(0) ;

		Scalar tmp(mp) ;

		tmp = 0.5 * lapse / rr ;
		// dlapse/dr = lapse/2/rah

		bc = 1. ;
		for (int j=0; j<nt; j++) {
		    for (int k=0; k<np; k++) {
		        bc.set(0,k,j,0) = tmp.val_grid_point(1,k,j,0) ;
		    }
		}


	    }

	}
	else {

	    if (first) {

	        bc = - 0.5 ;  // lapse = 0.5

	    }
	    else {

	        bc = 1./sqrt(2.) - 1. ;  // lapse = 1./sqrt(2.)

	    }

	}

    }

    bc.std_base_scal() ;
    return bc ;

}


const Valeur Black_hole::bc_shift_x(double omega_r) const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    const Mg3d* mg_angu = mg->get_angu() ;
    Valeur bc(mg_angu) ;

    Base_val** bases = mp.get_mg()->std_base_vect_cart() ;

    Scalar rr(mp) ;
    rr = mp.r ;
    rr.std_spectral_base() ;
    Scalar st(mp) ;
    st = mp.sint ;
    st.std_spectral_base() ;
    Scalar cp(mp) ;
    cp = mp.cosp ;
    cp.std_spectral_base() ;
    Scalar yy(mp) ;
    yy = mp.y ;
    yy.std_spectral_base() ;

    Scalar tmp(mp) ;

    if (kerrschild) {

        tmp = lapse_bh * (lapse / confo / confo) * st * cp
	  - omega_r * yy - shift_bh(1) ;

	//        tmp = lap_bh * lap_bh * st * cp - omega_r * yy ;
    }
    else {  // Isotropic coordinates with the maximal slicing

        tmp = (lapse / confo / confo) * st * cp - omega_r * yy ;

    }

    int nt = mg->get_nt(0) ;
    int np = mg->get_np(0) ;

    bc = 1. ;
    for (int j=0; j<nt; j++) {
        for (int k=0; k<np; k++) {
	    bc.set(0,k,j,0) = tmp.val_grid_point(1,k,j,0) ;
	}
    }

    bc.base = *bases[0] ;
    //    bc.std_base_scal() ;

    for (int i=0; i<3; i++)
      delete bases[i] ;

    delete [] bases ;

    return bc ;

}


const Valeur Black_hole::bc_shift_y(double omega_r) const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    const Mg3d* mg_angu = mg->get_angu() ;
    Valeur bc(mg_angu) ;

    Base_val** bases = mp.get_mg()->std_base_vect_cart() ;

    Scalar rr(mp) ;
    rr = mp.r ;
    rr.std_spectral_base() ;
    Scalar st(mp) ;
    st = mp.sint ;
    st.std_spectral_base() ;
    Scalar sp(mp) ;
    sp = mp.sinp ;
    sp.std_spectral_base() ;
    Scalar xx(mp) ;
    xx = mp.x ;
    xx.std_spectral_base() ;

    Scalar tmp(mp) ;

    if (kerrschild) {

        tmp = lapse_bh * (lapse / confo / confo) * st * sp
	  + omega_r * xx - shift_bh(2) ;

	//        tmp = lap_bh * lap_bh * st * sp + omega_r * xx ;
    }
    else {

        tmp = (lapse / confo / confo) * st * sp + omega_r * xx ;

    }

    int nt = mg->get_nt(0) ;
    int np = mg->get_np(0) ;

    bc = 1. ;
    for (int j=0; j<nt; j++) {
        for (int k=0; k<np; k++) {
	    bc.set(0,k,j,0) = tmp.val_grid_point(1,k,j,0) ;
	}
    }

    bc.base = *bases[1] ;
    //    bc.std_base_scal() ;

    for (int i=0; i<3; i++)
      delete bases[i] ;

    delete [] bases ;

    return bc ;

}


const Valeur Black_hole::bc_shift_z() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    const Mg3d* mg_angu = mg->get_angu() ;
    Valeur bc(mg_angu) ;

    Base_val** bases = mp.get_mg()->std_base_vect_cart() ;

    Scalar rr(mp) ;
    rr = mp.r ;
    rr.std_spectral_base() ;
    Scalar ct(mp) ;
    ct = mp.cost ;
    ct.std_spectral_base() ;

    Scalar tmp(mp) ;

    if (kerrschild) {

        tmp = lapse_bh * (lapse / confo / confo) * ct
	  - shift_bh(3) ;

        //        tmp = lap_bh * lap_bh * ct ;
    }
    else {

        tmp = (lapse / confo / confo) * ct ;

    }

    int nt = mg->get_nt(0) ;
    int np = mg->get_np(0) ;

    bc = 1. ;
    for (int j=0; j<nt; j++) {
        for (int k=0; k<np; k++) {
	    bc.set(0,k,j,0) = tmp.val_grid_point(1,k,j,0) ;
	}
    }

    bc.base = *bases[2] ;
    //    bc.std_base_scal() ;

    for (int i=0; i<3; i++)
      delete bases[i] ;

    delete [] bases ;

    return bc ;

}


const Valeur Black_hole::bc_confo() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    const Mg3d* mg_angu = mg->get_angu() ;
    Valeur bc(mg_angu) ;

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

    int nt = mg->get_nt(0) ;
    int np = mg->get_np(0) ;

    Scalar tmp(mp) ;

    if (kerrschild) { // Assumes that r_BH = 1.

        Scalar divshift(mp) ;
	divshift = shift_rs(1).deriv(1) + shift_rs(2).deriv(2)
	  + shift_rs(3).deriv(3) ;
	divshift.std_spectral_base() ;

	Scalar llshift(mp) ;
	llshift = st*cp*shift_rs(1) + st*sp*shift_rs(2) + ct*shift_rs(3) ;
	llshift.std_spectral_base() ;

	Scalar lldllsh = llshift.dsdr() ;
	lldllsh.std_spectral_base() ;

	Scalar tmp1 = divshift ;
	Scalar tmp2 = -3.*lldllsh ;

        Scalar tmp5 = 0.5*confo*(lapse_bh*confo*confo/lapse - 1.)/rr ;
        tmp1.set_dzpuis(tmp5.get_dzpuis()) ;
	tmp2.set_dzpuis(tmp5.get_dzpuis()) ;

        Scalar tmp3 = 2. * lapse_bh * lapse_bh * mass * llshift / rr / rr ;
	Scalar tmp4 = 4. * pow(lapse_bh,3.) * mass * (1.+3.*mass/rr)
	  * lapse_rs / rr / rr ;

	tmp3.set_dzpuis(tmp5.get_dzpuis()) ;
	tmp4.set_dzpuis(tmp5.get_dzpuis()) ;

        tmp = tmp5 + pow(confo,3.)*(tmp1+tmp2+tmp3+tmp4)/12./lapse/lapse_bh ;

	//	tmp = -0.5 * (1. - 2. * mass / rr) / rr ;

    }
    else {  // Isotropic coordinates with the maximal slicing

        Scalar divshift(mp) ;
	divshift = shift(1).deriv(1) + shift(2).deriv(2)
	  + shift(3).deriv(3) ;
	divshift.std_spectral_base() ;

	Scalar llshift(mp) ;
	llshift = st*cp*shift(1) + st*sp*shift(2) + ct*shift(3) ;
	llshift.std_spectral_base() ;

	Scalar lldllsh = llshift.dsdr() ;
	lldllsh.std_spectral_base() ;

	Scalar tmp1 = divshift ;
	Scalar tmp2 = -3.*lldllsh ;

        Scalar tmp5 = - 0.5 * confo / rr ;

	tmp1.set_dzpuis(tmp5.get_dzpuis()) ;
	tmp2.set_dzpuis(tmp5.get_dzpuis()) ;

	tmp = tmp5 + pow(confo, 3.) * (tmp1 + tmp2) / 12. / lapse ;

    }

    bc = 1. ;
    for (int j=0; j<nt; j++) {
        for (int k=0; k<np; k++) {
	    bc.set(0,k,j,0) = tmp.val_grid_point(1,k,j,0) ;
	}
    }

    bc.std_base_scal() ;
    return bc ;

}
