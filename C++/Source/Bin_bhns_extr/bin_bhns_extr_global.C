/*
 *  Methods of class Bin_bhns_extr to compute global quantities
 *
 *    (see file bin_bhns_extr.h for documentation).
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

char bin_bhns_extr_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:52:42  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2005/02/28 23:07:12  k_taniguchi
 * Suppression of the ADM mass and so on.
 *
 * Revision 1.1  2004/11/30 20:46:13  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include <cmath>

// Lorene headers
#include "bin_bhns_extr.h"
#include "coord.h"
#include "unites.h"

          //--------------------------------------------------//
          //          X coordinate of the barycenter          //
          //--------------------------------------------------//

namespace Lorene {
double Bin_bhns_extr::xa_barycenter_extr() const {

  using namespace Unites ;

    if (p_xa_barycenter_extr == 0x0) {    // a new computation is required

        p_xa_barycenter_extr = new double ;

	*p_xa_barycenter_extr = 0 ;

        const Map& mp = star.get_mp() ;
        Cmp xxa(mp) ;
	xxa = mp.xa ;	// Absolute X coordinate
	xxa.std_base_scal() ;

	if (star.in_kerrschild()) { // Kerr-Schild background metric

	    const Coord& xx = mp.x ;
	    const Coord& yy = mp.y ;
	    const Coord& zz = mp.z ;

	    Tenseur r_bh(mp) ;
	    r_bh.set_etat_qcq() ;
	    r_bh.set() = pow( (xx+separ)*(xx+separ) + yy*yy + zz*zz, 0.5) ;
	    r_bh.set_std_base() ;

	    Tenseur msr(mp) ;
	    msr = ggrav * mass_bh / r_bh ;
	    msr.set_std_base() ;

	    Cmp tmp = sqrt(1. + 2.*msr()) ;
	    tmp.std_base_scal() ;

	    Tenseur acar = star.get_a_car() ;
	    acar.set_std_base() ;

	    Tenseur g_euler = star.get_gam_euler() ;
	    g_euler.set_std_base() ;

	    Tenseur nbary = star.get_nbar() ;
	    nbary.set_std_base() ;

	    Cmp dens = acar() * sqrt(acar()) * g_euler() * nbary()
	      * tmp * xxa ;
	    dens.std_base_scal() ;

	    *p_xa_barycenter_extr = dens.integrale() / mass_b_extr() ;

	}
	else { // Conformally flat background metrci

	    Tenseur acar = star.get_a_car() ;
	    acar.set_std_base() ;

	    Tenseur g_euler = star.get_gam_euler() ;
	    g_euler.set_std_base() ;

	    Tenseur nbary = star.get_nbar() ;
	    nbary.set_std_base() ;

	    Cmp dens = acar() * sqrt(acar()) * g_euler() * nbary() * xxa ;
	    dens.std_base_scal() ;

	    *p_xa_barycenter_extr = dens.integrale() / mass_b_extr() ;

	}

    }

    return *p_xa_barycenter_extr ;

}

          //--------------------------------------------------//
          //          Y coordinate of the barycenter          //
          //--------------------------------------------------//

double Bin_bhns_extr::ya_barycenter_extr() const {

  using namespace Unites ;

    if (p_ya_barycenter_extr == 0x0) {    // a new computation is required

        p_ya_barycenter_extr = new double ;

	*p_ya_barycenter_extr = 0 ;

        const Map& mp = star.get_mp() ;
        Cmp yya(mp) ;
	yya = mp.ya ;	// Absolute Y coordinate
	yya.std_base_scal() ;

	if (star.in_kerrschild()) { // Kerr-Schild background metric

	    const Coord& xx = mp.x ;
	    const Coord& yy = mp.y ;
	    const Coord& zz = mp.z ;

	    Tenseur r_bh(mp) ;
	    r_bh.set_etat_qcq() ;
	    r_bh.set() = pow( (xx+separ)*(xx+separ) + yy*yy + zz*zz, 0.5) ;
	    r_bh.set_std_base() ;

	    Tenseur msr(mp) ;
	    msr = ggrav * mass_bh / r_bh ;
	    msr.set_std_base() ;

	    Cmp tmp = sqrt(1. + 2.*msr()) ;
	    tmp.std_base_scal() ;

	    Tenseur acar = star.get_a_car() ;
	    acar.set_std_base() ;

	    Tenseur g_euler = star.get_gam_euler() ;
	    g_euler.set_std_base() ;

	    Tenseur nbary = star.get_nbar() ;
	    nbary.set_std_base() ;

	    Cmp dens = acar() * sqrt(acar()) * g_euler() * nbary()
	      * tmp * yya ;
	    dens.std_base_scal() ;

	    *p_ya_barycenter_extr = dens.integrale() / mass_b_extr() ;

	}
	else { // Conformally flat background metric
	       // It should be zero !

	    Tenseur acar = star.get_a_car() ;
	    acar.set_std_base() ;

	    Tenseur g_euler = star.get_gam_euler() ;
	    g_euler.set_std_base() ;

	    Tenseur nbary = star.get_nbar() ;
	    nbary.set_std_base() ;

	    Cmp dens = acar() * sqrt(acar()) * g_euler() * nbary() * yya ;
	    dens.std_base_scal() ;

	    *p_ya_barycenter_extr = dens.integrale() / mass_b_extr() ;

	}

    }

    return *p_ya_barycenter_extr ;

}

          //-------------------------------//
          //          Baryon mass          //
          //-------------------------------//

double Bin_bhns_extr::mass_b_extr() const {

  using namespace Unites ;

    if (p_mass_b_extr == 0x0) {    // a new computation is required

        p_mass_b_extr = new double ;

        if (star.is_relativistic()) {  // Relativistic case

	    *p_mass_b_extr = 0 ;

	    if (star.in_kerrschild()) { // Kerr-Schild background metric

	        const Map& mp = star.get_mp() ;

		const Coord& xx = mp.x ;
		const Coord& yy = mp.y ;
		const Coord& zz = mp.z ;

		Tenseur r_bh(mp) ;
		r_bh.set_etat_qcq() ;
		r_bh.set() = pow( (xx+separ)*(xx+separ) + yy*yy + zz*zz, 0.5) ;
		r_bh.set_std_base() ;

		Tenseur msr(mp) ;
		msr = ggrav * mass_bh / r_bh ;
		msr.set_std_base() ;

		Cmp tmp = sqrt(1. + 2.*msr()) ;
		tmp.std_base_scal() ;

		Tenseur acar = star.get_a_car() ;
		acar.set_std_base() ;

		Tenseur g_euler = star.get_gam_euler() ;
		g_euler.set_std_base() ;

		Tenseur nbary = star.get_nbar() ;
		nbary.set_std_base() ;

		Cmp dens = acar() * sqrt(acar()) * g_euler() * nbary() * tmp ;
		dens.std_base_scal() ;

		*p_mass_b_extr = dens.integrale() ;

	    }
	    else { // Conformally flat background metric

	        Tenseur acar = star.get_a_car() ;
		acar.set_std_base() ;

		Tenseur g_euler = star.get_gam_euler() ;
		g_euler.set_std_base() ;

		Tenseur nbary = star.get_nbar() ;
		nbary.set_std_base() ;

		Cmp dens = acar() * sqrt(acar()) * g_euler() * nbary() ;
		dens.std_base_scal() ;

		*p_mass_b_extr = dens.integrale() ;

	    }

	}
	else {

	    cout << "BH-NS binary system should be relativistic!!!" << endl ;
	    abort() ;

	}
    }

    return *p_mass_b_extr ;

}
}
