/*
 *  Method of class Bin_bhns to compute the analytic shift vector
 *
 *    (see file bin_bhns.h for documentation).
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

char bin_bhns_shift_ana_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:52:41  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2007/06/22 01:11:08  k_taniguchi
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
#include "bin_bhns.h"
#include "utilitaires.h"
#include "unites.h"

namespace Lorene {
void Bin_bhns::shift_analytic(double reduce_shift_bh, double reduce_shift_ns)
{

    using namespace Unites ;

    double massbh = hole.get_mass_bh() ;
    double massns = star.mass_g_bhns() ;
    double mass_bh = ggrav * massbh ;
    double mass_ns = ggrav * massns ;

    double mass_tot = mass_bh + mass_ns ;

    double comb = mass_bh * mass_ns * omega * separ / mass_tot ;

    //-----------------------------------------//
    //     Shift vector for the black hole     //
    //-----------------------------------------//

    const Vector& shift_bh_old = hole.get_shift_auto_rs() ;

    const Map& mp_bh = hole.get_mp() ;
    Scalar xx_bh(mp_bh) ;
    xx_bh = mp_bh.x ;
    xx_bh.std_spectral_base() ;
    Scalar yy_bh(mp_bh) ;
    yy_bh = mp_bh.y ;
    yy_bh.std_spectral_base() ;
    Scalar zz_bh(mp_bh) ;
    zz_bh = mp_bh.z ;
    zz_bh.std_spectral_base() ;
    Scalar rr_bh(mp_bh) ;
    rr_bh = mp_bh.r ;
    rr_bh.std_spectral_base() ;
    Scalar st_bh(mp_bh) ;
    st_bh = mp_bh.sint ;
    st_bh.std_spectral_base() ;
    Scalar ct_bh(mp_bh) ;
    ct_bh = mp_bh.cost ;
    ct_bh.std_spectral_base() ;
    Scalar sp_bh(mp_bh) ;
    sp_bh = mp_bh.sinp ;
    sp_bh.std_spectral_base() ;
    Scalar cp_bh(mp_bh) ;
    cp_bh = mp_bh.cosp ;
    cp_bh.std_spectral_base() ;

    double rad_bh = rr_bh.val_grid_point(1, 0, 0, 0) ;

    Scalar x_bh_ex(mp_bh) ;
    Scalar y_bh_ex(mp_bh) ;
    Scalar z_bh_ex(mp_bh) ;

    if (hole.is_irrotational()) {

        // x component
        // -----------
	x_bh_ex = 0.2 * comb * rad_bh * rad_bh
	  * st_bh * st_bh * cp_bh * sp_bh / pow(rr_bh, 3.) ;
	x_bh_ex.annule_domain(0) ;
	x_bh_ex.std_spectral_base() ;

	(hole.set_shift_auto_rs()).set(1) = shift_bh_old(1)
	  + reduce_shift_bh * x_bh_ex ;

	// y component
	// -----------
	y_bh_ex = 0.5 * comb * (7. + 0.2*rad_bh*rad_bh/rr_bh/rr_bh) / rr_bh
	  + 0.5 * comb * pow(st_bh*sp_bh,2.)
	  * (1. - 0.6*rad_bh*rad_bh/rr_bh/rr_bh) / rr_bh ;
	y_bh_ex.annule_domain(0) ;
	y_bh_ex.std_spectral_base() ;

	(hole.set_shift_auto_rs()).set(2) = shift_bh_old(2)
	  + reduce_shift_bh * y_bh_ex ;

	// z component
	// -----------
	z_bh_ex = 0.5 * comb * st_bh * sp_bh * ct_bh
	  * (1.-0.6*rad_bh*rad_bh/rr_bh/rr_bh) / rr_bh ;
	z_bh_ex.annule_domain(0) ;
	z_bh_ex.std_spectral_base() ;

	(hole.set_shift_auto_rs()).set(3) = shift_bh_old(3)
	  + reduce_shift_bh * z_bh_ex ;

	(hole.set_shift_auto_rs()).std_spectral_base() ;

    }
    else { // Corotational

        // x component
        // -----------
	x_bh_ex = - 0.6 * mass_ns * omega * rad_bh * rad_bh
	  * st_bh * sp_bh / pow(rr_bh, 2.)
	  + 0.5 * comb * st_bh * st_bh * cp_bh * sp_bh
	  * (1. - 0.6*rad_bh*rad_bh/rr_bh/rr_bh) / rr_bh
	  - 0.6*mass_bh*omega*rad_bh*rad_bh*pow(st_bh,3.)*pow(cp_bh,2.)*sp_bh
	  /pow(rr_bh, 2.) ;
	x_bh_ex.annule_domain(0) ;
	x_bh_ex.std_spectral_base() ;

	(hole.set_shift_auto_rs()).set(1) = shift_bh_old(1)
	  + reduce_shift_bh * x_bh_ex ;

	// y component
	// -----------
	y_bh_ex = 0.5 * comb * (7. + 0.2*rad_bh*rad_bh/rr_bh/rr_bh) / rr_bh
	  - 0.6 * mass_bh * omega * rad_bh * rad_bh * st_bh * cp_bh
	  / pow(rr_bh, 2.)
	  + 0.5 * comb * pow(st_bh*sp_bh,2.)
	  * (1. - 0.6*rad_bh*rad_bh/rr_bh/rr_bh) / rr_bh
	  - 0.6*mass_bh*omega*rad_bh*rad_bh*pow(st_bh,3.)*cp_bh*pow(sp_bh,2.)
	  /pow(rr_bh, 2.) ;
	y_bh_ex.annule_domain(0) ;
	y_bh_ex.std_spectral_base() ;

	(hole.set_shift_auto_rs()).set(2) = shift_bh_old(2)
	  + reduce_shift_bh * y_bh_ex ;

	// z component
	// -----------
	z_bh_ex = 0.5 * comb * st_bh * cp_bh * ct_bh
	  * (1. - 0.6*rad_bh*rad_bh/rr_bh/rr_bh) / rr_bh
	  - 0.6*mass_bh*omega*rad_bh*rad_bh*st_bh*st_bh*cp_bh*sp_bh*ct_bh
	  / pow(rr_bh, 2.) ;
	z_bh_ex.annule_domain(0) ;
	z_bh_ex.std_spectral_base() ;

	(hole.set_shift_auto_rs()).set(3) = shift_bh_old(3)
	  + reduce_shift_bh * z_bh_ex ;

	(hole.set_shift_auto_rs()).std_spectral_base() ;

    }


    //-------------------------------------------//
    //     Shift vector for the neutron star     //
    //-------------------------------------------//
    int nzet = star.get_nzet() ;
    int nz_ns = (star.get_mp()).get_mg()->get_nzone() ;

    const Map& mp_ns = star.get_mp() ;
    Scalar xx_ns(mp_ns) ;
    xx_ns = mp_ns.x ;
    xx_ns.std_spectral_base() ;
    Scalar yy_ns(mp_ns) ;
    yy_ns = mp_ns.y ;
    yy_ns.std_spectral_base() ;
    Scalar zz_ns(mp_ns) ;
    zz_ns = mp_ns.z ;
    zz_ns.std_spectral_base() ;
    Scalar rr_ns(mp_ns) ;
    rr_ns = mp_ns.r ;
    rr_ns.std_spectral_base() ;
    Scalar st_ns(mp_ns) ;
    st_ns = mp_ns.sint ;
    st_ns.std_spectral_base() ;
    Scalar ct_ns(mp_ns) ;
    ct_ns = mp_ns.cost ;
    ct_ns.std_spectral_base() ;
    Scalar sp_ns(mp_ns) ;
    sp_ns = mp_ns.sinp ;
    sp_ns.std_spectral_base() ;
    Scalar cp_ns(mp_ns) ;
    cp_ns = mp_ns.cosp ;
    cp_ns.std_spectral_base() ;

    double rad_ns = rr_ns.val_grid_point(1, 0, 0, 0) ;

    Scalar x_ns_in(mp_ns) ;
    Scalar x_ns_ex(mp_ns) ;
    Scalar y_ns_in(mp_ns) ;
    Scalar y_ns_ex(mp_ns) ;
    Scalar z_ns_in(mp_ns) ;
    Scalar z_ns_ex(mp_ns) ;

    if (star.is_irrotational()) {

        // x component
        // -----------
        x_ns_in = - 0.2 * comb * xx_ns * yy_ns / pow(rad_ns, 3.) ;
	x_ns_in.annule(nzet, nz_ns-1) ;
	x_ns_in.std_spectral_base() ;

	x_ns_ex = - 0.2 * comb * rad_ns * rad_ns
	  * st_ns * st_ns * cp_ns * sp_ns / pow(rr_ns, 3.) ;
	x_ns_ex.annule(0, nzet-1) ;
	x_ns_ex.std_spectral_base() ;

	(star.set_shift_auto()).set(1) = reduce_shift_ns
	  * (x_ns_in + x_ns_ex) ;

	// y component
	// -----------
	y_ns_in = - 0.5 * comb * (11. - 3.8*rr_ns*rr_ns/rad_ns/rad_ns) / rad_ns
	  - 0.2 * comb * yy_ns * yy_ns / pow(rad_ns, 3.) ;
	y_ns_in.annule(nzet, nz_ns-1) ;
	y_ns_in.std_spectral_base() ;

	y_ns_ex = - 0.5 * comb * (7. + 0.2*rad_ns*rad_ns/rr_ns/rr_ns) / rr_ns
	  - 0.5 * comb * pow(st_ns*sp_ns,2.)
	  * (1. - 0.6*rad_ns*rad_ns/rr_ns/rr_ns) / rr_ns ;
	y_ns_ex.annule(0, nzet-1) ;
	y_ns_ex.std_spectral_base() ;

	(star.set_shift_auto()).set(2) = reduce_shift_ns
	  * (y_ns_in + y_ns_ex) ;

	// z component
	// -----------
	z_ns_in = - 0.2 * comb * yy_ns * zz_ns / pow(rad_ns, 3.) ;
	z_ns_in.annule(nzet, nz_ns-1) ;
	z_ns_in.std_spectral_base() ;

	z_ns_ex = - 0.5 * comb * st_ns * sp_ns * ct_ns
	  * (1.-0.6*rad_ns*rad_ns/rr_ns/rr_ns) / rr_ns ;
	z_ns_ex.annule(0, nzet-1) ;
	z_ns_ex.std_spectral_base() ;

	(star.set_shift_auto()).set(3) = reduce_shift_ns
	  * (z_ns_in + z_ns_ex) ;

    }
    else { // Corotational

        // x component
        // -----------
        x_ns_in = 1.5 * mass_ns * omega * yy_ns
	  * (1. - 0.6*rr_ns*rr_ns/rad_ns/rad_ns) / rad_ns
	  - 0.2 * comb * xx_ns * yy_ns / pow(rad_ns, 3.)
	  + 0.6 * mass_ns * omega * xx_ns * xx_ns * yy_ns / pow(rad_ns, 3.) ;
        x_ns_in.annule(nzet, nz_ns-1) ;
	x_ns_in.std_spectral_base() ;

	x_ns_ex = 0.6 * mass_ns * omega * rad_ns * rad_ns
	  * st_ns * sp_ns / pow(rr_ns, 2.)
	  - 0.5 * comb * st_ns * st_ns * cp_ns * sp_ns
	  * (1. - 0.6*rad_ns*rad_ns/rr_ns/rr_ns) / rr_ns
	  + 0.6*mass_ns*omega*rad_ns*rad_ns*pow(st_ns,3.)*pow(cp_ns,2.)*sp_ns
	  /pow(rr_ns, 2.) ;
	x_ns_ex.annule(0, nzet-1) ;
	x_ns_ex.std_spectral_base() ;

	(star.set_shift_auto()).set(1) = reduce_shift_ns
	  * (x_ns_in + x_ns_ex) ;

	// y component
	// -----------
	y_ns_in = - 0.5 * comb * (11. - 3.8*rr_ns*rr_ns/rad_ns/rad_ns) / rad_ns
	  + 1.5 * mass_ns * omega * xx_ns
	  * (1. - 0.6*rr_ns*rr_ns/rad_ns/rad_ns) / rad_ns
	  - 0.2 * comb * yy_ns * yy_ns / pow(rad_ns, 3.)
	  + 0.6 * mass_ns * omega * xx_ns * yy_ns * yy_ns / pow(rad_ns, 3.) ;
	y_ns_in.annule(nzet, nz_ns-1) ;
	y_ns_in.std_spectral_base() ;

	y_ns_ex = - 0.5 * comb * (7. + 0.2*rad_ns*rad_ns/rr_ns/rr_ns) / rr_ns
	  + 0.6 * mass_ns * omega * rad_ns * rad_ns * st_ns * cp_ns
	  / pow(rr_ns, 2.)
	  - 0.5 * comb * pow(st_ns*sp_ns,2.)
	  * (1. - 0.6*rad_ns*rad_ns/rr_ns/rr_ns) / rr_ns
	  + 0.6*mass_ns*omega*rad_ns*rad_ns*pow(st_ns,3.)*cp_ns*pow(sp_ns,2.)
	  /pow(rr_ns, 2.) ;
	y_ns_ex.annule(0, nzet-1) ;
	y_ns_ex.std_spectral_base() ;

	(star.set_shift_auto()).set(2) = reduce_shift_ns
	  * (y_ns_in + y_ns_ex) ;

	// z component
	// -----------
	z_ns_in = - 0.2 * comb * yy_ns * zz_ns / pow(rad_ns, 3.)
	  + 0.6 * mass_ns * omega * xx_ns * yy_ns * zz_ns / pow(rad_ns, 3.) ;
	z_ns_in.annule(nzet, nz_ns-1) ;
	z_ns_in.std_spectral_base() ;

	z_ns_ex = - 0.5 * comb * st_ns * cp_ns * ct_ns
	  * (1. - 0.6*rad_ns*rad_ns/rr_ns/rr_ns) / rr_ns
	  + 0.6*mass_ns*omega*rad_ns*rad_ns*st_ns*st_ns*cp_ns*sp_ns*ct_ns
	  / pow(rr_ns, 2.) ;
	z_ns_ex.annule(0, nzet-1) ;
	z_ns_ex.std_spectral_base() ;

	(star.set_shift_auto()).set(3) = reduce_shift_ns
	  * (z_ns_in + z_ns_ex) ;

    }

    (star.set_shift_auto()).std_spectral_base() ;

}
}
