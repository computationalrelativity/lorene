/*
 *  Method of class Bin_ns_bh to set some analytical form
 *   to the shift vector of the neutron star
 *
 *    (see file bin_ns_bh.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Keisuke Taniguchi
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

char bin_ns_bh_anashift_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:52:43  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:01  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2005/11/30 11:09:06  p_grandclement
 * Changes for the Bin_ns_bh project
 *
 * Revision 1.1  2004/06/09 06:28:21  k_taniguchi
 * First revision.
 *
 *
 * $Header$
 *
 */

// C headers
#include <cmath>

// Lorene headers
#include "bin_ns_bh.h"
#include "unites.h"

namespace Lorene {
void Bin_ns_bh::analytical_shift() {

    // NS-BH binary systems should be relativistic
    // -------------------------------------------
    if ( !star.is_relativistic() ) {
        abort() ;
    }

    using namespace Unites ;

    // Radius of the neutron star
    double a0 = star.ray_eq() ;

    // Mass ratio
    double p_mass = star.mass_g() / (hole.masse_adm_seul()/ggrav) ;

    // G M Omega R / (1+p)
    double www = ggrav * star.mass_g() * omega
      * separation() / (1. + p_mass) ;

    const Map& mp = star.get_mp() ;
    Cmp tmp(mp) ;
    Cmp tmp_ext(mp) ;
    int nzet = star.get_nzet() ;
    int nzm1 = mp.get_mg()->get_nzone() - 1 ;

    //-------------------
    // Irrotational case
    //-------------------
    // Since this formula is only an initial guess, we use it
    //  also for the corotating case.

    // Computation of w_shift
    // ----------------------
    star.set_w_shift().set_etat_qcq() ;

    // X component
    // -----------
    star.set_w_shift().set(0) = 0. ;

    // Y component
    // -----------

    // For the incompressible case :
    tmp = -6. * www / a0 * ( 1. - (mp.r)*(mp.r) / (3.*a0*a0) ) ;

    tmp.annule(nzet, nzm1) ;
    tmp_ext = -4. * www / mp.r ;
    tmp_ext.annule(0, nzet-1) ;

    star.set_w_shift().set(1) = - tmp - tmp_ext ;

    // Z component
    // -----------
    star.set_w_shift().set(2) = 0. ;

    // Sets the standard spectral bases for Cartesian components
    star.set_w_shift().set_std_base() ;

    // Computation of khi_shift
    //-------------------------

    tmp = 2. * www / a0 * (mp.y) * ( 1. - 3.*(mp.r)*(mp.r) / (5.*a0*a0) ) ;
    tmp.annule(nzet, nzm1) ;
    tmp_ext = 0.8 * www * a0 * a0 * (mp.sint) * (mp.sinp) / ((mp.r)*(mp.r)) ;
    tmp_ext.annule(0, nzet-1) ;

    star.set_khi_shift() = - tmp - tmp_ext ;

    // Sets the standard spectral bases for a scalar field
    star.set_khi_shift().set_std_base() ;

}
}
