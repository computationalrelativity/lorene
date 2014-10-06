/*
 *  Method of class Bin_bhns_extr to set some analytical form
 *
 *    (see file bin_bhns_extr.h for documentation).
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

char bin_bhns_extr_anashift_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:13:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/11/30 20:45:29  k_taniguchi
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
#include "unites.h"

void Bin_bhns_extr::analytical_shift() {

  using namespace Unites ;

    // BH-NS binary systems should be relativistic
    // -------------------------------------------

    if ( !star.is_relativistic() ) {

        cout << "BH-NS binary systems should be relativistic !!!" << endl ;
        abort() ;
    }

    // Radius of the neutron star
    double a0 = star.ray_eq() ;

    // G M Omega R
    double www = ggrav * star.mass_g() * omega * separ ;
    // Approximates the mass ratio -> 0

    const Map& mp = star.get_mp() ;
    Tenseur tmp(mp) ;
    Tenseur tmp_ext(mp) ;
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
    tmp.set_etat_qcq() ;
    tmp.set() = 6. * www / a0 * ( 1. - (mp.r)*(mp.r) / (3.*a0*a0) ) ;
    tmp.set().annule(nzet, nzm1) ;
    tmp.set_std_base() ;

    tmp_ext.set_etat_qcq() ;
    tmp_ext.set() = 4. * www / mp.r ;
    tmp_ext.set().annule(0, nzet-1) ;
    tmp_ext.set_std_base() ;

    star.set_w_shift().set(1) = tmp() + tmp_ext() ;

    // Z component
    // -----------
    star.set_w_shift().set(2) = 0. ;

    // Sets the standard spectral bases for Cartesian components
    star.set_w_shift().set_std_base() ;

    // Computation of khi_shift
    //-------------------------

    tmp.set() = 2. * www / a0 * (mp.y)
      * ( 1. - 3.*(mp.r)*(mp.r) / (5.*a0*a0) ) ;
    tmp.set().annule(nzet, nzm1) ;
    tmp.set_std_base() ;

    tmp_ext.set() = 0.8 * www * a0 * a0 * (mp.sint) * (mp.sinp)
      / ((mp.r)*(mp.r)) ;
    tmp_ext.set().annule(0, nzet-1) ;
    tmp_ext.set_std_base() ;

    star.set_khi_shift() = tmp + tmp_ext ;

    // Sets the standard spectral bases for a scalar field
    star.set_khi_shift().set_std_base() ;

}
