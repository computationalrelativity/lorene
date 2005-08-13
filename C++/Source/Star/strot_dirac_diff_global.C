/*
 *  Methods for computing global quantities within the class Star_rot_Dirac_diff
 *
 *    (see file star.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Motoyuki Saijo
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

char strot_dirac_diff_global_C[] = "$Header$" ;

/*
 * $Header$
 *
 */


// C headers
#include <math.h>

// Lorene headers
#include "star_rot_dirac_diff.h"


                     //---------------------//
                     //        T/W          //
                     //---------------------//

double Star_rot_Dirac_diff::tsw() const {

  if (p_tsw == 0x0) {    // a new computation is required

    Vector phi_kill(mp, CON, mp.get_bvect_spher()) ;

    phi_kill.set(1).set_etat_zero() ;
    phi_kill.set(2).set_etat_zero() ;
    phi_kill.set(3) = 1. ;
    phi_kill.set(3).std_spectral_base() ;
    phi_kill.set(3).mult_rsint() ;

    Scalar j_source = contract(contract(gamma.cov(), 0, j_euler, 0),
                               0, phi_kill, 0) ;

    Scalar dens = sqrt( gamma.determinant() ) * j_source * omega_field ;

    dens.std_spectral_base() ;

    double tcin = 0.5 * dens.integrale() ;

    Scalar dens2 = sqrt( gamma.determinant() ) * gam_euler * ener ;
    
    dens2.std_spectral_base() ;
    
    double mass_p = dens2.integrale() ;

    p_tsw = new double( tcin / ( mass_p + tcin - mass_g() ) ) ;

 }

  return *p_tsw ;

    cout << "T/W :             " << p_tsw << '\n' ;

}
