/*
 *  Methods for computing global quantities within the class Star_rot_Dirac
 *
 *    (see file star_rot_dirac.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Lap-Ming Lin & Jerome Novak
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

char strot_dirac_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2005/01/31 08:51:48  j_novak
 * New files for rotating stars in Dirac gauge (still under developement).
 *
 *
 *
 * $Header$
 *
 */


// C headers
#include <math.h>
#include <assert.h>

// Lorene headers
#include "star_rot_dirac.h"
#include "utilitaires.h" 


        //-------------------------------------------------//
        //                Baryonic mass                    //
        //                                                 //
        //  Note: In Lorene units, neutron mass is unity   //
        //-------------------------------------------------//


double Star_rot_Dirac::mass_b() const {

  if (p_mass_b == 0x0) {    // a new computation is required

    Scalar dens = sqrt( gamma.determinant() ) * gam_euler * nbar ;

    dens.std_spectral_base() ;

    p_mass_b = new double( dens.integrale() ) ;

  }

  return *p_mass_b ;

}

          //---------------------------------------------------------//
          //            Gravitational mass                           //   
          //                                                         //
          // Note: This is the Komar mass for stationary and         //
          //       asymptotically flat spacetime (see, eg, Wald)     //
          //---------------------------------------------------------//

double Star_rot_Dirac::mass_g() const {

  if (p_mass_g == 0x0) {    // a new computation is required

    Scalar j_source = 2.* contract(contract(gamma.cov(), 0, j_euler, 0), 
				   0, shift, 0) ;
           
    Scalar dens = sqrt( gamma.determinant() ) *  
	      ( nnn * (ener_euler + s_euler) - j_source ) ;

    dens.std_spectral_base() ;

    p_mass_g = new double( dens.integrale() ) ;

  }

  return *p_mass_g ;

}

                //--------------------------------------//
                //    Angular momentum                  //
                //                                      // 
                // Komar-type integral (see, eg, Wald)  //
                //--------------------------------------// 

double Star_rot_Dirac::angu_mom() const {

  if (p_angu_mom == 0x0) {    // a new computation is required

    Vector phi_kill(mp, CON, mp.get_bvect_spher()) ;

    phi_kill.set(1).set_etat_zero() ;
    phi_kill.set(2).set_etat_zero() ;
    phi_kill.set(3) = 1. ;
    phi_kill.set(3).std_spectral_base() ;
    phi_kill.set(3).mult_rsint() ;

    Scalar j_source = contract(contract(gamma.cov(), 0, j_euler, 0),
			       0, phi_kill, 0) ;

    Scalar dens = sqrt( gamma.determinant() ) * j_source  ;

    dens.std_spectral_base() ;

    p_angu_mom = new double( dens.integrale() ) ;


  }

  return *p_angu_mom ;

}






