/*
 *  Methods for computing global quantities within the class Star_rot_Dirac
 *
 *    (see file star.h for documentation).
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
 * Revision 1.3  2005/02/02 10:11:24  j_novak
 * Better calculation of the GRV3 identity.
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
#include "unites.h"
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

    // phi_kill = axial killing vector 

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


                  //-----------------------//
                  //        GRV2           //   
                  // ** still under development  //
                  //-----------------------//

double Star_rot_Dirac::grv2() const {

  using namespace Unites ;
  if (p_grv2 == 0x0) {    // a new computation is required


    Scalar u_square = contract(contract(gamma.cov(),0, u_euler, 0),
			     0, u_euler, 0) ;

    //**
    // sou_m = 8\pi T_{\mu\nu} m^{\mu}m^{\nu}
    // m^{\mu} = (0,0,0, r sint/M) in the spherical orthonormal basis. 
    //
    // => sou_m = (r sint)^2 [ (E+P) U^2 + P ], U=v_i v^i 
    //
    // GRV2 paper (cf. Bonazzola & Gourgoulhon CQG 11, 1775 (1994). 
    //**

    Scalar sou_m = 2 * qpig * ( (ener_euler + press)*u_square + press ) ;

    sou_m.std_spectral_base() ;

    sou_m.mult_rsint() ;
    
    sou_m.mult_rsint() ;


    // aa_quad = \tilde{A}_{ij} A^{ij} = K_{ij} K^{ij} (for trace(K)=0)

    Scalar sou_q = 1.5 * aa_quad ;

    // Here is the term \nu_{|| a}\nu^{|| a} in the GRV2 paper. 
    //

    Scalar sou_tmp = gamma.con()(1,1) * logn.dsdr() * logn.dsdr() ;
    
    Scalar term_2 = 2 * gamma.con()(1,2) * logn.dsdr() * logn.dsdt() ;
    term_2.div_r_dzpuis(4) ;

    Scalar term_3 = gamma.con()(2,2) * logn.dsdt() * logn.dsdt() ;
    term_3.div_r_dzpuis(4) ;
    term_3.div_r_dzpuis(4) ;

    sou_tmp += term_2 + term_3 ;

    sou_q -= sou_tmp ;


    p_grv2 = new double( double(1) - lambda_grv2(sou_m, sou_q) ) ;
    
  }

  return *p_grv2 ;

}


     //-------------------------------------------------------------//
     //                 GRV3                                        //
     // cf. Eq. (29) of Gourgoulhon & Bonazzola CQG, 11, 443 (1994) //
     //-------------------------------------------------------------//

double Star_rot_Dirac::grv3() const {

  using namespace Unites ;

  if (p_grv3 == 0x0) {    // a new computation is required

    // Gravitational term 
    // -------------------

    Scalar sou_q = 0.75*aa_quad - contract(logn.derive_cov(gamma), 0,
					   logn.derive_con(gamma), 0) ;


    Tensor t_tmp = contract(gamma.connect().get_delta(), 2, 
			    gamma.connect().get_delta(), 0) ;

    Scalar tmp_1 = 0.25* contract( gamma.con(), 0, 1, 
	    contract(t_tmp, 0, 3), 0, 1 ) ;

    Scalar tmp_2 = 0.25* contract( gamma.con(), 0, 1, 
	      contract( contract( gamma.connect().get_delta(), 0, 1), 
                       0, gamma.connect().get_delta(), 0), 0, 1)  ;

    sou_q = sou_q + tmp_1 - tmp_2 ;

    sou_q = sqrt( gamma.determinant() ) * sou_q ; 

    sou_q.std_spectral_base() ;

    double int_grav = sou_q.integrale() ;


    // Matter term 
    // --------------

    Scalar sou_m = qpig*s_euler ;

    sou_m = sqrt( gamma.determinant() ) * sou_m ;

    sou_m.std_spectral_base() ;

    double int_mat = sou_m.integrale() ;

    //    p_grv3 = new double( (int_grav + int_mat) / int_mat ) ;
    
    p_grv3 = new double( int_grav + int_mat ) ;


  }

  return *p_grv3 ;

}
