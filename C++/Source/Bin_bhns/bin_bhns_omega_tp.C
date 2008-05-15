/*
 *  Methods of class Bin_bhns to compute an orbital angular velocity
 *   from two points at the stellar surface
 *
 *    (see file bin_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2006-2007 Keisuke Taniguchi
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

char star_bhns_omega_tp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2008/05/15 19:00:27  k_taniguchi
 * Change of some parameters.
 *
 * Revision 1.1  2007/06/22 01:10:00  k_taniguchi
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
#include "bin_bhns.h"
#include "unites.h"
#include "utilitaires.h"

          //---------------------------------------------------//
          //     Orbaital angular velocity from two points     //
          //---------------------------------------------------//

double Bin_bhns::omega_two_points() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_omega_two_points == 0x0) {   // a new computation is required

        double omega_two ;

	const Scalar& lapconf = star.get_lapconf_tot() ;
	const Scalar& confo = star.get_confo_tot() ;
	const Scalar& psi4 = star.get_psi4() ;
	const Vector& shift = star.get_shift_tot() ;
	const Scalar& gam = star.get_gam() ;

	int ii = (star.get_mp()).get_mg()->get_nr(0) - 1 ;
	int jj = (star.get_mp()).get_mg()->get_nt(0) - 1 ;
	int ka = 0 ;
	int kb = (star.get_mp()).get_mg()->get_np(0) / 2 ;

	double psi4_a = psi4.val_grid_point(0,ka,jj,ii) ;
	double psi4_b = psi4.val_grid_point(0,kb,jj,ii) ;
	double con2_a = confo.val_grid_point(0,ka,jj,ii)
	  * confo.val_grid_point(0,ka,jj,ii) ;
	double con2_b = confo.val_grid_point(0,kb,jj,ii)
	  * confo.val_grid_point(0,kb,jj,ii) ;
	double gam2_a = gam.val_grid_point(0,ka,jj,ii)
	  * gam.val_grid_point(0,ka,jj,ii) ;
	double gam2_b = gam.val_grid_point(0,kb,jj,ii)
	  * gam.val_grid_point(0,kb,jj,ii) ;
	double lap2_a = lapconf.val_grid_point(0,ka,jj,ii)
	  * lapconf.val_grid_point(0,ka,jj,ii) ;
	double lap2_b = lapconf.val_grid_point(0,kb,jj,ii)
	  * lapconf.val_grid_point(0,kb,jj,ii) ;
	double shiftx_a = shift(1).val_grid_point(0,ka,jj,ii) ;
	double shiftx_b = shift(1).val_grid_point(0,kb,jj,ii) ;
	double shifty_a = shift(2).val_grid_point(0,ka,jj,ii) ;
	double shifty_b = shift(2).val_grid_point(0,kb,jj,ii) ;
	double shiftz_a = shift(3).val_grid_point(0,ka,jj,ii) ;
	double shiftz_b = shift(3).val_grid_point(0,kb,jj,ii) ;

	double xns_rot = (star.get_mp()).get_ori_x() - x_rot ;
	double yns_rot = (star.get_mp()).get_ori_y() - y_rot ;

	double ra = star.ray_eq() ;
	double rb = star.ray_eq_pi() ;

	if (hole.is_kerrschild()) {

	    cout << "!!!!! WARNING: Not yet available !!!!!" << endl ;
	    abort() ;
	    /*
	    double y_separ = (star.get_mp()).get_ori_y() ;
	    double xbh_rot = (hole.get_mp()).get_ori_x() - x_rot ;
	    double mass = ggrav * hole.get_mass_bh() ;
	    double rbh_a = sqrt( (ra+separ)*(ra+separ) + y_separ*y_separ ) ;
	    double rbh_b = sqrt( (-rb+separ)*(-rb+separ) + y_separ*y_separ ) ;

	    double msr_a = 2.*mass / pow(rbh_a, 3.) ;
	    double msr_b = 2.*mass / pow(rbh_b, 3.) ;

	    double sa = shiftx_a*shiftx_a+shifty_a*shifty_a+shiftz_a*shiftz_a
	      + msr_a * ((ra+separ)*shiftx_a + y_separ*shifty_a)
	      * ((ra+separ)*shiftx_a + y_separ*shifty_a) ;

	    double sb = shiftx_b*shiftx_b+shifty_b*shifty_b+shiftz_b*shiftz_b
	      + msr_b * ((-rb+separ)*shiftx_b + y_separ*shifty_b)
	      * ((-rb+separ)*shiftx_b + y_separ*shifty_b) ;

	    double ta = -shiftx_a*yns_rot + shifty_a*(ra+xns_rot)
	      + msr_a * ((ra+separ)*shiftx_a + y_separ*shifty_a)
	      * y_separ * xbh_rot ;

	    double tb = -shiftx_b*yns_rot + shifty_b*(-rb+xns_rot)
	      + msr_b * ((-rb+separ)*shiftx_b + y_separ*shifty_b)
	      * y_separ * xbh_rot ;

	    double ua = yns_rot*yns_rot + (ra+xns_rot)*(ra+xns_rot)
	      + msr_a * y_separ * y_separ * xbh_rot * xbh_rot ;

	    double ub = yns_rot*yns_rot + (-rb+xns_rot)*(-rb+xns_rot)
	      + msr_b * y_separ * y_separ * xbh_rot * xbh_rot ;

	    // Coefficients : Omega^2 * aaa + 2*Omega * bbb + ccc = 0
	    // ------------------------------------------------------

	    double aaa = psi4_a * gam2_a * ua - psi4_b * gam2_b * ub ;
	    double bbb = psi4_a * gam2_a * ta - psi4_b * gam2_b * tb ;
	    double ccc = psi4_a * gam2_a * sa - psi4_b * gam2_b * sb
	      - lap2_a * gam2_a + lap2_b * gam2_b ;

	    // Term inside the square root : ddd = bbb*bbb - aaa*ccc
	    // -----------------------------------------------------

	    double ddd = bbb*bbb - aaa*ccc ;

	    if ( ddd < 0 ) {
	        cout <<
		  "!!! WARNING : Omega (from two points) does not exist !!!"
		     << endl ;

		omega_two = 0. ;
	    }
	    else {

	        double omega_1 = (-bbb + sqrt(ddd)) / aaa ;
		double omega_2 = (-bbb - sqrt(ddd)) / aaa ;

		cout << "Bin_bhns::omega_two_points:" << endl ;
		cout << "   omega_1 : " << omega_1 * f_unit << " [rad/s]"
		     << endl ;
		cout << "   omega_2 : " << omega_2 * f_unit << " [rad/s]"
		     << endl ;

		omega_two = omega_1 ;

	    }
	    */
	}
	else {  // Isotropic coordinates with the maximal slicing

	    double sa = shiftx_a*shiftx_a+shifty_a*shifty_a+shiftz_a*shiftz_a ;
	    double sb = shiftx_b*shiftx_b+shifty_b*shifty_b+shiftz_b*shiftz_b ;

	    double ta = -shiftx_a*yns_rot + shifty_a*(ra+xns_rot) ;
	    double tb = -shiftx_b*yns_rot + shifty_b*(-rb+xns_rot) ;

	    double ua = yns_rot*yns_rot + (ra+xns_rot)*(ra+xns_rot) ;
	    double ub = yns_rot*yns_rot + (-rb+xns_rot)*(-rb+xns_rot) ;

	    // Coefficients : Omega^2 * aaa + 2*Omega * bbb + ccc = 0
	    // ------------------------------------------------------

	    double aaa = psi4_a * gam2_a * ua - psi4_b * gam2_b * ub ;
	    double bbb = psi4_a * gam2_a * ta - psi4_b * gam2_b * tb ;
	    double ccc = psi4_a * gam2_a * sa - psi4_b * gam2_b * sb
	      - lap2_a * gam2_a / con2_a + lap2_b * gam2_b / con2_b ;

	    // Term inside the square root : ddd = bbb*bbb - aaa*ccc
	    // -----------------------------------------------------

	    double ddd = bbb*bbb - aaa*ccc ;

	    if ( ddd < 0 ) {
	        cout <<
		  "!!! WARNING : Omega (from two points) does not exist !!!"
		     << endl ;

		omega_two = 0. ;
	    }
	    else {

	        double omega_1 = (-bbb + sqrt(ddd)) / aaa ;
		double omega_2 = (-bbb - sqrt(ddd)) / aaa ;

		cout << "Bin_bhns::omega_two_points:" << endl ;
		cout << "   omega_1 : " << omega_1 * f_unit << " [rad/s]"
		     << endl ;
		cout << "   omega_2 : " << omega_2 * f_unit << " [rad/s]"
		     << endl ;

		omega_two = omega_1 ;

	    }

	}

        p_omega_two_points = new double( omega_two ) ;

    }

    return *p_omega_two_points ;

}
