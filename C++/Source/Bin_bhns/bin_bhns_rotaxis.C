/*
 *  Methods of class Bin_bhns to compute the location of the rotation axis
 *
 *    (see file bin_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005-2006 Keisuke Taniguchi
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

char bin_bhns_rotaxis_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:13:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2007/06/22 01:10:47  k_taniguchi
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
#include "unites.h"

void Bin_bhns::rotation_axis_x(double rot_exp_x) {

    using namespace Unites ;

    double momunit = (hole.get_mass_bh()+star.mass_g_bhns())*omega*separ ;

    double error_y = line_mom_bhns()(1) / momunit ;

    if (error_y >= 1.) {
      cout << "Bin_bhns::rotation_axis:" << endl ;
      cout << "   !!! WARNING : error_y is larger than +1 !!!" << endl ;
      error_y *= 0.1 ;
    }

    // Sets X_BH and X_NS
    // ------------------

    double gg = pow( (2.-error_y)/(2.-2.*error_y), rot_exp_x) ;

    double xbh_old = hole.get_mp().get_ori_x() ;

    cout << "Bin_bhns::rotation_axis:" << endl ;
    cout << "   error_y : " << error_y << "   gg : " << gg << endl ;
    cout << "   old X_BH : " << hole.get_mp().get_ori_x() / km << " [km]"
	 << "   old X_NS : " << star.get_mp().get_ori_x() / km << " [km]"
	 << endl ;

    double xbh_new = xbh_old * gg ;
    double xns_new = xbh_new + separ ;

    cout << "   new X_BH : " << xbh_new / km << " [km]"
	 << "   new X_NS : " << xns_new / km << " [km]"
	 << endl ;

    double yns_old = star.get_mp().get_ori_y() ;

    (hole.set_mp()).set_ori(xbh_new, 0., 0.) ;
    (star.set_mp()).set_ori(xns_new, yns_old, 0.) ;

    set_x_rot() = 0. ;

}

void Bin_bhns::rotation_axis_y(double thres_rot, double rot_exp_y,
			       double fact) {

    using namespace Unites ;

    double momunit = (hole.get_mass_bh()+star.mass_g_bhns())*omega*separ ;

    double error_x = line_mom_bhns()(0) / momunit ;

    if (error_x <= -1.) {
      cout << "Bin_bhns::rotation_axis:" << endl ;
      cout << "   !!! WARNING : error_x is smaller than -1 !!!" << endl ;
      error_x *= 0.1 ;
    }

    // Sets Y_NS
    // ---------

    double ff = pow( (2.+error_x)/(2.+2.*error_x), rot_exp_y) ;

    double yns_old = star.get_mp().get_ori_y() ;

    if ( fabs(error_x) < thres_rot ) {
        cout << "Bin_bhns::rotation_axis:" << endl ;
	cout << "   ff is set to 1 because error_x is smaller than" << endl ;
	cout << "   the threshold value (" << thres_rot << ")" << endl ;

        ff = 1. ;
    }

    cout << "Local center of mass of NS:" << endl ;
    cout << "   X_CM : " << xa_barycenter() / km << " [km]"
	 << "   Y_CM : " << ya_barycenter() / km << " [km]" << endl ;

    cout << "Bin_bhns::rotation_axis:" << endl ;
    cout << "   error_x : " << error_x << "   ff : " << ff << endl ;
    cout << "   old Y_BH : " << hole.get_mp().get_ori_y() / km << " [km]"
	 << "   old Y_NS : " << star.get_mp().get_ori_y() / km << " [km]"
	 << endl ;

    double aa = fact * separ ;
    double yns_new = yns_old + aa * (1. - ff) ;

    cout << "   new Y_BH : " << 0. / km << " [km]"
	 << "   new Y_NS : " << yns_new / km << " [km]"
	 << endl ;

    double xbh_old = hole.get_mp().get_ori_x() ;
    double xns_old = star.get_mp().get_ori_x() ;

    (hole.set_mp()).set_ori(xbh_old, 0., 0.) ;
    (star.set_mp()).set_ori(xns_old, yns_new, 0.) ;

    set_y_rot() = 0. ;

}
