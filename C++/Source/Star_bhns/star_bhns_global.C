/*
 *  Methods of class Star_bhns to compute global quantities
 *
 *    (see file star_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005,2007 Keisuke Taniguchi
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:18:16  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:40  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2008/05/15 19:15:20  k_taniguchi
 * Change of a parameter.
 *
 * Revision 1.1  2007/06/22 01:31:24  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
//#include <>

// Lorene headers
#include "star_bhns.h"
#include "unites.h"
#include "utilitaires.h"

                    //-------------------------------//
                    //          Baryon mass          //
                    //-------------------------------//

namespace Lorene {
double Star_bhns::mass_b() const {

    using namespace Unites ;

    if (p_mass_b == 0x0) {    // a new computation is required

        Scalar dens = pow(confo_tot, 6.) * nbar * gam_euler ;
	dens.std_spectral_base() ;

	p_mass_b = new double( dens.integrale() ) ;
    }

    return *p_mass_b ;

    /*
    cout << "Star_bhns::mass_b() is not available !!!" << endl
	 << " --> Use Star_bhns::mass_b_bhns(kerrschild, mass_bh, sepa)"
	 << endl ;
    abort() ;
    */

}

double Star_bhns::mass_b_bhns(bool kerrschild, const double& mass_bh,
			      const double& sepa) const {

    using namespace Unites ;

    if (p_mass_b_bhns == 0x0) {    // a new computation is required

	Scalar tmp(mp) ;

	if (kerrschild) {

	    Scalar xx(mp) ;
	    xx = mp.x ;
	    xx.std_spectral_base() ;
	    Scalar yy(mp) ;
	    yy = mp.y ;
	    yy.std_spectral_base() ;
	    Scalar zz(mp) ;
	    zz = mp.z ;
	    zz.std_spectral_base() ;

	    double yns = mp.get_ori_y() ;

	    Scalar rr(mp) ;
	    rr = sqrt( (xx+sepa)*(xx+sepa) + (yy+yns)*(yy+yns) + zz*zz ) ;
	    rr.std_spectral_base() ;

	    tmp = sqrt(1. + 2.*ggrav*mass_bh/rr) ;

	}
	else { // Isotropic coordinates with the maximal slicing

	    tmp = 1. ;

	}
	tmp.std_spectral_base() ;

        Scalar dens = pow(confo_tot, 6.) * nbar * gam_euler * tmp ;
	dens.std_spectral_base() ;

	p_mass_b_bhns = new double( dens.integrale() ) ;
    }

    return *p_mass_b_bhns ;

}


                    //--------------------------------------//
                    //          Gravitational mass          //
                    //--------------------------------------//

double Star_bhns::mass_g() const {

    cout << "Star_bhns::mass_g() is not available !!!" << endl
	 << " --> Use Star_bhns::mass_g_bhns()"
	 << endl ;
    abort() ;

}

double Star_bhns::mass_g_bhns() const {

    // This mass is valid only for an isolated spherical star

    if (p_mass_g_bhns == 0x0) {    // a new computation is required

        Scalar dens = lapconf_tot * pow(confo_tot, 5.)
	  * (ener_euler + s_euler) ;

	dens.std_spectral_base() ;

	p_mass_g_bhns = new double( dens.integrale() ) ;

    }

    return *p_mass_g_bhns ;

}
}
