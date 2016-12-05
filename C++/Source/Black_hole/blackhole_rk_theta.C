/*
 *  Methods of class Black_hole to compute a forth-order Runge-Kutta
 *  integration to the theta direction for the solution of the Killing vectors
 *
 *    (see file blackhole.h for documentation).
 *
 */

/*
 *   Copyright (c) 2007 Keisuke Taniguchi
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
 * Revision 1.5  2016/12/05 16:17:48  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:52:46  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:03  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/07/02 20:44:19  k_taniguchi
 * Typos removed.
 *
 * Revision 1.1  2008/05/15 19:33:50  k_taniguchi
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
#include "blackhole.h"
#include "unites.h"
#include "utilitaires.h"

          //---------------------------------------------------//
          //     Forth-order Runge-Kutta to the polar angle    //
          //---------------------------------------------------//

namespace Lorene {
Tbl Black_hole::runge_kutta_theta_bh(const Tbl& xi_i, const double& theta_i,
				     const double& phi,
				     const int& nrk_theta) const {

    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    int nt = mg->get_nt(1) ;

    Tbl xi_f(3) ;  // xi_f(0)=xi_hat{theta}, xi_f(1)=xi_hat{phi}, xi_f(2)=L
    xi_f.set_etat_qcq() ;

    if (kerrschild) {

      cout << "Not yet prepared!!!" << endl ;
      abort() ;

    }
    else {  // Isotropic coordinates

      // Initial data at phi on the equator
      double xi_t0 = xi_i(0) ;  // xi_hat{theta}
      double xi_p0 = xi_i(1) ;  // xi_hat{phi}
      double xi_l0 = xi_i(2) ;  // L
      double theta0 = theta_i ;

      double dt = - 0.5 * M_PI / double(nt-1) / double(nrk_theta) ;
      // Compute from M_PI/2 to 0

      double rah = rad_ah() ;

      Scalar dlnconfo(mp) ;
      dlnconfo = confo.stdsdp() / confo ;
      dlnconfo.std_spectral_base() ;

      Scalar laplnconfo(mp) ;
      laplnconfo = confo.lapang() / confo ;
      laplnconfo.std_spectral_base() ;

      Scalar confo2(mp) ;
      confo2 = confo * confo ;
      confo2.std_spectral_base() ;

      double xi_t1, xi_t2, xi_t3, xi_t4, xi_tf ;
      double xi_p1, xi_p2, xi_p3, xi_p4, xi_pf ;
      double xi_l1, xi_l2, xi_l3, xi_l4, xi_lf ;
      double f1, f2, f3, f4 ;
      double g1, g2, g3, g4 ;
      double h1, h2, h3, h4 ;

      // Forth-order Runge-Kutta
      // (nrk_theta times steps between two collocation points)
      // ------------------------------------------------------

      for (int i=0; i<nrk_theta; i++) {

	// First
	f1 = -2. * xi_p0 * dlnconfo.val_point(rah, theta0, phi) ;
	g1 = xi_l0 * rah * confo2.val_point(rah, theta0, phi)
	  + 2. * xi_t0 * dlnconfo.val_point(rah, theta0, phi) ;
	h1 = - (1. - 2.*laplnconfo.val_point(rah, theta0, phi)) * xi_p0
	  / rah / confo2.val_point(rah, theta0, phi) ;

	xi_t1 = dt * f1 ;
	xi_p1 = dt * g1 ;
	xi_l1 = dt * h1 ;

	// Second
	f2 = -2. * (xi_p0+0.5*xi_p1)
	  * dlnconfo.val_point(rah, theta0+0.5*dt, phi) ;
	g2 = (xi_l0+0.5*xi_l1) * rah
	  * confo2.val_point(rah, theta0+0.5*dt, phi)
	  + 2. * (xi_t0+0.5*xi_t1)
	  * dlnconfo.val_point(rah, theta0+0.5*dt, phi) ;
	h2 = - (1. - 2.*laplnconfo.val_point(rah, theta0+0.5*dt, phi))
	  * (xi_p0+0.5*xi_p1) / rah
	  / confo2.val_point(rah, theta0+0.5*dt, phi) ;

	xi_t2 = dt * f2 ;
	xi_p2 = dt * g2 ;
	xi_l2 = dt * h2 ;

	// Third
	f3 = -2. * (xi_p0+0.5*xi_p2)
	  * dlnconfo.val_point(rah, theta0+0.5*dt, phi) ;
	g3 = (xi_l0+0.5*xi_l2) * rah
	  * confo2.val_point(rah, theta0+0.5*dt, phi)
	  + 2. * (xi_t0+0.5*xi_t2)
	  * dlnconfo.val_point(rah, theta0+0.5*dt, phi) ;
	h3 = - (1. - 2.*laplnconfo.val_point(rah, theta0+0.5*dt, phi))
	  * (xi_p0+0.5*xi_p2) / rah
	  / confo2.val_point(rah, theta0+0.5*dt, phi) ;

	xi_t3 = dt * f3 ;
	xi_p3 = dt * g3 ;
	xi_l3 = dt * h3 ;

	// Forth
	f4 = -2. * (xi_p0+xi_p3) * dlnconfo.val_point(rah, theta0+dt, phi) ;
	g4 = (xi_l0+xi_l3) * rah * confo2.val_point(rah, theta0+dt, phi)
	  + 2. * (xi_t0+xi_t3) * dlnconfo.val_point(rah, theta0+dt, phi) ;
	h4 = - (1. - 2.*laplnconfo.val_point(rah, theta0+dt, phi))
	  * (xi_p0+xi_p3) / rah / confo2.val_point(rah, theta0+dt, phi) ;

	xi_t4 = dt * f4 ;
	xi_p4 = dt * g4 ;
	xi_l4 = dt * h4 ;

	// Final results
	// -------------
	xi_tf = xi_t0 + (xi_t1 + 2.*xi_t2 + 2.*xi_t3 + xi_t4) / 6. ;
	xi_pf = xi_p0 + (xi_p1 + 2.*xi_p2 + 2.*xi_p3 + xi_p4) / 6. ;
	xi_lf = xi_l0 + (xi_l1 + 2.*xi_l2 + 2.*xi_l3 + xi_l4) / 6. ;

	// Final results are put into the initial data
	// in order for the next step
	// -------------------------------------------
	xi_t0 = xi_tf ;
	xi_p0 = xi_pf ;
	xi_l0 = xi_lf ;

      } // End of the loop

      xi_f.set(0) = xi_tf ;
      xi_f.set(1) = xi_pf ;
      xi_f.set(2) = xi_lf ;

    }

    return xi_f ;

}
}
