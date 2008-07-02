/*
 *  Methods of class Hole_bhns to compute a forth-order Runge-Kutta
 *  integration to the phi direction for the solution of the Killing vectors
 *
 *    (see file hole_bhns.h for documentation).
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

char hole_bhns_rk_phi_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2008/07/02 20:47:55  k_taniguchi
 * Typos removed.
 *
 * Revision 1.1  2008/05/15 19:10:31  k_taniguchi
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
#include "hole_bhns.h"
#include "unites.h"
#include "utilitaires.h"

          //--------------------------------------------------//
          //      Forth-order Runge-Kutta on the equator      //
          //--------------------------------------------------//

Tbl Hole_bhns::runge_kutta_phi(const Tbl& xi_i, const double& phi_i,
			       const int& nrk_phi) const {

    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    int np = mg->get_np(1) ;

    Tbl xi_f(3) ;  // xi_f(0)=xi_hat{theta}, xi_f(1)=xi_hat{phi}, xi_f(2)=L
    xi_f.set_etat_qcq() ;

    if (kerrschild) {

      cout << "Not yet prepared!!!" << endl ;
      abort() ;

    }
    else {  // Isotropic coordinates

      // Initial data at phi=0 on the equator
      double xi_t0 = xi_i(0) ;  // xi_hat{theta}
      double xi_p0 = xi_i(1) ;  // xi_hat{phi}
      double xi_l0 = xi_i(2) ;  // L
      double phi0 = phi_i ;

      double dp = 2. * M_PI / double(np) / double(nrk_phi) ;

      double rah = rad_ah() ;

      Scalar dlnconfo(mp) ;
      dlnconfo = confo_tot.dsdt() / confo_tot ;
      dlnconfo.std_spectral_base() ;

      Scalar laplnconfo(mp) ;
      laplnconfo = confo_tot.lapang() / confo_tot ;
      laplnconfo.std_spectral_base() ;

      Scalar confo2(mp) ;
      confo2 = confo_tot * confo_tot ;
      confo2.std_spectral_base() ;

      double xi_t1, xi_t2, xi_t3, xi_t4, xi_tf ;
      double xi_p1, xi_p2, xi_p3, xi_p4, xi_pf ;
      double xi_l1, xi_l2, xi_l3, xi_l4, xi_lf ;
      double f1, f2, f3, f4 ;
      double g1, g2, g3, g4 ;
      double h1, h2, h3, h4 ;

      // Forth-order Runge-Kutta
      // (nrk_phi times steps between two collocation points)
      // ----------------------------------------------------

      for (int i=0; i<nrk_phi; i++) {

	// First
	f1 = - xi_l0 * rah * confo2.val_point(rah, M_PI/2., phi0)
	  + 2. * xi_p0 * dlnconfo.val_point(rah, M_PI/2., phi0) ;
	g1 = -2. * xi_t0 * dlnconfo.val_point(rah, M_PI/2., phi0) ;
	h1 = (1. - 2.*laplnconfo.val_point(rah, M_PI/2., phi0)) * xi_t0
	  / rah / confo2.val_point(rah, M_PI/2., phi0) ;

	xi_t1 = dp * f1 ;
	xi_p1 = dp * g1 ;
	xi_l1 = dp * h1 ;

	// Second
	f2 = - (xi_l0+0.5*xi_l1) * rah
	  * confo2.val_point(rah, M_PI/2., phi0+0.5*dp)
	  + 2. * (xi_p0+0.5*xi_p1)
	  * dlnconfo.val_point(rah, M_PI/2., phi0+0.5*dp) ;
	g2 = -2. * (xi_t0+0.5*xi_t1)
	  * dlnconfo.val_point(rah, M_PI/2., phi0+0.5*dp) ;
	h2 = (1. - 2.*laplnconfo.val_point(rah, M_PI/2., phi0+0.5*dp))
	  * (xi_t0+0.5*xi_t1) / rah
	  / confo2.val_point(rah, M_PI/2., phi0+0.5*dp) ;

	xi_t2 = dp * f2 ;
	xi_p2 = dp * g2 ;
	xi_l2 = dp * h2 ;

	// Third
	f3 = - (xi_l0+0.5*xi_l2) * rah
	  * confo2.val_point(rah, M_PI/2., phi0+0.5*dp)
	  + 2. * (xi_p0+0.5*xi_p2)
	  * dlnconfo.val_point(rah, M_PI/2., phi0+0.5*dp) ;
	g3 = -2. * (xi_t0+0.5*xi_t2)
	  * dlnconfo.val_point(rah, M_PI/2., phi0+0.5*dp) ;
	h3 = (1. - 2.*laplnconfo.val_point(rah, M_PI/2., phi0+0.5*dp))
	  * (xi_t0+0.5*xi_t2) / rah
	  / confo2.val_point(rah, M_PI/2., phi0+0.5*dp) ;

	xi_t3 = dp * f3 ;
	xi_p3 = dp * g3 ;
	xi_l3 = dp * h3 ;

	// Forth
	f4 = - (xi_l0+xi_l3) * rah * confo2.val_point(rah, M_PI/2., phi0+dp)
	  + 2. * (xi_p0+xi_p3) * dlnconfo.val_point(rah, M_PI/2., phi0+dp) ;
	g4 = -2. * (xi_t0+xi_t3) * dlnconfo.val_point(rah, M_PI/2., phi0+dp) ;
	h4 = (1. - 2.*laplnconfo.val_point(rah, M_PI/2., phi0+dp))
	  * (xi_t0+xi_t3) / rah / confo2.val_point(rah, M_PI/2., phi0+dp) ;

	xi_t4 = dp * f4 ;
	xi_p4 = dp * g4 ;
	xi_l4 = dp * h4 ;

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
