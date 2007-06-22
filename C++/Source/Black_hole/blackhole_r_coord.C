/*
 *  Method of class Black_hole to express the radial coordinate
 *  in Kerr-Schild coordinates by that in spatially isotropic coordinates
 *
 *    (see file blackhole.h for documentation).
 *
 */

/*
 *   Copyright (c) 2006 Keisuke Taniguchi
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

char blackhole_r_coord_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:20:33  k_taniguchi
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
#include "blackhole.h"
#include "unites.h"
#include "utilitaires.h"

// Local function
double gg(double, const double) ;

const Scalar Black_hole::r_coord(bool neumann, bool first) const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    const Mg3d* mg = mp.get_mg() ;
    int nz = mg->get_nzone() ;          // total number of domains
    int nr = mg->get_nr(0) ;
    int nt = mg->get_nt(0) ;
    int np = mg->get_np(0) ;

    double mass = ggrav * mass_bh ;

    Scalar r_iso(mp) ;
    r_iso = mp.r ;
    r_iso.std_spectral_base() ;

    Scalar r_are(mp) ;
    r_are = r_iso ;  // Initialization
    r_are.std_spectral_base() ;

    // Sets C/M^2 for each case of the lapse boundary condition
    // --------------------------------------------------------
    double cc ;

    if (neumann) {  // Neumann boundary condition
        if (first) {  // First condition
	  // d\alpha/dr = 0
	  // --------------
	  cc = 2. ;
	}
	else {  // Second condition
	  // d\alpha/dr = \alpha/(2 rah)
	  // ---------------------------
	  cc = 0.5 * (sqrt(17.) - 1.) ;

	  //	  cc = 0.5 * (sqrt(17.) + 1.) ;
	}
    }
    else {  // Dirichlet boundary condition
       if (first) {  // First condition
	 // \alpha = 1/2
	 // ------------
	 cc = 2. ;
       }
       else {  // Second condition
	 // \alpha = 1/sqrt(2.)
	 // -------------------
	 cc = 2. * sqrt(2.) ;
       }
    }

    int ll ;
    double diff ;
    double ratio ;
    double precis = 1.e-15 ;
    double dp ;
    double tmp ;
    double tr ;

    int nn = 1000 ;
    assert(nn%4 == 0) ;
    int mm = nn/4 ;
    double x1, x2, x3, x4, x5 ;
    double hh, integ ;

    // Boole's Rule (Newton-Cotes Integral) for integration
    // ----------------------------------------------------

    for (int l=1; l<nz; l++) {

      for (int i=0; i<nr; i++) {

	ratio = 1. ;
	dp = 10. ;
	tr = r_iso.val_grid_point(l,0,0,i) ;

	while ( dp > precis ) {

	  diff = 1. ;  // Initialization
	  ll = 0 ;
	  dp = 0.1 * dp ;

	  while ( diff > precis ) {

	    ll++ ;
	    tmp = ratio + ll * dp ;

	    double r_max = 2.*mass/tmp/tr ;

	    hh = r_max / double(nn) ;
	    integ = 0. ;

	    for (int n=0; n<mm; n++) {

	      x1 = hh * double(4*n) ;
	      x2 = hh * double(4*n+1) ;
	      x3 = hh * double(4*n+2) ;
	      x4 = hh * double(4*n+3) ;
	      x5 = hh * double(4*n+4) ;

	      integ += (hh/45.) * (14.*gg(x1,cc) + 64.*gg(x2,cc)
				   + 24.*gg(x3,cc) + 64.*gg(x4,cc)
				   + 14.*gg(x5,cc)) ;

	    }

	    diff = -log( tmp ) - integ ;

	    //	    cout << "diff: " << diff << "  x: " << tmp << endl ;

	  }

	  ratio += (ll - 1) * dp ;

	}

	for (int j=0; j<nt; j++) {
	  for (int k=0; k<np; k++) {

	    r_are.set_grid_point(l,k,j,i) = ratio ;

	  }
	}

	//	arrete() ;

      }
    }

    r_are.std_spectral_base() ;
    r_are.annule_domain(0) ;
    r_are.raccord(1) ;

    /*
    cout << "r_are:" << endl ;
    for (int l=0; l<nz; l++) {
      cout << r_are.val_grid_point(l,0,0,0) << "  "
	   << r_are.val_grid_point(l,0,0,nr-1) << endl ;
    }
    */

    return r_are ;

}

//*****************************************************************

double gg(double xx, const double cc) {

    double tcc2 = cc*cc/16. ;
    double tmp = sqrt(1. - xx + tcc2*pow(xx, 4.)) ;

    double resu = (-1. + tcc2 * pow(xx, 3.)) / tmp / (1. + tmp) ;

    return resu ;

}
