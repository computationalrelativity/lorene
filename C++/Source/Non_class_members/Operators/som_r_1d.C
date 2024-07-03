/*
 *   Copyright (c) 2024 Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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
 * Functions for direct spectral summation in the radial direction and a 1D array.
 * 
 *   SYNOPSYS:
 *     double som_r_1d_XX(double* ti, int nr, double x)
 *
 *     ti : array of spectral coefficients
 *     nr : size of the above array
 *     x  : point where to compute the summation; x in [0, 1] or x in [-1, 1]
 */


// Headers C
#include <cstdlib>

#include "headcpp.h"
#include "proto.h"

namespace Lorene {
			//-----------------------------
			//--- Non-implemented cases ---
			//-----------------------------

  double som_r_1d_pas_prevu(const double*, int, double)
  {
    throw(invalid_argument("som_r_1d : r basis not implemented yet !")) ;
    return 0. ;
  }

			//--------------
			//--- R_CHEB ---
			//--------------

  double som_r_1d_cheb(const double* ti, int nr, double x)
  {
    int i;
    const double* pi = ti ;	    // pi points on the array of coefficients
    double resu = 0. ;
    
    double cheb0 = 1 ;
    resu += (*pi)*cheb0 ;
    pi++ ;
    double cheb1 = x ;
    if (nr > 1) {
      resu += (*pi)*cheb1 ;
      pi++ ;
    }
    for (i=2; i<nr; i++) {
      double cheb2 = 2*x* cheb1 - cheb0 ; // Chebyshev recursion formula
      resu += (*pi)*cheb2 ;
      pi++ ;
      cheb0 = cheb1 ;
      cheb1 = cheb2 ;
    }
    
    return resu ;
  }


			//---------------
			//--- R_CHEBP ---
			//---------------

  double som_r_1d_chebp(const double* ti, int nr, double x)
  {
    int i ;
    const double* pi = ti ;	    // pi points on the array of coefficients
    double resu = 0. ;
    
    double cheb0 = 1. ;
    resu += (*pi)*cheb0 ;
    pi++ ;
    double cheb1 = x ;
    
    for (i=2; i<2*nr-1; i++) {
      double cheb2 = 2*x*cheb1 - cheb0 ; // Chebyshev recursion formula
      cheb0 = cheb1 ;
      cheb1 = cheb2 ;
      if (i%2 == 0) {
	resu += (*pi)*cheb2 ;
	pi++ ;
      }
    }
    return resu ;
  }


			//---------------
			//--- R_CHEBI ---
			//---------------

  double som_r_1d_chebi(const double* ti, int nr, double x)
  {
    int i ;
    const double* pi = ti ;	 // pi points on the array of coefficients
    double resu = 0. ;
    
    double cheb0 = 1 ;
    double cheb1 = x ;
    resu += (*pi)*cheb1 ;
    pi++ ;
    
    for (i=2; i<2*nr; i++) {
      double cheb2 = 2*x*cheb1 - cheb0 ; // Chebyshev recursion formula
      cheb0 = cheb1 ;
      cheb1 = cheb2 ;
      if (i%2 == 1) {
	resu += (*pi)*cheb2 ;
	pi++ ;
      }
    }
    return resu ;
  }

                        //--------------
			//--- R_LEG  ---
			//--------------

  double som_r_1d_leg(const double* ti, const int nr, double x)
  {
    int i ;
    const double* pi = ti ;	   // pi points on the array of coefficients
    double resu = 0. ;
    
    double leg0 = 1. ;
    resu += (*pi)*leg0 ;
    pi++ ;
    double leg1 = x ;
    if (nr > 1)  {
      resu += (*pi)*leg1 ;
      pi++ ;
    }
    
    for (i=2; i<nr; i++) {
      // Legendre recursion formula
      double leg2 = (double(2*i-1)*x* leg1 - double(i-1)*leg0) / double(i) ;
      resu += (*pi)*leg2 ;
      pi++ ;
      leg0 = leg1 ;
      leg1 = leg2 ;
    }
    return resu ;
  }


			//--------------
			//--- R_LEGP ---
			//--------------

  double som_r_1d_legp(const double* ti, const int nr, double x)
  {
    int i ;
    const double* pi = ti ;	// pi points on the array of coefficients
    double resu = 0. ;
    
    double leg0 = 1. ;
    resu += (*pi)*leg0 ;
    pi++ ;
    double leg1 = x ;
    for (i=2; i<2*nr-1; i++) {
      // Legendre recursion formula
      double leg2 = (double(2*i-1)*x* leg1 - double(i-1)*leg0) / double(i) ;
      leg0 = leg1 ;
      leg1 = leg2 ;
      if (i%2 == 0) {
	resu += (*pi)*leg2 ;
	pi++ ;
      }
    }
    return resu ;
  }


			//--------------
			//--- R_LEGI ---
			//--------------

  double som_r_1d_legi(const double* ti, int nr, double x)
  {
    int i ;
    const double* pi = ti ;	 // pi points on the array of coefficients
    double resu = 0. ;
    
    double leg0 = 1. ;
    double leg1 = x ;
    resu += (*pi)*leg1 ;
    pi++ ;
    for (i=2; i<2*nr; i++)
      {
	// Legendre recursion formula
	double leg2 = (double(2*i-1)*x* leg1 - double(i-1)*leg0) / double(i) ;
	leg0 = leg1 ;
	leg1 = leg2 ;
	if (i%2 == 1) {
	  resu += (*pi)*leg2 ;
	  pi++ ;
	}
      }
    return resu ;
  }

			//----------------
			//--- R_JACO02 ---
			//----------------

  double som_r_1d_jaco02(const double* ti, int nr, double x)
  {
    int i;
    const double* pi = ti ;	 // pi points on the array of coefficients
    double resu = 0. ;
    
    // Jacobi(0,2) polynomial values at x, up to degree nr-1
    double* jaco = jacobi(nr-1,x) ;
    for (i=0 ; i<nr ; i++) {
      resu += (*pi) * jaco[i] ;
      pi++ ;  // R suivant
    }
    return resu ;
  }

} // end of namespace
