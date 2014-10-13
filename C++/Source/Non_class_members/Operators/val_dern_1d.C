/*
 *   Copyright (c) 2004 Jerome Novak
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


char val_dern_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/13 08:53:27  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2004/02/17 09:21:39  j_novak
 * New functions for calculating values of the derivatives of a function
 * using its Chebyshev coefficients.
 *
 *
 * $Header$
 *
 */

#include "type_parite.h"
#include "tbl.h"

/*
 * Functions computing value of f^(n) at boundaries of the interval [-1, 1],
 * using the Chebyshev expansion of f. Note: n=0 works too.
 * 
 * Input : 1-dimensional Tbl containing the Chebyshev coefficients of f.
 *	    int base : base of spectral expansion.
 *
 * Output : double : the value of the n-th derivative of f at x=+/- 1.
 * 
 */


namespace Lorene {
double val1_dern_1d(int n, const Tbl& tb, int base_r)
{

  //This function should be OK for any radial base
  assert ( (base_r == R_CHEB) || (base_r == R_CHEBI) || (base_r == R_CHEBP) ||
	   (base_r == R_CHEBU) ) ; 

  assert (n>=0) ;
  assert (tb.get_ndim() == 1) ;
  int nr = tb.get_dim(0) ;

  double resu = 0. ;

  int n_ini = ( (base_r == R_CHEBP) || (base_r == R_CHEBI) ) ? n / 2 : n ;

  double *tbi = &tb.t[n_ini] ;
  for (int i=n_ini; i<nr; i++) {
    double fact = 1. ;
    int ii = i ;
    if (base_r == R_CHEBP) ii *= 2 ;
    if (base_r == R_CHEBI) ii = 2*i + 1 ;
    for (int j=0; j<n; j++) 
      fact *= double(ii*ii - j*j)/double(2*j + 1) ;
    resu += fact * (*tbi) ;
    tbi++ ;
  }

  return resu ;
}

double valm1_dern_1d(int n, const Tbl& tb, int base_r)
{

  //This function should be OK for any radial base
  assert ( (base_r == R_CHEB) || (base_r == R_CHEBI) || (base_r == R_CHEBP) ||
	   (base_r == R_CHEBU) ) ; 

  assert (n>=0) ;
  assert (tb.get_ndim() == 1) ;
  int nr = tb.get_dim(0) ;

  double resu = 0. ;
  double parite, fac ;
  int n_ini ;
  switch (base_r) {
  case R_CHEBP:
    n_ini = n / 2 ;
    parite = 1 ;
    fac = (n%2 == 0 ? 1 : -1) ;
    break ;
  case R_CHEBI: 
    n_ini = n / 2 ;
    fac = (n%2 == 0 ? -1 : 1) ;
    parite = 1 ;
    break ;
  default:
    n_ini = n ;
    parite = -1 ;
    fac = 1 ;
    break ;
  }
  double *tbi = &tb.t[n_ini] ;
  
  for (int i=n_ini; i<nr; i++) {
    double fact = fac ;
    int ii = i ;
    if (base_r == R_CHEBP) ii *= 2 ;
    if (base_r == R_CHEBI) ii = 2*i + 1 ;
    for (int j=0; j<n; j++)
      fact *= double(ii*ii - j*j)/double(2*j + 1) ;
    resu += fact * (*tbi) ;
    fac *= parite ;
    tbi++ ;
  }

  return resu ;
}
}
