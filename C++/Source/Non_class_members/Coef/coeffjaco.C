/*
 *   Copyright (c) 2007 Jean-Louis Cornou
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
 * $Id$
 * $Log$
 * Revision 1.2  2007/12/21 12:41:46  j_novak
 * Removed the #include<fftw3.h> not needed here. Corrected headers.
 *
 * Revision 1.1  2007/12/11 15:42:21  jl_cornou
 * Premiere version des fonctions liees aux polynomes de Jacobi(0,2)
 *
 * $Header$
 *
 */
 
#include "tbl.h"

Tbl jacobipointsgl(int) ;

double* coeffjaco(int n , double* ff) {

  Tbl jj = jacobipointsgl(n) ;
  double som ;
  int i,k ;
  double* aa = new double[n+1] ;

  for (k = 0 ; k < n ; k++) {
    som = 3*ff[0]*jj(k,0)/(jj(n,0)*jj(n,0)) ;
    for (i = 1 ; i < n+1 ; i++) {
      som += ff[i]*jj(k,i)/(jj(n,i)*jj(n,i)) ;
    }
    aa[k] = (2*k+3)/double(n*(n+3))*som ;
  }
  som = 3*ff[0]/jj(n,0) ;
    for (i = 1 ; i < n+1 ; i++) {
    som += ff[i]/jj(n,i) ;
    }
    aa[n]=som/double(n+3) ;
  return aa ;
}
