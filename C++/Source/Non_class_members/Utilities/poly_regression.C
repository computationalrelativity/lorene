/*
 *  Functions to compute polynomial regression of 1 data.
 *
 */

/*
 *   Copyright (c) 2022 Jerome Novak
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
 * Revision 1.1  2022/04/15 13:39:25  j_novak
 * New class Eos_compose_fit to generate fitted EoSs from CompOSE tables.
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "matrice.h"

namespace Lorene {

  // Function to initialize the transformation matrix Chebyshev -> Taylor expansion
  //-------------------------------------------------------------------------------
  Tbl initialize_dd(int ncoef) {
    Tbl dd(ncoef, ncoef) ;
    dd.annule_hard() ;
    dd.set(0,0) = 1. ; dd.set(0, 1) = 0 ;
    dd.set(1, 0) = 0 ; dd.set(1, 1) = 1. ;
    for (int i=2; i<ncoef; i++) {
      dd.set(i,0) = 0 ;
      dd.set(i, 1) = 0 ;
    }
    for (int j=2; j<ncoef; j++) {
      dd.set(0, j) = -dd(0, j-2) ;
      for (int i=1; i<ncoef; i++) {
	dd.set(i, j) = 2.*dd(i-1, j-1) - dd(i, j-2) ; // Recursion formula
	// of Chebyshev polynomials
      }
    }
    
    return dd ;
  }
  
  
  // Computes the polynomial regression of degree n_poly to data (xx, yy)
  // Result is an array of Chebyshev coefficients.
  Tbl poly_regression(const Tbl& xx, const Tbl& yy, int n_poly) {
    
    assert((xx.get_ndim() == 1) && (yy.get_ndim() == 1)) ;
    assert( xx.get_taille() == yy.get_taille() ) ;
    
    int n_data = xx.get_dim(0) ;
    
    Tbl powxj(n_data, 2*n_poly+1) ; powxj.set_etat_qcq() ;
    for (int i=0; i<n_data; i++)
      powxj.set(i, 0) = 1. ;
    for (int j=1; j<2*n_poly+1; j++) {
      for (int i=0; i<n_data; i++) {
	powxj.set(i,j) = powxj(i, j-1)*xx(i) ;
      }
    }
    
    Tbl xpow(2*n_poly+1) ; xpow.annule_hard() ;
    for (int j=0; j<2*n_poly+1; j++) {
      for (int i=0; i<n_data; i++) {
	xpow.set(j) += powxj(i, j) ;
      }
    }
    
    Tbl rhs(n_poly+1) ; rhs.annule_hard() ;
    for (int j=0; j<n_poly+1; j++) {
      for (int i=0; i<n_data; i++) {
	rhs.set(j) += powxj(i, j)*yy(i) ;
      }
    }
    
    Matrice mat(n_poly+1, n_poly+1) ; mat.set_etat_qcq() ;
    mat.set(0,0) = n_data ;
    for (int i=0; i<=n_poly; i++) {
      for (int j=0; j<=n_poly; j++) {
	mat.set(i,j) = xpow(i+j) ;
      }
    }
    mat.set_lu() ;
    
    Tbl sol_Taylor =  mat.inverse(rhs) ; // Solution expressed as Taylor coefficients
    
    Matrice cheb = initialize_dd(n_poly+1) ;
    cheb.set_lu() ;
    
    return cheb.inverse(sol_Taylor) ;
  }

}
