/*
 *   Copyright (c) 2004 Philippe Grandclement
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

char ope_helmholtz_minus_pseudo_1d_solh_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_sf_bessel.h>

#include "proto.h"
#include "ope_elementary.h"

		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _solh_helmholtz_minus_pseudo_1d_pas_prevu (int, int, 
					       double, double, double) {

    cout << " Solution homogene pas prevue ..... : "<< endl ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	
	       	//-------------------
	       //--  R_CHEBU   -----
	      //-------------------
	
Tbl _solh_helmholtz_minus_pseudo_1d_r_chebu (int n, int l, 
					     double masse, 
					     double alpha, double) {
    
 
    Tbl res(n) ;
    res.set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;

    for (int i=0 ; i<n-1 ; i++) {
      double air = 1./(alpha*(-1-cos(M_PI*i/(n-1)))) ;
      if ((l !=0) && (l!=1))
	coloc[i] = gsl_sf_bessel_kl_scaled (l-1, masse*air) / exp(masse*air) * air;
      else
	coloc[i] = exp(-masse*air) ;
    }
    coloc[n-1] = 0 ;

    cfrcheb(deg, deg, coloc, deg, coloc) ;
    for (int i=0 ; i<n ;i++)
      res.set(i) = coloc[i] ;

    delete [] coloc ;
    delete [] deg ;
    
    return res ;
}


Tbl Ope_helmholtz_minus_pseudo_1d::get_solh () const {

  // Routines de derivation
  static Tbl (*solh_helmholtz_minus_pseudo_1d[MAX_BASE]) (int, int, double, double, double) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solh_helmholtz_minus_pseudo_1d[i] = _solh_helmholtz_minus_pseudo_1d_pas_prevu ;
    }
    // Les routines existantes
    solh_helmholtz_minus_pseudo_1d[R_CHEBU >> TRA_R] = _solh_helmholtz_minus_pseudo_1d_r_chebu ;
  }
  
  Tbl res(solh_helmholtz_minus_pseudo_1d[base_r](nr, l_quant, masse, alpha, beta)) ;

  // Un peu tricky...
  
  if (res.get_ndim() == 1) {
    Tbl val_lim (val_solp (res, alpha, base_r)) ;
    val_lim *= sqrt (double(2)) ;

    s_one_plus   = val_lim(0) ;
    s_one_minus  = val_lim(1) ; 
    ds_one_plus  = val_lim(2) ;
    ds_one_minus = val_lim(3) ;

  }
  else {
    Tbl auxi (nr) ;
    auxi.set_etat_qcq() ;
    for (int i=0 ; i<nr ; i++)
      auxi.set(i) = res(0,i) ;

    Tbl val_one  (val_solp (auxi, alpha, base_r)) ; 
    val_one *= sqrt (double(2)) ;

    s_one_plus   = val_one(0) ;
    s_one_minus  = val_one(1) ; 
    ds_one_plus  = val_one(2) ;
    ds_one_minus = val_one(3) ;

    for (int i=0 ; i<nr ; i++)
      auxi.set(i) = res(1,i) ;

    Tbl val_two  (val_solp (auxi, alpha, base_r)) ;
    val_two *= sqrt(double(2)) ;

    s_two_plus   = val_two(0) ;
    s_two_minus  = val_two(1) ; 
    ds_two_plus  = val_two(2) ;
    ds_two_minus = val_two(3) ;   

  }
  
  return res ;
}
