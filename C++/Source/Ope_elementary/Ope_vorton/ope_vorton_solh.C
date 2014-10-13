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

char ope_vorton_solh_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <cmath>
#include <cstdlib>

#include "proto.h"
#include "ope_elementary.h"

		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
namespace Lorene {
Tbl _solh_vorton_pas_prevu (int, int, double, double, Tbl&) {

    cout << " Solution homogene pas prevue ..... : "<< endl ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
		//-------------------
	       //--  R_CHEBU   -----
	      //-------------------

Tbl _solh_vorton_r_chebu (int n, int l, double alpha, double, Tbl& val_lim) {

  double l_one = -double(l) ;
  double rminus = -0.5/alpha ;

  Tbl res(n) ;
  res.set_etat_qcq() ;
  double* coloc = new double[n] ;
  
  int * deg = new int[3] ;
  deg[0] = 1 ; 
  deg[1] = 1 ;
  deg[2] = n ;
  
  //Construction de la premiere solution homogene :
  for (int i=0 ; i<n ; i++)
    coloc[i] = pow(1./alpha/(-cos(M_PI*i/(n-1))-1) , l_one) ;
	
  cfrcheb(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(i) = coloc[i] ;

  delete [] coloc ;
  delete [] deg ;
   
  val_lim.set(0,0) = pow(rminus, l_one) ;
  val_lim.set(0,1) = l_one*pow(rminus, l_one-1.) ;
  val_lim.set(0,2) = 0. ;
  val_lim.set(0,3) = 0. ;
  val_lim  /= sqrt(double(2)) ;

  return res ;
}


		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _solh_vorton_r_cheb (int n, int l, double alpha, double beta, Tbl& val_lim) {


  double l_one = double(l+1) ;
  double l_two = double(-l) ;
  double rminus = beta - alpha ;
  double rplus = beta + alpha ;

  Tbl res(2, n) ;
  res.set_etat_qcq() ;
  double* coloc = new double[n] ;
  
  int * deg = new int[3] ;
  deg[0] = 1 ; 
  deg[1] = 1 ;
  deg[2] = n ;
  
  //Construction de la premiere solution homogene :
  for (int i=0 ; i<n ; i++)
    coloc[i] = pow(alpha*(-cos(M_PI*i/(n-1))) + beta, l_one) ;

  cfrcheb(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(0, i) = coloc[i] ;
    
  // construction de la seconde solution homogene :
  for (int i=0 ; i<n ; i++) 
    coloc[i] = pow(alpha*(-cos(M_PI*i/(n-1))) + beta, l_two) ;

  cfrcheb(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(1, i) = coloc[i] ;	
        
  delete [] coloc ;
  delete [] deg ;
   
  val_lim.set(0,0) = pow(rminus, l_one) ;
  val_lim.set(0,1) = l_one*pow(rminus, l_one-1) ;
  val_lim.set(0,2) = pow(rplus, l_one) ;
  val_lim.set(0,3) = l_one*pow(rplus, l_one-1) ;

  val_lim.set(1,0) = pow(rminus, l_two) ;
  val_lim.set(1,1) = l_two*pow(rminus, l_two-1) ;
  val_lim.set(1,2) = pow(rplus, l_two) ;
  val_lim.set(1,3) = l_two*pow(rplus, l_two-1) ;
  val_lim  /= sqrt(double(2)) ;

  return res ;
}


Tbl Ope_vorton::get_solh () const {

  // Routines de derivation
  static Tbl (*solh_vorton[MAX_BASE]) (int, int, double, double, Tbl&) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solh_vorton[i] = _solh_vorton_pas_prevu ;
    }
    // Les routines existantes
    solh_vorton[R_CHEB >> TRA_R] = _solh_vorton_r_cheb ;
    solh_vorton[R_CHEBU >> TRA_R] = _solh_vorton_r_chebu ;
  }

  Tbl val_lim (2 ,4) ;
  val_lim.set_etat_qcq() ;
  Tbl res(solh_vorton[base_r](nr,l_quant, alpha, beta, val_lim)) ;


  s_one_minus  = val_lim(0,0) ;
  ds_one_minus = val_lim(0,1) ;
  s_one_plus   = val_lim(0,2) ;
  ds_one_plus  = val_lim(0,3) ;

  if (res.get_ndim()>1) {
     s_two_minus  = val_lim(1,0) ;
     ds_two_minus = val_lim(1,1) ;
     s_two_plus   = val_lim(1,2) ;
     ds_two_plus  = val_lim(1,3) ;
  }

  return res ;
}
}
