/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char solh_helmholtz_minusC[] = "$Header $" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/12/11 14:48:49  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * 
 * $Header$
 *
 */

//fichiers includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "proto.h"
#include "matrice.h"
#include "tbl.h"
#include "type_parite.h"
#include "indent.h"

	//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _solh_helmholtz_minus_pas_prevu (int, double, double, double) {

  cout << "Homogeneous solution not implemented in hemlholtz_minus : "<< endl ;
  abort() ;
  exit(-1) ;
  Tbl res(1) ;
  return res;
}
	

	
		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _solh_helmholtz_minus_r_cheb (int n, double alpha, double beta, 
				  double masse) {
  
  assert (masse > 0) ;
  
  Tbl res(2,n) ;
  res.set_etat_qcq() ;
  double* coloc = new double[n] ;
    
  int * deg = new int[3] ;
  deg[0] = 1 ; 
  deg[1] = 1 ;
  deg[2] = n ;
  
  // First SH
  for (int i=0 ; i<n ; i++){
    double air = alpha*(-cos(M_PI*i/(n-1))) + beta ;
    coloc[i] = exp(masse*air)/air ;
  }

  cfrcheb(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(0,i) = coloc[i] ;

  // Second SH
  for (int i=0 ; i<n ; i++){
    double air = alpha*(-cos(M_PI*i/(n-1))) + beta ;
    coloc[i] = exp(-masse*air)/air ;
  }
  
  cfrcheb(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(1,i) = coloc[i] ;

  delete [] deg ;
  delete [] coloc ;
  return res ;
}
	
		//-------------------
	       //--  R_CHEBU  ------
	      //-------------------

Tbl _solh_helmholtz_minus_r_chebu (int n, double alpha, double, double masse) {
  
  assert (masse > 0) ;
  
  Tbl res(n) ;
  res.set_etat_qcq() ;
  double* coloc = new double[n] ;
    
  int * deg = new int[3] ;
  deg[0] = 1 ; 
  deg[1] = 1 ;
  deg[2] = n ;
  
  for (int i=0 ; i<n-1 ; i++){
    double air = 1./(alpha*(-1-cos(M_PI*i/(n-1)))) ;
    coloc[i] = exp(-masse*air)/air ;
  }
  coloc[n-1] = 0 ;
  
  cfrcheb(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(i) = coloc[i] ;

  delete [] deg ;
  delete [] coloc ;
  return res ;
}


	      	//-------------------
	       //--  Fonction   ----
	      //-------------------
	      
	      
Tbl solh_helmholtz_minus (int n, double alpha, double beta, 
			  double masse, int base_r) {

  // Routines de derivation
  static Tbl (*solh_helmholtz_minus[MAX_BASE])(int, double, double, double) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solh_helmholtz_minus[i] = _solh_helmholtz_minus_pas_prevu ;
    }
    // Les routines existantes
    solh_helmholtz_minus[R_CHEB >> TRA_R] = _solh_helmholtz_minus_r_cheb ;
    solh_helmholtz_minus[R_CHEBU >> TRA_R] = _solh_helmholtz_minus_r_chebu ;
  }
    
  Tbl res(solh_helmholtz_minus[base_r](n, alpha, beta, masse)) ;
  return res ;
}
