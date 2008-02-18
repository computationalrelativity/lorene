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


char val_solh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2008/02/18 13:53:45  j_novak
 * Removal of special indentation instructions.
 *
 * Revision 1.2  2003/12/11 15:37:09  p_grandclement
 * sqrt(2) ----> sqrt(double(2))
 *
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
#include "type_parite.h"


		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _val_solh_pas_prevu (int, double, double) {

    cout << " Solution homogene pas prevue ..... : "<< endl ;
    abort() ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	
		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _val_solh_r_cheb (int l, double alpha, double beta) {
  
  double echelle = beta/alpha ;

  Tbl res(2, 4) ;
  res.set_etat_qcq() ;
  
  // Solution 1 : (x+echelle)^l
  res.set(0,0) = pow(1.+echelle, l) ;
  res.set(0,1) = pow(-1.+echelle, l) ;
  res.set(0,2) = pow(1.+echelle, l-1)*l/alpha ;
  res.set(0,3) = pow(-1.+echelle, l-1)*l/alpha ;
  
  // Solution 2 : 1./(x+echelle)^(l+1)
  res.set(1,0) = 1./pow(1.+echelle, l+1) ;
  res.set(1,1) = 1./pow(-1.+echelle, l+1) ;
  res.set(1,2) = -1./pow(1.+echelle, l+2)*(l+1)/alpha ;
  res.set(1,3) = -1./pow(-1.+echelle, l+2)*(l+1)/alpha ;

  res /= sqrt(double(2)) ;
  return res ;
}	
	
		//-------------------
	       //--  R_CHEBP  ------
	      //-------------------

Tbl _val_solh_r_chebp (int l, double alpha, double) {
  
  Tbl res(4) ;
  res.set_etat_qcq() ;
  
  // Solution : x^l
  res.set(0) = 1. ;
  res.set(1) = (l==0) ? 1. : 0. ;
  res.set(2) = 1./alpha*l ;
  res.set(3) = (l==1) ? 1 : 0 ;
  
  res /= sqrt(double(2)) ;
  return res ;
}
	
	
	      	//-------------------
	       //--  R_CHEBI   -----
	      //-------------------
	
Tbl _val_solh_r_chebi (int l, double alpha, double) {
        
  Tbl res(4) ;
  res.set_etat_qcq() ;
  
  // Solution : x^l
  res.set(0) = 1. ;
  res.set(1) = (l==0) ? 1. : 0 ;
  res.set(2) = 1./alpha*l ;
  res.set(3) = (l==1) ? 1 : 0. ;
  
  res /= sqrt(double(2)) ;
  return res ;
}
	
	
	
	       	//-------------------
	       //--  R_CHEBU   -----
	      //-------------------
	
Tbl _val_solh_r_chebu (int l, double alpha, double) {
   
  Tbl res(4) ;
  res.set_etat_qcq() ;
  
  // Solution : 1/(x-1)^(l+1)
  res.set(0) = 0 ; /// Not used
  res.set(1) = pow(-2., l+1)/sqrt(double(2)) ; 
  res.set(2) = 0. ; /// noy used
  res.set(3) = -alpha*(l+1)*pow(-2., l+2)/sqrt(double(2)) ;        

  return res ;
}
	
	
	
	
	      	//-------------------
	       //--  Fonction   ----
	      //-------------------
	      
	      
Tbl val_solh(int l, double alpha, double beta, int base_r) {

		// Routines de derivation
    static Tbl (*val_solh[MAX_BASE])(int, double, double) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    val_solh[i] = _val_solh_pas_prevu ;
	}
		// Les routines existantes
	val_solh[R_CHEB >> TRA_R] = _val_solh_r_cheb ;
	val_solh[R_CHEBU >> TRA_R] = _val_solh_r_chebu ;
	val_solh[R_CHEBP >> TRA_R] = _val_solh_r_chebp ;
	val_solh[R_CHEBI >> TRA_R] = _val_solh_r_chebi ;
    }
    
    Tbl res(val_solh[base_r](l, alpha, beta)) ;
    return res ;
}
