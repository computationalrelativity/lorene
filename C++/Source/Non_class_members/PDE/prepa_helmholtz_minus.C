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


char prepa_helmholtz_minus_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/01/15 09:15:37  p_grandclement
 * Modification and addition of the Helmholtz operators
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

#include "matrice.h"
#include "tbl.h"
#include "type_parite.h"
#include "indent.h"
#include "proto.h"




		//------------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

Matrice _prepa_helmholtz_minus_nondege_pas_prevu(const Matrice &so) {
  
    cout << "Unknown case for prepa_helmholtz_minus_nondege" << endl ;
    abort() ;
    exit(-1) ;
    return so;
}



	     	//-------------------
	       //--  R_CHEBP   -------
	      //--------------------

Matrice _prepa_helmholtz_minus_nondege_r_chebp (const Matrice &lap) {
    
 int n = lap.get_dim(0) ;

 int non_dege = 1 ;
 
 Matrice res(n-non_dege, n-non_dege) ;
 res.set_etat_qcq() ;
 for (int i=0 ; i<n-non_dege ; i++)
   for (int j=0 ; j<n-non_dege ; j++)
     res.set(i, j) = lap(i, j+non_dege) ;
 
 res.set_band (4,1) ;
 res.set_lu() ;
 return res ;
} 
  
	     	//-------------------
	       //--  R_CHEB   -------
	      //--------------------

Matrice _prepa_helmholtz_minus_nondege_r_cheb (const Matrice &lap) {
    
  int n = lap.get_dim(0) ;
  int non_dege = 2 ;
  
  Matrice res(n-non_dege, n-non_dege) ;
  res.set_etat_qcq() ;
  for (int i=0 ; i<n-non_dege ; i++)
    for (int j=0 ; j<n-non_dege ; j++)
      res.set(i, j) = lap(i, j+non_dege) ;
  
  res.set_band (4,4) ;
  res.set_lu() ;
  return res ;
} 

	     	//-------------------
	       //--  R_CHEBU  -------
	      //--------------------
Matrice _prepa_helmholtz_minus_nondege_r_chebu (const Matrice &lap) {
    
  int n = lap.get_dim(0) ;
  int non_dege = 1 ;
  
  Matrice res(n-non_dege, n-non_dege) ;
  res.set_etat_qcq() ;
  for (int i=0 ; i<n-non_dege ; i++)
    for (int j=0 ; j<n-non_dege ; j++)
      res.set(i, j) = lap(i, j+non_dege) ;
  
  res.set_band (5,3) ;
  res.set_lu() ;
  return res ;
} 

        	//-------------------
	       //--  Fonction   ----
	      //-------------------
	      
Matrice prepa_helmholtz_minus_nondege(const Matrice &ope, int base_r) {

  // Routines de derivation
  static Matrice (*prepa_helmholtz_minus_nondege[MAX_BASE])
    (const Matrice&) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      prepa_helmholtz_minus_nondege[i] = 
	_prepa_helmholtz_minus_nondege_pas_prevu ;
    }
    // Les routines existantes
    prepa_helmholtz_minus_nondege[R_CHEB >> TRA_R] = 
      _prepa_helmholtz_minus_nondege_r_cheb ;
    prepa_helmholtz_minus_nondege[R_CHEBU >> TRA_R] = 
      _prepa_helmholtz_minus_nondege_r_chebu ;
    prepa_helmholtz_minus_nondege[R_CHEBP >> TRA_R] = 
      _prepa_helmholtz_minus_nondege_r_chebp ;
  }
  
  Matrice res(prepa_helmholtz_minus_nondege[base_r](ope)) ;
  return res ;
}

