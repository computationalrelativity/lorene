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

char ope_helmholtz_minus_2d_non_dege_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <math.h>
#include <stdlib.h>

#include "proto.h"
#include "ope_elementary.h"


		//------------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

Matrice _helmholtz_minus_2d_non_dege_pas_prevu(const Matrice &lap, int) {
    cout << "Construction non degeneree pas prevue..." << endl ;
    abort() ;
    exit(-1) ;
    return lap ;
}



	     	//-------------------
	       //--  R_CHEB   -------
	      //--------------------

Matrice _helmholtz_minus_2d_non_dege_r_cheb (const Matrice &lap, int) {
    
    
  int n = lap.get_dim(0) ;
    
  Matrice res(n-2, n-2) ;
  res.set_etat_qcq() ;
  for (int i=0 ; i<n-2 ; i++)
    for (int j=0 ; j<n-2 ; j++)
      res.set(i, j) = lap(i, j+2) ;
  
  res.set_band(4, 4) ;
  res.set_lu() ;
  
  return res ;
}


	     	//-------------------
	       //--  R_CHEBU   -----
	      //-------------------
Matrice _helmholtz_minus_2d_non_dege_r_chebu_deux (const Matrice&) ;    
	      
Matrice _helmholtz_minus_2d_non_dege_r_chebu (const Matrice &lap,  int puis) {

    switch (puis) {
	case 2 :
	    return _helmholtz_minus_2d_non_dege_r_chebu_deux (lap) ;
	default :
	    abort() ;
	    exit(-1) ;
	    return Matrice(0, 0) ;
    }
}


// Cas dzpuis = 2 ;
Matrice _helmholtz_minus_2d_non_dege_r_chebu_deux (const Matrice &lap) {
    
   int n = lap.get_dim(0) ;
    
   Matrice res(n-1, n-1) ;
   res.set_etat_qcq() ;
   for (int i=0 ;i<n-1 ; i++)
     for (int j=0 ; j<n-1 ; j++)
       res.set(i, j) = lap(i, j+1) ;
   res.set_band(5, 3) ;
   res.set_lu() ;
   
   return res ;
}


void Ope_helmholtz_minus_2d::do_non_dege() const {
  if (ope_cl == 0x0)
    do_ope_cl() ;
  
  if (non_dege != 0x0)
    delete non_dege ;
  
  // Routines de derivation
  static Matrice (*helmholtz_minus_2d_non_dege[MAX_BASE])(const Matrice&, int);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      helmholtz_minus_2d_non_dege[i] = _helmholtz_minus_2d_non_dege_pas_prevu ;
    }
    // Les routines existantes
    helmholtz_minus_2d_non_dege[R_CHEB >> TRA_R] = _helmholtz_minus_2d_non_dege_r_cheb ;
    helmholtz_minus_2d_non_dege[R_CHEBU >> TRA_R] = _helmholtz_minus_2d_non_dege_r_chebu ;
  }
  non_dege = new Matrice(helmholtz_minus_2d_non_dege[base_r](*ope_cl, dzpuis)) ;
}
