/*
 *   Copyright (c) 2003 Philippe Grandclement
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

 

/*
 * $Id$
 * $Header$
 *
 */
#include <cmath>
#include <cstdlib>

#include "proto.h"
#include "ope_elementary.h"

		//-------------------
	       //-- Pas prevu   ----
	      //-------------------

namespace Lorene {
Matrice _sec_order_non_dege_pas_prevu (const Matrice& so) {
  cout << "Sec_order non dege : not implemented" << endl ;
  abort() ;
  exit(-1) ;
  return so;
}

		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Matrice _sec_order_non_dege_r_cheb (const Matrice& source) {

  int n = source.get_dim(0) ;
  int non_dege = 2 ;
  
  Matrice res(n-non_dege, n-non_dege) ;
  res.set_etat_qcq() ;
  for (int i=0 ; i<n-non_dege ; i++)
    for (int j=0 ; j<n-non_dege ; j++)
      res.set(i, j) = source(i, j+non_dege) ;
  
  res.set_band (2,2) ;
  res.set_lu() ;
  return res ;
} 


void Ope_sec_order::do_non_dege() const {
  if (ope_cl == 0x0)
    do_ope_cl() ;
  
  if (non_dege != 0x0)
    delete non_dege ;
  
  // Routines de derivation
  static Matrice (*sec_order_non_dege[MAX_BASE])(const Matrice&);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      sec_order_non_dege[i] = _sec_order_non_dege_pas_prevu ;
    }
    // Les routines existantes
    sec_order_non_dege[R_CHEB >> TRA_R] = _sec_order_non_dege_r_cheb ;
  }
  non_dege = new Matrice(sec_order_non_dege[base_r](*ope_cl)) ;
}
}
