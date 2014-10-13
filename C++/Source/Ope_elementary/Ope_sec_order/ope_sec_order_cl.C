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

char ope_sec_order_cl_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <cmath>
#include <cstdlib>

#include "proto.h"
#include "ope_elementary.h"

// Version Matrice --> Matrice

		//-------------------
	       //-- Pas prevu   ----
	      //-------------------

namespace Lorene {
Matrice _cl_sec_order_pas_prevu (const Matrice& so) {
  cout << "CL Sec_order not implemented" << endl ;
  abort() ;
  exit(-1) ;
  return so;
}

		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Matrice _cl_sec_order_r_cheb (const Matrice& source) {

  int n = source.get_dim(0) ;
  assert (n==source.get_dim(1)) ;
  
  Matrice barre(source) ;
  int dirac = 1 ;
  for (int i=0 ; i<n-2 ; i++) {
    for (int j=0 ; j<n ; j++)
      barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j))
	/(i+1) ;
    if (i==0) dirac = 0 ;
  }
  
  Matrice res(barre) ;
  for (int i=0 ; i<n-4 ; i++)
    for (int j=0 ; j<n ; j++)
      res.set(i, j) = barre(i, j)-barre(i+2, j) ;
  
  return res ;
} 


void Ope_sec_order::do_ope_cl() const {
  if (ope_mat == 0x0)
    do_ope_mat() ;
  
  if (ope_cl != 0x0)
    delete ope_cl ;
  
  // Routines de derivation
  static Matrice (*cl_sec_order[MAX_BASE])(const Matrice&);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      cl_sec_order[i] = _cl_sec_order_pas_prevu ;
    }
    // Les routines existantes
    cl_sec_order[R_CHEB >> TRA_R] = _cl_sec_order_r_cheb ;
  }
  ope_cl = new Matrice(cl_sec_order[base_r](*ope_mat)) ;
}


}
