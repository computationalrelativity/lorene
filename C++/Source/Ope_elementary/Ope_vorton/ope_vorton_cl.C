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

char ope_vorton_cl_C[] = "$Header$" ;

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

Matrice _cl_vorton_pas_prevu (const Matrice& so, int) {
  cout << "CL vorton not implemented" << endl ;
  abort() ;
  exit(-1) ;
  return so;
}

		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Matrice _cl_vorton_r_cheb (const Matrice& source, int) {

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
	for (int j=0 ; j< n ; j++)
	    res.set(i, j) = barre(i, j)-barre(i+2, j) ;
  return res ;
}

	        //-------------------
	       //--  R_CHEBU   -----
	      //-------------------

Matrice _cl_vorton_r_chebu_trois (const Matrice &source) {

    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ; 
 
    Matrice barre(source) ;
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=0 ; j<n ; j++)
	     barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j)) ;
	if (i==0) dirac = 0 ;
	}
    
    Matrice tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    tilde.set(i, j) = (barre(i, j)-barre(i+2, j)) ;
    
    Matrice res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = (tilde(i, j)+tilde(i+1, j)) ;

    return res ;
}
    
Matrice _cl_vorton_r_chebu (const Matrice &source, int puis) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
    Matrice res(n, n) ;
    res.set_etat_qcq() ;
    
    switch (puis) {
	case 3 :
	    res = _cl_vorton_r_chebu_trois(source) ;
	    break ;
	default :
	    abort() ;
	    exit(-1) ;
      }
return res ;
}


void Ope_vorton::do_ope_cl() const {
  if (ope_mat == 0x0)
    do_ope_mat() ;
  
  if (ope_cl != 0x0)
    delete ope_cl ;
  
  // Routines de derivation
  static Matrice (*cl_vorton[MAX_BASE])(const Matrice&, int);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      cl_vorton[i] = _cl_vorton_pas_prevu ;
    }
    // Les routines existantes
    cl_vorton[R_CHEB >> TRA_R] = _cl_vorton_r_cheb ;
    cl_vorton[R_CHEBU >> TRA_R] = _cl_vorton_r_chebu ;
  }
  ope_cl = new Matrice(cl_vorton[base_r](*ope_mat, dzpuis)) ;
}


