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

char ope_helmholtz_minus_2d_cl_C[] = "$Header$" ;

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
namespace Lorene {
Matrice _cl_helmholtz_minus_2d_pas_prevu (const Matrice & source, int) {
    cout << "Combinaison lineaire pas prevu..." << endl ;
    abort() ;
    exit(-1) ;
    return source;
}


		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Matrice _cl_helmholtz_minus_2d_r_cheb (const Matrice &source, int) {
   int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
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

		//-------------------
	       //--  R_CHEBU   -----
	      //-------------------

Matrice _cl_helmholtz_minus_2d_r_chebu_deux (const Matrice&) ;


Matrice _cl_helmholtz_minus_2d_r_chebu (const Matrice &source, int puis) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
    Matrice res(n, n) ;
    res.set_etat_qcq() ;
    
    switch (puis) {
	case 2 :
	    res = _cl_helmholtz_minus_2d_r_chebu_deux(source) ;
	    break ;
	default :
	    abort() ;
	    exit(-1) ;
    }
    
    return res ;
}


//Cas dzpuis == 2
Matrice _cl_helmholtz_minus_2d_r_chebu_deux (const Matrice &source) {

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
  
  Matrice bis(tilde) ;
  for (int i=0 ; i<n-4 ; i++)
    for (int j=0 ; j<n ; j++)
      bis.set(i, j) = (tilde(i, j)+tilde(i+1, j)) ;

  Matrice res (bis) ;
  for (int i=0 ; i<n-4 ; i++)
    for (int j=0 ; j<n ; j++)
      res.set(i, j) = (bis(i, j)-bis(i+1, j)) ;

  return res ;
}

void Ope_helmholtz_minus_2d::do_ope_cl() const {
  if (ope_mat == 0x0)
    do_ope_mat() ;

  if (ope_cl != 0x0)
    delete ope_cl ;
  
  // Routines de derivation
  static Matrice (*cl_helmholtz_minus_2d[MAX_BASE])(const Matrice&, int);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      cl_helmholtz_minus_2d[i] = _cl_helmholtz_minus_2d_pas_prevu ;
    }
    // Les routines existantes
    cl_helmholtz_minus_2d[R_CHEB >> TRA_R] = _cl_helmholtz_minus_2d_r_cheb ;
    cl_helmholtz_minus_2d[R_CHEBU >> TRA_R] = _cl_helmholtz_minus_2d_r_chebu ;
  }
  ope_cl = new Matrice(cl_helmholtz_minus_2d[base_r](*ope_mat, dzpuis)) ;
}


}
