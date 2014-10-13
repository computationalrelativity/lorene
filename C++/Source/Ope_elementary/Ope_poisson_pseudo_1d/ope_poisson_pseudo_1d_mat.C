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

char ope_poisson_pseudo_1d_mat_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <cmath>
#include <cstdlib>

#include "proto.h"
#include "ope_elementary.h"

		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

namespace Lorene {
Matrice _poisson_pseudo_1d_mat_pas_prevu(int, int, double, double) {
    cout << "Operator pas prevu..." << endl ;
    abort() ;
    exit(-1) ;
    Matrice res(1, 1) ;
    return res;
}


		   //-------------------------
		   //--   CAS R_CHEBP    -----
		   //--------------------------
		    

Matrice _poisson_pseudo_1d_mat_r_chebp (int n, int l, double, double) {

  Matrice dd(n, n) ;
  dd.set_etat_qcq() ;
  Matrice xx(n, n) ;
  xx.set_etat_qcq() ;
  
  double* vect  = new double[n] ;
    
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 1 ;
    d2sdx2_1d (n, &vect, R_CHEBP) ;
    
    for (int j=0 ; j<n ; j++)
      dd.set(j, i) = vect[j] ; 
  }
  
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 1 ;
    sx2_1d (n, &vect, R_CHEBP) ;
    for (int j=0 ; j<n ; j++)
      xx.set(j, i) = vect[j] ;
  }
  
  delete [] vect ;
  
  Matrice res(n, n) ;
  res = dd-l*(l-1)*xx ;
	
  return res ;
}


		   //------------------------
		   //--   CAS R_CHEBI    ----
		   //------------------------
		    

Matrice _poisson_pseudo_1d_mat_r_chebi (int n, int l,
					double, double) {
  
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEBI) ;  // appel dans le cas impair
	for (int j=0 ; j<n ; j++)
	    dd.set(j, i) = vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sx2_1d (n, &vect, R_CHEBI) ;
	for (int j=0 ; j<n ; j++)
	    xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd-l*(l-1)*xx ;
    return res ;
} 

		   //-------------------------
		   //--   CAS R_CHEB    -----
		   //-----------------------
		    

Matrice _poisson_pseudo_1d_mat_r_cheb (int n, int l, 
				       double alf, double bet) {
            
  double echelle = bet / alf ;

  
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
	for (int j=0 ; j<n ; j++)
	    dd.set(j, i) = vect[j]*echelle*echelle ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
	multx_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<(n>i+1 ? i+1 : n) ; j++)
	    dd.set(j, i) += 2*echelle*vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
	multx_1d (n, &vect, R_CHEB) ;
	multx_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<(n>i+1 ? i+1 : n) ; j++)
	    dd.set(j, i) += vect[j] ;
    }
	   
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sx2_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<n ; j++)
	    xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd-l*(l-1)*xx ;
    return res ;
} 


void Ope_poisson_pseudo_1d::do_ope_mat() const {
  if (ope_mat != 0x0) 
    delete ope_mat ;

  // Routines de derivation
  static Matrice (*poisson_pseudo_1d_mat[MAX_BASE])(int, int, double, double);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      poisson_pseudo_1d_mat[i] = _poisson_pseudo_1d_mat_pas_prevu ;
    }
    // Les routines existantes
    poisson_pseudo_1d_mat[R_CHEB >> TRA_R] = _poisson_pseudo_1d_mat_r_cheb ;
    poisson_pseudo_1d_mat[R_CHEBP >> TRA_R] = _poisson_pseudo_1d_mat_r_chebp ;
    poisson_pseudo_1d_mat[R_CHEBI >> TRA_R] = _poisson_pseudo_1d_mat_r_chebi ;
  }
  ope_mat = new Matrice(poisson_pseudo_1d_mat[base_r](nr, l_quant, 
						      alpha, beta)) ;
}
}
