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

char ope_helmholtz_minus_pseudo_1d_mat_C[] = "$Header$" ;

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

Matrice _helmholtz_minus_pseudo_1d_mat_pas_prevu(int, int, double, double, 
						 double, int) {
    cout << "Operateur pas prevu..." << endl ;
    abort() ;
    exit(-1) ;
    Matrice res(1, 1) ;
    return res;
}



		   //-------------------------
		   //--   CAS R_CHEBU    -----
		   //-------------------------

Matrice _helmholtz_minus_pseudo_1d_mat_r_chebu_deux(int,int,double, double) ;

Matrice _helmholtz_minus_pseudo_1d_mat_r_chebu( int n, int l, double masse, 
						double alpha, double, int puis) {
  Matrice res(n-2, n-2) ; 
  res.set_etat_qcq() ;
  switch (puis) {
  case 2 :
    res = _helmholtz_minus_pseudo_1d_mat_r_chebu_deux (n, l,masse, alpha) ;
    break ;
  default :
    abort() ;
    exit(-1) ;
  }
  return res ;
}

    //Cas ou dzpuis = 2
Matrice _helmholtz_minus_pseudo_1d_mat_r_chebu_deux (int n, int l, double masse, 
						     double alpha) {
        

  Matrice res(n-2, n-2) ;
  res.set_etat_qcq() ;
  double* vect = new double[n] ;
  double* vect_bis = new double[n] ;
  double* vect_dd = new double[n] ;
  double* vect_d = new double[n] ;
  
  for (int i=0 ; i<n-2 ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 2*i+3 ;
    vect[i+1] = -4*i-4 ;
    vect[i+2] = 2*i+1 ;

    // Der sec.
    for (int j=0 ; j<n ; j++)
      vect_bis[j] = vect[j] ;
    
    d2sdx2_1d (n, &vect_bis, R_CHEBU) ;  // appel dans le cas unsurr
    mult2_xm1_1d_cheb (n, vect_bis, vect_dd) ; // multiplication par (x-1)^2
    
    // Der simple
    for (int j=0 ; j<n ; j++)
      vect_bis[j] = vect[j] ;

    dsdx_1d (n, &vect_bis, R_CHEBU) ;  // appel dans le cas unsurr
    mult_xm1_1d_cheb (n, vect_bis, vect_d) ; // multiplication par (x-1)
    
    // Mass term
    for (int j=0 ; j<n ; j++)
      vect_bis[j] = vect[j] ;
    sx2_1d (n, &vect_bis, R_CHEBU) ;
    
    for (int j=0 ; j<n-2 ; j++)
      res.set(j,i) = vect_dd[j] + 2*vect_d[j] - l*(l-1)*vect[j] - masse*masse/alpha/alpha*vect_bis[j] ; 
  }

  delete [] vect ;
  delete [] vect_bis ;
  delete [] vect_dd ;
  
  return res ;
} 


	
void Ope_helmholtz_minus_pseudo_1d::do_ope_mat() const {
  if (ope_mat != 0x0) 
    delete ope_mat ;

  // Routines de derivation
  static Matrice (*helmholtz_minus_pseudo_1d_mat[MAX_BASE])(int, int, double, 
						     double, double, int);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      helmholtz_minus_pseudo_1d_mat[i] = _helmholtz_minus_pseudo_1d_mat_pas_prevu ;
    }
    // Les routines existantes
    helmholtz_minus_pseudo_1d_mat[R_CHEBU >> TRA_R] = _helmholtz_minus_pseudo_1d_mat_r_chebu ;
  }
  ope_mat = new Matrice(helmholtz_minus_pseudo_1d_mat[base_r](nr, l_quant, masse, 
						       alpha, beta, dzpuis)) ;
}
