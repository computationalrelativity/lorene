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

char ope_sec_order_mat_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <math.h>
#include <stdlib.h>

#include "proto.h"
#include "ope_elementary.h"


                      //-----------------------------------
                     // Routine pour les cas non prevus -- 
                     //-----------------------------------

Matrice _sec_order_mat_pas_prevu(int, double, double, double, double) {
  cout << "Sec_order : base not implemented..." << endl ;
  abort() ;
  exit(-1) ;
  Matrice res(1, 1) ;
  return res;
}



                    //-------------------------
		    //--   CAS R_CHEB   -----
		    //------------------------

Matrice _sec_order_mat_r_cheb (int n, double alpha, 
			       double a, double b, double c) {
  
  double* vect = new double[n] ;
  
  Matrice dd (n,n) ;
  dd.set_etat_qcq() ;
  Matrice df (n,n) ;
  df.set_etat_qcq() ;
  Matrice ff (n,n) ;
  ff.set_etat_qcq() ;


  // Calcul terme en d^2...
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 1 ;
  
    d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin

    for (int j=0 ; j<n ; j++)
      dd.set(j,i) = vect[j] ;
  }

  // Calcul terme en d...
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 1 ;
  
    sxdsdx_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
    
    for (int j=0 ; j<n ; j++)
      df.set(j,i) = vect[j] ;
  }

  // Calcul terme Id
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      ff.set(j,i) = 0 ;
    ff.set(i,i) = 1 ;
  }

  delete [] vect ;

  return a*dd/alpha/alpha+b*df/alpha+c*ff ;
}

void Ope_sec_order::do_ope_mat() const {
  if (ope_mat != 0x0) 
    delete ope_mat ;

  // Routines de derivation
  static Matrice (*sec_order_mat[MAX_BASE])(int, double, 
					    double, double, double);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      sec_order_mat[i] = _sec_order_mat_pas_prevu ;
    }
    // Les routines existantes
    sec_order_mat[R_CHEB >> TRA_R] = _sec_order_mat_r_cheb ;
  }
  ope_mat = new Matrice(sec_order_mat[base_r](nr, alpha, 
					      a_param, b_param, c_param)) ;
}
