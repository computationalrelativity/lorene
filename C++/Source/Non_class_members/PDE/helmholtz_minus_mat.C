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


char helmholtz_minus_mat_C[] = "$$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2004/08/24 09:14:44  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.3  2004/01/15 09:15:37  p_grandclement
 * Modification and addition of the Helmholtz operators
 *
 * Revision 1.2  2003/12/11 15:57:26  p_grandclement
 * include stdlib.h encore ...
 *
 * Revision 1.1  2003/12/11 14:48:49  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * 
 * $Header$
 *
 */
#include <stdlib.h>

#include "matrice.h"
#include "type_parite.h"
#include "proto.h"

                     //-----------------------------------
                     // Routine pour les cas non prevus -- 
                     //-----------------------------------

Matrice _helmholtz_minus_mat_pas_prevu(int, int, double, double, double) {
  cout << "Helmholtz minus : base not implemented..." << endl ;
  abort() ;
  exit(-1) ;
  Matrice res(1, 1) ;
  return res;
}

                    //-------------------------
		    //--   CAS R_CHEBU    -----
		    //------------------------

Matrice _helmholtz_minus_mat_r_chebu (int n, int lq, double alpha, 
				      double, double masse) {

  assert (masse > 0) ;

  Matrice res(n-2, n-2) ;
  res.set_etat_qcq() ;
  
  double* vect = new double[n] ;
  double* vect_bis = new double[n] ;
  double* vect_dd = new double[n] ;
  
  for (int i=0 ; i<n-2 ; i++) {
    
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 2*i+3 ;
    vect[i+1] = -4*i-4 ;
    vect[i+2] = 2*i+1 ;
    
    for (int j=0 ; j<n ; j++)
      vect_bis[j] = vect[j] ;
    
    d2sdx2_1d (n, &vect_bis, R_CHEBU) ;  // appel dans le cas unsurr
    mult2_xm1_1d_cheb (n, vect_bis, vect_dd) ; // multiplication par (x-1)^2
   
    // Mass term
    for (int j=0 ; j<n ; j++)
      vect_bis[j] = vect[j] ;
    sx2_1d (n, &vect_bis, R_CHEBU) ;
    
    for (int j=0 ; j<n-2 ; j++)
      res.set(j,i) = vect_dd[j] - lq*(lq+1)*vect[j] 
	- masse*masse*vect_bis[j]/alpha/alpha ;
  }
  
  delete [] vect ;
  delete [] vect_bis ;
  delete [] vect_dd ;

  return res ;
}


                    //-------------------------
		    //--   CAS R_CHEB   -----
		    //------------------------

Matrice _helmholtz_minus_mat_r_cheb (int n, int lq, double alpha, double beta, 
				     double masse) {

  assert (masse > 0) ;
       
  double echelle = beta / alpha ;
    
  Matrice dd(n, n) ;
  dd.set_etat_qcq() ;
  Matrice xd(n, n) ;
  xd.set_etat_qcq() ;  
  Matrice xx(n, n) ;
  xx.set_etat_qcq() ;
  
  double* vect = new double[n] ;
  
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 1 ;
    d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
    vect[i] -= masse*masse*alpha*alpha ;
    for (int j=0 ; j<n ; j++)
      dd.set(j, i) = vect[j]*echelle*echelle ;
  }
  
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 1 ;
    d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
    vect[i] -= masse*masse*alpha*alpha ;
    multx_1d (n, &vect, R_CHEB) ;
    for (int j=0 ; j<n ; j++)
      dd.set(j, i) += 2*echelle*vect[j] ;
  }
  
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 1 ;
    d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
    vect[i] -= masse*masse*alpha*alpha ;
    multx_1d (n, &vect, R_CHEB) ;
    multx_1d (n, &vect, R_CHEB) ;
    for (int j=0 ; j<n ; j++)
      dd.set(j, i) += vect[j] ;
  }
  
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 1 ;
    sxdsdx_1d (n, &vect, R_CHEB) ;
    for (int j=0 ; j<n ; j++)
      xd.set(j, i) = vect[j]*echelle ;
  }
  
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 1 ;
    sxdsdx_1d (n, &vect, R_CHEB) ;
    multx_1d (n, &vect, R_CHEB) ;
    for (int j=0 ; j<n ; j++)
      xd.set(j, i) += vect[j] ;
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
  res = dd+2*xd - lq*(lq+1)*xx; 
  
  return res ;
}

	
                //--------------------------
		//- La routine a appeler  ---
	        //----------------------------

Matrice helmholtz_minus_mat(int n, int lq, 
			    double alpha, double beta, double masse, 
			    int base_r)
{
  
  // Routines de derivation
  static Matrice (*helmholtz_minus_mat[MAX_BASE])(int, int, 
						  double, double, double);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      helmholtz_minus_mat[i] = _helmholtz_minus_mat_pas_prevu ;
    }
    // Les routines existantes
    helmholtz_minus_mat[R_CHEB >> TRA_R] = _helmholtz_minus_mat_r_cheb ;
    helmholtz_minus_mat[R_CHEBU >> TRA_R] = _helmholtz_minus_mat_r_chebu ;
  }
  
  Matrice res(helmholtz_minus_mat[base_r](n, lq, alpha, beta, masse)) ;
  return res ;
}

