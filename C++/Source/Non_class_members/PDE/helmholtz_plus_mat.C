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


 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2016/12/05 16:18:09  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:53:28  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:16:08  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2007/05/06 10:48:12  p_grandclement
 * Modification of a few operators for the vorton project
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
#include <cstdlib>

#include "matrice.h"
#include "type_parite.h"
#include "proto.h"

                     //-----------------------------------
                     // Routine pour les cas non prevus -- 
                     //-----------------------------------

namespace Lorene {
Matrice _helmholtz_plus_mat_pas_prevu(int, int, double, double, double) {
  cout << "Helmholtz plus : base not implemented..." << endl ;
  abort() ;
  exit(-1) ;
  Matrice res(1, 1) ;
  return res;
}



                    //-------------------------
                    //--   CAS R_CHEBP    -----
		   //--------------------------
		    

Matrice _helmholtz_plus_mat_r_chebp (int n, int lq, double alpha, double, 
				     double masse) {
  assert (masse > 0) ;
 
  Matrice dd(n, n) ;
  dd.set_etat_qcq() ;
  Matrice xd(n, n) ;
  xd.set_etat_qcq() ;
  Matrice xx(n, n) ;
  xx.set_etat_qcq() ;
  Matrice sx2(n, n) ;
  sx2.set_etat_qcq() ;

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
    sxdsdx_1d (n, &vect, R_CHEBP) ;
    for (int j=0 ; j<n ; j++)
      xd.set(j, i) = vect[j] ;  
  }
   
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 1 ;
    sx2_1d (n, &vect, R_CHEBP) ;
    for (int j=0 ; j<n ; j++)
      sx2.set(j, i) = vect[j] ;  
  }

  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      xx.set(i,j) = 0 ;
    xx.set(i, i) = 1 ;  
  }
  
  delete [] vect ;
    
  Matrice res(n, n) ;
  res = dd+2*xd-lq*(lq+1)*sx2+masse*masse*alpha*alpha*xx ;
  
  return res ;
}
    

                    //-------------------------
		    //--   CAS R_CHEB   -----
		    //------------------------

Matrice _helmholtz_plus_mat_r_cheb (int n, int lq, double alpha, double beta, 
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
    vect[i] += masse*masse*alpha*alpha ;
    for (int j=0 ; j<n ; j++)
      dd.set(j, i) = vect[j]*echelle*echelle ;
  }
  
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 1 ;
    d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
    vect[i] += masse*masse*alpha*alpha ;
    multx_1d (n, &vect, R_CHEB) ;
    for (int j=0 ; j<n ; j++)
      dd.set(j, i) += 2*echelle*vect[j] ;
  }
  
  for (int i=0 ; i<n ; i++) {
    for (int j=0 ; j<n ; j++)
      vect[j] = 0 ;
    vect[i] = 1 ;
    d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
    vect[i] += masse*masse*alpha*alpha ;
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

Matrice helmholtz_plus_mat(int n, int lq, double alpha, double beta, double masse, 
			    int base_r)
{
  
  // Routines de derivation
  static Matrice (*helmholtz_plus_mat[MAX_BASE])(int, int, double, double, double);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      helmholtz_plus_mat[i] = _helmholtz_plus_mat_pas_prevu ;
    }
    // Les routines existantes
    helmholtz_plus_mat[R_CHEB >> TRA_R] = _helmholtz_plus_mat_r_cheb ;
    helmholtz_plus_mat[R_CHEBP >> TRA_R] = _helmholtz_plus_mat_r_chebp ;
  }
  
  Matrice res(helmholtz_plus_mat[base_r](n, lq, alpha, beta, masse)) ;
  return res ;
}

}
