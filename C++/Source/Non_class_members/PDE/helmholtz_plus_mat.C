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


char helmholtz_plus_mat_C[] = "$$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/12/11 14:48:49  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * 
 * $Header$
 *
 */

#include "matrice.h"
#include "type_parite.h"
#include "proto.h"

                     //-----------------------------------
                     // Routine pour les cas non prevus -- 
                     //-----------------------------------

Matrice _helmholtz_plus_mat_pas_prevu(int, double, double, double) {
  cout << "Helmholtz plus : base not implemented..." << endl ;
  abort() ;
  exit(-1) ;
  Matrice res(1, 1) ;
  return res;
}


                    //-------------------------
		    //--   CAS R_CHEB   -----
		    //------------------------

Matrice _helmholtz_plus_mat_r_cheb (int n, double alpha, double beta, 
				     double masse) {

  

  assert (masse > 0) ;

  const int nmax = 10 ;// Nombre de Matrices stockees
  static Matrice* tab[nmax] ;  // les matrices calculees
  static int nb_dejafait = 0 ; // nbre de matrices calculees
  static double masse_dejafait[nmax] ;
    
  static double vieux_alpha = 0 ;
  static double vieux_beta = 0 ;
  static int vieux_n = 0;
   
  // Si on a change alpha ou n, on detruit tout :
  if ((vieux_alpha != alpha) || (vieux_n != n) || (vieux_beta != beta)) {
    for (int i=0 ; i<nb_dejafait ; i++) {
      masse_dejafait[i] = 0 ;
      delete tab[i] ;
    }
    
    nb_dejafait = 0 ;
    vieux_alpha = alpha ;
    vieux_beta = beta ;
    vieux_n = n ;
  }
   
  int indice = -1 ;
   
  // On determine si la matrice a deja ete calculee :
  for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if (masse_dejafait[conte] == masse)
      indice = conte ;
    
  // Calcul a faire : 
  if (indice  == -1) {
    if (nb_dejafait >= nmax) {
      cout << "_helmholtz_plus_mat_r_cheb : trop de matrices" << endl ;
      abort() ;
      exit (-1) ;
    }
       
    double echelle = beta / alpha ;

    masse_dejafait[nb_dejafait] = masse ;  
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xd(n, n) ;
    xd.set_etat_qcq() ;
  
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
    
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd+2*xd ; 

    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
  }
    
  // Cas ou le calcul a deja ete effectue :
  else
    return *tab[indice] ;  
}
	
                //--------------------------
		//- La routine a appeler  ---
	        //----------------------------

Matrice helmholtz_plus_mat(int n, double alpha, double beta, double masse, 
			    int base_r)
{
  
  // Routines de derivation
  static Matrice (*helmholtz_plus_mat[MAX_BASE])(int, double, double, double);
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      helmholtz_plus_mat[i] = _helmholtz_plus_mat_pas_prevu ;
    }
    // Les routines existantes
    helmholtz_plus_mat[R_CHEB >> TRA_R] = _helmholtz_plus_mat_r_cheb ;
  }
  
  Matrice res(helmholtz_plus_mat[base_r](n, alpha, beta, masse)) ;
  return res ;
}

