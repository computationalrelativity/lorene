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

char ope_poisson_pseudo_1d_solh_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <math.h>
#include <stdlib.h>

#include "proto.h"
#include "ope_elementary.h"


		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _solh_poisson_pseudo_1d_pas_prevu (int, int,double, double, Tbl&) {

    cout << " Solution homogene pas prevue ..... : "<< endl ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	
		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _solh_poisson_pseudo_1d_r_cheb (int n, int l, double alpha, double beta, 
				    Tbl& val_lim) {
                

  double echelle = beta / alpha ;
  int expo_un = l ;
  int expo_deux = -l+1 ;

  val_lim.set(0,0) = pow(echelle-1, double(expo_un)) ;
  val_lim.set(0,1) = double(expo_un) * pow(echelle-1, double(expo_un-1))/alpha ;
  val_lim.set(0,2) = pow(echelle+1, double(expo_un)) ;
  val_lim.set(0,3) = double(expo_un) * pow(echelle+1, double(expo_un-1))/alpha ;
    
    
  val_lim.set(1,0) = pow(echelle-1, double(expo_deux)) ;
  val_lim.set(1,1) = double(expo_deux) * pow(echelle-1, double(expo_deux-1))/alpha ;
  val_lim.set(1,2) = pow(echelle+1, double(expo_deux)) ;
  val_lim.set(1,3) = double(expo_deux) * pow(echelle+1, double(expo_deux-1))/alpha ;
  
 
  Tbl res(2, n) ;
  res.set_etat_qcq() ;
  double* coloc = new double[n] ;
  
  int * deg = new int[3] ;
  deg[0] = 1 ; 
  deg[1] = 1 ;
  deg[2] = n ;
  
  //Construction de la premiere solution homogene :
  
  for (int i=0 ; i<n ; i++)
    coloc[i] = pow(echelle-cos(M_PI*i/(n-1)), double(expo_un)) ;
    
  cfrcheb(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(0, i) = coloc[i] ;
  
  // construction de la seconde solution homogene :
  for (int i=0 ; i<n ; i++)
    coloc[i] = pow(echelle-cos(M_PI*i/(n-1)), double(expo_deux)) ;
  
  cfrcheb(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(1, i) = coloc[i] ;	
  
  
  delete [] coloc ;
  delete [] deg ;
   
  return res ;
}
    
 	
		//-------------------
	       //--  R_CHEBP  ------
	      //-------------------

Tbl _solh_poisson_pseudo_1d_r_chebp (int n, int l, double alpha, 
			      double, Tbl& val_lim) {
  

  val_lim.set(0,0) = (l!=0) ? 1 : 0 ;
  val_lim.set(0,1) = (l!=1) ? 0 : 1 ;
  val_lim.set(0,2) = 1. ;
  val_lim.set(0,3) = double(l)/alpha ;


  assert (div(l, 2).rem ==0) ;
  
  Tbl res(n) ;
  res.set_etat_qcq() ;
  double* coloc = new double[n] ;
  
  int * deg = new int[3] ;
  deg[0] = 1 ; 
  deg[1] = 1 ;
  deg[2] = n ;
  
  for (int i=0 ; i<n ; i++)
    coloc[i] = pow(sin(M_PI*i/2/(n-1)), double(l)) ;

  cfrchebp(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(i) = coloc[i] ;
  
  delete [] coloc ;
  delete [] deg ;
   
  return res ;
}	
	
	      	//-------------------
	       //--  R_CHEBI   -----
	      //-------------------
	
Tbl _solh_poisson_pseudo_1d_r_chebi (int n, int l, 
				     double alpha, double, Tbl& val_lim) {

  val_lim.set(0,0) = 0 ;
  val_lim.set(0,1) = (l!=1) ? 0 : 1 ;
  val_lim.set(0,2) = 1. ;
  val_lim.set(0,3) = double(l)/alpha ;

 
    
  assert (div(l, 2).rem == 1)  ;

  Tbl res(n) ;
  res.set_etat_qcq() ;
  double* coloc = new double[n] ;
    
  int * deg = new int[3] ;
  deg[0] = 1 ; 
  deg[1] = 1 ;
  deg[2] = n ;
    
  for (int i=0 ; i<n ; i++)
    coloc[i] = pow(sin(M_PI*i/2/(n-1)), double(l)) ;
    
  cfrchebi(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(i) = coloc[i] ;
	
  delete [] coloc ;
  delete [] deg ;
   
  return res ;
}
	

Tbl Ope_poisson_pseudo_1d::get_solh () const {

  // Routines de derivation
  static Tbl (*solh_poisson_pseudo_1d[MAX_BASE]) (int, int, double, double, Tbl&) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solh_poisson_pseudo_1d[i] = _solh_poisson_pseudo_1d_pas_prevu ;
    }
    // Les routines existantes
    solh_poisson_pseudo_1d[R_CHEB >> TRA_R] = _solh_poisson_pseudo_1d_r_cheb ;
    solh_poisson_pseudo_1d[R_CHEBP >> TRA_R] = _solh_poisson_pseudo_1d_r_chebp ;
    solh_poisson_pseudo_1d[R_CHEBI >> TRA_R] = _solh_poisson_pseudo_1d_r_chebi ;
  }
  
  Tbl val_lim (2,4) ;
  val_lim.set_etat_qcq() ;
  Tbl res(solh_poisson_pseudo_1d[base_r](nr,l_quant, alpha, beta, val_lim)) ;

  s_one_minus  = val_lim(0,0) ;
  ds_one_minus = val_lim(0,1) ;
  s_one_plus   = val_lim(0,2) ;
  ds_one_plus  = val_lim(0,3) ;

  s_two_minus  = val_lim(1,0) ;
  ds_two_minus = val_lim(1,1) ;
  s_two_plus   = val_lim(1,2) ;
  ds_two_plus  = val_lim(1,3) ;

  return res ;
}
