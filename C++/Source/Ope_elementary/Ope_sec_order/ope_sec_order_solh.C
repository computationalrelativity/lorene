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

char ope_sec_order_solh_C[] = "$Header$" ;

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
Tbl _solh_sec_order_pas_prevu (int, double, double,double,double,double,Tbl&) {

  cout << "Homogeneous solution not implemented in Sec_order : "<< endl ;
  abort() ;
  exit(-1) ;
  Tbl res(1) ;
  return res;
}

	
		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _solh_sec_order_r_cheb (int n, double alpha, double beta, 
			       double a, double b, double c, Tbl& val_lim) {

  // Stuff to compute the coefs...
  Tbl res(2,n) ;
  res.set_etat_qcq() ;
  double* coloc = new double[n] ;
    
  int * deg = new int[3] ;
  deg[0] = 1 ; 
  deg[1] = 1 ;
  deg[2] = n ;
  
  // Array on the radius 
  double* sigma = new double[n] ;
  for (int i=0 ; i<n ; i++)
    sigma[i] = alpha*(-cos(M_PI*i/(n-1))) + beta ;

  double sigma_minus = beta-alpha ;
  double sigma_plus  = beta+alpha ;

  // Determinant :
  double delta = b*b-4*a*c ;
  int signe ;
  if (delta > 0)
    signe = + 1 ;
  if (delta < 0)
    signe = -1 ;
  if (fabs(delta) < 1e-14)
    signe = 0 ;

  switch (signe) {
    
  case 1: {
    // Two real solutions
    double lambda_one = (-b + sqrt(delta)) / 2./a ;
    double lambda_two = (-b - sqrt(delta)) / 2./a ;

    // First SH
    for (int i=0 ; i<n ; i++)
      coloc[i] = exp(sigma[i]*lambda_one) ;
    cfrcheb(deg, deg, coloc, deg, coloc) ;
    for (int i=0 ; i<n ;i++)
      res.set(0,i) = coloc[i] ;

    // Second SH
    for (int i=0 ; i<n ; i++)
      coloc[i] = exp(sigma[i]*lambda_two) ;
    cfrcheb(deg, deg, coloc, deg, coloc) ;
    for (int i=0 ; i<n ;i++)
      res.set(1,i) = coloc[i] ;
    
    // Limit on the boundaries :
    val_lim.set(0,0) = exp(sigma_minus*lambda_one) ;
    val_lim.set(0,1) = lambda_one*exp(sigma_minus*lambda_one)
      /exp(sigma_minus) ;
    val_lim.set(0,2) = exp(sigma_plus*lambda_one) ;
    val_lim.set(0,3) = lambda_one*exp(sigma_plus*lambda_one)
      /exp(sigma_plus) ;

    val_lim.set(1,0) = exp(sigma_minus*lambda_two) ;
    val_lim.set(1,1) = lambda_two*exp(sigma_minus*lambda_two)/
      exp(sigma_minus);
    val_lim.set(1,2) = exp(sigma_plus*lambda_two) ;
    val_lim.set(1,3) = lambda_two*exp(sigma_plus*lambda_two)
      /exp(sigma_plus) ;
    val_lim /= sqrt(double(2)) ;
    break ;
  }
  case 0: {
    // Only one solution :
    double lambda = -b/2./a ;
    // First SH
    for (int i=0 ; i<n ; i++)
      coloc[i] = exp(sigma[i]*lambda) ;
    cfrcheb(deg, deg, coloc, deg, coloc) ;
    for (int i=0 ; i<n ;i++)
      res.set(0,i) = coloc[i] ;

    // Second SH
    for (int i=0 ; i<n ; i++)
      coloc[i] = sigma[i]*exp(sigma[i]*lambda) ;
    cfrcheb(deg, deg, coloc, deg, coloc) ;
    for (int i=0 ; i<n ;i++)
      res.set(1,i) = coloc[i] ; 

    // Limit on the boundaries :
    val_lim.set(0,0) = exp(sigma_minus*lambda) ;
    val_lim.set(0,1) = lambda*exp(sigma_minus*lambda)/exp(sigma_minus) ;
    val_lim.set(0,2) = exp(sigma_plus*lambda) ;
    val_lim.set(0,3) = lambda*exp(sigma_plus*lambda)/exp(sigma_plus) ;

    val_lim.set(1,0) = sigma_minus*exp(sigma_minus*lambda) ;
    val_lim.set(1,1) = exp(sigma_minus*lambda)* (lambda*sigma_minus+1)
      /exp(sigma_minus);
    val_lim.set(1,2) = sigma_plus*exp(sigma_plus*lambda) ;
    val_lim.set(1,3) =  exp(sigma_plus*lambda)* (lambda*sigma_plus+1)/ 
      exp(sigma_plus) ;
    val_lim /= sqrt(double(2)) ;
    break ;
  }
  case -1:{
    // Two imaginary solutions :
    double real_part = -b/2./a ;
    double imag_part = sqrt(-delta)/2./a ;

    // First SH
    for (int i=0 ; i<n ; i++)
      coloc[i] = exp(sigma[i]*real_part)*cos(imag_part*sigma[i]) ;
    cfrcheb(deg, deg, coloc, deg, coloc) ;
    for (int i=0 ; i<n ;i++)
      res.set(0,i) = coloc[i] ;
    
    // Second SH
    for (int i=0 ; i<n ; i++)
      coloc[i] = exp(sigma[i]*real_part)*sin(imag_part*sigma[i]) ;
    cfrcheb(deg, deg, coloc, deg, coloc) ;
    for (int i=0 ; i<n ;i++)
      res.set(1,i) = coloc[i] ;

    // Limit on the boundaries :
    val_lim.set(0,0) = exp(sigma_minus*real_part)*cos(imag_part*sigma_minus) ;
    val_lim.set(0,1) =  (real_part*cos(imag_part*sigma_minus) - 
			 imag_part*sin(imag_part*sigma_minus)) * exp(real_part*sigma_minus)
      /exp(sigma_minus);
    val_lim.set(0,2) = exp(sigma_plus*real_part)*cos(imag_part*sigma_plus) ;
    val_lim.set(0,3) = (real_part*cos(imag_part*sigma_plus) - 
			 imag_part*sin(imag_part*sigma_plus)) * exp(real_part*sigma_plus)
      /exp(sigma_plus) ;

    val_lim.set(1,0) = exp(sigma_minus*real_part)*sin(imag_part*sigma_minus) ;
    val_lim.set(1,1) = (real_part*sin(imag_part*sigma_minus) + 
			 imag_part*cos(imag_part*sigma_minus)) * exp(real_part*sigma_minus)
      /exp(sigma_minus);
    val_lim.set(1,2) = exp(sigma_plus*real_part)*sin(imag_part*sigma_plus);
    val_lim.set(1,3) = (real_part*sin(imag_part*sigma_plus) + 
			imag_part*cos(imag_part*sigma_plus)) * exp(real_part*sigma_plus)
      /exp(sigma_plus) ;
    val_lim /= sqrt(double(2)) ;
    break ;
  }
  default:
    cout << "What are you doing here ? Get out or I call the police !" << endl;
    abort() ;
    break ;
  }

  delete [] deg ;
  delete [] coloc ;
  delete [] sigma ;

  return res ;
}


Tbl Ope_sec_order::get_solh () const {

  // Routines de derivation
  static Tbl (*solh_sec_order[MAX_BASE]) (int, double, double,
					     double, double, double, Tbl&) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solh_sec_order[i] = _solh_sec_order_pas_prevu ;
    }
    // Les routines existantes
    solh_sec_order[R_CHEB >> TRA_R] = _solh_sec_order_r_cheb ;
  }
  
  Tbl val_lim (2,4) ;
  val_lim.set_etat_qcq() ;
  Tbl res(solh_sec_order[base_r](nr,alpha,beta,a_param,b_param,c_param, 
				    val_lim)) ;

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
