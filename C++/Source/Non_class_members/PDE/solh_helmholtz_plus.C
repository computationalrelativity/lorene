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


char solh_helmholtz_plusC[] = "$Header $" ;


/*
 * $Id$
 * $Log$
 * Revision 1.6  2014/10/13 08:53:31  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:16:10  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2008/02/18 13:53:43  j_novak
 * Removal of special indentation instructions.
 *
 * Revision 1.3  2007/05/06 10:48:12  p_grandclement
 * Modification of a few operators for the vorton project
 *
 * Revision 1.2  2004/01/15 09:15:37  p_grandclement
 * Modification and addition of the Helmholtz operators
 *
 * Revision 1.1  2003/12/11 14:48:49  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * 
 * $Header$
 *
 */

//fichiers includes
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>

#include "proto.h"
#include "matrice.h"
#include "type_parite.h"

                //------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
namespace Lorene {
Tbl _solh_helmholtz_plus_pas_prevu (int, int, double, double, double) {

  cout << "Homogeneous solution not implemented in hemlholtz_plus : "<< endl ;
  abort() ;
  exit(-1) ;
  Tbl res(1) ;
  return res;
}
	

	
		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _solh_helmholtz_plus_r_cheb (int n, int lq, double alpha, double beta, 
				  double masse) {
  
  assert (masse > 0) ;
  
  Tbl res(2,n) ;
  res.set_etat_qcq() ;
  double* coloc = new double[n] ;
  
  int * deg = new int[3] ;
  deg[0] = 1 ; 
  deg[1] = 1 ;
  deg[2] = n ;
  
  // First SH
  for (int i=0 ; i<n ; i++){
    double air = alpha*(-cos(M_PI*i/(n-1))) + beta ;
    coloc[i] = -gsl_sf_bessel_jl(lq, masse*air) ;
  }
  
  cfrcheb(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(0,i) = coloc[i] ;
  
  // Second SH
  for (int i=0 ; i<n ; i++){
    double air = alpha*(-cos(M_PI*i/(n-1))) + beta ;
    coloc[i] = -gsl_sf_bessel_yl (lq, masse*air) ;
  }
  
  cfrcheb(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(1,i) = coloc[i] ;
  
  delete [] deg ;
  delete [] coloc ;
  return res ;
}

	
		//-------------------
	       //--  R_CHEBP   -----
	      //-------------------

Tbl _solh_helmholtz_plus_r_chebp (int n, int lq, double alpha, double, 
				  double masse) {
  
  assert (masse > 0) ;
  
  Tbl res(n) ;
  res.set_etat_qcq() ;
  double* coloc = new double[n] ;
  
  int * deg = new int[3] ;
  deg[0] = 1 ; 
  deg[1] = 1 ;
  deg[2] = n ;
  
  // SH
  for (int i=0 ; i<n ; i++){
    double air = alpha*sin(M_PI*i/2/(n-1)) ;
    coloc[i] = -gsl_sf_bessel_jl(lq, masse*air) ;
  }

  cfrchebp(deg, deg, coloc, deg, coloc) ;
  for (int i=0 ; i<n ;i++)
    res.set(i) = coloc[i] ;
  
  delete [] deg ;
  delete [] coloc ;
  return res ;
}


	      	//-------------------
	       //--  Fonction   ----
	      //-------------------
	      
	      
Tbl solh_helmholtz_plus (int n, int lq, double alpha, double beta, 
			  double masse, int base_r) {

  // Routines de derivation
  static Tbl (*solh_helmholtz_plus[MAX_BASE])(int, int, double, double, double) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solh_helmholtz_plus[i] = _solh_helmholtz_plus_pas_prevu ;
    }
    // Les routines existantes
    solh_helmholtz_plus[R_CHEB >> TRA_R] = _solh_helmholtz_plus_r_cheb ;
    solh_helmholtz_plus[R_CHEBP >> TRA_R] = _solh_helmholtz_plus_r_chebp ;
  }
    
  Tbl res(solh_helmholtz_plus[base_r](n, lq, alpha, beta, masse)) ;
  return res ;
}
}
