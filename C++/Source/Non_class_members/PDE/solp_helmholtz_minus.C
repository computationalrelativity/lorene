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


char solp_helmholtz_minus_C[] = "$Header$" ;

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

//fichiers includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrice.h"
#include "tbl.h"
#include "type_parite.h"
#include "indent.h"
#include "proto.h"

		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _solp_helmholtz_minus_pas_prevu (const Matrice &, const Matrice &, 
				     const Tbl &, double, double) {
    cout << " Solution homogene pas prevue ..... : "<< endl ;
    abort() ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	
		//-------------------
	       //--  R_CHEBU   ------
	      //-------------------


Tbl _solp_helmholtz_minus_r_chebu (const Matrice &lap, const Matrice &nondege, 
				   const Tbl &source, double, double) {
  
  int n = lap.get_dim(0)+2 ;	  
  int dege = n-nondege.get_dim(0) ;
  assert (dege==3) ;
  
  Tbl source_cl (cl_helmholtz_minus(source, R_CHEBU)) ;
  
  Tbl so(n-dege) ;
  so.set_etat_qcq() ;
  for (int i=0 ; i<n-dege ; i++)
    so.set(i) = source_cl(i);
  
  Tbl sol (nondege.inverse(so)) ;
    
  Tbl res(n) ;
  res.annule_hard() ;
  for (int i=1 ; i<n-2 ; i++) {
    res.set(i) += sol(i-1)*(2*i+3) ;
    res.set(i+1) += -sol(i-1)*(4*i+4) ;
    res.set(i+2) += sol(i-1)*(2*i+1) ;
  }
  return res ;
}


          	//-------------------
	       //--  R_CHEB   -----
	      //-------------------
Tbl _solp_helmholtz_minus_r_cheb (const Matrice &lap, const Matrice &nondege, 
				const Tbl &source, double alpha, double beta) {
  
  int n = lap.get_dim(0) ;	  
  int dege = n-nondege.get_dim(0) ;
  assert (dege ==2) ;
  
  Tbl source_aux(source*alpha*alpha) ;
  Tbl xso(source_aux) ;
  Tbl xxso(source_aux) ;
  multx_1d(n, &xso.t, R_CHEB) ;
  multx_1d(n, &xxso.t, R_CHEB) ;
  multx_1d(n, &xxso.t, R_CHEB) ;
  source_aux = beta*beta/alpha/alpha*source_aux+2*beta/alpha*xso+xxso ;
  
  source_aux = cl_helmholtz_minus (source_aux, R_CHEB) ;
  
  Tbl so(n-dege) ;
  so.set_etat_qcq() ;
  for (int i=0 ; i<n-dege ; i++)
    so.set(i) = source_aux(i) ;
 
  Tbl auxi(nondege.inverse(so)) ;

  Tbl res(n) ;
  res.set_etat_qcq() ;
  for (int i=dege ; i<n ; i++)
    res.set(i) = auxi(i-dege) ;
  
  for (int i=0 ; i<dege ; i++)
    res.set(i) = 0 ;
  return res ;
}


	      	//-------------------
	       //--  Fonction   ----
	      //-------------------
	      
	      
Tbl solp_helmholtz_minus (const Matrice &lap, const Matrice &nondege, 
			  const Tbl &source, double alpha, double beta, 
			  int base_r) {

  // Routines de derivation
  static Tbl (*solp_helmholtz_minus[MAX_BASE]) (const Matrice&, const Matrice&,
						const Tbl&, double, double) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solp_helmholtz_minus[i] = _solp_helmholtz_minus_pas_prevu ;
    }
    // Les routines existantes
    solp_helmholtz_minus[R_CHEB >> TRA_R] = _solp_helmholtz_minus_r_cheb ;
    solp_helmholtz_minus[R_CHEBU >> TRA_R] = _solp_helmholtz_minus_r_chebu ;
  }
  
  Tbl res(solp_helmholtz_minus[base_r] (lap, nondege, source, alpha, beta)) ;
  return res ;
}
