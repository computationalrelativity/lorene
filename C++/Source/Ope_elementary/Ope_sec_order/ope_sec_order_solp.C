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

char ope_sec_order_solp_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <cmath>
#include <cstdlib>

#include "proto.h"
#include "ope_elementary.h"


                //------------------------------------
		// Cl version Tbl -> Tbl            --
		//------------------------------------
namespace Lorene {
Tbl _cl_sec_order_pas_prevu (const Tbl &so) {

  cout << "Linear combination for Sec_order not implemented..." << endl ;
  abort() ;
  exit(-1) ;
  return so;
}

               //-------------------
	       //--  R_CHEB  -------
	      //--------------------
Tbl _cl_sec_order_r_cheb (const Tbl& source) {
  
  int n = source.get_dim(0) ;
  
  Tbl barre(source) ;
  int dirac = 1 ;
  for (int i=0 ; i<n-2 ; i++) {
    barre.set(i) = ((1+dirac)*source(i)-source(i+2))
      /(i+1) ;
    if (i==0) dirac = 0 ;
  }
  
  Tbl res(barre) ;
  for (int i=0 ; i<n-4 ; i++)
    res.set(i) = barre(i)-barre(i+2) ;

  return res ;
}
              

		//----------------------------
	       //- Routine a appeler        ---
	      //------------------------------

Tbl cl_sec_order (const Tbl &source, int base_r) {
    
  // Routines de derivation
  static Tbl (*cl_sec_order[MAX_BASE])(const Tbl &) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      cl_sec_order[i] = _cl_sec_order_pas_prevu ;
    }
    // Les routines existantes
    cl_sec_order[R_CHEB >> TRA_R] = _cl_sec_order_r_cheb ;
  }
    
    Tbl res(cl_sec_order[base_r](source)) ;
    return res ;
}


                       //*******************************
                       //  CALCUL SP proprement parler
                       //*******************************

		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _solp_sec_order_pas_prevu (const Matrice &, const Matrice &, 
				     const Tbl &) {
    cout << " Solution particuliere pas prevue in sec_order..... : "<< endl ;
    abort() ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}

                //-------------------
	       //--  R_CHEB   -----
	      //-------------------

Tbl _solp_sec_order_r_cheb (const Matrice &lap, const Matrice &nondege, 
			       const Tbl &source) {
  
  int n = lap.get_dim(0) ;	  
  int dege = n-nondege.get_dim(0) ;
  assert (dege ==2) ;
  
  Tbl source_aux (cl_sec_order (source, R_CHEB)) ;
  
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



Tbl Ope_sec_order::get_solp (const Tbl& so) const {

  if (non_dege == 0x0)
    do_non_dege() ;

  // Routines de derivation
  static Tbl (*solp_sec_order[MAX_BASE]) (const Matrice&, const Matrice&,
					     const Tbl&) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solp_sec_order[i] = _solp_sec_order_pas_prevu ;
    }
    // Les routines existantes
    solp_sec_order[R_CHEB >> TRA_R] = _solp_sec_order_r_cheb ;
  }
  
  Tbl res(solp_sec_order[base_r] (*ope_mat, *non_dege, so)) ;
  
  Tbl valeurs (val_solp (res, alpha, base_r)) ;
  sp_plus = valeurs(0) ;
  sp_minus = valeurs(1) ;
  dsp_plus = valeurs(2) / exp(alpha+beta) ;
  dsp_minus = valeurs(3) / exp(beta-alpha) ;
  
  return res ;
}
}
