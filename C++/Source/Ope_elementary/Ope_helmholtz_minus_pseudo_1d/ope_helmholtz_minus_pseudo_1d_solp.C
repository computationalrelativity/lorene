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

 

/*
 * $Id$
 * $Header$
 *
 */
#include <cmath>
#include <cstdlib>

#include "proto.h"
#include "ope_elementary.h"
//--------------------------------------------------
// Version Tbl --> Tbl a 1D pour la source 
//--------------------------------------------------


namespace Lorene {
Tbl _cl_helmholtz_minus_pseudo_1d_pas_prevu (const Tbl & source, int) {
     cout << "Combinaison lineaire pas prevue..." << endl ;
    abort() ;
    exit(-1) ;
    return source;
}




		//-------------------
	       //--  R_CHEBU   -----
	      //-------------------
Tbl _cl_helmholtz_minus_pseudo_1d_r_chebu_deux(const Tbl&) ;

Tbl _cl_helmholtz_minus_pseudo_1d_r_chebu (const Tbl &source, int puis) {

    int n=source.get_dim(0) ;
    Tbl res(n) ;
    res.set_etat_qcq() ;
    
    switch(puis) {
	case 2 :
	    res = _cl_helmholtz_minus_pseudo_1d_r_chebu_deux(source) ;
	    break ;
	
	default :
	    abort() ;
	    exit(-1) ;    
    }
   return res ;
}

// Cas dzpuis = 2 ;
Tbl _cl_helmholtz_minus_pseudo_1d_r_chebu_deux (const Tbl &source) {

  Tbl barre(source) ;
  int n = source.get_dim(0) ;
    
  int dirac = 1 ;
  for (int i=0 ; i<n-2 ; i++) {
    barre.set(i) = ((1+dirac)*source(i)-source(i+2)) ;
    if (i==0) dirac = 0 ;
  }
  
  Tbl tilde(barre) ;
  for (int i=0 ; i<n-4 ; i++)
    tilde.set(i) = (barre(i)-barre(i+2)) ;
  
  Tbl bis(tilde) ;
  for (int i=0 ; i<n-4 ; i++)
    bis.set(i) = (tilde(i)+tilde(i+1)) ;
  
  Tbl res(bis) ;
  for (int i=0 ; i<n-4 ; i++)
    res.set(i) = (bis(i)-bis(i+1)) ;
  
  return res ;
}


		//----------------------------
	       //- Routine a appeler        ---
	      //------------------------------

Tbl cl_helmholtz_minus_pseudo_1d (const Tbl &source, int puis, int base_r) {
		// Routines de derivation
    static Tbl (*cl_helmholtz_minus_pseudo_1d[MAX_BASE])(const Tbl &, int) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    cl_helmholtz_minus_pseudo_1d[i] = _cl_helmholtz_minus_pseudo_1d_pas_prevu ;
	}
		// Les routines existantes
	cl_helmholtz_minus_pseudo_1d[R_CHEBU >> TRA_R] = _cl_helmholtz_minus_pseudo_1d_r_chebu ;

    }
    
    Tbl res(cl_helmholtz_minus_pseudo_1d[base_r](source, puis)) ;
    return res ;
}


		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _solp_helmholtz_minus_pseudo_1d_pas_prevu (const Matrice &, const Matrice &,
				double,  double, const Tbl &, int) {
    cout << " Solution homogene pas prevue ..... : "<< endl ;
    abort() ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	       	//-------------------
	       //--  R_CHEBU   -----
	      //-------------------
Tbl _solp_helmholtz_minus_pseudo_1d_r_chebu_deux (const Matrice&, const Matrice&, 
				   const Tbl&) ;

Tbl _solp_helmholtz_minus_pseudo_1d_r_chebu (const Matrice &lap, const Matrice &nondege, 
			     double, double, 
			     const Tbl &source, int puis) {
    int n = lap.get_dim(0) ;
    Tbl res(n+2) ;
    res.set_etat_qcq() ;
    
    switch (puis) {
	case 2 :
	    res = _solp_helmholtz_minus_pseudo_1d_r_chebu_deux 
	      (lap, nondege, source) ;
	    break ;
	default :
	    abort() ;
	    exit(-1) ;
	}
return res ;
}
	
// Cas dzpuis = 2 ;
Tbl _solp_helmholtz_minus_pseudo_1d_r_chebu_deux (const Matrice &lap, const Matrice &nondege, 
				   const Tbl &source) {
  
  int n = lap.get_dim(0)+2 ;	  
  int dege = n-nondege.get_dim(0) ;
  assert (dege == 3) ;

  Tbl source_cl (cl_helmholtz_minus_pseudo_1d(source, 2, R_CHEBU)) ;
  
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


Tbl Ope_helmholtz_minus_pseudo_1d::get_solp (const Tbl& so) const {

  if (non_dege == 0x0)
    do_non_dege() ;

  // Routines de derivation
  static Tbl (*solp_helmholtz_minus_pseudo_1d[MAX_BASE]) (const Matrice&, const Matrice&,
					   double, double,const Tbl&, int) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solp_helmholtz_minus_pseudo_1d[i] = _solp_helmholtz_minus_pseudo_1d_pas_prevu ;
    }
    // Les routines existantes
    solp_helmholtz_minus_pseudo_1d[R_CHEBU >> TRA_R] = _solp_helmholtz_minus_pseudo_1d_r_chebu ;
  }
  
  Tbl res(solp_helmholtz_minus_pseudo_1d[base_r] (*ope_cl, *non_dege, 
						  alpha, beta, so, dzpuis)) ;

  Tbl valeurs (val_solp (res, alpha, base_r)) ;
  valeurs *= sqrt(double(2)) ;
  
  sp_plus = valeurs(0) ;
  sp_minus = valeurs(1) ;
  dsp_plus = valeurs(2) ;
  dsp_minus = valeurs(3) ;

 
  return res ;
}
}
