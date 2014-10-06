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

char ope_vorton_solp_C[] = "$Header$" ;

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
Tbl _cl_vorton_pas_prevu (const Tbl &so, int) {

  cout << "Linear combination for vorton not implemented..." << endl ;
  abort() ;
  exit(-1) ;
  return so;
}

               //-------------------
	       //--  R_CHEB  -------
	      //--------------------
Tbl _cl_vorton_r_cheb (const Tbl& source, int) {
  
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

	        //-------------------
	       //--  R_CHEBU   -----
	      //-------------------

Tbl _cl_vorton_r_chebu_trois (const Tbl &source) {
    int n = source.get_dim(0) ;
    Tbl barre(source) ;
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	     barre.set(i) = ((1+dirac)*source(i)-source(i+2)) ;
	if (i==0) dirac = 0 ;
	}
    
    Tbl tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	    tilde.set(i) = (barre(i)-barre(i+2)) ;
    
    Tbl res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	    res.set(i) = (tilde(i)+tilde(i+1)) ;
    return res ;
}

Tbl _cl_vorton_r_chebu (const Tbl& source, int puis) {
  int n = source.get_dim(0) ;
  Tbl res (n) ;
  res.set_etat_qcq() ;
   
   switch (puis) {
	case 3 :
	    res = _cl_vorton_r_chebu_trois(source) ;
	    break ;
	default :
	    abort() ;
	    exit(-1) ;
      }
return res ;
}

		//----------------------------
	       //- Routine a appeler        ---
	      //------------------------------

Tbl cl_vorton (const Tbl &source, int puis, int base_r) {
    
  // Routines de derivation
  static Tbl (*cl_vorton[MAX_BASE])(const Tbl &, int) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      cl_vorton[i] = _cl_vorton_pas_prevu ;
    }
    // Les routines existantes
    cl_vorton[R_CHEB >> TRA_R] = _cl_vorton_r_cheb ;
    cl_vorton[R_CHEBU >> TRA_R] = _cl_vorton_r_chebu ;
  }
    
    Tbl res(cl_vorton[base_r](source, puis)) ;
    return res ;
}


                       //*******************************
                       //  CALCUL SP proprement parler
                       //*******************************

		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _solp_vorton_pas_prevu (const Matrice &, const Matrice &, 
				     const Tbl &, double, double, int) {
    cout << " Solution particuliere pas prevue in sec_order..... : "<< endl ;
    abort() ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
                //-------------------
	       //--  R_CHEBU   -----
	      //-------------------

Tbl _solp_vorton_r_chebu_trois (const Matrice &lap, const Matrice &nondege, 
			       const Tbl &source, double alpha) {
  
  int n = lap.get_dim(0) ;	  
  int dege = n-nondege.get_dim(0) ;;
  assert ((dege==2) || (dege==1)) ;

  Tbl source_aux (alpha*cl_vorton (source, 3, R_CHEBU)) ;

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
  
  if (dege==2) {
      double somme = 0 ;
      for (int i=0 ; i<n ; i++)
	  somme -= res(i) ;
  res.set(0) = somme ;
  }

  return res ;
}


Tbl _solp_vorton_r_chebu (const Matrice &lap, const Matrice &nondege, 
			       const Tbl &source,double alpha, double, int puis) {
  int n = source.get_dim(0) ;
  Tbl res (n) ;
  res.set_etat_qcq() ;
   
   switch (puis) {
	case 3 :
	    res = _solp_vorton_r_chebu_trois(lap, nondege, source, alpha) ;
	    break ;
	default :
	    abort() ;
	    exit(-1) ;
      }
return res ;
}

                //-------------------
	       //--  R_CHEB   -----
	      //-------------------

Tbl _solp_vorton_r_cheb (const Matrice &lap, const Matrice &nondege, 
			       const Tbl &source,double alpha, double beta, int dz) {
  
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
  source_aux = cl_vorton (source_aux, dz, R_CHEB) ;

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



Tbl Ope_vorton::get_solp (const Tbl& so) const {

  if (non_dege == 0x0)
    do_non_dege() ;

  // Routines de derivation
  static Tbl (*solp_vorton[MAX_BASE]) (const Matrice&, const Matrice&,
					     const Tbl&, double, double, int) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solp_vorton[i] = _solp_vorton_pas_prevu ;
    }
    // Les routines existantes
    solp_vorton[R_CHEB >> TRA_R] = _solp_vorton_r_cheb ; 
    solp_vorton[R_CHEBU >> TRA_R] = _solp_vorton_r_chebu ;
  }
  
  Tbl res(solp_vorton[base_r] (*ope_mat, *non_dege, so, alpha, beta, dzpuis)) ;
  Tbl valeurs (val_solp (res, alpha, base_r)) ;

  sp_plus = valeurs(0) ;
  sp_minus = valeurs(1) ;
  dsp_plus = valeurs(2) ;
  dsp_minus = valeurs(3) ;

  return res ;
}
