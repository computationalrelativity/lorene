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

char ope_poisson_pseudo_1d_solp_C[] = "$Header$" ;

/*
 * $Id$
 * $Header$
 *
 */
#include <math.h>
#include <stdlib.h>

#include "proto.h"
#include "ope_elementary.h"
//--------------------------------------------------
// Version Tbl --> Tbl a 1D pour la source 
//--------------------------------------------------


Tbl _cl_poisson_pseudo_1d_pas_prevu (const Tbl &source) {
     cout << "Combinaison lineaire pas prevue..." << endl ;
    abort() ;
    exit(-1) ;
    return source;
}



		//-------------------
	       //--  R_CHEB  -------
	      //--------------------

Tbl _cl_poisson_pseudo_1d_r_cheb (const Tbl &source) {
    Tbl barre(source) ;
    int n = source.get_dim(0) ;
    
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
	       //--  R_CHEBP   -----
	      //-------------------

Tbl _cl_poisson_pseudo_1d_r_chebp (const Tbl &source) {
    Tbl barre(source) ;
    int n = source.get_dim(0) ;
    
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	    barre.set(i) = (1+dirac)*source(i)-source(i+2) ;
	if (i==0) dirac = 0 ;
    }

    Tbl tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	    tilde.set(i) = barre(i)-barre(i+2) ;

    Tbl res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	    res.set(i) = tilde(i)-tilde(i+1) ;
	
   return res ;
}


		//-------------------
	       //--  R_CHEBI   -----
	      //-------------------

Tbl _cl_poisson_pseudo_1d_r_chebi (const Tbl &source) {
    Tbl barre(source) ;
    int n = source.get_dim(0) ;
    
    for (int i=0 ; i<n-2 ; i++)
	    barre.set(i) = source(i)-source(i+2) ;

    Tbl tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	    tilde.set(i) = barre(i)-barre(i+2) ;    

    Tbl res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	    res.set(i) = tilde(i)-tilde(i+1) ;
	
   return res ;
}




		//----------------------------
	       //- Routine a appeler        ---
	      //------------------------------

Tbl cl_poisson_pseudo_1d (const Tbl &source, int base_r) {
    
		// Routines de derivation
    static Tbl (*cl_poisson_pseudo_1d[MAX_BASE])(const Tbl &) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    cl_poisson_pseudo_1d[i] = _cl_poisson_pseudo_1d_pas_prevu ;
	}
		// Les routines existantes
	cl_poisson_pseudo_1d[R_CHEB >> TRA_R] = _cl_poisson_pseudo_1d_r_cheb ;
	cl_poisson_pseudo_1d[R_CHEBP >> TRA_R] = _cl_poisson_pseudo_1d_r_chebp ;
	cl_poisson_pseudo_1d[R_CHEBI >> TRA_R] = _cl_poisson_pseudo_1d_r_chebi ;
    }
    
    Tbl res(cl_poisson_pseudo_1d[base_r](source)) ;
    return res ;
}


		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _solp_poisson_pseudo_1d_pas_prevu (const Matrice &, const Matrice &,
				double,  double, const Tbl &) {
    cout << " Solution homogene pas prevue ..... : "<< endl ;
    abort() ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	
		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _solp_poisson_pseudo_1d_r_cheb (const Matrice &lap, 
				    const Matrice &nondege, 
				    double alpha, double beta, 
				    const Tbl &source) {
    
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
    source_aux = cl_poisson_pseudo_1d(source_aux, R_CHEB) ;
	
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
	       //--  R_CHEBP   -----
	      //-------------------

Tbl _solp_poisson_pseudo_1d_r_chebp (const Matrice &lap, 
				     const Matrice &nondege, 
				     double alpha, double , const Tbl &source) {
    
    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;
    assert ((dege==2) || (dege == 1)) ;
    Tbl source_aux(alpha*alpha*source) ;
    source_aux = cl_poisson_pseudo_1d(source_aux, R_CHEBP) ;
    
    Tbl so(n-dege) ;
    so.set_etat_qcq() ;
    for (int i=0 ; i<n-dege ; i++)
	so.set(i) = source_aux(i);
   
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
	    if (i%2 == 0)
		somme -= res(i) ;
	    else somme += res(i) ;
	res.set(0) = somme ;
	return res ;
    }
    else return res ;
}
	
	
	      	//-------------------
	       //--  R_CHEBI   -----
	      //-------------------
	
Tbl _solp_poisson_pseudo_1d_r_chebi (const Matrice &lap, 
				     const Matrice &nondege, 
				     double alpha, double, const Tbl &source) {
   

  int n = lap.get_dim(0) ;	  
  int dege = n-nondege.get_dim(0) ;
  assert ((dege==2) || (dege == 1)) ;
  Tbl source_aux(source*alpha*alpha) ;
  source_aux = cl_poisson_pseudo_1d(source_aux, R_CHEBI) ;
  
  Tbl so(n-dege) ;
  so.set_etat_qcq() ;
  for (int i=0 ; i<n-dege ; i++)
    so.set(i) = source_aux(i);
  
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
      if (i%2 == 0)
	somme -= (2*i+1)*res(i) ;
      else somme += (2*i+1)*res(i) ;
    res.set(0) = somme ;
    return res ;
  }
  else return res ;
}
	

Tbl Ope_poisson_pseudo_1d::get_solp (const Tbl& so) const {

  if (non_dege == 0x0)
    do_non_dege() ;

  // Routines de derivation
  static Tbl (*solp_poisson_pseudo_1d[MAX_BASE]) (const Matrice&, 
						  const Matrice&,
						  double, double,const Tbl&) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solp_poisson_pseudo_1d[i] = _solp_poisson_pseudo_1d_pas_prevu ;
    }
    // Les routines existantes
    solp_poisson_pseudo_1d[R_CHEB >> TRA_R] = _solp_poisson_pseudo_1d_r_cheb ;
    solp_poisson_pseudo_1d[R_CHEBP >> TRA_R] = _solp_poisson_pseudo_1d_r_chebp ;
    solp_poisson_pseudo_1d[R_CHEBI >> TRA_R] = _solp_poisson_pseudo_1d_r_chebi ;
  }
  
  Tbl res(solp_poisson_pseudo_1d[base_r] (*ope_mat, *non_dege, 
				   alpha, beta, so)) ;
  
  Tbl valeurs (val_solp (res, alpha, base_r)) ;
  valeurs *= sqrt(double(2)) ;
  sp_plus = valeurs(0) ;
  sp_minus = valeurs(1) ;
  dsp_plus = valeurs(2) ;
  dsp_minus = valeurs(3) ;
  
  return res ;
}
