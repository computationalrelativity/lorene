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

char ope_poisson_2d_solp_C[] = "$Header$" ;

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


Tbl _cl_poisson_2d_pas_prevu (const Tbl &source, int puis) {
     cout << "Combinaison lineaire pas prevue..." << endl ;
    cout << "source : " << &source << endl ;
    cout << "dzpuis : " << puis << endl ;
    abort() ;
    exit(-1) ;
    return source;
}



		//-------------------
	       //--  R_CHEB  -------
	      //--------------------

Tbl _cl_poisson_2d_r_cheb (const Tbl &source, int) {
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

Tbl _cl_poisson_2d_r_chebp (const Tbl &source, int) {
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

Tbl _cl_poisson_2d_r_chebi (const Tbl &source, int) {
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


		//-------------------
	       //--  R_CHEBU   -----
	      //-------------------
Tbl _cl_poisson_2d_r_chebu_quatre(const Tbl&) ;
Tbl _cl_poisson_2d_r_chebu_trois(const Tbl&) ;
Tbl _cl_poisson_2d_r_chebu_deux(const Tbl&) ;

Tbl _cl_poisson_2d_r_chebu (const Tbl &source, int puis) {
    
    int n=source.get_dim(0) ;
    Tbl res(n) ;
    res.set_etat_qcq() ;
    
    switch(puis) {
	case 4 :
	    res = _cl_poisson_2d_r_chebu_quatre(source) ;
	    break ;
	case 3 :
	    res = _cl_poisson_2d_r_chebu_trois (source) ;
	    break ;
	case 2 :
	    res = _cl_poisson_2d_r_chebu_deux(source) ;
	    break ;
	
	default :
	    abort() ;
	    exit(-1) ;    
    }
   return res ;
}

// Cas dzpuis = 4 ;
Tbl _cl_poisson_2d_r_chebu_quatre (const Tbl &source) {
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
	    
    Tbl prime(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	    prime.set(i) = (tilde(i)-tilde(i+1)) ;
    
     Tbl res(prime) ;
	for (int i=0 ; i<n-4 ; i++)
		res.set(i) = (prime(i)-prime(i+2)) ;
 
   return res ;
}
// cas dzpuis = 3
Tbl _cl_poisson_2d_r_chebu_trois (const Tbl &source) {
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
	    
    Tbl res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	    res.set(i) = (tilde(i)+tilde(i+1)) ;
 
   return res ;
}

// Cas dzpuis = 2 ;
Tbl _cl_poisson_2d_r_chebu_deux (const Tbl &source) {
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
	    
    Tbl res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	    res.set(i) = (tilde(i)+tilde(i+1)) ;
   return res ;
}


		//----------------------------
	       //- Routine a appeler        ---
	      //------------------------------

Tbl cl_poisson_2d (const Tbl &source, int puis, int base_r) {
    
		// Routines de derivation
    static Tbl (*cl_poisson_2d[MAX_BASE])(const Tbl &, int) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    cl_poisson_2d[i] = _cl_poisson_2d_pas_prevu ;
	}
		// Les routines existantes
	cl_poisson_2d[R_CHEB >> TRA_R] = _cl_poisson_2d_r_cheb ;
	cl_poisson_2d[R_CHEBU >> TRA_R] = _cl_poisson_2d_r_chebu ;
	cl_poisson_2d[R_CHEBP >> TRA_R] = _cl_poisson_2d_r_chebp ;
	cl_poisson_2d[R_CHEBI >> TRA_R] = _cl_poisson_2d_r_chebi ;
    }
    
    Tbl res(cl_poisson_2d[base_r](source, puis)) ;
    return res ;
}


		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _solp_poisson_2d_pas_prevu (const Matrice &, const Matrice &,
				double,  double, const Tbl &, int) {
    cout << " Solution homogene pas prevue ..... : "<< endl ;
    abort() ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	
		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _solp_poisson_2d_r_cheb (const Matrice &lap, const Matrice &nondege, 
			     double alpha, double beta, 
			     const Tbl &source, int) {
    
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
    source_aux = cl_poisson_2d(source_aux, 0, R_CHEB) ;
	
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

Tbl _solp_poisson_2d_r_chebp (const Matrice &lap, const Matrice &nondege, 
			      double alpha, double , const Tbl &source, int) {
    
    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;
    assert ((dege==2) || (dege == 1)) ;
    Tbl source_aux(alpha*alpha*source) ;
    source_aux = cl_poisson_2d(source_aux, 0, R_CHEBP) ;
    
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
	
Tbl _solp_poisson_2d_r_chebi (const Matrice &lap, const Matrice &nondege, 
			      double alpha, double, const Tbl &source, int) {
   

    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;
    assert ((dege == 2) || (dege ==1)) ;
    Tbl source_aux(source*alpha*alpha) ;
    source_aux = cl_poisson_2d(source_aux, 0, R_CHEBI) ;
    
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
	
	
	
	       	//-------------------
	       //--  R_CHEBU   -----
	      //-------------------
Tbl _solp_poisson_2d_r_chebu_quatre (const Matrice&, const Matrice&, 
				     double, const Tbl&) ;
Tbl _solp_poisson_2d_r_chebu_trois (const Matrice&, const Matrice&, 
				     double, const Tbl&) ;
Tbl _solp_poisson_2d_r_chebu_deux (const Matrice&, const Matrice&, 
				   const Tbl&) ;

Tbl _solp_poisson_2d_r_chebu (const Matrice &lap, const Matrice &nondege, 
			     double alpha, double, 
			     const Tbl &source, int puis) {
    int n = lap.get_dim(0) ;
    Tbl res(n) ;
    res.set_etat_qcq() ;
    
    switch (puis) {

	case 4 :
	    res = _solp_poisson_2d_r_chebu_quatre 
	      (lap, nondege, alpha, source) ;
	    break ;
	case 3 :
	    res = _solp_poisson_2d_r_chebu_trois 
	      (lap, nondege, alpha, source) ;
	    break ;
	case 2 :
	    res = _solp_poisson_2d_r_chebu_deux 
	      (lap, nondege, source) ;
	    break ;
	default :
	    abort() ;
	    exit(-1) ;
	}
return res ;
}


// Cas dzpuis = 4 ;
Tbl _solp_poisson_2d_r_chebu_quatre (const Matrice &lap, const Matrice &nondege, 
				     double alpha, const Tbl &source) {
    
    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;
    assert ((dege==3) || (dege ==2)) ;
    Tbl source_aux(source*alpha*alpha) ;
    source_aux = cl_poisson_2d(source_aux, 4, R_CHEBU) ;
    
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
      
    if (dege==3) {
	double somme = 0 ;
	for (int i=0 ; i<n ; i++)
	    somme += i*i*res(i) ;
	double somme_deux = somme ;
	for (int i=0 ; i<n ; i++)
	    somme_deux -= res(i) ;
	res.set(1) = -somme ;
	res.set(0) = somme_deux ;
	return res ;
    }
    else {
	double somme = 0 ;
	for (int i=0 ; i<n ; i++)
	    somme += res(i) ;
	res.set(0) = -somme ;
	return res ;
    }
}	

// Cas dzpuis = 3 ;
Tbl _solp_poisson_2d_r_chebu_trois (const Matrice &lap, const Matrice &nondege, 
				    double alpha, const Tbl &source) {
    
    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;
    assert (dege ==2) ;
    
    Tbl source_aux(source*alpha) ;
    source_aux = cl_poisson_2d(source_aux, 3, R_CHEBU) ;
    
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
      
    double somme = 0 ;
    for (int i=0 ; i<n ; i++)
	somme += res(i) ;
    res.set(0) = -somme ;
    return res ;
}

	
// Cas dzpuis = 2 ;
Tbl _solp_poisson_2d_r_chebu_deux (const Matrice &lap, const Matrice &nondege, 
				   const Tbl &source) {
    
    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;
    assert ((dege==1) || (dege ==2)) ;
    Tbl source_aux(cl_poisson_2d(source, 2, R_CHEBU)) ;
    
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
    
    if (dege == 2) {
	double somme = 0 ;
	for (int i=0 ; i<n ; i++)
	    somme+=res(i) ;
    
	res.set(0) = -somme ;
    }
    
    return res ;
}


Tbl Ope_poisson_2d::get_solp (const Tbl& so) const {

  if (non_dege == 0x0)
    do_non_dege() ;

  // Routines de derivation
  static Tbl (*solp_poisson_2d[MAX_BASE]) (const Matrice&, const Matrice&,
					   double, double,const Tbl&, int) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      solp_poisson_2d[i] = _solp_poisson_2d_pas_prevu ;
    }
    // Les routines existantes
    solp_poisson_2d[R_CHEB >> TRA_R] = _solp_poisson_2d_r_cheb ;
    solp_poisson_2d[R_CHEBP >> TRA_R] = _solp_poisson_2d_r_chebp ;
    solp_poisson_2d[R_CHEBI >> TRA_R] = _solp_poisson_2d_r_chebi ;
    solp_poisson_2d[R_CHEBU >> TRA_R] = _solp_poisson_2d_r_chebu ;
  }
  
  Tbl res(solp_poisson_2d[base_r] (*ope_mat, *non_dege, 
				   alpha, beta, so, dzpuis)) ;
  
  Tbl valeurs (val_solp (res, alpha, base_r)) ;
  valeurs *= sqrt(double(2)) ;
  sp_plus = valeurs(0) ;
  sp_minus = valeurs(1) ;
  dsp_plus = valeurs(2) ;
  dsp_minus = valeurs(3) ;
  
  return res ;
}
