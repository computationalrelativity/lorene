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


char solp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2007/12/12 12:30:49  jl_cornou
 * *** empty log message ***
 *
 * Revision 1.6  2004/10/05 15:44:21  j_novak
 * Minor speed enhancements.
 *
 * Revision 1.5  2004/02/20 10:55:23  j_novak
 * The versions dzpuis 5 -> 3 has been improved and polished. Should be
 * operational now...
 *
 * Revision 1.4  2004/02/09 08:55:31  j_novak
 * Corrected error in the arguments of _solp_r_chebu_cinq
 *
 * Revision 1.3  2004/02/06 10:53:54  j_novak
 * New dzpuis = 5 -> dzpuis = 3 case (not ready yet).
 *
 * Revision 1.2  2002/10/16 14:37:12  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.14  2000/09/07  12:54:42  phil
 * *** empty log message ***
 *
 * Revision 2.13  2000/05/22  13:41:50  phil
 * ajout du cas dzpuis == 3 ;
 *
 * Revision 2.12  1999/11/30  17:45:04  phil
 * changement prototypage
 *
 * Revision 2.11  1999/10/12  09:44:44  phil
 * *** empty log message ***
 *
 * Revision 2.10  1999/10/12  09:37:54  phil
 * passage en const Matreice&
 *
 * Revision 2.9  1999/07/08  09:54:58  phil
 * changement appel de multx_1d
 *
 * Revision 2.8  1999/07/07  10:05:20  phil
 * Passage aux operateurs 1d
 * /
 *
 * Revision 2.7  1999/07/02  15:04:48  phil
 * *** empty log message ***
 *
 * Revision 2.6  1999/06/23  16:21:59  phil
 * *** empty log message ***
 *
 * Revision 2.5  1999/06/23  12:35:06  phil
 * ajout de dzpuis = 2
 *
 * Revision 2.4  1999/04/07  15:07:18  phil
 * *** empty log message ***
 *
 * Revision 2.3  1999/04/07  15:06:14  phil
 * Changement de calcul pour (-1)^n
 *
 * Revision 2.2  1999/04/07  14:55:46  phil
 * Changement prototypage
 *
 * Revision 2.1  1999/04/07  14:36:38  phil
 * passage par reference
 *
 * Revision 2.0  1999/04/07  14:11:24  phil
 * *** empty log message ***
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

/*
 * Calcul une solution particuliere a : Lap X = Y
 * 
 * Entree :
 *	lap : matrice de l'operateur
 *	nondege : matrice non degeneree associee
 *	alpha et beta : voire mapping
 *	source : Tbl contenant les coefficients en r de Y
 *	puis : puissance dans la ZEC
 * Sortie :
 *	Tbl contenant les coefficients de X
 * 
 * 
 */
		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _solp_pas_prevu (const Matrice &lap, const Matrice &nondege, double alpha, 
		    double beta, const Tbl &source, int puis) {
    cout << " Solution particuliere pas prevue ..... : "<< endl ;
    cout << " Laplacien : " << endl << lap << endl ;
    cout << " Non degenere  : " << endl << nondege << endl ;
    cout << " Source : " << &source << endl ;
    cout << " Alpha : " << alpha << endl ;
    cout << " Beta : " << beta << endl ;
    cout << " Puiss : " << puis << endl ;
    abort() ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	
		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _solp_r_cheb (const Matrice &lap, const Matrice &nondege, double alpha, 
		  double beta, const Tbl &source, int) {
    
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
    source_aux = combinaison(source_aux, 0, R_CHEB) ;
	
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
	       //--  R_JACO02 ------
	      //-------------------

Tbl _solp_r_jaco02 (const Matrice &lap, const Matrice &nondege, double alpha, 
		  double beta, const Tbl &source, int) {
    
    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;
    assert (dege ==2) ;
    
    Tbl source_aux(source*alpha*alpha) ;
    source_aux = combinaison(source_aux, 0, R_JACO02) ;
	
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

Tbl _solp_r_chebp (const Matrice &lap, const Matrice &nondege, double alpha, 
		    double , const Tbl &source, int) {
    
    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;
    assert ((dege==2) || (dege == 1)) ;
    Tbl source_aux(alpha*alpha*source) ;
    source_aux = combinaison(source_aux, 0, R_CHEBP) ;
    
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
	
Tbl _solp_r_chebi (const Matrice &lap, const Matrice &nondege, double alpha, 
		    double,  const Tbl &source, int) {
   

    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;
    assert ((dege == 2) || (dege ==1)) ;
    Tbl source_aux(source*alpha*alpha) ;
    source_aux = combinaison(source_aux, 0, R_CHEBI) ;
    
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

Tbl _solp_r_chebu (const Matrice &lap, const Matrice &nondege, double alpha, 
		    double, const Tbl &source, int puis) {
    int n = lap.get_dim(0) ;
    Tbl res(n) ;
    res.set_etat_qcq() ;
    
    switch (puis) {
	case 5 :
	    res = _solp_r_chebu_cinq (lap, nondege, source) ;
	    break ;
	case 4 :
	    res = _solp_r_chebu_quatre (lap, nondege, alpha, source) ;
	    break ;
	case 3 :
	    res = _solp_r_chebu_trois (lap, nondege, alpha, source) ;
	    break ;
	case 2 :
	    res = _solp_r_chebu_deux (lap, nondege, source) ;
	    break ;
	default :
	    abort() ;
	    exit(-1) ;
	}
return res ;
}


// Cas dzpuis = 4 ;
Tbl _solp_r_chebu_quatre (const Matrice &lap, const Matrice &nondege, double alpha, 
			    const Tbl &source) {
    
    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;
    assert ((dege==3) || (dege ==2)) ;
    Tbl source_aux(source*alpha*alpha) ;
    source_aux = combinaison(source_aux, 4, R_CHEBU) ;
    
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
Tbl _solp_r_chebu_trois (const Matrice &lap, const Matrice &nondege, double alpha, 
			    const Tbl &source) {
    
    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;
    assert (dege ==2) ;
    
    Tbl source_aux(source*alpha) ;
    source_aux = combinaison(source_aux, 3, R_CHEBU) ;
    
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
Tbl _solp_r_chebu_deux (const Matrice &lap, const Matrice &nondege, 
			const Tbl &source) {
    
    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;
    assert ((dege==1) || (dege ==2)) ;
    Tbl source_aux(combinaison(source, 2, R_CHEBU)) ;
    
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

// Cas dzpuis = 5 ;
Tbl _solp_r_chebu_cinq (const Matrice &lap, const Matrice &nondege, 
			const Tbl &source) {
    
    int n = lap.get_dim(0) ;	  
    int dege = n-nondege.get_dim(0) ;

    Tbl source_aux(combinaison(source, 5, R_CHEBU)) ;
    
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


	      	//-------------------
	       //--  Fonction   ----
	      //-------------------
	      
	      
Tbl solp(const Matrice &lap, const Matrice &nondege, double alpha,
	    double beta, const Tbl &source, int puis, int base_r) {

		// Routines de derivation
    static Tbl (*solp[MAX_BASE])(const Matrice&, const Matrice&, double, double,
				     const Tbl&, int) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    solp[i] = _solp_pas_prevu ;
	}
		// Les routines existantes
	solp[R_CHEB >> TRA_R] = _solp_r_cheb ;
	solp[R_CHEBU >> TRA_R] = _solp_r_chebu ;
	solp[R_CHEBP >> TRA_R] = _solp_r_chebp ;
	solp[R_CHEBI >> TRA_R] = _solp_r_chebi ;
	solp[R_JACO02 >> TRA_R] = _solp_r_jaco02 ;
    }
    
    return solp[base_r](lap, nondege, alpha, beta, source, puis) ;
}
