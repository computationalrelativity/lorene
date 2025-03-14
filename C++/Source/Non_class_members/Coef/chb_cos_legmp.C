/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *                 2009 Jerome Novak
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


 

/*
 *  Calcule les coefficients du developpement (suivant theta) en fonctions
 *  associees de Legendre P_l^m(cos(theta)) a partir des coefficients du 
 *  developpement en cos(j*theta) 
 *  representant une fonction 3-D symetrique par le retournement 
 *  (x, y, z) --> (-x, -y, z). 
 * 
 * Entree:
 * -------
 *  const int* deg : tableau du nombre effectif de degres de liberte dans chacune 
 *		     des 3 dimensions: 
 *			deg[0] = np : nombre de points de collocation en phi
 *			deg[1] = nt : nombre de points de collocation en theta
 *			deg[2] = nr : nombre de points de collocation en r
 *
 *  const double* cfi :  tableau des coefficients c_j du develop. en cos defini
 *			  comme suit (a r et phi fixes)
 *
 *			f(theta) = som_{j=0}^{nt-1} c_j cos( j theta ) 
 *
 * 		    L'espace memoire correspondant au pointeur cfi doit etre 
 *	            nr*nt*(np+2) et doit avoir ete alloue avant 
 *		    l'appel a la routine.	 
 *		    Le coefficient c_j (0 <= j <= nt-1) doit etre stoke dans le 
 *		    tableau cfi comme suit
 *		          c_j = cfi[ nr*nt* k + i + nr* j ]
 *		    ou k et i sont les indices correspondant a
 *		    phi et r respectivement.
 *
 * Sortie:
 * -------
 *   double* cfo :  tableau des coefficients a_l du develop. en fonctions de
 *		    Legendre associees P_n^m:
 *
 *			f(theta) = 
 *			    som_{l=m}^{nt-1} a_l P_l^m( cos(theta) )
 *		
 *		    avec m pair :  m = 0, 2, ..., np.		  
 *
 *		    P_n^m(x) represente la fonction de Legendre associee
 *		       de degre n et d'ordre m normalisee de facon a ce que
 *
 *			int_0^pi [ P_n^m(cos(theta)) ]^2  sin(theta) dtheta = 1
 *
 * 		    L'espace memoire correspondant au pointeur cfo doit etre 
 *	            nr*nt*(np+2) et doit avoir ete alloue avant 
 *		    l'appel a la routine.	 
 *		    Le coefficient a_l (0 <= l <= nt-1) est stoke dans le 
 *		    tableau cfo comme suit
 *		          a_l = cfo[ nr*nt* k + i + nr* l ]
 *		    ou k et i sont les indices correspondant a phi et r 
 *		    respectivement: m = 2( k/2 ).
 *		    NB: pour l < m,  a_l = 0
 *
 * NB:
 * ---
 *  Il n'est pas possible d'avoir le pointeur cfo egal a cfi.
 */

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:18:00  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:10  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:15:59  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2009/10/23 12:54:47  j_novak
 * New base T_LEG_MI
 *
 * Revision 1.1  2009/10/13 13:49:36  j_novak
 * New base T_LEG_MP.
 *
 *
 * $Header$
 *
 */


// headers du C
#include <cassert>
#include <cstdlib>

// Prototypage
#include "headcpp.h"
#include "proto.h"

namespace Lorene {
//******************************************************************************

void chb_cos_legmp(const int* deg , const double* cfi, double* cfo) {

int k2, l, jmin, j, i, m ;
 
// Nombres de degres de liberte en phi et theta :
    int np = deg[0] ;
    int nt = deg[1] ;
    int nr = deg[2] ;

    assert(np < 4*nt) ;

    // Tableau de travail
    double* som = new double[nr] ;

// Recherche de la matrice de passage  cos --> Legendre 
    double* aa = mat_cos_legmp(np, nt) ;
    
// Increment en m pour la matrice aa :
    int maa = nt * nt  ;
   
// Pointeurs de travail :
    double* resu = cfo ;
    const double* cc = cfi ;

// Increment en phi :
    int ntnr = nt * nr ;

// Indice courant en phi :
    int k = 0 ;

//----------------------------------------------------------------
//		    Cas axisymetrique       
//----------------------------------------------------------------

    if (np == 1) {

	m = 0 ; 

// Boucle sur l'indice l du developpement en Legendre

// ... produit matriciel (parallelise sur r)
		for (l=m; l<nt; l++) {
		    for (i=0; i<nr; i++) {
			som[i] = 0 ; 
		    }
		    
		    jmin =  l  ;    // pour m=0, aa_lj = 0 pour j<l
		    for (j=jmin; j<nt; j++) {
			double amlj = aa[nt*l + j] ;
			for (i=0; i<nr; i++) {
			    som[i] += amlj * cc[nr*j + i] ;
			}
		    }
		    
		    for (i=0; i<nr; i++) {
			*resu = som[i]  ;
			resu++ ;  
		    }
		    
		}  // fin de la boucle sur l 
	
	// Mise a zero des coefficients k=1 et k=2 :
	// ---------------------------------------
	
	for (i=ntnr; i<3*ntnr; i++) {
	    cfo[i] = 0 ;		 
	}	    
	    

	// on sort
	delete [] som ;
	return ; 
	
    }	// fin du cas np=1


//----------------------------------------------------------------
//		    Cas 3-D     
//----------------------------------------------------------------


// Boucle sur phi  : 


    for (m=0; m < np + 1 ; m+=2) {	    
    	
	for (k2=0; k2 < 2; k2++) {  // k2=0 : cos(m phi)  ;   k2=1 : sin(m phi)

	    if ( (k == 1) || (k == np+1) ) {	// On met les coef de sin(0 phi)
						// et sin( np phi)  a zero 
		for (l=0; l<nt; l++) {
		    for (i=0; i<nr; i++) {
			*resu = 0 ;
			resu++ ; 
		    }		    
		}
	    }
	    else {

// Boucle sur l'indice l du developpement en Legendre

		int lmax = (m<nt-1 ? m : nt-1) ;
		for (l=0; l<lmax; l++) {
		    for (i=0; i<nr; i++) {
			*resu = 0 ;
			resu++ ; 
		    }
		}		
// ... produit matriciel (parallelise sur r)
		for (l=m; l<nt; l++) {
		    for (i=0; i<nr; i++) {
			som[i] = 0 ; 
		    }
		    
		    jmin = ( m == 0 ) ? l : 0 ;  // pour m=0, aa_lj = 0 pour j<l
		    for (j=jmin; j<nt; j++) {
			double amlj = aa[nt*l + j] ;
			for (i=0; i<nr; i++) {
			    som[i] += amlj * cc[nr*j + i] ;
			}
		    }
		    
		    for (i=0; i<nr; i++) {
			*resu = som[i]  ;
			resu++ ;  
		    }
		    
		}  // fin de la boucle sur l 

	    }	// fin du cas k != 1 et k!=np+1
	    
// On passe au phi suivant :
	    cc = cc + ntnr	; 
	    k++ ;
	    	    
	}   // fin de la boucle sur k2 
	
// On passe a l'harmonique en phi suivante :

	aa += maa ;	// pointeur sur la nouvelle matrice de passage
		
    }	// fin de la boucle (m) sur phi  

//## verif : 
    assert(resu == cfo + (np+2)*ntnr) ;

    // Menage
    delete [] som ;
    
}
}
