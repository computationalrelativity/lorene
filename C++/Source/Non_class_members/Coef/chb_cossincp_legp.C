/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char chb_cossincp_legp_C[] = "$Header$" ;

/*
 *  Calcule les coefficients du developpement (suivant theta) en fonctions
 *  associees de Legendre P_l^m(cos(theta)) a partir des coefficients du 
 *  developpement en cos(2*j*theta) [m pair] / sin( (2*j+1) * theta) [m impair]
 *  representant une fonction 3-D symetrique par rapport au plan equatorial
 *  z = 0.
 * 
 * Entree:
 * -------
 *  const int* deg : tableau du nombre effectif de degres de liberte dans chacune 
 *		     des 3 dimensions: 
 *			deg[0] = np : nombre de points de collocation en phi
 *			deg[1] = nt : nombre de points de collocation en theta
 *			deg[2] = nr : nombre de points de collocation en r
 *
 *  const double* cfi :  tableau des coefficients c_j du develop. en cos/sin definis
 *			  comme suit (a r et phi fixes)
 *
 *	pour m pair:	f(theta) = som_{j=0}^{nt-1} c_j cos( 2 j theta ) 
 *			  
 *	pour m impair:	f(theta) = som_{j=0}^{nt-2} c_j sin( (2 j+1) theta ) 
 *
 * 		    L'espace memoire correspondant au pointeur cfi doit etre 
 *	            nr*nt*(np+2) et doit avoir ete alloue avant 
 *		    l'appel a la routine.	 
 *		    Le coefficient c_j (0 <= j <= nt-1) doit etre stoke dans le 
 *		    tableau cfi comme suit
 *		          c_j = cfi[ nr*nt* k + i + nr* j ]
 *		    ou k et i sont les indices correspondant a
 *		    phi et r respectivement: m = k/2.
 *	            Pour m impair, c_0 = c_{nt-1} = 0.
 *
 * Sortie:
 * -------
 *   double* cfo :  tableau des coefficients a_l du develop. en fonctions de
 *		    Legendre associees P_n^m:
 *
 *	pour m pair:	f(theta) = 
 *			    som_{l=m/2}^{nt-1} a_l P_{2l}^m( cos(theta) )
 *			  
 *	pour m impair:  f(theta) = 
 *			    som_{l=(m-1)/2}^{nt-2} a_l P_{2l+1}^m( cos(theta) )
 *
 *		    ou P_n^m(x) represente la fonction de Legendre associee
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
 *		    respectivement: m = k/2.
 *		    NB: pour m pair et l < m/2,  a_l = 0
 *			pour m impair et l < (m-1)/2,  a_l = 0
 *
 * NB:
 * ---
 *  Il n'est pas possible d'avoir le pointeur cfo egal a cfi.
 */

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  1999/02/22  15:45:31  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers du C++
#include <iostream.h>

// headers du C
#include <assert.h>
#include <malloc.h>
#include <stdlib.h>

// Prototypage
#include "proto.h"

//******************************************************************************

void chb_cossincp_legp(const int* deg , const double* cfi, double* cfo) {

// Espace de travail realloue eventuellement a chaque appel :

int ip, k2, l, jmin, j, i, m ;
 
// Nombres de degres de liberte en phi et theta :
    int np = deg[0] ;
    int nt = deg[1] ;
    int nr = deg[2] ;

    assert(np < 4*nt) ;

    // Tableau de travail
    double* som = (double*)( malloc( nr*sizeof(double) ) ) ;

// Recherche de la matrice de passage  cos/sin --> Legendre 
    double* aa = mat_cossincp_legp(np, nt) ;
    
// Increment en m pour la matrice aa :
    int maa = nt * nt  ;

//## Test
//    double* aat = aa ;
//    for ( m=0; m < np/2+1 ; m++) {
//	cout << "---------------------------------------" << endl ; 
//	cout << " m = " << m << endl ; 
//	cout << " " << endl ; 
//
//	for (l=m/2; l<nt; l++) {	
//	    cout << " l = " << l << " : " ; 	    
//	    for (j=0; j<nt; j++) {
//		cout << aat[nt*l + j] << " " ;
//	    }
//	    cout << endl ; 
//	}
//	arrete() ; 
//	aat += maa ;	// pointeur sur la nouvelle matrice de passage	
//    }			
//##    
   
// Pointeurs de travail :
    double* resu = cfo ;
    const double* cc = cfi ;

// Increment en phi :
    int ntnr = nt * nr ;

// Indice courant en phi :
    int k = 0 ;

// Ordre des harmoniques du developpement de Fourier en phi :
    m = 0 ;	    
      
// --------------
// Boucle sur phi  : k =    4*ip    4*ip+1   4*ip+2    4*ip+3
// --------------    m =    2*ip     2*ip    2*ip+1    2*ip+1
//		     k2 =     0	      1	       0	 1

    for (ip=0; ip < np/4 + 1 ; ip++) {	    
    
//--------------------------------
//  Partie  m pair
//--------------------------------


	for (k2=0; k2 < 2; k2++) {  // k2=0 : cos(m phi)  ;   k2=1 : sin(m phi)
	
	    if ( (k == 1) || (k == np+1) ) {	// On met les coef de sin(0 phi)
						// et sin( np/2 phi)  a zero 
		for (l=0; l<nt; l++) {
		    for (i=0; i<nr; i++) {
			*resu = 0 ;
			resu++ ; 
		    }		    
		}
	    }
	    else {

// Boucle sur l'indice l du developpement en Legendre

		for (l=0; l<m/2; l++) {
		    for (i=0; i<nr; i++) {
			*resu = 0 ;
			resu++ ; 
		    }
		}		
// ... produit matriciel (parallelise sur r)
		for (l=m/2; l<nt; l++) {
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

	    }	// fin du cas k != 1 
	    
// On passe au phi suivant :
	    cc = cc + ntnr	; 
	    k++ ;
	    	    
	}   // fin de la boucle sur k2 
	
// On passe a l'harmonique en phi suivante :
	m++ ;
	aa += maa ;	// pointeur sur la nouvelle matrice de passage
	
//--------------------------------
//  Partie m impair
//--------------------------------

	for (k2=0; k2 < 2; k2++) {  // k2=0 : cos(m phi)  ;   k2=1 : sin(m phi)
	
	    if ( k == np+1 ) {			// On met les coef de 
						// sin( np/2 phi)  a zero 
		for (l=0; l<nt; l++) {
		    for (i=0; i<nr; i++) {
			*resu = 0 ;
			resu++ ; 
		    }		    
		}
	    }

	    if (k < np+1) {  

// Boucle sur l'indice l du developpement en Legendre
		for (l=0; l<(m-1)/2; l++) {
		    for (i=0; i<nr; i++) {
			*resu = 0 ;
			resu++ ; 
		    }
		}

// ... produit matriciel (parallelise sur r)
		for (l=(m-1)/2; l<nt-1; l++) {
		    for (i=0; i<nr; i++) {
			som[i] = 0 ; 
		    }
		   
		    jmin = ( m == 1 ) ? l : 0 ;  // pour m=1, aa_lj = 0 pour j<l

		    for (j=jmin; j<nt-1; j++) {
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

// Dernier coef en l=nt-1 mis a zero pour le cas m impair : 
		for (i=0; i<nr; i++) {
		    *resu = 0  ;
		    resu++ ;  
		}
		

// On passe au phi suivant :
		cc = cc + ntnr	; 
		k++ ;
		
	    }	// fin du cas k < np+1
	    
	}  //  fin de la boucle sur k2 
	

// On passe a l'harmonique en phi suivante :
	m++ ;
	aa += maa ;	// pointeur sur la nouvelle matrice de passage 
	
    }	// fin de la boucle (ip) sur phi  

// Mise a zero des coefficients de sin( np/2 phi ) (k=np+1)

//## verif : 
//    assert(resu == cfo + (np+2)*ntnr) ;

    // Menage
    free (som) ;
    
}
