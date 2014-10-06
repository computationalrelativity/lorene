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


char chb_cossinci_legi_C[] = "$Header$" ;

/*
 *  Calcule les coefficients du developpement (suivant theta) en fonctions
 *  associees de Legendre P_l^m(cos(theta)) (l-m impair) a partir des 
 *  coefficients du developpement en 
 *   cos[(2*j+1)*theta] [m pair] / sin(2*j* theta) [m impair]
 *  representant une fonction 3-D antisymetrique par rapport au plan equatorial
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
 *	pour m pair:	f(theta) = som_{j=0}^{nt-2} c_j cos( (2 j+1) theta ) 
 *			  
 *	pour m impair:	f(theta) = som_{j=1}^{nt-2} c_j sin( 2j theta ) 
 *
 * 		    L'espace memoire correspondant au pointeur cfi doit etre 
 *	            nr*nt*(np+2) et doit avoir ete alloue avant 
 *		    l'appel a la routine.	 
 *		    Le coefficient c_j (0 <= j <= nt-1) doit etre stoke dans le 
 *		    tableau cfi comme suit
 *		          c_j = cfi[ nr*nt* k + i + nr* j ]
 *		    ou k et i sont les indices correspondant a
 *		    phi et r respectivement: m = k/2.
 *	            Pour m pair, c_{nt-1} = 0.
 *	            Pour m impair, c_0 = c_{nt-1} = 0.
 *
 * Sortie:
 * -------
 *   double* cfo :  tableau des coefficients a_l du develop. en fonctions de
 *		    Legendre associees P_n^m impaires:
 *
 *	pour m pair:	f(theta) = 
 *			    som_{l=m/2}^{nt-2} a_l P_{2l+1}^m( cos(theta) )
 *			  
 *	pour m impair:  f(theta) = 
 *			    som_{l=(m+1)/2}^{nt-2} a_l P_{2l}^m( cos(theta) )
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
 *		    NB: pour m pair: si l < m/2 ou l = nt-1,  a_l = 0
 *			pour m impair: si l < (m+1)/2 ou l = nt-1,  a_l = 0
 *
 * NB:
 * ---
 *  Il n'est pas possible d'avoir le pointeur cfo egal a cfi.
 */

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/06 15:16:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2005/02/18 13:14:10  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
 * Revision 1.2  2002/10/16 14:36:52  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/02/22  15:45:41  hyc
 * *** empty log message ***
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

//******************************************************************************

void chb_cossinci_legi(const int* deg , const double* cfi, double* cfo) {

// Espace de travail realloue eventuellement a chaque appel :

int ip, k2, l, jmin, j, i, m ;
 
// Nombres de degres de liberte en phi et theta :
    int np = deg[0] ;
    int nt = deg[1] ;
    int nr = deg[2] ;

    assert(np < 4*nt) ;

    // Tableau de travail
    double* som = new double[nr] ; 

// Recherche de la matrice de passage  cos/sin --> Legendre 
    double* aa = mat_cossinci_legi(np, nt) ;
    
// Increment en m pour la matrice aa :
    int maa = nt * nt  ;
   
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

	 //...  0 <= l <= m/2 - 1 : a_l = 0 
		for (l=0; l<m/2; l++) {
		    for (i=0; i<nr; i++) {
			*resu = 0 ;
			resu++ ; 
		    }
		}		

	 //...  m/2 <= l <= nt-2 :  produit matriciel (parallelise sur r)
		for (l=m/2; l<nt-1; l++) {
		    for (i=0; i<nr; i++) {
			som[i] = 0 ; 
		    }
		    
		    //## jmin = ( m == 0 ) ? l : 0 ;  // pour m=0, aa_lj = 0 pour j<l
		    jmin = 0 ; 
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


	 //...  l = nt-1 : a_l = 0 
		for (i=0; i<nr; i++) {
		    *resu = 0 ;
		    resu++ ; 
		}

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

	 //...  0 <= l <= (m-1)/2  : a_l = 0 

		for (l=0; l<(m+1)/2; l++) {
		    for (i=0; i<nr; i++) {
			*resu = 0 ;
			resu++ ; 
		    }
		}

	 //...  (m+1)/2 <= l <= nt-2 :  produit matriciel (parallelise sur r)
		for (l=(m+1)/2; l<nt-1; l++) {
		    for (i=0; i<nr; i++) {
			som[i] = 0 ; 
		    }
		   
		    //## jmin = ( m == 1 ) ? l : 0 ;  // pour m=1, aa_lj = 0 pour j<l
		    jmin = 1 ;	// pour j=0, sin(2*j*theta)=0  

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

	 //...  l = nt-1 : a_l = 0 
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
    assert(resu == cfo + (np+2)*ntnr) ;

    // Menage
    delete [] som ;
    
}
