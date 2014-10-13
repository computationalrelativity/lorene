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


char chb_cosp_legpp_C[] = "$Header$" ;

/*
 *  Calcule les coefficients du developpement (suivant theta) en fonctions
 *  associees de Legendre P_l^m(cos(theta)) a partir des coefficients du 
 *  developpement en cos(2*j*theta) 
 *  representant une fonction 3-D symetrique par rapport au plan equatorial
 *  z = 0 et symetrique par le retournement (x, y, z) --> (-x, -y, z). 
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
 *			f(theta) = som_{j=0}^{nt-1} c_j cos( 2 j theta ) 
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
 *			    som_{l=m/2}^{nt-1} a_l P_{2l}^m( cos(theta) )
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
 *		    NB: pour l < m/2,  a_l = 0
 *
 * NB:
 * ---
 *  Il n'est pas possible d'avoir le pointeur cfo egal a cfi.
 */

/*
 * $Id$
 * $Log$
 * Revision 1.7  2014/10/13 08:53:10  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:16:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2005/02/18 13:14:10  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
 * Revision 1.4  2003/12/19 16:21:46  j_novak
 * Shadow hunt
 *
 * Revision 1.3  2003/01/31 10:31:23  e_gourgoulhon
 * Suppressed the directive #include <malloc.h> for malloc is defined
 * in <stdlib.h>
 *
 * Revision 1.2  2002/10/16 14:36:52  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  2000/09/29  16:06:31  eric
 * Mise a zero des coefficients k=1 et k=2.
 *
 * Revision 2.0  1999/02/22  15:45:54  hyc
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

namespace Lorene {
//******************************************************************************

void chb_cosp_legpp(const int* deg , const double* cfi, double* cfo) {

int k2, l, jmin, j, i, m ;
 
// Nombres de degres de liberte en phi et theta :
    int np = deg[0] ;
    int nt = deg[1] ;
    int nr = deg[2] ;

    assert(np < 4*nt) ;

    // Tableau de travail
    double* som = new double[nr] ;

// Recherche de la matrice de passage  cos --> Legendre 
    double* aa = mat_cosp_legpp(np, nt) ;
    
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
		for (l=m/2; l<nt; l++) {
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
