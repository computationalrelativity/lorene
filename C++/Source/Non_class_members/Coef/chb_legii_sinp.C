/*
 *   Copyright (c) 2003 Jerome Novak
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


char chb_legii_sinp_C[] = "$Header$" ;

/*
 *  Calcule les coefficients du developpement (suivant theta) 
 *  en sin(2j theta) 
 *  a partir des coefficients du developpement en fonctions
 *  associees de Legendre P_l^m(cos(theta)) (l pair et m impair)
 *  pour une une fonction 3-D antisymetrique par rapport au plan equatorial
 *  z = 0 et antisymetrique par le retournement (x, y, z) --> (-x, -y, z). 
 * 
 * Entree:
 * -------
 *  const int* deg : tableau du nombre effectif de degres de liberte dans chacune 
 *		     des 3 dimensions: 
 *			deg[0] = np : nombre de points de collocation en phi
 *			deg[1] = nt : nombre de points de collocation en theta
 *			deg[2] = nr : nombre de points de collocation en r
 *
 *  const double* cfi : tableau des coefficients a_j du develop. en fonctions de
 *		    Legendre associees P_n^m:
 *
 *		    f(theta) = 
 *			    som_{l=(m+1)/2}^{nt-2} a_j P_{2j}^m( cos(theta) )
 *			  
 *		(m impair) 
 *
 *		    ou P_l^m(x) represente la fonction de Legendre associee
 *		       de degre l et d'ordre m normalisee de facon a ce que
 *
 *			int_0^pi [ P_l^m(cos(theta)) ]^2  sin(theta) dtheta = 1
 *
 * 		    L'espace memoire correspondant au pointeur cfi doit etre 
 *	            nr*nt*(np+2) et doit avoir ete alloue avant 
 *		    l'appel a la routine.	 
 *		    Le coefficient a_j (0 <= j <= nt-1) doit etre stoke dans le 
 *		    tableau cfi comme suit
 *		          a_j = cfi[ nr*nt* k + i + nr* j ]
 *		    ou k et i sont les indices correspondant a phi et r 
 *		    respectivement: m = 2 (k/2).
 *		    NB: pour j<(m+1)/2,  a_j = 0
 *
 * Sortie:
 * -------
 *   double* cfo :  tableau des coefficients c_j du develop. en sin definis
 *			  comme suit (a r et phi fixes) :
 *
 *			f(theta) = som_{j=1}^{nt-2} c_j sin( 2j theta ) 
 *			  
 * 		    L'espace memoire correspondant au pointeur cfo doit etre 
 *	            nr*nt*(np+2) et doit avoir ete alloue avant 
 *		    l'appel a la routine.	 
 *		    Le coefficient c_j (0 <= j <= nt-1) est stoke dans le 
 *		    tableau cfo comme suit
 *		          c_j = cfo[ nr*nt* k + i + nr* j ]
 *		    ou k et i sont les indices correspondant a
 *		    phi et r respectivement.
 *	            NB:	c_0 =  c_{nt-1} = 0.
 *
 *
 * NB:
 * ---
 *  Il n'est pas possible d'avoir le pointeur cfo egal a cfi.
 */

/*
 * $Id$
 * $Log$
 * Revision 1.2  2005/02/18 13:14:10  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
 * Revision 1.1  2003/09/16 08:58:01  j_novak
 * New functions for the T_LEG_II base
 *
 *
 * $Header$
 *
 */

// headers du C
#include <stdlib.h>
#include <assert.h>

// Headers Lorene
#include "headcpp.h"
#include "proto.h"

//******************************************************************************

void chb_legii_sinp(const int* deg , const double* cfi, double* cfo) {

int k2, l, j, i, m ;
 
// Nombres de degres de liberte en phi et theta :
    int np = deg[0] ;
    int nt = deg[1] ;
    int nr = deg[2] ;

    assert(np < 4*nt) ;
    assert( cfi != cfo ) ; 

    // Tableau de travail
    double* som = new double[nr] ;
// Recherche de la matrice de passage  Legendre -->  cos/sin 
    double* bb = mat_legii_sinp(np, nt) ;
    
// Increment en m pour la matrice bb :
    int mbb = nt * nt  ;
   
// Pointeurs de travail :
    double* resu = cfo ;
    const double* cc = cfi ;

// Increment en phi :
    int ntnr = nt * nr ;

// Indice courant en phi :
    int k = 0 ;

    // Cas k=0 (m=1 : cos(phi)) 
    // ------------------------

    //... premier coef en j=0 mis a zero: 
    for (i=0; i<nr; i++) {
      *resu = 0  ;
      resu++ ;  
    }

    // Boucle sur l'indice j du developpement en sin( 2j theta) 

    for (j=1; j<nt-1; j++) {

	// ... produit matriciel (parallelise sur r)
	for (i=0; i<nr; i++) {
	    som[i] = 0 ; 
	}

	for (l=1; l<nt-1; l++) {
	    double bmjl = bb[nt*j + l] ;
	    for (i=0; i<nr; i++) {
		som[i] += bmjl * cc[nr*l + i] ;
	    }
	}
		    
	for (i=0; i<nr; i++) {
	    *resu = som[i]  ;
	    resu++ ;  
	}
		    
    }  // fin de la boucle sur j 

    // Dernier coef en j=nt-1 mis a zero pour le cas m impair : 
    for (i=0; i<nr; i++) {
	*resu = 0  ;
	resu++ ;  
    }
		
    // Special case np=1 (axisymmetry)
    // -------------------------------
    if (np==1) {
	for (i=0; i<2*ntnr; i++) {
	    *resu = 0 ;
	    resu++ ; 
	}
	delete []  som  ; 
	return ; 			    
    }
	
    // On passe au phi suivant :
    cc = cc + ntnr	; 
    k++ ;
		
    // Cas k=1 : tout est mis a zero
    // -----------------------------	

    for (l=0; l<nt; l++) {
	for (i=0; i<nr; i++) {
	    *resu = 0 ;
	    resu++ ; 
	}		    
    }

    // On passe au phi suivant :
    cc = cc + ntnr	; 
    k++ ;
		
    // Cas k=2 (m=1 : sin(phi)) 
    // ------------------------

    //... premier coef en j=0 mis a zero: 
    for (i=0; i<nr; i++) {
      *resu = 0  ;
      resu++ ;  
    }

    // Boucle sur l'indice j du developpement en sin( 2j theta) 

    for (j=1; j<nt-1; j++) {

	// ... produit matriciel (parallelise sur r)
	for (i=0; i<nr; i++) {
	    som[i] = 0 ; 
	}

	for (l=1; l<nt-1; l++) {
	    double bmjl = bb[nt*j + l] ;
	    for (i=0; i<nr; i++) {
		som[i] += bmjl * cc[nr*l + i] ;
	    }
	}
		    
	for (i=0; i<nr; i++) {
	    *resu = som[i]  ;
	    resu++ ;  
	}
		    
    }  // fin de la boucle sur j 

    // Dernier coef en j=nt-1 mis a zero pour le cas m impair : 
    for (i=0; i<nr; i++) {
	*resu = 0  ;
	resu++ ;  
    }
		
    // On passe au phi suivant :
    cc = cc + ntnr	; 
    k++ ;
		
    // On passe au m suivant :
    bb += mbb ;	// pointeur sur la nouvelle matrice de passage 

    // Cas k >= 3
    // ----------

    for (m=3; m < np ; m+=2) {	    
    	
      for (k2=0; k2 < 2; k2++) {  // k2=0 : cos(m phi)  ;   k2=1 : sin(m phi)
	
	// Boucle sur l'indice j du developpement en sin( 2j theta) 

	//... premier coef en j=0 mis a zero: 
	for (i=0; i<nr; i++) {
	  *resu = 0  ;
	  resu++ ;  
	}

	for (j=1; j<nt-1; j++) {

	  for (i=0; i<nr; i++) {
	    som[i] = 0 ; 
	  }

	  for (l=(m+1)/2; l<nt-1; l++) {
	    double bmjl = bb[nt*j + l] ;
	    for (i=0; i<nr; i++) {
	      som[i] += bmjl * cc[nr*l + i] ;
	    }
	  }
	    
	  for (i=0; i<nr; i++) {
	    *resu = som[i]  ;
	    resu++ ;  
	  }
	    
	}  // fin de la boucle sur j 
	  
	// Dernier coef en j=nt-1 mis a zero pour le cas m impair : 
	for (i=0; i<nr; i++) {
	  *resu = 0  ;
	  resu++ ;  
	}
	  
	// On passe au phi suivant :
	cc = cc + ntnr	; 
	k++ ;
	  
      }  //  fin de la boucle sur k2 
	
      // On passe a l'harmonique en phi suivante :
      bb += mbb ;	// pointeur sur la nouvelle matrice de passage 
      
    }	// fin de la boucle (m) sur phi  


    // Cas k=np+1 : tout est mis a zero
    // --------------------------------	

    for (l=0; l<nt; l++) {
	for (i=0; i<nr; i++) {
	    *resu = 0 ;
	    resu++ ; 
	}		    
    }


//## verif : 
    assert(resu == cfo + (np+2)*ntnr) ;

    // Menage
    delete [] som ;
    
}
