/*
 *   Copyright (c) 1999-2002 Eric Gourgoulhon
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
 * Transformation de Jacobi (cas fin) sur le troisieme indice (indice
 *  correspondant a r) d'un tableau 3-D.
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en r est  nr = deg[2] et doit etre de la forme
 * 			nr = 2*p + 1
 *   int* dimf	: tableau du nombre d'elements de ff dans chacune des trois
 *	          dimensions.
 *		  On doit avoir  dimf[2] >= deg[2] = nr.
 *		  NB: pour dimf[0] = 1 (un seul point en phi), la transformation
 *		      est bien effectuee.
 *		      pour dimf[0] > 1 (plus d'un point en phi), la
 *		      transformation n'est effectuee que pour les indices (en phi)
 *		      j != 1 et j != dimf[0]-1 (cf. commentaires sur borne_phi).
 *
 *   double* ff : tableau des valeurs de la fonction aux nr points de
 *                        de collocation
 *
 *			  x_i =  pointsgausslobatto(nr-1)
 *
 *		    Les valeurs de la fonction doivent etre stockees dans le 
 *		    tableau ff comme suit
 *			   f( x_i ) = ff[ dimf[1]*dimf[2] * j + dimf[2] * k + i ]
 *		    ou j et k sont les indices correspondant a phi et theta 
 *		    respectivement.
 * 		    L'espace memoire correspondant a ce pointeur doit etre 
 *		    dimf[0]*dimf[1]*dimf[2] et doit avoir ete alloue avant 
 *		    l'appel a la routine.	 
 *
 *   int* dimc	: tableau du nombre d'elements de cf dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dimc[2] >= deg[2] = nr. 
 *
 * Sortie:
 * -------
 *   double* cf	:   tableau des coefficients c_i de la fonction definis
 *		    comme suit (a theta et phi fixes)
 *
 *				    f(x) = som_{i=0}^{nr-1} c_i J_i(x) , 
 *
 *		    ou J_i(x) designe le polynome de Jacobi d'indice (0,2) 
 *                  de degre i. 	 
 *		    Les coefficients c_i (0 <= i <= nr-1) sont stockees dans
 *		    le tableau cf comme suit
 *			   c_i = cf[ dimc[1]*dimc[2] * j + dimc[2] * k + i ]
 *		    ou j et k sont les indices correspondant a phi et theta 
 *		    respectivement.
 * 		    L'espace memoire correspondant au pointeur cf doit etre 
 *		    dimc[0]*dimc[1]*dimc[2] et doit avoir ete alloue avant 
 *		    l'appel a la routine.	 
 *
 * NB: Si le pointeur ff est egal a cf, la routine ne travaille que sur un 
 *     seul tableau, qui constitue une entree/sortie.
 */

// headers du C
#include <cstdlib>
#include <cassert>

#include "tbl.h"

// Prototypage des sous-routines utilisees:
namespace Lorene {
double* coeffjaco(int , double* ) ;

//*****************************************************************************

void cfrjaco02(const int* deg, const int* dimf, double* ff, const int* dimc,
		double* cf)

{

int i, j, k ;

// Dimensions des tableaux ff et cf  :
    int n1f = dimf[0] ;
    int n2f = dimf[1] ;
    int n3f = dimf[2] ;
    int n1c = dimc[0] ;
    int n2c = dimc[1] ;
    int n3c = dimc[2] ;

// Nombres de degres de liberte en r :    
    int nr = deg[2] ;
    
// Tests de dimension:
    if (nr > n3f) {
	cout << "cfrjaco02: nr > n3f : nr = " << nr << " ,  n3f = " 
	<< n3f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (nr > n3c) {
	cout << "cfrjaco02: nr > n3c : nr = " << nr << " ,  n3c = " 
	<< n3c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n1f > n1c) {
	cout << "cfrjaco02: n1f > n1c : n1f = " << n1f << " ,  n1c = " 
	<< n1c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n2f > n2c) {
	cout << "cfrjaco02: n2f > n2c : n2f = " << n2f << " ,  n2c = " 
	<< n2c << endl ;
	abort () ;
	exit(-1) ;
    }

// Nombre de points en radial:
    int nm1 = nr - 1 ;
    
// boucle sur phi et theta

    int n2n3f = n2f * n3f ;
    int n2n3c = n2c * n3c ;

/*   
 * Borne de la boucle sur phi: 
 *    si n1f = 1, on effectue la boucle une fois seulement.
 *    si n1f > 1, on va jusqu'a j = n1f-2 en sautant j = 1 (les coefficients
 *	j=n1f-1 et j=0 ne sont pas consideres car nuls). 
 */
    int borne_phi = ( n1f > 1 ) ? n1f-1 : 1 ;

    for (j=0; j< borne_phi; j++) {
    
	if (j==1) continue ;	// on ne traite pas le terme en sin(0 phi)

	for (k=0; k<n2f; k++) {

	    int i0 = n2n3f * j + n3f * k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3c * j + n3c * k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau resultat

	    double* aa = coeffjaco(nm1,ff0) ;
	    for ( i = 0; i < nr ; i++ ) {
	      cf0[i] = aa[i];
	    }   // fin de la boucle sur r
	    delete [] aa ;
	} 	// fin de la boucle sur theta 
   }	// fin de la boucle sur phi

}
}
