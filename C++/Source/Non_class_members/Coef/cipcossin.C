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


char cipcossin_C[] = "$Header$" ;


/*
 * Transformation de Fourier inverse sur le premier indice d'un tableau 3-D
 *  par le biais de la routine FFT Fortran FFT991
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune 
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en phi est  np = deg[0] et doit etre de la forme
 * 			np = 2^p 3^q 5^r 
 *   int* dimc	: tableau du nombre d'elements dans chacune des trois 
 *    		  dimensions du tableau cf
 *		  On doit avoir  dimc[0] >= deg[0] + 2  = np + 2. 
 *
 *   int* dimf	: tableau du nombre d'elements dans chacune des trois 
 *    		  dimensions du tableau ff
 *		  On doit avoir  dimf[0] >= deg[0]  = np . 
 *
 *
 * Entree / sortie :
 * -----------------
 *   double* cf : entree: tableau des coefficients de la fonction f; 
 *			  L'espace memoire correspondant a ce
 *                        pointeur doit etre dimc[0]*dimc[1]*dimc[2] et doit 
 *			  avoir ete alloue avant l'appel a la routine.
 *			  La convention de stokage est la suivante:
 *		cf[ dimc[2]*dimc[1]*k + dimc[2]*j + i ] = c_k      0<= k <= np,
 *			  ou les indices j et i correspondent respectivement
 *			  a theta et a r et ou les c_k sont les coefficients 
 *			  du developpement de f en series de Fourier:
 *
 *			f(phi) = som_{l=0}^{np/2} c_{2l} cos( 2 pi/T l phi )
 *				 + c_{2l+1} sin( 2 pi/T l phi ),
 *
 *  			  ou T est la periode de f.
 *			  En particulier cf[1] = 0.
 *			  De plus, cf[np+1] n'est pas egal a c_{np+1}
 *			  mais a zero.
 *   !!!! CE TABLEAU EST DETRUIT EN SORTIE !!!!!
 *
 * Sortie:
 * -------
 *  double* ff : tableau des valeurs de la fonction aux points de
 *                        de collocation
 *
 *				phi_k = T k/np      0 <= k <= np-1
 *
 *			  suivant la convention de stokage:
 *	ff[ dimf[2]*dimf[1]*k + dimf[2]*j + i ] = f(phi_k)    0 <= k <= np-1,
 *			 les indices j et i correspondant respectivement
 *			  a theta et a r. 
 *			  L'espace memoire correspondant a ce
 *                        pointeur doit etre dimf[0]*dimf[1]*dimf[2] et doit 
 *			  avoir ete alloue avant l'appel a la routine.
 */

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  1999/02/22  15:43:58  hyc
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
#include <stdarg.h>
#include <malloc.h>
#include <unistd.h>
#include <stdlib.h>

// Prototypage des sous-routines utilisees:
int*	facto_ini(int ) ;
double*	trigo_ini(int ) ;
extern "C" void fft991_(double [], double [], double [], int [],
                int *, int *, int *, int *, int *) ;
//*****************************************************************************

void cipcossin(const int* deg, const int* dimc, const int* dimf, 
	       double* cf, double* ff)
{

int i, j, k, index ;

// Dimensions des tableaux ff et cf  :
    int n1f = dimf[0] ;
    int n2f = dimf[1] ;
    int n3f = dimf[2] ;
    int n1c = dimc[0] ;
    int n2c = dimc[1] ;
    int n3c = dimc[2] ;

// Nombres de degres de liberte en phi:    
    int np = deg[0] ;

// Tests de dimension:
    if (np+2 > n1c) {
	cout << "cipcossin: np+2 > n1c : np = " << np << " ,  n1c = " 
	<< n1c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (np > n1f) {
	cout << "cipcossin: np > n1f : np = " << np << " ,  n1f = " 
	<< n1f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n3f > n3c) {
	cout << "cipcossin: n3f > n3c : n3f = " << n3f << " ,  n3c = " 
	<< n3c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n2f > n2c) {
	cout << "cipcossin: n2f > n2c : n2f = " << n2f << " ,  n2c = " 
	<< n2c << endl ;
	abort () ;
	exit(-1) ;
    }

    // Recherche des tables
    int* facto = facto_ini(np) ;
    double* trigo = trigo_ini(np) ;
	
    // Tableau de travail
    double* t1 = (double*)( malloc( (np+2)*sizeof(double) ) ) ;

// Denormalisation des cosinus
    int n2n3c = n2c * n3c ;
    for (i=2; i<np; i += 2 )	{
    	for (j=0; j<n2c; j++) {
	    for (k=0; k<n3c; k++) {
		index = n2n3c * i + n3c * j + k ;
		cf[index] *=  .5 ;	
	   }
	}
    }

// Normalisation des sinus (les termes sin(0*phi) et sin(np/2 *phi) ne sont pas
//			    traites)
    for (i=3; i<np+1; i += 2 )	{ 
    	for (j=0; j<n2c; j++) {
	    for (k=0; k<n3c; k++) {
		index = n2n3c * i + n3c * j + k ;
		cf[index] *=  -.5 ;	
	   }
	}
    }

// Parametres pour la routine FFT991
    int jump = 1 ;
    int inc = n2n3c ;
    int lot = 1 ;
    int isign = 1 ;

// boucle sur r et theta

    for (j=0; j<n2c; j++) {
	for (k=0; k<n3c; k++) {
	
	    index = n3c * j + k ;

// FFT inverse
	    double* debut = cf + index ;

    	    fft991_( debut, t1, trigo, facto, &inc, &jump, &np,
		     &lot, &isign) ;
	} 	// fin de la boucle sur r 
   }	// fin de la boucle sur theta

// On recopie le resultat dans ff si besoin est (c'est-a-dire si le
//  pointeur ff est different de cf) :
   
   if ( ff != cf ) {
	int n2n3f = n2f * n3f ;
	for (i=0; i<np; i++) {
    	    for (j=0; j<n2f; j++) {
	    	for (k=0; k<n3f; k++) {
		    int indexc = n2n3c * i + n3c * j + k ;
		    int indexf = n2n3f * i + n3f * j + k ;
		    ff[indexf] = cf[indexc]  ;	
		}
	    }
	}	
	
   }

    // Menage
    free (t1) ;
    
}
