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


char cfpcossin_C[] = "$Header$" ;

/*
 * Transformation de Fourier sur le premier indice d'un tableau 3-D
 *  par le biais de la routine FFT Fortran FFT991
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune 
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en phi est  np = deg[0] et doit etre de la forme
 * 			np = 2^p 3^q 5^r 
 *   int* dim	: tableau du nombre d'elements de ff dans chacune des trois
 *	          dimensions.
 *		  On doit avoir  dim[0] >= deg[0] + 2  = np + 2. 
 *
 * Entree/Sortie : 
 * ---------------
 *   double* cf : entree: tableau des valeurs de la fonction f aux np points de
 *                        de collocation 
 *				phi_k = T k/np      0 <= k <= np-1
 *			  T etant la periode de f. La convention de stokage
 *			  utilisee est la suivante
 *	cf[ dim[2]*dim[1]*k + dim[2]*j + i ] = f(phi_k)    0 <= k <= np-1,
 *			 les indices j et i correspondant respectivement
 *			  a theta et a r. 
 *			  L'espace memoire correspondant au
 *                        pointeur cf doit etre dim[0]*dim[1]*dim[2] et doit 
 *			  etre alloue avant l'appel a la routine.	
 * 
 * 		  sortie: tableau des coefficients de la fonction suivant	 
 *			  la convention de stokage
 *		cf[ dim[2]*dim[1]*k + dim[2]*j + i ] = c_k      0<= k <= np,
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
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/10/16 14:36:43  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/09/09 13:00:39  e_gourgoulhon
 * Modification of declaration of Fortran 77 prototypes for
 * a better portability (in particular on IBM AIX systems):
 * All Fortran subroutine names are now written F77_* and are
 * defined in the new file C++/Include/proto_f77.h.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/02/22  15:48:58  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */
// headers du C
#include <assert.h>
#include <malloc.h>
#include <stdlib.h>

// Prototypes of F77 subroutines
#include "headcpp.h"
#include "proto_f77.h"

// Prototypage des sous-routines utilisees:
int*	facto_ini(int ) ;
double*	trigo_ini(int ) ;
//*****************************************************************************

void cfpcossin(const int* deg, const int* dim, double* cf)

{

int i, j, k, index ;

// Dimensions du tableau cf :
    int n1 = dim[0] ;
    int n2 = dim[1] ;
    int n3 = dim[2] ;

// Nombres de degres de liberte en phi :    
    int np = deg[0] ;
    
// Tests de dimension:
    if (np+2 > n1) {
	cout << "cfpcossin: np+2 > n1 : np+2 = " << np+2 << " ,  n1 = " 
	<< n1 << endl ;
	abort () ;
	exit(-1) ;
    }

// Recherche des tables

    int* facto = facto_ini(np) ;
    double* trigo = trigo_ini(np) ;
	
// Tableau de travail
    double* t1 = (double*)( malloc( (np+2)*sizeof(double) ) ) ;

// Parametres pour la routine FFT991
    int jump = 1 ;
    int inc = n2 * n3 ;
    int lot = 1 ;
    int isign = - 1 ;

// boucle sur theta et r  
    for (j=0; j<n2; j++) {
	for (k=0; k<n3; k++) {
	    index = n3 * j + k ;
// FFT
	    double* debut = cf + index ;

    	    F77_fft991( debut, t1, trigo, facto, &inc, &jump, &np,
		     &lot, &isign) ;
	} 	// fin de la boucle sur r 
   }	// fin de la boucle sur theta

// Normalisation des cosinus
    int n2n3 = n2 * n3 ;
    for (i=2; i<np; i += 2 )	{
    	for (j=0; j<n2; j++) {
	    for (k=0; k<n3; k++) {
		index = n2n3 * i + n3 * j + k ;
		cf[index] *=  2. ;	
	   }
	}
    }

// Normalisation des sinus
    for (i=1; i<np+1; i += 2 )	{
    	for (j=0; j<n2; j++) {
	    for (k=0; k<n3; k++) {
		index = n2n3 * i + n3 * j + k ;
		cf[index] *=  -2. ;	
	   }
	}
    }

    // Menage
    free(t1) ;

}
