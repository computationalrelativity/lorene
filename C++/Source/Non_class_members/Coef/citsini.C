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


char citsini_C[] = "$Header$" ;

/*
 * Transformation inverse  sin((2*l+1)*theta)  sur le deuxieme indice (theta)
 *  d'un tableau 3-D representant une fonction symetrique par rapport
 *  au plan z=0.
 *  Utilise la routine FFT Fortran FFT991
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune 
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en theta est  nt = deg[1] et doit etre de la forme
 * 			nt = 2^p 3^q 5^r + 1 
 *   int* dimc	: tableau du nombre d'elements de cf dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dimc[1] >= deg[1] = nt. 
 *
 *   double* cf	: 	tableau des coefficients c_l de la fonction definis
 *			  comme suit (a r et phi fixes)
 *
 *			f(theta) = som_{l=0}^{nt-2} c_l sin( (2 l+1) theta ) . 
 *
 * 			  L'espace memoire correspondant a ce
 *                        pointeur doit etre dimc[0]*dimc[1]*dimc[2] et doit 
 *			  avoir ete alloue avant l'appel a la routine.	 
 *			 Le coefficient c_l (0 <= l <= nt-2) doit etre stoke dans
 *			 le tableau cf comme suit
 *		          c_l = cf[ dimc[1]*dimc[2] * j + k + dimc[2] * l ]
 *			 ou j et k sont les indices correspondant a
 *			 phi et r respectivement.
 *
 *   int* dimf	: tableau du nombre d'elements de ff dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dimf[1] >= deg[1] = nt. 
 *
 * Sortie:
 * -------
 *   double* ff : tableau des valeurs de la fonction aux nt points de
 *                        de collocation
 *
 *			  theta_l =  pi/2 l/(nt-1)       0 <= l <= nt-1 
 *
 * 			  L'espace memoire correspondant a ce
 *                        pointeur doit etre dimf[0]*dimf[1]*dimf[2] et doit 
 *			  avoir ete alloue avant l'appel a la routine.	 
 *			  Les valeurs de la fonction sont stokees
 *			  dans le tableau ff comme suit
 *		    f( theta_l ) = ff[ dimf[1]*dimf[2] * j + k + dimf[2] * l ]
 *			 ou j et k sont les indices correspondant a
 *			 phi et r respectivement.
 *
 * NB: Si le pointeur cf est egal a ff, la routine ne travaille que sur un 
 *     seul tableau, qui constitue une entree/sortie.
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  1999/02/22  15:41:16  hyc
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
double* cheb_ini(const int) ;
double* chebimp_ini(const int ) ;
extern "C" void fft991_(double [], double [], double [], int [],
                int *, int *, int *, int *, int *) ;
//*****************************************************************************

void citsini(const int* deg, const int* dimc, double* cf, const int* dimf,
		   double* ff)
{

int i, j, k ;

// Dimensions des tableaux ff et cf  :
    int n1f = dimf[0] ;
    int n2f = dimf[1] ;
    int n3f = dimf[2] ;
    int n1c = dimc[0] ;
    int n2c = dimc[1] ;
    int n3c = dimc[2] ;

// Nombres de degres de liberte en theta :    
    int nt = deg[1] ;
    
// Tests de dimension:
    if (nt > n2f) {
	cout << "citsini: nt > n2f : nt = " << nt << " ,  n2f = " 
	<< n2f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (nt > n2c) {
	cout << "citsini: nt > n2c : nt = " << nt << " ,  n2c = " 
	<< n2c << endl ;
	abort () ;
	exit(-1) ;
    }
    if ( (n1f > 1) && (n1c > n1f) ) {
	cout << "citsini: n1c > n1f : n1c = " << n1c << " ,  n1f = " 
	<< n1f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n3c > n3f) {
	cout << "citsini: n3c > n3f : n3c = " << n3c << " ,  n3f = " 
	<< n3f << endl ;
	abort () ;
	exit(-1) ;
    }

// Nombre de points pour la FFT:
    int nm1 = nt - 1;
    int nm1s2 = nm1 / 2;

// Recherche des tables pour la FFT:
    int* facto = facto_ini(nm1) ;
    double* trigo = trigo_ini(nm1) ;

// Recherche de la table des sin(psi) :
    double* sinp = cheb_ini(nt);	
	
// Recherche de la table des sin( theta_l ) :
    double* sinth = chebimp_ini(nt);	

    // tableau de travail t1 et g
    //   (la dimension nm1+2 = nt+1 est exigee par la routine fft991)
    double* g =  (double*)( malloc( (nm1+2)*sizeof(double) ) ) ;	
    double* t1 = (double*)( malloc( (nm1+2)*sizeof(double) ) ) ;

// Parametres pour la routine FFT991
    int jump = 1 ;
    int inc = 1 ;
    int lot = 1 ;
    int isign = 1 ;
    
// boucle sur phi et r 

    int n2n3f = n2f * n3f ;
    int n2n3c = n2c * n3c ;

/*   
 * Borne de la boucle sur phi: 
 *    si n1f = 1, on effectue la boucle une fois seulement.
 *    si n1f > 1, on va jusqu'a j = n1c-2 en sautant j = 1 (les coefficients
 *	j=n1c-1 et j=0 ne sont pas consideres car nuls). 
 */

    int borne_phi =  n1c-1  ;
    if (n1f == 1) borne_phi = 1 ; 

    for (j=0; j< borne_phi; j++) {
    
	if (j==1) continue ;	// on ne traite pas le terme en sin(0 phi)

	for (k=0; k<n3c; k++) {

	    int i0 = n2n3c * j + k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat

// Calcul des coefficients du developpement en cos(2 l theta)
//  de la fonction h(theta) := f(theta) sin(theta)
// en fonction de ceux de f (le resultat est stoke dans le tableau t1) :
	t1[0] = .5 * cf0[0] ;
	for (i=1; i<nm1; i++) {
	    t1[i] = .5 * ( cf0[ n3c*i ] - cf0[ n3c*(i-1) ] ) ;
	}
	t1[nm1] = -.5 * cf0[ n3c*(nt-2) ] ;    

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a theta par  psi = 2 theta  et F(psi) = h(theta(psi)).  
 */
 
// Calcul des coefficients de Fourier de la fonction 
//   G(psi) = F+(psi) + F_(psi) sin(psi)
// en fonction des coefficients en cos(2l theta) de h:

// Coefficients impairs de G
//--------------------------
 
	    double c1 = t1[1] ;

    	    double som = 0;
	    ff0[n3f] = 0 ;
    	    for ( i = 3; i < nt; i += 2 ) {
	    	ff0[ n3f*i ] = t1[i] - c1 ;
		som += ff0[ n3f*i ] ;
    	    }	

// Valeur en psi=0 de la partie antisymetrique de F, F_ :
	    double fmoins0 = nm1s2 * c1 + som ;

// Coef. impairs de G
// NB: le facteur 0.25 est du a la normalisation de fft991; si fft991
//     donnait exactement les coef. des sinus, ce facteur serait -0.5.
    	    g[1] = 0 ;
    	    for ( i = 3; i < nt; i += 2 ) {
		g[i] = 0.25 * ( ff0[ n3f*i ] - ff0[ n3f*(i-2) ] ) ;
    	    }
    	    g[nt] = 0 ;	


// Coefficients pairs de G
//------------------------
//  Ces coefficients sont egaux aux coefficients pairs du developpement de
//   h.
// NB: le facteur 0.5 est du a la normalisation de fft991; si fft991
//     donnait exactement les coef. des cosinus, ce facteur serait 1.

	    g[0] = t1[0] ;
    	    for (i=2; i<nm1; i += 2 ) g[i] = 0.5 * t1[i] ;	
    	    g[nm1] = t1[nm1] ;

// Transformation de Fourier inverse de G 
//---------------------------------------

// FFT inverse
    	    fft991_( g, t1, trigo, facto, &inc, &jump, &nm1, &lot, &isign) ;

// Valeurs de f deduites de celles de G
//-------------------------------------

    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
	
		double fp = 0.5 * ( g[i] + g[isym] ) ;
		double fm = 0.5 * ( g[i] - g[isym] ) / sinp[i] ;
		ff0[ n3f*i ] = ( fp + fm ) / sinth[i] ;
		ff0[ n3f*isym ] = ( fp - fm ) / sinth[isym] ;
    	    }
	
//... cas particuliers:
	    ff0[0] = 0 ;
	    ff0[ n3f*nm1 ] = g[0] - fmoins0 ;
	    ff0[ n3f*nm1s2 ] = g[nm1s2] / sinth[nm1s2];


	} 	// fin de la boucle sur r 
   }	// fin de la boucle sur phi

    // Menage
    free (t1) ;
    free (g) ;
    
}
