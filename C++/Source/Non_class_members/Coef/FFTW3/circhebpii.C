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


char circhebpii_C[] = "$Header$" ;


/*
 * Transformation de Tchebyshev inverse (cas rare) sur le troisieme indice 
 * (indice correspondant a r) d'un tableau 3-D decrivant une fonction quelconque. 
 *  Utilise la bibliotheque fftw.
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune 
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en r est  nr = deg[2] et doit etre de la forme
 * 			nr = 2*p + 1 
 *   int* dimc	: tableau du nombre d'elements de cf dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dimc[2] >= deg[2] = nr. 
 *		  NB: pour dimc[0] = 1 (un seul point en phi), la transformation
 *		      est bien effectuee.
 *		      pour dimc[0] > 1 (plus d'un point en phi), la
 *		      transformation n'est effectuee que pour les indices (en phi)
 *		      j != 1 et j != dimc[0]-1 (cf. commentaires sur borne_phi).
 *
 *   double* cf	:   tableau des coefficients c_i de la fonction definis
 *		    comme suit (a theta et phi fixes)
 *
 *		   Si l impair f(x) = som_{i=0}^{nr-1} c_i T_{2i}(x) , 
 *                 Si l pair f(x) = som_{i=0}^{nr-1} c_i T_{2i+1}(x) ,
 *
 *		    ou T_{i}(x) designe le polynome de Tchebyshev de degre i. 	 
 *		    Les coefficients c_i (0 <= i <= nr-1) doivent etre stokes 
 *		    dans le tableau cf comme suit
 *			   c_i = cf[ dimc[1]*dimc[2] * j + dimc[2] * k + i ]
 *		    ou j et k sont les indices correspondant a phi et theta 
 *		    respectivement.
 * 		    L'espace memoire correspondant a ce pointeur doit etre 
 *		    dimc[0]*dimc[1]*dimc[2] et doit etre alloue avant l'appel a
 *		    la routine.	 
 *
 *   int* dimf	: tableau du nombre d'elements de ff dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dimf[2] >= deg[2] = nr. 
 *
 * Sortie:
 * -------
 *   double* ff : tableau des valeurs de la fonction aux nr points de
 *                        de collocation
 *
 *			  x_i = sin( pi/2 i/(nr-1) )      0 <= i <= nr-1 
 *
 *		    Les valeurs de la fonction sont stokees dans le 
 *		    tableau ff comme suit
 *			 f( x_i ) = ff[ dimf[1]*dimf[2] * j + dimf[2] * k + i ]
 *		    ou j et k sont les indices correspondant a phi et theta 
 *		    respectivement.
 * 		    L'espace memoire correspondant a ce pointeur doit etre 
 *		    dimf[0]*dimf[1]*dimf[2] et doit avoir ete alloue avant 
 *		    l'appel a la routine.	 
 *
 * NB: Si le pointeur cf est egal a ff, la routine ne travaille que sur un 
 *     seul tableau, qui constitue une entree/sortie.
 */

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:20  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:18:49  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/12/21 17:06:02  j_novak
 * Added all files for using fftw3.
 *
 * Revision 1.1  2004/11/23 15:13:50  m_forot
 * Added the bases for the cases without any equatorial symmetry
 * (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 *
 * $Header$
 *
 */

// headers du C
#include <cstdlib>
#include <fftw3.h>

//Lorene prototypes
#include "tbl.h"

// Prototypage des sous-routines utilisees:
namespace Lorene {
fftw_plan back_fft(int, Tbl*&) ;
double* cheb_ini(const int) ;
double* chebimp_ini(const int ) ;

//*****************************************************************************

void circhebpii(const int* deg, const int* dimc, double* cf, 
		    const int* dimf, double* ff)

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
    if (nr > n3c) {
	cout << "circhebpii: nr > n3c : nr = " << nr << " ,  n3c = " 
	<< n3c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (nr > n3f) {
	cout << "circhebpii: nr > n3f : nr = " << nr << " ,  n3f = " 
	<< n3f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n1c > n1f) {
	cout << "circhebpii: n1c > n1f : n1c = " << n1c << " ,  n1f = " 
	<< n1f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n2c > n2f) {
	cout << "circhebpii: n2c > n2f : n2c = " << n2c << " ,  n2f = " 
	<< n2f << endl ;
	abort () ;
	exit(-1) ;
    }

// Nombre de points pour la FFT:
    int nm1 = nr - 1;
    int nm1s2 = nm1 / 2;

// Recherche des tables pour la FFT:
    Tbl* pg = 0x0 ;
    fftw_plan p = back_fft(nm1, pg) ;
    Tbl& g = *pg ;
    double* t1 = new double[nr] ;

// Recherche de la table des sin(psi) :
    double* sinp = cheb_ini(nr);

// Recherche de la table des points de collocations x_k :
    double* x = chebimp_ini(nr);	
	
// boucle sur phi et theta

    int n2n3f = n2f * n3f ;
    int n2n3c = n2c * n3c ;

/*   
 * Borne de la boucle sur phi: 
 *    si n1c = 1, on effectue la boucle une fois seulement.
 *    si n1c > 1, on va jusqu'a j = n1c-2 en sautant j = 1 (les coefficients
 *	j=n1c-1 et j=0 ne sont pas consideres car nuls). 
 */
    int borne_phi = ( n1c > 1 ) ? n1c-1 : 1 ;

    for (j=0; j< borne_phi; j++) {
    
	if (j==1) continue ;	// on ne traite pas le terme en sin(0 phi)


	/************ Cas l impair **********/

	for (k=1; k<n2c; k+=2) {

	    int i0 = n2n3c * j + n3c * k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + n3f * k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat


// Calcul des coefficients de Fourier de la fonction 
//   G(psi) = F+(theta) + F_(theta) sin(theta)
// en fonction des coefficients de Tchebyshev de f:

// Coefficients impairs de G
//--------------------------
 
	    double c1 = cf0[1] ;

    	    double som = 0;
	    ff0[1] = 0 ;
    	    for ( i = 3; i < nr; i += 2 ) {
	    	ff0[i] = cf0[i] - c1 ;
		som += ff0[i] ;
    	    }	

// Valeur en theta=0 de la partie antisymetrique de F, F_ :
	    double fmoins0 = nm1s2 * c1 + som ;

// Coef. impairs de G
// NB: le facteur 0.25 est du a la normalisation de fftw; si fftw
//     donnait exactement les coef. des sinus, ce facteur serait -0.5.
    	    for ( i = 3; i < nr; i += 2 ) {
		g.set(nm1-i/2) = 0.25 * ( ff0[i] - ff0[i-2] ) ;
    	    }


// Coefficients pairs de G
//------------------------
//  Ces coefficients sont egaux aux coefficients pairs du developpement de
//   f.
// NB: le facteur 0.5 est du a la normalisation de fftw; si fftw
//     donnait exactement les coef. des cosinus, ce facteur serait 1.

	    g.set(0) = cf0[0] ;
    	    for (i=1; i<nm1s2; i++) g.set(i) = 0.5 * cf0[2*i] ;	
    	    g.set(nm1s2) = cf0[nm1] ;

// Transformation de Fourier inverse de G 
//---------------------------------------

// FFT inverse
	    fftw_execute(p) ;

// Valeurs de f deduites de celles de G
//-------------------------------------

    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice (dans le tableau g) du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
// ... indice (dans le tableau ff0) du point x correspondant a psi
		int ix = nm1 - i ;
// ... indice (dans le tableau ff0) du point x correspondant a sym(psi)
		int ixsym = nm1 -  isym ;

		double fp = .5 * ( g(i) + g(isym) ) ;
		double fm = .5 * ( g(i) - g(isym) ) / sinp[i] ;

		ff0[ix] = fp + fm ;
		ff0[ixsym] = fp - fm ;
    	    }
	
//... cas particuliers:
	    ff0[0] = g(0) - fmoins0 ;
	    ff0[nm1] = g(0) + fmoins0 ;
	    ff0[nm1s2] = g(nm1s2) ;
	} 	// fin de la boucle sur theta 
	
	/*********** Cas l pair **********/
	
	for (k=0; k<n2c; k+=2) {

	    int i0 = n2n3c * j + n3c * k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + n3f * k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat

// Calcul des coefficients du developpement en T_{2i}(x) de la fonction
//  h(x) := x f(x) a partir des coefficients de f (resultat stoke dans le
//  tableau t1 :
	t1[0] = .5 * cf0[0] ;
	for (i=1; i<nm1; i++) t1[i] = .5 * ( cf0[i] + cf0[i-1] ) ;
	t1[nm1] = .5 * cf0[nr-2] ;


// Calcul des coefficients de Fourier de la fonction 
//   G(psi) = F+(theta) + F_(theta) sin(theta)
// en fonction des coefficients de Tchebyshev de f:

// Coefficients impairs de G
//--------------------------
 
	    double c1 = t1[1] ;

    	    double som = 0;
	    ff0[1] = 0 ;
    	    for ( i = 3; i < nr; i += 2 ) {
	    	ff0[i] = t1[i] - c1 ;
		som += ff0[i] ;
    	    }	

// Valeur en theta=0 de la partie antisymetrique de F, F_ :
	    double fmoins0 = nm1s2 * c1 + som ;

// Coef. impairs de G
// NB: le facteur 0.25 est du a la normalisation de fftw; si fftw
//     donnait exactement les coef. des sinus, ce facteur serait -0.5.
    	    for ( i = 3; i < nr; i += 2 ) {
		g.set(nm1-i/2) = 0.25 * ( ff0[i] - ff0[i-2] ) ;
    	    }


// Coefficients pairs de G
//------------------------
//  Ces coefficients sont egaux aux coefficients pairs du developpement de
//   f.
// NB: le facteur 0.5 est du a la normalisation de fftw; si fftw
//     donnait exactement les coef. des cosinus, ce facteur serait 1.

	    g.set(0) = t1[0] ;
    	    for (i=1; i<nm1s2; i++ ) g.set(i) = 0.5 * t1[2*i] ;	
    	    g.set(nm1s2) = t1[nm1] ;

// Transformation de Fourier inverse de G 
//---------------------------------------

// FFT inverse
	    fftw_execute(p) ;

// Valeurs de f deduites de celles de G
//-------------------------------------

    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice (dans le tableau g) du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
// ... indice (dans le tableau ff0) du point x correspondant a psi
		int ix = nm1 - i ;
// ... indice (dans le tableau ff0) du point x correspondant a sym(psi)
		int ixsym = nm1 -  isym ;

		double fp = .5 * ( g(i) + g(isym) ) ;
		double fm = .5 * ( g(i) - g(isym) ) / sinp[i] ;

		ff0[ix] = ( fp + fm ) / x[ix];
		ff0[ixsym] = ( fp - fm ) / x[ixsym] ;
    	    }
	
//... cas particuliers:
	    ff0[0] = 0 ;
	    ff0[nm1] = g(0) + fmoins0 ;
	    ff0[nm1s2] = g(nm1s2) / x[nm1s2] ;

	} 	// fin de la boucle sur theta 


   }	// fin de la boucle sur phi

    delete [] t1 ;
}
}
