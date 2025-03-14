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


 

/*
 * Transformation de Tchebyshev inverse T_{2k+1}/T_{2k} (suivant la parite de 
 * l'indice m en phi) sur le troisieme indice 
 * (indice correspondant a r) d'un tableau 3-D decrivant par exemple
 * d/dr d'une fonction symetrique
 * par rapport au plan equatorial z = 0 et sans aucune autre symetrie,  
 *  cad que l'on a effectue
 *	1/ un developpement en polynomes de Tchebyshev impairs pour m pair 
 *	2/ un developpement en polynomes de Tchebyshev pairs pour m impair 
 *
 *
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
 *
 *   double* cf	:   tableau des coefficients c_i de la fonction definis
 *		    comme suit (a theta et phi fixes)
 *
 *		    -- pour m pair (i.e. j = 0, 1, 4, 5, 8, 9, ...) :
 *
 *			  f(x) = som_{i=0}^{nr-2} c_i T_{2i+1}(x)
 * 
 *		       ou T_{2i+1}(x) designe le polynome de Tchebyshev de 
 *		       degre 2i+1. 
 *
 *		    -- pour m impair (i.e. j = 2, 3, 6, 7, 10, 11, ...) :
 *
 *			    f(x) =  som_{i=0}^{nr-1} c_i T_{2i}(x) , 
 *
 *		       ou T_{2i}(x) designe le polynome de Tchebyshev de 
 *		       degre 2i. 
 *	 
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
 * Revision 1.5  2016/12/05 16:18:05  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:20  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:18:49  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2004/12/27 14:27:28  j_novak
 * Added forgotten "delete [] t1"
 *
 * Revision 1.1  2004/12/21 17:06:02  j_novak
 * Added all files for using fftw3.
 *
 * Revision 1.4  2003/01/31 10:31:23  e_gourgoulhon
 * Suppressed the directive #include <malloc.h> for malloc is defined
 * in <stdlib.h>
 *
 * Revision 1.3  2002/10/16 14:36:53  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/09/09 13:00:40  e_gourgoulhon
 * Modification of declaration of Fortran 77 prototypes for
 * a better portability (in particular on IBM AIX systems):
 * All Fortran subroutine names are now written F77_* and are
 * defined in the new file C++/Include/proto_f77.h.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/02/22  15:43:19  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */



// headers du C
#include <cassert>
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

void circhebpimi(const int* deg, const int* dimc, double* cf,
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
	cout << "circhebpimi: nr > n3c : nr = " << nr << " ,  n3c = " 
	<< n3c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (nr > n3f) {
	cout << "circhebpimi: nr > n3f : nr = " << nr << " ,  n3f = " 
	<< n3f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n1c > n1f) {
	cout << "circhebpimi: n1c > n1f : n1c = " << n1c << " ,  n1f = " 
	<< n1f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n2c > n2f) {
	cout << "circhebpimi: n2c > n2f : n2c = " << n2c << " ,  n2f = " 
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

//=======================================================================
//				Cas m pair
//=======================================================================

    j = 0 ;
    
    while (j<n1f-1) {   //le dernier coef en phi n'est pas considere
			// (car nul)


//------------------------------------------------------------------------
//  partie cos(m phi) avec m pair : developpement en T_{2i+1}(x)
//------------------------------------------------------------------------

	for (k=0; k<n2c; k++) {

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

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a x par  x = cos(psi/2)   et F(psi) = h(x(psi)).  
 */

// Calcul des coefficients de Fourier de la fonction 
//   G(psi) = F+(psi) + F_(psi) sin(psi)
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

// Valeur en psi=0 de la partie antisymetrique de F, F_ :
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

//------------------------------------------------------------------------
//  partie sin(m phi) avec m pair : developpement en T_{2i+1}(x)
//------------------------------------------------------------------------

	j++ ;

	if ( (j != 1) && (j != n1f-1) ) {  
//  on effectue le calcul seulement dans les cas ou les coef en phi ne sont 
//  pas nuls 

	for (k=0; k<n2c; k++) {

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

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a x par  x = cos(psi/2)   et F(psi) = h(x(psi)).  
 */

// Calcul des coefficients de Fourier de la fonction 
//   G(psi) = F+(psi) + F_(psi) sin(psi)
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

// Valeur en psi=0 de la partie antisymetrique de F, F_ :
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

	}    // fin du cas ou le calcul etait necessaire (i.e. ou les
	     // coef en phi n'etaient pas nuls)

// On passe au cas m pair suivant:
	j+=3 ;

   }	// fin de la boucle sur les cas m pair

    if (n1f<=3) {	    // cas m=0 seulement (symetrie axiale)
	delete [] t1 ;
	return ;
    }

//=======================================================================
//				Cas m impair
//=======================================================================

    j = 2 ;
    
    while (j<n1f-1) {   
    
//--------------------------------------------------------------------
//  partie cos(m phi) avec m impair : developpement en T_{2i}(x)
//--------------------------------------------------------------------

	for (k=0; k<n2c; k++) {

	    int i0 = n2n3c * j + n3c * k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + n3f * k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a x par  x = cos(psi/2)   et F(psi) = f(x(psi)).  
 */

// Calcul des coefficients de Fourier de la fonction 
//   G(psi) = F+(psi) + F_(psi) sin(psi)
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

// Valeur en psi=0 de la partie antisymetrique de F, F_ :
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


//--------------------------------------------------------------------
//  partie sin(m phi) avec m impair : developpement en T_{2i}(x)
//--------------------------------------------------------------------

	j++ ;

	if ( j != n1f-1 ) {  
//  on effectue le calcul seulement dans les cas ou les coef en phi ne sont 
//  pas nuls 

	for (k=0; k<n2c; k++) {

	    int i0 = n2n3c * j + n3c * k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + n3f * k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a x par  x = cos(psi/2)   et F(psi) = f(x(psi)).  
 */

// Calcul des coefficients de Fourier de la fonction 
//   G(psi) = F+(psi) + F_(psi) sin(psi)
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

// Valeur en psi=0 de la partie antisymetrique de F, F_ :
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

	}    // fin du cas ou le calcul etait necessaire (i.e. ou les
	     // coef en phi n'etaient pas nuls)

// On passe au cas m impair suivant:
	j+=3 ;

   }	// fin de la boucle sur les cas m impair

    delete [] t1 ;
}

}
