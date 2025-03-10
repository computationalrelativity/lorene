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
 * Transformation en sin(2*l*theta) ou cos((2*l+1)*theta) (suivant la parite
 *  de l'indice m en phi) sur le deuxieme indice (theta)
 *  d'un tableau 3-D representant une fonction antisymetrique par rapport
 *  au plan z=0.
 *  Utilise la routine FFT Fortran FFT991
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune 
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en theta est  nt = deg[1] et doit etre de la forme
 * 			nt = 2^p 3^q 5^r + 1 
 *   int* dimf	: tableau du nombre d'elements de ff dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dimf[1] >= deg[1] = nt. 
 *
 *   double* ff : tableau des valeurs de la fonction aux nt points de
 *                        de collocation
 *
 *			  theta_l =  pi/2 l/(nt-1)       0 <= l <= nt-1 
 *
 * 			  L'espace memoire correspondant a ce
 *                        pointeur doit etre dimf[0]*dimf[1]*dimf[2] et doit 
 *			  etre alloue avant l'appel a la routine.	 
 *			  Les valeurs de la fonction doivent etre stokees
 *			  dans le tableau ff comme suit
 *		    f( theta_l ) = ff[ dimf[1]*dimf[2] * m + k + dimf[2] * l ]
 *			 ou m et k sont les indices correspondant a
 *			 phi et r respectivement.
 *	NB: cette routine suppose que la transformation en phi a deja ete
 *	    effectuee: ainsi m est un indice de Fourier, non un indice de
 *	    point de collocation en phi.
 *
 *   int* dimc	: tableau du nombre d'elements de cf dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dimc[1] >= deg[1] = nt. 
 * Sortie:
 * -------
 *   double* cf	: 	tableau des coefficients c_l de la fonction definis
 *			  comme suit (a r et phi fixes)
 *
 *			  pour m pair:
 *			   f(theta) = som_{l=0}^{nt-2} c_l sin( 2 l theta ) .
 *			  pour m impair:
 *			   f(theta) = som_{l=0}^{nt-2} c_l cos( (2 l+1) theta ) . 
 *
 * 			  L'espace memoire correspondant a ce
 *                        pointeur doit etre dimc[0]*dimc[1]*dimc[2] et doit 
 *			  etre alloue avant l'appel a la routine.	 
 *			 Le coefficient c_l (0 <= l <= nt-1) est stoke dans
 *			 le tableau cf comme suit
 *		          c_l = cf[ dimc[1]*dimc[2] * m + k + dimc[2] * l ]
 *			 ou m et k sont les indices correspondant a
 *			 phi et r respectivement.
 *			 Pour m pair, c_0 = c_{nt-1} = 0.
 *			 Pour m impair, c_{nt-1} = 0.
 *
 * NB: Si le pointeur ff est egal a cf, la routine ne travaille que sur un 
 *     seul tableau, qui constitue une entree/sortie.
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:18:03  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/15 12:48:21  j_novak
 * Corrected namespace declaration.
 *
 * Revision 1.3  2014/10/13 08:53:16  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:18:45  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/12/21 17:06:01  j_novak
 * Added all files for using fftw3.
 *
 * Revision 1.4  2003/01/31 10:31:23  e_gourgoulhon
 * Suppressed the directive #include <malloc.h> for malloc is defined
 * in <stdlib.h>
 *
 * Revision 1.3  2002/10/16 14:36:51  j_novak
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
 * Revision 2.0  1999/02/22  15:47:12  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers du C
#include <cassert>
#include <cstdlib>

// Prototypes of F77 subroutines
#include "headcpp.h"
#include "proto_f77.h"

// Prototypage des sous-routines utilisees:
namespace Lorene {
int*	facto_ini(int ) ;
double*	trigo_ini(int ) ;
double* cheb_ini(const int) ;
double* chebimp_ini(const int ) ;
//*****************************************************************************

void cftcossinsp(const int* deg, const int* dimf, double* ff, const int* dimc,
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

// Nombre de degres de liberte en theta :    
    int nt = deg[1] ;
    
// Tests de dimension:
    if (nt > n2f) {
	cout << "cftcossinsp: nt > n2f : nt = " << nt << " ,  n2f = " 
	<< n2f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (nt > n2c) {
	cout << "cftcossinsp: nt > n2c : nt = " << nt << " ,  n2c = " 
	<< n2c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n1f > n1c) {
	cout << "cftcossinsp: n1f > n1c : n1f = " << n1f << " ,  n1c = " 
	<< n1c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n3f > n3c) {
	cout << "cftcossinsp: n3f > n3c : n3f = " << n3f << " ,  n3c = " 
	<< n3c << endl ;
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
	
// Recherche de la table des points de collocations x_k = cos(theta_{nt-1-k}) :
    double* x = chebimp_ini(nt);	

    // tableau de travail G et t1
    //   (la dimension nm1+2 = nt+1 est exigee par la routine fft991)
    double* g = (double*)( malloc( (nm1+2)*sizeof(double) ) );	
    double* t1 = (double*)( malloc( (nm1+2)*sizeof(double) ) ) ;

// Parametres pour la routine FFT991
    int jump = 1 ;
    int inc = 1 ;
    int lot = 1 ;
    int isign = - 1 ;

// boucle sur phi et r (les boucles vont resp. de 0 a dimf[0]-1
//			 et 0 a dimf[2])

    int n2n3f = n2f * n3f ;
    int n2n3c = n2c * n3c ;

//=======================================================================
//				Cas m pair
//=======================================================================

    j = 0 ;
    
    while (j<n1f-1) {   //le dernier coef en phi n'est pas considere
			// (car nul)

//------------------------------------------------------------------------
//  partie cos(m phi) avec m pair : transformation en sin((2 l) theta)
//------------------------------------------------------------------------

	for (k=0; k<n3f; k++) {

	    int i0 = n2n3f * j + k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3c * j + k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau resultat

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a theta par  psi = 2 theta   et F(psi) = f(theta(psi)).  
 */
 
// Fonction G(psi) = F+(psi) sin(psi) + F_(psi) 
//---------------------------------------------
    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice (dans le tableau g) du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
// ... indice (dans le tableau ff0) du point theta correspondant a psi
		int ix = n3f * i ;
// ... indice (dans le tableau ff0) du point theta correspondant a sym(psi)
		int ixsym = n3f * isym ;
// ... F+(psi) sin(psi)
		double fps = 0.5 * ( ff0[ix] + ff0[ixsym] ) * sinp[i] ;	
// ... F_(psi) 
		double fm = 0.5 * ( ff0[ix] - ff0[ixsym] ) ; 
		g[i] = fps + fm ;
		g[isym] = fps - fm ;
    	    }
//... cas particuliers:
    	    g[0] = 0.5 * ( ff0[0] - ff0[ n3f*nm1 ] ) ;
    	    g[nm1s2] = ff0[ n3f*nm1s2 ] ;

// Developpement de G en series de Fourier par une FFT
//----------------------------------------------------

    	    F77_fft991( g, t1, trigo, facto, &inc, &jump, &nm1, &lot, &isign) ;

// Coefficients pairs du developmt. sin(2l theta) de f
//----------------------------------------------------
//  Ces coefficients sont egaux aux coefficients en sinus du developpement
//  de G en series de Fourier (le facteur -2. vient de la normalisation
//  de fft991: si fft991 donnait reellement les coefficients en sinus,
//  il faudrait le remplacer par un +1) :

	    cf0[0] = 0. ;
    	    for (i=2; i<nm1; i += 2 ) cf0[n3c*i] = -2. * g[i+1] ;	
	    cf0[n3c*nm1] = 0. ;    

// Coefficients impairs du developmt. en sin(2l theta) de f
//---------------------------------------------------------
// Ces coefficients sont obtenus par une recurrence a partir des 
// coefficients en cosinus du developpement de G en series de Fourier
// (le facteur +4. vient de la normalisation
//  de fft991: si fft991 donnait reellement les coefficients en cosinus,
//  il faudrait le remplacer par un +2.)

    	    cf0[n3c] = 2.* g[0] ;
    	    for ( i = 3; i < nt; i += 2 ) {
		cf0[ n3c*i ] = cf0[ n3c*(i-2) ] + 4. * g[i-1] ;
    	    }

	} 	// fin de la boucle sur r 

//--------------------------------------------------------------------
//  partie sin(m phi) avec m pair : transformation en sin((2 l) theta)
//--------------------------------------------------------------------

	j++ ;

	if ( (j != 1) && (j != n1f-1 ) ) {  
//  on effectue le calcul seulement dans les cas ou les coef en phi ne sont 
//  pas nuls 

	for (k=0; k<n3f; k++) {

	    int i0 = n2n3f * j + k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3c * j + k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau resultat

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a theta par  psi = 2 theta   et F(psi) = f(theta(psi)).  
 */
 
// Fonction G(psi) = F+(psi) sin(psi) + F_(psi) 
//---------------------------------------------
    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice (dans le tableau g) du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
// ... indice (dans le tableau ff0) du point theta correspondant a psi
		int ix = n3f * i ;
// ... indice (dans le tableau ff0) du point theta correspondant a sym(psi)
		int ixsym = n3f * isym ;
// ... F+(psi) sin(psi)
		double fps = 0.5 * ( ff0[ix] + ff0[ixsym] ) * sinp[i] ;	
// ... F_(psi) 
		double fm = 0.5 * ( ff0[ix] - ff0[ixsym] ) ; 
		g[i] = fps + fm ;
		g[isym] = fps - fm ;
    	    }
//... cas particuliers:
    	    g[0] = 0.5 * ( ff0[0] - ff0[ n3f*nm1 ] ) ;
    	    g[nm1s2] = ff0[ n3f*nm1s2 ] ;

// Developpement de G en series de Fourier par une FFT
//----------------------------------------------------

    	    F77_fft991( g, t1, trigo, facto, &inc, &jump, &nm1, &lot, &isign) ;

// Coefficients pairs du developmt. sin(2l theta) de f
//----------------------------------------------------
//  Ces coefficients sont egaux aux coefficients en sinus du developpement
//  de G en series de Fourier (le facteur -2. vient de la normalisation
//  de fft991: si fft991 donnait reellement les coefficients en sinus,
//  il faudrait le remplacer par un +1) :

	    cf0[0] = 0. ;
    	    for (i=2; i<nm1; i += 2 ) cf0[n3c*i] = -2. * g[i+1] ;	
	    cf0[n3c*nm1] = 0. ;    

// Coefficients impairs du developmt. en sin(2l theta) de f
//---------------------------------------------------------
// Ces coefficients sont obtenus par une recurrence a partir des 
// coefficients en cosinus du developpement de G en series de Fourier
// (le facteur +4. vient de la normalisation
//  de fft991: si fft991 donnait reellement les coefficients en cosinus,
//  il faudrait le remplacer par un +2.)

    	    cf0[n3c] = 2.* g[0] ;
    	    for ( i = 3; i < nt; i += 2 ) {
		cf0[ n3c*i ] = cf0[ n3c*(i-2) ] + 4. * g[i-1] ;
    	    }

	} 	// fin de la boucle sur r 

	}    // fin du cas ou le calcul etait necessaire (i.e. ou les
	     // coef en phi n'etaient pas nuls)

// On passe au cas m pair suivant:
	j+=3 ;

   }	// fin de la boucle sur les cas m pair

    if (n1f<=3) {	    // cas m=0 seulement (symetrie axiale)
	free (t1) ;
	free (g) ;
	return ;
    }
    
//=======================================================================
//				Cas m impair
//=======================================================================

    j = 2 ;
    
    while (j<n1f-1) {   //le dernier coef en phi n'est pas considere
			// (car nul)

//--------------------------------------------------------------------
//  partie cos(m phi) avec m impair : transformation en cos((2 l+1) theta)
//--------------------------------------------------------------------

	for (k=0; k<n3f; k++) {

	    int i0 = n2n3f * j + k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3c * j + k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau resultat

// Multiplication de la fonction par x=cos(theta) (pour la rendre paire) 
// NB: dans les commentaires qui suivent, on note h(x) = x f(x).
// (Les valeurs de h dans l'espace des configurations sont stokees dans le
//  tableau cf0).
	    for (i=0; i<nt-1; i++) cf0[n3c*i] = x[nm1-i] * ff0[n3f*i] ;
	    cf0[n3c*nm1] = 0 ;

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a theta par  psi = 2 theta   et F(psi) = f(theta(psi)).  
 */
 
// Valeur en psi=0 de la partie antisymetrique de F, F_ :
    	    double fmoins0 = 0.5 * ( cf0[0] - cf0[ n3c*nm1 ] );

// Fonction G(psi) = F+(psi) + F_(psi) sin(psi) 
//---------------------------------------------
    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice (dans le tableau g) du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
// ... indice (dans le tableau cf0) du point theta correspondant a psi
		int ix = n3c * i ;
// ... indice (dans le tableau cf0) du point theta correspondant a sym(psi)
		int ixsym = n3c * isym ;
// ... F+(psi)
		double fp = 0.5 * ( cf0[ix] + cf0[ixsym] ) ;	
// ... F_(psi) sin(psi)
		double fms = 0.5 * ( cf0[ix] - cf0[ixsym] ) * sinp[i] ; 
		g[i] = fp + fms ;
		g[isym] = fp - fms ;
    	    }
//... cas particuliers:
    	    g[0] = 0.5 * ( cf0[0] + cf0[ n3c*nm1 ] );
    	    g[nm1s2] = cf0[ n3c*nm1s2 ];

// Developpement de G en series de Fourier par une FFT
//----------------------------------------------------

    	    F77_fft991( g, t1, trigo, facto, &inc, &jump, &nm1, &lot, &isign) ;

// Coefficients pairs du developmt. de Tchebyshev de h
//----------------------------------------------------
//  Ces coefficients sont egaux aux coefficients en cosinus du developpement
//  de G en series de Fourier (le facteur 2. vient de la normalisation
//  de fft991; si fft991 donnait reellement les coef. en cosinus, il faudrait le
//   remplacer par un +1.)  :

	    cf0[0] = g[0] ;
    	    for (i=2; i<nm1; i += 2 ) cf0[n3c*i] = 2.* g[i] ;	
	    cf0[n3c*nm1] = g[nm1] ;    

// Coefficients impairs du developmt. de Tchebyshev de h
//------------------------------------------------------
// 1. Coef. c'_k (recurrence amorcee a partir de zero)
//  Le +4. en facteur de g[i] est du a la normalisation de fft991
//  (si fft991 donnait reellement les coef. en sinus, il faudrait le
//   remplacer par un -2.) 
    	    cf0[n3c] = 0 ;
    	    double som = 0;
    	    for ( i = 3; i < nt; i += 2 ) {
		cf0[n3c*i] = cf0[n3c*(i-2)] + 4. * g[i] ;
	    	som += cf0[n3c*i] ;
    	    }

// 2. Calcul de c_1 :
    	    double c1 = ( fmoins0 - som ) / nm1s2 ;

// 3. Coef. c_k avec k impair:	
    	    cf0[n3c] = c1 ;
    	    for ( i = 3; i < nt; i += 2 ) cf0[n3c*i] += c1 ;


// Coefficients de f en fonction de ceux de h
//-------------------------------------------

	    cf0[0] = 2* cf0[0] ;
	    for (i=1; i<nm1; i++) {
		cf0[n3c*i] = 2 * cf0[n3c*i] - cf0[n3c*(i-1)] ;
	    }    
	    cf0[n3c*nm1] = 0 ;

	} 	// fin de la boucle sur r 

//------------------------------------------------------------------------
//  partie sin(m phi) avec m impair : transformation en cos((2 l+1) theta)
//------------------------------------------------------------------------

	j++ ;

	if ( j != n1f-1  ) {  
//  on effectue le calcul seulement dans les cas ou les coef en phi ne sont 
//  pas nuls 

	for (k=0; k<n3f; k++) {

	    int i0 = n2n3f * j + k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3c * j + k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau resultat

// Multiplication de la fonction par x=cos(theta) (pour la rendre paire) 
// NB: dans les commentaires qui suivent, on note h(x) = x f(x).
// (Les valeurs de h dans l'espace des configurations sont stokees dans le
//  tableau cf0).
	    for (i=0; i<nt-1; i++) cf0[n3c*i] = x[nm1-i] * ff0[n3f*i] ;
	    cf0[n3c*nm1] = 0 ;

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a theta par  psi = 2 theta   et F(psi) = f(theta(psi)).  
 */
 
// Valeur en psi=0 de la partie antisymetrique de F, F_ :
    	    double fmoins0 = 0.5 * ( cf0[0] - cf0[ n3c*nm1 ] );

// Fonction G(psi) = F+(psi) + F_(psi) sin(psi) 
//---------------------------------------------
    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice (dans le tableau g) du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
// ... indice (dans le tableau cf0) du point theta correspondant a psi
		int ix = n3c * i ;
// ... indice (dans le tableau cf0) du point theta correspondant a sym(psi)
		int ixsym = n3c * isym ;
// ... F+(psi)
		double fp = 0.5 * ( cf0[ix] + cf0[ixsym] ) ;	
// ... F_(psi) sin(psi)
		double fms = 0.5 * ( cf0[ix] - cf0[ixsym] ) * sinp[i] ; 
		g[i] = fp + fms ;
		g[isym] = fp - fms ;
    	    }
//... cas particuliers:
    	    g[0] = 0.5 * ( cf0[0] + cf0[ n3c*nm1 ] );
    	    g[nm1s2] = cf0[ n3c*nm1s2 ];

// Developpement de G en series de Fourier par une FFT
//----------------------------------------------------

    	    F77_fft991( g, t1, trigo, facto, &inc, &jump, &nm1, &lot, &isign) ;

// Coefficients pairs du developmt. de Tchebyshev de h
//----------------------------------------------------
//  Ces coefficients sont egaux aux coefficients en cosinus du developpement
//  de G en series de Fourier (le facteur 2. vient de la normalisation
//  de fft991; si fft991 donnait reellement les coef. en cosinus, il faudrait le
//   remplacer par un +1.)  :

	    cf0[0] = g[0] ;
    	    for (i=2; i<nm1; i += 2 ) cf0[n3c*i] = 2.* g[i] ;	
	    cf0[n3c*nm1] = g[nm1] ;    

// Coefficients impairs du developmt. de Tchebyshev de h
//------------------------------------------------------
// 1. Coef. c'_k (recurrence amorcee a partir de zero)
//  Le +4. en facteur de g[i] est du a la normalisation de fft991
//  (si fft991 donnait reellement les coef. en sinus, il faudrait le
//   remplacer par un -2.) 
    	    cf0[n3c] = 0 ;
    	    double som = 0;
    	    for ( i = 3; i < nt; i += 2 ) {
		cf0[n3c*i] = cf0[n3c*(i-2)] + 4. * g[i] ;
	    	som += cf0[n3c*i] ;
    	    }

// 2. Calcul de c_1 :
    	    double c1 = ( fmoins0 - som ) / nm1s2 ;

// 3. Coef. c_k avec k impair:	
    	    cf0[n3c] = c1 ;
    	    for ( i = 3; i < nt; i += 2 ) cf0[n3c*i] += c1 ;


// Coefficients de f en fonction de ceux de h
//-------------------------------------------

	    cf0[0] = 2* cf0[0] ;
	    for (i=1; i<nm1; i++) {
		cf0[n3c*i] = 2 * cf0[n3c*i] - cf0[n3c*(i-1)] ;
	    }    
	    cf0[n3c*nm1] = 0 ;

	} 	// fin de la boucle sur r 

	}    // fin du cas ou le calcul etait necessaire (i.e. ou les
	     // coef en phi n'etaient pas nuls)


// On passe au cas m impair suivant:
	j+=3 ;

   }	// fin de la boucle sur les cas m impair

    // Menage
    free (t1) ;
    free (g) ;

}
}
