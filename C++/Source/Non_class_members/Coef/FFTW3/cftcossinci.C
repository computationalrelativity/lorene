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


char cftcossinci_C[] = "$Header$" ;

/*
 * Transformation en  cos((2*l+1)*theta) ou sin(2*l*theta)(suivant la parite
 *  de l'indice m en phi) sur le deuxieme indice (theta)
 *  d'un tableau 3-D representant une fonction antisymetrique par rapport
 *  au plan z=0.
 *  Utilise la bibliotheque fftw.
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune 
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en theta est  nt = deg[1] et doit etre de la forme
 * 			nt = 2*p + 1 
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
 *			   f(theta) = som_{l=0}^{nt-2} c_l cos( (2 l+1) theta ) . 
 *			  pour m impair:
 *			   f(theta) = som_{l=0}^{nt-2} c_l sin( 2 l theta ) . 
 *
 * 			  L'espace memoire correspondant a ce
 *                        pointeur doit etre dimc[0]*dimc[1]*dimc[2] et doit 
 *			  etre alloue avant l'appel a la routine.	 
 *			 Le coefficient c_l (0 <= l <= nt-1) est stoke dans
 *			 le tableau cf comme suit
 *		          c_l = cf[ dimc[1]*dimc[2] * m + k + dimc[2] * l ]
 *			 ou m et k sont les indices correspondant a
 *			 phi et r respectivement.
 *			 Pour m pair, c_{nt-1} = 0.
 *			 Pour m impair, c_0 = c_{nt-1} = 0 .
 *
 * NB: Si le pointeur ff est egal a cf, la routine ne travaille que sur un 
 *     seul tableau, qui constitue une entree/sortie.
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:18:48  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/12/21 17:06:02  j_novak
 * Added all files for using fftw3.
 *
 * Revision 1.4  2003/01/31 10:31:23  e_gourgoulhon
 * Suppressed the directive #include <malloc.h> for malloc is defined
 * in <stdlib.h>
 *
 * Revision 1.3  2002/10/16 14:36:46  j_novak
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
 * Revision 2.0  1999/02/22  15:47:40  hyc
 * *** empty log message ***
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
fftw_plan prepare_fft(int, Tbl*&) ;
double* cheb_ini(const int) ;
double* chebimp_ini(const int ) ;
//*****************************************************************************

void cftcossinci(const int* deg, const int* dimf, double* ff, const int* dimc,
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
	cout << "cftcossinci: nt > n2f : nt = " << nt << " ,  n2f = " 
	<< n2f << endl ;
	abort () ;
    }
    if (nt > n2c) {
	cout << "cftcossinci: nt > n2c : nt = " << nt << " ,  n2c = " 
	<< n2c << endl ;
	abort () ;
    }
    if (n1f > n1c) {
	cout << "cftcossinci: n1f > n1c : n1f = " << n1f << " ,  n1c = " 
	<< n1c << endl ;
	abort () ;
    }
    if (n3f > n3c) {
	cout << "cftcossinci: n3f > n3c : n3f = " << n3f << " ,  n3c = " 
	<< n3c << endl ;
	abort () ;
    }

// Nombre de points pour la FFT:
    int nm1 = nt - 1;
    int nm1s2 = nm1 / 2;

// Recherche des tables pour la FFT:
    Tbl* pg = 0x0 ;
    fftw_plan p = prepare_fft(nm1, pg) ;
    Tbl& g = *pg ;

// Recherche de la table des sin(psi) :
    double* sinp = cheb_ini(nt);	
	
// Recherche de la table des points de collocations x_k = cos(theta_{nt-1-k}) :
    double* x = chebimp_ini(nt);	

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

//--------------------------------------------------------------------
//  partie cos(m phi) avec m pair : transformation en cos((2 l+1) theta)
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
		g.set(i) = fp + fms ;
		g.set(isym) = fp - fms ;
    	    }
//... cas particuliers:
    	    g.set(0) = 0.5 * ( cf0[0] + cf0[ n3c*nm1 ] );
    	    g.set(nm1s2) = cf0[ n3c*nm1s2 ];

// Developpement de G en series de Fourier par une FFT
//----------------------------------------------------

	    fftw_execute(p) ;

// Coefficients pairs du developmt. de Tchebyshev de h
//----------------------------------------------------
//  Ces coefficients sont egaux aux coefficients en cosinus du developpement
//  de G en series de Fourier (le facteur 2/nm1 vient de la normalisation
//  de fftw; si fftw donnait reellement les coef. en cosinus, il faudrait le
//   remplacer par un +1.)  :

	    double fac = 2./double(nm1) ;
	    cf0[0] = g(0) / double(nm1) ;
    	    for (i=2; i<nm1; i += 2 ) cf0[n3c*i] = fac * g(i/2) ;	
	    cf0[n3c*nm1] = g(nm1s2) ;    

// Coefficients impairs du developmt. de Tchebyshev de h
//------------------------------------------------------
// 1. Coef. c'_k (recurrence amorcee a partir de zero)
//  Le 4/nm1 en facteur de g(i) est du a la normalisation de fftw
//  (si fftw donnait reellement les coef. en sinus, il faudrait le
//   remplacer par un -2.) 
	    fac *= 2. ;
    	    cf0[n3c] = 0 ;
    	    double som = 0;
    	    for ( i = 3; i < nt; i += 2 ) {
		cf0[n3c*i] = cf0[n3c*(i-2)] + fac * g(nm1 - i/2) ;
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


//--------------------------------------------------------------------
//  partie sin(m phi) avec m pair : transformation en cos((2 l+1) theta)
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
		g.set(i) = fp + fms ;
		g.set(isym) = fp - fms ;
    	    }
//... cas particuliers:
    	    g.set(0) = 0.5 * ( cf0[0] + cf0[ n3c*nm1 ] );
    	    g.set(nm1s2) = cf0[ n3c*nm1s2 ];

// Developpement de G en series de Fourier par une FFT
//----------------------------------------------------

	    fftw_execute(p) ;

// Coefficients pairs du developmt. de Tchebyshev de h
//----------------------------------------------------
//  Ces coefficients sont egaux aux coefficients en cosinus du developpement
//  de G en series de Fourier (le facteur 2/nm1 vient de la normalisation
//  de fftw; si fftw donnait reellement les coef. en cosinus, il faudrait le
//   remplacer par un +1.)  :

	    double fac = 2./double(nm1) ;
	    cf0[0] = g(0) / double(nm1) ;
    	    for (i=2; i<nm1; i += 2 ) cf0[n3c*i] = fac * g(i/2) ;	
	    cf0[n3c*nm1] = g(nm1s2) ;    

// Coefficients impairs du developmt. de Tchebyshev de h
//------------------------------------------------------
// 1. Coef. c'_k (recurrence amorcee a partir de zero)
//  Le 4/nm1 en facteur de g(i) est du a la normalisation de fftw
//  (si fftw donnait reellement les coef. en sinus, il faudrait le
//   remplacer par un -2.) 
	    fac *= 2. ;
    	    cf0[n3c] = 0 ;
    	    double som = 0;
    	    for ( i = 3; i < nt; i += 2 ) {
		cf0[n3c*i] = cf0[n3c*(i-2)] + fac * g(nm1 - i/2) ;
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

// On passe au cas m pair suivant:
	j+=3 ;

   }	// fin de la boucle sur les cas m pair

    if (n1f<=3) 	    // cas m=0 seulement (symetrie axiale)
	return ;
    
//=======================================================================
//				Cas m impair
//=======================================================================

    j = 2 ;
    
    while (j<n1f-1) {   //le dernier coef en phi n'est pas considere
			// (car nul)

//--------------------------------------------------------------------
//  partie cos(m phi) avec m impair : transformation en sin((2 l) theta)
//--------------------------------------------------------------------

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
		g.set(i) = fps + fm ;
		g.set(isym) = fps - fm ;
    	    }
//... cas particuliers:
    	    g.set(0) = 0.5 * ( ff0[0] - ff0[ n3f*nm1 ] ) ;
    	    g.set(nm1s2) = ff0[ n3f*nm1s2 ] ;

// Developpement de G en series de Fourier par une FFT
//----------------------------------------------------
	    fftw_execute(p) ;

// Coefficients pairs du developmt. sin(2l theta) de f
//----------------------------------------------------
//  Ces coefficients sont egaux aux coefficients en sinus du developpement
//  de G en series de Fourier (le facteur -2/nm1 vient de la normalisation
//  de fftw: si fftw donnait reellement les coefficients en sinus,
//  il faudrait le remplacer par un +1) :
	    double fac = 2./double(nm1) ;
	    cf0[0] = 0. ;
    	    for (i=2; i<nm1; i += 2 ) cf0[n3c*i] = -fac * g(nm1 - i/2); 
	    cf0[n3c*nm1] = 0. ;    

// Coefficients impairs du developmt. en sin(2l theta) de f
//---------------------------------------------------------
// Ces coefficients sont obtenus par une recurrence a partir des 
// coefficients en cosinus du developpement de G en series de Fourier
// (le facteur 4/nm1 vient de la normalisation
//  de fftw: si fftw donnait reellement les coefficients en cosinus,
//  il faudrait le remplacer par un +2.)

    	    cf0[n3c] = fac * g(0) ;
	    fac *= 2. ;
    	    for ( i = 3; i < nt; i += 2 ) {
		cf0[ n3c*i ] = cf0[ n3c*(i-2) ] + fac * g(i/2) ;
    	    }

	} 	// fin de la boucle sur r 

//------------------------------------------------------------------------
//  partie sin(m phi) avec m impair : transformation en sin((2 l) theta)
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
		g.set(i) = fps + fm ;
		g.set(isym) = fps - fm ;
    	    }
//... cas particuliers:
    	    g.set(0) = 0.5 * ( ff0[0] - ff0[ n3f*nm1 ] ) ;
    	    g.set(nm1s2) = ff0[ n3f*nm1s2 ] ;

// Developpement de G en series de Fourier par une FFT
//----------------------------------------------------

	    fftw_execute(p) ;

// Coefficients pairs du developmt. sin(2l theta) de f
//----------------------------------------------------
//  Ces coefficients sont egaux aux coefficients en sinus du developpement
//  de G en series de Fourier (le facteur -2/nm1 vient de la normalisation
//  de fftw: si fftw donnait reellement les coefficients en sinus,
//  il faudrait le remplacer par un +1) :
	    double fac = 2./double(nm1) ;
	    cf0[0] = 0. ;
    	    for (i=2; i<nm1; i += 2 ) cf0[n3c*i] = -fac * g(nm1 - i/2); 
	    cf0[n3c*nm1] = 0. ;    

// Coefficients impairs du developmt. en sin(2l theta) de f
//---------------------------------------------------------
// Ces coefficients sont obtenus par une recurrence a partir des 
// coefficients en cosinus du developpement de G en series de Fourier
// (le facteur 4/nm1 vient de la normalisation
//  de fftw: si fftw donnait reellement les coefficients en cosinus,
//  il faudrait le remplacer par un +2.)

    	    cf0[n3c] = fac * g(0) ;
	    fac *= 2. ;
    	    for ( i = 3; i < nt; i += 2 ) {
		cf0[ n3c*i ] = cf0[ n3c*(i-2) ] + fac * g(i/2) ;
    	    }

	} 	// fin de la boucle sur r 

	}    // fin du cas ou le calcul etait necessaire (i.e. ou les
	     // coef en phi n'etaient pas nuls)


// On passe au cas m impair suivant:
	j+=3 ;

   }	// fin de la boucle sur les cas m impair


}
