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


char citcossinci_C[] = "$Header$" ;

/*
 * Transformation inverse cos((2*l+1)*theta) ou sin(2*l*theta) (suivant la 
 *  parite de l'indice m en phi) sur le deuxieme indice (theta)
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
 *   int* dimc	: tableau du nombre d'elements de cf dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dimc[1] >= deg[1] = nt. 
 *
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
 *			  avoir ete alloue avant l'appel a la routine.	 
 *			 Le coefficient c_l (0 <= l <= nt-1) doit etre stoke dans
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
 * Revision 1.2  2004/12/27 14:27:28  j_novak
 * Added forgotten "delete [] t1"
 *
 * Revision 1.1  2004/12/21 17:06:03  j_novak
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
 * Revision 2.0  1999/02/22  15:42:35  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers du C
#include <stdlib.h>
#include <fftw3.h>

//Lorene prototypes
#include "tbl.h"

// Prototypage des sous-routines utilisees:
fftw_plan back_fft(int, Tbl*&) ;
double* cheb_ini(const int) ;
double* chebimp_ini(const int ) ;
//*****************************************************************************

void citcossinci(const int* deg, const int* dimc, double* cf, const int* dimf,
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
	cout << "citcossinci: nt > n2f : nt = " << nt << " ,  n2f = " 
	<< n2f << endl ;
	abort () ;
    }
    if (nt > n2c) {
	cout << "citcossinci: nt > n2c : nt = " << nt << " ,  n2c = " 
	<< n2c << endl ;
	abort () ;
    }
    if (n1c > n1f) {
	cout << "citcossinci: n1c > n1f : n1c = " << n1c << " ,  n1f = " 
	<< n1f << endl ;
	abort () ;
    }
    if (n3c > n3f) {
	cout << "citcossinci: n3c > n3f : n3c = " << n3c << " ,  n3f = " 
	<< n3f << endl ;
	abort () ;
    }

// Nombre de points pour la FFT:
    int nm1 = nt - 1;
    int nm1s2 = nm1 / 2;

// Recherche des tables pour la FFT:
    Tbl* pg = 0x0 ;
    fftw_plan p = back_fft(nm1, pg) ;
    Tbl& g = *pg ;
    double* t1 = new double[nt] ;

// Recherche de la table des sin(psi) :
    double* sinp = cheb_ini(nt);	
	
// Recherche de la table des points de collocations x_k = cos(theta_{nt-1-k}):
    double* x = chebimp_ini(nt) ;	

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

//-----------------------------------------------------------------------
//  partie cos(m phi) avec m pair : transformation cos((2 l+1) theta) inverse
//-----------------------------------------------------------------------

	for (k=0; k<n3c; k++) {

	    int i0 = n2n3c * j + k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat
	     
// Calcul des coefficients du developpement en T_{2i}(x) de la fonction
//  h(x) := x f(x)  (x=cos(theta)) a partir des coefficients de f 
// (resultat stoke dans le tableau t1 :
	t1[0] = .5 * cf0[0] ; 
	for (i=1; i<nm1; i++) t1[i] = .5 * ( cf0[n3c*i] + cf0[n3c*(i-1)] ) ;
	t1[nm1] = .5 * cf0[n3c*(nt-2)] ;

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a theta par psi = 2 theta   et F(psi) = f(theta(psi)).  
 */
 
// Calcul des coefficients de Fourier de la fonction 
//   G(psi) = F+(psi) + F_(psi) sin(psi)
// en fonction des coefficients en cos(2l theta) de f:

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
// NB: le facteur 0.25 est du a la normalisation de fftw; si fftw
//     donnait exactement les coef. des sinus, ce facteur serait -0.5.
    	    for ( i = 3; i < nt; i += 2 ) {
		g.set(nm1-i/2) = 0.25 * ( ff0[ n3f*i ] - ff0[ n3f*(i-2) ] ) ;
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
// ... indice du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
	
		double fp = 0.5 * ( g(i) + g(isym) ) ;
		double fm = 0.5 * ( g(i) - g(isym) ) / sinp[i] ;
		ff0[ n3f*i ] =  ( fp + fm ) / x[isym] ;
		ff0[ n3f*isym ] = ( fp - fm ) / x[i] ;
    	    }
	
//... cas particuliers:
	    ff0[0] = g(0) + fmoins0 ;
	    ff0[ n3f*nm1 ] = 0 ;
	    ff0[ n3f*nm1s2 ] = g(nm1s2) / x[nm1s2] ;

	} 	// fin de la boucle sur r 

//-----------------------------------------------------------------------
//  partie sin(m phi) avec m pair : transformation cos((2 l+1) theta) inverse
//-----------------------------------------------------------------------

	j++ ;

	if ( (j != 1) && (j != n1f-1 ) ) {  
//  on effectue le calcul seulement dans les cas ou les coef en phi ne sont 
//  pas nuls 

	for (k=0; k<n3c; k++) {

	    int i0 = n2n3c * j + k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat
	     
// Calcul des coefficients du developpement en T_{2i}(x) de la fonction
//  h(x) := x f(x)  (x=cos(theta)) a partir des coefficients de f 
// (resultat stoke dans le tableau t1 :
	t1[0] = .5 * cf0[0] ;
	for (i=1; i<nm1; i++) t1[i] = .5 * ( cf0[n3c*i] + cf0[n3c*(i-1)] ) ;
	t1[nm1] = .5 * cf0[n3c*(nt-2)] ;

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a theta par psi = 2 theta   et F(psi) = f(theta(psi)).  
 */
 
// Calcul des coefficients de Fourier de la fonction 
//   G(psi) = F+(psi) + F_(psi) sin(psi)
// en fonction des coefficients en cos(2l theta) de f:

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
// NB: le facteur 0.25 est du a la normalisation de fftw; si fftw
//     donnait exactement les coef. des sinus, ce facteur serait -0.5.
    	    for ( i = 3; i < nt; i += 2 ) {
		g.set(nm1-i/2) = 0.25 * ( ff0[ n3f*i ] - ff0[ n3f*(i-2) ] ) ;
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
// ... indice du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
	
		double fp = 0.5 * ( g(i) + g(isym) ) ;
		double fm = 0.5 * ( g(i) - g(isym) ) / sinp[i] ;
		ff0[ n3f*i ] =  ( fp + fm ) / x[isym] ;
		ff0[ n3f*isym ] = ( fp - fm ) / x[i] ;
    	    }
	
//... cas particuliers:
	    ff0[0] = g(0) + fmoins0 ;
	    ff0[ n3f*nm1 ] = 0 ;
	    ff0[ n3f*nm1s2 ] = g(nm1s2) / x[nm1s2] ;


	} 	// fin de la boucle sur r 

	}    // fin du cas ou le calcul etait necessaire (i.e. ou les
	     // coef en phi n'etaient pas nuls)

// On passe au cas m pair suivant:
	j+=3 ;

   }	// fin de la boucle sur les cas m pair

//##
    if (n1f<=3) {
      delete [] t1 ;
      return ;
    }
    
//=======================================================================
//				Cas m impair
//=======================================================================

    j = 2 ;
    
    while (j<n1f-1) {   //le dernier coef en phi n'est pas considere
			// (car nul)

//--------------------------------------------------------------------------
//  partie cos(m phi) avec m impair : transformation sin((2 l) theta) inv.
//--------------------------------------------------------------------------

	for (k=0; k<n3c; k++) {

	    int i0 = n2n3c * j + k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat
	     

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a theta par psi = 2 theta  et F(psi) = f(theta(psi)).  
 */
 

// Calcul des coefficients de Fourier de la fonction 
//   G(psi) = F+(psi) + F_(psi) sin(psi)
// en fonction des coefficients en sin(2l theta) de f:

// Coefficients en sinus de G
//---------------------------
//  Ces coefficients sont egaux aux coefficients pairs du developmt. en
//  sin(2l theta) de f (le facteur -.5 vient de la normalisation
//  de fftw: si fftw donnait reellement les coefficients en sinus,
//  il faudrait le remplacer par un +1) :
    	    for (i=2; i<nm1; i += 2 ) g.set(nm1-i/2) = - .5 * cf0[n3c*i] ;	

// Coefficients en cosinus de G
//-----------------------------
// Ces coefficients se deduisent des coefficients impairs du developmt. 
// en sin(2l theta) de f (le facteur +.25 vient de la normalisation
//  de fftw: si fftw donnait reellement les coefficients en cosinus,
//  il faudrait le remplacer par un +.5)

	    g.set(0) = .5 * cf0[n3c] ;
    	    for ( i = 1; i < nm1s2; i++ ) {
		g.set(i) = .25 * ( cf0[ n3c*(2*i+1) ] - cf0[ n3c*(2*i-1) ] ) ;
    	    }
	    g.set(nm1s2) = - .5 * cf0[ n3c*(nt-2) ] ;


// Transformation de Fourier inverse de G 
//---------------------------------------

// FFT inverse
	    fftw_execute(p) ;

// Valeurs de f deduites de celles de G
//-------------------------------------

    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
	
		double fp = 0.5 * ( g(i) + g(isym) ) / sinp[i] ;
		double fm = 0.5 * ( g(i) - g(isym) ) ;
		ff0[ n3f*i ] = fp + fm ;
		ff0[ n3f*isym ] = fp - fm ;
    	    }
	
//... cas particuliers:
	    ff0[0] = 0. ;
	    ff0[ n3f*nm1 ] = -2. * g(0) ;
	    ff0[ n3f*nm1s2 ] = g(nm1s2) ;

	} 	// fin de la boucle sur r 


//--------------------------------------------------------------------------
//  partie sin(m phi) avec m impair : transformation sin( (2 l) theta) inv.
//--------------------------------------------------------------------------

	j++ ;

	if ( j != n1f-1  ) {  
//  on effectue le calcul seulement dans les cas ou les coef en phi ne sont 
//  pas nuls 

	for (k=0; k<n3c; k++) {

	    int i0 = n2n3c * j + k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat
	     

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a theta par psi = 2 theta  et F(psi) = f(theta(psi)).  
 */
 

// Calcul des coefficients de Fourier de la fonction 
//   G(psi) = F+(psi) + F_(psi) sin(psi)
// en fonction des coefficients en sin(2l theta) de f:

// Coefficients en sinus de G
//---------------------------
//  Ces coefficients sont egaux aux coefficients pairs du developmt. en
//  sin(2l theta) de f (le facteur -.5 vient de la normalisation
//  de fftw: si fftw donnait reellement les coefficients en sinus,
//  il faudrait le remplacer par un +1) :
    	    for (i=2; i<nm1; i += 2 ) g.set(nm1-i/2) = - .5 * cf0[n3c*i] ;	

// Coefficients en cosinus de G
//-----------------------------
// Ces coefficients se deduisent des coefficients impairs du developmt. 
// en sin(2l theta) de f (le facteur +.25 vient de la normalisation
//  de fftw: si fftw donnait reellement les coefficients en cosinus,
//  il faudrait le remplacer par un +.5)

	    g.set(0) = .5 * cf0[n3c] ;
    	    for ( i = 1; i < nm1s2; i++ ) {
		g.set(i) = .25 * ( cf0[ n3c*(2*i+1) ] - cf0[ n3c*(2*i-1) ] ) ;
    	    }
	    g.set(nm1s2) = - .5 * cf0[ n3c*(nt-2) ] ;


// Transformation de Fourier inverse de G 
//---------------------------------------

// FFT inverse
	    fftw_execute(p) ;

// Valeurs de f deduites de celles de G
//-------------------------------------

    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
	
		double fp = 0.5 * ( g(i) + g(isym) ) / sinp[i] ;
		double fm = 0.5 * ( g(i) - g(isym) ) ;
		ff0[ n3f*i ] = fp + fm ;
		ff0[ n3f*isym ] = fp - fm ;
    	    }
	
//... cas particuliers:
	    ff0[0] = 0. ;
	    ff0[ n3f*nm1 ] = -2. * g(0) ;
	    ff0[ n3f*nm1s2 ] = g(nm1s2) ;

	} 	// fin de la boucle sur r 


	}    // fin du cas ou le calcul etait necessaire (i.e. ou les
	     // coef en phi n'etaient pas nuls)

// On passe au cas m impair suivant:
	j+=3 ;

   }	// fin de la boucle sur les cas m impair

    delete [] t1 ;

}
