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
 * Transformation inverse cos(l*theta) ou sin(l*theta) (suivant la 
 *  parite de l'indice m en phi) sur le deuxieme indice (theta)
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
 *   int* dimc	: tableau du nombre d'elements de cc dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dimc[1] >= deg[1] = nt. 
 *
 *   double* cf	: 	tableau des coefficients c_l de la fonction definis
 *			  comme suit (a r et phi fixes)
 *
 *			  pour m pair:
 *			f(theta) = som_{l=0}^{nt-1} c_l sin( l theta ) . 
 *			  pour m impair:
 *			f(theta) = som_{l=0}^{nt-2} c_l cos( l theta ) . 
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
 *			  theta_l =  pi l/(nt-1)       0 <= l <= nt-1 
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
 * Revision 1.5  2016/12/05 16:18:04  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/15 12:48:22  j_novak
 * Corrected namespace declaration.
 *
 * Revision 1.3  2014/10/13 08:53:17  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:18:46  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/12/21 17:06:01  j_novak
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

void citcossins(const int* deg, const int* dimc, double* cf, const int* dimf,
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
	cout << "citcossins: nt > n2f : nt = " << nt << " ,  n2f = " 
	<< n2f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (nt > n2c) {
	cout << "citcossins: nt > n2c : nt = " << nt << " ,  n2c = " 
	<< n2c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n1c > n1f) {
	cout << "citcossins: n1c > n1f : n1c = " << n1c << " ,  n1f = " 
	<< n1f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n3c > n3f) {
	cout << "citcossins: n3c > n3f : n3c = " << n3c << " ,  n3f = " 
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
    double* g = (double*)( malloc( (nm1+2)*sizeof(double) ) ) ;	
    double* t1 = (double*)( malloc( (nm1+2)*sizeof(double) ) ) ;

// Parametres pour la routine FFT991
    int jump = 1 ;
    int inc = 1 ;
    int lot = 1 ;
    int isign = 1 ;
    
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

//--------------------------------------------------------------------------
//  partie cos(m phi) avec m pair : transformation sin(l theta) inv.
//--------------------------------------------------------------------------

      	for (k=0; k<n3c; k++) {

	    int i0 = n2n3c * j + k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat
	     
// Coefficients impairs de G
//--------------------------
 
	    g[1] = 0 ;
	    for (i=2; i<nm1; i += 2 ) g[i+1] = -0.5 * cf0[ n3c*i ] ;	
	    g[nt] = 0 ;	


// Coefficients pairs de G
//------------------------

	    g[0] = .5 * cf0[n3c] ;
    	    for ( i = 3; i < nt; i += 2 ) {
		g[i-1] = .25 * ( cf0[ n3c*i ] - cf0[ n3c*(i-2) ] ) ;
    	    }
	    g[nm1] = - .5 * cf0[ n3c*(nt-2) ] ;
	    
// Transformation de Fourier inverse de G 
//---------------------------------------

// FFT inverse
    	    F77_fft991( g, t1, trigo, facto, &inc, &jump, &nm1, &lot, &isign) ;

// Valeurs de f deduites de celles de G
//-------------------------------------

    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice du pt symetrique de theta par rapport a pi/2:
		int isym = nm1 - i ; 
	
		double fp = 0.5 * ( g[i] + g[isym] ) / sinp[i]  ;
		double fm = 0.5 * ( g[i] - g[isym] ) ;
		ff0[ n3f*i ] = fp + fm ;
		ff0[ n3f*isym ] = fp - fm ;
    	    }
	
//... cas particuliers:
	    ff0[0] = 0. ;
	    ff0[ n3f*nm1 ] = -2*g[0] ;
	    ff0[ n3f*nm1s2 ] = g[nm1s2] ;


	} 	// fin de la boucle sur r 

//--------------------------------------------------------------------------
//  partie sin(m phi) avec m pair : transformation sin(l theta) inv.
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
	     
// Coefficients impairs de G
//--------------------------
 
	    g[1] = 0 ;
	    for (i=2; i<nm1; i += 2 ) g[i+1] = -0.5 * cf0[ n3c*i ] ;	
	    g[nt] = 0 ;	


// Coefficients pairs de G
//------------------------

	    g[0] = .5 * cf0[n3c] ;
    	    for ( i = 3; i < nt; i += 2 ) {
		g[i-1] = .25 * ( cf0[ n3c*i ] - cf0[ n3c*(i-2) ] ) ;
    	    }
	    g[nm1] = - .5 * cf0[ n3c*(nt-2) ] ;
	    
// Transformation de Fourier inverse de G 
//---------------------------------------

// FFT inverse
    	    F77_fft991( g, t1, trigo, facto, &inc, &jump, &nm1, &lot, &isign) ;

// Valeurs de f deduites de celles de G
//-------------------------------------

    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
	
		double fp = 0.5 * ( g[i] + g[isym] ) / sinp[i]  ;
		double fm = 0.5 * ( g[i] - g[isym] ) ;
		ff0[ n3f*i ] = fp + fm ;
		ff0[ n3f*isym ] = fp - fm ;
    	    }
	
//... cas particuliers:
	    ff0[0] = 0. ;
	    ff0[ n3f*nm1 ] = -2*g[0] ;
	    ff0[ n3f*nm1s2 ] = g[nm1s2] ;


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

//-----------------------------------------------------------------------
//  partie cos(m phi) avec m impair : transformation cos( l theta) inverse
//-----------------------------------------------------------------------

	for (k=0; k<n3c; k++) {

	    int i0 = n2n3c * j + k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3f * j + k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau resultat
	     
 
// Coefficients impairs de G
//--------------------------
 
	    double c1 = cf0[n3c] ;

    	    double som = 0;
	    ff0[n3f] = 0 ;
    	    for ( i = 3; i < nt; i += 2 ) {
	    	ff0[ n3f*i ] = cf0[ n3c*i ] - c1 ;
		som += ff0[ n3f*i ] ;
    	    }	

// Valeur en theta=0 de la partie antisymetrique de F, F_ :
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
//   f.
// NB: le facteur 0.5 est du a la normalisation de fft991; si fft991
//     donnait exactement les coef. des cosinus, ce facteur serait 1.

	    g[0] = cf0[0] ;
    	    for (i=2; i<nm1; i += 2 ) g[i] = 0.5 * cf0[ n3c*i ] ;	
    	    g[nm1] = cf0[ n3c*nm1 ] ;

// Transformation de Fourier inverse de G 
//---------------------------------------

// FFT inverse
    	    F77_fft991( g, t1, trigo, facto, &inc, &jump, &nm1, &lot, &isign) ;

// Valeurs de f deduites de celles de G
//-------------------------------------

    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice du pt symetrique de theta par rapport a pi/2:
		int isym = nm1 - i ; 
	
		double fp = 0.5 * ( g[i] + g[isym] ) ;
		double fm = 0.5 * ( g[i] - g[isym] ) / sinp[i] ;
		ff0[ n3f*i ] = fp + fm ;
		ff0[ n3f*isym ] = fp - fm ;
    	    }
	
//... cas particuliers:
	    ff0[0] = g[0] + fmoins0 ;
	    ff0[ n3f*nm1 ] = g[0] - fmoins0 ;
	    ff0[ n3f*nm1s2 ] = g[nm1s2] ;


	} 	// fin de la boucle sur r 

//-----------------------------------------------------------------------
//  partie sin(m phi) avec m impair : transformation cos(l theta) inverse
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
	     
// Coefficients impairs de G
//--------------------------
 
	    double c1 = cf0[n3c] ;

    	    double som = 0;
	    ff0[n3f] = 0 ;
    	    for ( i = 3; i < nt; i += 2 ) {
	    	ff0[ n3f*i ] = cf0[ n3c*i ] - c1 ;
		som += ff0[ n3f*i ] ;
    	    }	

// Valeur en theta=0 de la partie antisymetrique de F, F_ :
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
//   f.
// NB: le facteur 0.5 est du a la normalisation de fft991; si fft991
//     donnait exactement les coef. des cosinus, ce facteur serait 1.

	    g[0] = cf0[0] ;
    	    for (i=2; i<nm1; i += 2 ) g[i] = 0.5 * cf0[ n3c*i ] ;	
    	    g[nm1] = cf0[ n3c*nm1 ] ;

// Transformation de Fourier inverse de G 
//---------------------------------------

// FFT inverse
    	    F77_fft991( g, t1, trigo, facto, &inc, &jump, &nm1, &lot, &isign) ;

// Valeurs de f deduites de celles de G
//-------------------------------------

    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice du pt symetrique de theta par rapport a pi/2:
		int isym = nm1 - i ; 
	
		double fp = 0.5 * ( g[i] + g[isym] ) ;
		double fm = 0.5 * ( g[i] - g[isym] ) / sinp[i] ;
		ff0[ n3f*i ] = fp + fm ;
		ff0[ n3f*isym ] = fp - fm ;
    	    }
	
//... cas particuliers:
	    ff0[0] = g[0] + fmoins0 ;
	    ff0[ n3f*nm1 ] = g[0] - fmoins0 ;
	    ff0[ n3f*nm1s2 ] = g[nm1s2] ;


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
