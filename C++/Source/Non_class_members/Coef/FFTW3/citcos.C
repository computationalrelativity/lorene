/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 2002 Jerome Novak
 *   
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
 * Transformation en cos(l*theta) inverse sur le deuxieme indice (theta)
 *  d'un tableau 3-D representant une fonction quelconque (theta variant de 0
 *  a pi). Utilise la bibliotheque fftw.
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune 
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en theta est  nt = deg[1] et doit etre de la forme
 * 			nt = 2*p + 1 
 *   int* dimc	: tableau du nombre d'elements de cc dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dimc[1] >= deg[1] = nt. 
 *		  NB: pour dimc[0] = 1 (un seul point en phi), la transformation
 *		      est bien effectuee.
 *		      pour dimc[0] > 1 (plus d'un point en phi), la
 *		      transformation n'est effectuee que pour les indices (en phi)
 *		      j != 1 et j != dimc[0]-1 (cf. commentaires sur borne_phi).
 *
 *   double* cf	: 	tableau des coefficients c_l de la fonction definis
 *			  comme suit (a r et phi fixes)
 *
 *			   f(theta) = som_{l=0}^{nt-1} c_l cos( l theta ) . 
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
 * Revision 1.4  2016/12/05 16:18:05  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:20  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:18:49  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/12/21 17:06:03  j_novak
 * Added all files for using fftw3.
 *
 * Revision 1.2  2003/01/31 10:31:23  e_gourgoulhon
 * Suppressed the directive #include <malloc.h> for malloc is defined
 * in <stdlib.h>
 *
 * Revision 1.1  2002/11/12 17:43:53  j_novak
 * Added transformatin functions for T_COS basis.
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
//*****************************************************************************

void citcos(const int* deg, const int* dimc, double* cf, const int* dimf,
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
	cout << "citcos: nt > n2f : nt = " << nt << " ,  n2f = " 
	<< n2f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (nt > n2c) {
	cout << "citcos: nt > n2c : nt = " << nt << " ,  n2c = " 
	<< n2c << endl ;
	abort () ;
	exit(-1) ;
    }
    if ( (n1f > 1) && (n1c > n1f) ) {
	cout << "citcos: n1c > n1f : n1c = " << n1c << " ,  n1f = " 
	<< n1f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n3c > n3f) {
	cout << "citcos: n3c > n3f : n3c = " << n3c << " ,  n3f = " 
	<< n3f << endl ;
	abort () ;
	exit(-1) ;
    }

// Nombre de points pour la FFT:
    int nm1 = nt - 1;
    int nm1s2 = nm1 / 2;

// Recherche des tables pour la FFT:
    Tbl* pg = 0x0 ;
    fftw_plan p = back_fft(nm1, pg) ;
    Tbl& g = *pg ;

// Recherche de la table des sin(psi) :
    double* sinp = cheb_ini(nt);	
	
// boucle sur phi et r (les boucles vont resp. de 0 a max(dimc[0]-2,0) et 
//			0 a dimc[2]-1 )

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
	     
/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a theta par psi = 2 theta   et F(psi) = f(theta(psi)).  
 */
 
// Calcul des coefficients de Fourier de la fonction 
//   G(psi) = F+(psi) + F_(psi) sin(psi)
// en fonction des coefficients en cos(2l theta) de f:

// Coefficients impairs de G
//--------------------------
 
	    double c1 = cf0[n3c] ;

    	    double som = 0;
	    ff0[n3f] = 0 ;
    	    for ( i = 3; i < nt; i += 2 ) {
	    	ff0[ n3f*i ] = cf0[ n3c*i ] - c1 ;
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

	    g.set(0) = cf0[0] ;
    	    for (i=1; i<nm1s2; i++ ) g.set(i) = 0.5 * cf0[ n3c*2*i ] ;	
    	    g.set(nm1s2) = cf0[ n3c*nm1 ] ;

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
		ff0[ n3f*i ] = fp + fm ;
		ff0[ n3f*isym ] = fp - fm ;
    	    }
	
//... cas particuliers:
	    ff0[0] = g(0) + fmoins0 ;
	    ff0[ n3f*nm1 ] = g(0) - fmoins0 ;
	    ff0[ n3f*nm1s2 ] = g(nm1s2) ;


	} 	// fin de la boucle sur r 
   }	// fin de la boucle sur phi

    
}
}
