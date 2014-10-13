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


char cfrchebp_C[] = "$Header$" ;


/*
 * Transformation de Tchebyshev (cas rare) sur le troisieme indice (indice 
 *  correspondant a r) d'un tableau 3-D decrivant une fonction paire. 
 *  Utilise la bibliotheque fftw.
 *
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune 
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en r est  nr = deg[2] et doit etre de la forme
 * 			nr = 2*p + 1 
 *   int* dimf	: tableau du nombre d'elements de ff dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dimf[2] >= deg[2] = nr. 
 *		  NB: pour dimf[0] = 1 (un seul point en phi), la transformation
 *		      est bien effectuee.
 *		      pour dimf[0] > 1 (plus d'un point en phi), la
 *		      transformation n'est effectuee que pour les indices (en phi)
 *		      j != 1 et j != dimf[0]-1 (cf. commentaires sur borne_phi).
 *
 *   double* ff : tableau des valeurs de la fonction aux nr points de
 *                        de collocation
 *
 *			  x_i = sin( pi/2 i/(nr-1) )      0 <= i <= nr-1 
 *
 *		    Les valeurs de la fonction doivent etre stokees dans le 
 *		    tableau ff comme suit
 *			   f( x_i ) = ff[ dimf[1]*dimf[2] * j + dimf[2] * k + i ]
 *		    ou j et k sont les indices correspondant a phi et theta 
 *		    respectivement.
 * 		    L'espace memoire correspondant a ce pointeur doit etre 
 *		    dimf[0]*dimf[1]*dimf[2] et doit etre alloue avant l'appel a 
 *		    la routine.	 
 *
 *   int* dimc	: tableau du nombre d'elements de cf dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dimc[2] >= deg[2] = nr. 
 *
 * Sortie:
 * -------
 *   double* cf	:   tableau des nr coefficients c_i de la fonction definis
 *		    comme suit (a theta et phi fixes)
 *
 *				    f(x) = som_{i=0}^{nr-1} c_i T_{2i}(x) , 
 *
 *		    ou T_{2i}(x) designe le polynome de Tchebyshev de degre 2i. 	 
 *		    Les coefficients c_i (0 <= i <= nr-1) sont stokes dans
 *		    le tableau cf comme suit
 *			   c_i = cf[ dimc[1]*dimc[2] * j + dimc[2] * k + i ]
 *		    ou j et k sont les indices correspondant a phi et theta 
 *		    respectivement.
 * 		    L'espace memoire correspondant a ce pointeur doit etre 
 *		    dimc[0]*dimc[1]*dimc[2] et doit avoir ete alloue avant 
 *		    l'appel a la routine.	 
 *
 * NB: Si le pointeur ff est egal a cf, la routine ne travaille que sur un 
 *     seul tableau, qui constitue une entree/sortie.
 */
 
/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:18  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:18:47  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/12/21 17:06:02  j_novak
 * Added all files for using fftw3.
 *
 * Revision 1.4  2003/01/31 10:31:23  e_gourgoulhon
 * Suppressed the directive #include <malloc.h> for malloc is defined
 * in <stdlib.h>
 *
 * Revision 1.3  2002/10/16 14:36:44  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/09/09 13:00:39  e_gourgoulhon
 * Modification of declaration of Fortran 77 prototypes for
 * a better portability (in particular on IBM AIX systems):
 * All Fortran subroutine names are now written F77_* and are
 * defined in the new file C++/Include/proto_f77.h.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/02/22  15:48:30  hyc
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
namespace Lorene {
fftw_plan prepare_fft(int, Tbl*&) ;
double* cheb_ini(const int) ;

//*****************************************************************************

void cfrchebp(const int* deg, const int* dimf, double* ff, const int* dimc,
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

// Nombres de degres de liberte en r :    
    int nr = deg[2] ;
    
// Tests de dimension:
    if (nr > n3f) {
	cout << "cfrchebp: nr > n3f : nr = " << nr << " ,  n3f = " 
	<< n3f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (nr > n3c) {
	cout << "cfrchebp: nr > n3c : nr = " << nr << " ,  n3c = " 
	<< n3c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n1f > n1c) {
	cout << "cfrchebp: n1f > n1c : n1f = " << n1f << " ,  n1c = " 
	<< n1c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n2f > n2c) {
	cout << "cfrchebp: n2f > n2c : n2f = " << n2f << " ,  n2c = " 
	<< n2c << endl ;
	abort () ;
	exit(-1) ;
    }

// Nombre de points pour la FFT:
    int nm1 = nr - 1;
    int nm1s2 = nm1 / 2;

// Recherche des tables pour la FFT:
    Tbl* pg = 0x0 ;
    fftw_plan p = prepare_fft(nm1, pg) ;
    Tbl& g = *pg ;

// Recherche de la table des sin(psi) :
    double* sinp = cheb_ini(nr);	
	
// boucle sur phi et theta

    int n2n3f = n2f * n3f ;
    int n2n3c = n2c * n3c ;

/*   
 * Borne de la boucle sur phi: 
 *    si n1f = 1, on effectue la boucle une fois seulement.
 *    si n1f > 1, on va jusqu'a j = n1f-2 en sautant j = 1 (les coefficients
 *	j=n1f-1 et j=0 ne sont pas consideres car nuls). 
 */
    int borne_phi = ( n1f > 1 ) ? n1f-1 : 1 ;

    for (j=0; j< borne_phi; j++) {
    
	if (j==1) continue ;	// on ne traite pas le terme en sin(0 phi)

	for (k=0; k<n2f; k++) {

	    int i0 = n2n3f * j + n3f * k ; // indice de depart 
	    double* ff0 = ff + i0 ;    // tableau des donnees a transformer

	    i0 = n2n3c * j + n3c * k ; // indice de depart 
	    double* cf0 = cf + i0 ;    // tableau resultat

/*
 * NB: dans les commentaires qui suivent, psi designe la variable de [0, pi]
 *     reliee a x par  x = cos(psi/2)   et F(psi) = f(x(psi)).  
 */
 
// Valeur en psi=0 de la partie antisymetrique de F, F_ :
    	    double fmoins0 = 0.5 * ( ff0[nm1] - ff0[0] );

// Fonction G(psi) = F+(psi) + F_(psi) sin(psi) 
//---------------------------------------------
    	    for ( i = 1; i < nm1s2 ; i++ ) {
// ... indice (dans le tableau g) du pt symetrique de psi par rapport a pi/2:
		int isym = nm1 - i ; 
// ... indice (dans le tableau ff0) du point x correspondant a psi
		int ix = nm1 - i ;
// ... indice (dans le tableau ff0) du point x correspondant a sym(psi)
		int ixsym = nm1 -  isym ;
	
// ... F+(psi)
		double fp = 0.5 * ( ff0[ix] + ff0[ixsym] ) ;	
// ... F_(psi) sin(psi)
		double fms = 0.5 * ( ff0[ix] - ff0[ixsym] ) * sinp[i] ; 
		g.set(i) = fp + fms ;
		g.set(isym) = fp - fms ;
    	    }
//... cas particuliers:
    	    g.set(0) = 0.5 * ( ff0[nm1] + ff0[0] );
    	    g.set(nm1s2) = ff0[nm1s2];

// Developpement de G en series de Fourier par une FFT
//----------------------------------------------------

	    fftw_execute(p) ;

// Coefficients pairs du developmt. de Tchebyshev de f
//----------------------------------------------------
//  Ces coefficients sont egaux aux coefficients en cosinus du developpement
//  de G en series de Fourier (le facteur 2/nm1 vient de la normalisation
//  de fftw) :

	    double fac = 2./double(nm1) ;
	    cf0[0] = g(0) / double(nm1) ;
	    for (i=2; i<nm1; i += 2) cf0[i] = fac*g(i/2) ;
	    cf0[nm1] = g(nm1s2) / double(nm1) ;

// Coefficients impairs du developmt. de Tchebyshev de f
//------------------------------------------------------
// 1. Coef. c'_k (recurrence amorcee a partir de zero)
//  Le 4/nm1 en facteur de g(i) est du a la normalisation de fftw
//  (si fftw donnait reellement les coef. en sinus, il faudrait le
//   remplacer par un -2.)
	    fac *= 2. ;
    	    cf0[1] = 0. ;
    	    double som = 0.;
    	    for ( i = 3; i < nr; i += 2 ) {
		cf0[i] = cf0[i-2] + fac * g(nm1-i/2) ;
	    	som += cf0[i] ;
    	    }

// 2. Calcul de c_1 :
    	    double c1 = ( fmoins0 - som ) / nm1s2 ;

// 3. Coef. c_k avec k impair:	
    	    cf0[1] = c1 ;
    	    for ( i = 3; i < nr; i += 2 ) cf0[i] += c1 ;

	} 	// fin de la boucle sur theta 
   }	// fin de la boucle sur phi


}
}
