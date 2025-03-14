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
 * Transformation de Fourier inverse sur le premier indice d'un tableau 3-D
 *  par le biais de la routine FFT Fortran FFT991
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune 
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en phi est  np = deg[0] et doit etre de la forme
 * 			np = 2*p
 *   int* dimc	: tableau du nombre d'elements dans chacune des trois 
 *    		  dimensions du tableau cf
 *		  On doit avoir  dimc[0] >= deg[0] + 2  = np + 2. 
 *
 *   int* dimf	: tableau du nombre d'elements dans chacune des trois 
 *    		  dimensions du tableau ff
 *		  On doit avoir  dimf[0] >= deg[0]  = np . 
 *
 *
 * Entree / sortie :
 * -----------------
 *   double* cf : entree: tableau des coefficients de la fonction f; 
 *			  L'espace memoire correspondant a ce
 *                        pointeur doit etre dimc[0]*dimc[1]*dimc[2] et doit 
 *			  avoir ete alloue avant l'appel a la routine.
 *			  La convention de stokage est la suivante:
 *		cf[ dimc[2]*dimc[1]*k + dimc[2]*j + i ] = c_k      0<= k <= np,
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
 *   !!!! CE TABLEAU EST DETRUIT EN SORTIE !!!!!
 *
 * Sortie:
 * -------
 *  double* ff : tableau des valeurs de la fonction aux points de
 *                        de collocation
 *
 *				phi_k = T k/np      0 <= k <= np-1
 *
 *			  suivant la convention de stokage:
 *	ff[ dimf[2]*dimf[1]*k + dimf[2]*j + i ] = f(phi_k)    0 <= k <= np-1,
 *			 les indices j et i correspondant respectivement
 *			  a theta et a r. 
 *			  L'espace memoire correspondant a ce
 *                        pointeur doit etre dimf[0]*dimf[1]*dimf[2] et doit 
 *			  avoir ete alloue avant l'appel a la routine.
 */

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:18:05  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:19  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:18:49  j_novak
 * Modified #include directives to use c++ syntax.
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
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/02/22  15:43:58  hyc
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
fftw_plan back_fft(int, Tbl*&) ;
//*****************************************************************************

void cipcossin(const int* deg, const int* dimc, const int* dimf, 
	       double* cf, double* ff)
{

// Dimensions des tableaux ff et cf  :
    int n1f = dimf[0] ;
    int n2f = dimf[1] ;
    int n3f = dimf[2] ;
    int n1c = dimc[0] ;
    int n2c = dimc[1] ;
    int n3c = dimc[2] ;

// Nombres de degres de liberte en phi:    
    int np = deg[0] ;

// Tests de dimension:
    if (np+2 > n1c) {
	cout << "cipcossin: np+2 > n1c : np = " << np << " ,  n1c = " 
	<< n1c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (np > n1f) {
	cout << "cipcossin: np > n1f : np = " << np << " ,  n1f = " 
	<< n1f << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n3f > n3c) {
	cout << "cipcossin: n3f > n3c : n3f = " << n3f << " ,  n3c = " 
	<< n3c << endl ;
	abort () ;
	exit(-1) ;
    }
    if (n2f > n2c) {
	cout << "cipcossin: n2f > n2c : n2f = " << n2f << " ,  n2c = " 
	<< n2c << endl ;
	abort () ;
	exit(-1) ;
    }

    // Recherche des tables
    Tbl* pg = 0x0 ;
    fftw_plan p = back_fft(np, pg) ;
    Tbl& g = *pg ;    
    int index ;
    int n2n3c = n2c*n3c ;
    int n2n3f = n2f*n3f ;

// boucle sur r et theta

    for (int j=0; j<n2c; j++) {
      for (int k=0; k<n3c; k++) {
	index = n3c * j + k ;
	double* debut = cf + index ;
	g.set(0) = *debut ;
	debut += 2*n2n3c ;
	for (int i=1; i<np/2; i++) {
	  int isym = np - i ;
	  g.set(i) = 0.5 * (*debut) ; debut += n2n3c ;
	  g.set(isym) = -0.5 * (*debut) ; debut += n2n3c ;
	}
	g.set(np/2) = *debut ;
	    
// FFT inverse

	fftw_execute(p) ;

// On range dans ff
	if ((j<n2f) && (k<n3f)) {
	  debut = ff +  n3f * j + k ;
	  for (int i=0; i<np; i++) {
	    *debut = g(i) ;
	    debut += n2n3f ;
	  }
	}
      } 	// fin de la boucle sur r 
    }	// fin de la boucle sur theta
    
}
}
