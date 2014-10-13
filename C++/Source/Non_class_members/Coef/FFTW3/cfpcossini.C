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


char cfpcossini_C[] = "$Header$" ;

/*
 * Transformation de Fourier sur le premier indice d'un tableau 3-D
 * Cas d'une fonction antisymetrique par la transformation phi --> phi + Pi
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune 
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en phi est  np = deg[0] et doit etre de la forme
 * 			np = 2*p
 *   int* dim	: tableau du nombre d'elements de ff dans chacune des trois 
 *	          dimensions.
 *		  On doit avoir  dim[0] >= deg[0] + 2  = np + 2. 
 *
 * Entree/Sortie : 
 * ---------------
 *   double* cf : entree: tableau des valeurs de la fonction f aux np points de
 *                        de collocation 
 *				phi_k = Pi k/np      0 <= k <= np-1
 *			  La convention de stokage utilisee est la suivante
 *	cf[ dim[2]*dim[1]*k + dim[2]*j + i ] = f(phi_k)    0 <= k <= np-1,
 *			 les indices j et i correspondant respectivement
 *			  a theta et a r. 
 *			  L'espace memoire correspondant au
 *                        pointeur cf doit etre dim[0]*dim[1]*dim[2] et doit 
 *			  etre alloue avant l'appel a la routine.	
 * 
 * 		  sortie: tableau des coefficients de la fonction suivant	 
 *			  la convention de stokage
 *		cf[ dim[2]*dim[1]*k + dim[2]*j + i ] = c_k      0<= k < np,
 *			  ou les indices j et i correspondent respectivement
 *			  a theta et a r et ou les c_k sont les coefficients 
 *			  du developpement de f en series de Fourier:
 *
 *			f(phi) = c_0 cos(phi) + c_2 sin(phi)
 *			    + som_{l=1}^{np/2-1} c_{2l+1} cos( (2l+1) phi )
 *				 + c_{2l+2} sin( (2l+1) phi ),
 *
 *			En particulier: c_1 = 0 et c_{np+1} = 0
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
 * Revision 1.2  2002/10/16 14:36:43  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  2000/09/08  15:55:04  eric
 * Premiere version testee.
 *
 * Revision 2.0  2000/09/07  15:15:06  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>

// Headers Lorene
#include "headcpp.h"
#include "proto.h"

namespace Lorene {
//*****************************************************************************

void cfpcossini(const int* deg, const int* dim, double* cf) {

    // Dimensions du tableau cf :
    int n1 = dim[0] ;
    int nt = dim[1] ;
    int nr = dim[2] ;
    int ntnr = nt * nr ; 

    // Nombres de degres de liberte en phi :    
    int np = deg[0] ;
    
    // Tests de dimension:
    if (np+2 > n1) {
	cout << "cfpcossini: np+2 > n1 : np+2 = " << np+2 << " ,  n1 = " 
	<< n1 << endl ;
	abort() ;
    }

    // Tableau contenant les points de 0 a 2 Pi (et non seulement de 0 a Pi)
    int np2 = 2*np ; 
    int deg2[] = {np2, nt, nr} ; 
    int dim2[] = {np2+2, nt, nr} ;
    
    double* cf2 = new double[(np2+2)*nt*nr] ; 
    
    // Recopie des valeurs pour phi dans [0, Pi[ :
    for (int k=0; k<np; k++) {
	for (int ij = 0; ij <ntnr; ij++) {
	    cf2[k*ntnr + ij] = cf[k*ntnr + ij] ; 
	}
    }
    
    // Valeurs pour phi dans [Pi, 2Pi[ obtenues par antisymetrie:
    int npntnr = np * ntnr ; 
    for (int k=0; k<np; k++) {
	for (int ij = 0; ij <ntnr; ij++) {
	    cf2[npntnr + k*ntnr + ij] = - cf[k*ntnr + ij] ; 
	}
    }
    
    // Transformation de Fourier sur cf2 : 
    cfpcossin(deg2, dim2, cf2) ;
    
    // Recopie des coefficients dans cf :
    
    // Terme k=0   cos(phi)
    for (int ij = 0; ij <ntnr; ij++) {
	    cf[ij] = cf2[2*ntnr + ij] ; 
    }
    
    // Terme k=1
    for (int ij = 0; ij <ntnr; ij++) {
	    cf[ntnr + ij] = 0 ; 
    }
    
    // Terme k=2   sin(phi)
    for (int ij = 0; ij <ntnr; ij++) {
	    cf[2*ntnr + ij] = cf2[3*ntnr + ij] ; 
    }
    
    // Termes suivants:
    for (int k=3; k<np; k+=2) {
	int k2 = 2*(k-1) + 2 ;
	// Terme en cosinus 
	for (int ij = 0; ij <ntnr; ij++) {
	    cf[k*ntnr + ij] = cf2[k2*ntnr + ij] ; 
	}
	// Terme en sinus
	k2++ ; 
	int k1 = k+1 ; 
	for (int ij = 0; ij <ntnr; ij++) {
	    cf[k1*ntnr + ij] = cf2[k2*ntnr + ij] ; 
	}	
    }
    
    // Terme k=np+1 
    for (int ij = 0; ij <ntnr; ij++) {
	    cf[(np+1)*ntnr + ij] = 0 ; 
    }
    
    
    
    // Menage
    delete [] cf2 ;    
 
}
}
