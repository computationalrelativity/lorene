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


char cipcossini_C[] = "$Header$" ;

/*
 * Transformation de Fourier inverse sur le premier indice d'un tableau 3-D
 * Cas d'une fonction antisymetrique par la transformation phi --> phi + Pi
 * 
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune 
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en phi est  np = deg[0] et doit etre de la forme
 * 			np = 2^p 3^q 5^r 
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
 *		cf[ dimc[2]*dimc[1]*k + dimc[2]*j + i ] = c_k      0<= k < np,
 *			  ou les indices j et i correspondent respectivement
 *			  a theta et a r et ou les c_k sont les coefficients 
 *			  du developpement de f en series de Fourier:
 *
 *			f(phi) = c_0 cos(phi) + c_2 sin(phi)
 *			    + som_{l=1}^{np/2-1} c_{2l+1} cos( (2l+1) phi )
 *				 + c_{2l+2} sin( (2l+1) phi ),
 *
 *			En particulier: c_1 = 0 et c_{np+1} = 0
 *
 *   !!!! CE TABLEAU EST DETRUIT EN SORTIE !!!!!
 *
 * Sortie:
 * -------
 *  double* ff : tableau des valeurs de la fonction aux points de
 *                        de collocation
 *
 *				phi_k = Pi k/np      0 <= k <= np-1
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
 * Revision 1.2  2014/10/06 15:18:49  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/12/21 17:06:02  j_novak
 * Added all files for using fftw3.
 *
 * Revision 1.2  2002/10/16 14:36:53  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  2000/09/08  15:55:21  eric
 * Premiere version testee.
 *
 * Revision 2.0  2000/09/07  15:15:17  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>

#include "headcpp.h"

// Headers Lorene
#include "proto.h"

//*****************************************************************************

void cipcossini(const int* deg, const int* dimc, const int* dimf, 
	       double* cf, double* ff)
{

    // Nombres de degres de liberte en phi:    
    int np = deg[0] ;
    
    // Tableaux pour echantillonnage sur [0,2 Pi[ :
    int np2 = 2*np ; 
    
    int deg2[] = {np2, deg[1], deg[2]} ; 
    int dimc2[] = {np2+2, dimc[1], dimc[2]} ; 
    int dimf2[] = {np2, dimf[1], dimf[2]} ; 
    
    double* cf2 = new double[ dimc2[0]*dimc2[1]*dimc2[2] ] ; 
    double* ff2 = new double[ dimf2[0]*dimf2[1]*dimf2[2] ] ; 
    
    // Recopie des coefficients dans cf2 :
    int ntnrc = dimc[1] * dimc[2] ; 

    // Harmonique m=0 (cosinus et sinus) mise a zero
    for (int ij = 0; ij <2*ntnrc; ij++) {
	     cf2[ij] = 0 ; 
    }
    
    // Harmonique m=1 (cosinus) 
    for (int ij = 0; ij <ntnrc; ij++) {
	     cf2[2*ntnrc + ij] = cf[ij] ; 
    }
    
    // Harmonique m=1 (sinus) 
    for (int ij = 0; ij <ntnrc; ij++) {
	     cf2[3*ntnrc + ij] = cf[2*ntnrc + ij] ; 
    }

    for (int k2=4; k2<np2; k2 +=4) {
	// Harmoniques paires : 
	for (int ij = 0; ij <2*ntnrc; ij++) {
	     cf2[k2*ntnrc + ij] = 0 ; 
	}

	// Harmoniques impaires (cosinus et sinus)
	int k = k2 / 2 + 1 ; 
	for (int ij = 0; ij <2*ntnrc; ij++) {
	     cf2[(k2+2)*ntnrc + ij] = cf[k*ntnrc + ij] ; 
	}
    }

    // Mise a zero des deux derniers coefficients:
    for (int ij = 0; ij <2*ntnrc; ij++) {
	cf2[np2*ntnrc + ij] = 0 ; 
    }
        
    // Transformation de Fourier inverse sur [0, 2 Pi]
    cipcossin(deg2, dimc2, dimf2, cf2, ff2) ; 
    
    // Recopie des valeurs de la fonction dans ff : 
    int ntnrf = dimf[1] * dimf[2] ; 
    for (int k=0; k<np; k++) {
	for (int ij = 0; ij <ntnrf; ij++) {
	     ff[k*ntnrf + ij] = ff2[k*ntnrf + ij] ; 
	}
    }
    
    
    // Menage
    delete [] cf2 ;    
    delete [] ff2 ;     
    
    
}
