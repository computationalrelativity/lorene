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


char legendre_C[] = "$Header$" ;

/*
 * Calcule les valeurs des fonctions de Legendre associees 
 *   P_l^m(cos(theta)) / (2m-1)!!
 * aux points 
 *	theta_j = pi/2 j/(nt-1)	    0 <= j <= nt-1 
 * qui echantillonnent uniformement l'intervalle [0, pi/2].
 *
 *
 * Entree:
 * -------
 * int m    : ordre de la fonction de Legendre associee P_l^m
 * int nt   : nombre de points en theta
 * 
 * Sortie (valeur de retour) :
 * -------------------------
 * double* legendre :	ensemble des (nt-m)*nt valeurs 
 *			    P_l^m(cos(theta))/(2m-1)!! 
 *			stokees comme suit:
 *
 *	    legendre[nt* (l-m) + j] = P_l^m( cos(theta_j) ) / (2m-1)!! 
 *
 *			avec   m <= l <= nt-1.
 *
 * NB: Cette routine effectue le calcul a chaque appel et ne renvoie pas
 *     un pointeur sur des valeurs precedemment calculees.
 */


/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  1999/02/22  15:37:13  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers du C++
#include <iostream.h>

// headers du C
#include <assert.h>
#include <math.h>
#include <malloc.h>

//******************************************************************************

double* legendre(int m, int nt) {

int i, j, l ;
   
    int lmax = nt - 1 ; 
    assert(m >= 0) ; 
    assert(m <= lmax) ; 

    double dt = M_PI / double(2*(nt-1)) ;

// Allocation memoire pour le tableau resultat 
//--------------------------------------------

    double* resu = (double *)(malloc( (lmax-m+1)*nt * sizeof(double) )) ;

    // Tableau de travail
    double* cost = (double*)( malloc( nt*sizeof(double) ) ) ;

//-----------------------
// 1/ Calcul de P_m^m
//-----------------------

    if (m==0) {
	for (j=0; j<nt; j++) {
	    resu[j] = 1. ;	    // P_0^0(x) = 1. 
	}
    }
    else {

//... P_m^m(x) = (-1)^m  (1-x^2)^{m/2}  <--- cette formule donne un P_m^m
//					     plus petit par un facteur
//					     (2m-1)!! que celui de la litterature

	for (j=0; j<nt; j++) {
	    double y = 1. ; 
	    double s = sin(j*dt) ;
		for (i=1 ; i<2*m; i+=2) {
		    y *= - s ;		
// NB: Pour obtenir le P_m^m de la litterature, il faudrait remplacer la ligne
//  ci-dessus par :  y *= - i*s ;
		}
	    resu[j] = y ;	    
//##	    resu[j] = pow(-s, double(m)) ; 
	}
    }	// fin du cas m != 0
    
    if (lmax==m) {
	free (cost) ;
	return resu ;
    }
    else {

//-----------------------
// 2/ Calcul de P_{m+1}^m
//-----------------------

//... Calcul des cos( theta_j ) :
	for (j=0; j<nt; j++) {
	    cost[j] = cos(j*dt)  ;	    
	}

	for (j=0; j<nt; j++) {
	    resu[nt+j] = cost[j] * (2.*m+1) * resu[j] ;	    
	}
	
//-----------------------
// 3/ Calcul de P_l^m  pour m+2 <= l <= lmax 
//-----------------------

	for (l=m+2; l < lmax+1 ; l++) {
	    int i_l = nt*(l-m) ;
	    int i_lm1 = nt*(l-1-m) ; 
	    int i_lm2 = nt*(l-2-m) ; 
	    int a = 2*l - 1 ;
	    int b = l + m - 1 ; 
	    int c = l - m ;

	    for (j=0; j<nt; j++) {
		resu[i_l+j] = ( cost[j] * a * resu[i_lm1+j] 
				- b * resu[i_lm2+j] ) / c ;
	    }
	}

	free (cost) ;
	return resu ; 
    
    } // fin du cas lmax > m 

}



