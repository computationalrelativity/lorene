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


char legendre_norm_C[] = "$Header$" ;

/*
 * Calcule les valeurs des fonctions de Legendre associees 
 *   P_l^m(cos(theta)) normalisees de facon a ce que
 * 
 *	 int_0^pi [ P_l^m(cos(theta)) ]^2  sin(theta) dtheta = 1
 * 
 * NB: Cette normalisation est differente de celle de la litterature
 *
 * Le calcul est effectue aux 2*nt-1 points  
 *	theta_j = pi/2 j/(2*nt-2)		0 <= j <= 2*nt-2 
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
 * double* legendre_norm : ensemble des (2*nt-1-m)*(2*nt-1) valeurs 
 *			    P_l^m(cos(theta))
 *			stokees comme suit:
 *
 *	    legendre_norm[(2*nt-1)* (l-m) + j] = P_l^m( cos(theta_j) ) 
 *
 *			avec   m <= l <= 2*nt-2.
 *
 * NB: Cette routine effectue le calcul a chaque appel et ne renvoie pas
 *     un pointeur sur des valeurs precedemment calculees.
 */


/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:54  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/02/22  15:37:00  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers du C
#include <assert.h>
#include <math.h>
#include <malloc.h>

// Prototypage
#include "headcpp.h"
#include "proto.h"

//******************************************************************************

double* legendre_norm(int m, int nt) {

int l, j ;

    int lmax = 2*nt - 2 ; 

// Sur-echantillonnage pour calculer les carres sans aliasing:
    int nt2 = 2*nt - 1 ; 
    int nt2m1 = nt2 - 1 ; 

    int deg[3] ;
    deg[0] = 1 ;
    deg[1] = 1 ;
    deg[2] = nt2 ;

// Tableau de travail
    double* yy = (double*)( malloc( nt2*sizeof(double) ) ) ;
    
// Recherche des fonctions de legendre associees non normalisees 
// -------------------------------------------------------------
//  (NB: elles different de celles de la litterature par un facteur (2m-1)!!) :

    double* leg = legendre(m, nt2) ;
    
// Normalisation 
// -------------
     for (l=m; l<lmax+1; l++) {
     
	int ml = (m+1)*(l+1) ; 
	
// On divise les fonctions de Legendre par (m+1)*(l+1) 
// pour obtenir des nombres pas trop grands:

	for (j=0; j<nt2; j++) {
	   leg[nt2*(l-m)+j] /= ml ;
	}   

// Carre : 
	for (j=0; j<nt2; j++) {	
	   double w = leg[nt2*(l-m)+j]  ;
	   yy[nt2m1-j] =  w * w  ;   // le rangement est celui qui convient
				    // a cfrchebp
	}   
	
// Developpement en polynomes de Tchebyshev pairs (x=cos(theta)) : 

	cfrchebp(deg, deg, yy, deg, yy) ;
    
// Integrale sur [0,Pi] = 2 fois l'integrale sur [0,1] pour x = cos(theta) : 
	double integ = 2.*int1d_chebp(nt2, yy) ;
    
// Facteur de normalisation
	double fact = 1. / sqrt(integ) ;
 
/* Test: Comparaison avec le resultat analytique
 * 
 *	double fact_test = ml* factorielle2(2*m-1) * sqrt( double(2*l+1)/2.
 *		* factorielle(l-m) / factorielle(l+m) ) ;
 *	double diff = (fact - fact_test) / fact_test ; 
 *
 *	cout << "m, l : "<< m << " " << l << " : " << fact << " " << fact_test 
 *	    << " " << diff << endl ;
 */
    
	for (j=0; j<nt2; j++) {
	   leg[nt2*(l-m)+j] *= fact ;
	}   
       
     }	    // fin de la boucle sur l
     
// Liberation espace memoire :
    free(yy) ; 

     return leg ; 

}



