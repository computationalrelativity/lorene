/*
 *   Copyright (c) 2003 Jerome Novak
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


char cftlegii_C[] = "$Header$" ;

/*
 * Transformation en fonctions de Legendre associees sur le deuxieme indice
 *  (theta) d'un tableau 3-D representant une fonction 3-D symetrique
 *  par rapport au plan equatorial z = 0 et antisymetrique par le retournement 
 *  (x, y, z) --> (-x, -y, z). 
 *
 *
 * Entree:
 * -------
 *   int* deg	: tableau du nombre effectif de degres de liberte dans chacune 
 *		  des 3 dimensions: le nombre de points de collocation
 *		  en theta est  nt = deg[1] et doit etre de la forme
 * 			nt = 2^p 3^q 5^r + 1 
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
 *		    f( theta_l ) = ff[ dimf[1]*dimf[2] * k + i + dimf[2] * l ]
 *			 ou k et i sont les indices correspondant a
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
 *   double* cf :  tableau des coefficients a_j du develop. en fonctions de
 *		    Legendre associees P_l^m (l pair, m impair)
 *
 *		    f(theta) = 
 *			    som_{l=(m-1)/2}^{nt-2} a_j P_{2j}^m( cos(theta) )
 *
 *		    avec m impair.		  
 *
 *		    P_l^m(x) represente la fonction de Legendre associee
 *		       de degre l et d'ordre m normalisee de facon a ce que
 *
 *			int_0^pi [ P_l^m(cos(theta)) ]^2  sin(theta) dtheta = 1
 *
 * 		    L'espace memoire correspondant au pointeur cfo doit etre 
 *	            nr*nt*(np+2) et doit avoir ete alloue avant 
 *		    l'appel a la routine.	 
 *		    Le coefficient a_j (0 <= j <= nt-1) est stoke dans le 
 *		    tableau cfo comme suit
 *		          a_j = cfo[ nr*nt* k + i + nr* j ]
 *		    ou k et i sont les indices correspondant a phi et r 
 *		    respectivement: m = 2( k/2 ).
 *		    NB: pour j<(m+1)/2,  a_j = 0
 *
 *
 * NB: Si le pointeur cf est egal a ff, la routine ne travaille que sur un 
 *     seul tableau, qui constitue une entree/sortie.
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/09/16 13:07:41  j_novak
 * New files for coefficient trnasformation to/from the T_LEG_II base.
 *
 *
 * $Header$
 *
 */

// headers du C
#include <assert.h>

// headers bien de chez nous
#include "headcpp.h"
#include "proto.h"
//*****************************************************************************

void cftlegii(const int* deg, const int* dimf, double* ff, const int* dimc,
		double* cf)
{

// Limitations de la routine:
    assert(dimc[0]==deg[0]+2) ;
    assert(dimc[1]==deg[1]) ;
    assert(dimc[2]==deg[2]) ;


    // Tableau de travail :
    int taille = dimc[0]*dimc[1]*dimc[2] ;
    double* cf_cs =  new double[taille] ;	

//--------------------------------------------------------------
// 1/ Transformation esp. des configurations --> sin(2j theta) 
//--------------------------------------------------------------

    cftsinp(deg, dimf, ff, dimc, cf_cs) ;

//--------------------------------------------------------------
// 2/ Transformation  sin(2j theta)  ---> Legendre 
//--------------------------------------------------------------

    chb_sinp_legii(deg , cf_cs, cf) ;

    // Menage
    delete [] cf_cs ;
}
