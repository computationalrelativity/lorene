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


char mat_legpi_sini_C[] = "$Header$" ;

/*
 * Fournit la matrice de passage pour la transformation des coefficients du
 * developpement en fonctions associees de Legendre 
 * P_l^m(cos(theta)) avec l impair et m impair
 * dans les coefficients du developpement 
 * en sin( (2j+1) theta ).
 * 
 * Cette routine n'effectue le calcul de la matrice que si celui-ci n'a pas 
 * deja ete fait, sinon elle renvoie le pointeur sur une valeur precedemment 
 * calculee.
 * 
 * Entree:
 * -------
 *  int np  :	Nombre de degres de liberte en phi 
 *  int nt  :	Nombre de degres de liberte en theta
 *
 * Sortie (valeur de retour) :
 * ---------------------------
 *  double* mat_legpi_sini : pointeur sur le tableau contenant l'ensemble
 *			        (pour les np/2 valeurs de m: m=1,3,...,np-1)) des 
 *				matrices de passage.
 *				La dimension du tableau est (np/2+1)*nt^2
 *				Le stokage est le suivant: 
 *
 *		mat_legpi_sini[ nt*nt* m/2 + nt*j + l] = B_{mjl} 
 *			    
 *				ou B_{mjl} (m impair) est defini par
 *
 *   P_{2l+1}^m( cos(theta) ) = som_{j=0}^{nt-2} B_{mjl} sin((2*j+1)*theta) 
 *
 *  ou P_n^m(x) represente la fonction de Legendre associee de degre n et 
 *     d'ordre m normalisee de facon a ce que
 *
 *	 int_0^pi [ P_n^m(cos(theta)) ]^2  sin(theta) dtheta = 1
 *
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.1  2000/11/14  15:12:36  eric
 * Traitement du cas np=1
 *
 * Revision 2.0  2000/09/29  16:09:41  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers du C++
#include <iostream.h>

// headers du C
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

// Headers Lorene
#include "proto.h"

// Variable de loch
int loch_mat_legpi_sini = 0 ;

//******************************************************************************

double* mat_legpi_sini(int np, int nt) {

#define NMAX	30		// Nombre maximun de couples(np,nt) differents 
static	double*	tab[NMAX] ;	// Tableau des pointeurs sur les tableaux 
static	int	nb_dejafait = 0 ;	// Nombre de tableaux deja initialises 
static	int	np_dejafait[NMAX] ;    // Valeurs de np pour lesquelles le
				       // calcul a deja ete fait
static	int	nt_dejafait[NMAX] ;    // Valeurs de np pour lesquelles le
				       // calcul a deja ete fait

int i, indice,  j,  j2,  m,  l ;

    #pragma critical (loch_mat_legpi_sini)
    {
// Les matrices B_{mjl} pour ce couple (np,nt) ont-elles deja ete calculees ? 

    indice = -1 ;
    for ( i=0 ; i < nb_dejafait ; i++ ) {
	if ( (np_dejafait[i] == np) && (nt_dejafait[i] == nt) ) indice = i ;
	}


// Si le calcul n'a pas deja ete fait, il faut le faire : 
    if (indice == -1) {		   
	if ( nb_dejafait >= NMAX ) {
	    cout << "mat_legpp_cosp: nb_dejafait >= NMAX : "
		<< nb_dejafait << " <-> " << NMAX << endl ;
	    abort () ;	
	    exit(-1) ;
	    }
	indice = nb_dejafait ; 
	nb_dejafait++ ; 
	np_dejafait[indice] = np ;
	nt_dejafait[indice] = nt ;

	tab[indice] = (double *) malloc( sizeof(double) * (np/2+1)*nt*nt ) ;

//-----------------------
// Preparation du calcul 
//-----------------------

// Sur-echantillonnage :
	int nt2 = 2*nt - 1 ; 

	int deg[3] ;
	deg[0] = 1 ;
	deg[1] = nt2 ;
	deg[2] = 1 ;

// Tableaux de travail
	double* yy = (double*)( malloc( nt2*sizeof(double) ) ) ;


//-------------------
// Boucle sur m
//-------------------

	int m_max = (np == 1) ? 1 : np-1 ; 

	for (m=1; m <= m_max ; m+=2) {

	    // Recherche des fonctions de Legendre associees d'ordre m :

	    double* leg = legendre_norm(m, nt) ;	
	

	    for (l=(m-1)/2; l<nt-1; l++) {	// boucle sur les P_{2l+1}^m
	    
		int ll = 2*l+1 ;	// degre des fonctions de Legendre
		    	
		for (j2=0; j2<nt2; j2++) {
		    yy[j2] = leg[nt2* (ll-m) + j2] ;
		}

		//....... transformation en sin((2j+1)*theta) : 

		cftsini(deg, deg, yy, deg, yy) ; 

		//....... le resultat fournit les elements de matrice : 
		for (j=0; j<nt-1; j++) {
		    tab[indice][ nt*nt* ((m-1)/2) + nt*j + l] = yy[j] ; 
		}
	    
	    }  // fin de la boucle sur l (indice de P_{2l+1}^m)
	    
	    free(leg) ;
	   	    
	}  // fin de la boucle sur m

// Liberation espace memoire
// -------------------------

	free(yy) ;

    } // fin du cas ou le calcul etait necessaire

    }	// Fin de zone critique
    
    return tab[indice] ;

}


