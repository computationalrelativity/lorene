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


char cheb_ini_C[] = "$Header$" ;

/* 
 * Routine de calcul des sin(psi) pour les transformations de Tchebyshev
 *   (psi decrivant uniformement l'intervalle [0, pi]). 
 *
 * Entree:
 *   n		nombre de degres de liberte
 * Sortie:
 *   chebf_ini	pointeur double* sur la table des sinus
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:53  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  1999/11/24  16:16:25  eric
 * Modif affichage.
 *
 * Revision 2.0  1999/02/22  15:44:22  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers du C
#include <math.h>
#include <stdlib.h>
#include <malloc.h>

#include "headcpp.h"

// Variables externes de loch
int loch_cheb_ini = 0 ;

//---------------------------------------------------------------------------------

double* cheb_ini(const int n )
{

// Variables locales statiques
// ---------------------------
#define NMAX	30			/* Nombre maximun de dimensions differentes */
static	double*	table_sin[NMAX] ;	/* Tableau des pointeurs sur les tableaux */
static	int	nwork = 0 ;		/* Nombre de tableaux deja initialises */
static	int	tbn[NMAX] ;		/* Tableau des points deja initialises */
int indice ;

// Cette routine est entierement critique
#pragma critical (loch_cheb_ini)
{

    // Ce nombre de points a-t-il deja ete utilise ?
    indice = -1 ;
    int i ;
    for ( i=0 ; i < nwork ; i++ ) {
	if ( tbn[i] == n ) indice = i ;
	}

    // Initialisation
    if (indice == -1) {		    /* Il faut une nouvelle initialisation */
	if ( nwork >= NMAX ) {
	    cout << "cheb_ini : nwork >= NMAX !" << endl ; 
	    abort() ; 
	}
	indice = nwork ; nwork++ ; tbn[indice] = n ;

	int nm1s2 = (n-1) / 2 ;  		
	table_sin[indice] = (double *) malloc( sizeof(double) * nm1s2 ) ;
	if ( table_sin[indice] == 0 ) {
	    cout << "cheb_ini : malloc error !" << endl ; 
	    abort() ; 
	}

	double xx = M_PI / double(n-1);
	for ( i = 0; i < nm1s2 ; i++ ) {
	    table_sin[indice][i] = sin( xx * i );
	    }
	}
    }   // Fin de la region critique

    // Valeurs de retour
    return table_sin[indice] ;
	
}
