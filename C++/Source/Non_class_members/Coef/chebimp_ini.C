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


char chebimp_ini_C[] = "$Header$" ;


/* 
 * Routine de calcul des points de collocations 
 *
 *			  x_i = sin( pi/2 i/(n-1) )    0 <= i <= n-1  
 *
 * des echantillonnages de [0, 1] rarefies en 0.  
 * Ces valeurs sont notamment utilisees pour les transformations de Tchebyshev 
 * impaires. 
 * 
 *
 * Entree:
 * -------
 *   const int n :	nombre de degres de liberte
 *
 * Sortie:
 * -------
 *   double* chebimp_ini :  pointeur sur la table des points de collocation
 *			    L'espace memoire correspondant a ce pointeur
 *			    est alloue par la routine.
 */

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.1  1999/11/24  16:18:29  eric
 * Modif affichage.
 *
 * Revision 2.0  1999/02/22  15:44:12  hyc
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
#include <stdlib.h>
#include <malloc.h>

// Variable externe de loch
int loch_chebimp_ini = 0 ;

//*****************************************************************************

double* chebimp_ini(const int n )
{

// Variables locales statiques
// ---------------------------
#define NMAX	30			/* Nombre maximun de dimensions differentes */
static	double*	table_x[NMAX] ;		/* Tableau des pointeurs sur les tableaux */
static	int	nwork = 0 ;		/* Nombre de tableaux deja initialises */
static	int	tbn[NMAX] ;		/* Tableau des points deja initialises */
int indice ;

    // Mise en zone critique de toute la routine
    #pragma critical (loch_chebimp_ini)
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
	    cout << "chebimp_ini: nwork > NMAX : "
		<< nwork << " <-> " << NMAX << endl ;
	    abort () ;	
	    }
	indice = nwork ; nwork++ ; tbn[indice] = n ;

	table_x[indice] = (double *) malloc( sizeof(double) * n ) ;
	if ( table_x[indice] == 0 ) {
	    cout << "chebimp_ini : malloc error !" << endl ; 
	    abort() ; 
	    }

	double xx = M_PI / double(2*(n-1));
	for ( i = 0; i < n ; i++ ) {
	    table_x[indice][i] = sin( xx * i );
	    }
	}

    }	    // Fin de la zone critique

    // Valeurs de retour
    return table_x[indice] ;

}
