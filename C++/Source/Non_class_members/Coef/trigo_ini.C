/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
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


char trigo_ini_C[] = "$Header$" ;

/* 
 * Routine d'initialisation des tables trigo
 * Version speciale FAX
 *
 * Entree:
 *   n		nombre de degres de liberte 
 * Sortie:
 *   trigo_ini	pointeur double* sur la table trigo 
 *
 * Doit etre en zone critique
 */

/*
 * $Id$
 * $Log$
 * Revision 1.5  2003/12/19 16:21:47  j_novak
 * Shadow hunt
 *
 * Revision 1.4  2003/01/31 10:31:24  e_gourgoulhon
 * Suppressed the directive #include <malloc.h> for malloc is defined
 * in <stdlib.h>
 *
 * Revision 1.3  2002/10/16 14:36:57  j_novak
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
 * Revision 2.1  1999/11/24  16:23:22  eric
 * Modif affichage.
 *
 * Revision 2.0  1999/02/22  15:29:33  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers du C
#include <math.h>
#include <stdlib.h>

// Prototypes of F77 subroutines
#include "headcpp.h"
#include "proto_f77.h"

// Variable de loch
int loch_trigo_ini = 0 ;

double *trigo_ini( int n )
{
// Variables locales statiques
// ---------------------------
#define NMAX	30			/* Nombre maximun de dimensions differentes */
static	double	*table_trigo[NMAX] ;	/* Tableau des pointeurs sur les tableaux */
static	int	nwork = 0 ;		/* Nombre de tableaux deja initialises */
static	int	tbn[NMAX] ;		/* Tableau des points deja initialises */
static	int trois = 3 ;
int indice ;

//#pragma critical (loch_trigo_ini)
{
    // Ce nombre de points a-t-il deja ete utilise ?
    indice = -1 ;
    int i ;
    for ( i=0 ; i < nwork ; i++ ) {
	if ( tbn[i] == n ) indice = i ;
    }

    // Initialisation
	if (indice == -1) {	/* Il faut une nouvelle initialisation */
		if ( nwork >= NMAX ) {
		    cout << "trigo_ini : nwork >= NMAX !" << endl ; 
		    abort() ; 
			}
		indice = nwork ; nwork++ ; tbn[indice] = n ;

//		table_trigo[indice] = new double[3*n/2 + 1] ;
		table_trigo[indice] = (double *) malloc( sizeof(double) * (3*n/2 + 1) ) ;
		if ( table_trigo[indice] == 0 ) {
		    cout << "trigo_ini : malloc error !" << endl ; 
		    abort() ; 
			}

		F77_fftrig( table_trigo[indice], &n, &trois ) ;
		}

    }	// Fin de zone critique
    
    // Valeurs de retour
    return table_trigo[indice] ;
}
