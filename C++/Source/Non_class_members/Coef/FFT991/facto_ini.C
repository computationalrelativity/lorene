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


 


/* 
 * Routine d'initialisation des tables de factorisation de n
 * Version speciale FAX
 *
 * Entree:
 *   n		nombre de degres de liberte
 * Sortie: 
 *   facto_ini	pointeur int* sur la table de factorisation de n 
 *
 * Cette routine doit etre entierement en zone critique
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.6  2016/12/05 16:18:04  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/15 12:48:22  j_novak
 * Corrected namespace declaration.
 *
 * Revision 1.4  2014/10/13 08:53:18  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:18:47  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2009/11/25 08:00:54  j_novak
 * Removed a #pragma directive
 *
 * Revision 1.1  2004/12/21 17:06:01  j_novak
 * Added all files for using fftw3.
 *
 * Revision 1.4  2003/01/31 10:31:24  e_gourgoulhon
 * Suppressed the directive #include <malloc.h> for malloc is defined
 * in <stdlib.h>
 *
 * Revision 1.3  2002/10/16 14:36:54  j_novak
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
 * Revision 2.1  1999/11/24  16:21:27  eric
 * Modif affichage.
 *
 * Revision 2.0  1999/02/22  15:38:52  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers du C
#include <cmath>
#include <cstdlib>

// Prototypes of F77 subroutines
#include "headcpp.h"
#include "proto_f77.h"

namespace Lorene {
int *facto_ini( int n )
{

// Variables locales statiques
// ---------------------------
#define NMAX	30			/* Nombre maximun de dimensions differentes */
static	int   	*table_facto[NMAX] ;	/* Tableau des pointeurs sur les tableaux */
static	int	nwork = 0 ;		/* Nombre de tableaux deja initialises */
static	int	tbn[NMAX] ;		/* Tableau des points deja initialises */
static	int trois = 3 ;			// On aurait aime ecrire: "const int trois=3 ;"
int indice ;

{
    // Ce nombre de points a-t-il deja ete utilise ?
    int i ;
    indice = -1 ;
    for ( i=0 ; i < nwork ; i++ ) {
	if ( tbn[i] == n ) indice = i ;
    }

// Initialisation si necessaire
	if (indice == -1) {	/* Il faut une nouvelle initialisation */
		if ( nwork >= NMAX ) {
		    cout << "facto_ini : nwork >= NMAX !" << endl ; 
		    abort() ; 
			}
		indice = nwork ; nwork++ ; tbn[indice] = n ;

		table_facto[indice] = (int *) malloc( sizeof(int) * 20 ) ;
//		table_facto[indice] = new int[20] ;
		if ( table_facto[indice] == 0 ) {
		    cout << "facto_ini : malloc error !" << endl ; 
		    abort() ; 
			}

		F77_fax( table_facto[indice], &n, &trois ) ;
		}

}   // Fin de zone critique

    // Valeurs de retour
    return table_facto[indice] ;

}
}
