/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char sxdsdx_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:37:00  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  1999/07/08  09:55:37  phil
 * correction gestion memoire
 *
 * Revision 2.0  1999/07/07  10:15:14  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Includes
#include <stdlib.h>

#include "tbl.h"
#include "type_parite.h"

/*
 * Operateur :
 *  -R_CHEB : ds/dx 
 *  -R_CHEBP ou R_CHEBI : (f'-f'(0))/x
 * 
 * 
 * 
 * Entree : coefficients de f dans tb
 *	    nr : nombre de points en r
 * Sortie : coeffieicient du resultat dans tb
 * 
 * 
 */


		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

void _sxdsdx_1d_pas_prevu(int nr, double* tb, double *res) {
    cout << "sxdsdx pas prevu..." << endl ;
    cout << " valeurs: " << tb << "   " << res << endl ;
    cout << "nr : " << nr << endl ;
    abort () ;
    exit(-1) ;
}



			//---------------
			// cas R_CHEB ---
			//---------------

void _dsdx_1d_r_cheb(int nr, double* tb, double *xo)
{

    double som ;
	    
    xo[nr-1] = 0 ;
    som = 2*(nr-1) * tb[nr-1] ;
    xo[nr-2] = som ;
    for (int i = nr-4 ; i >= 0 ; i -= 2 ) {
	som += 2*(i+1) * tb[i+1] ;
	xo[i] = som ;
    }	// Fin de la premiere boucle sur r
    som = 2*(nr-2) * tb[nr-2] ;
    xo[nr-3] = som ;
    for (int i = nr-5 ; i >= 0 ; i -= 2 ) {
	som += 2*(i+1) * tb[i+1] ;
	xo[i] = som ;
    }	// Fin de la deuxieme boucle sur r
    xo[0] *= .5 ;

}

			//----------------
			// cas R_CHEBP ---
			//----------------

void _sxdsdx_1d_r_chebp(int nr, double* tb, double *xo)
{

    double som ;
      
    xo[nr-1] = 0 ;
    som = 8 * (nr-1) * tb[nr-1] ;
    xo[nr-2] = som ;
    for (int i = nr-4 ; i >= 0 ; i -= 2 ) {
	som += 8 * (i+1) * tb[i+1] ;
	xo[i] = som ;
    }	// Fin de la premiere boucle sur r

    som = 8 * (nr-2) * tb[nr-2] ;
    xo[nr-3] = som ;
    for (int i = nr-5 ; i >= 0 ; i -= 2 ) {
	som += 8 * (i+1) * tb[i+1] ;
	xo[i] = som ;
    }	// Fin de la deuxieme boucle sur r
	    
    xo[0] *= .5 ;
	
}

			//----------------
			// cas R_CHEBI ---
			//----------------
			
void _sxdsdx_1d_r_chebi(int nr, double* tb, double *xo)
{

 
    double som ;

    xo[nr-1] = 0 ;
    som = 4 * (2*(nr-1)+1) * tb[nr-1] ;
    xo[nr-2] = som ;
    for (int i = nr-4 ; i >= 0 ; i -= 2 ) {
	som += 4 * (2*(i+1)+1) * tb[i+1] ;
	xo[i] = som ;
    }	// Fin de la premiere boucle sur r

    som = 4 * (2*(nr-2)+1) * tb[nr-2] ;
    xo[nr-3] = som ;
    for (int i = nr-5 ; i >= 0 ; i -= 2 ) {
	som += 4 * (2*(i+1)+1) * tb[i+1] ;
	xo[i] = som ;
    }	// Fin de la deuxieme boucle sur r

}


		    //----------------------------
		    // La routine a appeler   ----
		    //----------------------------
		    
		    
void sxdsdx_1d(int nr, double **tb, int base_r)	    // Version appliquee a this
{

// Routines de derivation
static void (*sxdsdx_1d[MAX_BASE])(int, double *, double *) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    sxdsdx_1d[i] = _sxdsdx_1d_pas_prevu ;
	}
	// Les routines existantes
	sxdsdx_1d[R_CHEB >> TRA_R] = _dsdx_1d_r_cheb ;
	sxdsdx_1d[R_CHEBP >> TRA_R] = _sxdsdx_1d_r_chebp ;
	sxdsdx_1d[R_CHEBI >> TRA_R] = _sxdsdx_1d_r_chebi ;
	}
    
    double *result = new double[nr] ;
    sxdsdx_1d[base_r](nr, *tb, result) ;
    
    delete [] (*tb) ;
    (*tb) = result ;
}
