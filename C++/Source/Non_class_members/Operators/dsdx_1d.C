/*
 *   Copyright (c) 2004 Jerome Novak
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


char dsdx_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2006/04/10 15:19:20  j_novak
 * New definition of 1D operators dsdx and sx in the nucleus (bases R_CHEBP and
 * R_CHEBI).
 *
 * Revision 1.2  2005/01/10 16:34:53  j_novak
 * New class for 1D mono-domain differential operators.
 *
 * Revision 1.1  2004/02/06 10:53:53  j_novak
 * New dzpuis = 5 -> dzpuis = 3 case (not ready yet).
 *
 *
 * $Header$
 *
 */

#include <stdlib.h>
#include "type_parite.h"
#include "headcpp.h"
#include "proto.h"

/*
 * Routine appliquant l'operateur dsdx.
 * 
 * Entree : tb contient les coefficients du developpement 
 *	    int nr : nombre de points en r. 
 *
 * Sortie : tb contient dsdx
 * 
 */

		//----------------------------------
		// Routine pour les cas non prevus --
		//----------------------------------

void _dsdx_1d_pas_prevu(int nr, double* tb, double *xo) {
    cout << "dsdx pas prevu..." << endl ;
    cout << "Nombre de points : " << nr << endl ;
    cout << "Valeurs : " << tb << "  " << xo <<endl ;
    abort() ;
    exit(-1) ;
}

			//----------------
			// cas R_CHEBU ---
			//----------------

void _dsdx_1d_r_chebu(int nr,  double* tb, double *xo)
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
			// cas R_CHEBI ---
			//----------------

void _dsdx_1d_r_chebi(int nr,  double* tb, double *xo)
{

    double som ;
	    
    xo[nr-1] = 0 ;
    som = 2*(2*nr-3) * tb[nr-2] ;
    xo[nr-2] = som ;
    for (int i = nr-3 ; i >= 0 ; i -- ) {
      som += 2*(2*i+1) * tb[i] ;
      xo[i] = som ;
    }	
    xo[0] *= .5 ;
}

			//----------------
			// cas R_CHEBP ---
			//----------------

void _dsdx_1d_r_chebp(int nr,  double* tb, double *xo)
{

    double som ;
	    
    xo[nr-1] = 0 ;
    som = 4*(nr-1) * tb[nr-1] ;
    xo[nr-2] = som ;
    for (int i = nr-3 ; i >= 0 ; i --) {
      som += 4*(i+1) * tb[i+1] ;
      xo[i] = som ;
    }	
}

		// ---------------------
		// La routine a appeler
		//----------------------
		
		
void dsdx_1d(int nr, double** tb, int base_r)
{

		// Routines de derivation
    static void (*dsdx_1d[MAX_BASE])(int, double*, double *) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    dsdx_1d[i] = _dsdx_1d_pas_prevu ;
	}
		// Les routines existantes
	dsdx_1d[R_CHEBU >> TRA_R] = _dsdx_1d_r_chebu ;
	dsdx_1d[R_CHEBP >> TRA_R] = _dsdx_1d_r_chebp ;
	dsdx_1d[R_CHEBI >> TRA_R] = _dsdx_1d_r_chebi ;
	dsdx_1d[R_CHEB >> TRA_R] = _dsdx_1d_r_cheb ;

    }
    
    double *result = new double[nr] ;
    
    dsdx_1d[base_r](nr, *tb, result) ;
    
    delete [] (*tb) ;
    (*tb) = result ;
}
