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


char sxp12_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:16:07  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2007/12/11 15:42:23  jl_cornou
 * Premiere version des fonctions liees aux polynomes de Jacobi(0,2)
 *
 * Revision 1.2  2002/10/16 14:36:58  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/10/11  09:53:08  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Includes
#include <cstdlib>
#include <cassert>

#include "headcpp.h"
#include "type_parite.h"
#include "proto.h"

void sxpun_1d(int, double **, int) ;

		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

void _sxp12_1d_pas_prevu(int nr, double* tb, double *res) {
    cout << "sxp12 pas prevu..." << tb << "    " << res << endl ;
    cout << "nr : " << nr << endl ;
    abort() ;
    exit(-1) ;
}

			//---------------
			// cas R_JACO02 -
			//---------------

void _sxp12_1d_r_jaco02(int nr, double* tb, double *xo) {
 
    assert (nr>2) ;
    sxpun_1d(nr, &tb, R_JACO02 >> TRA_R) ;
    sxpun_1d(nr, &tb, R_JACO02 >> TRA_R) ;
	for (int i = 0 ; i<nr ; i++) {
	xo[i] = tb[i] ;
	}
}


		    //----------------------
		    // La routine a appeler
		    //----------------------
		    
void sxp12_1d(int nr, double **tb, int base_r)	    // Version appliquee a this
{

// Routines de derivation
static void (*sxp12_1d[MAX_BASE])(int, double *, double*) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    sxp12_1d[i] = _sxp12_1d_pas_prevu ;
	}
	// Les routines existantes
	sxp12_1d[R_JACO02 >> TRA_R] = _sxp12_1d_r_jaco02 ;
    }
    
    double *result = new double[nr] ;
    sxp12_1d[base_r](nr, *tb, result) ;
    
    delete [] (*tb) ;
    (*tb) = result ;
}		
