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


char sxpundsdx_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:27  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:16:07  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2007/12/11 15:42:23  jl_cornou
 * Premiere version des fonctions liees aux polynomes de Jacobi(0,2)
 *
 * Revision 1.2  2002/10/16 14:37:11  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/10/11  09:55:46  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */
 
 // Includes
#include <cstdlib>

#include "headcpp.h"
#include "type_parite.h"

namespace Lorene {
void dsdx_1d(int, double**, int) ;
void sxpun_1d(int, double**, int) ;


		//----------------------------------
		// Routine pour les cas non prevus --
		//----------------------------------

void _sxpundsdx_1d_pas_prevu(int nr, double* tb, double *xo) {
    cout << "sxpundsdx pas prevu..." << endl ;
    cout << "Nombre de points : " << nr << endl ;
    cout << "Valeurs : " << tb << "  " << xo <<endl ;
    abort() ;
    exit(-1) ;
}


			//---------------
			// cas R_JACO02 -
			//---------------

void _sxpundsdx_1d_r_jaco02(int nr, double* tb, double *xo) {
    
	dsdx_1d(nr, &tb, R_JACO02 >> TRA_R) ;
	sxpun_1d(nr, &tb, R_JACO02 >> TRA_R) ;
	for (int j=0 ; j<nr ; j++) {
		xo[j]=tb[j] ;
	}
}

		// ---------------------
		// La routine a appeler
		//----------------------
		
		
void sxpundsdx_1d(int nr, double** tb, int base_r)
{

		// Routines de derivation
    static void (*sxpundsdx_1d[MAX_BASE])(int, double*, double *) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    sxpundsdx_1d[i] = _sxpundsdx_1d_pas_prevu ;
	}
		// Les routines existantes
	sxpundsdx_1d[R_JACO02 >> TRA_R] = _sxpundsdx_1d_r_jaco02 ;
    }
    
    double *result = new double[nr] ;
    
    sxpundsdx_1d[base_r](nr, *tb, result) ;
    
    delete [] (*tb) ;
    (*tb) = result ;
}
}
