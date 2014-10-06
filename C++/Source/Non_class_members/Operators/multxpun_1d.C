/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


char multxpun_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:16:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2007/12/11 15:42:22  jl_cornou
 * Premiere version des fonctions liees aux polynomes de Jacobi(0,2)
 *
 * Revision 1.2  2002/10/16 14:37:07  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/04/03  17:01:59  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


#include <cstdlib>
#include <cmath>

#include "tbl.h"
#include "type_parite.h"

/*
 * Operateur :
 *  -R_JACO02 : (f-f(-1))*(x+1)
 * 
 *
 * Entree : coefficients de f dans tb
 *	    nr : nombre de points en r
 *
 * Sortie : coefficient du resultat dans xo
 * 
 * 
 */


		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

void _multxpun_1d_pas_prevu(int nr, double* tb, double *res) {
    cout << "multxpun pas prevu..." << endl ;
    cout << " valeurs: " << tb << "   " << res << endl ;
    cout << "nr : " << nr << endl ;
    abort () ;
    exit(-1) ;
}

			//---------------
			// cas R_JACO02 -
			//---------------

void _multxpun_1d_r_jaco02 (int nr, double* tb, double *xo)
{
    
    xo[nr-1] = 0 ;
    for (int i = 1 ; i < nr-1 ; i++) {
	xo[i] = i*(i+2)/double((i+1)*(2*i+1))*tb[i-1] + (i*i+3*i+3)/double((i+1)*(i+2))*tb[i] + (i+1)*(i+3)/double((i+2)*(2*i+5))*tb[i+1] ;
    }
    xo[0] = 1.5*tb[0] + 0.3*tb[1] ;
    xo[nr-1] = (nr*nr-1)/double((nr)*(2*nr-1))*tb[nr-2] + (1+1/double((nr)*(nr+1)))*tb[nr-1] ;
}

		    //----------------------------
		    // La routine a appeler   ----
		    //----------------------------
		    
		    
void multxpun_1d(int nr, double **tb, int base_r)	    // Version appliquee a this
{

// Routines de derivation
static void (*multxpun_1d[MAX_BASE])(int, double *, double *) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    multxpun_1d[i] = _multxpun_1d_pas_prevu ;
	}
	// Les routines existantes
	multxpun_1d[R_JACO02 >> TRA_R] = _multxpun_1d_r_jaco02 ;
	}
    
    double *result = new double[nr] ;
    multxpun_1d[base_r](nr, *tb, result) ;
    
    delete [] (*tb) ;
    (*tb) = result ;
}
