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


char multx_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:16:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2002/10/16 14:36:58  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  1999/07/08  09:54:30  phil
 * correction gestion memoire
 *
 * Revision 2.0  1999/07/07  10:15:40  phil
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


/*
 * Operateurs de multiplication par x
 * 
 * Uniquement le cas R_CHEB
 * 
 * Entree :
 * 
 *	int nr : nombre de points en r
 *	tb contient les coefficients du developpement
 *	evenetuellement e dans echelle
 *	
 *   Sortie :
 *	tb contient les coefficients du resultat...
 */
 
 
		//--------------------------------------
		// Routine pour les cas non prevus -----
		//--------------------------------------

void _multx_1d_pas_prevu(int nr, double* tb, double *result) {
    cout << "multx pas prevu..." << endl ;
    cout << "Valeurs : " << tb << "  " << result << endl ;
    cout << " nr : " << nr << endl ;
    abort() ;
    exit(-1) ;
}


		//------------------
		// cas R_CHEB	 ---
		//------------------

void _multx_1d_r_cheb (int nr, double* tb, double* res) {
    
    res[0] = tb[1]/2. ;
    res[1] = (2*tb[0]+tb[2])/2. ;
    res[nr-1] = tb[nr-2]/2. ;
    
    for (int i=2 ; i<nr-1 ; i++)
	res[i] = (tb[i-1]+tb[i+1])/2. ;
	
}



		// ----------------------
		// La routine a appeler 
		//-----------------------
		
void multx_1d(int nr,  double **tb, int base_r)
{

		// Routines de derivation
    static void (*multx_1d[MAX_BASE])(int, double *, double *) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    multx_1d[i] = _multx_1d_pas_prevu ;
	}
		// Les routines existantes
	multx_1d[R_CHEB >> TRA_R] = _multx_1d_r_cheb ;
    }
    
    
    double *result = new double[nr] ;
    multx_1d[base_r](nr, *tb, result) ;
    delete [] (*tb) ;
    (*tb) = result ;
}

