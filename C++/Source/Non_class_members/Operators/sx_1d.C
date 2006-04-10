/*
 *   Copyright (c) 2006 Jerome Novak
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


char sx_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2006/04/10 15:19:20  j_novak
 * New definition of 1D operators dsdx and sx in the nucleus (bases R_CHEBP and
 * R_CHEBI).
 *
 *
 * $Header$
 *
 */
 
 
 // Includes
#include <stdlib.h>

#include "headcpp.h"
#include "type_parite.h"
#include "proto.h"


/*
 * 1/x operator (division by x)
 * 
 * Only for the bases R_CHEBP and R_CHEBI
 * 
 * Input :
 * 
 *	int nr : number of coefficients
 *	tb     : the array of coefficients of the input function
 *	
 * Output :
 *	tb     : the array of coefficients of the result
 */
 
 
		//----------------------------
		// Bases not implemented -----
		//----------------------------

void _sx_1d_pas_prevu(int , double* , double *) {
    cout << "sx_1d : base not implemented..." << endl ;
    abort() ;
    exit(-1) ;
}


		//-------------------
		// case R_CHEBP	 ---
		//-------------------

void _sx_1d_r_chebp (int nr, double* tb, double* res) {
   
    double som ;
    int sign = 1 ; 
    res[nr-1] = 0 ;
    som = 2 * sign * tb[nr-1] ;
    res[nr-2] = som ;
    for (int i=nr-3 ; i>=0 ; i--) {
	sign = - sign ;
	som += 2 * sign * tb[i+1] ;
	res[i] = (i%2 == 0 ? -1 : 1) * som ;
    }
 
}


		//-------------------
		// case R_CHEBI	 ---
		//-------------------

void _sx_1d_r_chebi (int nr, double* tb, double* res) {
   
    double som ;
    int sign = 1 ; 
    res[nr-1] = 0 ;
    som = 2 * sign * tb[nr-2] ;
    res[nr-2] = som ;
    for (int i=nr-3 ; i>=0 ; i--) {
	sign = - sign ;
	som += 2 * sign * tb[i] ;
	res[i] = (i%2 == 0 ? -1 : 1) * som ;
    }
    res[0] *= 0.5 ;
}

		//-------------------
		// case R_CHEBU	 ---
		//-------------------

void _sx_1d_r_chebu (int nr, double* tb, double* res) {
   
    for (int i=0; i<nr; i++)
	res[i] = tb[i] ;

    sxm1_1d_cheb(nr, res) ;

}



		// ----------------------
		// The calling function 
		//-----------------------
		
void sx_1d(int nr,  double **tb, int base_r)
{

    static void (*sx_1d[MAX_BASE])(int, double *, double *) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    sx_1d[i] = _sx_1d_pas_prevu ;
	}
		// Les routines existantes
	sx_1d[R_CHEBP >> TRA_R] = _sx_1d_r_chebp ;
	sx_1d[R_CHEBI >> TRA_R] = _sx_1d_r_chebi ;
	sx_1d[R_CHEBU >> TRA_R] = _sx_1d_r_chebu ;
    }
    
    
    double *result = new double[nr] ;
    sx_1d[base_r](nr, *tb, result) ;
    delete [] (*tb) ;
    (*tb) = result ;
}

