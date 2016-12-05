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


 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:18:09  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2015/03/05 08:49:32  j_novak
 * Implemented operators with Legendre bases.
 *
 * Revision 1.3  2014/10/13 08:53:27  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:16:07  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2006/04/10 15:19:20  j_novak
 * New definition of 1D operators dsdx and sx in the nucleus (bases R_CHEBP and
 * R_CHEBI).
 *
 *
 * $Header$
 *
 */
 
 
 // Includes
#include <cstdlib>

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

namespace Lorene {
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

		//------------------
		// case R_LEGP	 ---
		//------------------

void _sx_1d_r_legp (int nr, double* tb, double* res) {
   
    double term = 0 ;
    res[nr-1] = 0 ;
    for (int i=nr-2 ; i>=0 ; i--) {
	term += tb[i+1] ;
	res[i] = double(4*i+3)/double(2*i+2)*term ;
	term *= -double(2*i+1)/double(2*i+2) ;
    }
 
}


		//------------------
		// case R_LEGI	 ---
		//------------------

void _sx_1d_r_legi (int nr, double* tb, double* res) {
 
  double term = 0 ;
  res[nr-1] = 0 ;
  for (int i=nr-2; i>=0; i--) {
    term += tb[i] ;
    res[i] = double(4*i+1)/double(2*i+1)*term ;
    term *= -double(2*i)/double(2*i+1) ;
  }
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
	sx_1d[R_LEGP >> TRA_R] = _sx_1d_r_legp ;
	sx_1d[R_LEGI >> TRA_R] = _sx_1d_r_legi ;
	sx_1d[R_CHEBU >> TRA_R] = _sx_1d_r_chebu ;
    }
    
    
    double *result = new double[nr] ;
    sx_1d[base_r](nr, *tb, result) ;
    delete [] (*tb) ;
    (*tb) = result ;
}

}
