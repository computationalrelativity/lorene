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


char multx2_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2015/03/05 08:49:32  j_novak
 * Implemented operators with Legendre bases.
 *
 * Revision 1.4  2014/10/13 08:53:24  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
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

		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

namespace Lorene {
void _multx2_1d_pas_prevu(int nr, double* tb, double *res) {
    cout << "multx2 pas prevu..." << tb << "    " << res << endl ;
    cout << "nr : " << nr << endl ;
    abort() ;
    exit(-1) ;
}

			//---------------
			// cas R_CHEBP --
			//---------------

void _multx2_1d_r_chebp(int nr, double* tb, double *xo) {
    
    assert (nr>2) ;
    
    xo[0] = (2*tb[0]+tb[1])/4 ;
    xo[1] = (2*tb[0]+2*tb[1]+tb[2])/4 ;
    
    for (int i=2 ; i<nr-1 ; i++)
	xo[i] = (tb[i-1]+2*tb[i]+tb[i+1])/4 ;

    xo[nr-1] = (tb[nr-2]+2*tb[nr-1])/4 ;
}


			//---------------
			// cas R_CHEBI --
			//---------------

void _multx2_1d_r_chebi(int nr, double* tb, double *xo){

    assert(nr>2) ;

    xo[0] = (3*tb[0]+tb[1])/4 ;

    for (int i=1 ; i<nr-1 ; i++)
	xo[i] = (tb[i-1]+2*tb[i]+tb[i+1])/4 ;

    xo[nr-1] = 0. ;//(tb[nr-2]+2*tb[nr-1])/4 ;
}

			//---------------
			// cas R_CHEB --
			//---------------

void _multx2_1d_r_cheb(int nr, double* tb, double *xo){

    assert(nr>3) ;

    xo[0] = (2*tb[0]+tb[2])/4 ;
    xo[1] = (3*tb[1]+tb[3])/4 ;
    xo[2] = (2*tb[0]+2*tb[2]+tb[4])/4 ;

    for (int i=3 ; i<nr-2 ; i++)
	xo[i] = (tb[i-2]+2*tb[i]+tb[i+2])/4 ;

    for (int i=nr-2 ; i<nr ; i++)
	xo[i] = (tb[i-2]+2*tb[i])/4 ;
}


			//--------------
			// cas R_LEGP --
			//--------------

void _multx2_1d_r_legp(int nr, double* tb, double *xo) {
   
  assert (nr>2) ;
  
  xo[0] = tb[0]/3. + 2.*tb[1]/15. ;
  
  for (int i=1 ; i<nr-1 ; i++)
    xo[i] = double(2*i*(2*i-1))*tb[i-1]/double((4*i-1)*(4*i-3))
      + (double(4*i*i)/double((4*i-1)*(4*i+1)) 
	 + double((2*i+1)*(2*i+1))/double((4*i+1)*(4*i+3)))*tb[i]
      + double((2*i+1)*(2*i+2))*tb[i+1]/double((4*i+3)*(4*i+5)) ;
  xo[nr-1] = double((2*nr-2)*(2*nr-3))*tb[nr-2]/double((4*nr-5)*(4*nr-7))
    + (double(4*(nr-1)*(nr-1))/double((4*nr-5)*(4*nr-3))
       + double((2*nr-1)*(2*nr-1))/double((4*nr-3)*(4*nr-1)))*tb[nr-1] ;
}


			//--------------
			// cas R_LEGI --
			//--------------

void _multx2_1d_r_legi(int nr, double* tb, double *xo){

  assert(nr>2) ;

  xo[0] = 0.6*tb[0] + 6.*tb[1]/35. ;

    for (int i=1 ; i<nr-1 ; i++)
      xo[i] = double(2*i*(2*i+1))*tb[i-1]/double((4*i+1)*(4*i-1))
	+ (double((2*i+1)*(2*i+1))/double((4*i+1)*(4*i+3))
	   + double((2*i+2)*(2*i+2))/double((4*i+3)*(4*i+5)) )*tb[i]
	+ double((2*i+2)*(2*i+3))*tb[i+1]/double((4*i+7)*(4*i+5)) ;
    xo[nr-1] = 0. ;
}

			//-------------
			// cas R_LEG --
			//-------------

void _multx2_1d_r_leg(int nr, double* tb, double *xo){

    assert(nr>3) ;

    xo[0] = tb[0]/3. + 2.*tb[2]/15. ;
    xo[1] = 0.6*tb[1] + 6.*tb[3]/35. ;
    for (int i=2 ; i<nr-2 ; i++)
      xo[i] = double(i*(i-1))*tb[i-2]/double((2*i-1)*(2*i-3))
	+ (double(i*i)/double((2*i-1)*(2*i+1)) 
	   + double((i+1)*(i+1))/double((2*i+1)*(2*i+3)))*tb[i]
	+ double((i+1)*(i+2))*tb[i+2]/double((2*i+3)*(2*i+5)) ;
    for (int i=nr-2 ; i<nr ; i++)
	xo[i] = double(i*(i-1))*tb[i-2]/double((2*i-1)*(2*i-3))
	+ (double(i*i)/double((2*i-1)*(2*i+1)) 
	   + double((i+1)*(i+1))/double((2*i+1)*(2*i+3)))*tb[i] ;
}




		    //----------------------
		    // La routine a appeler
		    //----------------------
		    
void multx2_1d(int nr, double **tb, int base_r)	    // Version appliquee a this
{

// Routines de derivation
static void (*multx2_1d[MAX_BASE])(int, double *, double*) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    multx2_1d[i] = _multx2_1d_pas_prevu ;
	}
	// Les routines existantes
	multx2_1d[R_CHEB >> TRA_R] = _multx2_1d_r_cheb ;
	multx2_1d[R_CHEBP >> TRA_R] = _multx2_1d_r_chebp ;
	multx2_1d[R_CHEBI >> TRA_R] = _multx2_1d_r_chebi ;
	multx2_1d[R_LEG >> TRA_R] = _multx2_1d_r_leg ;
	multx2_1d[R_LEGP >> TRA_R] = _multx2_1d_r_legp ;
	multx2_1d[R_LEGI >> TRA_R] = _multx2_1d_r_legi ;
    }
    
    double *result = new double[nr] ;
    multx2_1d[base_r](nr, *tb, result) ;
    
    delete [] (*tb) ;
    (*tb) = result ;
}		
}
