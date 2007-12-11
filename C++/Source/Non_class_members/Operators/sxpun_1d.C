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


char sxpun_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2007/12/11 15:28:21  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
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


#include <stdlib.h>
#include <math.h>

#include "tbl.h"
#include "type_parite.h"

/*
 * Operateur :
 *  -R_CHEB : (f-f(-1))/(x+1)
 *  -R_JACO02 : (f-f(-1))/(x+1)
 * 
 *
 * Entree : coefficients de f dans tb
 *	    nr : nombre de points en r
 * Sortie : coefficient du resultat dans tb
 * 
 * 
 */


		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

void _sxpun_1d_pas_prevu(int nr, double* tb, double *res) {
    cout << "sxpun pas prevu..." << endl ;
    cout << " valeurs: " << tb << "   " << res << endl ;
    cout << "nr : " << nr << endl ;
    abort () ;
    exit(-1) ;
}



			//---------------
			// cas R_CHEB ---
			//---------------

void _sxpun_1d_r_cheb (int nr, double* tb, double *xo)
{
    
    assert (nr >= 3) ; 
    
    xo[nr-1] = 0 ;
    xo[nr-2] = 2*(tb[nr-1]-xo[nr-1]) ;
    
    for (int i=nr-3 ; i>0 ; i--)
	xo[i] = 2*tb[i+1]-2*xo[i+1]-xo[i+2] ;
    
    double somme = 0 ;
    for (int i=0 ; i<nr ; i++)
	if (i%2 == 0)
	    somme += tb[i] ;
	else
	    somme -= tb[i] ;
    
    xo[0] = tb[0]-xo[1]/2.-somme ;
}


			//---------------
			// cas R_JACO02 -
			//---------------

void _sxpun_1d_r_jaco02 (int nr, double* tb, double *xo)
{
    
    xo[nr-1] = 0 ;
    double somme ;
    for (int i = 0 ; i < nr-1 ; i++) {
	somme = 0 ;
	for (int j = i+1 ; j < nr ; j++) {
	somme += pow((-1),(j-1-i))*((j+1)*(j+2)/double((i+1)*(i+2))-(i+1)*(i+2)/double((j+1)*(j+2)))*tb[j] ;
	}
	xo[i] = (2*i+3)/double(4)*somme ;
    }
}

		    //----------------------------
		    // La routine a appeler   ----
		    //----------------------------
		    
		    
void sxpun_1d(int nr, double **tb, int base_r)	    // Version appliquee a this
{

// Routines de derivation
static void (*sxpun_1d[MAX_BASE])(int, double *, double *) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    sxpun_1d[i] = _sxpun_1d_pas_prevu ;
	}
	// Les routines existantes
	sxpun_1d[R_CHEB >> TRA_R] = _sxpun_1d_r_cheb ;
	sxpun_1d[R_JACO02 >> TRA_R] = _sxpun_1d_r_jaco02 ;
	}
    
    double *result = new double[nr] ;
    sxpun_1d[base_r](nr, *tb, result) ;
    
    delete [] (*tb) ;
    (*tb) = result ;
}
