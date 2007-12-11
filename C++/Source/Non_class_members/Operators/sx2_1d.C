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


char sx2_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2007/12/11 15:28:18  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.2  2002/10/16 14:36:59  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  1999/07/08  09:55:21  phil
 * correction gestion memoire
 *
 * Revision 2.0  1999/07/07  10:15:03  phil
 * *** empty log message ***
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
 * Operateurs non singuliers :
 *	R_CHEB : Identite
 *	R_CHEBP, R_CHEBI : (f-f(0)-f'(0)x)/x^2
 *	R_CHEBU : (f-f(1)-(x-1)*f'(1))/(x-1)^2
 * 
 *Entree : tb contient coefficients de developpement en r
 *         nr : nombre de points en r.	
 *Sortie : tb contient les coefficients du resultat
 * 
 */

		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

void _sx2_1d_pas_prevu(int nr, double* tb, double *res) {
    cout << "sx2 pas prevu..." << tb << "    " << res << endl ;
    cout << "nr : " << nr << endl ;
    abort() ;
    exit(-1) ;
}


			//-------------
			// Identite --
			//------------

void _sx2_1d_identite(int nr, double*tb ,double *res) {
    for (int i=0 ; i<nr ; i++)
	res[i] = tb[i] ;
}


			//---------------
			// cas R_CHEBP --
			//---------------

void _sx2_1d_r_chebp(int nr, double* tb, double *xo)
{
   
   
   
    double somp, somn ;
    int sgn = 1 ;
	    
    xo[nr-1] = 0 ;
    somp = 4 * sgn * (nr-1) * tb[nr-1] ;
    somn = 2 * sgn * tb[nr-1] ;
    xo[nr-2] = somp - 2*(nr-2)*somn ;
    for (int i = nr-3 ; i >= 0 ; i-- ) {
	sgn = - sgn ;
	somp += 4 * (i+1) * sgn * tb[i+1] ;
	somn += 2 * sgn * tb[i+1] ;
	xo[i] = somp - 2*i * somn ;
    }	// Fin de la premiere boucle sur r
    for (int i=0 ; i<nr ; i+=2) {
	xo[i] = -xo[i] ;
    }	// Fin de la deuxieme boucle sur r
    xo[0] *= .5 ;

}


			//---------------
			// cas R_CHEBI --
			//---------------

void _sx2_1d_r_chebi(int nr, double* tb, double *xo)
{
    double somp, somn ;
    int sgn = 1 ;
	    
    xo[nr-1] = 0 ;
    somp = 2 * sgn * (2*(nr-1)+1) * tb[nr-1] ;
    somn = 2 * sgn * tb[nr-1] ;
    xo[nr-2] = somp - (2*(nr-2)+1)*somn ;
    for (int i = nr-3 ; i >= 0 ; i-- ) {
	sgn = - sgn ;
	somp += 2 * (2*(i+1)+1) * sgn * tb[i+1] ;
	somn += 2 * sgn * tb[i+1] ;
	xo[i] = somp - (2*i+1) * somn ;
    }	// Fin de la premiere boucle sur r
    for (int i=0 ; i<nr ; i+=2) {
	xo[i] = -xo[i] ;
    }	// Fin de la deuxieme boucle sur r

}
			//----------------
			// cas R_CHEBU --
			//---------------

void _sxm12_1d_r_chebu(int nr, double* tb, double *xo) {

    // On appelle juste deux fois la routine existante :
    sxm1_1d_cheb(nr, tb) ;
    sxm1_1d_cheb(nr, tb) ;
    for (int i=0 ; i<nr ; i++)
	xo[i] = tb[i] ;
}  

		    //----------------------
		    // La routine a appeler
		    //----------------------
		    
void sx2_1d(int nr, double **tb, int base_r)	    // Version appliquee a this
{

// Routines de derivation
static void (*sx2_1d[MAX_BASE])(int, double *, double*) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    sx2_1d[i] = _sx2_1d_pas_prevu ;
	}
	// Les routines existantes
	sx2_1d[R_CHEB >> TRA_R] = _sx2_1d_identite ;
	sx2_1d[R_CHEBP >> TRA_R] = _sx2_1d_r_chebp ;
	sx2_1d[R_CHEBI >> TRA_R] = _sx2_1d_r_chebi ;
	sx2_1d[R_CHEBU >> TRA_R] = _sxm12_1d_r_chebu ;
    }
    
    double *result = new double[nr] ;
    sx2_1d[base_r](nr, *tb, result) ;
    
    delete [] (*tb) ;
    (*tb) = result ;
}		
