/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 1999-2002 Eric Gourgoulhon
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


char som_phi_C[] = "$Header$" ;

/*
 * Ensemble des routine pour la sommation directe en phi
 * 
 *   SYNOPSYS:
 *     double som_phi_XX
 *	(double* ti, int np, double phi, double* xo)
 *
 *   ATTENTION: np est le nombre reel de points.
 *		on suppose que ti contient les n+2 
 *		avec les 0 qu'il faut.
 * 
 */

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/10/16 14:36:58  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/05/05 16:21:28  e_gourgoulhon
 * Error message (for unknown basis) in English.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2000/09/08  16:07:02  eric
 * Ajout de la base P_COSSIN_I
 *
 * Revision 2.1  2000/03/06  09:34:58  eric
 * Suppression des #include inutiles.
 *
 * Revision 2.0  1999/04/12  15:43:21  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <math.h>

#include "headcpp.h"

void som_phi_pas_prevu
    (double*, const int, const double, double*) {
	cout << "Mtbl_cf::val_point: phi basis not implemented yet ! "
	     << endl ;
	abort () ;
}

void som_phi_cossin
    (double* ti, const int np, const double phi, double* xo) {
    
    *xo = ti[0] ;   // premier element

    // Sommation sur les cosinus et les sinus
    for (int k=2 ; k<np-1 ; k +=2 ) {
	int m = k/2 ;
	*xo += ti[k] * cos(m * phi) ;
	*xo += ti[k+1] * sin(m * phi) ;
    }
    *xo += ti[np] * cos(np/2 * phi) ;
}

void som_phi_cossin_p
    (double* ti, const int np, const double phi, double* xo) {
    
    *xo = ti[0] ;   // premier element

    // Sommation sur les cosinus et les sinus
    for (int k=2 ; k<np-1 ; k +=2 ) {
	int m = 2*(k/2) ;
	*xo += ti[k] * cos(m * phi) ;
	*xo += ti[k+1] * sin(m * phi) ;
    }
    *xo += ti[np] * cos(np * phi) ;
}

void som_phi_cossin_i
    (double* ti, const int np, const double phi, double* xo) {
    
    *xo = ti[0] * cos(phi) + ti[2] * sin(phi) ; 

    // Sommation sur les harmoniques d'ordre m >= 3 : 
    for (int k=3 ; k<np ; k +=2 ) {
	int m = k ;
	*xo += ti[k] * cos(m * phi) ;
	*xo += ti[k+1] * sin(m * phi) ;
    }

}
