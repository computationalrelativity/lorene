/*
 * Function Mtbl_cf::lapang for the computation of the angular Laplacian:
 *
 *  d^2/dtheta^2 + cos(theta)/sin(theta) d/dtheta + 1/sin(theta) d^2/dphi^2
 *
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char mtbl_cf_lapang_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/09/16 12:11:59  j_novak
 * Added the base T_LEG_II.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2000/10/04  14:55:45  eric
 * Ajout des bases T_LEG_IP et T_LEG_PI.
 *
 * Revision 2.2  1999/10/18  13:41:58  eric
 * Suppression de l'argument base dans les routines de derivation des mtbl_cf.
 *
 * Revision 2.1  1999/09/30  12:54:31  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/04/26  16:42:17  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Headers Lorene
#include "mtbl_cf.h"
#include "base_val.h"
#include "type_parite.h"


// Prototypage des fonctions utilisees:
void _lapang_pas_prevu(Mtbl_cf *, int) ;
void _lapang_t_leg_p(Mtbl_cf *, int) ;
void _lapang_t_leg_i(Mtbl_cf *, int) ;
void _lapang_t_leg_pp(Mtbl_cf *, int) ;
void _lapang_t_leg_ip(Mtbl_cf *, int) ;
void _lapang_t_leg_pi(Mtbl_cf *, int) ;
void _lapang_t_leg_ii(Mtbl_cf *, int) ;

//*****************************************************************************

void Mtbl_cf::lapang()		    // Version appliquee a this
{

// Routines de derivation
static void (*_lapang[MAX_BASE])(Mtbl_cf *, int) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    _lapang[i] = _lapang_pas_prevu ;
	}
	// Les routines existantes
	_lapang[T_LEG_P >> TRA_T] = _lapang_t_leg_p ;
	_lapang[T_LEG_PP >> TRA_T] = _lapang_t_leg_pp ;
	_lapang[T_LEG_I >> TRA_T] = _lapang_t_leg_i ;
	_lapang[T_LEG_IP >> TRA_T] = _lapang_t_leg_ip ;
	_lapang[T_LEG_PI >> TRA_T] = _lapang_t_leg_pi ;
	_lapang[T_LEG_II >> TRA_T] = _lapang_t_leg_ii ;
    }

    // Boucle sur les zones
    for (int l=0 ; l<get_mg()->get_nzone() ; l++) {
	int base_t = (base.b[l] & MSQ_T) >> TRA_T ;
	_lapang[base_t](this, l) ;
    }
    
}
