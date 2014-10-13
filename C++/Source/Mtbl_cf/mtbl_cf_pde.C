/*
 * Method Mtbl_cf::poisson_angu().
 *
 *  (see file mtbl_cf.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003-2005 Eric Gourgoulhon & Jerome Novak
 *   Copyright (c) 2005 Michael Forot
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


char mtbl_cf_pde_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2014/10/13 08:53:08  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2009/10/23 12:56:19  j_novak
 * New base T_LEG_MI
 *
 * Revision 1.4  2009/10/13 19:44:41  j_novak
 * New base T_LEG_MP.
 *
 * Revision 1.3  2005/04/04 21:32:13  e_gourgoulhon
 * Added argument lambda to method poisson_angu
 * to deal with the generalized angular Poisson equation:
 *     Lap_ang u + lambda u = source.
 *
 * Revision 1.2  2004/12/17 13:35:03  m_forot
 * Add the case T_LEG
 *
 * Revision 1.1  2003/10/15 21:12:22  e_gourgoulhon
 * First version.
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
namespace Lorene {
void _poisangu_pas_prevu(Mtbl_cf *, int, double) ;
void _poisangu_t_leg_p(Mtbl_cf *, int, double) ;
void _poisangu_t_leg_i(Mtbl_cf *, int, double) ;
void _poisangu_t_leg_pp(Mtbl_cf *, int, double) ;
void _poisangu_t_leg_ip(Mtbl_cf *, int, double) ;
void _poisangu_t_leg_pi(Mtbl_cf *, int, double) ;
void _poisangu_t_leg_ii(Mtbl_cf *, int, double) ;
void _poisangu_t_leg_mp(Mtbl_cf *, int, double) ;
void _poisangu_t_leg_mi(Mtbl_cf *, int, double) ;
void _poisangu_t_leg(Mtbl_cf *, int, double) ;

//*****************************************************************************

void Mtbl_cf::poisson_angu(double lambda) {

	// Routines de derivation
	static void (*poisangu[MAX_BASE])(Mtbl_cf *, int, double) ;
	static int nap = 0 ;

    // Premier appel
    if (nap==0) {
		nap = 1 ;
		for (int i=0 ; i<MAX_BASE ; i++) {
	    	poisangu[i] = _poisangu_pas_prevu ;
		}
		// Les routines existantes
		poisangu[T_LEG_P >> TRA_T] = _poisangu_t_leg_p ;
		poisangu[T_LEG_PP >> TRA_T] = _poisangu_t_leg_pp ;
		poisangu[T_LEG_I >> TRA_T] = _poisangu_t_leg_i ;
		poisangu[T_LEG_IP >> TRA_T] = _poisangu_t_leg_ip ;
		poisangu[T_LEG_PI >> TRA_T] = _poisangu_t_leg_pi ;
		poisangu[T_LEG_MP >> TRA_T] = _poisangu_t_leg_mp ;
		poisangu[T_LEG_MI >> TRA_T] = _poisangu_t_leg_mi ;
		poisangu[T_LEG >> TRA_T] = _poisangu_t_leg ;
    }

    // Boucle sur les zones
    for (int l=0 ; l<get_mg()->get_nzone() ; l++) {
		int base_t = (base.b[l] & MSQ_T) >> TRA_T ;
		poisangu[base_t](this, l, lambda) ;
    }
    
}
}
