/*
 *   Copyright (c) 2000-2001 Jerome Novak
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


char map_et_dalembert_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/01/03 15:30:28  j_novak
 * Some comments modified.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2001/10/16  10:07:52  novak
 * *** empty log message ***
 *
 * Revision 1.2  2001/07/19 14:13:55  novak
 * new list of arguments for Map_et::dalembert
 *
 * Revision 1.1  2000/10/19 15:41:15  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "map.h"
#include "cmp.h"
#include "param.h"

Mtbl_cf sol_dalembert(Param&, const Map_af&, const Mtbl_cf&) ;

//*****************************************************************************

void Map_et::dalembert(Param& par, Cmp& fJp1, const Cmp& fJ, const Cmp& fJm1,
		       const Cmp& source) const {
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp()->get_mg() == mg) ; 
    assert(fJ.get_etat() != ETATNONDEF) ; 
    assert(fJ.get_mp()->get_mg() == mg) ; 
    assert(fJm1.get_etat() != ETATNONDEF) ; 
    assert(fJm1.get_mp()->get_mg() == mg) ; 
    assert(fJp1.get_mp()->get_mg() == mg) ; 

    assert(par.get_n_double() >= 1) ;
    cout << "Not implemented" << endl ;
    cout << par.get_n_double() << fJp1 << fJ << fJm1 << source ;
    abort() ;

    
}


