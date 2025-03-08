/*
 *  Method of the class Map_et for computing the integral of a Cmp over
 *  all space.
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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
 * Revision 1.6  2016/12/05 16:17:57  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:53:05  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:13  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2004/07/26 16:02:23  j_novak
 * Added a flag to specify whether the primitive should be zero either at r=0
 * or at r going to infinity.
 *
 * Revision 1.2  2004/06/14 15:27:35  e_gourgoulhon
 * Added method primr (not ready yet !).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2000/08/16  12:12:12  eric
 * Ajout de ciaff.set_dzpuis( ci.get_dzpuis() ) pour tenir compte de
 *  la suppression de Mtbl::dzpuis.
 *
 * Revision 2.1  2000/01/17  12:40:32  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/01/17  11:17:06  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Headers C
#include <cstdlib>


// Headers Lorene
#include "map.h"
#include "cmp.h"

namespace Lorene {
Tbl* Map_et::integrale(const Cmp& ci) const {

    assert(ci.get_etat() != ETATNONDEF) ; 

    int nz = mg->get_nzone() ; 
    
    if (ci.get_etat() == ETATZERO) {
	Tbl* resu = new Tbl(nz) ;
	resu->annule_hard() ; 
	return resu ; 
    }
    
    assert( ci.get_etat() == ETATQCQ ) ; 
            
    // Construction of an affine mapping to call Map_af::integrale   
    Map_af mpaff(*this) ; 

    Cmp ciaff(mpaff) ; 
    
    // Multiplication by the reducted Jacobian of the mapping
    
    Base_val sauve_base = (ci.va).base ; 
    
    ciaff = (ci.va) * rsx2drdx ;
    
    (ciaff.va).set_base(sauve_base) ; 
	
    ciaff.set_dzpuis( ci.get_dzpuis() ) ; 
        
    // Call to the Map_af version :
    
    return mpaff.integrale(ciaff) ;
    
}


void Map_et::primr(const Scalar& , Scalar& , bool) const {

    cout << "Map_et::primr : not ready yet !" << endl ; 
    abort() ; 
}
}
