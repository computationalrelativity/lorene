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


 

/*
 * $Id$
 * $Log$
 * Revision 1.6  2016/12/05 16:17:57  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:53:03  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2004/03/01 09:57:03  j_novak
 * the wave equation is solved with Scalars. It now accepts a grid with a
 * compactified external domain, which the solver ignores and where it copies
 * the values of the field from one time-step to the next.
 *
 * Revision 1.3  2003/06/18 08:45:27  j_novak
 * In class Mg3d: added the member get_radial, returning only a radial grid
 * For dAlembert solver: the way the coefficients of the operator are defined has been changed.
 *
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
#include "tensor.h"
#include "param.h"

namespace Lorene {
Mtbl_cf sol_dalembert(Param&, const Map_af&, const Mtbl_cf&) ;

//*****************************************************************************

void Map_et::dalembert(Param& , Scalar& fJp1, const Scalar& fJ, const Scalar& fJm1,
		       const Scalar& source) const {
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp().get_mg() == mg) ; 
    assert(fJ.get_etat() != ETATNONDEF) ; 
    assert(fJ.get_mp().get_mg() == mg) ; 
    assert(fJm1.get_etat() != ETATNONDEF) ; 
    assert(fJm1.get_mp().get_mg() == mg) ; 
    assert(fJp1.get_mp().get_mg() == mg) ; 

    cout << "Map_et_dalembert:" << endl ;
    cout << "Not implemented" << endl ;
    cout << fJp1 << fJ << fJm1 << source ;
    abort() ;

    
}


}
