/*
 * Output computational coordinates \f$\xi\f$ 
 *
 * for:
 *   - Valeur
 *   - Mtbl_cf
 */

/*
 *   Copyright (c) 1999-2001 Jerome Novak
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


 


// Headers C
#include <cassert>

// Headers Lorene
#include "mtbl.h"
#include "valeur.h"

namespace Lorene {
void Valeur::va_x() {

    // Protection
    assert(etat != ETATNONDEF) ; 

    for(int l=0; l<mg->get_nzone(); l++)
     for(int k=0;k<mg->get_np(l);k++)
      for(int j=0;j<mg->get_nt(l);j++)
       for(int i=0;i<mg->get_nr(l);i++)
	    set(l,k,j,i) = mg->get_grille3d(l)->x[i];

}
}
