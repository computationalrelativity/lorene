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


char tenseur_inv_pois_vect_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2005/08/29 09:41:45  p_grandclement
 * Minor modif
 *
 * Revision 1.2  2002/08/07 16:14:11  j_novak
 * class Tenseur can now also handle tensor densities, this should be transparent to older codes
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/10/19  09:49:47  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//Standard
#include <stdlib.h>

//Lorene
#include "tenseur.h"

// Inversion de Poisson vectoriel :
Tenseur Tenseur::inverse_poisson_vect (double lambda) const {
    
    assert (valence == 1) ;
    assert (etat != ETATNONDEF) ;
    if (etat == ETATZERO)
	return (*this) ;

    Tenseur inverse (*mp, 1, CON, *get_triad(), metric, poids) ;
    Tenseur grad (contract(this->gradient(), 0, 1)) ;
    grad.dec2_dzpuis() ;
    Tenseur grad_shift (grad.gradient()) ;
    grad_shift.inc2_dzpuis() ;
    inverse.set_etat_qcq() ;
    for (int i=0 ; i<3 ; i++)
	inverse.set(i) = (*this)(i).laplacien(4)+lambda*grad_shift(i) ;
    
    return inverse ;
}
