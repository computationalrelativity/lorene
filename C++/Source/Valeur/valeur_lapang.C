/*
 * Function Valeur::lapang for the computation of the angular Laplacian:
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


char valeur_lapang_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/13 08:53:50  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  1999/11/30  12:44:22  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.3  1999/11/23  16:17:37  eric
 * Reorganisation du calcul dans le cas ETATZERO.
 *
 * Revision 2.2  1999/11/19  09:32:40  eric
 * La valeur de retour est desormais const Valeur &.
 *
 * Revision 2.1  1999/10/18  13:42:09  eric
 * Suppression de l'argument base dans les routines de derivation des mtbl_cf.
 *
 * Revision 2.0  1999/04/26  16:43:22  phil
 * *** empty log message ***
 *
 * Revision 1.1  1999/04/26  16:42:30  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */


// Headers Lorene
#include "mtbl_cf.h"
#include "valeur.h"


namespace Lorene {
const Valeur& Valeur::lapang() const {
    
    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_lapang != 0x0) {
	return *p_lapang ;
    }
    
    // ... si, il faut bosser

    p_lapang = new Valeur(mg) ;

    if (etat == ETATZERO) {
	p_lapang->set_etat_zero() ; 
    }
    else {
	assert(etat == ETATQCQ) ; 
	p_lapang->set_etat_cf_qcq() ;
	Mtbl_cf* cfp = p_lapang->c_cf ; // Pointeur sur le Mtbl_cf qui vient d'etre
					// cree par le set_etat_cf_qcq()

	// Initialisation de *cfp : recopie des coef. de la fonction
	if (c_cf == 0x0) {
	    coef() ;
	}
	*cfp = *c_cf ;	
 
	cfp->lapang() ;	// calcul 
    
	p_lapang->base = cfp->base ; // On remonte la base de sortie au niveau Valeur
    }

    // Termine
    return *p_lapang ;    
    
}
}
