/*
 *  Computes 1/sin(theta) d/dphi   of a Valeur
 */

/*
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


char valeur_stdsdp_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:51  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:24  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  1999/11/23  16:18:43  eric
 * Reorganisation du calcul dans le cas ETATZERO.
 *
 * Revision 2.0  1999/11/19  11:22:30  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cassert>

// Headers Lorene
#include "valeur.h"

namespace Lorene {
const Valeur& Valeur::stdsdp() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_stdsdp != 0x0) {
	return *p_stdsdp ;
    }
    
    // ... si, il faut bosser

    p_stdsdp = new Valeur( dsdp() ) ;
 
    *p_stdsdp = p_stdsdp->ssint() ;
    
    // Termine
    return *p_stdsdp ;
}
}
