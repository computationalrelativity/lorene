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


char map_poisson_vect_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/07/09 16:46:23  p_grandclement
 * The Param in the case of an affine mapping is now 0x0 and not deleted
 * (I wonder why it was working before)
 *
 * Revision 1.2  2002/05/07 07:10:45  e_gourgoulhon
 * Compatibilty with xlC compiler on IBM SP2:
 * 	suppressed the parenthesis around argument of instruction new:
 * 	e.g.   aa = new (Tbl*[nzone])  --->  aa = new Tbl*[nzone]
 * 		result = new (Param)   --->  result = new Param
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2000/03/07  16:53:51  eric
 * niter est desormais passe en parametres.
 *
 * Revision 2.2  2000/02/15  10:25:49  phil
 * suppression des fonctions poisson_vect et poisson_vect_oohara
 *
 * Revision 2.1  2000/02/09  10:01:46  phil
 * ajout version oohara
 *
 * Revision 2.0  2000/01/21  12:59:10  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


#include "map.h"
#include "tenseur.h"
#include "param.h"

Param* Map_af::donne_para_poisson_vect(Param& para, int) const {
    return 0x0 ;
}


Param* Map_et::donne_para_poisson_vect(Param& para, int i) const {
    assert ((i>=0) && (i<4)) ;
    
    Param* result ;
    result = new Param ;
    result->add_int(para.get_int()) ;	//nbre max iterations
    result->add_double(para.get_double(0), 0) ;	// relaxation
    result->add_double(para.get_double(1), 1) ; // precision
    
    if (i!=3)
	result->add_cmp_mod(para.get_tenseur_mod().set(i)) ; //source au pas precedent.
    else
	result->add_cmp_mod(para.get_cmp_mod()) ; //la scalaire...
    

    result->add_int_mod(para.get_int_mod()) ;
    return result ;
}
