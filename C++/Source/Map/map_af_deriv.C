/*
 * Computations of Cmp partial derivatives for a Map_af mapping
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


char map_af_deriv_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.12  2000/02/25  08:59:51  eric
 * Remplacement de ci.get_dzpuis() == 0  par ci.check_dzpuis(0).
 * Suppression de l'affectation des dzpuis Mtbl/Mtnl_cf a la fin car
 *   c'est fait par Cmp::set_dzpuis.
 *
 * Revision 2.11  2000/01/26  13:09:18  eric
 * Reprototypage complet des routines de derivation:
 * le resultat est desormais suppose alloue a l'exterieur de la routine
 * et est passe en argument (Cmp& resu), si bien que le prototypage
 * complet devient:
 *             void DERIV(const Cmp& ci, Cmp& resu) const
 *
 * Revision 2.10  1999/11/30  12:51:32  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.9  1999/11/26  14:23:55  eric
 * Traitement dzpuis des Cmp.
 *
 * Revision 2.8  1999/11/26  10:58:02  eric
 * Traitement dzpuis.
 *
 * Revision 2.7  1999/11/25  16:29:29  eric
 * Reorganisation complete du calcul des derivees partielles.
 *
 * Revision 2.6  1999/10/27  15:45:23  eric
 * Suppression du membre Cmp::c.
 *
 * Revision 2.5  1999/10/27  08:47:03  eric
 * Introduction de Cmp::va a la place de *(Cmp::c).
 *
 * Revision 2.4  1999/10/22  08:16:21  eric
 * const Map*.
 *
 * Revision 2.3  1999/10/14  14:27:17  eric
 * Methodes const.
 *
 * Revision 2.2  1999/10/13  15:54:40  eric
 * Mg3d* -> const Mg3d*
 *
 * Revision 2.1  1999/09/17  10:01:09  phil
 * correction pour deriv_x et deriv_y
 *
 * Revision 2.0  1999/09/14  16:37:06  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */
 
// Header Lorene
#include "map.h"
#include "cmp.h"


			//---------------------//
			//	d/dr	       //
			//---------------------//
			
// In the ZEC : 

void Map_af::dsdr(const Cmp& ci, Cmp& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp()->get_mg() == mg) ; 

    
    if (ci.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
    }
    else {    
	assert( ci.get_etat() == ETATQCQ ) ; 
	assert( ci.check_dzpuis(0) ) ; 

	(ci.va).coef() ;    // (ci.va).c_cf is up to date
	
	resu = (ci.va).dsdx() * dxdr ;     //  dxi/dR, - dxi/dU (ZEC)
	
	(resu.va).base = (ci.va).dsdx().base ;	// same basis as d/dxi
	
	int nz = mg->get_nzone() ; 
	if (mg->get_type_r(nz-1) == UNSURR) {
	    resu.set_dzpuis(2) ;	    // r^2 d/dr has been computed in the
					    // external domain
	}

    }
    
}

			//------------------------//
			//	1/r d/dtheta      //
			//------------------------//

void Map_af::srdsdt(const Cmp& ci, Cmp& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp()->get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
    }
    else {

	assert( ci.get_etat() == ETATQCQ ) ; 
	assert( ci.check_dzpuis(0) ) ; 

	(ci.va).coef() ;    // (ci.va).c_cf is up to date

	Valeur tmp = ci.va ; 
	
	tmp = tmp.dsdt() ;	// d/dtheta
	tmp = tmp.sx() ;	// 1/xi, Id, 1/(xi-1)
	
	Base_val sauve_base( tmp.base ) ; 
	
	tmp = tmp * xsr ;	// xi/R, 1/R, (xi-1)/U
	
	tmp.base = sauve_base ;   // The above operation does not the basis

	resu = tmp ;

	int nz = mg->get_nzone() ; 
	if (mg->get_type_r(nz-1) == UNSURR) {
	    resu.set_dzpuis(2) ;	    // r d/dtheta has been computed in
					    // the external domain
	}

    }
    
}


			//------------------------------------//
			//	1/(r sin(theta))  d/dphi      //
			//------------------------------------//

void Map_af::srstdsdp(const Cmp& ci, Cmp& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp()->get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
    }
    else {

	assert( ci.get_etat() == ETATQCQ) ; 
	assert( ci.check_dzpuis(0) ) ; 

	(ci.va).coef() ;    // (ci.va).c_cf is up to date

	Valeur tmp = ci.va ; 
	
	tmp = tmp.dsdp() ;	// d/dphi
	tmp = tmp.ssint() ;	// 1/sin(theta)
	tmp = tmp.sx() ;	// 1/xi, Id, 1/(xi-1)
	
	Base_val sauve_base( tmp.base ) ; 
	
	tmp = tmp * xsr ;	// xi/R, 1/R, (xi-1)/U
	
	tmp.base = sauve_base ;   // The above operation does not change the basis

	resu = tmp ;

	int nz = mg->get_nzone() ; 
	if (mg->get_type_r(nz-1) == UNSURR) {
	    resu.set_dzpuis(2) ;	    // r/sin(theta) d/dphi has been 
					    // computed in the external domain
	}

    }
    
}

