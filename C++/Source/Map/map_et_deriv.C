/*
 * Computations of Cmp partial derivatives for a Map_et mapping
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


char map_et_deriv_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.3  2000/02/25  09:01:28  eric
 * Remplacement de ci.get_dzpuis() == 0  par ci.check_dzpuis(0).
 * Suppression de l'affectation des dzpuis Mtbl/Mtnl_cf a la fin car
 *   c'est fait par Cmp::set_dzpuis.
 *
 * Revision 1.2  2000/01/26  13:09:52  eric
 * Reprototypage complet des routines de derivation:
 * le resultat est desormais suppose alloue a l'exterieur de la routine
 * et est passe en argument (Cmp& resu), si bien que le prototypage
 * complet devient:
 *            void DERIV(const Cmp& ci, Cmp& resu) const
 *
 * Revision 1.1  1999/12/17  12:59:29  eric
 * Initial revision
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
			
void Map_et::dsdr(const Cmp& ci, Cmp& resu) const {

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

void Map_et::srdsdt(const Cmp& ci, Cmp& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp()->get_mg() == mg) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
    }
    else {

	assert( ci.get_etat() == ETATQCQ ) ; 
	assert( ci.check_dzpuis(0) ) ; 

	(ci.va).coef() ;    // (ci.va).c_cf is up to date

	// Computation of 1/R df/dtheta'   ---> srdfdt
	// ----------------------------
	Valeur srdfdt = ci.va ; 
	
	srdfdt = srdfdt.dsdt() ;	// d/dtheta'
	srdfdt = srdfdt.sx() ;		// 1/xi, Id, 1/(xi-1)
	
	Base_val sauve_base( srdfdt.base ) ; 
	
	srdfdt = srdfdt * xsr ;	// xi/R, 1/R, (xi-1)/U
	
	srdfdt.base = sauve_base ;   // The above operation does not change the basis

	// Computation of 1/(dR/dx) 1/R dR/dtheta' df/dx   ----> adfdx
	// ----------------------------------------------

	Valeur adfdx = ci.va ; 

	adfdx = adfdx.dsdx()  ;	    // df/dx 
		    
	sauve_base = adfdx.base ; 
	adfdx = adfdx * dxdr * srdrdt ;  // 1/(dR/dx) 1/R dR/dtheta' df/dx
	adfdx.base = sauve_base ; 

	// Final result 
	// ------------

	resu = srdfdt - adfdx ;

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

void Map_et::srstdsdp(const Cmp& ci, Cmp& resu) const {

    assert (ci.get_etat() != ETATNONDEF) ; 
    assert (ci.get_mp()->get_mg() == mg) ; 

    if (ci.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
    }
    else {

	assert( ci.get_etat() == ETATQCQ) ; 
	assert( ci.check_dzpuis(0) ) ; 

	(ci.va).coef() ;    // (ci.va).c_cf is up to date

	// Computation of 1/(R sin(theta')) df/dphi'   ---> srstdfdp
	// -----------------------------------------

	Valeur srstdfdp = ci.va ; 
	
	srstdfdp = srstdfdp.dsdp() ;	// d/dphi
	srstdfdp = srstdfdp.ssint() ;	// 1/sin(theta)
	srstdfdp = srstdfdp.sx() ;	// 1/xi, Id, 1/(xi-1)
	
	Base_val sauve_base( srstdfdp.base ) ; 
	
	srstdfdp = srstdfdp * xsr ;	// xi/R, 1/R, (xi-1)/U
	
	srstdfdp.base = sauve_base ;   // The above operation does not change the basis

	// Computation of 1/(dR/dx) 1/(R sin(theta') dR/dphi' df/dx   --> bdfdx
	// --------------------------------------------------------
	Valeur bdfdx = ci.va ; 

	bdfdx = bdfdx.dsdx()  ;	    // df/dx 
		    
	sauve_base = bdfdx.base ; 
	bdfdx = bdfdx * dxdr * srstdrdp  ;  
	bdfdx.base = sauve_base ; 

	// Final result 
	// ------------

	resu = srstdfdp - bdfdx ;

	int nz = mg->get_nzone() ; 
	if (mg->get_type_r(nz-1) == UNSURR) {
	    resu.set_dzpuis(2) ;	    // r/sin(theta) d/dphi has been 
					    // computed in the external domain
	}

    }
    
}

