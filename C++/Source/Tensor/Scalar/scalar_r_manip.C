/*
 *  Member functions of the class Scalar for various r manipulations
 *
 *    See file scalar.h for documentation. 
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Eric Gourgoulhon  (for a preceding Cmp version)
 *   Copyright (c) 1999-2001 Philippe Grandclement  (for a preceding Cmp version)
 *   Copyright (c) 2001 Jerome Novak (for a preceding Cmp version)
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


char scalar_r_manip_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/10/08 12:26:03  j_novak
 * Second part of the bug (sorry!)
 *
 * Revision 1.3  2003/10/08 12:19:12  j_novak
 * Bug corrected, thanks to purify
 *
 * Revision 1.2  2003/10/05 21:16:41  e_gourgoulhon
 * Added methods div_r_ced() and div_rsint_ced().
 *
 * Revision 1.1  2003/09/25 09:12:01  e_gourgoulhon
 * First version (uses Cmp as intermediate variable).
 *
 *
 * $Header$
 *
 */

#include "tensor.h" 
#include "cmp.h" 


			//-------------------//
			//	    div_r        //
			//-------------------//

void Scalar::div_r() {
    
	Cmp cuu(*this) ; 
	
    mp->div_r(cuu) ;   // Call of the appropriate routine of the mapping
    
	operator=(cuu) ; 
	
    del_deriv() ;   // Delete the derived members

}


			//---------------------//
			//	    div_r_ced      //
			//---------------------//


void Scalar::div_r_ced() {
    
	if (etat == ETATZERO) {
		dzpuis += 2 ; 
		return ; 
	}
	
	assert(etat == ETATQCQ) ; 
	
	int nzm1 = mp->get_mg()->get_nzone() - 1 ; // index of the CED
	
	// Copy of the CED part of *this into uu_ext and multiplication by r
	Scalar uu_ext(*mp) ; 
	uu_ext.allocate_all() ;
	uu_ext.annule(0,nzm1-1) ; // zero in all domains but the CED
	uu_ext.set_domain(nzm1) = domain(nzm1) ; 
	uu_ext.set_spectral_base(va.base) ; 
	uu_ext.mult_r_zec() ; // multiplication by r in the CED

	// Division by r in all domains but the CED
	annule(nzm1, nzm1) ; 	// zero in the CED
	div_r() ; 
	
	// Add the CED part
	set_domain(nzm1) = uu_ext.domain(nzm1) ; 
	
	dzpuis += 2 ; 
	
    del_deriv() ;   // Delete the derived members

}

			//---------------------------//
			//	    mult_r	     //
			//---------------------------//

void Scalar::mult_r() {
    
	Cmp cuu(*this) ; 

    mp->mult_r(cuu) ;   // Call of the appropriate routine of the mapping
    
	operator=(cuu) ; 
	
    del_deriv() ;   // Delete the derived members
    
}

			//---------------------------//
			//	    mult_r_zec	     //
			//---------------------------//

void Scalar::mult_r_zec() {
    
	Cmp cuu(*this) ; 

    mp->mult_r_zec(cuu) ;   // Call of the appropriate routine of the mapping
    
	operator=(cuu) ; 

    del_deriv() ;   // Delete the derived members

}

			//---------------------------//
			//	    mult_rsint	     //
			//---------------------------//

void Scalar::mult_rsint() {
    
	Cmp cuu(*this) ; 

    mp->mult_rsint(cuu) ;   // Call of the appropriate routine of the mapping 
    
	operator=(cuu) ; 

    del_deriv() ;   // Delete the derived members

}

			//---------------------------//
			//	    div_rsint	     //
			//---------------------------//

void Scalar::div_rsint() {
    
	Cmp cuu(*this) ; 

    mp->div_rsint(cuu) ;   // Call of the appropriate routine of the mapping
    
	operator=(cuu) ; 

    del_deriv() ;   // Delete the derived members

}


			//-------------------------//
			//	    div_rsint_ced      //
			//-------------------------//


void Scalar::div_rsint_ced() {
    
	if (etat == ETATZERO) {
		dzpuis += 2 ; 
		return ; 
	}
	
	assert(etat == ETATQCQ) ; 
	
	int nzm1 = mp->get_mg()->get_nzone() - 1 ; // index of the CED
	
	// Copy of the CED part of *this into uu_ext and multiplication by r
	Scalar uu_ext(*mp) ; 
	uu_ext.allocate_all() ;
	uu_ext.annule(0,nzm1-1) ; // zero in all domains but the CED
	uu_ext.set_domain(nzm1) = domain(nzm1) ; 
	uu_ext.set_spectral_base(va.base) ; 
	uu_ext.mult_r_zec() ; // multiplication by r in the CED

	// Division by sin(theta) in the CED :

	// what follows does not apply if the mapping is not radial:
	assert( dynamic_cast<const Map_radial*>(mp) != 0x0 ) ; 
	Valeur vtmp = (uu_ext.get_spectral_va()).ssint() ;
	uu_ext.set_spectral_va() = vtmp ;  

	// Division by r sin(theta) in all domains but the CED
	annule(nzm1, nzm1) ; 	// zero in the CED
	div_rsint() ; 
	
	// Add the CED part
	set_domain(nzm1) = uu_ext.domain(nzm1) ; 
	
	dzpuis += 2 ; 
	
    del_deriv() ;   // Delete the derived members

}




			//---------------------------//
			//	    dec_dzpuis	     //
			//---------------------------//

void Scalar::dec_dzpuis() {
    
	Cmp cuu(*this) ; 

    mp->dec_dzpuis(cuu) ;   // Call of the appropriate routine of the mapping
    
	operator=(cuu) ; 

}

			//---------------------------//
			//	    inc_dzpuis	     //
			//---------------------------//

void Scalar::inc_dzpuis() {
    
	Cmp cuu(*this) ; 

    mp->inc_dzpuis(cuu) ;   // Call of the appropriate routine of the mapping
    
	operator=(cuu) ; 
}



			//---------------------------//
			//	    dec2_dzpuis	     //
			//---------------------------//

void Scalar::dec2_dzpuis() {
    
	Cmp cuu(*this) ; 

    mp->dec2_dzpuis(cuu) ;   // Call of the appropriate routine of the mapping 
    
	operator=(cuu) ; 

}

			//---------------------------//
			//	    inc2_dzpuis	     //
			//---------------------------//

void Scalar::inc2_dzpuis() {
    
	Cmp cuu(*this) ; 

    mp->inc2_dzpuis(cuu) ;   // Call of the appropriate routine of the mapping 
    
	operator=(cuu) ; 
}


