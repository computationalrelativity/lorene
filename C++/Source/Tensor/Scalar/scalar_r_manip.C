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
 * Revision 1.1  2003/09/25 09:12:01  e_gourgoulhon
 * First version (uses Cmp as intermediate variable).
 *
 *
 * $Header$
 *
 */

#include "tensor.h" 
#include "cmp.h" 


			//---------------------------//
			//	    div_r	     //
			//---------------------------//

void Scalar::div_r() {
    
	Cmp cuu(*this) ; 
	
    mp->div_r(cuu) ;   // Call of the appropriate routine of the mapping
    
	operator=(cuu) ; 
	
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


