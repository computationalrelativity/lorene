/*
 *  Member functions of the class Cmp for various r manipulations
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 2001 Jerome Novak
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


char cmp_r_manip_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.5  2001/10/29  15:37:06  novak
 * Ajout de Cmp::div_r()
 *
 * Revision 1.4  2000/08/31 13:04:46  eric
 * Ajout des fonctions mult_rsint et div_rsint.
 *
 * Revision 1.3  2000/05/22  14:39:52  phil
 * ajout de inc_dzpuis et dec_dzpuis
 *
 * Revision 1.2  1999/12/10  16:33:48  eric
 * Appel de del_deriv().
 *
 * Revision 1.1  1999/11/30  14:22:54  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

#include "cmp.h" 


			//---------------------------//
			//	    div_r	     //
			//---------------------------//

void Cmp::div_r() {
    
    mp->div_r(*this) ;   // Call of the appropriate routine of the mapping
    
    del_deriv() ;   // Delete the derived members

}
			//---------------------------//
			//	    mult_r	     //
			//---------------------------//

void Cmp::mult_r() {
    
    mp->mult_r(*this) ;   // Call of the appropriate routine of the mapping
    
    del_deriv() ;   // Delete the derived members
    
}

			//---------------------------//
			//	    mult_r_zec	     //
			//---------------------------//

void Cmp::mult_r_zec() {
    
    mp->mult_r_zec(*this) ;   // Call of the appropriate routine of the mapping
    
    del_deriv() ;   // Delete the derived members

}

			//---------------------------//
			//	    mult_rsint	     //
			//---------------------------//

void Cmp::mult_rsint() {
    
    mp->mult_rsint(*this) ;   // Call of the appropriate routine of the mapping 
    
    del_deriv() ;   // Delete the derived members

}

			//---------------------------//
			//	    div_rsint	     //
			//---------------------------//

void Cmp::div_rsint() {
    
    mp->div_rsint(*this) ;   // Call of the appropriate routine of the mapping
    
    del_deriv() ;   // Delete the derived members

}

			//---------------------------//
			//	    dec_dzpuis	     //
			//---------------------------//

void Cmp::dec_dzpuis() {
    
    mp->dec_dzpuis(*this) ;   // Call of the appropriate routine of the mapping
    
}

			//---------------------------//
			//	    inc_dzpuis	     //
			//---------------------------//

void Cmp::inc_dzpuis() {
    
    mp->inc_dzpuis(*this) ;   // Call of the appropriate routine of the mapping
    
}



			//---------------------------//
			//	    dec2_dzpuis	     //
			//---------------------------//

void Cmp::dec2_dzpuis() {
    
    mp->dec2_dzpuis(*this) ;   // Call of the appropriate routine of the mapping 
    
}

			//---------------------------//
			//	    inc2_dzpuis	     //
			//---------------------------//

void Cmp::inc2_dzpuis() {
    
    mp->inc2_dzpuis(*this) ;   // Call of the appropriate routine of the mapping 
    
}


