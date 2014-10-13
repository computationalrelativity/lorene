/*
 *  Member functions of the class Scalar for various theta manipulations
 *
 *    See file scalar.h for documentation. 
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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


char scalar_th_manip_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:47  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2006/05/26 09:00:12  j_novak
 * New members for multiplication or division by cos(theta).
 *
 * Revision 1.1  2003/11/04 23:00:59  e_gourgoulhon
 * First version.
 *
 *
 *
 * $Header$
 *
 */

#include "tensor.h" 


			//-------------------//
			//	    mult_cost    //
			//-------------------//


namespace Lorene {
void Scalar::mult_cost() {

    mp->mult_cost(*this) ;   // Call of the appropriate routine of the mapping 
    
    del_deriv() ;   // Delete the derived members

}
			//-------------------//
			//	    div_cost     //
			//-------------------//


void Scalar::div_cost() {

    mp->div_cost(*this) ;   // Call of the appropriate routine of the mapping 
    
    del_deriv() ;   // Delete the derived members

}


			//-------------------//
			//	    mult_sint    //
			//-------------------//


void Scalar::mult_sint() {

    mp->mult_sint(*this) ;   // Call of the appropriate routine of the mapping 
    
    del_deriv() ;   // Delete the derived members

}


			//-------------------//
			//	    div_sint     //
			//-------------------//


void Scalar::div_sint() {

    mp->div_sint(*this) ;   // Call of the appropriate routine of the mapping 
    
    del_deriv() ;   // Delete the derived members

}

			//-------------------//
			//	    div_tant     //
			//-------------------//


void Scalar::div_tant() {

    mp->div_tant(*this) ;   // Call of the appropriate routine of the mapping 
    
    del_deriv() ;   // Delete the derived members

}




}
