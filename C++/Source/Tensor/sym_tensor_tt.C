/*
 *  Methods of class Sym_tensor_tt
 *
 *   (see file sym_tensor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
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


char sym_tensor_tt_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/03/03 13:16:21  j_novak
 * New potential khi (p_khi) and the functions manipulating it.
 *
 * Revision 1.2  2004/01/04 20:52:45  e_gourgoulhon
 * Added assignement (operator=) to a Tensor_sym.
 *
 * Revision 1.1  2003/10/27 10:50:54  e_gourgoulhon
 * First version.
 *
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>

// Headers Lorene
#include "tensor.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
Sym_tensor_tt::Sym_tensor_tt(const Map& map, const Base_vect& triad_i,
		const Metric& met) 
	: Sym_tensor_trans(map, triad_i, met ) {
				
	set_der_0x0() ;

}

// Copy constructor
// ----------------
Sym_tensor_tt::Sym_tensor_tt (const Sym_tensor_tt& source)
	: Sym_tensor_trans(source) {
    
	set_der_0x0() ;

	if (source.p_khi != 0x0) p_khi = new Scalar( *(source.p_khi) ) ; 
	if (source.p_eta != 0x0) p_eta = new Scalar( *(source.p_eta) ) ; 
	if (source.p_mu != 0x0) p_mu = new Scalar( *(source.p_mu) ) ; 
	
}   


// Constructor from a file
// -----------------------
Sym_tensor_tt::Sym_tensor_tt(const Map& mapping, const Base_vect& triad_i, 
	const Metric& met, FILE* fd) 
	: Sym_tensor_trans(mapping, triad_i, met, fd) {

	set_der_0x0() ;
}

			//--------------//
			//  Destructor  //
			//--------------//

Sym_tensor_tt::~Sym_tensor_tt() {

  Sym_tensor_tt::del_deriv() ;	// in order not to follow the virtual aspect
  									// of del_deriv()

}



			//-------------------//
			// Memory managment  //
			//-------------------//

void Sym_tensor_tt::del_deriv() const {

	if (p_khi != 0x0) delete p_khi ; 
	if (p_eta != 0x0) delete p_eta ; 
	if (p_mu != 0x0) delete p_mu ; 
	
	set_der_0x0() ;
	
	Sym_tensor_trans::del_deriv() ;

}

void Sym_tensor_tt::set_der_0x0() const {

  p_khi = 0x0 ;
  p_eta = 0x0 ; 
  p_mu = 0x0 ; 

}


			//--------------//
			//  Assignment  //
			//--------------//

void Sym_tensor_tt::operator=(const Sym_tensor_tt& source) {
    
    // Assignment of quantities common to all derived classes of Sym_tensor_trans
	Sym_tensor_trans::operator=(source) ; 
	
	del_deriv() ; 	
}


void Sym_tensor_tt::operator=(const Sym_tensor_trans& source) {
    
    // Assignment of quantities common to all derived classes of Sym_tensor_trans
	Sym_tensor_trans::operator=(source) ; 

	del_deriv() ; 	
}



void Sym_tensor_tt::operator=(const Sym_tensor& source) {
    
    // Assignment of quantities common to all derived classes of Sym_tensor_trans
	Sym_tensor_trans::operator=(source) ; 

	del_deriv() ; 	
}


void Sym_tensor_tt::operator=(const Tensor_sym& source) {
    
    // Assignment of quantities common to all derived classes of Sym_tensor_trans
    Sym_tensor_trans::operator=(source) ; 
	
    del_deriv() ; 	
}


void Sym_tensor_tt::operator=(const Tensor& source) {
    
    // Assignment of quantities common to all derived classes of Sym_tensor_trans
	Sym_tensor_trans::operator=(source) ; 
	
	del_deriv() ; 	
}




