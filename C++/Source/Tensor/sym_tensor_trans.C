/*
 *  Methods of class Sym_tensor_trans
 *
 *   (see file sym_tensor.h for documentation)
 *
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


char sym_tensor_trans_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/10/28 21:24:52  e_gourgoulhon
 * Added new methods trace() and tt_part().
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
#include "metric.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
Sym_tensor_trans::Sym_tensor_trans(const Map& map, const Base_vect& triad_i,
		const Metric& met) 
	: Sym_tensor(map, CON, triad_i),
	  met_div(&met) {
				
	set_der_0x0() ;

}

// Copy constructor
// ----------------
Sym_tensor_trans::Sym_tensor_trans (const Sym_tensor_trans& source)
	: Sym_tensor(source), 
	  met_div(source.met_div) {
    
	set_der_0x0() ;

	if (source.p_trace != 0x0) p_trace = new Scalar( *(source.p_trace) ) ; 
	if (source.p_tt != 0x0) p_tt = new Sym_tensor_tt( *(source.p_tt) ) ; 
	
}   


// Constructor from a file
// -----------------------
Sym_tensor_trans::Sym_tensor_trans(const Map& mapping, const Base_vect& triad_i, 
	const Metric& met, FILE* fd) 
	: Sym_tensor(mapping, triad_i, fd),
	  met_div(&met) {

	set_der_0x0() ;
}

			//--------------//
			//  Destructor  //
			//--------------//

Sym_tensor_trans::~Sym_tensor_trans() {

  Sym_tensor_trans::del_deriv() ;	// in order not to follow the virtual aspect
  									// of del_deriv()

}



			//-------------------//
			// Memory managment  //
			//-------------------//

void Sym_tensor_trans::del_deriv() const {

	if (p_trace != 0x0) delete p_trace ; 
	if (p_tt != 0x0) delete p_tt ; 
	
	set_der_0x0() ;
	Sym_tensor::del_deriv() ;

}

void Sym_tensor_trans::set_der_0x0() const {

	p_trace = 0x0 ; 
	p_tt = 0x0 ; 

}


			//--------------//
			//  Assignment  //
			//--------------//

void Sym_tensor_trans::operator=(const Sym_tensor_trans& source) {
    
    // Assignment of quantities common to all the derived classes of Sym_tensor
	Sym_tensor::operator=(source) ; 

    // Assignment of proper quantities of class Sym_tensor_trans
	assert(met_div == source.met_div) ; 
	
	del_deriv() ; 	
}


void Sym_tensor_trans::operator=(const Sym_tensor& source) {
    
    // Assignment of quantities common to all the derived classes of Sym_tensor
	Sym_tensor::operator=(source) ; 

  	// The metric which was set by the constructor is kept
	
	del_deriv() ; 	
}


void Sym_tensor_trans::operator=(const Tensor& source) {
    
    // Assignment of quantities common to all the derived classes of Sym_tensor
	Sym_tensor::operator=(source) ; 

  	// The metric which was set by the constructor is kept
	
	del_deriv() ; 	
}


			//-----------------------------//
			//    Computational methods    //
			//-----------------------------//

const Scalar& Sym_tensor_trans::trace() const {

	if (p_trace == 0x0) {   // a new computation is necessary

		assert( (type_indice(0)==CON) && (type_indice(1)==CON) ) ; 
		
		Tensor tmp = contract( met_div->cov(), 0, *this, 0 ) ;
		 
		p_trace = new Scalar( tmp.scontract(0,1) ) ; 
		
	}
	
	return *p_trace ; 

}


const Sym_tensor_tt& Sym_tensor_trans::tt_part() const {

	if (p_tt == 0x0) {   // a new computation is necessary

		cout << "Sym_tensor_trans::tt_part() : not implemented yet ! " << endl ; 
		abort() ; 
		
	}
	
	return *p_tt ; 

}






