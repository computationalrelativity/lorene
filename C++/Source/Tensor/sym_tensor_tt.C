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


 

/*
 * $Id$
 * $Log$
 * Revision 1.8  2016/12/05 16:18:17  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2014/10/13 08:53:43  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:13:19  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2005/04/01 14:28:32  j_novak
 * Members p_eta and p_mu are now defined in class Sym_tensor.
 *
 * Revision 1.4  2004/03/30 14:01:19  j_novak
 * Copy constructors and operator= now copy the "derived" members.
 *
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
#include <cstdlib>

// Headers Lorene
#include "tensor.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
namespace Lorene {
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
	
	set_der_0x0() ;
	
	Sym_tensor_trans::del_deriv() ;

}

void Sym_tensor_tt::set_der_0x0() const {

  p_khi = 0x0 ;
}


			//--------------//
			//  Assignment  //
			//--------------//

void Sym_tensor_tt::operator=(const Sym_tensor_tt& source) {
    
    // Assignment of quantities common to all derived classes of Sym_tensor_trans
	Sym_tensor_trans::operator=(source) ; 
	
	del_deriv() ; 
	
	if (source.p_khi != 0x0) p_khi = new Scalar( *(source.p_khi) ) ; 
	
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




}
