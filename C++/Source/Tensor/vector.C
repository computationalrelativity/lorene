/*
 *  Methods of class Vector
 *
 *   (see file vector.h for documentation)
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


char vector_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2003/10/06 13:58:48  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.4  2003/10/05 21:14:20  e_gourgoulhon
 * Added method std_spectral_base().
 *
 * Revision 1.3  2003/10/03 14:10:32  e_gourgoulhon
 * Added constructor from Tensor.
 *
 * Revision 1.2  2003/10/03 14:08:46  j_novak
 * Removed old change_trid...
 *
 * Revision 1.1  2003/09/26 08:05:31  j_novak
 * New class Vector.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <assert.h>

// Headers Lorene
#include "tensor.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
Vector::Vector(const Map& map, int tipe, const Base_vect& triad_i) 
		: Tensor(map, 1, tipe, triad_i) {
		
  set_der_0x0() ;

}

// Standard constructor with the triad passed as a pointer
// -------------------------------------------------------
Vector::Vector(const Map& map, int tipe, const Base_vect* triad_i) 
		: Tensor(map, 1, tipe, *triad_i) {
		
  set_der_0x0() ;
}
	
// Copy constructor
// ----------------
Vector::Vector (const Vector& source) : 
    Tensor(source) {
  
  assert(valence == 1) ;
  set_der_0x0() ;

}   


// Constructor from a {\tt Tensor}.
//--------------------------------
Vector::Vector(const Tensor& uu) : Tensor(uu) {

  assert(valence == 1) ;
  set_der_0x0() ;

}


// Constructor from a file
// -----------------------
Vector::Vector(const Map& mapping, const Base_vect& triad_i, FILE* fd) : 
  Tensor(mapping, triad_i, fd) {
   
  assert ( (valence == 1) && (n_comp == 3) ) ;
  set_der_0x0() ;

}


			//--------------//
			//  Destructor  //
			//--------------//


Vector::~Vector () {
  del_deriv() ;
}


			//-------------------//
			// Memory managment  //
			//-------------------//

void Vector::del_deriv() const {

  set_der_0x0() ;
  Tensor::del_deriv() ;

}

void Vector::set_der_0x0() const {

}

void Vector::operator=(const Tensor& t) {
    
    assert (t.valence == 1) ;

    del_deriv() ;

    triad = t.triad ; 

    assert(t.type_indice(0) == type_indice(0)) ;

    for (int i=0 ; i<3 ; i++) {
      *cmp[i] = *t.cmp[i] ;
    }
    del_deriv() ;
}



// Affectation d'une composante :
Scalar& Vector::set(int index) {
    
  assert ( (index>=1) && (index<=3) ) ;

  del_deriv() ;
  
  return *cmp[index - 1] ;
}

const Scalar& Vector::operator()(int index) const {
    
  assert ( (index>=1) && (index<=3) ) ;

  return *cmp[index - 1] ;

}


// Sets the standard spectal bases of decomposition for each component

void Vector::std_spectral_base() {

	Base_val** bases = 0x0 ;

	if ( triad->identify() == (mp->get_bvect_cart()).identify() ) {

		// Cartesian case
		bases = mp->get_mg()->std_base_vect_cart() ;

	}
	else {
		// Spherical case
		assert( triad->identify() == (mp->get_bvect_spher()).identify()) ;
		bases = mp->get_mg()->std_base_vect_spher() ;
	}
	    
	for (int i=0 ; i<3 ; i++) {
		cmp[i]->set_spectral_base( *bases[i] ) ; 
	}
		
	for (int i=0 ; i<3 ; i++) {
		delete bases[i] ;
	}
	delete [] bases ;


}

