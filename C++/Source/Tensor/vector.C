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
 * Revision 1.8  2003/10/20 09:32:12  j_novak
 * Members p_potential and p_div_free of the Helmholtz decomposition
 * + the method decompose_div(Metric).
 *
 * Revision 1.7  2003/10/16 14:21:37  j_novak
 * The calculation of the divergence of a Tensor is now possible.
 *
 * Revision 1.6  2003/10/13 13:52:40  j_novak
 * Better managment of derived quantities.
 *
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
#include "metric.h"


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

  Vector::del_deriv() ;

}


			//-------------------//
			// Memory managment  //
			//-------------------//

void Vector::del_deriv() const {

  for (int i=0; i<N_MET_MAX; i++) 
    del_derive_met(i) ;
  
  set_der_0x0() ;
  Tensor::del_deriv() ;

}

void Vector::set_der_0x0() const {

  for (int i=0; i<N_MET_MAX; i++) 
    set_der_met_0x0(i) ;

}

void Vector::del_derive_met(int j) const {

  assert( (j>=0) && (j<N_MET_MAX) ) ;

  if (met_depend[j] != 0x0) {
    if (p_potential[j] != 0x0)
      delete p_potential[j] ;
    if (p_div_free[j] != 0x0)
      delete p_div_free[j] ;
    
    set_der_met_0x0(j) ;
    
    Tensor::del_derive_met(j) ;
  }
}

void Vector::set_der_met_0x0(int i) const {

  assert( (i>=0) && (i<N_MET_MAX) ) ;

  p_potential[i] = 0x0 ;
  p_div_free[i] = 0x0 ;

  Tensor::set_der_met_0x0(i) ;

}
void Vector::operator=(const Vector& t) {
    
    triad = t.triad ; 

    assert(t.type_indice(0) == type_indice(0)) ;

    for (int i=0 ; i<3 ; i++) {
      *cmp[i] = *t.cmp[i] ;
    }
    del_deriv() ;
}

void Vector::operator=(const Tensor& t) {
    
    assert (t.valence == 1) ;

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

const Scalar& Vector::divergence(const Metric& metre) const {
  
  set_dependance(metre) ;
  int j = get_place_met(metre) ;
  assert ((j>=0) && (j<N_MET_MAX)) ;
  if (p_divergence[j] == 0x0) {
    p_divergence[j] = metre.get_connect().p_divergence(*this) ;
  }

  const Scalar* pscal = dynamic_cast<const Scalar*>(p_divergence[j]) ;

  assert(pscal != 0x0) ;

  return *pscal ;
}

const Scalar& Vector::potential(const Metric& metre) const {

  set_dependance(metre) ;
  int j = get_place_met(metre) ;
  assert ((j>=0) && (j<N_MET_MAX)) ;
  if (p_potential[j] == 0x0) {
    decompose_div(metre) ;
  }

  return *p_potential[j] ;
}

const Vector_divfree& Vector::div_free(const Metric& metre) const {

  set_dependance(metre) ;
  int j = get_place_met(metre) ;
  assert ((j>=0) && (j<N_MET_MAX)) ;
  if (p_div_free[j] == 0x0) {
    decompose_div(metre) ;
  }

  return *p_div_free[j] ;
}

void Vector::decompose_div(const Metric& metre) const {

  assert( type_indice(0) == CON ) ; //Only for contravariant vectors...

  set_dependance(metre) ;
  int j =  get_place_met(metre) ;
  assert ((j>=0) && (j<N_MET_MAX)) ;

  if ( (p_potential[j] != 0x0) && (p_div_free[j] != 0x0) ) 
    return ; // Nothing to do ...

  else {
    if (p_div_free[j] != 0x0)
      delete p_div_free[j] ;  
    if (p_potential[j] != 0x0) 
      delete p_potential[j] ;

    p_potential[j] = new Scalar( (divergence(metre)).poisson() ) ;

    p_div_free[j] = new Vector_divfree(*mp, *triad, metre) ;

    *p_div_free[j] = ( *this - p_potential[j]->derive_con(metre) ) ;

  }
  
}

