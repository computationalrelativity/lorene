/*
 *  Methods of class Vector_divfree
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


char vector_divfree_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/10/22 13:08:06  j_novak
 * Better handling of dzpuis flags
 *
 * Revision 1.3  2003/10/17 16:34:32  e_gourgoulhon
 * Added new methods set_vr_eta_mu and set_vr_mu.
 *
 * Revision 1.2  2003/10/15 13:52:57  j_novak
 * Initialization of met_div in the copy constructor
 *
 * Revision 1.1  2003/10/13 20:53:53  e_gourgoulhon
 * First version
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
Vector_divfree::Vector_divfree(const Map& map, const Base_vect& triad_i,
		const Metric& met) 
		: Vector(map, CON, triad_i),
			met_div(&met) {
		
  set_der_0x0() ;

}

	
// Copy constructor
// ----------------
Vector_divfree::Vector_divfree (const Vector_divfree& source) : 
    Vector(source), met_div(source.met_div) {
  
	set_der_0x0() ;
	if (source.p_eta != 0x0) p_eta = new Scalar(*(source.p_eta)) ;
	if (source.p_mu != 0x0) p_mu = new Scalar(*(source.p_mu)) ;

}   


// Constructor from a file
// -----------------------
Vector_divfree::Vector_divfree(const Map& mapping, const Base_vect& triad_i, 
	const Metric& met, FILE* fd) : Vector(mapping, triad_i, fd),
	met_div(&met) {
   
  set_der_0x0() ;

}


			//--------------//
			//  Destructor  //
			//--------------//


Vector_divfree::~Vector_divfree () {

  Vector_divfree::del_deriv() ;

}


			//-------------------//
			// Memory managment  //
			//-------------------//

void Vector_divfree::del_deriv() const {

  if (p_eta != 0x0) delete p_eta ; 
  if (p_mu != 0x0) delete p_mu ; 
  set_der_0x0() ;
  Vector::del_deriv() ;

}

void Vector_divfree::set_der_0x0() const {

	p_eta = 0x0 ; 
	p_mu = 0x0 ; 

}

		//-------------------------//
		//  Mutators / assignment  //
		//-------------------------//

void Vector_divfree::operator=(const Vector_divfree& source) {
    
    // Assignment of quantities common to all the derived classes of Vector
	Vector::operator=(source) ; 

    // Assignment of proper quantities of class Vector_divfree
	assert(met_div == source.met_div) ; 
	
	del_deriv() ; 	
	
}

void Vector_divfree::operator=(const Vector& source) {
    
  // Assignment of quantities common to all the derived classes of Vector
  Vector::operator=(source) ; 

  //The metric which was set by the constructor is kept

  del_deriv() ;
}

void Vector_divfree::operator=(const Tensor& source) {
    
  // Assignment of quantities common to all the derived classes of Vector
  Vector::operator=(source) ; 
  
  //The metric which was set by the constructor is kept

  del_deriv() ;
}


void Vector_divfree::set_vr_eta_mu(const Scalar& vr_i, const Scalar eta_i,
		const Scalar mu_i) {
		
		// All this has a meaning only for spherical components:
		#ifndef NDEBUG 
		const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
		assert(bvs != 0x0) ; 
		#endif
		
		del_deriv() ; // delete previous p_eta and p_mu, as well as 
					  //  derived quantities
		
		*cmp[0] = vr_i ; 	// V^r
		
		p_eta = new Scalar( eta_i ) ; 	// eta

		p_mu = new Scalar( mu_i ) ; 	// mu 
		
		update_vtvp() ; // V^theta and V^phi
		

}


void Vector_divfree::set_vr_mu(const Scalar& vr_i, const Scalar mu_i) {
		
		// All this has a meaning only for spherical components:
		#ifndef NDEBUG 
		const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
		assert(bvs != 0x0) ; 
		#endif
		
		del_deriv() ; // delete previous p_eta and p_mu, as well as 
					  //  derived quantities
		
		*cmp[0] = vr_i ; 	// V^r
		
		p_mu = new Scalar( mu_i ) ; 	// mu 
		
		eta() ; // computes eta form the divergence-free condition

		update_vtvp() ; // V^theta and V^phi		

}



