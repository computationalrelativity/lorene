/*
 *  Methods of class Connection_fcart.
 *
 *	(see file connection.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003	Eric Gourgoulhon & Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char connection_fcart_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2003/10/16 15:26:03  e_gourgoulhon
 * Suppressed unsued variable
 *
 * Revision 1.5  2003/10/16 14:21:36  j_novak
 * The calculation of the divergence of a Tensor is now possible.
 *
 * Revision 1.4  2003/10/11 16:45:43  e_gourgoulhon
 * Suppressed the call to Itbl::set_etat_qcq() after
 * the construction of the Itbl's.
 *
 * Revision 1.3  2003/10/11 14:39:50  e_gourgoulhon
 * Suppressed declaration of unusued arguments in some methods.
 *
 * Revision 1.2  2003/10/06 13:58:46  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.1  2003/10/03 14:11:48  e_gourgoulhon
 * Methods of class Connection_fcart.
 *
 * 
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>

// Lorene headers
#include "connection.h"


//------------------------------//
//	       Constructors     //
//------------------------------//



// Contructor from a Cartesian flat-metric-orthonormal basis

Connection_fcart::Connection_fcart(const Map& mpi, const Base_vect_cart& bi) 
  : Connection_flat(mpi, bi) {

}		

// Copy constructor
Connection_fcart::Connection_fcart(const Connection_fcart& ci) 
  : Connection_flat(ci) {

}		

	
//----------------------------//
//	       Destructor     //
//----------------------------//


Connection_fcart::~Connection_fcart(){
	
}


//-----------------------------//
//     Mutators / assignment   //
//-----------------------------//


void Connection_fcart::operator=(const Connection_fcart& ) {
	
  cout << "Connection_fcart::operator= : not implemented yet !" << endl ; 
  abort() ; 

}	



//-----------------------------//
//    Computational methods    //
//-----------------------------//

// Covariant derivative, returning a value.
//-----------------------------------------

Tensor Connection_fcart::derive_cov(const Tensor& uu) const {

  int valence0 = uu.get_valence() ; 
  int ncomp0 = uu.get_n_comp() ;
	
  // Protections
  // -----------
  if (valence0 >= 1) {
    assert(uu.get_triad() == triad) ; 
  }

  // Indices of the result
  // ---------------------
  Itbl tipe(valence0+1) ; 
  tipe.set(0) = COV ; 
  const Itbl tipeuu = uu.get_index_type() ;  
  for (int id = 1; id<=valence0; id++) {
    tipe.set(id) = tipeuu(id-1) ; 
  }

  // Creation of the result tensor
  // -----------------------------
  Tensor resu(*mp, valence0+1, tipe, *triad) ;
	
  Itbl ind1(valence0+1) ; // working Itbl to store the indices of resu
	
  Itbl ind0(valence0) ; // working Itbl to store the indices of uu
	

  // Derivation index = x
  // --------------------
  int k = 1 ; 	

  // Loop on all the components of the input tensor
  for (int ic=0; ic<ncomp0; ic++) {
	
    // indices corresponding to the component no. ic in the input tensor
    ind0 = uu.indices(ic) ; 
		
    // indices (k,ind0) in the output tensor
    ind1.set(0) = k ; 
    for (int id = 1; id<=valence0; id++) {
      ind1.set(id) = ind0(id-1) ; 
    }
		
    Scalar& cresu = resu.set(ind1) ; 
		
    cresu = (uu(ind0)).dsdx() ; 	// d/dx
		
  }



  // Derivation index = y
  // ---------------------
  k = 2 ; 	

  // Loop on all the components of the input tensor
  for (int ic=0; ic<ncomp0; ic++) {
	
    // indices corresponding to the component no. ic in the input tensor
    ind0 = uu.indices(ic) ; 
		
    // indices (k,ind0) in the output tensor
    ind1.set(0) = k ; 
    for (int id = 1; id<=valence0; id++) {
      ind1.set(id) = ind0(id-1) ; 
    }
		
    Scalar& cresu = resu.set(ind1) ; 
		
    cresu = (uu(ind0)).dsdy() ;  	// d/dy	
		
  }


  // Derivation index = z
  // --------------------
  k = 3 ; 	

  // Loop on all the components of the input tensor
  for (int ic=0; ic<ncomp0; ic++) {
	
    // indices corresponding to the component no. ic in the input tensor
    ind0 = uu.indices(ic) ; 
		
    // indices (k,ind0) in the output tensor
    ind1.set(0) = k ; 
    for (int id = 1; id<=valence0; id++) {
      ind1.set(id) = ind0(id-1) ; 
    }
		
    Scalar& cresu = resu.set(ind1) ; 
		
    cresu = (uu(ind0)).dsdz() ;  	// d/dz

  }


  // C'est fini !
  // -----------
  return resu ; 

}


// Covariant derivative, returning a pointer.
//-------------------------------------------

Tensor* Connection_fcart::p_derive_cov(const Tensor& uu) const {

  int valence0 = uu.get_valence() ; 
  int ncomp0 = uu.get_n_comp() ;
	
  // Protections
  // -----------
  if (valence0 >= 1) {
    assert(uu.get_triad() == triad) ; 
  }

  // Indices of the result
  // ---------------------
  Itbl tipe(valence0+1) ; 

  // Creation of the result pointer
  // ------------------------------
  Tensor* resu ;

  // If u is a Scalar, the result is a vector
  //----------------------------------------
  if (valence0 == 0) 
    resu = new Vector(*mp, COV, triad) ;
  else {
    tipe.set(0) = COV ; 
    const Itbl tipeuu = uu.get_index_type() ;  
    for (int id = 1; id<=valence0; id++) {
      tipe.set(id) = tipeuu(id-1) ; 
    }
    const Sym_tensor* stuu 
      = dynamic_cast<const Sym_tensor*>(&uu) ;
    if (stuu != 0x0) { //Then the type Tensor_delta reduces the storage
      resu = new Tensor_delta(*mp, tipe, *triad) ;
    }
    else { //Most general case...
      resu = new Tensor(*mp, valence0+1, tipe, *triad) ;
    }
  }
	
  Itbl ind1(valence0+1) ; // working Itbl to store the indices of resu
	
  Itbl ind0(valence0) ; // working Itbl to store the indices of uu
	

  // Derivation index = x
  // --------------------
  int k = 1 ; 	

  // Loop on all the components of the input tensor
  for (int ic=0; ic<ncomp0; ic++) {
	
    // indices corresponding to the component no. ic in the input tensor
    ind0 = uu.indices(ic) ; 
		
    // indices (k,ind0) in the output tensor
    ind1.set(0) = k ; 
    for (int id = 1; id<=valence0; id++) {
      ind1.set(id) = ind0(id-1) ; 
    }
		
    Scalar& cresu = resu->set(ind1) ; 
		
    cresu = (uu(ind0)).dsdx() ; 	// d/dx
		
  }



  // Derivation index = y
  // ---------------------
  k = 2 ; 	

  // Loop on all the components of the input tensor
  for (int ic=0; ic<ncomp0; ic++) {
	
    // indices corresponding to the component no. ic in the input tensor
    ind0 = uu.indices(ic) ; 
		
    // indices (k,ind0) in the output tensor
    ind1.set(0) = k ; 
    for (int id = 1; id<=valence0; id++) {
      ind1.set(id) = ind0(id-1) ; 
    }
		
    Scalar& cresu = resu->set(ind1) ; 
		
    cresu = (uu(ind0)).dsdy() ;  	// d/dy	
		
  }


  // Derivation index = z
  // --------------------
  k = 3 ; 	

  // Loop on all the components of the input tensor
  for (int ic=0; ic<ncomp0; ic++) {
	
    // indices corresponding to the component no. ic in the input tensor
    ind0 = uu.indices(ic) ; 
		
    // indices (k,ind0) in the output tensor
    ind1.set(0) = k ; 
    for (int id = 1; id<=valence0; id++) {
      ind1.set(id) = ind0(id-1) ; 
    }
		
    Scalar& cresu = resu->set(ind1) ; 
		
    cresu = (uu(ind0)).dsdz() ;  	// d/dz

  }


  // C'est fini !
  // -----------
  return resu ; 

}

// Divergence, returning a pointer.
//---------------------------------

Tensor* Connection_fcart::p_divergence(const Tensor& uu) const {

  int valence0 = uu.get_valence() ; 
	
  // Protections
  // -----------
  assert (valence0 >= 1) ;
  assert (uu.get_triad() == triad) ; 
  assert (uu.get_index_type(0) == CON) ;


  // Indices of the result
  // ---------------------
  Itbl tipe(valence0-1) ; 

  // Creation of the pointer on the result tensor
  // --------------------------------------------
  Tensor* resu ;

  // If u is a Vector, the result is a Scalar
  //----------------------------------------
  if (valence0 == 1) 
    resu = new Scalar(*mp) ;
  else {
    const Itbl tipeuu = uu.get_index_type() ;  
    for (int id = 0; id<valence0-1; id++) {
      tipe.set(id) = tipeuu(id+1) ; 
    }
    if (valence0 == 2) {
      resu = new Vector(*mp, tipe(0), *triad) ;
    }
    else {
      const Tensor_delta* del_uu 
	= dynamic_cast<const Tensor_delta*>(&uu) ;
      if (del_uu != 0x0) { //Then the type Sym_tensor reduces the storage
	resu = new Sym_tensor(*mp, tipe, *triad) ;
      }
      else { //Most general case...
	resu = new Tensor(*mp, valence0-1, tipe, *triad) ;
      }
    }
  }

  int ncomp1 = resu->get_n_comp() ;
	
  Itbl ind0(valence0) ; // working Itbl to store the indices of uu
	
  Itbl ind1(valence0-1) ; // working Itbl to store the indices of resu
	
  // Loop on all the components of the output tensor
  for (int ic=0; ic<ncomp1; ic++) {
	
    ind1 = resu->indices(ic) ; 
    Scalar& cresu = resu->set(ind1) ;
    cresu.set_etat_zero() ;

    for (int k=1; k<=3; k++) {
      
      // indices (k,ind1) in the input tensor
      ind0.set(0) = k ; 
      for (int id = 1; id<valence0; id++) {
	ind0.set(id) = ind1(id-1) ; 
      }
      cresu += uu(ind0).deriv(k) ; //Addition of dT^i/dx^i
    }

  }

  // C'est fini !
  // -----------
  return resu ; 

}

