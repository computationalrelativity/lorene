/*
 *  Methods of class Connection_fspher.
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

char connection_fspher_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.11  2003/11/03 10:58:30  j_novak
 * Treatment of the general case for divergence.
 *
 * Revision 1.10  2003/10/22 13:08:03  j_novak
 * Better handling of dzpuis flags
 *
 * Revision 1.9  2003/10/16 15:26:48  e_gourgoulhon
 * Name of method Scalar::div_r_ced() changed to Scalar::div_r_inc2().
 *
 * Revision 1.8  2003/10/16 14:21:36  j_novak
 * The calculation of the divergence of a Tensor is now possible.
 *
 * Revision 1.7  2003/10/15 10:46:18  e_gourgoulhon
 * Introduced call to the new method Scalar::div_tant to perform
 * division by tan(theta) in derive_cov.
 *
 * Revision 1.6  2003/10/11 16:45:43  e_gourgoulhon
 * Suppressed the call to Itbl::set_etat_qcq() after
 * the construction of the Itbl's.
 *
 * Revision 1.5  2003/10/11 14:39:50  e_gourgoulhon
 * Suppressed declaration of unusued arguments in some methods.
 *
 * Revision 1.4  2003/10/06 13:58:47  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.3  2003/10/05 21:09:23  e_gourgoulhon
 * Method derive_cov: multiplication by r^2 in the CED.
 *
 * Revision 1.2  2003/10/01 21:49:45  e_gourgoulhon
 * First version of derive_cov --- not tested yet.
 *
 * Revision 1.1  2003/10/01 15:42:49  e_gourgoulhon
 * still ongoing...
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
//	       Constructors         //
//------------------------------//



// Contructor from a spherical flat-metric-orthonormal basis

Connection_fspher::Connection_fspher(const Map& mpi, const Base_vect_spher& bi) 
  : Connection_flat(mpi, bi) {

}		

// Copy constructor
Connection_fspher::Connection_fspher(const Connection_fspher& ci) 
  : Connection_flat(ci) {

}		

	
//----------------------------//
//	       Destructor         //
//----------------------------//


Connection_fspher::~Connection_fspher(){
	
}


//-----------------------------//
//     Mutators / assignment   //
//-----------------------------//


void Connection_fspher::operator=(const Connection_fspher& ) {
	
  cout << "Connection_fspher::operator= : not implemented yet !" << endl ; 
  abort() ; 

}	



//-----------------------------//
//    Computational methods    //
//-----------------------------//


// Covariant derivative, returning a value.
//-----------------------------------------

Tensor Connection_fspher::derive_cov(const Tensor& uu) const {

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
	
  Itbl ind(valence0) ; // working Itbl to store the indices of uu
	
  Scalar tmp(*mp) ;	// working scalar

	
  // Derivation index = r
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
		
    cresu = (uu(ind0)).dsdr() ; 	// d/dr
		
    // all the connection coefficients Gamma^i_{jk} are zero for k=1
  }



  // Derivation index = theta
  // ------------------------
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
    bool dz4 = (cresu.get_dzpuis() == 4) ;
		
    cresu = (uu(ind0)).srdsdt() ;  // 1/r d/dtheta 	
		
    // Loop on all the indices of uu
    for (int id=0; id<valence0; id++) {
		
      switch ( ind0(id) ) {
				
      case 1 : {	// Gamma^r_{l theta} V^l 
	// or -Gamma^l_{r theta} V_l 
	ind = ind0 ; 
	ind.set(id) = 2 ;   // l = theta

	// Division by r in all domains but the CED (where a
	//  multiplication by r is performed instead)
	tmp = uu(ind) ; 
	dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	cresu -= tmp ; 
	break ; 
      }
				
      case 2 : {	// Gamma^theta_{l theta} V^l 
	// or -Gamma^l_{theta theta} V_l
	ind = ind0 ; 
	ind.set(id) = 1 ;   // l = r
	tmp = uu(ind) ; 
	dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	cresu += tmp ; 
	break ; 
      }
				
      case 3 : {	// Gamma^phi_{l theta} V^l 
	// or -Gamma^l_{phi theta} V_l
	break ; 
      }
				
      default : {
	cout << "Connection_fspher::derive_cov : index problem ! "
	     << endl ; 
	abort() ;  
      }
      }

    }


  }


  // Derivation index = phi
  // ----------------------
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
    bool dz4 = (cresu.get_dzpuis() == 4) ;
		
    cresu = (uu(ind0)).srstdsdp() ;  // 1/(r sin(theta)) d/dphi 	
		
    // Loop on all the indices of uu
    for (int id=0; id<valence0; id++) {
		
      switch ( ind0(id) ) {
				
      case 1 : {	// Gamma^r_{l phi} V^l 
	// or -Gamma^l_{r phi} V_l 
	ind = ind0 ; 
	ind.set(id) = 3 ;   // l = phi
	tmp = uu(ind) ; 
	dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	cresu -= tmp ; 
	break ; 
      }
				
      case 2 : {	// Gamma^theta_{l phi} V^l 
	// or -Gamma^l_{theta phi} V_l
	ind = ind0 ; 
	ind.set(id) = 3 ;   // l = phi
	tmp = uu(ind) ; 
	dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	tmp.div_tant() ; 	// division by tan(theta)
										
	cresu -= tmp ; 
	break ; 
      }
				
      case 3 : {	// Gamma^phi_{l phi} V^l 
	// or -Gamma^l_{phi phi} V_l
							
	ind = ind0 ; 

	ind.set(id) = 1 ;   // l = r
	tmp = uu(ind) ; 
	dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	cresu += tmp ; 

	ind.set(id) = 2 ;   // l = theta
	tmp = uu(ind) ; 
	dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	tmp.div_tant() ; 	// division by tan(theta)

	cresu += tmp ; 
	break ; 
      }
				
      default : {
	cout << "Connection_fspher::derive_cov : index problem ! "
	     << endl ; 
	abort() ;  
      }
      }

    }


  }



  // C'est fini !
  // -----------
  return resu ; 

}



// Covariant derivative, returning a pointer.
//-------------------------------------------

Tensor* Connection_fspher::p_derive_cov(const Tensor& uu) const {

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

  // Creation of the pointer on the result tensor
  // --------------------------------------------
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
	
  Itbl ind(valence0) ; // working Itbl to store the indices of uu
	
  Scalar tmp(*mp) ;	// working scalar

	
  // Derivation index = r
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
		
    cresu = (uu(ind0)).dsdr() ; 	// d/dr
		
    // all the connection coefficients Gamma^i_{jk} are zero for k=1
  }



  // Derivation index = theta
  // ------------------------
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
		
    cresu = (uu(ind0)).srdsdt() ;  // 1/r d/dtheta 
    bool dz4 = (cresu.get_dzpuis() == 4) ;
		
    // Loop on all the indices of uu
    for (int id=0; id<valence0; id++) {
		
      switch ( ind0(id) ) {
				
      case 1 : {	// Gamma^r_{l theta} V^l 
	// or -Gamma^l_{r theta} V_l 
	ind = ind0 ; 
	ind.set(id) = 2 ;   // l = theta

	// Division by r in all domains but the CED (where a
	//  multiplication by r is performed instead)
	tmp = uu(ind) ; 
	dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	cresu -= tmp ; 
	break ; 
      }
				
      case 2 : {	// Gamma^theta_{l theta} V^l 
	// or -Gamma^l_{theta theta} V_l
	ind = ind0 ; 
	ind.set(id) = 1 ;   // l = r
	tmp = uu(ind) ; 
	dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	cresu += tmp ; 
	break ; 
      }
				
      case 3 : {	// Gamma^phi_{l theta} V^l 
	// or -Gamma^l_{phi theta} V_l
	break ; 
      }
				
      default : {
	cout << "Connection_fspher::derive_cov : index problem ! "
	     << endl ; 
	abort() ;  
      }
      }

    }


  }


  // Derivation index = phi
  // ----------------------
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
		
    cresu = (uu(ind0)).srstdsdp() ;  // 1/(r sin(theta)) d/dphi 	
    bool dz4 = (cresu.get_dzpuis() == 4) ;
		
    // Loop on all the indices of uu
    for (int id=0; id<valence0; id++) {
		
      switch ( ind0(id) ) {
				
      case 1 : {	// Gamma^r_{l phi} V^l 
	// or -Gamma^l_{r phi} V_l 
	ind = ind0 ; 
	ind.set(id) = 3 ;   // l = phi
	tmp = uu(ind) ; 
	dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	cresu -= tmp ; 
	break ; 
      }
				
      case 2 : {	// Gamma^theta_{l phi} V^l 
	// or -Gamma^l_{theta phi} V_l
	ind = ind0 ; 
	ind.set(id) = 3 ;   // l = phi
	tmp = uu(ind) ; 
	dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	tmp.div_tant() ; 	// division by tan(theta)
					
	cresu -= tmp ; 
	break ; 
      }
				
      case 3 : {	// Gamma^phi_{l phi} V^l 
	// or -Gamma^l_{phi phi} V_l
							
	ind = ind0 ; 

	ind.set(id) = 1 ;   // l = r
	tmp = uu(ind) ; 
	dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	cresu += tmp ; 

	ind.set(id) = 2 ;   // l = theta
	tmp = uu(ind) ; 
	dz4 ? tmp.div_r() : tmp.div_r_inc2() ; 
	tmp.div_tant() ; 	// division by tan(theta)

	cresu += tmp ; 
	break ; 
      }
				
      default : {
	cout << "Connection_fspher::derive_cov : index problem ! "
	     << endl ; 
	abort() ;  
      }
      }

    }


  }



  // C'est fini !
  // -----------
  return resu ; 

}

// Divergence, returning a pointer.
//---------------------------------

Tensor* Connection_fspher::p_divergence(const Tensor& uu) const {

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

  // Booleans to know the dynamical type (better than typeid...)
  //------------------------------------------------------------
  //##  bool bvec = false ;
  //##  int sym_flag = 1 ;

  // If u is a Vector, the result is a Scalar
  //----------------------------------------
  if (valence0 == 1) {
    //    bvec = true ;
    resu = new Scalar(*mp) ;
  }
  else {
    const Itbl tipeuu = uu.get_index_type() ;  
    for (int id = 0; id<valence0-1; id++) {
      tipe.set(id) = tipeuu(id+1) ; 
    }
    if (valence0 == 2) {
      resu = new Vector(*mp, tipe(0), *triad) ;
//       const Sym_tensor* sym_uu 
// 	= dynamic_cast<const Sym_tensor*>(&uu) ;
//       if (sym_uu != 0x0) sym_flag = 2 ;
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
	
  Itbl ind(valence0) ; // working Itbl to store the indices of uu
	
  Scalar tmp1(*mp) ;	// working scalar
  Scalar tmp2(*mp) ;	// working scalar

	
  // Loop on all the components of the output tensor
  for (int ic=0; ic<ncomp1; ic++) {
	
    ind1 = resu->indices(ic) ; 
    Scalar& cresu = resu->set(ind1) ;

    // Derivation index = r
    // --------------------
    int k = 1 ; 	

    // indices (k,ind1) in the input tensor
    ind0.set(0) = k ; 
    for (int id = 1; id<valence0; id++) {
      ind0.set(id) = ind1(id-1) ; 
    }

    cresu = uu(ind0).dsdr() ; //dT^{r l}/dr

  // Derivation index = theta
  // ------------------------
    k = 2 ; 	

    // indices (k,ind1) in the input tensor
    ind0.set(0) = k ; 
    for (int id = 1; id<valence0; id++) {
      ind0.set(id) = ind1(id-1) ; 
    }
		
    tmp1 = uu(ind0).dsdt() ; //dT^{theta l} /dtheta

    ind = ind0 ;
    ind.set(0) = 1 ;
    tmp1 += uu(ind) ;//##Gamma^theta_{r theta}T^{r l} (div_r later)+sym_flag
    

    // Loop on all the indices of uu
    for (int id=1; id<valence0; id++) {
		
      switch ( ind0(id) ) {
      case 1 : {	// Gamma^r_{l theta} V^l 
	// or -Gamma^l_{r theta} V_l 
	ind = ind0 ; 
	ind.set(id) = 2 ;   // l = theta
	tmp1 -= uu(ind) ; 
	break ; 
      }
				
      case 2 : {	// Gamma^theta_{l theta} V^l 
	// or -Gamma^l_{theta theta} V_l
	ind = ind0 ; 
	ind.set(id) = 1 ;   // l = r
	tmp1 += uu(ind) ; 
	break ; 
      }
				
      case 3 : {	// Gamma^phi_{l theta} V^l 
	// or -Gamma^l_{phi theta} V_l
	break ; 
      }
				
      default : {
	cout << "Connection_fspher::divergence : index problem ! "
	     << endl ; 
	abort() ;  
      }
      }
    }

  // Derivation index = phi
  // ----------------------
    k = 3 ; 			
    // indices (k,ind1) in the input tensor
    ind0.set(0) = k ; 
    for (int id = 1; id<valence0; id++) {
      ind0.set(id) = ind1(id-1) ; 
    }
    
    tmp1 += uu(ind0).stdsdp() ; // 1/sin(theta) dT^phi / dphi
    
    ind = ind0 ;
    ind.set(0) = 1 ;
    tmp1 += uu(ind) ;//##Gamma^phi_{r phi}T^{r l} (div_r later)+sym_flag
    ind.set(0) = 2 ;
    tmp2 = uu(ind) ;//##Gamma^phi_{theta phi}T^{theta l} (div_r later)+sym_flag

    // Loop on all the indices of uu
    for (int id=1; id<valence0; id++) {
      
      switch ( ind0(id) ) {
      case 1 : {	// Gamma^r_{l phi} V^l 
	// or -Gamma^l_{r phi} V_l 
	ind = ind0 ; 
	ind.set(id) = 3 ;   // l = phi
	tmp1 -= uu(ind) ; 
	break ; 
      }
				
      case 2 : {	// Gamma^theta_{l phi} V^l 
	// or -Gamma^l_{theta phi} V_l
	ind = ind0 ; 
	ind.set(id) = 3 ;   // l = phi
	tmp2 -= uu(ind) ; 
	break ; 
      }
				
      case 3 : {	// Gamma^phi_{l phi} V^l 
	// or -Gamma^l_{phi phi} V_l
							
	ind = ind0 ; 

	ind.set(id) = 1 ;   // l = r
	tmp1 += uu(ind) ; 

	ind.set(id) = 2 ;   // l = theta
	tmp2 += uu(ind) ; 
	break ; 
      }
				
      default : {
	cout << "Connection_fspher::divergence : index problem ! "
	     << endl ; 
	abort() ;  
      }
      }
    }
    // There remains a division by tan(theta) and r:
    //----------------------------------------------
    tmp2.div_tant() ;
    tmp1 += tmp2 ;
    tmp1.div_r_inc2() ;

    cresu += tmp1 ; // the d/dr term...

  }

  // C'est fini !
  // -----------
  return resu ; 

}







