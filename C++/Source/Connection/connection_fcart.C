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
					//	       Constructors         //
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
					//	       Destructor         //
					//----------------------------//


Connection_fcart::~Connection_fcart(){
	
}


					//-----------------------------//
    				//     Mutators / assignment   //
					//-----------------------------//


void Connection_fcart::operator=(const Connection_fcart& ci) {
	
	cout << "Connection_fcart::operator= : not implemented yet !" << endl ; 
	abort() ; 

}	



					//-----------------------------//
					//    Computational methods    //
					//-----------------------------//



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
	tipe.set_etat_qcq() ; 
	tipe.set(0) = COV ; 
	const Itbl tipeuu = uu.get_index_type() ;  
	for (int id = 1; id<=valence0; id++) {
		tipe.set(id) = tipeuu(id-1) ; 
	}

	// Creation of the result tensor
	// -----------------------------
	Tensor resu(*mp, valence0+1, tipe, *triad) ;
	
	Itbl ind1(valence0+1) ; // working Itbl to store the indices of resu
	ind1.set_etat_qcq() ; 
	
	Itbl ind0(valence0) ; // working Itbl to store the indices of uu
	ind0.set_etat_qcq() ; 
	

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









