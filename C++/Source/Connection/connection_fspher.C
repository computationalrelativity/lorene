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


void Connection_fspher::operator=(const Connection_fspher& ci) {
	
	cout << "Connection_fspher::operator= : not implemented yet !" << endl ; 
	abort() ; 

}	



					//-----------------------------//
					//    Computational methods    //
					//-----------------------------//



Tensor Connection_fspher::derive_cov(const Tensor&) const {

	cout << "Connection_fspher::derive_cov: not implemented yet !" << endl ; 
	abort() ; 

}









