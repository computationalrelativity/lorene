/*
 *  Methods of class Connection.
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

char connection_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/10/01 15:42:49  e_gourgoulhon
 * still ongoing...
 *
 * Revision 1.1  2003/09/29 21:13:08  e_gourgoulhon
 * First version --- not ready yet.
 *
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


// Constructor ab initio

Connection::Connection(const Tensor_delta& delta_i) : mp(&(delta_i.get_mp())),
		triad(delta_i.get_triad()),
		delta(delta_i), 
		assoc_metric(false),
		flat_conn(0x0) {	//## flat_conn((delta_i.get_triad())->connect(delta_i.get_mp())
	
	set_der_0x0() ; 

	
}		


// Copy constructor

Connection::Connection(const Connection& conn_i) : mp(conn_i.mp),
		triad(conn_i.triad),
		delta(conn_i.delta), 
		assoc_metric(conn_i.assoc_metric),
		flat_conn(conn_i.flat_conn) {
	
	set_der_0x0() ; //## a voir
	
}		


// Constructor for derived classes

Connection::Connection(const Map& mpi, const Base_vect& bi) : mp(&mpi),
		triad(&bi),
		delta(mpi, CON, COV, COV, bi),
		assoc_metric(false),
		flat_conn(0x0){
		
	
}		


	
					//----------------------------//
					//	       Destructor         //
					//----------------------------//


Connection::~Connection(){

	del_deriv() ; 
	
}

					//-----------------------------//
					//        Memory management    //
					//-----------------------------//

void Connection::del_deriv() const {

	if (p_ricci != 0x0) delete p_ricci ; 
	
	set_der_0x0() ; 
	
}

void Connection::set_der_0x0() const {

	p_ricci = 0x0 ; 
	
}


					//-----------------------------//
    				//     Mutators / assignment   //
					//-----------------------------//


void Connection::operator=(const Connection& ci) {
	
	assert( triad == ci.triad ) ; 
	delta = ci.delta ; 
	
	del_deriv() ; 

}	



					//-----------------------------//
					//    Computational methods    //
					//-----------------------------//


Tensor Connection::derive_cov(const Tensor&) const {

	cout << "Connection::derive_cov : not implemented yet !" << endl ; 
	abort() ; 

} 



const Tensor& Connection::ricci() const {

	if (p_ricci == 0x0) {
		compute_ricci() ; 
	}
	
	return *p_ricci ; 
	
}


void Connection::compute_ricci() const {

	cout << "Connection::compute_ricci() : not implemented yet !" << endl ; 
	abort() ; 

}












