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
 * Revision 1.6  2003/10/11 14:39:49  e_gourgoulhon
 * Suppressed declaration of unusued arguments in some methods.
 *
 * Revision 1.5  2003/10/06 13:58:46  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.4  2003/10/03 14:16:04  e_gourgoulhon
 * Added set_der_0x0 in some constructors.
 *
 * Revision 1.3  2003/10/02 21:32:06  e_gourgoulhon
 * Added constructor from Metric.
 * Added functions fait_delta and update.
 *
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
#include "metric.h"


					//------------------------------//
					//	       Constructors         //
					//------------------------------//


// Constructor ab initio

Connection::Connection(const Tensor_delta& delta_i) : mp(&(delta_i.get_mp())),
		triad(delta_i.get_triad()),
		delta(delta_i), 
		assoc_metric(false),
		flat_conn(0x0) {
	
	const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ; 
	if (bvs != 0x0) {
		flat_conn = new Connection_fspher(*mp, *bvs) ; 
	}
	else{
		const Base_vect_cart* bvc = dynamic_cast<const Base_vect_cart*>(triad) ;
		if (bvc == 0x0) {
			cout << "Connection::Connection : unknown type of triad !" << endl ;
			abort() ; 			
		}
		//##
		// flat_conn = new Connection_fcart(*mp, *bvc) ; 
	}
	
	
	set_der_0x0() ; 

	
}		


// Standard constructor from a metric. 

Connection::Connection(const Metric& met) : mp(&(met.get_mp())),
		triad(met.cov().get_triad()),
		delta(*mp, CON, COV, COV, *triad),
		assoc_metric(true),
		flat_conn(0x0) {
		
	const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ; 
	if (bvs != 0x0) {
		flat_conn = new Connection_fspher(*mp, *bvs) ; 
	}
	else{
		const Base_vect_cart* bvc = dynamic_cast<const Base_vect_cart*>(triad) ;
		if (bvc == 0x0) {
			cout << "Connection::Connection : unknown type of triad !" << endl ;
			abort() ; 			
		}
		//##
		// flat_conn = new Connection_fcart(*mp, *bvc) ; 
	}

	fait_delta(met) ; 	// Computes delta

	set_der_0x0() ; 
}


// Copy constructor

Connection::Connection(const Connection& conn_i) : mp(conn_i.mp),
		triad(conn_i.triad),
		delta(conn_i.delta), 
		assoc_metric(conn_i.assoc_metric),
		flat_conn(0x0) {
		
	if (conn_i.flat_conn != 0x0) {
		const Connection_fspher* cfs = 
			dynamic_cast<const Connection_fspher*>(conn_i.flat_conn) ; 
		if (cfs != 0x0) {
			flat_conn = new Connection_fspher(*cfs) ; 			
		}
		//## cas Connection_fcart
	}
	
	set_der_0x0() ; 
	
}		


// Constructor for derived classes

Connection::Connection(const Map& mpi, const Base_vect& bi) : mp(&mpi),
		triad(&bi),
		delta(mpi, CON, COV, COV, bi),
		assoc_metric(false),
		flat_conn(0x0){
		
	set_der_0x0() ; 
	
}		


	
					//----------------------------//
					//	       Destructor         //
					//----------------------------//


Connection::~Connection(){

	if (flat_conn != 0x0) delete flat_conn ; 

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
	if (flat_conn != 0x0) delete flat_conn ; 
	flat_conn = 0x0 ; 

	if (ci.flat_conn != 0x0) {
		const Connection_fspher* cfs = 
			dynamic_cast<const Connection_fspher*>(ci.flat_conn) ; 
		if (cfs != 0x0) {
			flat_conn = new Connection_fspher(*cfs) ; 			
		}
		//## cas Connection_fcart
	}
	
	del_deriv() ; 

}	



					//-----------------------------//
					//    Computational methods    //
					//-----------------------------//


Tensor Connection::derive_cov(const Tensor& ti) const {

	cout << "Connection::derive_cov : not implemented yet !" << endl ; 
	abort() ; 
	return ti ;

} 

Tensor* Connection::p_derive_cov(const Tensor& ) const {

	cout << "Connection::p_derive_cov : not implemented yet !" << endl ; 
	abort() ; 
	return 0x0 ;
} 


const Tensor& Connection::ricci() const {

	if (p_ricci == 0x0) {
		compute_ricci() ; 
	}
	
	return *p_ricci ; 
	
}


void Connection::update(const Tensor_delta& delta_i) {

	assert(assoc_metric == false) ;
	
	assert(flat_conn != 0x0) ; // to guarantee we are not in a derived class
	
	delta = delta_i ; 
	
	del_deriv() ; 
	
}


void Connection::update(const Metric& met) {

	assert(assoc_metric == true) ;
	
	assert(flat_conn != 0x0) ; // to guarantee we are not in a derived class
	
	fait_delta(met) ; 
	
	del_deriv() ; 
	
}

void Connection::compute_ricci() const {

	cout << "Connection::compute_ricci() : not implemented yet !" << endl ; 
	abort() ; 

}


void Connection::fait_delta(const Metric& ) {

	cout << "Connection::fait_delta : not implemented yet !" << endl ; 
	abort() ; 

}  












