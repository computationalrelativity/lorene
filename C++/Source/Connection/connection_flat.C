/*
 *  Methods of class Connection_flat.
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

char connection_flat_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2003/10/11 14:39:50  e_gourgoulhon
 * Suppressed declaration of unusued arguments in some methods.
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


					//------------------------------//
					//	       Constructors         //
					//------------------------------//



// Constructor for derived classes

Connection_flat::Connection_flat(const Map& mpi, const Base_vect& bi) 
	: Connection(mpi, bi) {

	assoc_metric = true ;
	
	delta.set_etat_zero() ; 

}		

// Copy constructor
Connection_flat::Connection_flat(const Connection_flat& ci) 
	: Connection(ci) {

}		

	
					//----------------------------//
					//	       Destructor         //
					//----------------------------//


Connection_flat::~Connection_flat(){

}


					//-----------------------------//
    				//     Mutators / assignment   //
					//-----------------------------//


void Connection_flat::operator=(const Connection_flat& ) {
	
	cout << "Connection_flat::operator= : not implemented yet !" << endl ; 
	abort() ; 

}	



					//-----------------------------//
					//    Computational methods    //
					//-----------------------------//


void Connection_flat::compute_ricci() const {

	assert(p_ricci == 0x0) ;
	
	p_ricci = new Sym_tensor(*mp, COV, *triad) ;
	
	p_ricci->set_etat_zero() ; 
	
}












