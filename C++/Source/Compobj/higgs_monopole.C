/*
 *  Methods of the class HiggsMonopole
 *
 *    (see file compobj.h for documentation).
 *
 */

/*
 *   Copyright (c) 2014 Marie Leroy,  Eric Gourgoulhon
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

char HiggsMonopole_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2014/01/29 16:31:42  e_gourgoulhon
 * New class HiggsMonopole
 *
 *
 * $Header$
 *
 */


// C headers
#include <cassert>

// Lorene headers
#include "compobj.h"
#include "unites.h"

                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
// --------------------
HiggsMonopole::HiggsMonopole(Map& mpi, const char* file_name) :
			 Compobj(mpi) , 
			 hh(mpi)
{

    //ifstream file(file_name) ; 
    //if ( !file.good() ) {
        //cerr << "Problem in opening the file " << file_name << endl ;
        //abort() ;
    //}
    
    //file.getline(description1, 256) ;
    //file.getline(description2, 256) ;
    //cout << "description1 : " << description1 << endl ; 
    //cout << "description2 : " << description2 << endl ; 
//     description1[0] = " " ; 
//     description2[0] = " " ; 
//     cout << "description1 : " << description1 << endl ; 
//     cout << "description2 : " << description2 << endl ; 

}

// Copy constructor
// --------------------
HiggsMonopole::HiggsMonopole(const HiggsMonopole& other) :
			 Compobj(other),
			 hh(other.hh)
{
    // Pointers of derived quantities initialized to zero : 
    // set_der_0x0() ;
}

			    //------------//
			    // Destructor //
			    //------------//

HiggsMonopole::~HiggsMonopole(){

    // del_deriv() ; 

}

// Printing
// --------

ostream& HiggsMonopole::operator>>(ostream& ost) const {

 	using namespace Unites ;
	
    ost << endl << "Higgs monopole" << endl ; 

	Compobj::operator>>(ost) ; 

    return ost ; 
      
}



