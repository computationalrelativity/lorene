/*
 *  Methods of the class Boson_star
 *
 *    (see file boson_star.h for documentation).
 */

/*
 *   Copyright (c) 2012 Claire Some, Eric Gourgoulhon
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

char boson_star_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2012/11/22 16:04:12  c_some
 * New class Boson_star
 *
 *
 * $Header$
 *
 */


// C headers
#include <cassert>

// Lorene headers
#include "boson_star.h"

                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
// --------------------
Boson_star::Boson_star(Map& mpi) :
			 Star_QI(mpi) ,
			 rphi(mpi), 
			 iphi(mpi)
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;

    // Initialization of the scalar field to zero
    rphi = 0 ;   
    iphi = 0 ;   
    
}

// Copy constructor
// --------------------
Boson_star::Boson_star(const Boson_star& st) :
			 Star_QI(st), 
			 rphi(st.rphi), 
			 iphi(st.iphi)		 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}


// Constructor from a file
// -----------------------
Boson_star::Boson_star(Map& mpi, FILE* fich) :
			 Star_QI(mpi, fich) , 
			 rphi(mpi, *(mpi.get_mg()), fich) ,
			 iphi(mpi, *(mpi.get_mg()), fich)
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
        
}

			    //------------//
			    // Destructor //
			    //------------//

Boson_star::~Boson_star(){

    del_deriv() ; 

}


			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Boson_star::del_deriv() const {

    Star_QI::del_deriv() ; 

    Boson_star::set_der_0x0() ; 
}			    


void Boson_star::set_der_0x0() const {
 	 
}			    

			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Boson_star
// --------------------------------
void Boson_star::operator=(const Boson_star& st) {

    // Assignment of quantities common to all the derived classes of Star_QI
    Star_QI::operator=(st) ;	    
    
    rphi = st.rphi ; 
    iphi = st.iphi ; 

    del_deriv() ;  // Deletes all derived quantities
}	

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Boson_star::sauve(FILE* fich) const {

	Star_QI::sauve(fich) ; 
    rphi.sauve(fich) ; 
    iphi.sauve(fich) ; 
}

// Printing
// --------

ostream& Boson_star::operator>>(ostream& ost) const {
   
    Star_QI::operator>>(ost) ; 
    
    ost << endl << "Axisymmetric stationary boson star in quasi-isotropic coordinates (class Boson_star) " << endl ;
    
    ost << "Central value of the scalar field : " << rphi.val_grid_point(0,0,0,0) << " + i " << iphi.val_grid_point(0,0,0,0)<< endl ;  

    ost << "Real part of the scalar field : " << rphi << endl ; 
    ost << "Imaginary part of the scalar field : " << iphi << endl ; 
    	
    return ost ; 
      
}



