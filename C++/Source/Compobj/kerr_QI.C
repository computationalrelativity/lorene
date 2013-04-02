/*
 *  Methods of the class Kerr_QI
 *
 *    (see file compobj.h for documentation).
 *
 */

/*
 *   Copyright (c) 2013 Claire Some, Eric Gourgoulhon
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

char kerr_QI_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2013/04/02 23:17:18  e_gourgoulhon
 * New class Kerr_QI
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
Kerr_QI::Kerr_QI(Map& mpi, double mass, double a_over_m) :
			 Compobj_QI(mpi)
{
    mm = mass ; 
    aa = a_over_m * mm ; 
    double mm2 = mm*mm ; 
    double aa2 = aa*aa ; 
    double hh = sqrt(mm2-aa2) ; // Eq. (98)
    double hh2 = hh*hh ; 
    double r_hor = hh / double(2) ; // Eq. (10)

    const Coord& r = mp.r ;        // r field 
    const Coord& cost = mp.cost ;  // cos(theta) field
    const Coord& sint = mp.sint ;  // sin(theta) field
    Mtbl r2 = r*r ; 
    Mtbl cost2 = cost*cost ; 
    Mtbl sint2 = sint*sint ; 

    // A^2
    a_car = 1 + 2*mm/r + (3*mm2 + aa2*(2*cost2-1))/(2*r2) + hh2*mm/(2*r*r2) + 
            hh2*hh2/(16*r2*r2) ;  // Eq. (121)
    a_car.set_domain(0) = 1 ; 
    a_car.std_spectral_base() ;
    
    // Boyer-Lindquist radial coordinate:
    Mtbl rBL = r + hh2/(4*r) + mm ;  // Eq. (110)
    Mtbl rBL2 = rBL*rBL ; 
    Mtbl sigma = rBL2 + aa2*cost2 ;  // Eq. (93)
    
    // B^2
    b_car = (rBL2 + aa2 + 2*aa2*mm*rBL*sint2 / sigma) / r2 ;  // Eq. (125)
    b_car.set_domain(0) = 1 ; 
    b_car.std_spectral_base() ;
    bbb = sqrt(b_car) ; 
    bbb.std_spectral_base() ;
    
    // N
    nn = (1 - r_hor*r_hor / r2) / bbb ;  // Eq. (81) + (10)
    nn.set_domain(0) = 1 ; 
    nn.std_spectral_base() ;

    // N^phi
    nphi = 2*aa*mm / (sigma*(rBL+aa2/rBL) + 2*aa2*mm*sint2) ;  // Eq. (126)
    nphi.set_domain(0) = 0 ; 
    nphi.std_spectral_base() ;
    
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}

// Copy constructor
// --------------------
Kerr_QI::Kerr_QI(const Kerr_QI& other) :
			 Compobj_QI(other),
             mm(other.mm),
             aa(other.aa)		 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}


// Constructor from a file
// -----------------------
Kerr_QI::Kerr_QI(Map& mpi, FILE* fich) :
			 Compobj_QI(mpi)
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
    
    // Read of the saved fields:
    // ------------------------

    
}

			    //------------//
			    // Destructor //
			    //------------//

Kerr_QI::~Kerr_QI(){

    del_deriv() ; 

}


			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Kerr_QI::del_deriv() const {

    Compobj_QI::del_deriv() ; 


    Kerr_QI::set_der_0x0() ; 
}			    


void Kerr_QI::set_der_0x0() const {
 	 
}			    

			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Kerr_QI
// --------------------------------
void Kerr_QI::operator=(const Kerr_QI& other) {

    // Assignment of quantities common to all the derived classes of Compobj_QI
    Compobj_QI::operator=(other) ;	    
    
    mm = other.mm ; 
    aa = other.aa ; 

    del_deriv() ;  // Deletes all derived quantities
}	

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Kerr_QI::sauve(FILE* fich) const {

    
}

// Printing
// --------

ostream& Kerr_QI::operator>>(ostream& ost) const {

 	using namespace Unites ;
	
	Compobj_QI::operator>>(ost) ; 
    
    ost << endl << "Kerr spacetime in quasi-isotropic coordinates (class Kerr_QI) " << endl ; 

    ost << "M = " << mm << "  a = " << aa << endl ; 
    
    return ost ; 
      
}

			    //-------------------------//
			    //	Computational methods  //
			    //-------------------------//
			    
