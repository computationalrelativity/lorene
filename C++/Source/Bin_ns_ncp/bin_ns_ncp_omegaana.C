/*
 * Methods of class Bin_ns_ncp to set analytical value to omega
 *
 */

/*
 *   Copyright (c) 2003 Francois Limousin
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


 

/*
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "bin_ns_ncp.h"
#include "unites.h"


namespace Lorene {
void Bin_ns_ncp::analytical_omega() {
  
  using namespace Unites ;
    
    double rr = separation() ;
    double mtot = star1.mass_g() + star2.mass_g() ; 

    // Compacity factor
    double compact = ggrav * mtot / rr ; 

    // The compacity factor is set to zero in the Newtonian case
    if ( !star1.is_relativistic() ) {
	assert( !star2.is_relativistic() ) ; 
	compact = 0 ; 
    }
    
	
    double omega2 ; 
    
    if ( star1.is_irrotational() ) {
    
	    // Irrotational case
	    // -----------------
	
	assert( star2.is_irrotational() ) ; 
	
	omega2 = ggrav * mtot / pow(rr, 3) 
		* (1. - 2.75 * compact + 8.625 * compact*compact ) ; 
	
    }
    else{   // Corotating case
	    // ---------------
	
	assert( !star2.is_irrotational() ) ; 
    
	// a0/R
	double a0sr = star1.ray_eq() / rr ; 

	// Rescaled moment of inertia 5 I / (2 M a0^2)
	double ired = double(5)/double(3) * ( 1. - double(6) / M_PI / M_PI ) ; 
	omega2 =  ggrav * mtot / pow(rr, 3)
	    * (1. - compact * ( 2.75 + 2.*a0sr*a0sr * ired 
		  - 0.48*pow(a0sr, 4) * ired*ired ) 
		  + compact*compact * ( 8.625 + 2.75*a0sr*a0sr * ired
					+ 2.*pow(a0sr, 4) * ired*ired ) ) ;
    
    }    
    
    omega = sqrt( omega2 ) ; 
    
    // The derived quantities are obsolete:
    del_deriv() ; 
    
}
}
