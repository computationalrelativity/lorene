/*
 * Methods for computing global quantities within the class Star_bin
 *
 * (see file star.h for documentation)
 */

/*
 *   Copyright (c) 2004 Francois Limousin
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


char star_bin_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/01/20 15:18:17  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// Headers C

// Headers Lorene
#include "star.h"

			//--------------------------//
			//	Baryon mass	    //
			//--------------------------//

double Star_bin::mass_b() const {

    if (p_mass_b == 0x0) {    // a new computation is required
	

	// Works only if gamma is conformally flat
	Scalar a_car(gamma.cov()(1,1)) ;

	Scalar dens = pow(a_car, 3./2.) * gam_euler * nbar ;
	
	dens.std_spectral_base() ; 
	
	p_mass_b = new double( dens.integrale() ) ;
	
    }
    
    return *p_mass_b ; 

} 



			//----------------------------//
			//	Gravitational mass    //
			//----------------------------//

double Star_bin::mass_g() const {

    if (p_mass_g == 0x0) {    // a new computation is required
	
	// Works only if gamma is conformally flat
	Scalar a_car(gamma.cov()(1,1)) ;

	Scalar dens = pow(a_car, 3./2.) * nnn
	    * ( ener_euler + s_euler ) ;
	
	dens.std_spectral_base() ; 
	
	p_mass_g = new double( dens.integrale() ) ;

    }
    
    return *p_mass_g ; 

} 

		
			//----------------------------------//
			//  X coordinate of the barycenter  //
			//----------------------------------//


double Star_bin::xa_barycenter() const {

    if (p_xa_barycenter == 0x0) {    // a new computation is required
	
	Scalar xxa(mp) ; 
	xxa = mp.xa ;	// Absolute X coordinate
	xxa.std_spectral_base() ;

	// Works only if gamma is conformally flat
	Scalar a_car(gamma.cov()(1,1)) ;

	Scalar dens = pow(a_car, 3./2.) * gam_euler * nbar * xxa ; 
	
	dens.std_spectral_base() ; 

	p_xa_barycenter = new double( dens.integrale() / mass_b() ) ;
	
    }
    
    return *p_xa_barycenter ; 

}

