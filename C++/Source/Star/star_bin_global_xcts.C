/*
 * Methods for computing global quantities within the class Star_bin_xcts
 * (see file star.h for documentation)
 */

/*
 *   Copyright (c) 2010 Michal Bejger
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

char star_bin_global_xcts_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/13 08:53:38  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2010/05/04 07:51:05  m_bejger
 * Initial version
 *
 * $Header$
 *
 */

// Headers Lorene
#include "star.h"
#include "utilitaires.h"

			//--------------------------//
			//	Baryon mass	            //
			//--------------------------//

namespace Lorene {
double Star_bin_xcts::mass_b() const {

    if (p_mass_b == 0x0) {    // a new computation is required

	Scalar det_gamma = gamma.determinant() ;
	
	Scalar dens = sqrt(det_gamma) * gam_euler * nbar ;
	
	dens.std_spectral_base() ; 
	
	p_mass_b = new double( dens.integrale() ) ;
	
    }
    
    return *p_mass_b ; 

} 



			//----------------------------//
			//	Gravitational mass        //
			//----------------------------//

double Star_bin_xcts::mass_g() const {

    if (p_mass_g == 0x0) {    // a new computation is required
	
	Scalar det_gamma = gamma.determinant() ;

	Scalar dens = sqrt(det_gamma) * nn
	    * ( ener_euler + s_euler ) ;
	
	dens.std_spectral_base() ; 
	
	p_mass_g = new double( dens.integrale() ) ;

    }
    
    return *p_mass_g ; 

} 

		
			//----------------------------------//
			//  X coordinate of the barycenter  //
			//----------------------------------//


double Star_bin_xcts::xa_barycenter() const {

    if (p_xa_barycenter == 0x0) {    // a new computation is required
	
	Scalar xxa(mp) ; 
	xxa = mp.xa ;	// Absolute X coordinate
	xxa.std_spectral_base() ;

	Scalar det_gamma = gamma.determinant() ;

	Scalar dens = sqrt(det_gamma) * gam_euler * nbar * xxa ; 
	
	int nzone = mp.get_mg()->get_nzone() ;
	dens.annule_domain(nzone - 1) ;

	dens.std_spectral_base() ; 

	p_xa_barycenter = new double( dens.integrale() / mass_b() ) ;
	
    }
    
    return *p_xa_barycenter ; 

}
}
