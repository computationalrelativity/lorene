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


 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2016/12/05 16:18:14  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:53:38  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2005/09/13 19:38:31  f_limousin
 * Reintroduction of the resolution of the equations in cartesian coordinates.
 *
 * Revision 1.4  2005/02/17 17:33:25  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.3  2004/02/27 09:54:24  f_limousin
 * Generalisation of the formulas for mass_b and mass_g for non
 * conformally flat metrics.
 *
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
#include "utilitaires.h"

			//--------------------------//
			//	Baryon mass	    //
			//--------------------------//

namespace Lorene {
double Star_bin::mass_b() const {

    if (p_mass_b == 0x0) {    // a new computation is required

	Scalar det_gamma = gamma.determinant() ;
	
	Scalar dens = sqrt(det_gamma) * gam_euler * nbar ;

//	Scalar dens = psi4 * sqrt(psi4) * gam_euler * nbar ;
	
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


double Star_bin::xa_barycenter() const {

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
