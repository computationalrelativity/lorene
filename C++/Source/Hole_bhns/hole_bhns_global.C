/*
 *  Methods of class Hole_bhns to compute global quantities
 *
 *    (see file hole_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Keisuke Taniguchi
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

char hole_bhns_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/06/22 01:25:15  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
#include <math.h>

// Lorene headers
#include "hole_bhns.h"
#include "unites.h"
#include "utilitaires.h"

                    //-----------------------------------------//
                    //          Irreducible mass of BH         //
                    //-----------------------------------------//

double Hole_bhns::mass_irr_bhns() const {

    // Fundamental constants and units
    // -------------------------------
    using namespace Unites ;

    if (p_mass_irr_bhns == 0x0) {   // a new computation is required

        Scalar psi4(mp) ;
	psi4 = pow(confo_tot, 4.) ;
	psi4.std_spectral_base() ;
	psi4.annule_domain(0) ;
	psi4.raccord(1) ;

	double radius_ah = mp.val_r(1,-1.,M_PI/2.,0.) ;

	Map_af& mp_aff= dynamic_cast<Map_af&>(mp) ;

	double a_ah = mp_aff.integrale_surface(psi4, radius_ah) ;
	double mirr = sqrt(a_ah/16./M_PI) / ggrav ;

	p_mass_irr_bhns = new double( mirr ) ;

    }

    return *p_mass_irr_bhns ;

}
