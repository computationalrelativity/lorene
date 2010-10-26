/*
 * Method Star_bin_xcts::kinematics
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

char star_bin_kinema_xcts_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2010/10/26 19:57:02  m_bejger
 * Various cleanups
 *
 * Revision 1.2  2010/06/04 19:57:52  m_bejger
 * Added condition for mp.get_rot_phi
 *
 * Revision 1.1  2010/05/04 07:51:05  m_bejger
 * Initial version
 *
 * $Header$
 *
 */

#include <math.h>

// Headers Lorene
#include "star.h"

void Star_bin_xcts::kinematics(double omega, double x_axe) {

    int nzm1 = mp.get_mg()->get_nzone() - 1 ; 

    // --------------------
    // Computation of B^i/N
    // --------------------
    
    //  1/ Computation of omega m^i

    const Coord& xa = mp.xa ; 
    const Coord& ya = mp.ya ; 

    bsn.change_triad(mp.get_bvect_cart()) ;

    if (fabs(mp.get_rot_phi()) < 1e-10){ 
      bsn.set(1) = - omega * ya ;
      bsn.set(2) = omega * (xa - x_axe) ;
      bsn.set(3) = 0 ;
    }
    else {
		
      bsn.set(1) = omega * ya ;
      bsn.set(2) = - omega * (xa - x_axe) ;
      bsn.set(3) = 0 ;
    }

    bsn.std_spectral_base() ;
    bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC
    
    // Addition of beta and division by lapse
    // Eq. 47 from Gourgoulhon et al. (2001)
    // New convention : l = Nn + B ==> B = \beta + \Omega d\phi
    //---------------------------------------------------------

    bsn = ( bsn + beta ) / nn ; 

    bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC
        
    //----------------------
    // Centrifugal potential
    //----------------------
    
    // Lorentz factor between the co-orbiting observer 
    // and the Eulerian one
    // Eq. 23 from Gourgoulhon et al. (2001)

    Sym_tensor gamma_cov (gamma.cov()) ;
    gamma_cov.change_triad(mp.get_bvect_cart()) ;

    Scalar gam0 = 1 / sqrt(1 - contract(gamma_cov, 0, 1, bsn * bsn, 0, 1)) ;
    pot_centri = - log( gam0 ) ;

    pot_centri.annule(nzm1, nzm1) ;	   	// set to zero in the external domain
    pot_centri.std_spectral_base() ;	// set the bases for spectral expansions
    
      // The derived quantities are obsolete
      // -----------------------------------
      
      del_deriv() ;                
      
}
