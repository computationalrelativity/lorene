/*
 * Method Et_bin_ncp::kinematics
 *
 * (see file et_bin_ncp.h for documentation)
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char et_bin_ncp_kinema_C[] = "$Header$" ;

/*
 * $Header$
 *
 */

// Headers Lorene
#include "et_bin_ncp.h"

void Et_bin_ncp::kinematics(double omega, double x_axe) {

    int nz = mp.get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 
    
    // --------------------
    // Computation of B^i/N
    // --------------------
    
    //  1/ Computation of  - omega m^i

    const Coord& xa = mp.xa ; 
    const Coord& ya = mp.ya ; 

    bsn.set_etat_qcq() ; 

    if (fabs(mp.get_rot_phi()) < 1e-10){ 
      bsn.set(0) =  omega * ya ;
      bsn.set(1) = - omega * (xa - x_axe) ;
      bsn.set(2) = 0 ;
    }
    else {
      bsn.set(0) = - omega * ya ;
      bsn.set(1) = omega * (xa - x_axe) ;
      bsn.set(2) = 0 ;
    }


    bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC
    
    //	2/ Addition of shift and division by lapse
    // See Eq (47) from Gourgoulhon et al. (2001)

    bsn = ( bsn + shift ) / nnn ; 
  
    bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC
    bsn.set_std_base() ;   // set the bases for spectral expansions
        
    //-------------------------
    // Centrifugal potentatial
    //-------------------------

    if (relativistic) {

	// Lorentz factor between the co-orbiting observer and the Eulerian one
      // See Eq (23) from Gourgoulhon et al. (2001)

      Tenseur gam0 = 1 / sqrt( 1-sprod(bsn, bsn) ) ;
	
      pot_centri = - log( gam0 ) ;

    }
    else {

	pot_centri.set_etat_qcq() ; 
	

	// See Eq (40) from Gourgoulhon et al. (2001)
	pot_centri.set() = - 0.5 * omega * omega * (
			    (xa - x_axe) * (xa - x_axe) + ya * ya ) ; 

    }
    

     
    pot_centri.annule(nzm1, nzm1) ;	// set to zero in the external domain
    pot_centri.set_std_base() ;   // set the bases for spectral expansions
    
    // The derived quantities are obsolete
    // -----------------------------------
    
    del_deriv() ;                
    
}
