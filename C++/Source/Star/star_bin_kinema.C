/*
 * Method Star_bin::kinematics
 *
 * (see file star.h for documentation)
 *
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


char star_bin_kinema_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/01/20 15:18:45  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

#include <math.h>

// Headers Lorene
#include "star.h"

void Star_bin::kinematics(double omega, double x_axe) {

    int nz = mp.get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 
    
    // --------------------
    // Computation of B^i/N
    // --------------------
    
    //  1/ Computation of  - omega m^i

    const Coord& xa = mp.xa ; 
    const Coord& ya = mp.ya ; 

    const Mtbl& ra = sqrt(xa * xa + ya * ya) ;

    bsn.set_etat_qcq() ; 

    if (fabs(mp.get_rot_phi()) < 1e-10){ 
      bsn.set(1) =  0 ;
      bsn.set(2) = - omega * ra ;
      bsn.set(3) = 0 ;
    }
    else {
      bsn.set(1) = 0 ;
      bsn.set(2) = omega * ra ;
      bsn.set(3) = 0 ;
    }


    bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC
    
    //	2/ Addition of shift and division by lapse
    // See Eq (47) from Gourgoulhon et al. (2001)

    bsn = ( bsn + shift ) / nnn ; 
  
    bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC
    bsn.std_spectral_base() ;   // set the bases for spectral expansions
        
    //-------------------------
    // Centrifugal potential
    //-------------------------

    // Lorentz factor between the co-orbiting observer and the Eulerian one
    // See Eq (23) from Gourgoulhon et al. (2001)

      Scalar gam0 = 1 / sqrt( 1-sprod(bsn, bsn) ) ;
	
      pot_centri = - log( gam0 ) ;

      pot_centri.annule(nzm1, nzm1) ;	// set to zero in the external domain
      pot_centri.std_spectral_base() ;   // set the bases for spectral expansions
      
      // The derived quantities are obsolete
      // -----------------------------------
      
      del_deriv() ;                
      
}
