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


char et_bin_kinema_nsbhC[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:52:56  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.2  2005/10/18 13:12:33  p_grandclement
 * update of the mixted binary codes
 *
 * Revision 1.1  2005/08/31 09:32:50  p_grandclement
 * Forgot this file
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "et_bin_nsbh.h"
#include "graphique.h"

namespace Lorene {
void Et_bin_nsbh::kinematics(double omega, double) {

    int nz = mp.get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 
    
    // --------------------
    // Computation of B^i/N
    // --------------------
    
    //  1/ Computation of  - omega m^i

    const Coord& xa = mp.xa ; 
    const Coord& ya = mp.ya ; 

    bsn.set_etat_qcq() ; 
    
    bsn.set(0) = -omega * ya ;
    bsn.set(1) = omega * xa ;
    bsn.set(2) = 0 ;
    
    bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC
   
    //	2/ Addition of shift and division by lapse
    // See Eq (47) from Gourgoulhon et al. (2001)

    bsn = -( bsn + shift ) / nnn ;
    
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
			    xa *xa + ya * ya ) ; 

    }
    
    pot_centri.annule(nzm1, nzm1) ;	// set to zero in the external domain
    pot_centri.set_std_base() ;   // set the bases for spectral expansions
    
    // The derived quantities are obsolete
    // -----------------------------------
    
    del_deriv() ;                
    
    
}
}
