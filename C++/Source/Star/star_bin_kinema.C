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
 * Revision 1.6  2005/02/17 17:33:54  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.5  2004/06/22 12:51:59  f_limousin
 * Simplify the computation of gam_pot to improve the convergence of the code.
 *
 * Revision 1.4  2004/05/25 14:21:26  f_limousin
 * New method to compute pot_centri to improve the convergence of the code.
 *
 * Revision 1.3  2004/02/27 09:56:10  f_limousin
 * Correction of an error on the computation of bsn.
 *
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

    bsn.change_triad(mp.get_bvect_cart()) ;

    if (fabs(mp.get_rot_phi()) < 1e-10){ 
      bsn.set(1) =  omega * ya ;
      bsn.set(2) = - omega * (xa - x_axe) ;
      bsn.set(3) = 0 ;
    }
    else {
      bsn.set(1) = - omega * ya ;
      bsn.set(2) = omega * (xa - x_axe) ;
      bsn.set(3) = 0 ;
    }

    bsn.std_spectral_base() ;
 
    bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC
    
    //	2/ Addition of shift and division by lapse
    // See Eq (47) from Gourgoulhon et al. (2001)
 
    bsn.change_triad(mp.get_bvect_spher()) ;
    bsn = ( bsn - beta ) / nn ; 
    bsn.change_triad(mp.get_bvect_cart()) ;

    bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC
        
    //-------------------------
    // Centrifugal potential
    //-------------------------
    
    // Lorentz factor between the co-orbiting observer and the Eulerian one
    // See Eq (23) from Gourgoulhon et al. (2001)

    Sym_tensor flat_cov = flat.cov() * psi4 ;
    flat_cov.change_triad(mp.get_bvect_cart()) ;
    Sym_tensor gamma_cov (gamma.cov()) ;
    gamma_cov.change_triad(mp.get_bvect_cart()) ;

    //## For the convergence of the code, we introduce a flat scalar 
    // product to compute gam_pot and so pot_centri. 
    Scalar gam_pot = 1 / sqrt( 1 - contract(flat_cov, 0, 1, bsn * bsn, 0, 1)) ;

    Scalar gam0 = 1 / sqrt(1 - contract(gamma_cov, 0, 1, bsn * bsn, 0, 1)) ;
    
    // Relative error make when we take gam_pot instead of gam0 for 
    // the computation of pot_centri. 

    cout << "gam_pot" << endl << norme(gam_pot) << endl ;
    cout << "gam0" << endl << norme(gam0) << endl ;
    cout << "Relative difference between gam0 and gam_pot : " << endl ; 
    cout << diffrel(gam0, gam_pot) << endl ;

    pot_centri = - log( gam_pot ) ;
    
    pot_centri.annule(nzm1, nzm1) ;	// set to zero in the external domain
    pot_centri.std_spectral_base() ;   // set the bases for spectral expansions
    
      // The derived quantities are obsolete
      // -----------------------------------
      
      del_deriv() ;                
      
}
