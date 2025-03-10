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


 

/*
 * $Id$
 * $Log$
 * Revision 1.11  2016/12/05 16:18:14  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.10  2014/10/13 08:53:38  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.9  2014/10/06 15:13:17  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2006/04/11 14:24:44  f_limousin
 * New version of the code : improvement of the computation of some
 * critical sources, estimation of the dirac gauge, helical symmetry...
 *
 * Revision 1.7  2005/09/13 19:38:31  f_limousin
 * Reintroduction of the resolution of the equations in cartesian coordinates.
 *
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

#include <cmath>

// Headers Lorene
#include "star.h"

namespace Lorene {
void Star_bin::kinematics(double omega, double x_axe) {

    int nz = mp.get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 
    
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
    
    //	2/ Addition of beta and division by lapse
    // See Eq (47) from Gourgoulhon et al. (2001)
    // New convention : l = Nn + B ==> B = \beta + \Omega d\phi
 
    bsn = ( bsn + beta ) / nn ; 

    bsn.annule(nzm1, nzm1) ;	// set to zero in the ZEC
        
    //-------------------------
    // Centrifugal potential
    //-------------------------
    
    // Lorentz factor between the co-orbiting observer and the Eulerian one
    // See Eq (23) from Gourgoulhon et al. (2001)

    Sym_tensor gamma_cov (gamma.cov()) ;
    gamma_cov.change_triad(mp.get_bvect_cart()) ;

    Scalar gam0 = 1 / sqrt(1 - contract(gamma_cov, 0, 1, bsn * bsn, 0, 1)) ;
    pot_centri = - log( gam0 ) ;

    pot_centri.annule(nzm1, nzm1) ;	// set to zero in the external domain
    pot_centri.std_spectral_base() ;   // set the bases for spectral expansions
    
      // The derived quantities are obsolete
      // -----------------------------------
      
      del_deriv() ;                
      
}
}
