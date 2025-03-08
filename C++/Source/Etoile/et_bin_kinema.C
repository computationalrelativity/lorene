/*
 * Method Etoile_bin::kinematics
 *
 * (see file etoile.h for documentation)
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2016/12/05 16:17:53  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:52:56  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2005/08/29 15:10:16  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.4  2003/03/03 19:18:12  f_limousin
 * Suppression of bsn.set_triad(ref_triad).
 *
 * Revision 1.3  2003/01/17 13:35:48  f_limousin
 * Add comments and replace A^2*flat_scalar_prod by sprod
 *
 * Revision 1.2  2002/12/10 15:56:43  k_taniguchi
 * Change the multiplication "*" to "%".
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/02/17  18:51:22  eric
 * pot_centri is set to zero only in the external compactified domain.
 *
 * Revision 2.3  2000/02/12  18:37:17  eric
 * Appel de set_std_base sur les quantites calculees.
 *
 * Revision 2.2  2000/02/08  19:29:27  eric
 * Calcul sur les tenseurs.
 *
 * Revision 2.1  2000/02/04  17:14:27  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/02/04  16:37:51  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "etoile.h"
#include "graphique.h"

namespace Lorene {
void Etoile_bin::kinematics(double omega, double x_axe) {

    int nz = mp.get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 
    
    // --------------------
    // Computation of B^i/N
    // --------------------
    
    //  1/ Computation of  - omega m^i

    const Coord& xa = mp.xa ; 
    const Coord& ya = mp.ya ; 

    bsn.set_etat_qcq() ; 
    
    bsn.set(0) = omega * ya ;
    bsn.set(1) = - omega * (xa - x_axe) ;
    bsn.set(2) = 0 ;
    
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
}
