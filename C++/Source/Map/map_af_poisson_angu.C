/*
 *  Resolution of the angular Poisson equation. 
 *
 * (see file map.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003-2005 Eric Gourgoulhon & Jerome Novak
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

char map_af_poisson_angu_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2005/04/04 21:31:31  e_gourgoulhon
 *  Added argument lambda to method poisson_angu
 *  to deal with the generalized angular Poisson equation:
 *     Lap_ang u + lambda u = source.
 *
 * Revision 1.2  2003/10/16 08:49:23  j_novak
 * Added a flag to decide wether the output is in the Ylm or in the standard base.
 *
 * Revision 1.1  2003/10/15 21:11:26  e_gourgoulhon
 * Added method poisson_angu.
 *
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "tensor.h"
#include "param.h"

void Map_af::poisson_angu(const Scalar& source, Param& par, Scalar& uu,
             double lambda) const {

    assert(source.get_etat() != ETATNONDEF) ; 
	
	assert(&(source.get_mp()) == this ) ;
	assert(&(uu.get_mp()) == this ) ;

    // Spherical harmonic expansion of the source
    // ------------------------------------------
    
    const Valeur& sourva = source.get_spectral_va() ; 

    if (sourva.get_etat() == ETATZERO) {
		uu.set_etat_zero() ;
		return ;  
    }

    // Spectral coefficients of the source
    assert(sourva.get_etat() == ETATQCQ) ; 
    sourva.coef() ; 
    
    Valeur resu(mg) ; 
    resu = *(sourva.c_cf) ;	// copy of the coefficients of the source
    
    resu.ylm() ;			// spherical harmonic transform 
        
    // Call to the Mtbl_cf version
    // ---------------------------
    (resu.c_cf)->poisson_angu(lambda) ; 
	
    if (par.get_n_int() == 0) resu.ylm_i() ; // Back to standard bases 
                                             //in the case of no flag present
                                             // in the Param

    // Final result returned as a Scalar
    // ---------------------------------
    
    uu = resu ;
	
    uu.set_dzpuis( source.get_dzpuis() ) ;  // dzpuis unchanged
}

