/*
 *  Resolution of the angular Poisson equation. 
 *
 * (see file map.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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

void Map_af::poisson_angu(const Scalar& source, Param& , Scalar& uu) const {

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
    (resu.c_cf)->poisson_angu() ; 
	
	resu.ylm_i() ; // Back to standard bases 

    // Final result returned as a Scalar
    // ---------------------------------
    
    uu.set_etat_zero() ;  // to call Scalar::del_t().

    uu.set_etat_qcq() ; 
    
    uu.set_spectral_va() = resu ;
	
	uu.set_dzpuis( source.get_dzpuis() ) ;  // dzpuis unchanged
}

