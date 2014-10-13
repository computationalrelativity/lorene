/*
 *  Method of the class Map_af for the resolution of the scalar Poisson
 *   equation with a multipole falloff condition at the outer boundary
 *
 *    (see file map.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Joshua A. Faber
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

char map_af_poisson_ylm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:03  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2005/02/18 13:14:08  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
 * Revision 1.1  2004/12/30 15:55:57  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "map.h"
#include "cmp.h"

namespace Lorene {
Mtbl_cf sol_poisson_ylm(const Map_af&, const Mtbl_cf&, const int, const double*) ;
//*****************************************************************************

void Map_af::poisson_ylm(const Cmp& source, Param& , Cmp& pot, int nylm, double* intvec) const {
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp()->get_mg() == mg) ; 
    assert(pot.get_mp()->get_mg() == mg) ; 

    // Spherical harmonic expansion of the source
    // ------------------------------------------
    
    const Valeur& sourva = source.va ; 

    if (sourva.get_etat() == ETATZERO) {
	pot.set_etat_zero() ;
	return ;  
    }

    // Spectral coefficients of the source
    assert(sourva.get_etat() == ETATQCQ) ; 
    
    Valeur rho(sourva.get_mg()) ; 
    sourva.coef() ; 
    rho = *(sourva.c_cf) ;	// copy of the coefficients of the source
    
    rho.ylm() ;			// spherical harmonic transforms 
        
    // Call to the Mtbl_cf version
    // ---------------------------
    Mtbl_cf resu = sol_poisson_ylm(*this, *(rho.c_cf), nylm, intvec) ;
    
    // Final result returned as a Cmp
    // ------------------------------
    
    pot.set_etat_zero() ;  // to call Cmp::del_t().

    pot.set_etat_qcq() ; 
    
    pot.va = resu ;
    (pot.va).ylm_i() ; // On repasse en base standard.	    

}

}
