/*
 * Method of the class Map_af for the resolution of the scalar 2D-Poisson
 *  equation
 */

/*
 *   Copyright (c) 2004 Philippe Grandclement
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
 * Revision 1.3  2016/12/05 16:17:56  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2014/10/13 08:53:02  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2004/08/24 09:14:42  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * 
 * $Header$
 *
 */

// Header Lorene:
#include "map.h"
#include "scalar.h"
#include "param_elliptic.h"

//*****************************************************************************

namespace Lorene {

void Map_af::sol_elliptic_2d(Param_elliptic& ope_var, 
			     const Scalar& source, Scalar& pot) const {
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp().get_mg() == mg) ; 
    assert(pot.get_mp().get_mg() == mg) ; 

    assert( source.check_dzpuis(2) || source.check_dzpuis(4) 
	    || source.check_dzpuis(3)) ; 

    // Spherical harmonic expansion of the source
    // ------------------------------------------
    
    const Valeur& sourva = source.get_spectral_va() ; 

    if (sourva.get_etat() == ETATZERO) {
	pot.set_etat_zero() ;
	return ;  
    }

    // Spectral coefficients of the source
    assert(sourva.get_etat() == ETATQCQ) ; 
    
    Valeur rho(sourva.get_mg()) ; 
    sourva.coef() ; 
    rho = *(sourva.c_cf) ;	// copy of the coefficients of the source
      
    // On calcule les coefs radiaux de F et G 
    ope_var.var_F.set_spectral_va().coef() ;
    ope_var.var_G.set_spectral_va().coef() ;

    // Call to the Mtbl_cf version
    // ---------------------------
    Mtbl_cf resu = elliptic_solver (ope_var, *(rho.c_cf)) ;

    pot.set_etat_zero() ; 
    pot.set_etat_qcq() ; 
    pot.set_spectral_va() = resu ;
    pot.set_dzpuis(0) ;
}


}
