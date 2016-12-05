/*
 * Method of the class Map_af for the resolution of the scalar Poisson
 *  equation
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
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
 * Revision 1.7  2016/12/05 16:17:57  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:53:02  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2005/08/25 12:14:09  p_grandclement
 * Addition of a new method to solve the scalar Poisson equation, based on a multi-domain Tau-method
 *
 * Revision 1.4  2004/05/06 15:25:39  e_gourgoulhon
 * The case dzpuis=5 with null value in CED is well treated now.
 *
 * Revision 1.3  2004/02/20 10:55:23  j_novak
 * The versions dzpuis 5 -> 3 has been improved and polished. Should be
 * operational now...
 *
 * Revision 1.2  2004/02/06 10:53:52  j_novak
 * New dzpuis = 5 -> dzpuis = 3 case (not ready yet).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.9  2000/05/22  13:46:48  phil
 * ajout du cas dzpuis = 3
 *
 * Revision 1.8  2000/02/09  14:44:24  eric
 * Traitement de dzpuis ameliore.
 *
 * Revision 1.7  1999/12/22  16:37:10  eric
 * Ajout de pot.set_dzpuis(0) a la fin.
 *
 * Revision 1.6  1999/12/22  15:11:03  eric
 * Remplacement du test source.get_mp() == this  par
 *  source.get_mp()->get_mg() == mg
 * (idem pour pot),
 * afin de permettre l'appel par Map_et::poisson.
 *
 * Revision 1.5  1999/12/21  13:02:37  eric
 * Changement de prototype de la routine poisson : la solution est
 *  desormais passee en argument (et non plus en valeur de retour)
 *  pour s'adapter au prototype general de la fonction virtuelle
 *   Map::poisson.
 *
 * Revision 1.4  1999/12/21  10:06:29  eric
 * Ajout de l'argument (muet) Param&.
 *
 * Revision 1.3  1999/12/07  16:48:50  phil
 * On fait ylm_i avant de quitter
 *
 * Revision 1.2  1999/12/02  16:12:22  eric
 * *** empty log message ***
 *
 * Revision 1.1  1999/12/02  14:30:07  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "map.h"
#include "cmp.h"

namespace Lorene {
Mtbl_cf sol_poisson(const Map_af&, const Mtbl_cf&, int, bool match = true) ;
Mtbl_cf sol_poisson_tau(const Map_af&, const Mtbl_cf&, int) ;
//*****************************************************************************

void Map_af::poisson(const Cmp& source, Param& , Cmp& pot) const {
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp()->get_mg() == mg) ; 
    assert(pot.get_mp()->get_mg() == mg) ; 

    assert( source.check_dzpuis(2) || source.check_dzpuis(4) 
	    || source.check_dzpuis(3) || source.check_dzpuis(5) ) ; 
    
    bool match = true ;

    int dzpuis ; 
    
    if ( (source.dz_nonzero()) || (source.get_dzpuis() > 3)) { //##awkward??
	dzpuis = source.get_dzpuis() ; 
    }
    else{
	dzpuis = 4 ; 
    }

    match = !(dzpuis == 5) ;

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
    Mtbl_cf resu = sol_poisson(*this, *(rho.c_cf), dzpuis, match) ;
    
    // Final result returned as a Cmp
    // ------------------------------
    
    pot.set_etat_zero() ;  // to call Cmp::del_t().

    pot.set_etat_qcq() ; 
    
    pot.va = resu ;
    (pot.va).ylm_i() ; // On repasse en base standard.	    

    (dzpuis == 5) ? pot.set_dzpuis(3) : pot.set_dzpuis(0) ; 
    
}


		//----------------------
		// Tau version method
		//---------------------


void Map_af::poisson_tau(const Cmp& source, Param& , Cmp& pot) const {
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp()->get_mg() == mg) ; 
    assert(pot.get_mp()->get_mg() == mg) ; 

    assert( source.check_dzpuis(2) || source.check_dzpuis(4) 
	    || source.check_dzpuis(3)) ; 

    int dzpuis ; 
    
    if ( (source.dz_nonzero()) || (source.get_dzpuis() > 3)) { //##awkward??
	dzpuis = source.get_dzpuis() ; 
    }
    else{
	dzpuis = 4 ; 
    }

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
    
    Mtbl_cf resu = sol_poisson_tau(*this, *(rho.c_cf), dzpuis) ;
    
    // Final result returned as a Cmp
    // ------------------------------
    
    pot.set_etat_zero() ;  // to call Cmp::del_t().

    pot.set_etat_qcq() ; 
    
    pot.va = resu ;
    (pot.va).ylm_i() ; // On repasse en base standard.	    
    pot.set_dzpuis(0) ;
}

}
