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


char map_af_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
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

Mtbl_cf sol_poisson(const Map_af&, const Mtbl_cf&, int) ;
//*****************************************************************************

void Map_af::poisson(const Cmp& source, Param&, Cmp& pot) const {
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp()->get_mg() == mg) ; 
    assert(pot.get_mp()->get_mg() == mg) ; 

    assert( source.check_dzpuis(2) || source.check_dzpuis(4) 
	    || source.check_dzpuis(3)) ; 
    
    int dzpuis ; 
    
    if (source.dz_nonzero()){
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
    Mtbl_cf resu = sol_poisson(*this, *(rho.c_cf), dzpuis) ;
    
    // Final result returned as a Cmp
    // ------------------------------
    
    pot.set_etat_zero() ;  // to call Cmp::del_t().

    pot.set_etat_qcq() ; 
    
    pot.va = resu ;
    (pot.va).ylm_i() ; // On repasse en base standard.	    

    pot.set_dzpuis(0) ; 
    
}


