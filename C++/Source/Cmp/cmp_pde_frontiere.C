/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


char cmp_pde_frontiere_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/10/03 15:58:45  j_novak
 * Cleaning of some headers
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2000/05/22  16:07:03  phil
 * *** empty log message ***
 *
 * Revision 2.5  2000/05/22  16:03:48  phil
 * ajout du cas dzpuis = 3
 *
 * Revision 2.4  2000/04/27  15:18:27  phil
 * ajout des procedures relatives a la resolution dans une seule zone avec deux conditions limites.
 *
 * Revision 2.3  2000/03/31  15:59:54  phil
 * gestion des cas ou la source est nulle.
 *
 * Revision 2.2  2000/03/20  13:08:53  phil
 * *** empty log message ***
 *
 * Revision 2.1  2000/03/17  17:33:05  phil
 * *** empty log message ***
 *
 * Revision 2.0  2000/03/17  17:25:08  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "cmp.h"

Mtbl_cf sol_poisson_frontiere(const Map_af&, const Mtbl_cf&, const Mtbl_cf&,
				    int, int, int) ;

Mtbl_cf sol_poisson_frontiere_double (const Map_af&, const Mtbl_cf&, const Mtbl_cf&,
				    const Mtbl_cf&, int) ;

Cmp Cmp::poisson_dirichlet (const Valeur& limite, int num_front) const {
    
    Cmp resu(*mp) ;
    mp->poisson_frontiere (*this, limite, 1, num_front, resu) ; 
    return resu ;          
}


Cmp Cmp::poisson_neumann (const Valeur& limite, int num_front) const {
    
    Cmp resu(*mp) ;
    mp->poisson_frontiere (*this, limite, 2, num_front, resu) ; 
    return resu ;    
}

Cmp Cmp::poisson_frontiere_double (const Valeur& lim_func, const Valeur& lim_der, 
				    int num_zone) const {
    Cmp resu(*mp) ;
    mp->poisson_frontiere_double (*this, lim_func, lim_der, num_zone, resu) ; 
    return resu ;    
}		

void Map_et::poisson_frontiere(const Cmp&, const Valeur&, int, int, Cmp&) const {
    cout << "Procedure non implantee ! " << endl ;
    abort() ;
}

void Map_et::poisson_frontiere_double (const Cmp&, const Valeur&, const Valeur&,
				int, Cmp&) const {
    cout << "Procedure non implantee ! " << endl ;
    abort() ;
}

void Map_af::poisson_frontiere(const Cmp& source, const Valeur& limite, int type_raccord, 
			int num_front, Cmp& pot) const {
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp()->get_mg() == mg) ; 
    assert(pot.get_mp()->get_mg() == mg) ; 
    assert (source.get_mp()->get_mg()->get_angu() == limite.get_mg()) ;
    
    assert( source.check_dzpuis(2) || source.check_dzpuis(4) 
	    || source.check_dzpuis(3)) ; 
    assert ((type_raccord == 1) || (type_raccord==2)) ;
    int dzpuis ; 
    
    if (source.dz_nonzero()){
	dzpuis = source.get_dzpuis() ; 
    }
    else{
	dzpuis = 4 ; 
    }

    // Spherical harmonic expansion of the source
    // ------------------------------------------
    
    Valeur sourva = source.va ;
    sourva.coef() ;
    sourva.ylm() ;
    
    // Pour gerer le cas ou source est dans ETAT_ZERO...
    Mtbl_cf so_cf (sourva.get_mg(), sourva.base) ;
    if (sourva.get_etat() == ETATZERO) {
	so_cf.set_etat_zero() ;
	}
    else
	so_cf = *sourva.c_cf ;
    
    
    Valeur conditions (limite) ;
    conditions.coef() ;
    conditions.ylm() ; // spherical harmonic transforms 
    
    // Pour gerer le cas ou condition est dans ETAT_ZERO...
    Mtbl_cf auxiliaire (conditions.get_mg(), conditions.base) ;
    if (conditions.get_etat() == ETATZERO)
	auxiliaire.set_etat_zero() ;
    else
	auxiliaire = *conditions.c_cf ;
	
    // Call to the Mtbl_cf version
    // ---------------------------
    
    Mtbl_cf resu = sol_poisson_frontiere(*this, so_cf, auxiliaire, 
					type_raccord, num_front, dzpuis) ;
    
    // Final result returned as a Cmp
    // ------------------------------
    
    pot.set_etat_zero() ;  // to call Cmp::del_t().
    pot.set_etat_qcq() ;  
    pot.va = resu ;
    (pot.va).ylm_i() ; // On repasse en base standard.	    
    pot.set_dzpuis(0) ; 
}


void Map_af::poisson_frontiere_double (const Cmp& source, const Valeur& lim_func, 
			    const Valeur& lim_der, int num_zone, Cmp& pot) const {
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp()->get_mg() == mg) ; 
    assert(pot.get_mp()->get_mg() == mg) ; 
    assert (source.get_mp()->get_mg()->get_angu() == lim_func.get_mg()) ;
    assert (source.get_mp()->get_mg()->get_angu() == lim_der.get_mg()) ;
    
    // Spherical harmonic expansion of the source
    // ------------------------------------------
    
    Valeur sourva = source.va ;
    sourva.coef() ;
    sourva.ylm() ;
    
    // Pour gerer le cas ou source est dans ETAT_ZERO...
    Mtbl_cf so_cf (sourva.get_mg(), sourva.base) ;
    if (sourva.get_etat() == ETATZERO) {
	so_cf.set_etat_zero() ;
	}
    else
	so_cf = *sourva.c_cf ;
    
    
    Valeur cond_func (lim_func) ;
    cond_func.coef() ;
    cond_func.ylm() ; // spherical harmonic transforms 
    
    // Pour gerer le cas ou condition est dans ETAT_ZERO...
    Mtbl_cf auxi_func (cond_func.get_mg(), cond_func.base) ;
    if (cond_func.get_etat() == ETATZERO)
	auxi_func.set_etat_zero() ;
    else
	auxi_func = *cond_func.c_cf ;
    
    Valeur cond_der (lim_der) ;
    cond_der.coef() ;
    cond_der.ylm() ; // spherical harmonic transforms 
    
    // Pour gerer le cas ou condition est dans ETAT_ZERO...
    Mtbl_cf auxi_der (cond_der.get_mg(), cond_der.base) ;
    if (cond_der.get_etat() == ETATZERO)
	auxi_der.set_etat_zero() ;
    else
	auxi_der = *cond_der.c_cf ;
    
    
    
    // Call to the Mtbl_cf version
    // ---------------------------
    
    Mtbl_cf resu = sol_poisson_frontiere_double (*this, so_cf, auxi_func,
				    auxi_der, num_zone) ;
    
    // Final result returned as a Cmp
    // ------------------------------
    
    pot.set_etat_zero() ;  // to call Cmp::del_t().
    pot.set_etat_qcq() ;  
    pot.va = resu ;
    (pot.va).ylm_i() ; // On repasse en base standard.
}
