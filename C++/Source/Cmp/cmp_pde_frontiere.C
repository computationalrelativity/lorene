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
 * Revision 1.5  2004/11/23 12:49:58  f_limousin
 * Intoduce function poisson_dir_neu(...) to solve a scalar poisson
 * equation with a mixed boundary condition (Dirichlet + Neumann).
 *
 * Revision 1.4  2004/05/20 07:04:02  k_taniguchi
 * Recovery of "->get_angu()" in the assertion of Map_af::poisson_frontiere
 * because "limite" is the boundary value.
 *
 * Revision 1.3  2004/03/31 11:18:42  f_limousin
 * Methods Map_et::poisson_interne, Map_af::poisson_interne and
 * Cmp::poisson_neumann_interne have been implemented to solve the
 * continuity equation for strange stars.
 *
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
#include "scalar.h" 
#include "param.h" 
#include "cmp.h"

Mtbl_cf sol_poisson_frontiere(const Map_af&, const Mtbl_cf&, const Mtbl_cf&,
			      int, int, int, double = 0.,
			      double = 0.) ;

Mtbl_cf sol_poisson_frontiere_double (const Map_af&, const Mtbl_cf&, const Mtbl_cf&,
				    const Mtbl_cf&, int) ;

Mtbl_cf sol_poisson_interne(const Map_af&, const Mtbl_cf&, const Mtbl_cf&) ;

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

Cmp Cmp::poisson_neumann_interne (const Valeur& limite, 
				  Param& par, Cmp& resu) const {
    
    mp->poisson_interne (*this, limite, par, resu) ; 
    return resu ;    
}

Cmp Cmp::poisson_frontiere_double (const Valeur& lim_func, const Valeur& lim_der, 
				    int num_zone) const {
    Cmp resu(*mp) ;
    mp->poisson_frontiere_double (*this, lim_func, lim_der, num_zone, resu) ; 
    return resu ;    
}		

void Map_et::poisson_frontiere(const Cmp&, const Valeur&, int, int, Cmp&, double, double) const {
    cout << "Procedure non implantee ! " << endl ;
    abort() ;
}

void Map_et::poisson_frontiere_double (const Cmp&, const Valeur&, const Valeur&,
				int, Cmp&) const {
    cout << "Procedure non implantee ! " << endl ;
    abort() ;
}

void Map_af::poisson_frontiere(const Cmp& source, const Valeur& limite, 
			       int type_raccord, int num_front, 
			       Cmp& pot, double fact_dir, double fact_neu) const {
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp()->get_mg() == mg) ; 
    assert(pot.get_mp()->get_mg() == mg) ; 
    
    assert( source.check_dzpuis(2) || source.check_dzpuis(4) 
	    || source.check_dzpuis(3)) ; 
    assert ((type_raccord == 1) || (type_raccord==2)|| (type_raccord==3)) ;
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
    
    Mtbl_cf resu (sourva.get_mg(), sourva.base) ;

    if (type_raccord == 3){

    // Call to the Mtbl_cf version
    // ---------------------------
	resu = sol_poisson_frontiere(*this, so_cf, auxiliaire, 
					type_raccord, num_front, dzpuis,
					 fact_dir, fact_neu) ;
    }
    else{
	resu = sol_poisson_frontiere(*this, so_cf, auxiliaire, 
					     type_raccord, num_front, dzpuis) ;
    }

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



void Map_et::poisson_interne(const Cmp& source, const Valeur& limite,
			     Param& par, Cmp& uu) const {

    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp() == this) ;

    assert(uu.get_mp() == this) ; 

    int nz = mg->get_nzone() ; 
     
    //-------------------------------
    // Computation of the prefactor a  ---> Cmp apre
    //-------------------------------

    Mtbl unjj = 1 + srdrdt*srdrdt + srstdrdp*srstdrdp ;

    Mtbl apre1(*mg) ; 
    apre1.set_etat_qcq() ; 
    for (int l=0; l<nz; l++) {
	*(apre1.t[l]) = alpha[l]*alpha[l] ; 
    }

    apre1 = apre1 * dxdr * dxdr * unjj ;

    Cmp apre(*this) ; 
    apre = apre1 ; 
    
    Tbl amax0 = max(apre1) ;	// maximum values in each domain

    // The maximum values of a in each domain are put in a Mtbl
    Mtbl amax1(*mg) ; 
    amax1.set_etat_qcq() ; 
    for (int l=0; l<nz; l++) {
	*(amax1.t[l]) = amax0(l) ; 
    }
    
    Cmp amax(*this) ; 
    amax = amax1 ; 

    
    //-------------------
    //  Initializations 
    //-------------------
    
    int nitermax = par.get_int() ; 
    int& niter = par.get_int_mod() ; 
    double lambda = par.get_double() ; 
    double unmlambda = 1. - lambda ; 
    double precis = par.get_double(1) ;     
    
    Cmp& ssj = par.get_cmp_mod() ; 
    
    Cmp ssjm1 = ssj ; 
    Cmp ssjm2 = ssjm1 ; 

    Valeur& vuu = uu.va ; 

    Valeur vuujm1(*mg) ;
    if (uu.get_etat() == ETATZERO) {
	vuujm1 = 1 ;	// to take relative differences
	vuujm1.set_base( vuu.base ) ; 
    }
    else{
	vuujm1 = vuu ; 
    }
    
    // Affine mapping for the Laplacian-tilde

    Map_af mpaff(*this) ; 
    Param par_nul ; 

    cout << "Map_et::poisson : relat. diff. u^J <-> u^{J-1} : " << endl ;
    
//==========================================================================
//==========================================================================
//			    Start of iteration 
//==========================================================================
//==========================================================================

    Tbl tdiff(nz) ; 
    niter = 0 ; 
    
    do {

    //====================================================================
    //		Computation of R(u)    (the result is put in uu)
    //====================================================================


    //-----------------------
    // First operations on uu
    //-----------------------
    
    Valeur duudx = (uu.va).dsdx() ;	    // d/dx 

    const Valeur& d2uudtdx = duudx.dsdt() ;	    // d^2/dxdtheta

    const Valeur& std2uudpdx = duudx.stdsdp() ;    // 1/sin(theta) d^2/dxdphi   

    //------------------
    // Angular Laplacian 
    //------------------
    
    Valeur sxlapang = uu.va ; 

    sxlapang.ylm() ; 
    
    sxlapang = sxlapang.lapang() ;    
    
    sxlapang = sxlapang.sx() ;	  //  Lap_ang(uu) /x      in the nucleus
				    
    //---------------------------------------------------------------
    //  Computation of 
    // [ 2 /(dRdx) ( A - 1 ) duu/dx + 1/R (B - 1) Lap_ang(uu) ] / x
    //
    // with A = 1/(dRdx) R/(x+beta_l/alpha_l) unjj 
    //	    B = [1/(dRdx) R/(x+beta_l/alpha_l)]^2 unjj 
    // 
    //  The result is put in uu (via vuu)
    //---------------------------------------------------------------

    Valeur varduudx = duudx ; 

    uu.set_etat_qcq() ; 
    
    Base_val sauve_base = varduudx.base ; 
    
    vuu = 2. * dxdr * ( rsxdxdr * unjj - 1.) * varduudx 
		+ ( rsxdxdr*rsxdxdr * unjj - 1.) * xsr * sxlapang ; 

    vuu.set_base(sauve_base) ; 

    vuu = vuu.sx() ; 

    //---------------------------------------
    // Computation of  R(u)  
    //
    //  The result is put in uu (via vuu)
    //---------------------------------------


    sauve_base = vuu.base ; 

    vuu =  xsr * vuu 
		+ 2. * dxdr * (	sr2drdt * d2uudtdx 
			      + sr2stdrdp * std2uudpdx ) ;

    vuu += dxdr * ( lapr_tp + dxdr * ( 
		dxdr* unjj * d2rdx2 
		- 2. * ( sr2drdt * d2rdtdx  + sr2stdrdp * sstd2rdpdx ) ) 
				 ) * duudx ;		    

    vuu.set_base(sauve_base) ; 

    //====================================================================
    //   Computation of the effective source s^J of the ``affine''
    //     Poisson equation 
    //====================================================================
    
    ssj = lambda * ssjm1 + unmlambda * ssjm2 ; 
    
    ssj = ( source + uu + (amax - apre) * ssj ) / amax ; 

    (ssj.va).set_base((source.va).base) ; 
    
    //====================================================================
    //   Resolution of the ``affine'' Poisson equation 
    //====================================================================
    
    mpaff.poisson_interne(ssj, limite, par_nul, uu) ; 

    tdiff = diffrel(vuu, vuujm1) ; 

    cout << "  step " << niter << " :  " ; 
    cout << tdiff(0) << "  " ;  
 
    cout << endl ; 

    //=================================
    //  Updates for the next iteration
    //=================================
    
    ssjm2 = ssjm1 ; 
    ssjm1 = ssj ; 
    vuujm1 = vuu ; 
    
    niter++ ; 
    
    }	// End of iteration 
    while ( (tdiff(0) > precis) && (niter < nitermax) ) ;
    
//==========================================================================
//==========================================================================
//			    End of iteration 
//==========================================================================
//==========================================================================

}


void Map_af::poisson_interne(const Cmp& source, const Valeur& limite, 
			     Param& par, Cmp& pot) const {
    
    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp()->get_mg() == mg) ; 
    assert(pot.get_mp()->get_mg() == mg) ; 
    assert (source.get_mp()->get_mg()->get_angu() == limite.get_mg()) ;
 
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
    
    Mtbl_cf resu = sol_poisson_interne(*this, so_cf, auxiliaire) ;
    
    // Final result returned as a Cmp
    // ------------------------------
    
    pot.set_etat_zero() ;  // to call Cmp::del_t().
    pot.set_etat_qcq() ;  
    pot.va = resu ;
    (pot.va).ylm_i() ; // On repasse en base standard.	    
    pot.set_dzpuis(0) ; 
}

