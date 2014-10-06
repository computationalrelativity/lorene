/*
 * Method of the class Map_et for the resolution of the 2-D Poisson
 *  equation.
 *
 * (see file map.h for documentation).
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


char map_et_poisson2d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:13:13  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2002/02/07 14:55:58  e_gourgoulhon
 * Corrected a bug when the source is known only in the coefficient
 * space.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/11/07  14:21:03  eric
 * Correction d'une erreur dans le cas T_SIN_I (calcul de R(u)).
 *
 * Revision 2.3  2000/10/26  15:58:00  eric
 * Correction cas T_COS_P : l'import de saff_q se fait par copie du Tbl.
 *
 * Revision 2.2  2000/10/12  15:37:43  eric
 * Traitement des bases spectrales dans le cas T_COS_P.
 *
 * Revision 2.1  2000/10/11  15:15:43  eric
 * 1ere version operationnelle.
 *
 * Revision 2.0  2000/10/09  13:47:17  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene:
#include "map.h"
#include "cmp.h"
#include "param.h"

//*****************************************************************************

void Map_et::poisson2d(const Cmp& source_mat, const Cmp& source_quad, 
		       Param& par, Cmp& uu) const {
    
    assert(source_mat.get_etat() != ETATNONDEF) ; 
    assert(source_quad.get_etat() != ETATNONDEF) ; 
    assert(source_mat.get_mp()->get_mg() == mg) ; 
    assert(source_quad.get_mp()->get_mg() == mg) ; 
    assert(uu.get_mp()->get_mg() == mg) ; 

    assert( source_quad.check_dzpuis(4) ) ; 
    
    double& lambda = par.get_double_mod(0) ; 
    int nz = mg->get_nzone() ; 
    int nzm1 = nz-1 ; 

    // Special case of a vanishing source 
    // ----------------------------------

    if (    (source_mat.get_etat() == ETATZERO) 
	 && (source_quad.get_etat() == ETATZERO) ) {
	
	uu = 0 ; 
	lambda = 1 ; 
	return ; 
    }

    int base_t = ((source_mat.va).base).get_base_t(0) ; 

    switch (base_t) {
    
    //==================================================================
    //		case T_COS_P
    //==================================================================

	case T_COS_P : {

    // Construction of a Map_af which coincides with *this on the equator
    // ------------------------------------------------------------------
    
    double theta0 = M_PI / 2 ;	    // Equator
    double phi0 = 0 ; 

    Map_af mpaff(*this) ; 

    for (int l=0 ; l<nz ; l++) {
	double rmax = val_r(l, double(1), theta0, phi0) ;
	switch ( mg->get_type_r(l) ) {
	    case RARE:	{
		double rmin = val_r(l, double(0), theta0, phi0) ;
		mpaff.set_alpha(rmax - rmin, l) ;
		mpaff.set_beta(rmin, l) ;
		break ; 
	    }
	    
	    case FIN:	{
		double rmin = val_r(l, double(-1), theta0, phi0) ;
		mpaff.set_alpha( double(.5) * (rmax - rmin), l ) ;
		mpaff.set_beta( double(.5) * (rmax + rmin), l) ;
		break ;
	    }
	    
	    case UNSURR: {
		double rmin = val_r(l, double(-1), theta0, phi0) ;
		double umax = double(1) / rmin ;
		double umin = double(1) / rmax ;
		mpaff.set_alpha( double(.5) * (umin - umax),  l) ; 
		mpaff.set_beta( double(.5) * (umin + umax), l) ; 
		break ;
	    }
	    
	    default: {
		cout << "Map_et::poisson2d: unknown type_r ! " << endl ;
		abort () ;
		break ;
	    }
	    
	}
    }

    // Importation of source_mat and source_quad of the affine mapping
    // ---------------------------------------------------------------
    Cmp saff_m(mpaff) ; 
    saff_m.import( nzm1, source_mat ) ; 
    (saff_m.va).set_base( (source_mat.va).base ) ; 
    
    Cmp saff_q(mpaff) ; 
    
    // In order to use Cmp::import with dzpuis != 0 :
    Cmp set_q = source_quad ; 
    set_q.set_dzpuis(0) ;	// dzpuis artificially set to 0

    saff_q.import( nzm1, set_q ) ; 
    (saff_q.va).set_base( (set_q.va).base ) ; 

    // Copy in the external domain :
    if ( (set_q.va).get_etat() == ETATQCQ) {
        (set_q.va).coef_i() ; // the values in configuration space are required
        assert(   (set_q.va).c->get_etat() == ETATQCQ ) ;
        assert(  (saff_q.va).c->get_etat() == ETATQCQ ) ;
	*( (saff_q.va).c->t[nzm1] ) = *( (set_q.va).c->t[nzm1] ) ;
    }

    // the true dzpuis is restored : 
    saff_q.set_dzpuis( source_quad.get_dzpuis() ) ; 

    // Resolution of the 2-D Poisson equation on the spherical domains
    // ---------------------------------------------------------------
    
    Cmp uaff(mpaff) ; 
    
    mpaff.poisson2d(saff_m, saff_q, par, uaff) ;
    
    // Importation of the solution on the Map_et mapping *this
    // -------------------------------------------------------
    
    uu.import(uaff) ; 
    
    uu.va.set_base( uaff.va.base ) ;	// same spectral bases 
    
    break ;
    }
    
    //==================================================================
    //		case T_SIN_I
    //==================================================================

	case T_SIN_I : {
	     
    //-------------------------------
    // Computation of the prefactor a  ---> Cmp apre
    //-------------------------------

    Mtbl unjj = 1 + srdrdt*srdrdt  ;

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
    double lambda_relax = par.get_double() ; 
    double unmlambda_relax = 1. - lambda_relax ; 
    double precis = par.get_double(1) ;     
    
    Cmp& ssj = par.get_cmp_mod() ; 
    
    Cmp ssjm1 = ssj ; 
    Cmp ssjm2 = ssjm1 ; 

    Cmp ssj_q(*this) ;
    ssj_q = 0 ;  
    
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

    cout << "Map_et::poisson2d : relat. diff. u^J <-> u^{J-1} : " << endl ;
    
//==========================================================================
//==========================================================================
//			    Start of iteration 
//==========================================================================
//==========================================================================

    Tbl tdiff(nz) ; 
    double diff ; 
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

    //-------------------
    // 1/x d^2uu/dtheta^2
    //-------------------
    
    Valeur sxlapang = uu.va ; 

    sxlapang = sxlapang.d2sdt2() ;    
    
    sxlapang = sxlapang.sx() ;	  //  d^2(uu)/dth^2 /x      in the nucleus
				  //  d^2(uu)/dth^2         in the shells
				  //  d^2(uu)/dth^2 /(x-1)  in the ZEC
    
    //---------------------------------------------------------------
    //  Computation of 
    // [ (dR/dx)^{-1} ( A - 1 ) duu/dx + 1/R (B - 1) d^2uu/dth^2 ] / x
    //
    // with A = 1/(dRdx) R/(x+beta_l/alpha_l) unjj 
    //	    B = [1/(dRdx) R/(x+beta_l/alpha_l)]^2 unjj 
    // 
    //  The result is put in uu (via vuu)
    //---------------------------------------------------------------

    // Intermediate quantity jac which value is
    //     (dR/dx)^{-1}   in the nucleus and the shells
    //     +(dU/dx)^{-1}  in the ZEC

    Mtbl jac = dxdr ; 
    if (mg->get_type_r(nzm1) == UNSURR) {
	jac.annule(nzm1, nzm1) ; 
	Mtbl jac_ext = dxdr ; 
	jac_ext.annule(0, nzm1-1) ;
	jac_ext = - jac_ext ; 
	jac = jac + jac_ext ; 
    }

    uu.set_etat_qcq() ; 
    
    Base_val sauve_base = duudx.base ; 
    
    vuu = jac * ( rsxdxdr * unjj - 1.) * duudx 
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
		+ 2. * dxdr * sr2drdt * d2uudtdx  ;

    vuu += dxdr * ( sr2d2rdt2 + dxdr * ( 
		dxdr* unjj * d2rdx2 
		- 2. * sr2drdt * d2rdtdx  ) 
				 ) * duudx ;		    

    vuu.set_base(sauve_base) ; 

    // Since the assignment is performed on vuu (uu.va), the treatment
    //  of uu.dzpuis must be performed by hand:
    
    uu.set_dzpuis(4) ; 

    //====================================================================
    //   Computation of the effective source s^J of the ``affine''
    //     Poisson equation 
    //====================================================================
    
    ssj = lambda_relax * ssjm1 + unmlambda_relax * ssjm2 ; 
    
    ssj = ( source_mat + source_quad + uu + (amax - apre) * ssj ) / amax ; 

    (ssj.va).set_base((source_mat.va).base) ; 
    
    //====================================================================
    //   Resolution of the ``affine'' Poisson equation 
    //====================================================================
    
    assert( uu.check_dzpuis( ssj.get_dzpuis() ) ) ; 
    					               
    mpaff.poisson2d(ssj, ssj_q, par, uu) ; 

    tdiff = diffrel(vuu, vuujm1) ; 

    diff = max(tdiff) ;
    

    cout << "  step " << niter << " :  " ; 
    for (int l=0; l<nz; l++) {
	cout << tdiff(l) << "  " ;  
    }
    cout << endl ; 

    //=================================
    //  Updates for the next iteration
    //=================================
    
    ssjm2 = ssjm1 ; 
    ssjm1 = ssj ; 
    vuujm1 = vuu ; 
    
    niter++ ; 
    
    }	// End of iteration 
    while ( (diff > precis) && (niter < nitermax) ) ;
    
//==========================================================================
//==========================================================================
//			    End of iteration 
//==========================================================================
//==========================================================================
	

    break ; 
    }
    
	default : {
	    cout << "Map_et::poisson2d : unkown theta basis !" << endl ; 
	    cout << "  basis : " << hex << base_t << endl ; 
	    abort() ; 
	    break ; 
	}
    }	// End of switch on base_t
}


