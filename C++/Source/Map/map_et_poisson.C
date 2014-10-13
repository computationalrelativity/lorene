/*
 * Method of the class Map_et for the (iterative) resolution of the scalar
 *  Poisson equation.
 *
 * (see file map.h for the documentation).
 *
 */

/*
 *   Copyright (c) 2004 Francois Limousin
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
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


char map_et_poisson_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.7  2014/10/13 08:53:05  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2005/08/25 12:14:09  p_grandclement
 * Addition of a new method to solve the scalar Poisson equation, based on a multi-domain Tau-method
 *
 * Revision 1.5  2005/04/04 21:31:31  e_gourgoulhon
 *  Added argument lambda to method poisson_angu
 *  to deal with the generalized angular Poisson equation:
 *     Lap_ang u + lambda u = source.
 *
 * Revision 1.4  2004/06/22 12:20:17  j_novak
 * *** empty log message ***
 *
 * Revision 1.3  2004/05/25 14:28:01  f_limousin
 * First version of method Map_et::poisson_angu().
 *
 * Revision 1.2  2003/10/15 21:11:26  e_gourgoulhon
 * Added method poisson_angu.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.7  2000/05/22  14:55:30  phil
 * ajout du cas dzpuis = 3
 *
 * Revision 1.6  2000/03/30  09:18:37  eric
 * Modifs affichage.
 *
 * Revision 1.5  2000/03/29  12:01:38  eric
 * *** empty log message ***
 *
 * Revision 1.4  2000/03/29  11:48:09  eric
 * Modifs affichage.
 *
 * Revision 1.3  2000/03/10  15:48:25  eric
 * MODIFS IMPORTANTES:
 *   ssj est desormais traitee comme un Cmp (et non plus une Valeur) ce qui
 *     permet un traitement automatique du dzpuis associe.
 *   Traitement de dzpuis.
 *
 * Revision 1.2  2000/03/07  16:50:57  eric
 * Possibilite d'avoir une source avec dzpuis = 2.
 *
 * Revision 1.1  1999/12/22  17:11:24  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "map.h"
#include "cmp.h"
#include "scalar.h"
#include "param.h"
#include "graphique.h"

//*****************************************************************************

namespace Lorene {

void Map_et::poisson(const Cmp& source, Param& par, Cmp& uu) const {

    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp() == this) ;

    assert( source.check_dzpuis(2) || source.check_dzpuis(4) 
	    || source.check_dzpuis(3)) ; 
     
    assert(uu.get_mp() == this) ; 
    assert(uu.check_dzpuis(0)) ; 

    int nz = mg->get_nzone() ; 
    int nzm1 = nz - 1 ; 
    
    // Indicator of existence of a compactified external domain
    bool zec = false ; 		
    if (mg->get_type_r(nzm1) == UNSURR) {
	zec = true ;
    }
     
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

    const Valeur& std2uudpdx = duudx.stdsdp() ;    // 1/sin(theta) d^2/dxdphi   

    //------------------
    // Angular Laplacian 
    //------------------
    
    Valeur sxlapang = uu.va ; 

    sxlapang.ylm() ; 
    
    sxlapang = sxlapang.lapang() ;    
    
    sxlapang = sxlapang.sx() ;	  //  Lap_ang(uu) /x      in the nucleus
				  //  Lap_ang(uu)         in the shells
				  //  Lap_ang(uu) /(x-1)  in the ZEC
    
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

    if (zec) {
	varduudx.annule(nzm1) ;	    // term in d/dx set to zero in the ZEC
    }

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

    // Since the assignment is performed on vuu (uu.va), the treatment
    //  of uu.dzpuis must be performed by hand:
    
    uu.set_dzpuis(4) ; 

    if (source.get_dzpuis() == 2) {
	uu.dec2_dzpuis() ;		    // uu.dzpuis: 4 -> 2
    }
    
    if (source.get_dzpuis() == 3) {
	uu.dec_dzpuis() ;		//uu.dzpuis 4 -> 3
    }
    
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
    
    if ( source.get_dzpuis() == 0 ){
	ssj.set_dzpuis( 4 ) ;
    }
    else {
	ssj.set_dzpuis( source.get_dzpuis() ) ;    // Choice of the resolution 
					       //  dzpuis = 2, 3 or 4
    }

    assert( uu.check_dzpuis( ssj.get_dzpuis() ) ) ; 
    		
    mpaff.poisson(ssj, par_nul, uu) ; 

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

}



//*****************************************************************************
// VERSION WITH A TAU METHOD
//*****************************************************************************

void Map_et::poisson_tau(const Cmp& source, Param& par, Cmp& uu) const {

    assert(source.get_etat() != ETATNONDEF) ; 
    assert(source.get_mp() == this) ;

    assert( source.check_dzpuis(2) || source.check_dzpuis(4) 
	    || source.check_dzpuis(3)) ; 
     
    assert(uu.get_mp() == this) ; 
    assert(uu.check_dzpuis(0)) ; 

    int nz = mg->get_nzone() ; 
    int nzm1 = nz - 1 ; 
    
    // Indicator of existence of a compactified external domain
    bool zec = false ; 		
    if (mg->get_type_r(nzm1) == UNSURR) {
	zec = true ;
    }
     
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

    cout << "Map_et::poisson_tau : relat. diff. u^J <-> u^{J-1} : " << endl ;
    
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

    const Valeur& std2uudpdx = duudx.stdsdp() ;    // 1/sin(theta) d^2/dxdphi   

    //------------------
    // Angular Laplacian 
    //------------------
    
    Valeur sxlapang = uu.va ; 

    sxlapang.ylm() ; 
    
    sxlapang = sxlapang.lapang() ;    
    
    sxlapang = sxlapang.sx() ;	  //  Lap_ang(uu) /x      in the nucleus
				  //  Lap_ang(uu)         in the shells
				  //  Lap_ang(uu) /(x-1)  in the ZEC
    
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

    if (zec) {
	varduudx.annule(nzm1) ;	    // term in d/dx set to zero in the ZEC
    }

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

    // Since the assignment is performed on vuu (uu.va), the treatment
    //  of uu.dzpuis must be performed by hand:
    
    uu.set_dzpuis(4) ; 

    if (source.get_dzpuis() == 2) {
	uu.dec2_dzpuis() ;		    // uu.dzpuis: 4 -> 2
    }
    
    if (source.get_dzpuis() == 3) {
	uu.dec_dzpuis() ;		//uu.dzpuis 4 -> 3
    }
    
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
    
    if ( source.get_dzpuis() == 0 ){
	ssj.set_dzpuis( 4 ) ;
    }
    else {
	ssj.set_dzpuis( source.get_dzpuis() ) ;    // Choice of the resolution 
					       //  dzpuis = 2, 3 or 4
    }

    assert( uu.check_dzpuis( ssj.get_dzpuis() ) ) ; 
    		
    mpaff.poisson_tau(ssj, par_nul, uu) ; 

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
}

void Map_et::poisson_angu(const Scalar& source, Param& par, Scalar& uu,
    double lambda) const {
    
    if (lambda != double(0)) {
        cout << 
        "Map_et::poisson_angu : the case lambda != 0 is not treated yet !"
        << endl ; 
        abort() ; 
    }

    assert(source.get_mp() == *this) ;
    assert(uu.get_mp() == *this) ; 

    int nz = mg->get_nzone() ; 
    int nzm1 = nz - 1 ; 
 
    int* nrm6 = new int[nz];
    for (int l=0; l<=nzm1; l++) 
	nrm6[l] = mg->get_nr(l) - 6 ; 
 
//##     // Indicator of existence of a compactified external domain
//     bool zec = false ; 		
//     if (mg->get_type_r(nzm1) == UNSURR) {
// 	zec = true ;
//     }

    //-------------------
    //  Initializations 
    //-------------------
    
    int nitermax = par.get_int() ; 
    int& niter = par.get_int_mod() ; 
    double relax = par.get_double() ; 
    double precis = par.get_double(1) ;     
    
    Cmp& ssjcmp = par.get_cmp_mod() ; 
    
    Scalar ssj (ssjcmp) ;
    Scalar ssjm1 (ssj) ;
 
    int dzpuis = source.get_dzpuis() ;
    ssj.set_dzpuis(dzpuis) ;
    uu.set_dzpuis(dzpuis) ;
    ssjm1.set_dzpuis(dzpuis) ;

    Valeur& vuu = uu.set_spectral_va() ;

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

    cout << "Map_et::poisson angu : relat. diff. u^J <-> u^{J-1} : " << endl ;
    
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
    
    Valeur duudx = (uu.set_spectral_va()).dsdx() ;	    // d/dx 

    const Valeur& d2uudxdx = duudx.dsdx() ;	    // d^2/dxdx


    const Valeur& d2uudtdx = duudx.dsdt() ;	    // d^2/dxdtheta

    const Valeur& std2uudpdx = duudx.stdsdp() ;    // 1/sin(theta) d^2/dxdphi  

    //---------------------------------------
    // Computation of  R(u)  
    //
    //  The result is put in uu (via vuu)
    //---------------------------------------

    Mtbl unjj = srdrdt*srdrdt + srstdrdp*srstdrdp ;

    Base_val sauve_base = vuu.base ; 

    vuu =  - d2uudxdx * dxdr * dxdr * unjj
		+ 2. * dxdr * ( sr2drdt * d2uudtdx 
			      + sr2stdrdp * std2uudpdx ) ;

    vuu.set_base(sauve_base) ; 

    vuu += dxdr * ( lapr_tp + dxdr * ( 
		dxdr * unjj * d2rdx2 
		- 2. * ( sr2drdt * d2rdtdx  + sr2stdrdp * sstd2rdpdx ) ) 
				 ) * duudx ;		    

    vuu.set_base(sauve_base) ; 

    uu.mult_r() ;
    uu.mult_r() ;

    //====================================================================
    //   Computation of the effective source s^J of the ``affine''
    //     Poisson equation 
    //====================================================================

    uu.filtre_r(nrm6) ;
//    uu.filtre_phi(1) ;
//    uu.filtre_theta(1) ;
 
    ssj = source + uu ; 

    ssj = (1-relax) * ssj + relax * ssjm1 ; 
 
    (ssj.set_spectral_va()).set_base((source.get_spectral_va()).base) ; 


    //====================================================================
    //   Resolution of the ``affine'' Poisson equation 
    //====================================================================

    mpaff.poisson_angu(ssj, par_nul, uu) ; 
    
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
    
    vuujm1 = vuu ; 
    ssjm1 = ssj ;
 
    niter++ ; 
    
    }	// End of iteration 
    while ( (diff > precis) && (niter < nitermax) ) ;

//==========================================================================
//==========================================================================
//			    End of iteration 
//==========================================================================
//==========================================================================

    uu.set_dzpuis( source.get_dzpuis() ) ;  // dzpuis unchanged

}



}
