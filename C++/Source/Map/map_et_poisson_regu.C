/*
 * Method of the class Map_et for the (iterative) resolution of the scalar
 *  Poisson equation by using regularized source.
 *
 * (see file map.h for the documentation).
 *
 */

/*
 *   Copyright (c) 2000-2001 Keisuke Taniguchi
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
 * Revision 1.3  2016/12/05 16:17:58  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2014/10/13 08:53:05  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.8  2000/09/27  14:07:14  keisuke
 * Traitement des bases spectrales de d_logn_auto_div.
 *
 * Revision 2.7  2000/09/26  15:41:20  keisuke
 * Correction erreur: la triade de duu_div doit etre celle de *this et
 *  non celle de l'objet temporaire mpaff.
 *
 * Revision 2.6  2000/09/25  15:03:34  keisuke
 * Correct the derivative duu_div.
 *
 * Revision 2.5  2000/09/11  14:00:20  keisuke
 * Suppress "uu = uu_regu + uu_div" because of double setting (in poisson_regular).
 *
 * Revision 2.4  2000/09/07  15:51:29  keisuke
 * Minor change.
 *
 * Revision 2.3  2000/09/07  15:30:07  keisuke
 * Add a new argument Cmp& uu.
 *
 * Revision 2.2  2000/09/04  15:56:15  keisuke
 * Change the argumant Cmp& duu_div_r into Tenseur& duu_div.
 *
 * Revision 2.1  2000/09/04  14:52:17  keisuke
 * Change the scheme of code into that of map_et_poisson.C.
 *
 * Revision 2.0  2000/09/01  09:55:33  keisuke
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "map.h"
#include "cmp.h"
#include "tenseur.h"
#include "param.h"

//*****************************************************************************

namespace Lorene {

void Map_et::poisson_regular(const Cmp& source, int k_div, int nzet,
			     double unsgam1, Param& par, Cmp& uu,
			     Cmp& uu_regu, Cmp& uu_div, Tenseur& duu_div,
			     Cmp& source_regu, Cmp& source_div) const {


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

    cout << "Map_et::poisson_regular : relat. diff. u^J <-> u^{J-1} : "
	 << endl ;

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


    //------------------------
    // First operations on uu
    //------------------------

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

    //------------------------------------------------------------------
    //  Computation of
    // [ 2 /(dRdx) ( A - 1 ) duu/dx + 1/R (B - 1) Lap_ang(uu) ] / x
    //
    // with A = 1/(dRdx) R/(x+beta_l/alpha_l) unjj
    //	    B = [1/(dRdx) R/(x+beta_l/alpha_l)]^2 unjj
    //
    //  The result is put in uu (via vuu)
    //------------------------------------------------------------------

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

    //----------------------------------------
    // Computation of R(u)
    //
    //  The result is put in uu (via vuu)
    //----------------------------------------

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
	ssj.set_dzpuis( source.get_dzpuis() ) ;
	                                       // Choice of the resolution
					       //  dzpuis = 2, 3 or 4
    }

    assert( uu.check_dzpuis( ssj.get_dzpuis() ) ) ;

    mpaff.poisson_regular(ssj, k_div, nzet, unsgam1, par_nul, uu,
                          uu_regu, uu_div, duu_div,
                          source_regu, source_div) ;
			  
    //======================================
    //  Gradient of the diverging part (from that computed on the Map_af)
    //======================================

    Valeur& dr_uu_div = duu_div.set(0).va ; 
    Valeur& dt_uu_div = duu_div.set(1).va ; 
    Valeur& dp_uu_div = duu_div.set(2).va ; 

    Base_val bv = dr_uu_div.base ; 
    dr_uu_div = alpha[0] * dr_uu_div * dxdr ;    
    dr_uu_div.set_base( bv ) ; 
    
    bv = dt_uu_div.base ; 
    dt_uu_div = alpha[0] * dt_uu_div * xsr - srdrdt * dr_uu_div ;
    dt_uu_div.set_base( bv ) ; 

    bv = dp_uu_div.base ; 
    dp_uu_div = alpha[0] * dp_uu_div * xsr - srstdrdp * dr_uu_div ;
    dp_uu_div.set_base( bv ) ; 
    
    duu_div.set_triad( this->get_bvect_spher() ) ;
    
    
    //========================================
    // Relative difference with previous step
    //========================================

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
}
