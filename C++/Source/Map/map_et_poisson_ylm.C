/*
 *  Method of the class Map_et for the (iterative) resolution of the scalar
 *   Poisson equation with a multipole falloff condition at the outer boundary
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.3  2016/12/05 16:17:58  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2014/10/13 08:53:05  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2004/12/30 15:56:42  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "map.h"
#include "cmp.h"
#include "param.h"

//*****************************************************************************

namespace Lorene {

void Map_et::poisson_ylm(const Cmp& source, Param& par, Cmp& uu, int nylm, double* intvec) const {

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
    lambda=1.0;
    double unmlambda = 1. - lambda ; 
    double precis = par.get_double(1) ;     

    cout <<"relax in map_et:"<<lambda<<" "<<unmlambda<<endl;
    
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
    
    // *****************************************************************

    mpaff.poisson_ylm(ssj, par_nul, uu, nylm, intvec) ; 

    // *****************************************************************

    tdiff = diffrel(vuu, vuujm1) ; 

    diff = max(tdiff) ;
    

        cout << "  iter: " << niter << " :  " ; 
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
