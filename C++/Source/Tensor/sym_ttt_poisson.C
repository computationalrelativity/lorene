/*
 *  Resolution of the TT tensor Poisson equation
 *
 *    (see file sym_tensor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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

char sym_ttt_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:53:44  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2004/12/28 10:37:24  j_novak
 * Better way of enforcing zero divergence.
 *
 * Revision 1.3  2004/12/27 14:33:12  j_novak
 * New algorithm for the tensor Poisson eq.
 *
 * Revision 1.2  2004/03/04 09:50:41  e_gourgoulhon
 * Method poisson: use of new methods khi() and set_khi_mu.
 *
 * Revision 1.1  2003/11/07 16:53:52  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

// C headers
//#include <>

// Lorene headers
#include "tensor.h"
#include "param_elliptic.h"


namespace Lorene {
Sym_tensor_tt Sym_tensor_tt::poisson(int dzfin) const {

    // All this has a meaning only for spherical components...
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 
    //## ... and affine mapping, for the moment!
    assert(dynamic_cast<const Map_af*>(mp) != 0x0) ;
    assert( (dzfin == 0) || (dzfin == 2) ) ;
    Sym_tensor_tt resu(*mp, *triad, *met_div) ; 

    // Solution for the rr-component
    // ----------------------------

    const Scalar& source_rr = operator()(1,1) ;
    Scalar h_rr(*mp) ;
    int nz = mp->get_mg()->get_nzone() ;
   
    if (source_rr.get_etat() != ETATZERO) {

	//------------------------
	// The elliptic operator
	//------------------------
	  
	Param_elliptic param_hr(source_rr) ;
	for (int lz=0; lz<nz; lz++) 
	    param_hr.set_poisson_tens_rr(lz) ;
	  
	h_rr = source_rr.sol_elliptic(param_hr) ;
    }
    else
	h_rr.set_etat_zero() ;

    h_rr.inc_dzpuis(dzfin) ; //## can we improve here?
    resu.set(1,1) = h_rr ;

    // Solution for (eta / r) 
    //-----------------------
//     Scalar source_eta = - source_rr ;
//     source_eta.mult_r_dzpuis(3) ;
//     source_eta.mult_r_dzpuis(2) ;
//     h_rr.set_spectral_va().ylm() ;
//     Scalar tmp = 2*h_rr + h_rr.lapang() ;
//     if (dzfin == 0) 
// 	tmp.inc_dzpuis(2) ;
//     source_eta += tmp ;
//     source_eta = source_eta.primr() ;

//     source_eta.div_r_dzpuis(dzfin) ;

//     Scalar etasurr = (h_rr+source_eta).poisson_angu() ;

    Scalar source_eta = -resu(1,1).dsdr() ;
    source_eta.mult_r_dzpuis(dzfin) ;
    source_eta -= 3.*resu(1,1) ;
    Scalar etasurr = source_eta.poisson_angu() ;
		
    // Solution for mu
    // ---------------
	
    Scalar musurr = mu().poisson() ;
    musurr.div_r_dzpuis(dzfin) ;

    resu.set(1,1).set_spectral_va().ylm_i() ;
    
    Scalar** rcmp = resu.cmp ;

    Itbl idx(2) ;
    idx.set(0) = 1 ;	// r index
	
    // h^{r theta} : 
    // ------------
    idx.set(1) = 2 ;	// theta index
    *rcmp[position(idx)] = etasurr.dsdt() - musurr.stdsdp() ; 
    
    // h^{r phi} :
    // ------------
    idx.set(1) = 3 ;	// phi index
    *rcmp[position(idx)] = etasurr.stdsdp() + musurr.dsdt() ; 

    // h^{theta phi} and h^{phi phi}
    // -----------------------------
	
    //--------------  Computation of T^theta   --> taut : 
    
    Scalar tautst = resu(1,2).dsdr() ; 

    // dhrr contains  dh^{rt}/dr in all domains but the CED,    
    // in the CED:    r^2 dh^{rt}/dr        if dzfin = 0          (1)
    //                r^3 dh^{rt}/dr        if dzfin = 2          (2)
                                                    
    // Multiplication by r of dh^{rt}/dr (with the same dzpuis than h^{rt})
    tautst.mult_r_dzpuis(dzfin) ; 	
    
    // Addition of the remaining parts :	
    tautst += 3 * resu(1,2) - resu(1,1).dsdt() ; 
    tautst.mult_sint() ; 
	
    Scalar tmp = resu(1,1) ;
    tmp.mult_cost() ; 		// h^{rr} cos(th)
	
    tautst -= tmp ; 	// T^th / sin(th)
	
    Scalar taut = tautst ; 
    taut.mult_sint() ; 	// T^th
	

    //----------- Computation of T^phi   --> taup : 
    
    Scalar taupst = - resu(1,3).dsdr() ; 

    // dhrr contains  - dh^{rp}/dr in all domains but the CED,  
    // in the CED:    - r^2 dh^{rp}/dr        if dzfin = 0          (3)
    //                - r^3 dh^{rp}/dr        if dzfin = 2          (4)
                                                    	        
    // Multiplication by r of -dh^{rp}/dr  (with the same dzpuis than h^{rp})
    taupst.mult_r_dzpuis(dzfin) ; 	
                          
    // Addition of the remaining part :	
	
    taupst -= 3 * resu(1,3) ; 
    taupst.mult_sint() ; 	// T^ph / sin(th)
	
    Scalar taup = taupst ; 
    taup.mult_sint() ; 		// T^ph 
	
    //------------------- Computation of F and h^[ph ph}
    
    tmp = tautst ; 
    tmp.mult_cost() ; 	// T^th / tan(th)
	
    // dT^th/dth + T^th / tan(th) + 1/sin(th) dT^ph/dph :
    tmp = taut.dsdt() + tmp + taup.stdsdp() ;

    Scalar tmp2 (*mp) ;		
    tmp2 = tmp.poisson_angu() ;  // F
    tmp2.div_sint() ; 
    tmp2.div_sint() ; // h^{ph ph}
    
    idx.set(0) = 3 ;	// phi index
    idx.set(1) = 3 ;	// phi index
    *rcmp[position(idx)] = tmp2 ; 		// h^{ph ph} is updated
	
    
    //------------------- Computation of G and h^[th ph}
    
    tmp = taupst ; 
    tmp.mult_cost() ; // T^ph / tan(th)
    
    // - 1/sin(th) dT^th/dph + dT^ph/dth + T^ph / tan(th) :
    tmp = - taut.stdsdp() + taup.dsdt() + tmp ; 
    
    tmp2 = tmp.poisson_angu() ;  // G
    tmp2.div_sint() ; 
    tmp2.div_sint() ; // h^{th ph}
    
    idx.set(0) = 2 ;	// theta index
    idx.set(1) = 3 ;	// phi index
    *rcmp[position(idx)] = tmp2 ; 		// h^{th ph} is updated
    
    // h^{th th}  (from the trace-free condition)
    // ---------
    idx.set(1) = 2 ;	// theta index
    *rcmp[position(idx)] = - resu(1,1) - resu(3,3) ; 
    
     return resu ;   
}
}
