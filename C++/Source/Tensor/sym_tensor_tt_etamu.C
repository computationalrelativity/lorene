/*
 *  Methods of class Sym_tensor_tt related to eta and mu
 *
 *   (see file sym_tensor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
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


char sym_tensor_tt_etamu_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2004/03/03 13:16:21  j_novak
 * New potential khi (p_khi) and the functions manipulating it.
 *
 * Revision 1.7  2004/02/05 13:44:50  e_gourgoulhon
 * Major modif. of methods eta(), mu() and update() to treat
 * any value of dzpuis, thanks to the new definitions of
 * Scalar::mult_r(), Scalar::dsdr(), etc...
 *
 * Revision 1.6  2004/01/28 13:25:41  j_novak
 * The ced_mult_r arguments have been suppressed from the Scalar::*dsd* methods.
 * In the div/mult _r_dzpuis, there is no more default value.
 *
 * Revision 1.5  2003/11/05 15:28:31  e_gourgoulhon
 * Corrected error in update.
 *
 * Revision 1.4  2003/11/04 23:03:34  e_gourgoulhon
 * First full version of method update().
 * Add method set_rr_mu.
 * Method set_eta_mu ---> set_rr_eta_mu.
 *
 * Revision 1.3  2003/11/04 09:35:27  e_gourgoulhon
 * First operational version of update_tp().
 *
 * Revision 1.2  2003/11/03 22:33:36  e_gourgoulhon
 * Added methods update_tp and set_eta_mu.
 *
 * Revision 1.1  2003/11/03 17:08:37  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <assert.h>

// Headers Lorene
#include "tensor.h"

			//--------------//
			//     khi      //
			//--------------//

const Scalar& Sym_tensor_tt::khi() const {

  if (p_khi == 0x0) {   // a new computation is necessary
		
    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    // khi is computed from $h^{rr}$ component

    p_khi = new Scalar(operator()(1,1)) ;
    p_khi->mult_r() ;
    p_khi->mult_r() ;
  }

  return *p_khi ; 

}


			//--------------//
			//     eta      //
			//--------------//
			
			
const Scalar& Sym_tensor_tt::eta() const {


	if (p_eta == 0x0) {   // a new computation is necessary
		
		// All this has a meaning only for spherical components:
		assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

		// eta is computed from the divergence-free condition:

        int dzp = operator()(1,1).get_dzpuis() ; 
        
		Scalar dhrr = - operator()(1,1).dsdr() ; 	
        
        // dhrr contains - dh^{rr}/dr in all domains but the CED,                                           
        // in the CED:   - r^2 dh^{rr}/dr        if dzp = 0          (1)
        //               - r^(dzp+1) dh^{rr}/dr  if dzp > 0          (2)
                                                    
		        
        // CED value stored in dhrr_ced : 
		int nzm1 = mp->get_mg()->get_nzone() - 1 ; // index of the CED
        Tbl dhrr_ced =  dhrr.domain(nzm1) ;
        
		// Multiplication by r of (-d h^{rr}/dr) in all domains but the CED:
		dhrr.annule_domain(nzm1) ; 
		dhrr.mult_r() ; 

        // The CED value is restored:
        dhrr.set_domain(nzm1) = dhrr_ced ; 
        
        // Special treatment for dzpuis=2 (cf. Eqs. (1) and (2) above) : 
        if ( dhrr.check_dzpuis(2) )  {
            dhrr.dec_dzpuis() ;   // dhrr contains now -r dh^{rr}/dr in the CED
        }
        
        dhrr.set_dzpuis(dzp) ;  // dhrr is now equal to - r dh^{rr}/dr with the
                                //   same dzpuis than h^{rr}
                          

        // Substraction of the h^rr part and multiplication by r :
        dhrr -= 3. * operator()(1,1) ;                          

        dhrr.mult_r() ; 
        
		
		// Resolution of the angular Poisson equation for eta
		// --------------------------------------------------
		p_eta = new Scalar( dhrr.poisson_angu() ) ; 
	
	}

	return *p_eta ; 

}


			//--------------//
			//     mu       //
			//--------------//
			

const Scalar& Sym_tensor_tt::mu() const {

	if (p_mu == 0x0) {   // a new computation is necessary
		
		// All this has a meaning only for spherical components:
		assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

		Scalar tmp = operator()(1,3) ; 	// h^{r ph}
		tmp.div_tant() ; 		// h^{r ph} / tan(th)
		
		// dh^{r ph}/dth + h^{r ph}/tan(th) - 1/sin(th) dh^{r th}/dphi 
		tmp = operator()(1,3).dsdt() + tmp - operator()(1,2).stdsdp() ; 
		
		// Multiplication by r
		tmp.mult_r() ; 
		
		// Resolution of the angular Poisson equation for mu
		// --------------------------------------------------
		p_mu = new Scalar( tmp.poisson_angu() ) ;  

	}

	return *p_mu ; 

}

			//-------------------//
			//  set_rr_eta_mu    //
			//-------------------//
			

void Sym_tensor_tt::set_rr_eta_mu(const Scalar& hrr, const Scalar& eta_i, 
		const Scalar& mu_i) {

		// All this has a meaning only for spherical components:
		assert( dynamic_cast<const Base_vect_spher*>(triad) != 0x0 ) ; 
						
		set(1,1) = hrr ; 	// h^{rr}
							// calls del_deriv() and therefore delete previous
							// p_eta and p_mu
		
		p_eta = new Scalar( eta_i ) ; 	// eta

		p_mu = new Scalar( mu_i ) ; 	// mu 
		
		update() ; // all h^{ij}, except for h^{rr}
		
}
			
			//---------------//
			//  set_rr_mu    //
			//---------------//
			

void Sym_tensor_tt::set_rr_mu(const Scalar& hrr, const Scalar& mu_i) {

		// All this has a meaning only for spherical components:
		assert( dynamic_cast<const Base_vect_spher*>(triad) != 0x0 ) ; 
						
		set(1,1) = hrr ; 	// h^{rr}
							// calls del_deriv() and therefore delete previous
							// p_eta and p_mu
		
		p_mu = new Scalar( mu_i ) ; 	// mu 
		
		eta() ; // computes eta form the divergence-free condition
		
		update() ; // all h^{ij}, except for h^{rr}
		
}
			
			//-------------------//
			//  set_khi_eta_mu    //
			//-------------------//
			

void Sym_tensor_tt::set_khi_eta_mu(const Scalar& khi_i, const Scalar& eta_i, 
		const Scalar& mu_i) {

  // All this has a meaning only for spherical components:
  assert( dynamic_cast<const Base_vect_spher*>(triad) != 0x0 ) ; 
			
  set(1,1) = khi_i ;
  set(1,1).div_r() ;
  set(1,1).div_r() ;     // h^{rr}

  // calls del_deriv() and therefore delete previous
  // p_khi, p_eta and p_mu
		
  p_khi = new Scalar( khi_i ) ;        // khi

  p_eta = new Scalar( eta_i ) ; 	// eta
  
  p_mu = new Scalar( mu_i ) ; 	// mu 
  
  update() ; // all h^{ij}, except for h^{rr}
		
}
			
			//---------------//
			//  set_khi_mu    //
			//---------------//
			

void Sym_tensor_tt::set_khi_mu(const Scalar& khi_i, const Scalar& mu_i) {

  // All this has a meaning only for spherical components:
  assert( dynamic_cast<const Base_vect_spher*>(triad) != 0x0 ) ; 
						
  set(1,1) = khi_i ; 
                        // calls del_deriv() and therefore delete previous
                        // p_eta and p_mu
  set(1,1).div_r() ;
  set(1,1).div_r() ;	// h^{rr}
		
  p_khi = new Scalar ( khi_i ) ;  // khi

  p_mu = new Scalar( mu_i ) ; 	// mu 
		
  eta() ; // computes eta form the divergence-free condition
  
  update() ; // all h^{ij}, except for h^{rr}
		
}
			

			//-------------//
			//   update    //
			//-------------//
			

void Sym_tensor_tt::update() {

	// All this has a meaning only for spherical components:
	assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

	assert( (p_eta != 0x0) && (p_mu != 0x0) ) ; 
	
	
    Itbl idx(2) ;
    idx.set(0) = 1 ;	// r index
	
	// h^{r theta} : 
	// ------------
	idx.set(1) = 2 ;	// theta index
	*cmp[position(idx)] = p_eta->srdsdt() - p_mu->srstdsdp() ; 
    int dec_dzp = (p_eta->get_dzpuis() == 0) ? 2 : 1 ; 
	cmp[position(idx)]->dec_dzpuis(dec_dzp) ;
	
	// h^{r phi} :
	// ------------
	idx.set(1) = 3 ;	// phi index
	*cmp[position(idx)] = p_eta->srstdsdp() + p_mu->srdsdt() ; 
	cmp[position(idx)]->dec_dzpuis(dec_dzp) ;
	
	
	// h^{theta phi} and h^{phi phi}
	// -----------------------------
	
	//--------------  Computation of T^theta   --> taut : 
    
    Scalar tautst = operator()(1,2).dsdr() ; 

    // dhrr contains  dh^{rt}/dr in all domains but the CED,                                           
    // in the CED:    r^2 dh^{rt}/dr        if dzp = 0          (1)
    //                r^(dzp+1) dh^{rt}/dr  if dzp > 0          (2)
                                                    
		        
    // CED value stored in tmp_ced : 
    int nzm1 = mp->get_mg()->get_nzone() - 1 ; // index of the CED
    Tbl tmp_ced = tautst.domain(nzm1) ;
        
	// Multiplication by r of dh^{rt}/dr in all domains but the CED:
	tautst.annule_domain(nzm1) ; 
    tautst.mult_r() ; 

    // The CED value is restored:
    tautst.set_domain(nzm1) = tmp_ced ; 
        
    // Special treatment for dzpuis=2 (cf. Eqs. (1) and (2) above) : 
    if ( tautst.check_dzpuis(2) )  {
        tautst.dec_dzpuis() ;   // tautst contains now  r dh^{rt}/dr in the CED
    }
        
    tautst.set_dzpuis( operator()(1,2).get_dzpuis() ) ;  
    // taust is now equal to  r dh^{rt}/dr with the same dzpuis than h^{rt}
                          
    // Addition of the remaining parts :	
	tautst += 3 * operator()(1,2) - operator()(1,1).dsdt() ; 
	tautst.mult_sint() ; 
	
	Scalar tmp = operator()(1,1) ;
	tmp.mult_cost() ; 		// h^{rr} cos(th)
	
	tautst -= tmp ; 	// T^th / sin(th)
	
	Scalar taut = tautst ; 
	taut.mult_sint() ; 	// T^th
	

	//----------- Computation of T^phi   --> taup : 
    
	Scalar taupst = - operator()(1,3).dsdr() ; 

    // dhrr contains  - dh^{rp}/dr in all domains but the CED,                                           
    // in the CED:    - r^2 dh^{rp}/dr        if dzp = 0          (3)
    //                - r^(dzp+1) dh^{rp}/dr  if dzp > 0          (4)
                                                    	        
    tmp_ced = taupst.domain(nzm1) ;  // CED value stored in tmp_ced 
        
	// Multiplication by r of -dh^{rp}/dr in all domains but the CED:
	taupst.annule_domain(nzm1) ; 
    taupst.mult_r() ; 

    // The CED value is restored:
    taupst.set_domain(nzm1) = tmp_ced ; 
        
    // Special treatment for dzpuis=2 (cf. Eqs. (3) and (4) above) : 
    if ( taupst.check_dzpuis(2) )  {
        taupst.dec_dzpuis() ;   // taupst contains now  -r dh^{rp}/dr in the CED
    }
        
    taupst.set_dzpuis( operator()(1,3).get_dzpuis() ) ;  
    // taupst is now equal to  -r dh^{rp}/dr with the same dzpuis than h^{rp}
                          
    // Addition of the remaining part :	
	
	taupst -= 3 * operator()(1,3) ; 
	taupst.mult_sint() ; 	// T^ph / sin(th)
	
	Scalar taup = taupst ; 
	taup.mult_sint() ; 		// T^ph 
	
    
    //------------------- Computation of F and h^[ph ph}
    
	tmp = tautst ; 
	tmp.mult_cost() ; 	// T^th / tan(th)
	
	// dT^th/dth + T^th / tan(th) + 1/sin(th) dT^ph/dph :
	tmp = taut.dsdt() + tmp + taup.stdsdp() ;
		
	Scalar tmp2 = tmp.poisson_angu() ;  // F
	
	tmp2.div_sint() ; 
	tmp2.div_sint() ; // h^{ph ph}
	
	idx.set(0) = 3 ;	// phi index
	idx.set(1) = 3 ;	// phi index
	*cmp[position(idx)] = tmp2 ; 		// h^{ph ph} is updated
	
    
    //------------------- Computation of G and h^[th ph}
    
	tmp = taupst ; 
	tmp.mult_cost() ; // T^ph / tan(th)
	
	// - 1/sin(th) dT^th/dph + dT^ph/dth + T^ph / tan(th) :
	tmp = - taut.stdsdp() + taup.dsdt() + tmp ; 
	
	tmp2 = tmp.poisson_angu() ; 	// G
	
	tmp2.div_sint() ; 
	tmp2.div_sint() ; // h^{th ph}
	
	idx.set(0) = 2 ;	// theta index
	idx.set(1) = 3 ;	// phi index
	*cmp[position(idx)] = tmp2 ; 		// h^{th ph} is updated
	
	// h^{th th}  (from the trace-free condition)
	// ---------
	idx.set(1) = 2 ;	// theta index
	*cmp[position(idx)] = - operator()(1,1) - operator()(3,3) ; 
	

	Sym_tensor_trans::del_deriv() ; //## in order not to delete p_eta and p_mu
	


}			






