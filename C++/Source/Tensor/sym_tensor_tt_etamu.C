/*
 *  Methods of class Sym_tensor_tt related to eta and mu
 *
 *   (see file sym_tensor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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
			//     eta      //
			//--------------//
			
			
const Scalar& Sym_tensor_tt::eta() const {


	if (p_eta == 0x0) {   // a new computation is necessary
		
		// All this has a meaning only for spherical components:
		assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

		// eta is computed from the divergence-free condition:

		Scalar dhrr = - operator()(1,1).dsdr() ; 	// - dh^{rr}/dr  
		                                            // ( - r^2 dh^{rr}/dr in the CED)
		
		// treatment of the CED : 
		int nzm1 = mp->get_mg()->get_nzone() - 1 ; // index of the CED
		Scalar dhrr_ext(*mp) ;
		dhrr_ext.allocate_all() ; 
		dhrr_ext.annule(0,nzm1-1) ; // zero in all domains but the CED
		dhrr_ext.set_domain(nzm1) = dhrr.domain(nzm1) ; // equal to dhrr in the CED
		dhrr_ext.set_spectral_base( (dhrr.get_spectral_va()).get_base() ) ; 
		
		// Multiplication by r^2 in all domains but the CED:
		dhrr.annule_domain(nzm1) ; 
		dhrr.mult_r() ; 
		dhrr.mult_r() ; 
		
		// Adding the CED part
		dhrr.set_domain(nzm1) = dhrr_ext.domain(nzm1) ; 
		dhrr.set_dzpuis(0) ; 

		// dhrr_ext now used to store r h^{rr}
		dhrr_ext = operator()(1,1) ;
		dhrr_ext.mult_r() ; 
		dhrr_ext.inc_dzpuis() ; // since mult_r() did nothing but dzpuis -= 1.
		   
		// Final result for the h^rr source for eta:
		dhrr -= 3. * dhrr_ext ; 
		
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
		tmp.inc_dzpuis() ; // since mult_r() did nothing but dzpuis -= 1.
		
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
	*cmp[position(idx)] = p_eta->srdsdt(0) - p_mu->srstdsdp(0) ; 
	
	// h^{r phi} :
	// ------------
	idx.set(1) = 3 ;	// phi index
	*cmp[position(idx)] = p_eta->srstdsdp(0) + p_mu->srdsdt(0) ; 
	
	
	// h^{theta phi} and h^{phi phi}
	// -----------------------------
	
	Scalar tautst = operator()(1,2).dsdr() ; // dh^{rt}/dr (r^2 dh^{rt}/dr in the CED)
	tautst.dec_dzpuis(1) ; // r dh^{rt}/dr in the CED
	tautst.mult_r() ; // r dh^{rt}/dr  
				   // NB: mult_r() does nothing in the CED but dzpuis--
	
	tautst += 3 * operator()(1,2) - operator()(1,1).dsdt() ; 
	tautst.mult_sint() ; 
	
	Scalar tmp = operator()(1,1) ;
	tmp.mult_cost() ; 		// h^{rr} cos(th)
	
	tautst -= tmp ; 	// T^th / sin(th)
	
	Scalar taut = tautst ; 
	taut.mult_sint() ; 	// T^th
	
	Scalar taupst = - operator()(1,3).dsdr() ; // - dh^{rp}/dr (- r^2 dh^{rp}/dr in the CED)
	taupst.dec_dzpuis(1) ; // - r dh^{rp}/dr in the CED
	taupst.mult_r() ; // - r dh^{rp}/dr  
				   // NB: mult_r() does nothing in the CED but dzpuis--
	
	taupst -= 3 * operator()(1,3) ; 
	taupst.mult_sint() ; 	// T^ph / sin(th)
	
	Scalar taup = taupst ; 
	taup.mult_sint() ; 		// T^ph 
	
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
	
	tmp = taupst ; 
	tmp.mult_cost() ; // T^ph / tan(th)
	
	// 1/sin(th) dT^th/dph - dT^ph/dth - T^ph / tan(th) :
	tmp = taut.stdsdp() - taup.dsdt() - tmp ; 
	
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






