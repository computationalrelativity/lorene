/*
 *  Methods of class Sym_tensor linked with auxiliary members (eta, mu, W, X...)
 *
 *   (see file sym_tensor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2005 Jerome Novak
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


char sym_tensor__aux_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2005/04/05 09:22:15  j_novak
 * Use of the right formula with poisson_angu(2) for the determination of W and
 * X.
 *
 * Revision 1.2  2005/04/04 15:25:24  j_novak
 * Added new members www, xxx, ttt and the associated methods.
 *
 * Revision 1.1  2005/04/01 14:39:57  j_novak
 * Methods dealing with auxilliary derived members of the symmetric
 * tensor (eta, mu, W, X, etc ...).
 *
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// Headers Lorene
#include "metric.h"

    
			//--------------//
			//     eta      //
			//--------------//
			
			
const Scalar& Sym_tensor::eta(Param* par) const {

  if (p_eta == 0x0) {   // a new computation is necessary
	
    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    // eta is computed from its definition (148) of Bonazzola et al. (2004)
    int dzp = operator()(1,1).get_dzpuis() ; 
    int dzp_resu = ((dzp == 0) ? 0 : dzp-1) ;

    Scalar source_eta = operator()(1,2) ;
    source_eta.div_tant() ;
    source_eta += operator()(1,2).dsdt() + operator()(1,3).stdsdp() ;
    source_eta.mult_r_dzpuis(dzp_resu) ;
    
    // Resolution of the angular Poisson equation for eta
    // --------------------------------------------------
    if (dynamic_cast<const Map_af*>(mp) != 0x0) {
      p_eta = new Scalar( source_eta.poisson_angu() ) ; 
    }
    else {
      Scalar resu (*mp) ;
      resu = 0. ;
      mp->poisson_angu(source_eta, *par, resu) ;
      p_eta = new Scalar( resu ) ;  	    
    }
	
  }

  return *p_eta ; 

}

			
			//--------------//
			//     mu       //
			//--------------//
			

const Scalar& Sym_tensor::mu(Param* par) const {

  if (p_mu == 0x0) {   // a new computation is necessary
		
    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    Scalar source_mu = operator()(1,3) ; 	// h^{r ph}
    source_mu.div_tant() ; 		// h^{r ph} / tan(th)
    
    // dh^{r ph}/dth + h^{r ph}/tan(th) - 1/sin(th) dh^{r th}/dphi 
    source_mu += operator()(1,3).dsdt() - operator()(1,2).stdsdp() ; 
    
    // Multiplication by r
    int dzp = operator()(1,2).get_dzpuis() ; 
    int dzp_resu = ((dzp == 0) ? 0 : dzp-1) ;
    source_mu.mult_r_dzpuis(dzp_resu) ; 
    
    // Resolution of the angular Poisson equation for mu
    // --------------------------------------------------
    if (dynamic_cast<const Map_af*>(mp) != 0x0) {
      p_mu = new Scalar( source_mu.poisson_angu() ) ;  
    }
    else {
      Scalar resu (*mp) ;
      resu = 0. ;
      mp->poisson_angu(source_mu, *par, resu) ;
      p_mu = new Scalar( resu ) ;  	    
    }
  }
  return *p_mu ; 

}

			//-------------//
			//     T       //
			//-------------//
			

const Scalar& Sym_tensor::ttt() const {
  
  if (p_ttt == 0x0) { // a new computation is necessary

    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    p_ttt = new Scalar( operator()(2,2) + operator()(3,3) ) ;

  }
  return *p_ttt ;

}

			//------------//
			//     W      //
			//------------//
			
			
const Scalar& Sym_tensor::www() const {

  if (p_www == 0x0) {   // a new computation is necessary
	
    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    Scalar ppp = 0.5*(operator()(2,2) - operator()(3,3)) ;
    Scalar tmp = ppp ;
    tmp.mult_sint() ;
    Scalar source_w = tmp.dsdt().dsdt() ;
    tmp = operator()(2,3) ;
    tmp.mult_sint() ;
    source_w += 2*tmp.dsdt().stdsdp() ;
    source_w.div_sint() ;
    source_w -= ppp.lapang() ;
    
    // Resolution of the angular Poisson equation for W
    p_www = new Scalar( source_w.poisson_angu().poisson_angu(2) ) ; 

  }

  return *p_www ; 

}

			
			//------------//
			//     X      //
			//------------//
			
			
const Scalar& Sym_tensor::xxx() const {

  if (p_xxx == 0x0) {   // a new computation is necessary
	
    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    Scalar tmp = operator()(2,3) ;
    tmp.mult_sint() ;
    Scalar source_x = tmp.dsdt().dsdt() ;
    tmp = 0.5*(operator()(2,2) - operator()(3,3)) ;
    tmp.mult_sint() ;
    source_x -= 2*tmp.dsdt().stdsdp() ;
    source_x.div_sint() ;
    source_x -= operator()(2,3).lapang() ;
    
    // Resolution of the angular Poisson equation for W
    p_xxx = new Scalar( source_x.poisson_angu().poisson_angu(2) ) ;

  }

  return *p_xxx ; 

}

void Sym_tensor::set_auxiliary(const Scalar& trr, const Scalar& eta_in, 
			       const Scalar& mu_in, const Scalar& w_in, 
			       const Scalar& x_in, const Scalar& t_in ) {

    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 
    int dzp = trr.get_dzpuis() ;
#ifndef NDEBUG
    int dzeta = (dzp == 0 ? 0 : dzp - 1) ;
    assert(eta_in.check_dzpuis(dzeta)) ;
    assert(mu_in.check_dzpuis(dzeta)) ;
    assert(w_in.check_dzpuis(dzp)) ;
    assert(x_in.check_dzpuis(dzp)) ;
    assert(t_in.check_dzpuis(dzp)) ;
#endif

    set(1,1) = trr ;
    set(1,2) = eta_in.dsdt() - mu_in.stdsdp() ;
    set(1,2).div_r_dzpuis(dzp) ;
    set(1,3) = eta_in.stdsdp() + mu_in.dsdt() ;
    set(1,3).div_r_dzpuis(dzp) ;
    Scalar ppp = 2*w_in.dsdt().dsdt() - w_in.lapang() - 2*x_in.stdsdp().dsdt() ;
    set(2,3) = 2*x_in.dsdt().dsdt() - x_in.lapang() + 2*w_in.stdsdp().dsdt() ;
    set(2,2) = 0.5*t_in + ppp ;
    set(3,3) = 0.5*t_in - ppp ;

    // Deleting old derived quantities ...
    del_deriv() ; 

    // .. and affecting new ones.
    p_eta = new Scalar(eta_in) ;
    p_mu = new Scalar(mu_in) ;
    p_www = new Scalar(w_in) ;
    p_xxx = new Scalar(x_in) ;
    p_ttt = new Scalar(t_in) ;

}
