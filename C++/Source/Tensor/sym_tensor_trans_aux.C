/*
 *  Manipulation of auxiliary potentials for Sym_tensor_trans objects.
 *
 *    (see file sym_tensor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005  Jerome Novak
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

char sym_tensor_trans_aux_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2005/04/06 15:49:46  j_novak
 * Error corrected.
 *
 * Revision 1.1  2005/04/06 15:43:59  j_novak
 * New method Sym_tensor_trans::T_from_det_one(...).
 *
 *
 * $Header$
 *
 */

// C headers
#include <assert.h>

// Lorene headers
#include "tensor.h"

void Sym_tensor_trans::T_from_det_one( const Scalar& hrr, const Scalar& eta_in,
				       const Scalar& mu_in, const Scalar& w_in,
				       const Scalar& x_in) {

    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 
    int dzp = hrr.get_dzpuis() ;
#ifndef NDEBUG
    int dzeta = (dzp == 0 ? 0 : dzp - 1) ;
    assert(eta_in.check_dzpuis(dzeta)) ;
    assert(mu_in.check_dzpuis(dzeta)) ;
    assert(w_in.check_dzpuis(dzp)) ;
    assert(x_in.check_dzpuis(dzp)) ;
#endif
    assert(&eta_in != p_eta) ;
    assert(&mu_in != p_mu) ;
    assert(&w_in != p_www) ;
    assert(&x_in != p_xxx) ;
 
    set(1,1) = hrr ;
    set(1,2) = eta_in.dsdt() - mu_in.stdsdp() ;
    set(1,2).div_r_dzpuis(dzp) ;
    set(1,3) = eta_in.stdsdp() + mu_in.dsdt() ;
    set(1,3).div_r_dzpuis(dzp) ;
    Scalar tmp = w_in.lapang() ;
    tmp.set_spectral_va().ylm_i() ;
    Scalar ppp = 2*w_in.dsdt().dsdt() -  tmp - 2*x_in.stdsdp().dsdt() ;
    tmp = x_in.lapang() ;
    tmp.set_spectral_va().ylm_i() ;
    set(2,3) = 2*x_in.dsdt().dsdt() - tmp + 2*w_in.stdsdp().dsdt() ;

    Scalar A = 0.25*(hrr + 1) ;

    Scalar B = -0.5*( operator()(1,2) * operator()(1,2)
		      + operator()(1,3) * operator()(1,3) ) + hrr + 1. ;

    Scalar C = (hrr + 1)*(ppp*ppp + operator()(2,3)*operator()(2,3)) +
	(ppp + 1)*(operator()(1,3)*operator()(1,3)) +
	(1 - ppp)*(operator()(1,2)*operator()(1,2)) -
	2*operator()(1,2)*operator()(1,3)*operator()(2,3) + hrr ;

    Scalar t_new = (-2*C)/(B + sqrt(B*B - 4*A*C))  ;
    t_new.std_spectral_base() ;

    set(2,2) = 0.5*ttt() + ppp ;
    set(3,3) = 0.5*ttt() - ppp ;

   // Deleting old derived quantities ...
    del_deriv() ; 

    // .. and affecting new ones.
    p_eta = new Scalar(eta_in) ;
    p_mu = new Scalar(mu_in) ;
    p_www = new Scalar(w_in) ;
    p_xxx = new Scalar(x_in) ;
    p_ttt = new Scalar(t_new) ;
    
    

}
