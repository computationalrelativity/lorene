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
 * Revision 1.5  2005/09/07 16:47:43  j_novak
 * Removed method Sym_tensor_trans::T_from_det_one
 * Modified Sym_tensor::set_auxiliary, so that it takes eta/r and mu/r as
 * arguments.
 * Modified Sym_tensor_trans::set_hrr_mu.
 * Added new protected method Sym_tensor_trans::solve_hrr
 *
 * Revision 1.4  2005/04/08 08:22:04  j_novak
 * New methods set_hrr_mu_det_one() and set_WX_det_one(). Not tested yet...
 *
 * Revision 1.3  2005/04/07 07:56:22  j_novak
 * Better handling of dzpuis (first try).
 *
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
#include "graphique.h"


void Sym_tensor_trans::set_hrr_mu_det_one(const Scalar& hrr, const Scalar& mu_in,
					  double precis, int it_max ) {

    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 
    assert(hrr.check_dzpuis(0)) ;
    assert(mu_in.check_dzpuis(0)) ;
    assert(&mu_in != p_mu) ;
    assert( (precis > 0.) && (it_max > 0) ) ;

    Sym_tensor_tt tens_tt(*mp, *triad, *met_div) ;
    tens_tt.set_rr_mu(hrr, mu_in) ;
    tens_tt.inc_dzpuis(2) ;
    trace_from_det_one(tens_tt, precis, it_max) ;
    dec_dzpuis(2) ;

    return ;

}

void Sym_tensor_trans::set_WX_det_one(const Scalar& w_in, const Scalar& x_in,
					  double precis, int it_max ) {

    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 
    assert(w_in.check_dzpuis(0)) ;
    assert(x_in.check_dzpuis(0)) ;
    assert(&w_in != p_www) ;
    assert(&x_in != p_xxx) ;
    assert( (precis > 0.) && (it_max > 0) ) ;

    //Computation of mu from X
    //------------------------
    Scalar source_mu = - x_in.lapang() ;
    source_mu.set_spectral_va().ylm_i() ;
    source_mu -= 2*x_in ;
    source_mu.div_r_dzpuis(1) ;
    Scalar mu_over_r = source_mu.sol_divergence(3) ;
    mu_over_r.annule_l(0,0) ;
    //    mu_over_r.mult_r() ;

    // Preparation for the iteration
    //------------------------------
    Scalar h_old(*mp) ;
    h_old.set_etat_zero() ;
    double lambda = 1. ;
    Scalar w_lapang =  w_in.lapang() ;
    Scalar w_ylm = w_in ;
    w_ylm.set_spectral_va().ylm() ;

    for (int it=0; it<=it_max; it++) {
	
	w_lapang.set_spectral_va().ylm() ;
	Scalar sou_hrr = w_lapang.lapang() + 2*w_lapang ;
	Scalar tmp = h_old.dsdr() ;
	tmp.mult_r_dzpuis(0) ;
	tmp += 3*h_old ;
	tmp.set_spectral_va().ylm() ;
	tmp += 0.5*h_old.lapang() ;
	sou_hrr += tmp ;

	Scalar hrr_new(*mp) ;
	solve_hrr(sou_hrr, hrr_new) ;
	Scalar t_new = h_old - hrr_new ;

//## 	tmp = -2*w_in - 0.5*h_old ;
// 	tmp.set_spectral_va().ylm() ;
// 	tmp -= w_lapang ;
// 	Scalar sou_eta = tmp.dsdr() ;
// 	sou_eta.mult_r_dzpuis(0) ;
// 	sou_eta -= 3*w_lapang + 6*w_ylm  ;
// 	tmp = h_old ;
// 	tmp.set_spectral_va().ylm() ;
// 	sou_eta -= tmp ;
// 	sou_eta.annule_l(0,0, true) ;
// 	Scalar eta_new(*mp) ;
//## 	solve_hrr(sou_eta, eta_new) ;

    	Scalar sou_eta = -hrr_new.dsdr() ;
    	sou_eta.mult_r_dzpuis(0) ;
    	sou_eta -= 2*hrr_new - t_new;
    	Scalar eta_new = sou_eta.poisson_angu() ;

	set_auxiliary(hrr_new, eta_new, mu_over_r, w_in, x_in, t_new) ;

	const Sym_tensor_trans& hij = *this ;
	Scalar h_new = hij(1,1) * hij(2,3) * hij(2,3) 
	    + hij(2,2) * hij(1,3) * hij(1,3) + hij(3,3) * hij(1,2) * hij(1,2)
	    - 2.* hij(1,2) * hij(1,3) * hij(2,3) 
            - hij(1,1) * hij(2,2) * hij(3,3)
	    + hij(1,2) * hij(1,2) + hij(1,3) * hij(1,3) 
	    + hij(2,3) * hij(2,3) - hij(1,1) * hij(2,2) 
	    - hij(1,1) * hij(3,3) - hij(2,2) * hij(3,3) ;
	h_new.set_spectral_base(hrr_new.get_spectral_base()) ;

	Tbl tdif = max(abs(h_new - h_old)) ;
	double diff = max(tdif) ;
        cout << "Sym_tensor_trans::set_WX_det_one : " 
	     << "iteration : " << it << " convergence on h: " 
	     << diff << endl ; 
        if (diff < precis) break ;
        else h_old = lambda*h_new +(1-lambda)*h_old ;

        if (it == it_max) {
            cout << "Sym_tensor_trans:::set_WX_det_one : convergence not reached \n" ;
            cout << "  for the required accuracy (" << precis << ") ! " 
		 << endl ;
            abort() ;
	}
    }

    return ;

}
