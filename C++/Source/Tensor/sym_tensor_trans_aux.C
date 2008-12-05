/*
 *  Manipulation of auxiliary potentials for Sym_tensor_trans objects.
 *
 *    (see file sym_tensor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005-2006  Jerome Novak
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
 * Revision 1.18  2008/12/05 08:46:19  j_novak
 * New method Sym_tensor_trans_aux::set_tt_part_det_one.
 *
 * Revision 1.17  2007/04/27 13:52:55  j_novak
 * In method Sym_tensor_trans::set_AtBtt_det_one(), the member p_tt (the TT-part)
 * is updated.
 *
 * Revision 1.16  2007/03/20 12:20:56  j_novak
 * In Sym_tensor_trans::set_AtBtt_det_one(), the trace is stored in the resulting
 * object.
 *
 * Revision 1.15  2006/10/24 13:03:19  j_novak
 * New methods for the solution of the tensor wave equation. Perhaps, first
 * operational version...
 *
 * Revision 1.14  2006/09/15 08:48:01  j_novak
 * Suppression of some messages in the optimized version.
 *
 * Revision 1.13  2006/08/31 12:13:22  j_novak
 * Added an argument of type Param to Sym_tensor_trans::sol_Dirac_A().
 *
 * Revision 1.12  2006/06/26 09:28:13  j_novak
 * Added a forgotten initialisation in set_AtB_trace_zero().
 *
 * Revision 1.11  2006/06/20 12:07:15  j_novak
 * Improved execution speed for sol_Dirac_tildeB...
 *
 * Revision 1.10  2006/06/14 10:04:21  j_novak
 * New methods sol_Dirac_l01, set_AtB_det_one and set_AtB_trace_zero.
 *
 * Revision 1.9  2006/01/17 15:50:53  j_novak
 * Slight re-arrangment of the det=1 formula.
 *
 * Revision 1.8  2005/11/28 14:45:17  j_novak
 * Improved solution of the Poisson tensor equation in the case of a transverse
 * tensor.
 *
 * Revision 1.7  2005/09/16 13:34:44  j_novak
 * Back to dzpuis 1 for the source for mu. eta is computed the same way as hrr.
 *
 * Revision 1.6  2005/09/08 16:00:23  j_novak
 * Change of dzpuis for source for mu (temporary?).
 *
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
#include "metric.h"
#include "param.h"

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
				      const Scalar* h_prev, double precis, int it_max ) {

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

    // Preparation for the iteration
    //------------------------------
    Scalar h_old(*mp) ;
    if (h_prev != 0x0) 
	h_old = *h_prev ;
    else
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

 	tmp = -2*w_in - 0.5*h_old ;
 	tmp.set_spectral_va().ylm() ;
 	tmp -= w_lapang ;
 	Scalar sou_eta = tmp.dsdr() ;
 	sou_eta.mult_r_dzpuis(0) ;
 	sou_eta -= 3*w_lapang + 6*w_ylm  ;
 	tmp = h_old ;
 	tmp.set_spectral_va().ylm() ;
 	sou_eta -= tmp ;
 	sou_eta.annule_l(0,0, true) ;
 	Scalar eta_new(*mp) ;
 	solve_hrr(sou_eta, eta_new) ;

	set_auxiliary(hrr_new, eta_new, mu_over_r, w_in, x_in, t_new) ;

	const Sym_tensor_trans& hij = *this ;
	Scalar h_new = (1 + hij(1,1))*( hij(2,3)*hij(2,3) - hij(2,2)*hij(3,3) )
	    + hij(1,2)*hij(1,2)*(1 + hij(3,3)) 
	    + hij(1,3)*hij(1,3)*(1 + hij(2,2)) 
	    - hij(1,1)*(hij(2,2) + hij(3,3)) - 2*hij(1,2)*hij(1,3)*hij(2,3) ;
	h_new.set_spectral_base(hrr_new.get_spectral_base()) ;

	Tbl tdif = max(abs(h_new - h_old)) ;
	double diff = max(tdif) ;
#ifndef NDEBUG
        cout << "Sym_tensor_trans::set_WX_det_one : " 
	     << "iteration : " << it << " convergence on h: " 
	     << diff << endl ; 
#endif
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

void Sym_tensor_trans::set_AtBtt_det_one(const Scalar& a_in, const Scalar& tbtt_in,
					 const Scalar* h_prev, Param* par_bc, 
					 Param* par_mat, double precis, 
					 int it_max ) {
    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 
    assert(a_in.check_dzpuis(2)) ;
    assert(tbtt_in.check_dzpuis(2)) ;
    assert(&a_in != p_aaa) ;
    assert(&tbtt_in != p_tilde_b) ;
    assert( (precis > 0.) && (it_max > 0) ) ;

    //Computation of mu and X from A
    //-------------------------------
    Scalar mu_over_r(*mp) ;
    Scalar x_new(*mp) ;
    sol_Dirac_A(a_in, mu_over_r, x_new, par_bc) ;

    // Preparation for the iteration
    //------------------------------
    Scalar h_old(*mp) ;
    if (h_prev != 0x0) 
	h_old = *h_prev ;
    else 
	h_old.set_etat_zero() ;
    double lambda = 1. ;
    Scalar zero(*mp) ;
    zero.set_etat_zero() ;

    Scalar hrr_tt(*mp) ;
    Scalar eta_sr_tt(*mp) ;
    Scalar w_tt(*mp) ;
    sol_Dirac_tilde_B(tbtt_in, zero, hrr_tt, eta_sr_tt, w_tt, par_bc, par_mat) ;
    Sym_tensor_tt hijtt(*mp, *triad, *met_div) ;
    hijtt.set_auxiliary(hrr_tt, eta_sr_tt, mu_over_r, w_tt, x_new, zero) ;

    Scalar hrr_new(*mp) ;
    Scalar eta_over_r(*mp) ;
    Scalar w_new(*mp) ;

    for (int it=0; it<=it_max; it++) {
	Scalar tilde_B = get_tilde_B_from_TT_trace(zero, h_old) ;
	sol_Dirac_tilde_B(tilde_B, h_old, hrr_new, eta_over_r, w_new, 0x0, par_mat) ;
	
	set_auxiliary(hrr_new+hrr_tt, eta_over_r+eta_sr_tt, mu_over_r, 
		      w_new+w_tt, x_new, h_old - hrr_new-hrr_tt) ;

	const Sym_tensor_trans& hij = *this ;
	Scalar h_new = (1 + hij(1,1))*( hij(2,3)*hij(2,3) - hij(2,2)*hij(3,3) )
	    + hij(1,2)*hij(1,2)*(1 + hij(3,3)) 
	    + hij(1,3)*hij(1,3)*(1 + hij(2,2)) 
	    - hij(1,1)*(hij(2,2) + hij(3,3)) - 2*hij(1,2)*hij(1,3)*hij(2,3) ;
	h_new.set_spectral_base(hrr_tt.get_spectral_base()) ;

	double diff = max(max(abs(h_new - h_old))) ;
#ifndef NDEBUG
        cout << "Sym_tensor_trans::set_AtB_det_one : " 
	     << "iteration : " << it << " convergence on h: " 
	     << diff << endl ; 
#endif
        if (diff < precis) {
  	    set_auxiliary(hrr_new+hrr_tt, eta_over_r+eta_sr_tt, mu_over_r, 
		      w_new+w_tt, x_new, h_old - hrr_new-hrr_tt) ;
	    if (p_aaa != 0x0) delete p_aaa ;
	    p_aaa = new Scalar(a_in) ;
	    if (p_tilde_b != 0x0) delete p_tilde_b ;
	    p_tilde_b = new Scalar(tilde_B + tbtt_in) ;
	    if (p_trace != 0x0) delete p_trace ;
	    p_trace = new Scalar(h_old) ;
	    if (p_tt != 0x0) delete p_tt ;
	    p_tt = new Sym_tensor_tt(hijtt) ;
	    p_tt->inc_dzpuis(2) ;
	    break ;
	}
        else {
	    h_old = lambda*h_new +(1-lambda)*h_old ;
	}

        if (it == it_max) {
            cout << "Sym_tensor_trans:::set_AtBtt_det_one : convergence not reached \n" ;
            cout << "  for the required accuracy (" << precis << ") ! " 
		 << endl ;
            abort() ;
	}
    }
    return ;

}

void Sym_tensor_trans::set_tt_part_det_one(const Sym_tensor_tt& hijtt, const 
					   Scalar* h_prev, Param* par_mat, 
					   double precis, int it_max ) {
    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 
    assert( (precis > 0.) && (it_max > 0) ) ;

    Scalar mu_over_r = hijtt.mu() ;
    mu_over_r.div_r() ;
    const Scalar& x_new = hijtt.xxx() ;

    // Preparation for the iteration
    //------------------------------
    Scalar h_old(*mp) ;
    if (h_prev != 0x0) 
	h_old = *h_prev ;
    else 
	h_old.set_etat_zero() ;
    double lambda = 1. ;
    Scalar zero(*mp) ;
    zero.set_etat_zero() ;

    const Scalar& hrr_tt = hijtt( 1, 1 ) ;
    Scalar eta_sr_tt = hijtt.eta() ;
    eta_sr_tt.div_r() ;
    const Scalar w_tt = hijtt.www() ;

    Scalar hrr_new(*mp) ;
    Scalar eta_over_r(*mp) ;
    Scalar w_new(*mp) ;

    for (int it=0; it<=it_max; it++) {
	Scalar tilde_B = get_tilde_B_from_TT_trace(zero, h_old) ;
	sol_Dirac_tilde_B(tilde_B, h_old, hrr_new, eta_over_r, w_new, 0x0, par_mat) ;
	
	set_auxiliary(hrr_new+hrr_tt, eta_over_r+eta_sr_tt, mu_over_r, 
		      w_new+w_tt, x_new, h_old - hrr_new-hrr_tt) ;

	const Sym_tensor_trans& hij = *this ;
	Scalar h_new = (1 + hij(1,1))*( hij(2,3)*hij(2,3) - hij(2,2)*hij(3,3) )
	    + hij(1,2)*hij(1,2)*(1 + hij(3,3)) 
	    + hij(1,3)*hij(1,3)*(1 + hij(2,2)) 
	    - hij(1,1)*(hij(2,2) + hij(3,3)) - 2*hij(1,2)*hij(1,3)*hij(2,3) ;
	h_new.set_spectral_base(hrr_tt.get_spectral_base()) ;

	double diff = max(max(abs(h_new - h_old))) ;
#ifndef NDEBUG
        cout << "Sym_tensor_trans::set_tt_part_det_one : " 
	     << "iteration : " << it << " convergence on h: " 
	     << diff << endl ; 
#endif
        if (diff < precis) {
  	    set_auxiliary(hrr_new+hrr_tt, eta_over_r+eta_sr_tt, mu_over_r, 
		      w_new+w_tt, x_new, h_old - hrr_new-hrr_tt) ;
	    if (p_trace != 0x0) delete p_trace ;
	    p_trace = new Scalar(h_old) ;
	    if (p_tt != 0x0) delete p_tt ;
	    p_tt = new Sym_tensor_tt(hijtt) ;
	    p_tt->inc_dzpuis(2) ;
	    break ;
	}
        else {
	    h_old = lambda*h_new +(1-lambda)*h_old ;
	}

        if (it == it_max) {
            cout << "Sym_tensor_trans:::set_AtBtt_det_one : convergence not reached \n" ;
            cout << "  for the required accuracy (" << precis << ") ! " 
		 << endl ;
            abort() ;
	}
    }
    return ;

}

void Sym_tensor_trans::set_AtB_trace(const Scalar& a_in, const Scalar& tb_in,
				     const Scalar& hh, Param* par_bc, Param* par_mat) {

    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 
    assert(a_in.check_dzpuis(2)) ;
    assert(tb_in.check_dzpuis(2)) ;
    assert(hh.check_dzpuis(0)) ;
    assert(&a_in != p_aaa) ;
    assert(&tb_in != p_tilde_b) ;

    //Computation of mu and X from A
    //-------------------------------
    Scalar mu_over_r(*mp) ;
    Scalar x_new(*mp) ;
    sol_Dirac_A(a_in, mu_over_r, x_new, par_bc) ;

    // Computation of the other potentials
    //------------------------------------
    Scalar hrr_new(*mp) ;
    Scalar eta_over_r(*mp) ;
    Scalar w_new(*mp) ;

    sol_Dirac_tilde_B(tb_in, hh, hrr_new, eta_over_r, w_new, par_bc, par_mat) ;

    set_auxiliary(hrr_new, eta_over_r, mu_over_r, w_new, x_new, hh - hrr_new) ;
    if (p_aaa != 0x0) delete p_aaa ;
    p_aaa = new Scalar(a_in) ;
    if (p_tilde_b != 0x0) delete p_tilde_b ;
    p_tilde_b = new Scalar(tb_in) ;

    return ;

}
