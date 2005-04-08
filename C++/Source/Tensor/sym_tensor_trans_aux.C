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

    Scalar t_new(*mp) ;

    if (dzp == 0) {
	Scalar A = 0.25*(hrr + 1) ;
	
	Scalar B = -0.5*( operator()(1,2) * operator()(1,2)
		   + operator()(1,3) * operator()(1,3) ) + hrr + 1. ;
	
	Scalar C = (hrr + 1)*(ppp*ppp + operator()(2,3)*operator()(2,3)) +
	    (ppp + 1)*(operator()(1,3)*operator()(1,3)) +
	    (1 - ppp)*(operator()(1,2)*operator()(1,2)) -
	    2*operator()(1,2)*operator()(1,3)*operator()(2,3) + hrr ;
	t_new = (-2*C)/(B + sqrt(B*B - 4*A*C))  ;
    }
    else { //trying to be careful with the dzpuis ...
	Scalar hrr0 = hrr ;
	hrr0.dec_dzpuis(dzp) ;
	Scalar hrt0 = operator()(1,2) ;
	hrt0.dec_dzpuis(dzp) ;
	Scalar hrp0 = operator()(1,3) ;
	hrp0.dec_dzpuis(dzp) ;
	Scalar htp0 = operator()(2,3) ;
	htp0.dec_dzpuis(dzp) ;
	Scalar ppp0 = ppp ;
	ppp0.dec_dzpuis(dzp) ;
	Scalar A = 0.25*(hrr0 + 1) ;
	
	Scalar B = -0.5*( hrt0*hrt0 + hrp0*hrp0 ) + hrr0 + 1. ;
	
	Scalar C = (hrr0 + 1)*(ppp*ppp0 + operator()(2,3)*htp0) +
	    (ppp0 + 1)*(operator()(1,3)*hrp0) +
	    (1 - ppp0)*(operator()(1,2)*hrt0) -
	    2*hrt0*hrp0*operator()(2,3) + hrr ;
	Scalar C0 = (hrr0 + 1)*(ppp0*ppp0 + htp0*htp0) +
	    (ppp0 + 1)*(hrp0*hrp0) + (1 - ppp0)*(hrt0*hrt0) -
	    2*hrt0*hrp0*operator()(2,3) + hrr0 ;
	t_new = (-2*C)/(B + sqrt(B*B - 4*A*C0))  ;
    }

    assert(t_new.check_dzpuis(dzp)) ;
    t_new.std_spectral_base() ;

    set(2,2) = 0.5*t_new + ppp ;
    set(3,3) = 0.5*t_new - ppp ;

   // Deleting old derived quantities ...
    del_deriv() ; 

    // .. and affecting new ones.
    p_eta = new Scalar(eta_in) ;
    p_mu = new Scalar(mu_in) ;
    p_www = new Scalar(w_in) ;
    p_xxx = new Scalar(x_in) ;
    p_ttt = new Scalar(t_new) ;   

}

void Sym_tensor_trans::set_hrr_mu_det_one(const Scalar& hrr, const Scalar& mu_in,
					  double precis, int it_max ) {

    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 
    int dzp = hrr.get_dzpuis() ;
    int dzeta = (dzp == 0 ? 0 : dzp - 1) ;
    assert(mu_in.check_dzpuis(dzeta)) ;
    assert(&mu_in != p_mu) ;
    assert( (precis > 0.) && (it_max > 0) ) ;

    //Computation of X from mu
    //------------------------
    Scalar source_x = (-2.)*mu_in ;
    source_x.div_r_dzpuis(dzp) ;
    Scalar tmp = mu_in.dsdr() ;
    if (dzp == 0)
	tmp.dec_dzpuis(2) ;
    source_x -= tmp ;
    Scalar x_new = source_x.poisson_angu(2) ;

    // Preparation for the iteration
    //------------------------------
    Scalar T_old = - hrr ;
    Scalar dhrr = hrr.dsdr() ;
    dhrr.mult_r_dzpuis(dzp) ;
    dhrr += 2*hrr ;
    dhrr.mult_r_dzpuis(dzeta) ;

    for (int it=0; it<=it_max; it++) {
	
	Scalar source_eta = T_old ;
	source_eta.mult_r_dzpuis(dzeta) ;
	source_eta -= dhrr ;
	Scalar eta_new = source_eta.poisson_angu() ;

	Scalar source_w = (-2.)*eta_new ;
	source_w.div_r_dzpuis(dzp) ;
	tmp = eta_new.dsdr() ;
	if (dzp == 0)
	    tmp.dec_dzpuis(2) ;
	source_w -= tmp + 0.5*T_old ;
	Scalar w_new = source_w.poisson_angu(2) ;

	T_from_det_one(hrr, eta_new, mu_in, w_new, x_new) ;

	double diff = max(max(abs(ttt() - T_old))) ;
        cout << "Sym_tensor_trans::set_hrr_mu_det_one : " 
	     << "iteration : " << it << " convergence on T: " << diff << endl ;
        if (diff < precis) break ;
        else T_old = ttt() ;

        if (it == it_max) {
            cout << "Sym_tensor_trans:::set_hrr_mu_det_one : convergence not reached \n" ;
            cout << "  for the required accuracy (" << precis << ") ! " << endl ;
            abort() ;
	}
    }


}

void Sym_tensor_trans::set_WX_det_one(const Scalar& w_in, const Scalar& x_in,
					  double precis, int it_max ) {

    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 
    int dzp = w_in.get_dzpuis() ;
    int dzeta = (dzp == 0 ? 0 : dzp - 1) ;
    assert(x_in.check_dzpuis(dzp)) ;
    assert(&w_in != p_www) ;
    assert(&x_in != p_xxx) ;
    assert( (precis > 0.) && (it_max > 0) ) ;

    //Computation of mu from X
    //------------------------
    Scalar source_mu = 2*x_in + x_in.lapang() ;
    int dz_tmp = (dzp <2 ? 1 : dzp - 1) ; 
    source_mu.mult_r_dzpuis(dz_tmp) ;
    source_mu.mult_r_dzpuis(2) ;
    Scalar mu_new = source_mu.primr() ;
    mu_new.div_r_dzpuis((dzeta+1)/2) ;
    mu_new.div_r_dzpuis(dzeta) ;

    // Preparation for the iteration
    //------------------------------
    Scalar T_old(*mp) ;
    T_old.set_etat_zero() ;

    for (int it=0; it<=it_max; it++) {
	
	Scalar source_eta = 0.5*T_old + 2*w_in + w_in.lapang() ;
	dz_tmp = (dzp <2 ? 1 : dzp - 1) ; 
	source_eta.mult_r_dzpuis(dz_tmp) ;
	source_eta.mult_r_dzpuis(2) ;
	Scalar eta_new = source_eta.primr() ;
	eta_new.div_r_dzpuis((dzeta+1)/2) ;
	eta_new.div_r_dzpuis(dzeta) ;

	Scalar source_hrr = T_old ;
	source_hrr.mult_r_dzpuis(dzeta) ;
	source_hrr -= eta_new.lapang() ;
	dzeta >2 ? source_hrr.dec_dzpuis(dzeta-2) : source_hrr.inc_dzpuis(2-dzeta) ;
	Scalar hrr_new = source_hrr.primr() ;
	hrr_new.div_r_dzpuis((dzp+1)/2) ;
	hrr_new.div_r_dzpuis(dzp) ;

	T_from_det_one(hrr_new, eta_new, mu_new, w_in, x_in) ;

	double diff = max(max(abs(ttt() - T_old))) ;
        cout << "Sym_tensor_trans::set_WX_det_one : " 
	     << "iteration : " << it << " convergence on T: " << diff << endl ;
        if (diff < precis) break ;
        else T_old = ttt() ;

        if (it == it_max) {
            cout << "Sym_tensor_trans:::set_WX_det_one : convergence not reached \n" ;
            cout << "  for the required accuracy (" << precis << ") ! " << endl ;
            abort() ;
	}
    }


}
