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
 * Revision 1.15  2014/10/13 08:53:43  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.14  2014/10/06 15:13:19  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.13  2007/11/27 15:49:51  n_vasset
 * new function compute_tilde_C in class sym_tensor
 *
 * Revision 1.12  2006/10/24 14:52:38  j_novak
 * The Mtbl corresponding to the physical space is destroyed after the
 * calculation of tilde_B(tt), to get the updated result.
 *
 * Revision 1.11  2006/10/24 13:03:19  j_novak
 * New methods for the solution of the tensor wave equation. Perhaps, first
 * operational version...
 *
 * Revision 1.10  2006/08/31 12:12:43  j_novak
 * Correction of a small mistake in a phi loop.
 *
 * Revision 1.9  2006/06/28 07:48:26  j_novak
 * Better treatment of some null cases.
 *
 * Revision 1.8  2006/06/12 11:40:07  j_novak
 * Better memory and spectral base handling for Sym_tensor::compute_tilde_B.
 *
 * Revision 1.7  2006/06/12 07:42:29  j_novak
 * Fields A and tilde{B} are defined only for l>1.
 *
 * Revision 1.6  2006/06/12 07:27:20  j_novak
 * New members concerning A and tilde{B}, dealing with the transverse part of the
 * Sym_tensor.
 *
 * Revision 1.5  2005/09/07 16:47:43  j_novak
 * Removed method Sym_tensor_trans::T_from_det_one
 * Modified Sym_tensor::set_auxiliary, so that it takes eta/r and mu/r as
 * arguments.
 * Modified Sym_tensor_trans::set_hrr_mu.
 * Added new protected method Sym_tensor_trans::solve_hrr
 *
 * Revision 1.4  2005/04/05 15:38:08  j_novak
 * Changed the formulas for W and X. There was an error before...
 *
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
#include <cstdlib>
#include <cassert>
#include <cmath>

// Headers Lorene
#include "metric.h"
#include "proto.h"

    
			//--------------//
			//     eta      //
			//--------------//
			
			
namespace Lorene {
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
    Scalar tmp = ppp.dsdt() ;
    tmp.div_tant() ;
    Scalar source_w = ppp.dsdt().dsdt() +3*tmp - ppp.stdsdp().stdsdp() -2*ppp ;
    tmp = operator()(2,3) ;
    tmp.div_tant() ;
    tmp += operator()(2,3).dsdt() ;
    source_w += 2*tmp.stdsdp() ;
    
    // Resolution of the angular Poisson equation for W
    p_www = new Scalar( source_w.poisson_angu(2).poisson_angu() ) ; 

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

    const Scalar& htp = operator()(2,3) ;
    Scalar tmp = htp.dsdt() ;
    tmp.div_tant() ;
    Scalar source_x = htp.dsdt().dsdt() + 3*tmp - htp.stdsdp().stdsdp() -2*htp ;
    Scalar ppp = 0.5*(operator()(2,2) - operator()(3,3)) ;
    tmp = ppp ;
    tmp.div_tant() ;
    tmp += ppp.dsdt() ;
    source_x -= 2*tmp.stdsdp() ;
    
    // Resolution of the angular Poisson equation for W
    p_xxx = new Scalar( source_x.poisson_angu(2).poisson_angu() ) ;

  }

  return *p_xxx ; 

}

void Sym_tensor::set_auxiliary(const Scalar& trr, const Scalar& eta_over_r, 
			       const Scalar& mu_over_r, const Scalar& w_in, 
			       const Scalar& x_in, const Scalar& t_in ) {

    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 
    int dzp = trr.get_dzpuis() ;
    int dzeta = (dzp == 0 ? 0 : dzp - 1) ;

    assert(eta_over_r.check_dzpuis(dzp)) ;
    assert(mu_over_r.check_dzpuis(dzp)) ;
    assert(w_in.check_dzpuis(dzp)) ;
    assert(x_in.check_dzpuis(dzp)) ;
    assert(t_in.check_dzpuis(dzp)) ;

    assert(&w_in != p_www) ;
    assert(&x_in != p_xxx) ;
    assert(&t_in != p_ttt) ;

    set(1,1) = trr ;
    set(1,2) = eta_over_r.dsdt() - mu_over_r.stdsdp() ;
    //   set(1,2).div_r_dzpuis(dzp) ;
    set(1,3) = eta_over_r.stdsdp() + mu_over_r.dsdt() ;
    //    set(1,3).div_r_dzpuis(dzp) ;
    Scalar tmp = w_in.lapang() ;
    tmp.set_spectral_va().ylm_i() ;
    Scalar ppp = 2*w_in.dsdt().dsdt() -  tmp - 2*x_in.stdsdp().dsdt() ;
    tmp = x_in.lapang() ;
    tmp.set_spectral_va().ylm_i() ;
    set(2,3) = 2*x_in.dsdt().dsdt() - tmp + 2*w_in.stdsdp().dsdt() ;
    set(2,2) = 0.5*t_in + ppp ;
    set(3,3) = 0.5*t_in - ppp ;

    // Deleting old derived quantities ...
    del_deriv() ; 

    // .. and affecting new ones.
    p_eta = new Scalar(eta_over_r) ; p_eta->mult_r_dzpuis(dzeta) ;
    p_mu = new Scalar(mu_over_r) ; p_mu->mult_r_dzpuis(dzeta) ;
    p_www = new Scalar(w_in) ;
    p_xxx = new Scalar(x_in) ;
    p_ttt = new Scalar(t_in) ;

}

			//------------//
			//     A      //
			//------------//
			
			
const Scalar& Sym_tensor::compute_A(bool output_ylm, Param* par) const {

  if (p_aaa == 0x0) {   // a new computation is necessary
	
    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    int dzp = xxx().get_dzpuis() ; 
    int dzp_resu = ((dzp == 0) ? 2 : dzp+1) ;
 
    Scalar source_mu = operator()(1,3) ; 	// h^{r ph}
    source_mu.div_tant() ; 		// h^{r ph} / tan(th)
    
    // dh^{r ph}/dth + h^{r ph}/tan(th) - 1/sin(th) dh^{r th}/dphi 
    source_mu += operator()(1,3).dsdt() - operator()(1,2).stdsdp() ; 
    
    // Resolution of the angular Poisson equation for mu / r
    // -----------------------------------------------------
    Scalar tilde_mu(*mp) ;
    tilde_mu = 0 ;

    if (dynamic_cast<const Map_af*>(mp) != 0x0) {
	tilde_mu = source_mu.poisson_angu()  ;  
    }
    else {
	assert(par != 0x0) ;
	mp->poisson_angu(source_mu, *par, tilde_mu) ;
    }
    tilde_mu.div_r_dzpuis(dzp_resu) ;
    
    p_aaa = new Scalar( xxx().dsdr() - tilde_mu) ;
    p_aaa->annule_l(0, 1, output_ylm) ;
  }

  if (output_ylm) p_aaa->set_spectral_va().ylm() ;
  else  p_aaa->set_spectral_va().ylm_i() ;

  return *p_aaa ; 

}

			//--------------//
			//   tilde(B)   //
			//--------------//
			
			
const Scalar& Sym_tensor::compute_tilde_B(bool output_ylm, Param* par) const {

  if (p_tilde_b == 0x0) {   // a new computation is necessary
	
    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    int dzp = operator()(1,1).get_dzpuis() ; 
    int dzp_resu = ((dzp == 0) ? 2 : dzp+1) ;

    p_tilde_b = new Scalar(*mp) ;
    p_tilde_b->set_etat_qcq() ;

    Scalar source_eta = operator()(1,2) ;
    source_eta.div_tant() ;
    source_eta += operator()(1,2).dsdt() + operator()(1,3).stdsdp() ;
    
    // Resolution of the angular Poisson equation for eta / r
    // ------------------------------------------------------
    Scalar tilde_eta(*mp) ;
    tilde_eta = 0 ;

    if (dynamic_cast<const Map_af*>(mp) != 0x0) {
	tilde_eta = source_eta.poisson_angu() ; 
    }
    else {
	assert (par != 0x0) ;
	mp->poisson_angu(source_eta, *par, tilde_eta) ;
    }
 
    Scalar dwdr = www().dsdr() ;
    dwdr.set_spectral_va().ylm() ;
    Scalar wsr = www() ;
    wsr.div_r_dzpuis(dzp_resu) ;
    wsr.set_spectral_va().ylm() ;
    Scalar etasr2 = tilde_eta ;
    etasr2.div_r_dzpuis(dzp_resu) ;
    etasr2.set_spectral_va().ylm() ;
    Scalar dtdr = ttt().dsdr() ;
    dtdr.set_spectral_va().ylm() ;
    Scalar tsr = ttt() ;
    tsr.div_r_dzpuis(dzp_resu) ;
    tsr.set_spectral_va().ylm() ;
    Scalar hrrsr = operator()(1,1) ;
    hrrsr.div_r_dzpuis(dzp_resu) ;
    hrrsr.set_spectral_va().ylm() ;

    int nz = mp->get_mg()->get_nzone() ;

    Base_val base(nz) ;

    if (wsr.get_etat() != ETATZERO) {
	base = wsr.get_spectral_base() ;
    }
    else {
	if (etasr2.get_etat() != ETATZERO) {
	    base = etasr2.get_spectral_base() ;
	}
	else {
	    if (tsr.get_etat() != ETATZERO) {
		base = tsr.get_spectral_base() ;
	    }
	    else {
		if (hrrsr.get_etat() != ETATZERO) {
		    base = hrrsr.get_spectral_base() ;
		}
		else { //Everything is zero...
		    p_tilde_b->set_etat_zero() ;
		    return *p_tilde_b ;
		}
	    }
	}
    }
    
    p_tilde_b->set_spectral_base(base) ;
    p_tilde_b->set_spectral_va().set_etat_cf_qcq() ;
    p_tilde_b->set_spectral_va().c_cf->annule_hard() ;

    if (dwdr.get_etat() == ETATZERO) dwdr.annule_hard() ;
    if (wsr.get_etat() == ETATZERO) wsr.annule_hard() ;
    if (etasr2.get_etat() == ETATZERO) etasr2.annule_hard() ;
    if (dtdr.get_etat() == ETATZERO) dtdr.annule_hard() ;
    if (tsr.get_etat() == ETATZERO) tsr.annule_hard() ;
    if (hrrsr.get_etat() == ETATZERO) hrrsr.annule_hard() ;

    int m_q, l_q, base_r ;
    for (int lz=0; lz<nz; lz++) {
	int np = mp->get_mg()->get_np(lz) ;
	int nt = mp->get_mg()->get_nt(lz) ;
	int nr = mp->get_mg()->get_nr(lz) ;
	for (int k=0; k<np+1; k++)
	    for (int j=0; j<nt; j++) {
		base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
		if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1))
		{
		    for (int i=0; i<nr; i++) 
			p_tilde_b->set_spectral_va().c_cf->set(lz, k, j, i)
			    = (l_q+2)*(*dwdr.get_spectral_va().c_cf)(lz, k, j, i)
			    + l_q*(l_q+2)*(*wsr.get_spectral_va().c_cf)(lz, k, j, i)
			    - 2*(*etasr2.get_spectral_va().c_cf)(lz, k, j, i)
	    + 0.5*double(l_q+2)/double(l_q+1)*(*tsr.get_spectral_va().c_cf)(lz, k, j, i)
	    + 0.5/double(l_q+1)*(*dtdr.get_spectral_va().c_cf)(lz, k, j, i)
	    - 1./double(l_q+1)*(*hrrsr.get_spectral_va().c_cf)(lz, k, j, i) ;
		}
	    }
    }
    p_tilde_b->set_dzpuis(dzp_resu) ;
  } //End of p_tilde_b == 0x0    

  if (output_ylm) p_tilde_b->set_spectral_va().ylm() ;
  else  p_tilde_b->set_spectral_va().ylm_i() ;

  return *p_tilde_b ; 

}

Scalar Sym_tensor::compute_tilde_B_tt(bool output_ylm, Param* par) const {

    Scalar resu = compute_tilde_B(true, par) ;
    if (resu.get_etat() == ETATZERO) return resu ;
    else {
	assert(resu.get_etat() == ETATQCQ) ;
	assert(resu.get_spectral_va().c_cf != 0x0) ;
	int dzp = operator()(1,1).get_dzpuis() ; 
	int dzp_resu = ((dzp == 0) ? 2 : dzp+1) ;
	
	Scalar hsr = operator()(1,1) + ttt() ;
	if (hsr.get_etat() == ETATZERO) return resu ;
	Scalar dhdr = hsr.dsdr() ;
	dhdr.set_spectral_va().ylm() ;
	hsr.div_r_dzpuis(dzp_resu) ;
	hsr.set_spectral_va().ylm() ;
	
	int nz = mp->get_mg()->get_nzone() ;

	const Base_val& base = resu.get_spectral_base() ;
    
	if (dhdr.get_etat() == ETATZERO) dhdr.annule_hard() ;

	int m_q, l_q, base_r ;
	for (int lz=0; lz<nz; lz++) {
	    int np = mp->get_mg()->get_np(lz) ;
	int nt = mp->get_mg()->get_nt(lz) ;
	int nr = mp->get_mg()->get_nr(lz) ;
	for (int k=0; k<np+1; k++)
	    for (int j=0; j<nt; j++) {
		base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
		if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1))
		{
		    for (int i=0; i<nr; i++) 
			resu.set_spectral_va().c_cf->set(lz, k, j, i)
	    -= 	0.5*(*hsr.get_spectral_va().c_cf)(lz, k, j, i)
	      + 0.5/double(l_q+1)*(*dhdr.get_spectral_va().c_cf)(lz, k, j, i);
		}
	    }
	}
	resu.set_dzpuis(dzp_resu) ;
	if (resu.set_spectral_va().c != 0x0) 
	    delete resu.set_spectral_va().c ;
	resu.set_spectral_va().c = 0x0 ;
    } //End of resu != 0    

  if (output_ylm) resu.set_spectral_va().ylm() ;
  else  resu.set_spectral_va().ylm_i() ;

  return resu ; 

}

Scalar Sym_tensor::get_tilde_B_from_TT_trace(const Scalar& tbtt, const Scalar&
    tras) const {
    
    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    Scalar resu = tbtt ;
    if (resu.get_etat() == ETATZERO) {
	if (tras.get_etat() == ETATZERO) return resu ;
	else {
	    resu.annule_hard() ;
	    Base_val base = tras.get_spectral_base() ;
	    base.mult_x() ;
	    resu.set_spectral_base(base) ;
	}
    }
    resu.set_spectral_va().ylm() ;
    int dzp = operator()(1,1).get_dzpuis() ; 
    int dzp_resu = ((dzp == 0) ? 2 : dzp+1) ;
    
    Scalar hsr = tras ;
    if (hsr.get_etat() == ETATZERO) return resu ;
    else {
	Scalar dhdr = hsr.dsdr() ;
	if (dhdr.get_etat() == ETATZERO) dhdr.annule_hard() ;
	dhdr.set_spectral_va().ylm() ;
	hsr.div_r_dzpuis(dzp_resu) ;
	hsr.set_spectral_va().ylm() ;
	
	int nz = mp->get_mg()->get_nzone() ;	
	const Base_val& base = resu.get_spectral_base() ;    

	int m_q, l_q, base_r ;
	for (int lz=0; lz<nz; lz++) {
	    if ((*resu.get_spectral_va().c_cf)(lz).get_etat() == ETATZERO)
		resu.get_spectral_va().c_cf->set(lz).annule_hard() ;
	    int np = mp->get_mg()->get_np(lz) ;
	    int nt = mp->get_mg()->get_nt(lz) ;
	    int nr = mp->get_mg()->get_nr(lz) ;
	    for (int k=0; k<np+1; k++)
		for (int j=0; j<nt; j++) {
		    base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
		    if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1))
		    {
			for (int i=0; i<nr; i++) 
			    resu.set_spectral_va().c_cf->set(lz, k, j, i)
	    += 	0.5*(*hsr.get_spectral_va().c_cf)(lz, k, j, i)
	      + 0.5/double(l_q+1)*(*dhdr.get_spectral_va().c_cf)(lz, k, j, i);
		    }
		}
	}
	resu.set_dzpuis(dzp_resu) ;
	if (resu.set_spectral_va().c != 0x0) 
	    delete resu.set_spectral_va().c ;
	resu.set_spectral_va().c = 0x0 ;
    } //End of trace != 0    
    return resu ;
}


			//--------------//
			//   tilde(C)   //
			//--------------//
			
			
const Scalar& Sym_tensor::compute_tilde_C(bool output_ylm, Param* par) const {

  if (p_tilde_c == 0x0) {   // a new computation is necessary
	
    // All this has a meaning only for spherical components:
    assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    int dzp = operator()(1,1).get_dzpuis() ; 
    int dzp_resu = ((dzp == 0) ? 2 : dzp+1) ;

    p_tilde_c = new Scalar(*mp) ;
    p_tilde_c->set_etat_qcq() ;

    Scalar source_eta = operator()(1,2) ;
    source_eta.div_tant() ;
    source_eta += operator()(1,2).dsdt() + operator()(1,3).stdsdp() ;
    
    // Resolution of the angular Poisson equation for eta / r
    // ------------------------------------------------------
    Scalar tilde_eta(*mp) ;
    tilde_eta = 0 ;

    if (dynamic_cast<const Map_af*>(mp) != 0x0) {
	tilde_eta = source_eta.poisson_angu() ; 
    }
    else {
	assert (par != 0x0) ;
	mp->poisson_angu(source_eta, *par, tilde_eta) ;
    }
 
    Scalar dwdr = www().dsdr() ;
    dwdr.set_spectral_va().ylm() ;
    Scalar wsr = www() ;
    wsr.div_r_dzpuis(dzp_resu) ;
    wsr.set_spectral_va().ylm() ;
    Scalar etasr2 = tilde_eta ;
    etasr2.div_r_dzpuis(dzp_resu) ;
    etasr2.set_spectral_va().ylm() ;
    Scalar dtdr = ttt().dsdr() ;
    dtdr.set_spectral_va().ylm() ;
    Scalar tsr = ttt() ;
    tsr.div_r_dzpuis(dzp_resu) ;
    tsr.set_spectral_va().ylm() ;
    Scalar hrrsr = operator()(1,1) ;
    hrrsr.div_r_dzpuis(dzp_resu) ;
    hrrsr.set_spectral_va().ylm() ;

    int nz = mp->get_mg()->get_nzone() ;

    Base_val base(nz) ;

    if (wsr.get_etat() != ETATZERO) {
	base = wsr.get_spectral_base() ;
    }
    else {
	if (etasr2.get_etat() != ETATZERO) {
	    base = etasr2.get_spectral_base() ;
	}
	else {
	    if (tsr.get_etat() != ETATZERO) {
		base = tsr.get_spectral_base() ;
	    }
	    else {
		if (hrrsr.get_etat() != ETATZERO) {
		    base = hrrsr.get_spectral_base() ;
		}
		else { //Everything is zero...
		    p_tilde_c->set_etat_zero() ;
		    return *p_tilde_c ;
		}
	    }
	}
    }
    
    p_tilde_c->set_spectral_base(base) ;
    p_tilde_c->set_spectral_va().set_etat_cf_qcq() ;
    p_tilde_c->set_spectral_va().c_cf->annule_hard() ;

    if (dwdr.get_etat() == ETATZERO) dwdr.annule_hard() ;
    if (wsr.get_etat() == ETATZERO) wsr.annule_hard() ;
    if (etasr2.get_etat() == ETATZERO) etasr2.annule_hard() ;
    if (dtdr.get_etat() == ETATZERO) dtdr.annule_hard() ;
    if (tsr.get_etat() == ETATZERO) tsr.annule_hard() ;
    if (hrrsr.get_etat() == ETATZERO) hrrsr.annule_hard() ;

    int m_q, l_q, base_r ;
    for (int lz=0; lz<nz; lz++) {
	int np = mp->get_mg()->get_np(lz) ;
	int nt = mp->get_mg()->get_nt(lz) ;
	int nr = mp->get_mg()->get_nr(lz) ;
	for (int k=0; k<np+1; k++)
	    for (int j=0; j<nt; j++) {
		base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
		if ( (nullite_plm(j, nt, k, np, base) == 1) && (l_q > 1))
		{
		    for (int i=0; i<nr; i++) 
			p_tilde_c->set_spectral_va().c_cf->set(lz, k, j, i)
			    = - (l_q - 1)*(*dwdr.get_spectral_va().c_cf)(lz, k, j, i)
			    +  (l_q + 1)*(l_q - 1)*(*wsr.get_spectral_va().c_cf)(lz, k, j, i)
			    - 2*(*etasr2.get_spectral_va().c_cf)(lz, k, j, i)
	    + 0.5*double(l_q-1)/double(l_q)*(*tsr.get_spectral_va().c_cf)(lz, k, j, i)
	    - 0.5/double(l_q)*(*dtdr.get_spectral_va().c_cf)(lz, k, j, i)
	    + 1./double(l_q)*(*hrrsr.get_spectral_va().c_cf)(lz, k, j, i) ;
		}
	    }
    }
    p_tilde_c->set_dzpuis(dzp_resu) ;
  } //End of p_tilde_c == 0x0    

  if (output_ylm) p_tilde_c->set_spectral_va().ylm() ;
  else  p_tilde_c->set_spectral_va().ylm_i() ;

  return *p_tilde_c ; 

}
}
