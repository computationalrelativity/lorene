/*
 * Methods Star_bin::update_metric_der_comp
 *
 * (see file star.h for documentation)
 *
 */

/*
 *   Copyright (c) 2004 Francois Limousin
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


char star_bin_upmetr_der_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.12  2005/02/24 16:07:23  f_limousin
 * Improve the computation of dlogn, dlnq and dlnpsi. We import the part
 * coming from the companion and add the auto part.
 *
 * Revision 1.11  2005/02/18 13:14:18  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
 * Revision 1.10  2005/02/17 17:34:28  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.9  2004/06/22 12:52:47  f_limousin
 * Change qq, qq_auto and qq_comp to beta, beta_auto and beta_comp.
 *
 * Revision 1.8  2004/06/07 16:25:14  f_limousin
 * Minor modif.
 *
 * Revision 1.7  2004/04/08 16:33:32  f_limousin
 * The new variable is ln(Q) instead of Q=psi^2*N. It improves the
 * convergence of the code.
 *
 * Revision 1.6  2004/03/23 10:00:09  f_limousin
 * We now make the derivation with respect to the metric tilde
 * instead of the flat metric for the computation of dshift_comp.
 *
 * Revision 1.5  2004/02/27 09:53:14  f_limousin
 * Correction of an error on the computation of kcar_comp.
 *
 * Revision 1.4  2004/02/18 18:47:01  e_gourgoulhon
 * divshift_comp now computed via Tensor::divergence, the
 * method Tensor::scontract having disappeared.
 *
 * Revision 1.3  2004/01/20 15:20:23  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "star.h"
#include "utilitaires.h"
#include "graphique.h"

void Star_bin::update_metric_der_comp(const Star_bin& ) {
  

  // Derivatives of metric coefficients
  // ----------------------------------
  
    // dlogn
    //--------
    
    Vector dlogn_comp (mp, COV, mp.get_bvect_cart()) ;
    dlogn_comp.set_etat_qcq() ;
    Vector auxi (comp.logn_auto.derive_cov(comp.flat)) ;
    auxi.dec_dzpuis(2) ;
    auxi.change_triad(auxi.get_mp().get_bvect_cart()) ;
    auxi.change_triad(mp.get_bvect_cart()) ;
    assert ( *(auxi.get_triad()) == *(dlogn_comp.get_triad())) ;

    dlogn_comp.set(1).import_symy(auxi(1)) ;
    dlogn_comp.set(2).import_asymy(auxi(2)) ;
    dlogn_comp.set(3).import_symy(auxi(3)) ;
    dlogn_comp.std_spectral_base() ;
    dlogn_comp.inc_dzpuis(2) ;
    dlogn_comp.change_triad(mp.get_bvect_spher()) ;
  
    dlogn = logn_auto.derive_cov(flat) + dlogn_comp ;

    // dlnpsi
    //--------

    Scalar lnpsi_auto  = 0.5 * (lnq_auto - logn_auto) ;
    
    Vector dlnpsi_comp (mp, COV, mp.get_bvect_cart()) ;
    dlnpsi_comp.set_etat_qcq() ;
    auxi = 0.5*(comp.lnq_auto - comp.logn_auto).derive_cov(comp.flat) ;
    auxi.dec_dzpuis(2) ;
    auxi.change_triad(auxi.get_mp().get_bvect_cart()) ;
    auxi.change_triad(mp.get_bvect_cart()) ;
    assert ( *(auxi.get_triad()) == *(dlnpsi_comp.get_triad())) ;
    
    dlnpsi_comp.set(1).import_symy(auxi(1)) ;
    dlnpsi_comp.set(2).import_asymy(auxi(2)) ;
    dlnpsi_comp.set(3).import_symy(auxi(3)) ;
    dlnpsi_comp.std_spectral_base() ;
    dlnpsi_comp.inc_dzpuis(2) ;
    dlnpsi_comp.change_triad(mp.get_bvect_spher()) ;
    
    dlnpsi = lnpsi_auto.derive_cov(flat) + dlnpsi_comp ;
    
    // dlnq
    //-------

    Vector dlnq_comp (mp, COV, mp.get_bvect_cart()) ;
    dlnq_comp.set_etat_qcq() ;
    auxi = comp.lnq_auto.derive_cov(comp.flat) ;
    auxi.dec_dzpuis(2) ;
    auxi.change_triad(auxi.get_mp().get_bvect_cart()) ;
    auxi.change_triad(mp.get_bvect_cart()) ;
    assert ( *(auxi.get_triad()) == *(dlnq_comp.get_triad())) ;
    
    dlnq_comp.set(1).import_symy(auxi(1)) ;
    dlnq_comp.set(2).import_asymy(auxi(2)) ;
    dlnq_comp.set(3).import_symy(auxi(3)) ;
    dlnq_comp.std_spectral_base() ;
    dlnq_comp.inc_dzpuis(2) ;
    dlnq_comp.change_triad(mp.get_bvect_spher()) ;
    
    dlnq = lnq_auto.derive_cov(flat) + dlnq_comp ;

    // New value of hh_auto and hh_comp
    // ----------------------------------

    // The old hij_auto and hij_comp are TT but hij is not.
    hh_auto = hh_auto + (hh - hh_auto - hh_comp) * decouple ;
    hh_comp = hh_comp + (hh - hh_auto - hh_comp) * (1-decouple) ;


    // Computation of A^{ij}_comp
    // --------------------------
    
    aa_comp = beta_comp.ope_killing_conf(gtilde) ;
      
    aa_comp = 0.5 * aa_comp / nn ;
      
    // Computation of aa_quad_comp
    // ------------------------

    Tensor aa_auto_dd = aa_auto.down(0, gtilde).down(1, gtilde) ;

    aa_quad_comp = contract(aa_auto_dd, 0, 1, aa_comp, 0, 1, true) ; 

}      

