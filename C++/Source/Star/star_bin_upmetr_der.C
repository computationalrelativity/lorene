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

void Star_bin::update_metric_der_comp(const Star_bin& comp) {
  

  // Derivatives of metric coefficients
  // ----------------------------------
  
    dcov_logn = logn.derive_cov(flat) ;
    dcon_logn = logn.derive_con(flat) ;

    Scalar lnpsi = log (psi4) / 4. ;
    lnpsi.std_spectral_base() ;
    dcov_lnpsi = lnpsi.derive_cov(flat) ;
    dcon_lnpsi = lnpsi.derive_con(flat) ;

 
  // Computation of tkij_comp
  // ------------------------
    
    // Gradient tilde (partial derivatives with respect to
    //           the Spherical coordinates of the mapping)
    // D~^j beta^i
    
    const Tensor& dshift_comp = shift_comp.derive_con(flat) ;
    
    // Trace of D~_j beta^i  :
    Scalar divshift_comp = shift_comp.divergence(flat) ;
    
    // Computation of K^{ij}
    // -------------------------
      
    for (int i=1; i<=3; i++) 
	for (int j=i; j<=3; j++) {

	  tkij_comp.set(i, j) = dshift_comp(i, j) + dshift_comp(j, i) - 
	    double(2) /double(3) * divshift_comp * (gtilde.con())(i,j) ; 
	}
      
      tkij_comp = 0.5 * tkij_comp / nnn ;
      tkij_comp.std_spectral_base() ;
      
      // Computation of kcar_comp
      // ------------------------

      Tensor tkij_auto_cov = tkij_auto.down(0, gtilde).down(1, gtilde) ;

      kcar_comp = contract(tkij_auto_cov, 0, 1, tkij_comp, 0, 1, true) ; 
        
      kcar_comp.std_spectral_base() ;

}      

