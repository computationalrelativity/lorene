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
  
    dcov_logn = logn.derive_cov(flat) ;
    dcon_logn = logn.derive_con(flat) ;

    Scalar lnpsi = 0.5 * (lnq - logn) ;
    dcov_lnpsi = lnpsi.derive_cov(flat) ;
    dcon_lnpsi = lnpsi.derive_con(flat) ;

    // New value of hh_auto and hh_comp
    // ----------------------------------

    // The old hij_auto and hij_comp are TT but hij is not.
    hh_auto = hh_auto + (hh - hh_auto - hh_comp) * decouple ;
    hh_comp = hh_comp + (hh - hh_auto - hh_comp) * (1-decouple) ;


    // Computation of aa_comp
    // ------------------------
    
    // Gradient tilde (partial derivatives with respect to
    //           the Spherical coordinates of the mapping)
    // D~^j beta^i
    
    const Tensor& dbeta_comp = beta_comp.derive_con(gtilde) ;
    
    // Trace of D~_j beta^i  :
    Scalar divbeta_comp = beta_comp.divergence(gtilde) ;
    
    // Computation of A^{ij}
    // -------------------------
      
    for (int i=1; i<=3; i++) 
	for (int j=i; j<=3; j++) {

	  aa_comp.set(i, j) = dbeta_comp(i, j) + dbeta_comp(j, i) - 
	    double(2) /double(3) * divbeta_comp * (gtilde.con())(i,j) ; 
	}
      
      aa_comp = 0.5 * aa_comp / nn ;
      
      // Computation of aa_quad_comp
      // ------------------------

      Tensor aa_auto_cov = aa_auto.down(0, gtilde).down(1, gtilde) ;

      aa_quad_comp = contract(aa_auto_cov, 0, 1, aa_comp, 0, 1, true) ; 

}      

