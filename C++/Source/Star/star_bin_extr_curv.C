/*
 * Method of class Star_bin to compute the extrinsic curvature tensor
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


char star_bin_extr_curv_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.10  2014/10/13 08:53:38  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.9  2014/10/06 15:13:16  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2005/09/13 19:38:31  f_limousin
 * Reintroduction of the resolution of the equations in cartesian coordinates.
 *
 * Revision 1.7  2005/02/24 16:04:44  f_limousin
 * Change the name of some variables (for instance dcov_logn --> dlogn).
 *
 * Revision 1.6  2005/02/17 17:33:38  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.5  2004/05/25 14:19:01  f_limousin
 * Correction of an error : kcar_comp was computed instead
 * ok kcar_auto.
 *
 * Revision 1.4  2004/03/23 09:57:57  f_limousin
 * We now make the derivation with respect to the metric tilde
 * instead of the flat metric for the computation of dshift.
 *
 * Revision 1.3  2004/02/27 09:52:41  f_limousin
 * Correction of an error on the computation of kcar_auto.
 *
 * Revision 1.2  2004/01/20 15:18:00  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// C headers
#include <cmath>

// Headers Lorene
#include "star.h"

namespace Lorene {
void Star_bin::extrinsic_curvature(double om){
    
  // Construction of Omega d/dphi
  // ----------------------------
  
  const Mg3d* mg = mp.get_mg() ; 
  int nz = mg->get_nzone() ;	    // total number of domains
  Vector omdsdp (mp, CON, mp.get_bvect_cart()) ;
  Scalar yya (mp) ;
  yya = mp.ya ;
  Scalar xxa (mp) ;
  xxa = mp.xa ;
  
  if (fabs(mp.get_rot_phi()) < 1e-10){ 
    omdsdp.set(1) = - om * yya ;
    omdsdp.set(2) = om * xxa ;
    omdsdp.set(3).annule_hard() ;
  }
  else{
    omdsdp.set(1) = om * yya ;
    omdsdp.set(2) = - om * xxa ;
    omdsdp.set(3).annule_hard() ;
  }
  
  omdsdp.set(1).set_spectral_va()
    .set_base(*(mp.get_mg()->std_base_vect_cart()[0])) ;
  omdsdp.set(2).set_spectral_va()
    .set_base(*(mp.get_mg()->std_base_vect_cart()[1])) ;
  omdsdp.set(3).set_spectral_va()
    .set_base(*(mp.get_mg()->std_base_vect_cart()[2])) ;
  
  omdsdp.annule_domain(nz-1) ;


  // Gradient tilde (with respect to the cartesian coordinates
  //           of the mapping)
  // D~_j beta^i 
  
  const Tensor& dbeta = beta_auto.derive_con(gtilde) ;
                           
  // Trace of D~_j beta^i : 
  Scalar div_beta = beta_auto.divergence(gtilde) ;

  // Computation of K^{ij}
  // See Eq (49) from Gourgoulhon et al. (2001)
  // ------------------------------------------
  
  for (int i=1; i<=3; i++) 
    for (int j=1; j<=i; j++) {
      tkij_auto.set(i, j) = dbeta(i, j) + dbeta(j, i) - 
	double(2) /double(3) * div_beta * (gtilde.con())(i,j) ; 
    }
  

  // Addition (or not !) of u^{ij}
  tkij_auto = tkij_auto - 0*hij_auto.derive_lie(omdsdp) ;

  tkij_auto = 0.5 * tkij_auto / nn ;   
  
  // Computation of K_{ij} K^{ij}
  // ----------------------------
  
  Sym_tensor tkij_auto_cov = tkij_auto.up_down(gtilde) ;
  
  kcar_auto = contract(tkij_auto_cov, 0, 1, tkij_auto, 0, 1,true) ; 
      
  // The derived quantities are obsolete
  // -----------------------------------
  
  del_deriv() ;
}
}
