/*
 * Method of class Star_bin_xcts to compute 
 * the extrinsic curvature tensor 
 * (see file star.h for documentation)
 */

/*
 *   Copyright (c) 2010 Michal Bejger
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:18:14  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:38  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:17  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2010/06/15 08:10:29  m_bejger
 * *** empty log message ***
 *
 * Revision 1.1  2010/05/04 07:51:05  m_bejger
 * Initial version
 *
 * $Header$
 *
 */

// C headers
#include <cmath>

// Headers Lorene
#include "star.h"

namespace Lorene {
void Star_bin_xcts::extrinsic_curvature() {
     
  // D~_j beta^i 
  const Tensor& dbeta = beta_auto.derive_con(flat) ;
                           
  // Trace of D~_j beta^i : 
  Scalar div_beta = beta_auto.divergence(flat) ;

  // Computation of \hat{A}^{ij}, Eq. 8.130 of arXiv:gr-qc/0703035 
  // -------------------------------------------------------------
  
    for (int i=1; i<=3; i++) 
    for (int j=1; j<=i; j++) {
		
      haij_auto.set(i, j) = dbeta(i, j) + dbeta(j, i) - 
		double(2) /double(3) * div_beta * (flat.con())(i,j) ; 
    }

  haij_auto = 0.5 * pow(Psi, 7.) * haij_auto / chi ; 
  //## for comparison: old formulation
  //haij_auto = 0.5 * haij_auto / nn ; 
  
  haij_auto.std_spectral_base() ; 
  
  // Computation of (\hat{A}_{ij}\hat{A}^{ij})_{auto}
  // ------------------------------------------------
  
  Sym_tensor haij_auto_cov = haij_auto.up_down(flat) ;
  
  hacar_auto = contract(haij_auto_cov, 0, 1, haij_auto, 0, 1, true) ; 
      
  // The derived quantities are obsolete
  // -----------------------------------
  
  del_deriv() ;

}
}
