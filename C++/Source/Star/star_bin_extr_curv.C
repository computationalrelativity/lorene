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

// Headers Lorene
#include "star.h"

void Star_bin::extrinsic_curvature(){
    
    // Computation of K^{ij}
    // --------------------

    aa_auto = beta_auto.ope_killing_conf(gtilde) ;
    
    aa_auto = 0.5 * aa_auto / nn ;   
     
    // Computation of K_{ij} K^{ij}
    // ----------------------------
    
    Tensor aa_auto_dd = aa_auto.up_down(gtilde) ;
  
    aa_quad_auto = contract(aa_auto_dd, 0, 1, aa_auto, 0, 1, true) ; 
      
}
