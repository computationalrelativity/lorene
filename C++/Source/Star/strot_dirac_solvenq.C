/*
 *  Solution of the two scalar Poisson equations for rotating stars 
 *  in Dirac gauge.
 *
 *    (see file star_rot_Dirac.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Lap-Ming Lin & Jerome Novak
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

char strot_dirac_solvenq_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:53:40  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:18  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2005/02/17 17:31:29  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.2  2005/02/02 09:23:20  lm_lin
 *
 * Correct a mistake in the calculation of log(N).
 *
 *
 * $Header$
 *
 */

// C headers
#include <cmath>
#include <cassert>

// Lorene headers
#include "star_rot_dirac.h"
#include "unites.h"

namespace Lorene {
void Star_rot_Dirac::solve_logn_f(Scalar& ln_f_new) const {

    using namespace Unites ;
    
  //================================================
  // Source for log(n) containing matter field terms
  //================================================
    
    Scalar source_logn = qpig* psi4* (ener_euler + s_euler) ;
    
    ln_f_new = source_logn.poisson() ;

}

void Star_rot_Dirac::solve_logn_q(Scalar& ln_q_new) const {

    const Metric_flat& mets = mp.flat_met_spher() ;
    const Base_vect_spher& bspher = mp.get_bvect_spher() ;
    
    const Vector& dln_psi = ln_psi.derive_cov(mets) ; // D_i ln(Psi)
    const Vector& dln = logn.derive_cov(mets) ;         // D_i ln(N)



  //=============================================
  // Source for log(n) containing quadratic terms
  //=============================================

    Scalar source_logn = psi4* aa_quad 
      - contract(dln.up_down(mets), 0, dln, 0)
      - 2.* contract(dln_psi, 0, logn.derive_con(tgamma), 0) ;

         
    Tensor_sym tmp(mp, 2, COV, bspher, 0, 1) ;
    tmp = dln.derive_cov(mets) ;
    tmp.inc_dzpuis() ; //dzpuis 3 -> 4

    source_logn -= contract( hh, 0, 1, tmp+dln*dln, 0, 1 ) ;
        
    ln_q_new = source_logn.poisson() ;
    
}

void Star_rot_Dirac::solve_qqq(Scalar& q_new) const {

    using namespace Unites ;

    const Metric_flat& mets = mp.flat_met_spher() ;

    // Source for Q
    // ------------
        
    const Vector& dln_psi = ln_psi.derive_cov(mets) ; // D_i ln(Psi)
    const Vector& dqq = qqq.derive_cov(mets) ;           // D_i Q
    const Tensor_sym& dhh = hh.derive_cov(mets) ;       // D_k h^{ij}
    const Tensor_sym& dtgam = tgamma.cov().derive_cov(mets) ;    
                                                    // D_k {\tilde \gamma}_{ij}
    Scalar source_qq = psi4 * qqq * ( qpig* s_euler + 0.75* aa_quad ) ;
        
    Scalar tmp = contract( hh, 0, 1, dqq.derive_cov(mets), 0, 1 ) ; 
    tmp.inc_dzpuis() ; 
        
    source_qq -= tmp ;  
                        

    // tmp = \tilde{R}_{*}/4 + 2 D_i ln(Psi) \tilde{D}^i ln(Psi)
    tmp = 0.0625 * contract( dhh, 0, 1, dtgam, 0, 1 ).trace(tgamma) 
          - 0.125 * contract( dhh, 0, 1, dtgam, 0, 2 ).trace(tgamma) 
          + 2.* contract( dln_psi, 0, ln_psi.derive_con(tgamma), 0) ; 
     
    source_qq += psi2 * ( nn * tmp 
                + 2*contract(dln_psi, 0, nn.derive_con(tgamma), 0) ) ; 
                             
               
    q_new = source_qq.poisson() + 1. ; 

    if (q_new.get_etat() == ETATUN) q_new.std_spectral_base() ; 

}
}
