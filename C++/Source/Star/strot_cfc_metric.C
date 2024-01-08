/*
 *  Solution of the metric Poisson equations for rotating stars in CFC. 
 *
 *    (see file star_rot_cfc.h for documentation).
 *
 */

/*
 *   Copyright (c) 2024 Jerome Novak
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

 
// Lorene headers
#include "star_rot_cfc.h"
#include "unites.h"

namespace Lorene {

  //===================================================
  // Equation for log(n) containing matter field terms
  //===================================================
   
  void Star_rot_CFC::solve_logn_f(Scalar& ln_f_new) const {

    using namespace Unites ;
    
    Scalar source_logn = qpig* psi4* (ener_euler + s_euler) ;
    
    ln_f_new = source_logn.poisson() ;

}

  //================================================
  // Equation for log(n) containing quadratic terms
  //================================================

void Star_rot_CFC::solve_logn_q(Scalar& ln_q_new) const {

  Scalar ln_psi = log(psi) ;
  ln_psi.std_spectral_base() ;
  const Vector& dln_psi = ln_psi.derive_cov(flat) ; // D_i ln(Psi)
  const Vector& dln = logn.derive_cov(flat) ;         // D_i ln(N)

  Scalar source_logn = hatA_quad/(psi4*psi4) 
    - contract(dln.up_down(flat), 0, dln, 0)
    - 2.* contract(dln_psi, 0, logn.derive_con(flat), 0) ;
  
  ln_q_new = source_logn.poisson() ;
  
}

  //========================================
  // Equation for the conformal factor psi
  //========================================
  
  void Star_rot_CFC::solve_psi(Scalar& psi_new) {
    
    using namespace Unites ;
    
    // Source for the equation for X^i
    //---------------------------------
    
    Vector sou_Xi = 2*qpig*psi4*psi4*psi2*j_euler ;
    double lambda = 1./3. ;
    Vector Xi_new = sou_Xi.poisson(lambda) ;
    //    Xi_new.std_spectral_base() ;
    
    // Computation of \hat{A}^{ij} and its norm
    //------------------------------------------
    hatA = Xi_new.ope_killing_conf(flat) ;
    hatA_quad = contract(hatA, 0, 1, hatA.up_down(flat), 0, 1) ;
    
    
    // Source for conformal factor psi
    // --------------------------------
    Scalar source_psi = -0.5*qpig*psi4*psi*ener_euler
      - 0.125*hatA_quad/(psi4*psi2*psi) ;
    source_psi.std_spectral_base() ;
    
    psi_new = source_psi.poisson() + 1. ; 
    
    if (psi_new.get_etat() == ETATUN) psi_new.std_spectral_base() ; 
    
  }

  //========================
  // Equation for the shift
  //========================
  
  void Star_rot_CFC::solve_shift(Vector& beta_new) const {
    
    using namespace Unites ;
    
    double lambda = 1./3. ;

    Sym_tensor hatA_tmp = 2*nn*hatA/(psi4*psi2) ;
    Vector sou_shift = hatA_tmp.divergence(flat) ;
    sou_shift.inc_dzpuis(1) ;

    beta_new = sou_shift.poisson(lambda) ;
    beta_new.set(1) = 0 ; //## these components are null
    beta_new.set(2) = 0 ; //## in axial symmetry
   
  }

  //=========================================
  // Update of the metric derived quantities
  //=========================================
  void Star_rot_CFC::update_metric() {

    nn = exp(logn) ;
    nn.std_spectral_base() ;

    psi2 = psi*psi ;
    psi4 = psi2*psi2 ;

    gamma = psi4*flat.cov() ;

    del_deriv() ;
  }
  
}
