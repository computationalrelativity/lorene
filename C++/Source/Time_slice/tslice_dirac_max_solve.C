/*
 *  Methods of class Tslice_dirac_max for solving Einstein equations
 *
 *    (see file time_slice.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Eric Gourgoulhon & Jerome Novak
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

char tslice_dirac_max_solve_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/04/30 10:52:14  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdlib.h>
#include <assert.h>

// Lorene headers
#include "time_slice.h"
#include "metric.h"
#include "evolution.h"
#include "unites.h"

                    //--------------------------//
                    //      Equation for N      //
                    //--------------------------//

Scalar Tslice_dirac_max::solve_n(const Scalar* p_ener_dens,
                                 const Scalar* p_trace_stress) const {

    using namespace Unites ;

    const Map& map = nn().get_mp() ; 
    Scalar ener_dens(map) ; 
    if (p_ener_dens != 0x0) ener_dens = *(p_ener_dens) ; 
    else ener_dens.set_etat_zero() ; 
    
    Scalar trace_stress(map) ; 
    if (p_trace_stress != 0x0) trace_stress = *(p_trace_stress) ; 
    else trace_stress.set_etat_zero() ; 
    
    // Source for N 
    // ------------
        
    const Vector& dln_psi = ln_psi().derive_cov(ff) ; // D_i ln(Psi)
    const Vector& dnn = nn().derive_cov(ff) ;         // D_i N

    Sym_tensor taa = aa().up_down(tgam()) ; 
        
    Scalar aa_quad = contract(taa, 0, 1, aa(), 0, 1) ; 
        
    Scalar source_nn = psi4()* nn()*( qpig* (ener_dens + trace_stress) 
                                + aa_quad )
                      - 2.* contract(dln_psi, 0, nn().derive_con(tgam()), 0) ; 
         
    Scalar tmp = - contract( hh(), 0, 1, dnn.derive_cov(ff), 0, 1 ) ;
        
    tmp.inc_dzpuis() ; // dzpuis: 3 -> 4
        
    source_nn += tmp ;
    
    // Resolution of the Poisson equation for the lapse
    // ------------------------------------------------
        
    Scalar nn_new = source_nn.poisson() + 1. ; 

    // Test:
    maxabs(nn_new.laplacian() - source_nn,
                "Absolute error in the resolution of the equation for N") ;  

    return nn_new ; 

}

                
                    //--------------------------//
                    //      Equation for Q      //
                    //--------------------------//


Scalar Tslice_dirac_max::solve_q(const Scalar* p_trace_stress) const {

    using namespace Unites ;

    Scalar trace_stress(qq().get_mp()) ; 
    if (p_trace_stress != 0x0) trace_stress = *(p_trace_stress) ; 
    else trace_stress.set_etat_zero() ; 
    
    // Source for Q
    // ------------
        
    const Vector& dqq = qq().derive_cov(ff) ;           // D_i Q
    const Vector& dln_psi = ln_psi().derive_cov(ff) ;   // D_i ln(Psi)
    const Tensor_sym& dhh = hh().derive_cov(ff) ;       // D_k h^{ij}
    const Tensor_sym& dtgam = tgam().cov().derive_cov(ff) ;    
                                                    // D_k {\tilde \gamma}_{ij}

    Sym_tensor taa = aa().up_down(tgam()) ; 
        
    Scalar aa_quad = contract(taa, 0, 1, aa(), 0, 1) ; 
        
    Scalar source_qq = psi4() * qq() * ( qpig* trace_stress + 0.75* aa_quad ) ;
        
    Scalar tmp = contract( hh(), 0, 1, dqq.derive_cov(ff), 0, 1 ) ;             
    tmp.inc_dzpuis() ; 
        
    source_qq -= tmp ;  
                        
    tmp = 0.0625 * contract( dhh, 0, 1, dtgam, 0, 1 ).trace(tgam()) 
          - 0.125 * contract( dhh, 0, 1, dtgam, 0, 2 ).trace(tgam()) 
          + 2.* contract( dln_psi, 0, ln_psi().derive_con(tgam()), 0) ; 
     
    source_qq += psi()*psi() * ( nn() * tmp 
                + 2*contract(dln_psi, 0, nn().derive_con(tgam()), 0) ) ; 
                             
               
    // Resolution of the Poisson equation for Q
    // -----------------------------------------
        
    Scalar qq_new = source_qq.poisson() + 1. ; 

    // Test:
    maxabs(qq_new.laplacian() - source_qq,
                "Absolute error in the resolution of the equation for Q") ;  

    return qq_new ; 

}
                


                    //--------------------------//
                    //      Equation for beta   //
                    //--------------------------//


Vector Tslice_dirac_max::solve_beta(const Vector* p_mom_dens) const {

    using namespace Unites ;

    Vector mom_dens(beta().get_mp(), CON, beta().get_triad()) ; 
    if (p_mom_dens != 0x0) mom_dens = *(p_mom_dens) ; 
    else mom_dens.set_etat_zero() ; 

    // Source for beta
    // ---------------
        
    const Vector& dln_psi = ln_psi().derive_cov(ff) ; // D_i ln(Psi)
    const Vector& dnn = nn().derive_cov(ff) ;         // D_i N
 
    Vector source_beta = 2.* contract(aa(), 1, 
                                   dnn - 6.*nn() * dln_psi, 0) ;
                
    source_beta += 2.* nn() * ( 2.*qpig* psi4() * mom_dens 
                        - contract(tgam().connect().get_delta(), 1, 2, 
                                   aa(), 0, 1) ) ;
            
    Vector vtmp = contract(hh(), 0, 1, 
                           beta().derive_cov(ff).derive_cov(ff), 1, 2)
                + 0.3333333333333333*
                  contract(hh(), 1, beta().divergence(ff).derive_cov(ff), 0) ; 
    vtmp.inc_dzpuis() ; // dzpuis: 3 -> 4
                    
    source_beta -= vtmp ; 


    // Resolution of the vector Poisson equation 
    //------------------------------------------
    
    int method = 1 ;  // method used to solve the vector Poisson equation
    
    Vector beta_new = source_beta.poisson(0.3333333333333333, method) ; 
        
    // Test:
    Vector test_beta = (beta_new.derive_con(ff)).divergence(ff)
            +  0.3333333333333333 * (beta_new.divergence(ff)).derive_con(ff) ;
    test_beta.inc_dzpuis() ;  
    maxabs(test_beta - source_beta,
                "Absolute error in the resolution of the equation for beta") ; 

    return beta_new ; 

}
                
