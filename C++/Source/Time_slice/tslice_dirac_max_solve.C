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
 * Revision 1.17  2014/10/06 15:13:22  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.16  2010/10/20 07:58:10  j_novak
 * Better implementation of the explicit time-integration. Not fully-tested yet.
 *
 * Revision 1.15  2008/12/02 15:02:22  j_novak
 * Implementation of the new constrained formalism, following Cordero et al. 2009
 * paper. The evolution eqs. are solved as a first-order system. Not tested yet!
 *
 * Revision 1.14  2004/07/08 12:29:01  j_novak
 * use of new method Tensor::annule_extern_cn
 *
 * Revision 1.13  2004/06/30 08:02:40  j_novak
 * Added filtering in l of khi_new and mu_new. ki_source is forced to go to
 * zero at least as r^2.
 *
 * Revision 1.12  2004/06/17 07:07:11  e_gourgoulhon
 * Method solve_hij:
 *   -- replaced the attenuation of khi_source with tempo by a call to
 *      the new method Tensor::annule_extern_c2
 *   -- suppressed filtre_r on khi_source and khi_new
 *   -- added graphs of W^i and LW^{ij} (transverse decomp of S^ij).
 *
 * Revision 1.11  2004/06/15 09:43:36  j_novak
 * Attenuation of the source for khi in the last shell (temporary?).
 *
 * Revision 1.10  2004/06/14 20:47:31  e_gourgoulhon
 * Added argument method_poisson to method solve_hij.
 *
 * Revision 1.9  2004/06/03 10:02:45  j_novak
 * Some filtering is done on source_khi and khi_new.
 *
 * Revision 1.8  2004/05/24 21:00:44  e_gourgoulhon
 * Method solve_hij: khi and mu.smooth_decay(2,2) --> smooth_decay(2,1) ;
 *   added exponential_decay(khi) and exponential_decay(mu) after the
 *   call to smooth_decay. Method exponential_decay is provisory defined
 *   in this file.
 *
 * Revision 1.7  2004/05/17 19:56:25  e_gourgoulhon
 * -- Method solve_beta: added argument method
 * -- Method solve_hij: added argument graph_device
 *
 * Revision 1.6  2004/05/12 15:24:20  e_gourgoulhon
 * Reorganized the #include 's, taking into account that
 * time_slice.h contains now an #include "metric.h".
 *
 * Revision 1.5  2004/05/05 14:47:05  e_gourgoulhon
 * Modified text and graphical outputs.
 *
 * Revision 1.4  2004/05/03 15:06:27  e_gourgoulhon
 * Added matter source in solve_hij.
 *
 * Revision 1.3  2004/05/03 14:50:00  e_gourgoulhon
 * Finished the implementation of method solve_hij.
 *
 * Revision 1.2  2004/04/30 14:36:15  j_novak
 * Added the method Tslice_dirac_max::solve_hij(...)
 * NOT READY YET!!!
 *
 * Revision 1.1  2004/04/30 10:52:14  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

// C headers
#include <cassert>

// Lorene headers
#include "time_slice.h"
#include "unites.h"
#include "graphique.h"
#include "proto.h"

                    //----------------------------//
                    //    Equation for N\Psi     //
                    //----------------------------//

Scalar Tslice_dirac_max::solve_npsi(const Scalar* p_ener_dens,
                                 const Scalar* p_trace_stress) const {

    using namespace Unites ;

    const Map& map = npsi().get_mp() ;
    Scalar ener_dens(map) ; 
    if (p_ener_dens != 0x0) ener_dens = *(p_ener_dens) ; 
    else ener_dens.set_etat_zero() ; 
    
    Scalar trace_stress(map) ; 
    if (p_trace_stress != 0x0) trace_stress = *(p_trace_stress) ; 
    else trace_stress.set_etat_zero() ; 
    
    // Source for N\Psi 
    // ----------------
        
    const Vector& dnpsi = npsi().derive_cov(ff) ; // D_i N\Psi
    const Tensor_sym& dhh = hh().derive_cov(ff) ;       // D_k h^{ij}
    const Tensor_sym& dtgam = tgam().cov().derive_cov(ff) ;    
                                                    // D_k {\tilde \gamma}_{ij}
    Sym_tensor taa = aa().up_down(tgam()) ; 
        
    Scalar aa_quad = contract(taa, 0, 1, aa(), 0, 1) ; 
    Scalar tildeR = 0.25 * contract( dhh, 0, 1, dtgam, 0, 1 ).trace(tgam()) 
          - 0.5 * contract( dhh, 0, 1, dtgam, 0, 2 ).trace(tgam()) ; 
        
    Scalar source_npsi = -  contract( hh(), 0, 1, dnpsi.derive_cov(ff), 0, 1 ) ;
    source_npsi.inc_dzpuis() ;
    source_npsi += npsi()*( 0.5*qpig*(ener_dens + 2.*trace_stress )/( psi()*psi() ) 
			    + 0.875*aa_quad/( psi4()*psi4() ) + 0.125*tildeR ) ;
        
    // Resolution of the Poisson equation for N\Psi
    // --------------------------------------------
        
    Scalar npsi_new = source_npsi.poisson() + 1. ; 

    if (npsi_new.get_etat() == ETATUN) npsi_new.std_spectral_base() ; 

#ifndef NDEBUG
    // Test:
    maxabs(npsi_new.laplacian() - source_npsi,
                "Absolute error in the resolution of the equation for N") ;  
#endif
    return npsi_new ; 

}

                
                    //--------------------------//
                    //     Equation for \Psi    //
                    //--------------------------//


Scalar Tslice_dirac_max::solve_psi(const Scalar* p_ener_dens) const {

    using namespace Unites ;

    const Map& map = psi().get_mp() ;
    Scalar ener_dens(map) ; 
    if (p_ener_dens != 0x0) ener_dens = *(p_ener_dens) ; 
    else ener_dens.set_etat_zero() ; 
    
    // Source for \Psi
    // ---------------
        
    const Vector& dpsi = psi().derive_cov(ff) ;           // D_i Psi
    const Tensor_sym& dhh = hh().derive_cov(ff) ;       // D_k h^{ij}
    const Tensor_sym& dtgam = tgam().cov().derive_cov(ff) ;    
                                                    // D_k {\tilde \gamma}_{ij}

    Sym_tensor taa = hata().up_down(tgam()) ; 
        
    Scalar aa_quad = contract(taa, 0, 1, hata(), 0, 1) ; 
        
    Scalar tildeR = 0.25 * contract( dhh, 0, 1, dtgam, 0, 1 ).trace(tgam()) 
          - 0.5 * contract( dhh, 0, 1, dtgam, 0, 2 ).trace(tgam()) ; 

    Scalar source_psi = -contract( hh(), 0, 1, dpsi.derive_cov(ff), 0, 1 ) ;
    source_psi.inc_dzpuis() ;
    source_psi -= 0.5*qpig*ener_dens/psi()  
	+ 0.125*( aa_quad*pow(psi(), -7) - tildeR*psi() ) ; 
               
    // Resolution of the Poisson equation for Psi
    // ------------------------------------------
        
    Scalar psi_new = source_psi.poisson() + 1. ; 

    if (psi_new.get_etat() == ETATUN) psi_new.std_spectral_base() ; 

#ifndef NDEBUG
    // Test:
    maxabs(psi_new.laplacian() - source_psi,
                "Absolute error in the resolution of the equation for Psi") ;  
#endif

    return psi_new ; 
}
                


                    //--------------------------//
                    //      Equation for beta   //
                    //--------------------------//


Vector Tslice_dirac_max::solve_beta(int method) 
    const {

    // Source for beta
    // ---------------
    Vector source_beta = 
	- contract(hh(), 0, 1, beta().derive_cov(ff).derive_cov(ff), 1, 2)
	- 0.3333333333333333*contract(hh(), 1, beta().divergence(ff).derive_cov(ff), 0) ; 

    Sym_tensor sym_tmp = 2*nn()*aa() ;
    source_beta += sym_tmp.divergence(ff) ;
    source_beta.inc_dzpuis() ; // dzpuis: 3 -> 4

    // Resolution of the vector Poisson equation 
    //------------------------------------------
    
    Vector beta_new = source_beta.poisson(0.3333333333333333, ff, method) ; 
        
#ifndef NDEBUG
    // Test:
    Vector test_beta = (beta_new.derive_con(ff)).divergence(ff)
            +  0.3333333333333333 * (beta_new.divergence(ff)).derive_con(ff) ;
    test_beta.inc_dzpuis() ;  
    maxabs(test_beta - source_beta,
                "Absolute error in the resolution of the equation for beta") ; 
#endif
    return beta_new ; 

}
                

