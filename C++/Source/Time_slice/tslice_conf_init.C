/*
 *  Method of class Time_slice_conf to compute valid initial data
 *
 *    (see file time_slice.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon & Jerome Novak
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

char tslice_conf_init_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/04/07 07:58:21  e_gourgoulhon
 * Constructor as Minkowski slice: added call to std_spectral_base()
 * after setting the lapse to 1.
 *
 * Revision 1.1  2004/04/05 21:25:37  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>
#include <assert.h>

// Lorene headers
#include "time_slice.h"
#include "metric.h"
#include "evolution.h"
#include "unites.h"
#include "graphique.h"
#include "utilitaires.h"

void Time_slice_conf::initial_data_cts(const Sym_tensor& hh_in, 
                const Sym_tensor& uu, const Scalar& trk_in, 
                const Scalar& trk_point, const Scalar* p_ener_dens,
                const Vector* p_mom_dens, const Scalar* p_trace_stress) {

    using namespace Unites ;

    // Verifications
    assert(trk_in.check_dzpuis(2)) ; 
    assert(trk_point.check_dzpuis(4)) ; 

    // Initialisations
    double ttime = the_time[jtime] ; 
    
    hh_evol.update(hh_in, jtime, ttime) ; 

    // Reset of quantities depending on h^{ij}:
    if (p_tgamma != 0x0) {
        delete p_tgamma ;
        p_tgamma = 0x0 ; 
    } 
    if (p_hdirac != 0x0) {
        delete p_hdirac ; 
        p_hdirac = 0x0 ; 
    }
    if (p_gamma != 0x0) {
        delete p_gamma ; 
        p_gamma = 0x0 ;
    }
    gam_dd_evol.downdate(jtime) ; 
    gam_uu_evol.downdate(jtime) ; 
     
    trk_evol.update(trk_in, jtime, ttime) ; 

    // Reset of quantities depending on K:
    k_dd_evol.downdate(jtime) ; 
    k_uu_evol.downdate(jtime) ; 
   
    aa_evol.update(uu / (2.* nn()), jtime, ttime) ; 

    // Reset of quantities depending on A^{ij}:
    k_dd_evol.downdate(jtime) ; 
    k_uu_evol.downdate(jtime) ; 
    
    
    
    // Iteration
    int imax = 100 ; 
    double precis = 1.e-12 ; 
    const Map& map = hh_in.get_mp() ; 
    const Base_vect& triad = *(hh_in.get_triad()) ;

    Scalar ener_dens(map) ; 
    if (p_ener_dens != 0x0) ener_dens = *(p_ener_dens) ; 
    else ener_dens.set_etat_zero() ; 
    
    Vector mom_dens(map, CON, triad) ; 
    if (p_mom_dens != 0x0) mom_dens = *(p_mom_dens) ; 
    else mom_dens.set_etat_zero() ; 
    
    Scalar trace_stress(map) ; 
    if (p_trace_stress != 0x0) trace_stress = *(p_trace_stress) ; 
    else trace_stress.set_etat_zero() ; 
    
    Scalar tmp(map) ; 
    Scalar source_psi(map) ; 
    Scalar source_nn(map) ; 
    Vector source_beta(map, CON, triad) ; 
    
    for (int i=0; i<imax; i++) {
    
        //===============================================
        //  Computations of sources 
        //===============================================
    
        const Vector& dpsi = psi().derive_cov(ff) ;       // D_i Psi
        const Vector& dln_psi = ln_psi().derive_cov(ff) ; // D_i ln(Psi)
        const Vector& dnn = nn().derive_cov(ff) ;         // D_i N
        
        Sym_tensor taa = aa().up_down(tgam()) ; 
        
        Scalar aa_quad = contract(taa, 0, 1, aa(), 0, 1) ; 
        
        // Source for Psi 
        // --------------
        tmp = 0.125* psi() * tgam().ricci_scal() 
                - contract(hh(), 0, 1, dpsi.derive_cov(ff), 0, 1 ) ;
        tmp.inc_dzpuis() ; // dzpuis : 3 -> 4

        tmp -= contract(hdirac(), 0, dpsi, 0) ;  
                
        source_psi = tmp - psi()*psi4()* ( 0.5*qpig* ener_dens 
                        + 0.125* aa_quad 
                        - 8.33333333333333e-2* trk()*trk() ) ;  
                        
        // Source for N 
        // ------------
        
        source_nn = psi4()*( nn()*( qpig* (ener_dens + trace_stress) + aa_quad
                                    - 0.3333333333333333* trk()*trk() )
                             - trk_point ) 
                    - 2.* contract(dln_psi, 0, nn().derive_con(tgam()), 0)  
                    - contract(hdirac(), 0, dnn, 0) ; 
        
        tmp = psi4()* contract(beta(), 0, trk().derive_cov(ff), 0)
                - contract( hh(), 0, 1, dnn.derive_cov(ff), 0, 1 ) ;
        
        tmp.inc_dzpuis() ; // dzpuis: 3 -> 4
        
        source_nn += tmp ;


        // Source for beta 
        // ---------------

        source_beta = 2.* contract(aa(), 1, 
                                   dnn - 6.*nn() * dln_psi, 0) ;
                
        source_beta += 2.* nn() * ( 2.*qpig* psi4() * mom_dens 
            + 0.66666666666666666* trk().derive_con(tgam()) 
            - contract(tgam().connect().get_delta(), 1, 2, 
                                  aa(), 0, 1) ) ;
            
        Vector vtmp = contract(hh(), 0, 1, 
                           beta().derive_cov(ff).derive_cov(ff), 1, 2)
                + 0.3333333333333333*
                  contract(hh(), 1, beta().divergence(ff).derive_cov(ff), 0) 
                - hdirac().derive_lie(beta()) 
                + uu.divergence(ff) ; 
        vtmp.inc_dzpuis() ; // dzpuis: 3 -> 4
                    
        source_beta -= vtmp ; 
        
        source_beta += 0.66666666666666666* beta().divergence(ff) * hdirac() ;
        

        //=============================================
        // Resolution of elliptic equations
        //=============================================
        
        // Resolution of the Poisson equation for Psi
        // ------------------------------------------
        
        Scalar psi_jp1 = source_psi.poisson() + 1. ; 

        // Test:
        diffrel(psi_jp1.laplacian(), source_psi,
                "Relative error in the resolution of the equation for Psi") ;  

        des_meridian(psi_jp1, 0., 5., "Psi", 1) ; 

        // Resolution of the Poisson equation for the lapse
        // ------------------------------------------------
        
        Scalar nn_jp1 = source_nn.poisson() + 1. ; 

        // Test:
        diffrel(nn_jp1.laplacian(), source_nn,
                "Relative error in the resolution of the equation for N") ;  

        des_meridian(nn_jp1, 0., 5., "N", 2) ; 
        
        // Resolution of the vector Poisson equation for the shift
        //---------------------------------------------------------
        
        Vector beta_jp1 = source_beta.poisson(0.3333333333333333, 1) ; 
        
        des_meridian(beta_jp1(1), 0., 5., "\\gb\\ur\\d", 3) ; 
        des_meridian(beta_jp1(2), 0., 5., "\\gb\\u\\gh\\d", 4) ; 
        des_meridian(beta_jp1(3), 0., 5., "\\gb\\u\\gf\\d", 5) ; 
        
        // Test:
        Vector test_beta = (beta_jp1.derive_con(ff)).divergence(ff)
            +  0.3333333333333333 * (beta_jp1.divergence(ff)).derive_con(ff) ;
        test_beta.inc_dzpuis() ;  
        diffrel(test_beta, source_beta,
                "Relative error (L^1) in the resolution for beta") ; 
        diffrelmax(test_beta, source_beta,
                "Relative error (L^infty) in the resolution for beta") ; 

        //===========================================
        //      Convergence control
        //===========================================
    
        double diff_psi = max( diffrel(psi(), psi_jp1) ) ; 
        double diff_nn = max( diffrel(nn(), nn_jp1) ) ; 
        double diff_beta = max( diffrel(beta(), beta_jp1) ) ; 
        
        cout << "step = " << i << " :  diff_psi = " << diff_psi 
             << "  diff_nn = " << diff_nn
             << "  diff_beta = " << diff_beta << endl ; 
        if ( (diff_psi < precis) && (diff_nn < precis) && (diff_beta < precis) )
            break ; 

        //=============================================
        //      Updates for next step 
        //=============================================

        psi_evol.update(psi_jp1, jtime, ttime) ; 

        // Reset of quantities depending on Psi:
        qq_evol.downdate(jtime) ; 
        if (p_psi4 != 0x0) {
            delete p_psi4 ; 
            p_psi4 = 0x0 ; 
        }
        if (p_ln_psi != 0x0) {
            delete p_ln_psi ; 
            p_ln_psi = 0x0 ; 
        }
        if (p_gamma != 0x0) {
            delete p_gamma ; 
            p_gamma = 0x0 ; 
        }
        gam_dd_evol.downdate(jtime) ; 
        gam_uu_evol.downdate(jtime) ; 
     
        n_evol.update(nn_jp1, jtime, ttime) ; 

        // Reset of quantities depending on N:
        qq_evol.downdate(jtime) ; 

        beta_evol.update(beta_jp1, jtime, ttime) ; 

        // New value of A^{ij}:
        Sym_tensor aa_jp1 = ( beta().ope_killing_conf(tgam()) + uu ) 
                                / (2.* nn()) ; 
        
        //## Alternative formula:
        // Sym_tensor aa_jp1 = ( beta().ope_killing_conf(ff) 
        //                      - hh().derive_lie(beta())
        //                      - 0.6666666666666666 * beta.divergence() * hh()
        //                      + uu ) / (2.* nn()) ; 
                             
        aa_evol.update(aa_jp1, jtime, ttime) ; 

        // Reset of quantities depending on A^{ij}:
        k_dd_evol.downdate(jtime) ; 
        k_uu_evol.downdate(jtime) ; 

        // arrete() ; 

    }

} 













