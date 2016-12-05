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

 

/*
 * $Id$
 * $Log$
 * Revision 1.14  2016/12/05 16:18:19  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.13  2014/10/13 08:53:48  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.12  2014/10/06 15:13:22  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.11  2010/10/20 07:58:09  j_novak
 * Better implementation of the explicit time-integration. Not fully-tested yet.
 *
 * Revision 1.10  2008/12/04 18:22:49  j_novak
 * Enhancement of the dzpuis treatment + various bug fixes.
 *
 * Revision 1.9  2008/12/02 15:02:22  j_novak
 * Implementation of the new constrained formalism, following Cordero et al. 2009
 * paper. The evolution eqs. are solved as a first-order system. Not tested yet!
 *
 * Revision 1.8  2004/05/17 19:53:13  e_gourgoulhon
 * Added arguments graph_device and  method_poisson_vect.
 *
 * Revision 1.7  2004/05/12 15:24:20  e_gourgoulhon
 * Reorganized the #include 's, taking into account that
 * time_slice.h contains now an #include "metric.h".
 *
 * Revision 1.6  2004/05/10 09:12:01  e_gourgoulhon
 * Added a call to del_deriv() at the end.
 *
 * Revision 1.5  2004/05/03 14:48:48  e_gourgoulhon
 * Treatment of special cases nn_jp1.etat = ETATUN and psi_jp1.etat = ETATUN.
 *
 * Revision 1.4  2004/04/29 17:10:36  e_gourgoulhon
 * Added argument pdt and update of depth slices at the end,
 * taking into account the known time derivatives.
 *
 * Revision 1.3  2004/04/08 16:45:11  e_gourgoulhon
 * Use of new methods set_*.
 *
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

// C headers
#include <cassert>

// Lorene headers
#include "time_slice.h"
#include "unites.h"
#include "graphique.h"
#include "utilitaires.h"

namespace Lorene {
void Time_slice_conf::initial_data_cts(const Sym_tensor& uu, 
                const Scalar& trk_in, const Scalar& trk_point, 
                double pdt, double precis, int method_poisson_vect,
                const char* graph_device, const Scalar* p_ener_dens, 
                const Vector* p_mom_dens, const Scalar* p_trace_stress) {

    using namespace Unites ;

    // Verifications
    // -------------
    double tr_uu = max(maxabs(uu.trace(tgam()), "trace tgam_{ij} u^{ij}")) ; 
    if (tr_uu > 1.e-7) {
        cerr << 
        "Time_slice_conf::initial_data_cts : the trace of u^{ij} with respect\n"
        << "  to the conformal metric tgam_{ij} is not zero !\n" 
        << "  error = " << tr_uu << endl ; 
        abort() ; 
    }

    assert(trk_in.check_dzpuis(2)) ; 
    assert(trk_point.check_dzpuis(4)) ; 

    // Initialisations
    // ---------------
    double ttime = the_time[jtime] ; 
         
    trk_evol.update(trk_in, jtime, ttime) ; 

    // Reset of quantities depending on K:
    k_dd_evol.downdate(jtime) ; 
    k_uu_evol.downdate(jtime) ; 
   
    set_hata(psi4()*psi()*psi()* uu / (2.* nn()) ) ; 

    const Map& map = uu.get_mp() ; 
    const Base_vect& triad = *(uu.get_triad()) ;
    
    // For graphical outputs:
    int ngraph0 = 10 ;  // index of the first graphic device to be used
    int nz = map.get_mg()->get_nzone() ; 
    double ray_des = 1.25 * map.val_r(nz-2, 1., 0., 0.) ; // outermost radius
                                                          // for plots

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
    
    // Iteration
    // ---------
    int imax = 100 ; 
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

        if (psi_jp1.get_etat() == ETATUN) psi_jp1.std_spectral_base() ; 

        // Test:
        maxabs(psi_jp1.laplacian() - source_psi,
                "Absolute error in the resolution of the equation for Psi") ;  

        des_meridian(psi_jp1, 0., ray_des, "Psi", ngraph0, graph_device) ; 

        // Resolution of the Poisson equation for the lapse
        // ------------------------------------------------
        
        Scalar nn_jp1 = source_nn.poisson() + 1. ; 

        if (nn_jp1.get_etat() == ETATUN) nn_jp1.std_spectral_base() ; 

        // Test:
        maxabs(nn_jp1.laplacian() - source_nn,
                "Absolute error in the resolution of the equation for N") ;  

        des_meridian(nn_jp1, 0., ray_des, "N", ngraph0+1, graph_device) ; 
        
        // Resolution of the vector Poisson equation for the shift
        //---------------------------------------------------------
        
        Vector beta_jp1 = source_beta.poisson(0.3333333333333333, ff, 
                                              method_poisson_vect) ; 
        
        des_meridian(beta_jp1(1), 0., ray_des, "\\gb\\ur\\d", ngraph0+2, 
                     graph_device) ; 
        des_meridian(beta_jp1(2), 0., ray_des, "\\gb\\u\\gh\\d", ngraph0+3, 
                     graph_device) ; 
        des_meridian(beta_jp1(3), 0., ray_des, "\\gb\\u\\gf\\d", ngraph0+4, 
                     graph_device) ; 
        
        // Test:
        Vector test_beta = (beta_jp1.derive_con(ff)).divergence(ff)
            +  0.3333333333333333 * (beta_jp1.divergence(ff)).derive_con(ff) ;
        test_beta.inc_dzpuis() ;  
        maxabs(test_beta - source_beta,
                "Absolute error in the resolution for beta") ; 

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

        set_psi_del_npsi(psi_jp1) ; 
     
        n_evol.update(nn_jp1, jtime, ttime) ; 

        beta_evol.update(beta_jp1, jtime, ttime) ; 

        // New value of A^{ij}:
        Sym_tensor aa_jp1 = ( beta().ope_killing_conf(tgam()) + uu ) 
                                / (2.* nn()) ; 
        
        set_hata( aa_jp1 / (psi4()*psi()*psi()) ) ; 

    }
    
    //==================================================================
    // End of iteration 
    //===================================================================

    npsi_evol.update( n_evol[jtime]*psi_evol[jtime], jtime, ttime ) ;
    A_hata() ;
    B_hata() ;

    // Push forward in time to enable the computation of time derivatives
    // ------------------------------------------------------------------
    
    double ttime1 = ttime ; 
    int jtime1 = jtime ; 
    for (int j=1; j < depth; j++) {
        jtime1++ ; 
        ttime1 += pdt ; 
        psi_evol.update(psi_evol[jtime], jtime1, ttime1) ;  
        npsi_evol.update(npsi_evol[jtime], jtime1, ttime1) ;  
        n_evol.update(n_evol[jtime], jtime1, ttime1) ;  
        beta_evol.update(beta_evol[jtime], jtime1, ttime1) ;  
        hh_evol.update(hh_evol[jtime], jtime1, ttime1) ;
        hata_evol.update(hata_evol[jtime], jtime1, ttime1) ;
	A_hata_evol.update(A_hata_evol[jtime], jtime1, ttime1) ;
	B_hata_evol.update(B_hata_evol[jtime], jtime1, ttime1) ;
        trk_evol.update(trk_evol[jtime], jtime1, ttime1) ;
        the_time.update(ttime1, jtime1, ttime1) ;         
    } 
    jtime += depth - 1 ; 
    
    // Taking into account the time derivative of h^{ij} and K : 
    // ---------------------------------------------------------
    Sym_tensor uu0 = uu ; 
    uu0.dec_dzpuis(2) ; // dzpuis: 2 --> 0
    
    for (int j=1; j < depth; j++) {
        hh_evol.update(hh_evol[jtime] - j*pdt* uu0, 
                       jtime-j, the_time[jtime-j]) ;
                       
        trk_evol.update(trk_evol[jtime] - j*pdt* trk_point, 
                       jtime-j, the_time[jtime-j]) ;
                       
    } 
    
    // Reset of derived quantities (at the new time step jtime)
    // ---------------------------
    del_deriv() ; 
    
} 
}
