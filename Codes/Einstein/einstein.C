/*
 *  Main code for time evolving Einstein equations 
 *   in Dirac gauge.
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

char einstein_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.14  2004/04/05 14:46:33  e_gourgoulhon
 * The graphical functions are now declared in graphique.h and
 * are included in Lorene's library.
 *
 * Revision 1.13  2004/04/05 07:00:49  e_gourgoulhon
 * Corrected minor error in the computation of taa: ff --> tgam.
 *
 * Revision 1.12  2004/03/31 20:28:57  e_gourgoulhon
 * Update to take into account the new prototypes of
 * Evolution_std<TyT>::Evolution_std<TyT> and
 * Evolution_std<TyT>::update.
 * Test of shift equation.
 *
 * Revision 1.11  2004/03/11 12:09:32  e_gourgoulhon
 * Use of new method Scalar::visu_section_anim to produce outputs for
 * movies.
 *
 * Revision 1.10  2004/03/08 00:36:27  e_gourgoulhon
 * Added output for OpenDX (visu_section).
 * Initial amplitude = 1e-3.
 *
 * Revision 1.9  2004/03/06 21:16:38  e_gourgoulhon
 * First version with all equations implemented, with full sources,
 * including the time derivatives.
 * The only missing part is the trace of h.
 *
 * Revision 1.8  2004/03/05 15:11:18  e_gourgoulhon
 * Use of new method Scalar::smooth_decay on khi_jp1.
 *
 * Revision 1.7  2004/03/04 22:19:25  e_gourgoulhon
 * All sources completed (except for the time derivatives and
 * shift terms in the source for h).
 *
 * Revision 1.6  2004/03/04 16:25:56  e_gourgoulhon
 * Still in progress...
 *
 * Revision 1.5  2004/03/03 11:35:25  e_gourgoulhon
 * First version with Evolution_std's and d'Alembert.
 *
 * Revision 1.4  2004/03/02 14:54:17  e_gourgoulhon
 * Started to encode source for h from new equations.
 *
 * Revision 1.3  2004/02/27 21:17:26  e_gourgoulhon
 * Still in progress...
 *
 * Revision 1.2  2004/02/19 22:16:42  e_gourgoulhon
 * Sources of equations for Q, N and beta completed.
 *
 * Revision 1.1  2004/02/18 19:16:28  e_gourgoulhon
 * First version: c'est loin d'etre pret tout ca !!!
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Lorene headers
#include "metric.h"
#include "evolution.h"
#include "param.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"

int main() {

    //======================================================================
    //      Parameters of the computation
    //======================================================================

    double pdt = 0.01 ; 
    int jmax = 1000 ; 
    int jstop = jmax ; 
    bool compute_source = true ; 

    double relativistic_init = 0. ;     // 0 = flat space
    double ampli_h_init = 0.001 ;     // 0 = flat space
        

    //======================================================================
    //      Construction and initialization of the various objects
    //======================================================================

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 4 ; 	// Number of domains
    int nr = 17 ; 	// Number of collocation points in r in each domain
    int nt = 9 ; 	// Number of collocation points in theta in each domain
    int np = 8 ; 	// Number of collocation points in phi in each domain
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = SYM ; // no symmetry in phi
    bool compact = true ; // external domain is compactified
  
    // Multi-domain grid construction:
    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
    cout << mgrid << endl ; 

  
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------

    // radial boundaries of each domain:
    double r_limits[] = {0., 1., 2., 4., __infinity} ; 
    assert( nz == 4 ) ;  // since the above array described only 3 domains
  
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    cout << map << endl ;  
    
    // Flat metric f
    // -------------

    const Metric_flat& ff = map.flat_met_spher() ; 
    
    // Triad orthonormal with respect to the flat metric f
    // ----------------------------------------------------

    const Base_vect_spher& otriad = map.get_bvect_spher() ;
    
    
    // Set up of tensor h
    // ------------------
    
    Sym_tensor_trans hh_init(map, otriad, ff) ;  // hh is a transverse tensor
                                            // with respect to the flat metric
                                            // thanks to Dirac gauge
    
    // Test with the tensor h^{ij} = D^i D^j Phi  with Lap(Phi) = 0
    
    const Coord& x = map.x ; 
    const Coord& y = map.y ; 
    // const Coord& z = map.z ; 
    const Coord& r = map.r ; 
    //const Coord& cost = map.cost ; 
    //const Coord& sint = map.sint ; 
    //const Coord& cosp = map.cosp ; 
    // const Coord& sinp = map.sinp ; 
    
    Scalar khi_init(map) ; 

    khi_init = ampli_h_init * exp( - r*r ) * x*y ;
    
    //khi_init = ampli_h_init * (3*cost*cost-1) / 
    //   ( (r*r + 1./(r*r)) * sqrt(1.+r*r) ) ; 
    khi_init.std_spectral_base() ; 

    khi_init.smooth_decay(2, 1) ; 
    
    //##
    // des_meridian(khi_init, 0., 3., "khi_init before", 1) ; 
    // arrete() ; 
    // khi_init.smooth_decay(3, 4) ; 
    // khi_init.spectral_display("khi_init") ;   
    //  des_meridian(khi_init, 0., 3., "khi_init after", 2) ; 
    //  arrete() ; 
    //## 
                 
    Scalar mu_init(map) ; 
    mu_init = 0. * ampli_h_init / (1+r*r*r*r*r*r) ; 
    mu_init.std_spectral_base() ; 
    mu_init.mult_r() ; 
    mu_init.mult_r() ; 
    mu_init.mult_r() ; 
    mu_init.mult_cost() ; 
    
    //##
    // des_meridian(mu_init, 0., 3., "mu_init before", 1) ; 
    // arrete() ; 
    // mu_init.smooth_decay(3, 4) ; 
    // mu_init.spectral_display("mu_init") ;   
    // des_meridian(mu_init, 0., 3., "mu_init after", 2) ; 
    // arrete() ; 
    //## 
                  

    Sym_tensor_tt htt_init(map, otriad, ff) ;  // htt is the TT part of hh
        
    htt_init.set_khi_mu(khi_init, mu_init) ; 
    
    hh_init = htt_init ; 
        
    // des_meridian(hh_init, 0., 5., "hh_init") ; 
    maxabs( hh_init.divergence(ff), "Divergence of hh_init") ; 
    maxabs( hh_init.trace(ff), "Trace of hh_init") ; 

    arrete() ; 

    // Set up of field Q = Psi^2 N
    // ---------------------------
    
    Scalar qq_init(map) ; 
    Scalar tmp(map) ; 
    
    qq_init = 1. - relativistic_init / (1. + r*r) ;  
     
    qq_init.std_spectral_base() ;    // sets standard spectral bases


    // Set up of conformal metric gamma_tilde
    // --------------------------------------
    
    Metric tgam( ff.con() ) ;   // construction from the flat metric

    tgam = ff.con() + hh_init ;      // initialization  [ Eq. (51) ]
    

    // Set up of shift vector beta
    // ---------------------------    

    Vector beta_init(map, CON, otriad ) ; 
    beta_init.set_etat_zero() ; 

    // Set up of lapse function N
    // --------------------------
    
    Scalar nn_init(map) ; 

    nn_init = 1. - relativistic_init / sqrt(1. + r*r) ; 

    nn_init.std_spectral_base() ;    // sets standard spectral bases

    // Working stuff
    // -------------
    
    Scalar tmp0(map) ; 
    Sym_tensor sym_tmp(map, CON, otriad) ; 

    //======================================================================
    //                  Start of time evolution
    //======================================================================
    
    double ttime = 0. ; 
    int jtime = 0 ; 

    Evolution_std<Scalar> nn_time(nn_init, 3) ; 
    Evolution_std<Vector> beta_time(beta_init, 3) ; 
    Evolution_std<Scalar> qq_time(qq_init, 3) ; 
    Evolution_std<Sym_tensor_trans> hh_time(hh_init, 3) ; 
    Evolution_std<Scalar> khi_time(khi_init, 3) ; 
    Evolution_std<Scalar> mu_time(mu_init, 3) ; 
    
    ttime += pdt ; 
    jtime++ ; 
    nn_time.update(nn_init, jtime, ttime) ; 
    beta_time.update(beta_init, jtime, ttime) ; 
    qq_time.update(qq_init, jtime, ttime) ; 
    hh_time.update(hh_init, jtime, ttime) ; 
    khi_time.update(khi_init, jtime, ttime) ; 
    mu_time.update(mu_init, jtime, ttime) ; 
    
    ttime += pdt ; 
    jtime++ ; 
    nn_time.update(nn_init, jtime, ttime) ; 
    beta_time.update(beta_init, jtime, ttime) ; 
    qq_time.update(qq_init, jtime, ttime) ; 
    hh_time.update(hh_init, jtime, ttime) ; 
    khi_time.update(khi_init, jtime, ttime) ; 
    mu_time.update(mu_init, jtime, ttime) ; 
    
    
    // Parameters for the d'Alembert equations
    // ----------------------------------------
    int bc = 2 ;    // type of boundary condition : 2 = Bayliss & Turkel outgoing wave
 
    Param par_khi ; 
    par_khi.add_double(pdt) ; 
    par_khi.add_int(bc) ; 
    int *workflag_khi = new int(0) ; // working flag 
    par_khi.add_int_mod(*workflag_khi) ; 
    
    Param par_mu ; 
    par_mu.add_double(pdt) ; 
    par_mu.add_int(bc) ; 
    int *workflag_mu = new int(0) ; // working flag 
    par_mu.add_int_mod(*workflag_mu) ; 

    
    for (jtime = 2; jtime <= jmax; jtime++) {
    
        cout << 
        "==============================================================\n"
        << "  step: " << jtime << "   time = " << ttime << endl  
        << "==============================================================\n" ;
    
        // Values at time step jtime: 
        const Scalar& nn = nn_time[jtime] ; 
        const Vector& beta = beta_time[jtime] ; 
        const Scalar& qq = qq_time[jtime] ; 
        const Sym_tensor_trans& hh = hh_time[jtime] ;
        
        // Time derivatives:
        Scalar nn_point = nn_time.time_derive(jtime) ; 
        nn_point.inc_dzpuis(2) ; // dzpuis : 0 -> 2
        
        Vector beta_point = beta_time.time_derive(jtime) ; 
        beta_point.inc_dzpuis(2) ; // dzpuis : 0 -> 2
        
        Sym_tensor_trans hh_point = hh_time.time_derive(jtime) ; 
        hh_point.inc_dzpuis(2) ; // dzpuis : 0 -> 2
        
        // Sources of the Einstein equations:
        Scalar source_nn(map) ; 
        Vector source_beta(map, CON, otriad) ;
        Scalar source_qq(map) ; 
        Sym_tensor source_hh(map, CON, otriad) ;  
        
        if (compute_source) {
        //==============================================
        //  Definition of references on derivatives: 
        //   the source objects should not be modified
        //   in this scope
        //==============================================  

        Scalar psi = sqrt(qq / nn) ;   
        psi.std_spectral_base() ;    

        Scalar ln_psi = log( psi ) ;  
        ln_psi.std_spectral_base() ;    
        Scalar psi2 = psi * psi ; 
        Scalar psi4 = psi2 * psi2 ; 

        const Sym_tensor& tgam_dd = tgam.cov() ;    // {\tilde \gamma}_{ij}
        const Sym_tensor& tgam_uu = tgam.con() ;    // {\tilde \gamma}^{ij}
        const Tensor_sym& dtgam = tgam_dd.derive_cov(ff) ;    
                                                    // D_k {\tilde \gamma}_{ij}
        const Tensor_sym& dhh = hh.derive_cov(ff) ; // D_k h^{ij}
        const Vector& dln_psi = ln_psi.derive_cov(ff) ; // D_i ln(Psi)
        const Vector& tdln_psi_u = ln_psi.derive_con(tgam) ; // tD^i ln(Psi)
        const Vector& dnn = nn.derive_cov(ff) ;         // D_i N
        const Vector& tdnn_u = nn.derive_con(tgam) ;       // tD^i N
        const Vector& dqq = qq.derive_cov(ff) ;         // D_i Q
        const Scalar& div_beta = beta.divergence(ff) ;  // D_k beta^k
        
        // Conformal Killing operator applied to beta
        // ------------------------------------------
        Sym_tensor l_beta(map, CON, otriad) ;   // (L beta)^{ij}
        for (int i=1; i<=3; i++) {
            for (int j=1; j<=i; j++) {
                l_beta.set(i,j) = beta.derive_con(ff)(i,j)  
                                + beta.derive_con(ff)(j,i) 
                                - 0.6666666666666666* div_beta * ff.con()(i,j) ;
            }
        }

        // Conformal extrinsic curvature A
        // -------------------------------
        Sym_tensor aa(map, CON, otriad) ;   // A^{ij}

        aa = ( hh_point - hh.derive_lie(beta) + l_beta 
                - 0.6666666666666666 * div_beta * hh ) / (2.*nn) ; 
    
        Sym_tensor taa(map, COV, otriad) ;  // {\tilde A}_{ij}
        taa = aa.up_down(tgam) ;

        des_meridian(aa(2,3), 0., 5., "A\\u\\gh\\gf\\d", 20) ; 
        des_meridian(aa(3,3), 0., 5., "A\\u\\gf\\gf\\d", 21) ; 
        

       // Source for Q  [ Eq. (76) ]
        // ------------
        
        Scalar aa_quad = contract(taa, 0, 1, aa, 0, 1) ; 
        
        source_qq = 0.75 * psi4 * qq * aa_quad ;
        
        tmp = contract( hh, 0, 1, dqq.derive_cov(ff), 0, 1 ) ;             
        tmp.inc_dzpuis() ; 
        
        source_qq -= tmp ;  
                        
        tmp = 0.0625 * contract( dhh, 0, 1, dtgam, 0, 1 ).trace(tgam) 
             - 0.125 * contract( dhh, 0, 1, dtgam, 0, 2 ).trace(tgam) 
             + 2.* contract( contract( tgam_uu, 0, dln_psi, 0), 0,
                             dln_psi, 0 ) ;
     
        tmp0 = 2. * contract( tgam_uu, 0, 1, 
                              dln_psi * dnn, 0, 1) ;
        
        source_qq += psi2 * ( nn * tmp + tmp0 ) ; 
                             
        // source_qq.spectral_display("source_qq") ; 


        // Source for N  [ Eq. (80) ]
        // ------------
        
        source_nn = psi4 * nn * aa_quad ;
        
        tmp = contract( hh, 0, 1, dnn.derive_cov(ff), 0, 1 ) ;
        tmp.inc_dzpuis() ; 
        
        source_nn -= tmp + tmp0 ; 
                    
        // source_nn.spectral_display("source_nn") ; 

        // Source for beta [ Eq. (79) ]
        // ---------------

        source_beta = 2. * ( contract(aa, 1, 
                        dnn - 6.*nn * dln_psi, 0)
                - nn * contract(tgam.connect().get_delta(), 1, 2, aa, 0, 1) ) ;
                
        Vector vtmp = contract(hh, 0, 1, 
                           beta.derive_cov(ff).derive_cov(ff), 1, 2)
                + 0.3333333333333333 *
                  contract(hh, 1, div_beta.derive_cov(ff), 0) ; 
        vtmp.inc_dzpuis() ; // dzpuis: 3 -> 4
                    
        source_beta -= vtmp ; 
        
        // source_beta.spectral_display("source_beta") ; 
        

        // Quadratic part of the Ricci tensor of gam_tilde 
        // ------------------------------------------------
        
        Sym_tensor ricci_star(map, CON, otriad) ; 
        
        ricci_star = contract(hh, 0, 1, dhh.derive_cov(ff), 2, 3) ; 

        ricci_star.inc_dzpuis() ;   // dzpuis : 3 --> 4


        for (int i=1; i<=3; i++) {
            for (int j=1; j<=i; j++) {
                tmp = 0 ; 
                for (int k=1; k<=3; k++) {
                    for (int l=1; l<=3; l++) {
                        tmp += dhh(i,k,l) * dhh(j,l,k) ; 
                    }
                }
                sym_tmp.set(i,j) = tmp ; 
            }
        }
        ricci_star -= sym_tmp ;

        for (int i=1; i<=3; i++) {
            for (int j=1; j<=i; j++) {
                tmp = 0 ; 
                for (int k=1; k<=3; k++) {
                    for (int l=1; l<=3; l++) {
                        for (int m=1; m<=3; m++) {
                            for (int n=1; n<=3; n++) {
                            
        tmp += 0.5 * tgam_uu(i,k)*tgam_uu(j,l) * dhh(m,n,k) * dtgam(m,n,l)
                + tgam_dd(n,l) * dhh(m,n,k) * ( tgam_uu(i,k) * dhh(j,l,m)
                                                + tgam_uu(j,k) *  dhh(i,l,m) )
                - tgam_dd(k,l) * tgam_uu(m,n) * dhh(i,k,m) * dhh(j,l,n) ;
                            }
                        } 
                    }
                }
                sym_tmp.set(i,j) = tmp ; 
            }
        }
        ricci_star += sym_tmp ;

        ricci_star = 0.5 * ricci_star ; 
        
        // des_meridian(ricci_star(1,1), 0., 4., "Ricci_star^rr", 12) ; 
        
        // Curvature scalar of conformal metric :
        // -------------------------------------
        
        Scalar tricci_scal = 
            0.25 * contract( tgam_uu, 0, 1,
                             contract(dhh, 0, 1, dtgam, 0, 1), 0, 1 ) 
          - 0.5  * contract( tgam_uu, 0, 1,
                             contract(dhh, 0, 1, dtgam, 0, 2), 0, 1 ) ;  
                                                       
        // Full quadratic part of source for h : S^{ij}
        // --------------------------------------------
        
        Sym_tensor ss(map, CON, otriad) ; 
        
        sym_tmp = nn * (ricci_star + 8.* tdln_psi_u * tdln_psi_u)
                + 4.*( tdln_psi_u * tdnn_u + tdnn_u * tdln_psi_u ) 
                - 0.3333333333333333 * ( 
                    nn * (tricci_scal 
                            + 8.* contract(dln_psi, 0, tdln_psi_u, 0) )
                    + 8.* contract(dln_psi, 0, tdnn_u, 0) ) * tgam_uu ;

        ss = sym_tmp / psi4  ;
        
        sym_tmp = contract( tgam_uu, 1, 
                            contract(tgam_uu, 1, dqq.derive_cov(ff), 0), 1) ;
                            
        sym_tmp.inc_dzpuis() ; // dzpuis : 3 --> 4
        
        for (int i=1; i<=3; i++) {
            for (int j=1; j<=i; j++) {
                tmp = 0 ; 
                for (int k=1; k<=3; k++) {
                    for (int l=1; l<=3; l++) {
                        tmp += ( hh(i,k)*dhh(l,j,k) + hh(k,j)*dhh(i,l,k)
                                    - hh(k,l)*dhh(i,j,k) ) * dqq(l) ; 
                    }
                }
                sym_tmp.set(i,j) += 0.5 * tmp ; 
            }
        }
        
        tmp = qq.derive_con(tgam).divergence(tgam) ; 
        tmp.inc_dzpuis() ; // dzpuis : 3 --> 4
        
        sym_tmp -= 0.3333333333333333 * tmp * tgam_uu ; 
                    
        ss -= sym_tmp / (psi4*psi2) ; 

        for (int i=1; i<=3; i++) {
            for (int j=1; j<=i; j++) {
                tmp = 0 ; 
                for (int k=1; k<=3; k++) {
                    for (int l=1; l<=3; l++) {
                        tmp += tgam_dd(k,l) * aa(i,k) * aa(j,l) ; 
                    }
                }
                sym_tmp.set(i,j) = 2. * tmp ; 
            }
        }
        
        ss += nn * sym_tmp ; 

        // Source for h 
        // ------------
                 
        Sym_tensor lbh = hh.derive_lie(beta) ; 

        source_hh = (nn*nn/psi4 - 1.) * hh.derive_con(ff).divergence(ff) 
            + 2.* hh_point.derive_lie(beta) - lbh.derive_lie(beta) ;
        source_hh.inc_dzpuis() ; 
        
        source_hh += 2.* nn * ss ;
              
        //## Provisory: waiting for the Lie derivative to allow
        //  derivation with respect to a vector with dzpuis != 0
        vtmp = beta_point ; 
        vtmp.dec_dzpuis(2) ; 
        sym_tmp = hh.derive_lie(vtmp) ; 
        sym_tmp.inc_dzpuis(2) ;             

        source_hh += sym_tmp 
              + 1.3333333333333333 * div_beta* (hh_point - lbh)
              + 2. * (nn_point - nn.derive_lie(beta)) * aa  ;
              

        for (int i=1; i<=3; i++) {
            for (int j=1; j<=i; j++) {
                tmp = 0 ; 
                for (int k=1; k<=3; k++) {
                    tmp += ( hh.derive_con(ff)(k,j,i) 
                           + hh.derive_con(ff)(i,k,j) 
                           - hh.derive_con(ff)(i,j,k) ) * dqq(k) ;
                }
                sym_tmp.set(i,j) = tmp ; 
            }
        }
            
        source_hh -= nn / (psi4*psi2) * sym_tmp ; 
         
        tmp =  beta_point.divergence(ff) - div_beta.derive_lie(beta) ; 
        tmp.inc_dzpuis() ; 
        source_hh += 0.6666666666666666* ( tmp
             - 0.6666666666666666* div_beta * div_beta ) * hh ; 
               
        
        // (L beta_point)^ij --> sym_tmp
        for (int i=1; i<=3; i++) {
            for (int j=1; j<=i; j++) {
                sym_tmp.set(i,j) = beta_point.derive_con(ff)(i,j)  
                                 + beta_point.derive_con(ff)(j,i) 
                                 - 0.6666666666666666* 
                                    beta_point.divergence(ff) * ff.con()(i,j) ;
            }
        }
        sym_tmp -= l_beta.derive_lie(beta) ;
        sym_tmp.inc_dzpuis() ; 
        
        source_hh += 0.6666666666666666* div_beta * l_beta - sym_tmp ; 
           
        // des_meridian(source_hh(1,1), 0., 2., "source_hh^rr", 10) ; 

        source_hh.spectral_display("source_hh") ; 
        maxabs(source_hh, "Maxabs source_hh") ; 

        maxabs( source_hh.divergence(ff), "Divergence of source_hh") ; 
        maxabs( source_hh.transverse(ff).divergence(ff), 
                "Divergence of source_hh_transverse") ; 
        maxabs( source_hh.transverse(ff).trace(ff), 
                "Trace of source_hh_transverse") ; 

        arrete(jtime%jstop) ;  

        }
        else{
            source_hh.set_etat_zero() ; 
        }
        //==============================================
        //  End of scope for references on derivatives
        //==============================================  

        //=============================================
        // Resolution of elliptic equations
        //=============================================
        
        // Resolution of the Poisson equation for the lapse
        // ------------------------------------------------
        
        Scalar nn_jp1 = source_nn.poisson() + 1. ; 

        //des_meridian(nn_jp1, 0., 5., "N", 40) ; 

        // Resolution of the Poisson equation for Q
        // -----------------------------------------
        
        Scalar qq_jp1 = source_qq.poisson() + 1. ; 

        // des_meridian(qq_jp1, 0., 5., "Q", 41) ; 
        
        // Resolution of the vector Poisson equation for the shift
        //---------------------------------------------------------
        
        Vector beta_jp1 = source_beta.poisson(0.3333333333333333, 0) ; 
        des_meridian(beta_jp1(1), 0., 5., "\\gb\\ur\\d", 42) ; 
        des_meridian(beta_jp1(2), 0., 5., "\\gb\\u\\gh\\d", 43) ; 
        des_meridian(beta_jp1(3), 0., 5., "\\gb\\u\\gf\\d", 44) ; 
        
        // Test:
        Vector test_beta = (beta_jp1.derive_con(ff)).divergence(ff)
            +  0.3333333333333333 * (beta_jp1.divergence(ff)).derive_con(ff) ;
        test_beta.inc_dzpuis() ;  
        cout << "Relative error in the resolution of the shift equation: \n" ;
        diffrel(test_beta, source_beta) ; 
        diffrelmax(test_beta, source_beta) ; 
        // arrete() ; 
        
        
        //=============================================
        // Resolution of wave equation for h
        //=============================================
    
        const Sym_tensor_tt& source_htt = source_hh.transverse(ff).tt_part() ;
  
        //## Sym_tensor_tt source_htt(map, otriad, ff) ; 
        //## source_htt.set_etat_zero() ; 
        
        maxabs( source_htt.divergence(ff), "Divergence of source_htt") ; 
        maxabs( source_htt.trace(ff), "Trace of source_hhtt") ; 

        const Scalar& khi_source = source_htt.khi() ; 
        
        const Scalar& mu_source = source_htt.mu() ; 
       
        // des_meridian(khi_source, 0., 2., "khi_source", 11) ; 
                     
        Scalar khi_jp1 = khi_time[jtime].avance_dalembert(par_khi,
                                         khi_time[jtime-1], khi_source) ;
        khi_jp1.smooth_decay(2,2) ; 
        
        // des_meridian(khi_jp1, 0., 5., "khi_jp1", 12) ; 
        
        // const Scalar* khides[] = {&khi_jp1, &(khi_time[jtime]), 
        //                           &(khi_time[jtime-1])} ; 
        // double thetades[] = {0., 0., 0.} ; 
        // double phides[] = {0., 0., 0.} ; 
        // des_profile_mult(khides, 3, 0., 4., thetades, phides, 1, false, "khi") ;
        
        maxabs(khi_jp1 - khi_time[jtime], "Variation of khi") ;  
        
        Scalar mu_jp1 = mu_time[jtime].avance_dalembert(par_mu,  
                                         mu_time[jtime-1], mu_source) ;
        mu_jp1.smooth_decay(2,2) ; 
                                        
        // des_meridian(mu_jp1, 0., 5., "mu_jp1", 14) ; 
        
                
        Sym_tensor_tt htt_jp1(map, otriad, ff) ;
        
        htt_jp1.set_khi_mu(khi_jp1, mu_jp1) ;
                 
        Sym_tensor_trans hh_jp1 = htt_jp1 ;     //## the trace should be added
    
        // des_meridian(hh_jp1, 0., 5., "hh") ; 
        des_meridian(hh(2,3), 0., 5., "h\\u\\gh\\gf\\d", 30) ; 
        des_meridian(hh(3,3), 0., 5., "h\\u\\gf\\gf\\d", 31) ; 

        // hh(2,3).visu_section('z', 0., -4., 4., -4., 4., "h^tp", "h_tp") ;
        
        hh(3,3).visu_section_anim('z', 0., -4., 4., -4., 4., jtime, ttime, 4,
                                  "h^pp", "hpp") ;

        cout << "Next step : " << jtime + 1 << endl ; 
        arrete(jtime%jstop) ;  

        // Next time step         
        // --------------
        
        ttime += pdt ; 
        
        nn_time.update(nn_jp1, jtime+1, ttime) ; 
        beta_time.update(beta_jp1, jtime+1, ttime) ; 
        qq_time.update(qq_jp1, jtime+1, ttime) ; 
        hh_time.update(hh_jp1, jtime+1, ttime) ; 
        khi_time.update(khi_jp1, jtime+1, ttime) ; 
        mu_time.update(mu_jp1, jtime+1, ttime) ; 
        
    }


    par_khi.clean_all() ; 
    par_mu.clean_all() ; 

    return EXIT_SUCCESS ; 
}

