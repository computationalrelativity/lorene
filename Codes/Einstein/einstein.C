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

// Lorene headers
#include "metric.h"
#include "evolution.h"
#include "param.h"
#include "nbr_spx.h"
#include "utilitaires.h"

int main() {


    //======================================================================
    //      Construction and initialization of the various objects
    //======================================================================

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 3 ; 	// Number of domains
    int nzm1 = nz - 1 ;
    int nr = 17 ; 	// Number of collocation points in r in each domain
    int nt = 5 ; 	// Number of collocation points in theta in each domain
    int np = 8 ; 	// Number of collocation points in phi in each domain
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = NONSYM ; // no symmetry in phi
    bool compact = true ; // external domain is compactified
  
    // Multi-domain grid construction:
    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
    cout << mgrid << endl ; 

  
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------

    // radial boundaries of each domain:
    double r_limits[] = {0., 2., 3., __infinity} ; 
    assert( nz == 3 ) ;  // since the above array described only 3 domains
  
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    cout << map << endl ;  
    
    // Flat metric f
    // -------------

    const Metric_flat& ff = map.flat_met_spher() ; 
    
    // Triad orthonormal with respect to the flat metric f
    // ----------------------------------------------------

    const Base_vect_spher& otriad = map.get_bvect_spher() ;
    
    // Parameter for the initial data 
    //-------------------------------
    
    double relativistic_init = 0.1 ;     // 0 = flat space
    
    
    // Set up of tensor h
    // ------------------
    
    Sym_tensor_trans hh_init(map, otriad, ff) ;  // hh is a transverse tensor
                                            // with respect to the flat metric
                                            // thanks to Dirac gauge
    
    // Test with the tensor h^{ij} = D^i D^j Phi  with Lap(Phi) = 0
    
    const Coord& x = map.x ; 
    const Coord& y = map.y ; 
    const Coord& z = map.z ; 
    const Coord& r = map.r ; 
    const Coord& cost = map.cost ; 
    const Coord& sint = map.sint ; 
    const Coord& cosp = map.cosp ; 
    const Coord& sinp = map.sinp ; 
    const Coord& phi = map.phi ; 
    
    Sym_tensor_tt htt(map, otriad, ff) ;  // htt is the TT part of hh
    
    Scalar htt_rr(map) ; 
    htt_rr = relativistic_init * (3*cost*cost-1) / (r*r*r*r + 1./(r*r)) ; 
    htt_rr.std_spectral_base() ; 
    
    Scalar htt_mu(map) ; 
    htt_mu = relativistic_init / (1+r*r*r*r*r*r) ; 
    htt_mu.std_spectral_base() ; 
    htt_mu.mult_r() ; 
    htt_mu.mult_cost() ; 
    
    htt.set_rr_mu(htt_rr, htt_mu) ; 
    
    hh_init = htt ; 
    
    hh_init.annule_domain(nzm1) ;           // h set to zero in the CED     


    // Set up of field Q = Psi^2 N
    // ---------------------------
    
    Scalar qq_init(map) ; 
    Scalar tmp(map) ; 
    
    qq_init = 1. + relativistic_init * r*r ; 
    tmp = 1. + relativistic_init / (r*r) ; 
    qq_init.set_domain(nzm1) = tmp.domain(nzm1) ; 
     
    qq_init.std_spectral_base() ;    // sets standard spectral bases


    // Set up of conformal metric gamma_tilde
    // --------------------------------------
    
    Metric tgam( ff.con() ) ;   // construction from the flat metric

    tgam = ff.con() + hh_init ;      // initialization  [ Eq. (51) ]
    

    // Set up of shift vector beta
    // ---------------------------    

    Vector beta_init(map, CON, otriad ) ; 
    tmp =  sint * cosp / ( r + 1 / r ) ;   
    tmp.std_spectral_base() ;    // sets standard spectral bases
    beta_init = relativistic_init * tmp.derive_con(ff) ;   

    // Set up of lapse function N
    // --------------------------
    
    Scalar nn_init(map) ; 

    nn_init = 1. - relativistic_init * r*r ; 
    tmp = 1. - relativistic_init / (r*r) ; 
    nn_init.set_domain(nzm1) = tmp.domain(nzm1) ; 

    nn_init.std_spectral_base() ;    // sets standard spectral bases

    // Working stuff
    // -------------
    
    Scalar tmp0(map) ; 
    Sym_tensor sym_tmp(map, CON, otriad) ; 

    //======================================================================
    //                  Start of time evolution
    //======================================================================
    
    double pdt = 0.01 ; 
    double ttime = 0 ; 
    

    Evolution_std<Scalar> nn_time(nn_init, ttime, 3) ; 
    Evolution_std<Vector> beta_time(beta_init, ttime, 3) ; 
    Evolution_std<Scalar> qq_time(qq_init, ttime, 3) ; 
    Evolution_std<Sym_tensor_trans> hh_time(hh_init, ttime, 3) ; 
    
    ttime += pdt ; 
    nn_time.update(nn_init, ttime) ; 
    beta_time.update(beta_init, ttime) ; 
    qq_time.update(qq_init, ttime) ; 
    hh_time.update(hh_init, ttime) ; 
    
    ttime += pdt ; 
    nn_time.update(nn_init, ttime) ; 
    beta_time.update(beta_init, ttime) ; 
    qq_time.update(qq_init, ttime) ; 
    hh_time.update(hh_init, ttime) ; 
    
    // khi and mu at previous time step for the resolution 
    // of the wave equation for h^{ij}_TT
    //-----------------------------------------------------
    Scalar khi_prev = ((hh_time[1]).tt_part())(1,1) ;  
    khi_prev.annule_domain(nzm1) ; 
    khi_prev.mult_r() ; 
    khi_prev.mult_r() ; 
    
    
    Scalar mu_prev = ((hh_time[1]).tt_part()).mu() ;  
    mu_prev.annule_domain(nzm1) ; 

    
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

    int jmax = 3 ; 
    
    for (int jtime = 2; jtime <= jmax; jtime++) {
    
        
        const Scalar& nn = nn_time[jtime] ; 
        const Vector& beta = beta_time[jtime] ; 
        const Scalar& qq = qq_time[jtime] ; 
        const Sym_tensor_trans& hh = hh_time[jtime] ;
        
        Scalar source_nn(map) ; 
        Vector source_beta(map, CON, otriad) ;
        Scalar source_qq(map) ; 
        Sym_tensor_trans source_hh(map, otriad, ff) ;  
        

        {
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

        // Conformal extrinsic curvature A
        // -----------------------------------------
    
        Sym_tensor aa(map, CON, otriad) ;   // A^{ij}
        aa.set_etat_zero() ;    //  initialization to zero
    
        Sym_tensor taa(map, COV, otriad) ;  // {\tilde A}_{ij}
        taa = aa.up_down(ff) ;
        

       // Source for Q  [ Eq. (76) ]
        // ------------
        
        Scalar aa_quad = contract(taa, 0, 1, aa, 0, 1) ; 
        
        source_qq = 0.75 * psi4 * qq * aa_quad  
            - contract( hh, 0, 1, dqq.derive_cov(ff), 0, 1 ) ; 
            
        source_qq.inc_dzpuis() ;  
            
        tmp = 0.0625 * contract( dhh, 0, 1, dtgam, 0, 1 ).trace(tgam) 
             - 0.125 * contract( dhh, 0, 1, dtgam, 0, 2 ).trace(tgam) 
             + 2.* contract( contract( tgam_uu, 0, dln_psi, 0), 0,
                             dln_psi, 0 ) ;
     
        tmp0 = 2. * contract( tgam_uu, 0, 1, 
                              dln_psi * dnn, 0, 1) ;
        
        source_qq += psi2 * ( nn * tmp + tmp0 ) ; 
                             
        source_qq.spectral_display("source_qq") ; 


        // Source for N  [ Eq. (80) ]
        // ------------
        
        source_nn = psi4 * nn * aa_quad ;
        
        tmp = contract( hh, 0, 1, dnn.derive_cov(ff), 0, 1 ) ;
        tmp.inc_dzpuis() ; 
        
        source_nn -= tmp + tmp0 ; 
                    
        source_nn.spectral_display("source_nn") ; 

        // Source for beta [ Eq. (79) ]
        // ---------------

        source_beta = 2. * ( contract(aa, 1, 
                        dnn - 6.*nn * dln_psi, 0)
                - nn * contract(tgam.connect().get_delta(), 1, 2, aa, 0, 1) )
                - contract(hh, 0, 1, 
                           beta.derive_cov(ff).derive_cov(ff), 1, 2)
                - 0.3333333333333333 *
                  contract(hh, 1, beta.divergence(ff).derive_cov(ff), 0) ; 
                    
        source_beta.spectral_display("source_beta") ; 
        

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
        
        ss = nn * (ricci_star + 8.* tdln_psi_u * tdln_psi_u)
                + 4.*( tdln_psi_u * tdnn_u + tdnn_u * tdln_psi_u ) 
                - 0.3333333333333333 * ( 
                    nn * (tricci_scal 
                            + 8.* contract(dln_psi, 0, tdln_psi_u, 0) )
                    + 8.* contract(dln_psi, 0, tdnn_u, 0) ) * tgam_uu ;

        ss = ss / psi4  ;


        // Source for h 
        // ------------
                 
        source_hh = (nn*nn/psi4 - 1.) * hh.derive_con(ff).divergence(ff) ;
        
        source_hh.inc_dzpuis() ; 
        

        source_hh += 2. * nn * ss ;       
        

        }
        //==============================================
        //  End of scope for references on derivatives
        //==============================================  

        //=============================================
        // Resolution of wave equation for h
        //=============================================
    
        const Sym_tensor_tt& source_htt = source_hh.tt_part() ;
        
        Scalar khi_source = source_htt(1,1) ; 
        khi_source.mult_r() ;  
        khi_source.mult_r() ;  
        
        const Scalar& mu_source = source_htt.mu() ; 
        
        Scalar khi_j = hh.tt_part()(1,1) ; 
        khi_j.mult_r() ; 
        khi_j.mult_r() ; 
             
        Scalar khi_jp1 = khi_j.avance_dalembert(par_khi, khi_prev, khi_source) ;
        
        Scalar mu_jp1 = hh.tt_part().mu().avance_dalembert(par_mu, mu_prev, 
                                                            mu_source) ;
        
        Sym_tensor_tt htt_jp1(map, otriad, ff) ;
        
	    khi_jp1.div_r() ;   		// division 
	    khi_jp1.div_r() ;   		// by r^2 --> khi_jp1 now contains h^{rr}

        htt_jp1.set_rr_mu(khi_jp1, mu_jp1) ;
         
        htt_jp1.annule_domain(nzm1) ; 
        
        Sym_tensor_trans hh_jp1 = htt_jp1 ;    

        cout << "Next step ?" << endl ; 
        arrete() ;  

        // Next time step         
        // --------------
        
        ttime += pdt ; 
        
        khi_prev = khi_j ; 
        khi_prev.annule_domain(nzm1) ; 

        mu_prev = hh.tt_part().mu() ; 
        mu_prev.annule_domain(nzm1) ; 
                

        nn_time.update(nn, ttime) ; 
        beta_time.update(beta, ttime) ; 
        qq_time.update(qq, ttime) ; 
        hh_time.update(hh_jp1, ttime) ; 
        
    }


    par_khi.clean_all() ; 
    par_mu.clean_all() ; 

    return EXIT_SUCCESS ; 
}
