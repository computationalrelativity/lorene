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
    
    Sym_tensor_trans hh(map, otriad, ff) ;  // hh is a transverse tensor
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

    Scalar pot_hh(map) ; 
    pot_hh =  // 1                    // P_0^0
           + x                  // P_1^1 cos(p)
           + y                  // P_1^1 sin(p)
           + (3*z*z - r*r)      // P_2^0
           + (x*x - y*y )       // P_2^2 cos(2p)
           + x*y ;                // P_2^2 sin(2p)
     //      + x*(5*z*z-r*r)      // P_3^1 cos(p)
     //      + y*(5*z*z-r*r)      // P_3^1 sin(p)
     //      + x*(x*x-3*y*y)      // P_3^3 cos(3p)
     //      + y*(y*y-3*x*x) ;    // P_3^3 sin(3p)
           
    pot_hh.annule_domain(nzm1) ; 

    Mtbl potced = // 1 / r         // P_0^0
           + sint*cosp / pow(r,2)   // P_1^1 cos(p)
           + sint*sinp / pow(r,2)   // P_1^1 sin(p)
           + (3*cost*cost - 1) / pow(r,3)           // P_2^0
           + sint*sint*cos(2*phi) / pow(r,3)        // P_2^2 cos(2p)
           + sint*sint*sin(2*phi) / pow(r,3)  ;      // P_2^2 sin(2p)
         //  + sint*(15*cost*cost-3)*cosp / pow(r,4)      // P_3^1 cos(p)
         //  + sint*(15*cost*cost-3)*sinp / pow(r,4)      // P_3^1 sin(p)
         //  + pow(sint,3) * cos(3*phi) / pow(r,4)        // P_3^3 cos(3p)
         //  + pow(sint,3) * sin(3*phi) / pow(r,4) ;      // P_3^3 sin(3p)
           
    pot_hh.set_domain(nzm1) = potced(nzm1) ; 
    pot_hh.std_spectral_base() ; 
    
    pot_hh = relativistic_init * pot_hh ; 
    
    maxabs(pot_hh.laplacian(), "Laplacian of potential") ; 
    
    hh = pot_hh.derive_con(ff).derive_con(ff) ; 
    hh.dec_dzpuis(3) ; 

    // hh.set_etat_zero() ;    //  initialization to zero

    // Set up of field Q = Psi^2 N
    // ---------------------------
    
    Scalar qq(map) ; 
    
    qq = 1. + relativistic_init * r*r ; 
    potced = 1. + relativistic_init / (r*r) ; 
    qq.set_domain(nzm1) = potced(nzm1) ; 
     
    qq.std_spectral_base() ;    // sets standard spectral bases


    // Set up of conformal metric gamma_tilde
    // --------------------------------------
    
    Metric tgam( ff.con() ) ;   // construction from the flat metric

    tgam = ff.con() + hh ;      // initialization  [ Eq. (51) ]
    
    // Some shortcuts:
    const Sym_tensor& tgam_dd = tgam.cov() ;    
    const Sym_tensor& tgam_uu = tgam.con() ;    
    const Tensor_sym& delta = tgam.connect().get_delta() ;   // Eq. (26)


    // Set up of shift vector beta
    // ---------------------------    

    Vector beta(map, CON, otriad ) ;   
    beta = relativistic_init * pot_hh.derive_con(ff) ;   

    // Set up of lapse function N
    // --------------------------
    
    Scalar nn(map) ; 

    nn = 1. - relativistic_init * r*r ; 
    potced = 1. - relativistic_init / (r*r) ; 
    nn.set_domain(nzm1) = potced(nzm1) ; 

    nn.std_spectral_base() ;    // sets standard spectral bases

    // Set up of conformal factor Psi
    // ------------------------------
    
    Scalar psi(map) ;           // Psi
    Scalar psi2(map) ;          // Psi^2
    Scalar psi4(map) ;          // Psi^4
    Scalar ln_psi(map) ;           // ln(Psi)
    
    psi = sqrt(qq / nn) ;   
    psi.std_spectral_base() ;    // sets standard spectral bases
    psi2 = psi * psi ; 
    psi4 = psi2 * psi2 ; 
    ln_psi = log( psi ) ;  
    ln_psi.std_spectral_base() ;    // sets standard spectral bases

    // Set up of conformal extrinsic curvature A
    // -----------------------------------------
    
    Sym_tensor aa(map, CON, otriad) ;   // A^{ij}
    aa.set_etat_zero() ;    //  initialization to zero
    
    Sym_tensor taa(map, COV, otriad) ;  // {\tilde A}_{ij}
    taa = aa.up_down(ff) ;
    
    // Working stuff
    // -------------
    
    Scalar tmp(map) ; 
    Scalar tmp0(map) ; 

    //======================================================================
    //                  Start of time evolution
    //======================================================================
    
    int jmax = 1 ; 
    
    for (int jtime = 0; jtime <= jmax; jtime++) {

        // Source for Q  [ Eq. (76) ]
        // ------------
        
        Scalar aa_quad = contract(taa, 0, 1, aa, 0, 1) ; 
        
        Scalar source_qq = 0.75 * psi4 * qq * aa_quad  
            - contract( hh, 0, 1, qq.derive_cov(ff).derive_cov(ff), 0, 1 ) ; 
            
        source_qq.inc_dzpuis() ;  
            
        tmp = 0.0625 * contract( hh.derive_cov(ff), 0, 1, 
                        tgam_dd.derive_cov(ff), 0, 1 ).trace(tgam) 
             - 0.125 * contract( hh.derive_cov(ff), 0, 1,      
                        tgam_dd.derive_cov(ff), 0, 2 ).trace(tgam) 
             + 2.* contract( contract( tgam_uu, 0, ln_psi.derive_cov(ff), 0), 0,
                                 ln_psi.derive_cov(ff), 0 ) ;
     
        tmp0 = 2. * contract( tgam_uu, 0, 1, 
                              ln_psi.derive_cov(ff) * nn.derive_cov(ff), 0, 1) ;
        
        source_qq += psi2 * ( nn * tmp + tmp0 ) ; 
                             
        source_qq.spectral_display("source_qq") ; 


        // Source for N  [ Eq. (80) ]
        // ------------
        
        Scalar source_nn = psi4 * nn * aa_quad ;
        
        tmp = contract( hh, 0, 1, nn.derive_cov(ff).derive_cov(ff), 0, 1 ) ;
        tmp.inc_dzpuis() ; 
        
        source_nn -= tmp + tmp0 ; 
                    
        source_nn.spectral_display("source_nn") ; 

        // Source for beta [ Eq. (79) ]
        // ---------------

        Vector source_beta = 2. * ( contract(aa, 1, 
                        nn.derive_cov(ff) - 6.*nn * ln_psi.derive_cov(ff), 0)
                    - nn * contract(delta, 1, 2, aa, 0, 1) )
                - contract(hh, 0, 1, 
                           beta.derive_cov(ff).derive_cov(ff), 1, 2)
                - 0.3333333333333333 *
                  contract(hh, 1, beta.divergence(ff).derive_cov(ff), 0) ; 
                    
        source_beta.spectral_display("source_beta") ; 
        
        // Source for h [ Eq. (93) ]
        // -------------
                 
        Sym_tensor_trans source_hh(map, otriad, ff) ;  
        
        source_hh = (1. - nn*nn/psi4) * hh.derive_con(ff).divergence(ff) ;
        
        source_hh.inc_dzpuis() ; 
        
        // Computation of D^i h^{j k} + D^j h^{ik} - D^k h^{ij}
        Tensor_sym tstmp(map, 3, CON, otriad, 0, 1) ; 
        for (int i=1; i<=3; i++) {
            for (int j=1; j<=i; j++) {
                for (int k=1; k<=3; k++) {
                    tstmp.set(i,j,k) = hh.derive_con(ff)(j,k,i)
                                     + hh.derive_con(ff)(i,k,j) 
                                     - hh.derive_con(ff)(i,j,k) ;
                }
            }
        }
        

        source_hh += nn / psi4 * contract( tstmp, 2, nn.derive_cov(ff)
                                    + 2.*nn* ln_psi.derive_cov(ff), 0) ; 
       
        Sym_tensor tmpsym = beta * beta ; 
        source_hh += contract(tmpsym, 0, 1, 
                              hh.derive_cov(ff).derive_cov(ff), 2, 3) ; 
        
        
        cout << "Next step ?" << endl ; 
        arrete() ;  
    }


    return EXIT_SUCCESS ; 
}
