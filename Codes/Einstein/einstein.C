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

int main() {


    //======================================================================
    //      Construction and initialization of the various objects
    //======================================================================

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 3 ; 	// Number of domains
    int nr = 7 ; 	// Number of collocation points in r in each domain
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
    
    // Set up of tensor h
    // ------------------
    
    Sym_tensor_trans hh(map, otriad, ff) ;  // hh is a transverse tensor
                                            // with respect to the flat metric
                                            // thanks to Dirac gauge
                                            
    hh.set_etat_zero() ;    //  initialization to zero

    // Set up of field Q = Psi^2 N
    // ---------------------------
    
    Scalar qq(map) ; 
    qq = 1 ; 
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
    beta.set_etat_zero() ;    //  initialization to zero

    // Set up of lapse function N
    // --------------------------
    
    Scalar nn(map) ; 
    nn = 1 ;    // initial value

    // Set up of conformal factor Psi
    // ------------------------------
    
    Scalar psi(map) ;           // Psi
    Scalar psi2(map) ;          // Psi^2
    Scalar psi4(map) ;          // Psi^4
    Scalar ln_psi(map) ;           // ln(Psi)
    
    psi = sqrt(qq / nn) ;   
    psi2 = psi * psi ; 
    psi4 = psi2 * psi2 ; 
    ln_psi = log( psi ) ;  

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
        
        Scalar source_nn = psi4 * nn * aa_quad 
              - contract( hh, 0, 1, nn.derive_cov(ff).derive_cov(ff), 0, 1 ) 
              - tmp0 ;  
                    
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
                 
    }


    return EXIT_SUCCESS ; 
}
