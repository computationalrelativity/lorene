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


    // Set up of conformal metric gamma_tilde
    // --------------------------------------
    
    Metric tgam( ff.con() ) ;   // construction from the flat metric
    tgam = ff.con() + hh ;  // initialization  
    
    const Sym_tensor& tgam_dd = tgam.cov() ;    // abreviation

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

    //======================================================================
    //                  Start of time evolution
    //======================================================================
    
    int jmax = 10 ; 
    
    for (int jtime = 0; jtime <= jmax; jtime++) {

        // Hamiltonian constraint combined with the lapse equation (Eq. for Q)
        // -------------------------------------------------------------------
        
        Scalar aa_quad = contract(taa, 0, 1, aa, 0, 1) ; 
        
        Scalar source_qq = 0.75 * psi4 * qq * aa_quad  
            - contract( hh, 0, 1, qq.derive_cov(ff).derive_cov(ff), 0, 1 ) ;  
            
        tmp = contract( hh.derive_cov(ff), 0, 1, 
                        tgam_dd.derive_cov(ff), 0, 1 ) ; 
                    

    }


    return EXIT_SUCCESS ; 
}
