/*
 *  Test of Metric and Connection classes through the Schwarzschild metric
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

char test_schwarzschild_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/01/27 14:46:22  e_gourgoulhon
 * First version of test_schwarzschild.
 *
 *
 * $Header$
 *
 */

// C++ headers
#include <headcpp.h>

// C headers
#include <stdlib.h>
#include <math.h>

// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "graphique.h"
#include "utilitaires.h"
#include "proto.h"

int main() {

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 3 ; 	// Number of domains
    int nr = 17 ; 	// Number of collocation points in r in each domain
    int nt = 5 ; 	// Number of collocation points in theta in each domain
    int np = 4 ; 	// Number of collocation points in phi in each domain
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
    assert( nz == 3 ) ;  // since the above array describes only 3 domains
  
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    cout << map << endl ;  
    
    // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = map.r ;        // r field 
    
    // Schwarzschild metric in isotropic coordinates
    // ----------------------------------------------
    
    // Total mass = twice the radius at the horizon 
    double mm = 2. * map.val_r(0, 1., 0., 0.) ; 
        
    // Conformal factor
    Scalar psi4(map) ; 
    psi4 = pow( 1 + mm / (2*r), 4) ; 
    psi4.set_domain(0) = 1 ; 
    psi4.std_spectral_base() ;

    
    // Spatial metric

    const Metric_flat& fmet = map.flat_met_spher() ; 
    
    Sym_tensor gij_spher = psi4 * fmet.cov() ; 
    
    Metric gam(gij_spher) ; // construction from the covariant components 
   
    cout << gam << endl ; 
    arrete() ; 
    


    // Test: covariant derivative of the metric / flat metric:
    
    const Tensor& dg_cov = gam.cov().derive_cov( fmet ) ; 
    
    Vector dpsi4 = psi4.derive_cov( fmet ) ; 
    Tensor_sym d_gij_schw = fmet.cov() * dpsi4 ;
    
    Tensor diff_dg = dg_cov - d_gij_schw ; 
    cout << "Error on the covariant derivative of the metric / flat metric:" << endl ; 
    // diff_dg.spectral_display() ; 
    maxabs(diff_dg) ; 
    arrete() ; 

    // Test: Connection symbols Delta
    // ------------------------------
    const Tensor_sym& delta =  gam.connect().get_delta() ;     
    cout << "Connection (delta) : " << endl ; 
    delta.spectral_display() ; 
    maxabs(delta) ; 

    Scalar diff = delta(1,1,1) - 0.5 * dpsi4(1) / psi4 ; 
    cout << "Error on Delta^r_rr: \n " ; maxabs(diff) ; 

    diff = delta(1,2,1) - 0 ; 
    cout << "Error on Delta^r_rt: \n " ; maxabs(diff) ; 
    diff = delta(1,3,1) - 0 ; 
    cout << "Error on Delta^r_rp: \n " ; maxabs(diff) ; 
    
    diff = delta(1,2,2) + 0.5 * dpsi4(1) / psi4 ; 
    cout << "Error on Delta^r_tt: \n " ; maxabs(diff) ; 

    diff = delta(1,3,2) - 0  ; 
    cout << "Error on Delta^r_tp: \n " ; maxabs(diff) ; 

    diff = delta(2,1,1) - 0 ; 
    cout << "Error on Delta^t_rr: \n " ; maxabs(diff) ; 

    diff = delta(2,2,1) - 0.5 * dpsi4(1) / psi4 ; 
    cout << "Error on Delta^t_rt: \n " ; maxabs(diff) ; 

    diff = delta(2,3,1) - 0 ; 
    cout << "Error on Delta^t_rp: \n " ; maxabs(diff) ; 

    diff = delta(2,2,2) - 0 ; 
    cout << "Error on Delta^t_tt: \n " ; maxabs(diff) ; 

    diff = delta(2,3,2) - 0 ; 
    cout << "Error on Delta^t_tp: \n " ; maxabs(diff) ; 

    diff = delta(2,3,3) - 0 ; 
    cout << "Error on Delta^t_pp: \n " ; maxabs(diff) ; 

    diff = delta(3,1,1) - 0 ; 
    cout << "Error on Delta^p_rr: \n " ; maxabs(diff) ; 

    diff = delta(3,2,1) - 0 ; 
    cout << "Error on Delta^p_rt: \n " ; maxabs(diff) ; 

    diff = delta(3,3,1) - 0.5 * dpsi4(1) / psi4 ; 
    cout << "Error on Delta^p_rp: \n " ; maxabs(diff) ; 

    arrete() ; 
    
    // Test: covariant derivative of the metric / itself:
    // --------------------------------------------------
    
    const Tensor& dg_auto = gam.cov().derive_cov( gam ) ; 
    
    cout << "Error on the covariant derivative of the metric / itself:" << endl ; 
    maxabs(dg_auto) ; 
    arrete() ; 

    // Lapse function
    // --------------

    Scalar nn(map) ;
    
    nn = (1 - mm / (2*r) ) / (1 + mm / (2*r) ) ;
    
    nn.set_domain(0) = 1 ; 
    
    nn.std_spectral_base() ;
    
    cout << "Lapse N : " << nn << endl ; 
    arrete() ; 

    // Hamiltonian constraint
    // ----------------------
    
    const Tensor& ricci = gam.ricci() ; 
    
    const Scalar& ricci_scal = gam.ricci_scal() ; 
    
    cout << "Ricci scalar : " << endl ; 
    ricci_scal.spectral_display() ; 
    maxabs(ricci_scal) ; 

    // Dynamical Einstein equations
    //-----------------------------
    
    Sym_tensor dyn1 = - (nn.derive_cov(gam)).derive_cov(gam) ;
    
    Sym_tensor dyn2 = nn * ricci ; 
        
    Sym_tensor dyn_einstein = dyn1 + dyn2 ; 
    
    cout << "Dynamical Einstein equations:" << endl ; 
    dyn_einstein.spectral_display() ;     
    maxabs(dyn_einstein) ; 
    
    cout << "maxabs(dyn1) : " << endl ; 
    maxabs(dyn1) ; 
    cout << "maxabs(dyn2) : " << endl ; 
    maxabs(dyn2) ; 
    
    cout << "Relative error:" << endl ; 
    diffrel(dyn2, -dyn1) ; 

    return EXIT_SUCCESS ; 

}
