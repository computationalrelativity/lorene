/*
 *  Test of Metric and Connection classes through the Kerr metric
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
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

char test_kerr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/01/04 21:03:12  e_gourgoulhon
 * Still some improvements...
 * More to come...
 *
 * Revision 1.1  2003/12/30 23:11:28  e_gourgoulhon
 * First version: not ready yet.
 *
 *
 * $Header$
 *
 */

// C++ headers
#include <headcpp.h>

// C headers
#include <stdlib.h>

// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"

int main() {

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 3 ; 	// Number of domains
    int nr = 7 ; 	// Number of collocation points in r in each domain
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
    const Coord& cost = map.cost ;  // cos(theta) field
    const Coord& sint = map.sint ;  // sin(theta) field
    
    // Kerr metric in quasi-isotropic coordinates
    // ------------------------------------------
    
    // Twice the radius at the horizon:
    double hh = 2 * map.val_r(0, 1., 0., 0.) ; 
    
    double aasm = 0.5 ; // Kerr parameter a/M = J/M^2
    
    double mm = hh / sqrt(1. - aasm*aasm) ; // total mass M
    
    double aa = aasm * mm ;  // angular momentum parameter a = J / M
    

    // R (Boyer-Lindquist radial coordinate):
    Mtbl rs = r + mm + (mm*mm - aa*aa) / (4*r) ;

    Mtbl erre = 1. + mm/r + (mm*mm - aa*aa) / (4*r*r) ; // ratio R / r :
    
    Mtbl sigmasr = rs + (aa*aa* cost*cost)/rs ;  // Sigma / R
    
    Scalar a2(map) ; Scalar b2(map) ;
    a2 = erre*erre + (aa*aa*cost*cost)/(r*r) ;  // A^2 = Sigma^2 / r^2
    
    b2 = erre*erre + aa*aa/(r*r) + 2*aa*aa*sint*sint*mm/(sigmasr*r*r) ; // B^2
    
    Scalar tmp(map) ; 
    tmp = 1 ; 
    a2.set_domain(0) = tmp.domain(0) ; 
    b2.set_domain(0) = tmp.domain(0) ; 
    a2.std_spectral_base() ;
    b2.std_spectral_base() ;


    // Spatial metric
    // --------------
    
    Sym_tensor gij_spher(map, COV,  map.get_bvect_spher()) ; 
    gij_spher.set(1,1) = a2 ; 
    gij_spher.set(1,2) = 0 ; 
    gij_spher.set(1,3) = 0 ; 
    gij_spher.set(2,2) = a2 ; 
    gij_spher.set(2,3) = 0 ; 
    gij_spher.set(3,3) = b2 ;
    
    Metric gam(gij_spher) ; // construction from the covariant components 
   
    cout << gam << endl ; 
    
    cout << "Contravariant components of the metric: " << endl ; 
    cout << gam.con() << endl ; 
    arrete() ; 

    // Lapse function
    // --------------
    
    
    tmp = 1 - 2*mm / sigmasr + 4*aa*aa*mm*mm* sint* sint / 
            ( sigmasr* (sigmasr * ( rs*rs + aa*aa)
              + 2*aa*aa*mm* sint*sint ) ) ; 
    tmp.annule_domain(0) ; 
    
    
    Scalar nn = sqrt(abs(tmp)) ; 
    
    tmp = 1 ; 
    nn.set_domain(0) = tmp.domain(0) ; 
    nn.std_spectral_base() ;
    
    cout << "Lapse N : " << nn << endl ; 
    arrete() ; 

    // Shift vector
    // ------------
    
    Scalar nphi(map) ; 
    nphi = ( 2*aa*mm / (sigmasr * ( erre*rs + aa*aa/r) 
                                + 2*aa*aa*mm* sint*sint/r ) ) * sint ; 
    nphi.annule_domain(0) ; 
                                    
    Vector beta(map, CON, map.get_bvect_spher() ) ;     
    beta.set(1) = 0 ; 
    beta.set(2) = 0 ; 
    beta.set(3) = nphi ; 
    beta.std_spectral_base() ;

    cout << "Shift vector beta: " << beta << endl ; 
    arrete() ; 
        
    // Extrinsic curvature
    // -------------------
    
    const Tensor_sym& delta =  gam.connect().get_delta() ; 
    
    cout << "Connection (delta) : " << endl ; 
    delta.spectral_display() ; 
    arrete() ; 
    
    Vector beta_cov = beta.down(0, gam) ; 
    cout << "beta_cov : " << endl ; 
    maxabs(beta_cov) ; 
    arrete() ; 
    
    cout << "Dbeta: " << endl ; 
    beta_cov.derive_cov(gam).spectral_display() ; 
    maxabs(beta_cov.derive_cov(gam)) ; 
    arrete() ; 
    
    Sym_tensor kk(map, COV, map.get_bvect_spher()) ; 
    
    for (int i=1; i<=3; i++) {
        for (int j=1; j<=i; j++) {
            kk.set(i,j) = 0.5 * ( beta_cov.derive_cov(gam)(i,j)
                + beta_cov.derive_cov(gam)(j,i) ) / nn ; 
        }
    } 
    
    kk.dec_dzpuis(2) ; 
    
    cout << "Extrinsic curvature : " << endl ;  
    kk.spectral_display() ; 
    maxabs(kk) ; 
    
    Tensor kkup = kk.up(1, gam) ; 
    Scalar trkk = kkup.scontract(0,1) ; 
    cout << "Trace of K:" << endl ; 
    trkk.spectral_display() ;
    maxabs(trkk) ; 
    arrete() ; 
    
    const Tensor& dkk = kkup.derive_cov(gam) ; 
    
    Vector mom_constr = dkk.scontract(1,2) - trkk.derive_cov(gam) ; 
    cout << "Momentum constraint : " << endl ; 
    mom_constr.spectral_display() ; 
    
    maxabs(mom_constr) ; 


    return EXIT_SUCCESS ; 

}
