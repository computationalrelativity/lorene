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
 * Revision 1.5  2004/01/23 13:28:14  e_gourgoulhon
 * Scalar::set_val_hor --> Scalar::set_inner_boundary.
 *
 * Revision 1.4  2004/01/23 08:01:36  e_gourgoulhon
 * All Einstein equations are now verified.
 *
 * Revision 1.3  2004/01/19 16:58:49  e_gourgoulhon
 * Momemtum constraint OK !
 * Next step: Hamiltonian constraint.
 *
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
#include <math.h>

// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "graphique.h"
#include "utilitaires.h"
#include "cmp.h"
#include "proto.h"

int main() {

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 3 ; 	// Number of domains
    int nzm1 = nz - 1 ; // Index of outermost domain
    int nr = 17 ; 	// Number of collocation points in r in each domain
    int nt = 9 ; 	// Number of collocation points in theta in each domain
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
    double hh = 2. * map.val_r(0, 1., 0., 0.) ; 
    
    double aasm = 0.5 ; // Kerr parameter a/M = J/M^2
    
    double mm = hh / sqrt(1. - aasm*aasm) ; // total mass M
    
    double aa = aasm * mm ;  // angular momentum parameter a = J / M
    

    // R (Boyer-Lindquist radial coordinate):
    Mtbl rs = r + mm + (mm*mm - aa*aa) / (4*r) ;

    Mtbl erre = 1 + mm/r + (mm*mm - aa*aa) / (4*r*r) ; // ratio R / r :
    
    Mtbl sigmasr = rs + (aa*aa* cost*cost)/rs ;  // Sigma / R
    
    Scalar a2(map) ; Scalar b2(map) ;
    a2 = erre*erre + (aa*aa*cost*cost)/(r*r) ;  // A^2 = Sigma^2 / r^2
    
    b2 = erre*erre + aa*aa/(r*r) + 2*aa*aa*sint*sint*mm/(sigmasr*r*r) ; // B^2
    
    a2.set_domain(0) = 1. ; 
    b2.set_domain(0) = 1. ; 
    a2.std_spectral_base() ;
    b2.std_spectral_base() ;

    /*
    des_coef_xi(a2.get_spectral_va(), 1, 0, 0, 1.e-14, "a2") ; 
    des_coef_xi(b2.get_spectral_va(), 1, 0, 0, 1.e-14, "b2") ; 
    des_coef_theta(a2.get_spectral_va(), 1, 0, 0, 1.e-14, "a2") ; 
    des_coef_theta(b2.get_spectral_va(), 1, 0, 0, 1.e-14, "b2") ; 
    */
    
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
    
//    cout << "Contravariant components of the metric: " << endl ; 
//    cout << gam.con() << endl ; 

    // Flat metric :
    // -----------
    const Metric_flat& fmet = map.flat_met_spher() ; 


    // Schwartzchild metric
    // --------------------
    
    Scalar psi4(map) ; 
    psi4 = pow( 1 + mm / (2*r), 4) ; 
    psi4.set_domain(0) = 1 ; 
    psi4.std_spectral_base() ;

    Vector dpsi4 = psi4.derive_cov( fmet ) ; 
    Scalar tmp(map) ; 
    tmp = - 2 * mm / (r*r) * pow( 1 + mm / (2*r), 3) ;
    tmp.set_domain(0) = 0 ;  
    Vector vtmp(map, COV, map.get_bvect_spher()) ; 
    vtmp.set(1) = tmp ;
    tmp = - 2 * mm * pow( 1 + mm / (2*r), 3) ; 
    vtmp.set(1).set_domain(nzm1) = tmp.domain(nzm1) ; 
    vtmp.set(1).set_dzpuis(2) ; 
    vtmp.set(2) = 0 ; 
    vtmp.set(3) = 0 ; 
    vtmp.std_spectral_base() ;
    vtmp -= dpsi4 ; 
    cout << "Error on Grad(Psi^4) : " << endl ; 
    vtmp.spectral_display() ; 
    maxabs(vtmp) ; 
    arrete() ; 
    
    Sym_tensor gij_schw = psi4 * fmet.cov() ; 
    
    Sym_tensor diff_schw = gam.cov() - gij_schw ; 
    cout << "Comparison with Schwarzschild: \n" ; 
    maxabs(diff_schw) ; 
    arrete() ; 

    // Test: covariant derivative of the metric / flat metric:
    
    const Tensor& dg_cov = gam.cov().derive_cov( fmet ) ; 
    
    Tensor_sym d_gij_schw = fmet.cov() * dpsi4 ;
    
    Tensor diff_dg = dg_cov - d_gij_schw ; 
    cout << "Error on the covariant derivative of the metric / flat metric:" << endl ; 
    diff_dg.spectral_display() ; 
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
    
    
    tmp = 1 - 2*mm / sigmasr 
            + 4*aa*aa*mm*mm* sint* sint / 
            ( sigmasr* (sigmasr * ( rs*rs + aa*aa)
              + 2*aa*aa*mm* sint*sint ) ) ; 
    tmp.set_domain(0) = 1 ; 
    
    Scalar nn = sqrt(abs(tmp)) ; 
    nn.set_inner_boundary(1, 0.) ; 

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

    /*
	des_coef_xi(beta(3).get_spectral_va(), 1, 0, 0, 1.e-14, "beta^phi") ; 
    des_coef_xi(beta(3).get_spectral_va(), 1, 0, 1, 1.e-14, "beta^phi") ; 
    des_coef_theta(beta(3).get_spectral_va(), 1, 0, 0, 1.e-14, "beta^phi") ; 
	*/

    cout << "Shift vector beta: " << beta << endl ; 
	arrete() ; 
	
    
        
    // Extrinsic curvature
    // -------------------
    
    Vector beta_cov = beta.down(0, gam) ; 
    cout << "beta_cov : " << endl ; 
    maxabs(beta_cov) ; 
    arrete() ; 
    
    cout << "Dbeta: " << endl ; 
    beta_cov.derive_cov(gam).spectral_display() ; 
    maxabs(beta_cov.derive_cov(gam)) ; 
    arrete() ; 
    
    // N / (xi+1):
    Scalar nn_xip1 = Scalar( division_xpun(Cmp(nn), 0) ) ; 
    
    Sym_tensor kk(map, COV, map.get_bvect_spher()) ; 
    
    for (int i=1; i<=3; i++) {
        for (int j=1; j<=i; j++) {

            tmp = 0.5 * ( beta_cov.derive_cov(gam)(i,j)
                        + beta_cov.derive_cov(gam)(j,i) ) ; 

            tmp.set_inner_boundary(1, 0.) ;

            cout << "######### Component " << i << " " << j << endl ; 
            for (int k=0; k<mgrid.get_np(1); k++) {
                for (int jj=0; jj<mgrid.get_nt(1); jj++) {
                    cout << "  " << tmp.get_spectral_va()(1,k,jj,0) ; 
                }
                cout << endl ; 
            }
            
           
            kk.set(i,j) = Scalar( division_xpun(Cmp(tmp), 0) ) / nn_xip1 ;  
            
   /*
    des_coef_xi(kk(i,j).get_spectral_va(), 1, 0, 0, 1.e-14, "K^ij") ; 
    des_coef_xi(kk(i,j).get_spectral_va(), 1, 0, 1, 1.e-14, "K^ij") ; 
    des_coef_theta(kk(i,j).get_spectral_va(), 1, 0, 0, 1.e-14, "K^ij") ; 
    */
  
        }
    } 
    
    cout << "Extrinsic curvature : " << endl ;  
    kk.spectral_display() ; 
    maxabs(kk) ; 
    
    // Momentum constraint
    // -------------------
    
    Tensor kk_du = kk.up(1, gam) ; 
    Scalar trkk = kk_du.scontract(0,1) ; 
    cout << "Trace of K:" << endl ; 
    trkk.spectral_display() ;
    maxabs(trkk) ; 
    arrete() ; 
    
    const Tensor& dkk = kk_du.derive_cov(gam) ; 
    
    Vector mom_constr = dkk.scontract(1,2) - trkk.derive_cov(gam) ; 
    cout << "Momentum constraint : " << endl ; 
    mom_constr.spectral_display() ;     
    maxabs(mom_constr) ; 

    // Hamiltonian constraint
    // ----------------------
    
    const Tensor& ricci = gam.ricci() ; 
    
    const Scalar& ricci_scal = gam.ricci_scal() ; 
    
    cout << "Ricci scalar : " << endl ; 
    ricci_scal.spectral_display() ; 
    maxabs(ricci_scal) ; 

    Tensor kk_uu = kk_du.up(0, gam) ; 
    
    tmp = trkk * trkk - contract( contract(kk, 1, kk_uu, 1), 0, 1) ; 
    tmp.dec_dzpuis() ; 
    Scalar ham_constr = ricci_scal + tmp ;
        
    cout << "Hamiltonian constraint : " << endl ; 
    ham_constr.spectral_display() ;     
    maxabs(ham_constr) ; 

   // ham_constr.visu_section('y', 0., 0., 5., 0., 5., "Ham. constr.") ; 
    
    // Dynamical Einstein equations
    //-----------------------------
    
    Sym_tensor dyn1 = - (nn.derive_cov(gam)).derive_cov(gam) ;
    
    Sym_tensor dyn2 = nn * ricci ; 
        
    Sym_tensor dyn3 = nn * ( trkk * kk - 2 * contract(kk_du, 1, kk, 0) ) ; 
    dyn3.dec_dzpuis() ; 
    
    // Lie derivative of K along beta:
    Sym_tensor dyn4 = contract( beta.derive_cov(gam), 0, kk, 0)
                + contract( kk, 0, beta.derive_cov(gam), 0) ;
    dyn4.dec_dzpuis() ; 
    dyn4 += contract(beta, 0, kk.derive_cov(gam), 2)  ;  

    Sym_tensor dyn_einstein = dyn1 + dyn2 + dyn3 + dyn4; 
    
    cout << "Dynamical Einstein equations:" << endl ; 
    dyn_einstein.spectral_display() ;     
    maxabs(dyn_einstein) ; 
    
    cout << "maxabs(dyn1) : " << endl ; 
    maxabs(dyn1) ; 
    cout << "maxabs(dyn2) : " << endl ; 
    maxabs(dyn2) ; 
    cout << "maxabs(dyn3) : " << endl ; 
    maxabs(dyn3) ; 
    cout << "maxabs(dyn4) : " << endl ; 
    maxabs(dyn4) ; 
    
    cout << "Relative error:" << endl ; 
    diffrel(dyn2+dyn3+dyn4, -dyn1) ; 

    return EXIT_SUCCESS ; 

}
