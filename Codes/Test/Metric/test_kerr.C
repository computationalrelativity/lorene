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
 * Revision 1.9  2004/02/04 15:47:07  e_gourgoulhon
 * Changed declaration of ricci and ricci_c from "const Tensor&" to
 *  "const Sym_tensor&".
 *
 * Revision 1.8  2004/01/29 16:07:58  e_gourgoulhon
 * More polishing.
 *
 * Revision 1.7  2004/01/29 15:22:47  e_gourgoulhon
 * Introduced call to method Tensor::divergence in momentum constraint.
 * Polishing.
 *
 * Revision 1.6  2004/01/28 15:28:27  e_gourgoulhon
 * Added tests with Cartesian components.
 *
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

    cout << "Value of the Kerr parameter a/M ?" << endl ; 
    double aasm ; 
    cin >> aasm ; 

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 4 ; 	// Number of domains
    int nzm1 = nz - 1 ; // Index of outermost domain
    int nr = 25 ; 	// Number of collocation points in r in each domain
    int nt = 17 ; 	// Number of collocation points in theta in each domain
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
    double r_limits[] = {0., 2., 3., 6., __infinity} ; 
    assert( nz == 4 ) ;  // since the above array describes only 3 domains
  
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    cout << map << endl ;  
    
    // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = map.r ;        // r field 
    const Coord& cost = map.cost ;  // cos(theta) field
    const Coord& sint = map.sint ;  // sin(theta) field
    
    cout << "*********************************************************" << endl 
     << "                 Tests in spherical components"  << endl 
     <<  "*********************************************************"
     << endl ; 

    // Analytic form of the Kerr metric in quasi-isotropic coordinates
    // ---------------------------------------------------------------
    
    // Twice the radius at the horizon:
    double hh = 2. * map.val_r(0, 1., 0., 0.) ; 
    
    double mm = hh / sqrt(1. - aasm*aasm) ; // total mass M
    
    double aa = aasm * mm ;  // angular momentum parameter a = J / M
    

    // R (Boyer-Lindquist radial coordinate):
    Mtbl r_bl = r + mm + (mm*mm - aa*aa) / (4*r) ;

    Mtbl r_blsr = 1 + mm/r + (mm*mm - aa*aa) / (4*r*r) ; // ratio R / r :
    
    Mtbl sigmasr = r_bl + (aa*aa* cost*cost) / r_bl ;  // Sigma / R
    
    Scalar a2(map) ; Scalar b2(map) ;
    a2 = r_blsr * r_blsr + (aa*aa*cost*cost)/(r*r) ;  // A^2 = Sigma^2 / r^2
    
    b2 = r_blsr * r_blsr + aa*aa/(r*r) 
                    + 2*aa*aa*sint*sint*mm/(sigmasr*r*r) ; // B^2
    
    a2.set_domain(0) = 1. ; // Metric set to 1 in the nucleus 
    b2.set_domain(0) = 1. ; //   ("inside" the horizon)
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
   
    cout << 
    "Minimum val. of the covariant components of the metric in each domain: " 
    << endl ; 
    min(gam.cov()) ;
    cout << 
    "Maximum val. of the covariant components of the metric in each domain: " 
    << endl ; 
    max(gam.cov()) ;
    arrete() ; 
    
    Scalar tmp(map) ; // working scalar

    // Comparison with Schwarzschild metric in the non-rotating case
    // -------------------------------------------------------------
    if (aasm == double(0)) {    
    
        // Flat metric :
        const Metric_flat& fmet = map.flat_met_spher() ; 

        // Schwarzchild metric : 
    
        Scalar psi4(map) ; 
        psi4 = pow( 1 + mm / (2*r), 4) ; // conformal factor
        psi4.set_domain(0) = 1 ; 
        psi4.std_spectral_base() ;

        Vector dpsi4 = psi4.derive_cov( fmet ) ; 
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
        cout << "Error on Grad(Psi^4) (max absolute value in each domain) : \n " ; 
        maxabs(vtmp) ; 
        arrete() ; 
    
        Sym_tensor gij_schw = psi4 * fmet.cov() ; 
    
        Sym_tensor diff_schw = gam.cov() - gij_schw ; 
        cout << 
        "Comparison with the Schwarzschild metric (max absolute error): \n" ; 
        maxabs(diff_schw) ; 
        arrete() ; 

        // Test: covariant derivative of the metric / flat metric:

        const Tensor& dg_cov = gam.cov().derive_cov( fmet ) ; 
    
        Tensor_sym d_gij_schw = fmet.cov() * dpsi4 ;
    
        Tensor diff_dg = dg_cov - d_gij_schw ; 
        cout << "Error on the covariant derivative of the metric / flat metric:" << endl ; 
        diff_dg.spectral_display() ; 
        cout << "Error on the covariant derivative of the metric / flat metric:" << endl ; 
        maxabs(diff_dg) ; 
        arrete() ; 

    }

        
    // Test: covariant derivative of the metric / itself:
    // --------------------------------------------------
    
    const Tensor& dg_auto = gam.cov().derive_cov( gam ) ; 
    
    cout << "\n Error on the covariant derivative of the metric / itself \n" ;
    cout << "  (maximum absolute value in each domain) : \n" ; 
    maxabs(dg_auto) ; 
    arrete() ; 


    // Lapse function
    // --------------
    
    tmp = 1 - 2*mm / sigmasr 
            + 4*aa*aa*mm*mm* sint* sint / 
            ( sigmasr* (sigmasr * ( r_bl * r_bl + aa*aa)
              + 2*aa*aa*mm* sint*sint ) ) ; 
    tmp.set_domain(0) = 1 ; 
    
    Scalar nn = sqrt(abs(tmp)) ; 
    nn.set_inner_boundary(1, 0.) ; 

    nn.std_spectral_base() ;
    
    cout << "Minimum value of the lapse N in each domain : " << endl ; 
    min(nn, cout) ;
    cout << "Maximum value of the lapse N in each domain : " << endl ; 
    max(nn, cout) ;
    arrete() ;  

    // Shift vector
    // ------------
    
    Scalar nphi(map) ; 
    nphi = ( 2*aa*mm / (sigmasr * ( r_blsr * r_bl + aa*aa/r) 
                                + 2*aa*aa*mm* sint*sint/r ) ) * sint ; 
    nphi.annule_domain(0) ; 
                                    
    Vector beta(map, CON, map.get_bvect_spher() ) ;     
    beta.set(1) = 0 ; 
    beta.set(2) = 0 ; 
    beta.set(3) = - nphi ; 
    beta.std_spectral_base() ;

    cout << "Minimum value of the shift vector beta^i in each domain : \n" ; 
    min(beta) ;
    cout << "Maximum value of the shift vector beta^i in each domain : \n" ; 
    max(beta) ;
    cout << endl ; 

    Vector beta_cov = beta.down(0, gam) ; 
    cout << "Minimum value of the covariant comp. beta_i in each domain : \n" ; 
    min(beta_cov) ;
    cout << "Maximum value of the covariant comp. beta_i in each domain : \n" ; 
    max(beta_cov) ;
    arrete() ; 
        
    // Extrinsic curvature
    // -------------------
    
    cout << "D_j beta_i : " << endl ; 
    beta_cov.derive_cov(gam).spectral_display() ; 
    cout << "D_j beta_i (max. absolute value in each domain) : " << endl ; 
    maxabs(beta_cov.derive_cov(gam)) ; 
    arrete() ; 
    
    // Division of the lapse by (xi+1) (where xi is the grid radial 
    // variable) in the first shell to get a non-vanishing quantity
    // at the horizon (horizon <=> xi=-1) :
    
    Scalar nn_xip1 = Scalar( division_xpun(Cmp(nn), 0) ) ; 
    
    Sym_tensor kk(map, COV, map.get_bvect_spher()) ; 
    
    for (int i=1; i<=3; i++) {
        for (int j=1; j<=i; j++) {

            tmp = 0.5 * ( beta_cov.derive_cov(gam)(i,j)
                        + beta_cov.derive_cov(gam)(j,i) ) ; 

            cout << "## Horizon values of D_(i beta_j) component " << i 
                 << " " << j << endl ; 
            for (int k=0; k<mgrid.get_np(1); k++) {
                for (int jj=0; jj<mgrid.get_nt(1); jj++) {
                    cout << "  " << tmp.get_spectral_va()(1,k,jj,0) ; 
                }
                cout << endl ; 
            }
            
            tmp.set_inner_boundary(1, 0.) ;
           
            // Division by (xi+1) (where xi is the grid radial 
            // variable) in the first shell to get a non-vanishing quantity
            // at the horizon (horizon <=> xi=-1) 
            // and then division by N/(xi+1) to get K_{ij] : 
            
            kk.set(i,j) = Scalar( division_xpun(Cmp(tmp), 0) ) / nn_xip1 ;  
            
        }
    } 
    arrete() ; 
        
    cout << "Extrinsic curvature : " << endl ;  
    kk.spectral_display() ; 
    cout << 
    "Extrinsic curvature K_{ij} (max. absolute value in each domain) : \n" ;  
    maxabs(kk) ; 
    
    // Trace of K
    // ----------
    
    Tensor kk_du = kk.up(1, gam) ; 
    Scalar trkk = kk_du.scontract(0,1) ; 
    cout << "Trace of K  (max. absolute value in each domain) :" << endl ; 
    maxabs(trkk) ; 
    
    // Test:
    Scalar div_beta = beta.divergence(gam) ; 
    cout << "Divergence of beta  (max. absolute value in each domain) :\n" ; 
    maxabs(div_beta) ; 
    
    tmp = nn * trkk - div_beta ; 
    cout << "N K - div(beta)  (max. absolute value in each domain) : \n" ;
    maxabs(tmp) ; 

    arrete() ; 

    // Momentum constraint
    // -------------------
    
    Vector mom_constr = kk_du.divergence(gam) - trkk.derive_cov(gam) ; 
    
    cout << "Momentum constraint : " << endl ; 
    mom_constr.spectral_display() ;     
    cout << "Momentum constraint (max absolute error in each domain): " << endl ; 
    maxabs(mom_constr) ; 

    const Tensor& dkk = kk_du.derive_cov(gam) ; 
    mom_constr = dkk.scontract(1,2) - trkk.derive_cov(gam) ; 
    cout << "Momentum constraint by direct computation of divergence\n "
         << "   (max absolute error in each domain) : " << endl ; 
    maxabs(mom_constr) ; 


    arrete() ; 

    // Hamiltonian constraint
    // ----------------------
        
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
    cout << "Hamiltonian constraint (max absolute error in each domain): " << endl ; 
    maxabs(ham_constr) ; 
    arrete() ; 

   // ham_constr.visu_section('y', 0., 0., 5., 0., 5., "Ham. constr.") ; 
    
    // Dynamical Einstein equations
    //-----------------------------
    
    const Sym_tensor& ricci = gam.ricci() ; 

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
    cout << "Dynamical Einstein equations (max absolute error in each domain): " << endl ; 
    maxabs(dyn_einstein) ; 
    
    cout << "maxabs(dyn1) : " << endl ; 
    maxabs(dyn1) ; 
    cout << "maxabs(dyn2) : " << endl ; 
    maxabs(dyn2) ; 
    cout << "maxabs(dyn3) : " << endl ; 
    maxabs(dyn3) ; 
    cout << "maxabs(dyn4) : " << endl ; 
    maxabs(dyn4) ; 
    
    cout << "Dynamical Einstein equations (relative error in each domain): " << endl ; 
    diffrel(dyn2+dyn3+dyn4, -dyn1) ; 
    arrete() ; 
    
    //==============================================================//
    //                  Tests in Cartesian components               //
    //==============================================================//

    cout << "*********************************************************" << endl 
     << "                 Tests in Cartesian components"  << endl 
     <<  "*********************************************************"
     << endl ; 

    Sym_tensor gij_cart = gam.cov() ; 
    gij_cart.change_triad( map.get_bvect_cart() ) ; 
    
    Metric gam_c(gij_cart) ; // construction from the covariant components 
   
    cout << "Metric in Cartesian components : \n" ; 
    gam_c.cov().spectral_display() ; 
    arrete() ; 
    
    Vector beta_c = beta ; 
    beta_c.change_triad( map.get_bvect_cart() ) ; 

    cout << "Shift vector in Cartesian components : \n" ; 
    beta_c.spectral_display() ; 
    arrete() ; 
    
    // Extrinsic curvature in Cartesian components
    // -------------------------------------------
    
    Vector beta_c_cov = beta_c.down(0, gam_c) ; 
    cout << "beta_c_cov : " << endl ; 
    maxabs(beta_c_cov) ; 
    arrete() ; 
        
    Sym_tensor kk_c(map, COV, map.get_bvect_cart()) ; 
    
    for (int i=1; i<=3; i++) {
        for (int j=1; j<=i; j++) {

            tmp = 0.5 * ( beta_c_cov.derive_cov(gam_c)(i,j)
                        + beta_c_cov.derive_cov(gam_c)(j,i) ) ; 


            cout << "## Horizon values of D_(i beta_j) component " << i 
                 << " " << j << endl ; 

            for (int k=0; k<mgrid.get_np(1); k++) {
                for (int jj=0; jj<mgrid.get_nt(1); jj++) {
                    cout << "  " << tmp.get_spectral_va()(1,k,jj,0) ; 
                }
                cout << endl ; 
            }
            
            // tmp.set_inner_boundary(1, 0.) ;
          
            kk_c.set(i,j) = Scalar( division_xpun(Cmp(tmp), 0) ) / nn_xip1 ;  
            
        }
    } 
    arrete() ; 
    
    cout << "Extrinsic curvature (Cartesian components) : " << endl ;  
    kk_c.spectral_display() ; 
    cout << "Extrinsic curvature K_{ij} \n" <<
    "  (max. absolute value of Cartesian components in each domain) : \n" ;  
    maxabs(kk_c) ; 
    
    // Test :
    Sym_tensor kk_test = kk ; 
    kk_test.change_triad( map.get_bvect_cart() ) ;
    cout << "Relative error on the Cartesian components of K_{ij}: \n" ;
    diffrelmax(kk_c, kk_test) ;  
    arrete() ; 
    
    // Trace of K
    // ----------
    
    Tensor kk_du_c = kk_c.up(1, gam_c) ; 
    Scalar trkk_c = kk_du_c.scontract(0,1) ; 
    cout << "Trace of K:" << endl ; 
    maxabs(trkk_c) ; 

    // Test:
    Scalar div_beta_c = beta_c.divergence(gam_c) ; 
    cout << "Divergence of beta:\n" ; 
    maxabs(div_beta_c) ; 
    
    tmp = nn * trkk_c - div_beta_c ; 
    cout << "N K - div(beta): \n" ;
    maxabs(tmp) ; 

    arrete() ; 

    // Momentum constraint
    // -------------------
        
    Vector mom_constr_c = kk_du_c.divergence(gam_c) - trkk_c.derive_cov(gam_c) ; 

    cout << "Momentum constraint (Cart.) : " << endl ; 
    mom_constr_c.spectral_display() ;     
    cout << "Momentum constraint (Cart.) (max absolute error in each domain): " << endl ; 
    maxabs(mom_constr_c) ; 

    const Tensor& dkk_c = kk_du_c.derive_cov(gam_c) ; 
    mom_constr_c = dkk_c.scontract(1,2) - trkk_c.derive_cov(gam_c) ; 
    
    cout << "Momentum constraint (Cart.) by direct computation of divergence \n"
        << "  (max absolute error in each domain) : " << endl ; 
    maxabs(mom_constr_c) ; 
    arrete() ; 

    // Hamiltonian constraint
    // ----------------------
        
    const Scalar& ricci_scal_c = gam_c.ricci_scal() ; 
    
    cout << "Ricci scalar (Cart)  : " << endl ; 
    ricci_scal_c.spectral_display() ; 
    maxabs(ricci_scal_c) ; 

    Tensor kk_uu_c = kk_du_c.up(0, gam_c) ; 
    
    tmp = trkk_c * trkk_c - contract( contract(kk_c, 1, kk_uu_c, 1), 0, 1) ; 
    tmp.dec_dzpuis() ; 
    Scalar ham_constr_c = ricci_scal_c + tmp ;
        
    cout << "Hamiltonian constraint (Cart) : " << endl ; 
    ham_constr_c.spectral_display() ;     
    cout << "Hamiltonian constraint (max absolute error in each domain): " << endl ; 
    maxabs(ham_constr_c) ; 
    arrete() ; 

    // Dynamical Einstein equations
    //-----------------------------
    
    const Sym_tensor& ricci_c = gam_c.ricci() ; 
    
    // Test :
    Sym_tensor ricci_test = ricci ; 
    ricci_test.change_triad( map.get_bvect_cart() ) ;
    cout << "Relative error on the Cartesian components of Ricci: \n" ;
    diffrelmax(ricci_c, ricci_test) ;  
    arrete() ; 

    Sym_tensor dyn1_c = - (nn.derive_cov(gam_c)).derive_cov(gam_c) ;
    
    Sym_tensor dyn2_c = nn * ricci_c ; 
        
    Sym_tensor dyn3_c = nn * ( trkk_c * kk_c 
        - 2 * contract(kk_du_c, 1, kk_c, 0) ) ; 
    dyn3_c.dec_dzpuis() ; 
    
    // Lie derivative of K along beta:
    Sym_tensor dyn4_c = contract( beta_c.derive_cov(gam_c), 0, kk_c, 0)
                + contract( kk_c, 0, beta_c.derive_cov(gam_c), 0) ;
    dyn4_c.dec_dzpuis() ; 
    dyn4_c += contract(beta_c, 0, kk_c.derive_cov(gam_c), 2)  ;  

    Sym_tensor dyn_einstein_c = dyn1_c + dyn2_c + dyn3_c + dyn4_c ; 
    
    cout << "Dynamical Einstein equations:" << endl ; 
    dyn_einstein_c.spectral_display() ;     
    cout << "Dynamical Einstein equations (max absolute error in each domain): " << endl ; 
    maxabs(dyn_einstein_c) ; 
    
    cout << "maxabs(dyn1_c) : " << endl ; 
    maxabs(dyn1_c) ; 
    cout << "maxabs(dyn2_c) : " << endl ; 
    maxabs(dyn2_c) ; 
    cout << "maxabs(dyn3_c) : " << endl ; 
    maxabs(dyn3_c) ; 
    cout << "maxabs(dyn4_c) : " << endl ; 
    maxabs(dyn4_c) ; 
    
    cout << "Dynamical Einstein equations (relative error in each domain): " << endl ; 
    diffrel(dyn2_c+dyn3_c+dyn4_c, -dyn1_c) ; 
    
    return EXIT_SUCCESS ; 

}
