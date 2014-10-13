/*
 *  Simple tensor manipulations
 *
 */

/*
 *   Copyright (c) 2005 Eric Gourgoulhon
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
 * Revision 1.6  2014/10/13 08:54:07  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:09:48  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2006/05/30 13:22:20  j_novak
 * Back to previous versions...
 *
 * Revision 1.2  2005/11/16 17:11:59  e_gourgoulhon
 * Added fij.std_spectral_base() (although not necessary in the
 *  present case)
 * Modified some comments.
 *
 * Revision 1.1  2005/11/16 09:46:30  e_gourgoulhon
 * Added demo_tensor.C
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>
#include <cassert>
#include <cmath>

// Lorene headers
#include "headcpp.h"    // standard input/output C++ headers (iostream, fstream)
#include "metric.h"     // classes Metric, Tensor, etc...     
#include "nbr_spx.h"    // defines infinity as an ordinary number: __infinity
#include "graphique.h"  // for graphical outputs
#include "utilitaires.h"    // utilities 

using namespace Lorene ;

int main() {

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 3 ; 	// Number of domains
    int nr = 17 ; 	// Number of collocation points in r in each domain
    int nt = 9 ; 	// Number of collocation points in theta in each domain
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
    double r_limits[] = {0., 1., 2., __infinity} ; 
    assert( nz == 3 ) ;  // since the above array describes only 3 domains
  
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    cout << map << endl ;  
    
    // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = map.r ;      
    const Coord& x = map.x ;      
    const Coord& y = map.y ;      
    
    // Some scalar field to be used as a conformal factor
    // --------------------------------------------------
    
    Scalar psi4(map) ; 
    
    psi4 = 1 + 5*x*y*exp(-r*r) ; 
    
    psi4.set_outer_boundary(nz-1, 1.) ;     // 1 at spatial infinity
                                            // (instead of NaN !)
    
    psi4.std_spectral_base() ;  // Standar polynomial bases will be used 
                                // to perform the spectral expansions

    // Graphical outputs:
    // -----------------
    
    // 1D view via PGPLOT
    des_profile(psi4, 0., 4., 1, M_PI/4, M_PI/4, "r", "\\gq\\u4") ;

    // 2D view of the slice z=0 via PGPLOT
    des_coupe_z(psi4, 0., -3., 3., -3., 3., "\\gq\\u4") ;

    // 3D view of the same slice via OpenDX
    psi4.visu_section('z', 0., -3., 3., -3., 3.) ; 
    
    cout << "Coefficients of the spectral expansion of Psi^4:" << endl ; 
    psi4.spectral_display() ; 
    
    arrete() ; // pause (waiting for return) 
    
    // Components of the flat metric in an orthonormal
    // spherical frame
    
    Sym_tensor fij(map, COV, map.get_bvect_spher()) ; 
    fij.set(1,1) = 1 ; 
    fij.set(1,2) = 0 ; 
    fij.set(1,3) = 0 ; 
    fij.set(2,2) = 1 ; 
    fij.set(2,3) = 0 ; 
    fij.set(3,3) = 1 ; 

    fij.std_spectral_base() ;  // Standar polynomial bases will be used 
                                // to perform the spectral expansions

     
    // Components of the physical metric in an orthonormal
    // spherical frame

    Sym_tensor gij = psi4 * fij ; 

    // Construction of the metric from the covariant components:     
    
    Metric gam(gij) ; 
    
    // Construction of a Vector : V^i = (Psi^4)^{;i}
    
    Vector vv = psi4.derive_con(gam) ;  // this is spherical components
                                        // (same triad as gam)

    vv.dec_dzpuis(2) ; // the dzpuis flag (power of r in the CED)
                       // is set to 0 (= 2 - 2) 
        
    // Cartesian components of the vector :
    Vector vv_cart = vv ; 
    vv_cart.change_triad( map.get_bvect_cart() ) ; 
    
    // Plot of the vector field :
    
    des_coupe_vect_z(vv_cart, 0., -4., 1., -2., 2., -2., 2., "Vector V") ;

    // A symmetric tensor of valence 2 : the Ricci tensor associated
    // with the metric gam :
    
    Sym_tensor tens1 = gam.ricci() ; 

    const Sym_tensor& tens2 = gam.ricci() ; // same as before except that 
                                            // no memory is allocated for a
                                            //    new tensor: tens2 is merely 
                                            //    a non-modifiable reference to
                                            //    the Ricci tensor of gam
                                            
    // Plot of tens1
    
    des_meridian(tens1, 0., 4., "Ricci (x r\\u3\\d in last domain)", 10) ; 

    // Another valence 2 tensor : the covariant derivative of V with respect to 
    //       the metric gam :
    
    Tensor tens3 = vv.derive_cov(gam) ; 
    
    const Tensor& tens4 = vv.derive_cov(gam) ; 
    
    // the reference tens4 is preferable over the new object tens3 if you do
    // not intend to modify tens4 or vv, because it does not perform any
    // memory allocation for a tensor. 

    // Raising an index with the metric gam :
    
    Tensor tens5 = tens3.up(1, gam) ; // 1 = second index 
                                      // (index j in the covariant derivative 
                                      //   V^i_{;j})
                  
    Tensor diff1 = tens5 - vv.derive_con(gam) ; // this should be zero

    // Check:
    cout << "Maximum value of diff1 in each domain : " << endl ; 
    Tbl tdiff1 =  max(diff1) ;  
      
    // Another valence 2 tensor : Lie_V R_{ij}
    
    Sym_tensor tens6 = tens1.derive_lie(vv) ; 
    
    // Contracting two tensors :
    
    Tensor tens7 = contract(tens1, 1, tens5, 0) ; // contracting the last index
                                                  // of tens1 with the first one 
                                                  // of tens5
    
    // self contraction of a tensor :
    
    Scalar scal1 = contract(tens3, 0, 1) ; // 0 = first index, 1 = second index
    
    // Each of these fields should be zero:
    
    Scalar diff2 = scal1 - vv.divergence(gam) ;   // divergence
    
    Scalar diff3 = scal1 - tens3.trace() ;        // trace                                      

    // Check : 
    cout << "Maximum value of diff2 in each domain : " 
        << max(abs(diff2)) << endl ;  
    
    cout << "Maximum value of diff3 in each domain : " 
        << max(abs(diff3)) << endl ;  
        
    arrete() ; 
    
    // Tensorial product : 
    
    Tensor_sym tens8 = tens1 * tens3 ;  // tens1 = R_{ij}
                                        // tens3 = V^k_{;l}
                                        // tens8 = (T8)_{ij}^k_l = R_ij V^k_{;l}
    
    cout << "Valence of tens8 : " << tens8.get_valence() << endl ; 
    
    arrete() ; 
    
    cout << "Spectral coefficients of the component (2,3,1,1) of tens8 : "
         << endl ; 
         
    tens8(2,3,1,1).spectral_display() ; 
    
    ////////////////////////////////////////////////////////////////////////
    //                                                                    //
    // To see more functions, please have a look to Lorene documentation  // 
    // at http://www.lorene.obspm.fr/Refguide/                            //
    //                                                                    //
    ////////////////////////////////////////////////////////////////////////
                                               
    return EXIT_SUCCESS ; 

}
