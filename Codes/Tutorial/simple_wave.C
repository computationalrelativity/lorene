/*
 *  Simple wave equation.
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon
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

char simple_wave_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.14  2005/03/25 20:39:27  e_gourgoulhon
 * Call to the new function des_coupe_z for Scalar's.
 *
 * Revision 1.13  2004/05/11 20:17:04  e_gourgoulhon
 * New prototype of des_evol.
 *
 * Revision 1.12  2004/04/06 08:26:21  j_novak
 * Update of the list of arguments for Evolution constructors.
 *
 * Revision 1.11  2004/02/25 16:44:57  j_novak
 * workflag is now allocated dynamically.
 *
 * Revision 1.10  2004/02/21 17:06:43  e_gourgoulhon
 * Method Scalar::point renamed Scalar::val_grid_point.
 *
 * Revision 1.9  2004/02/18 18:59:29  e_gourgoulhon
 * Suppressed unnecessary #include's.
 *
 * Revision 1.8  2004/02/17 22:20:49  e_gourgoulhon
 * Much better plots thanks to the new function des_profile_mult.
 * Added evolution of the central value of the field and its plot
 * through the new function des_evol.
 *
 * Revision 1.7  2004/02/17 09:19:08  j_novak
 * Cleaning of Param at the end of the code (memory leaks).
 *
 * Revision 1.6  2004/02/16 13:00:00  e_gourgoulhon
 * Added the C headers (not required by GNU g++ !!!).
 *
 * Revision 1.5  2004/02/15 22:08:16  e_gourgoulhon
 * The example with Poisson equation is now in file simple_poisson.C.
 * simple_wave.C contains now an example of resolution of d'Alembert
 * equation. The time evolution is managed thanks to the new
 * class Evolution_std.
 *
 * Revision 1.4  2003/12/16 06:33:31  e_gourgoulhon
 * Added call to method Scalar::visu_box.
 *
 * Revision 1.3  2003/12/14 21:53:26  e_gourgoulhon
 * Added 3D visualization of vector field through Vector::visu_arrows.
 *
 * Revision 1.2  2003/12/11 16:21:05  e_gourgoulhon
 * Use simplified version of Scalar::visu_section.
 *
 * Revision 1.1  2003/12/11 11:30:02  e_gourgoulhon
 * First version.
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include "stdlib.h"
#include "assert.h"
#include "math.h"

// Lorene headers
#include "metric.h"
#include "graphique.h"
#include "param.h"
#include "evolution.h"
#include "utilitaires.h"


int main() {

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 2 ; 	// Number of domains
    int nr = 17 ; 	// Number of collocation points in r in each domain
    int nt = 9 ; 	// Number of collocation points in theta in each domain
    int np = 8 ; 	// Number of collocation points in phi in each domain
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = NONSYM ; // no symmetry in phi
    bool compact = false ; // external domain is not compactified
  
    // Multi-domain grid construction:
    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
    cout << mgrid << endl ; 

  
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------

    // radial boundaries of each domain:
    double r_limits[] = {0., 2., 4.} ; 
    assert( nz == 2 ) ;  // since the above array described only 2 domains
  
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    cout << map << endl ;  
        
    // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = map.r ;        // r field 
    const Coord& x = map.x ;        // x field
    const Coord& y = map.y ;        // y field
    
    // Setup of a scalar field (initial value for the d'Alembert equation)
    // -------------------------------------------------------------------

    Scalar uu0(map) ;  // construction of an object of Lorene class Scalar
    
    uu0 = 2* exp( - r*r ) * (1 + x + x*y) ; 
        
    uu0.std_spectral_base() ; // sets the bases for the spectral expansions
                                 // to the standard ones for a scalar field

    cout << uu0 << endl ;    // prints to screen 
    
    uu0.spectral_display() ;     // prints the spectral expansions
    
    // 2-D visualization via PGPLOT
    // ----------------------------

    des_coupe_z( uu0, 0., 1, "Field U") ; 
    

    // 3-D visualization via OpenDX
    // ----------------------------
    
    // double z0 = 0 ;     // section plane : z = z0
  
    // uu0.visu_section('z', z0, -2., 2., -1.5, 1.5, "Example of section vis.") ;

    // uu0.visu_box(-2., 2., -1.5, 1.5, -1., 1., "Example of volume rendering", 0x0) ;
    
    // Time evolution : d'Alembert equation 
    // ------------------------------------
    
    Scalar source(map) ; // source of d'Alembert equation 
    source = 0 ; 

    double dt = 0.02 ;  // time step 
    int bc = 2 ;    // type of boundary condition : 2 = Bayliss & Turkel outgoing wave
    int *workflag = new int(0) ; // working flag 
 
    Param par ; 
    par.add_double(dt) ; 
    par.add_int(bc) ; 
    par.add_int_mod(*workflag) ; 
    
    
    double t = 0 ; 

    int j_min = 0 ;

    Evolution_std<Scalar> uu(uu0, 3, j_min, t) ; // Time evolution of U    
    // Time evolution of the central value of U : 
    Evolution_full<double> uu_c(uu0.val_grid_point(0,0,0,0), j_min, t) ; 
    j_min++ ;
    t += dt ;

    uu.update(uu0, j_min, t) ;
    uu_c.update(uu0.val_grid_point(0,0,0,0), j_min, t) ;

    int j_max = 500 ; 
    
    for (int j = j_min; j <= j_max ; j++) {
    
        Scalar uu_jp1 = uu[j].avance_dalembert(par, uu[j-1], source) ; 
    
	t += dt ;
        uu.update(uu_jp1, j+1, t) ; 
        
        uu_c.update(uu_jp1.val_grid_point(0,0,0,0), j+1, t) ; 
   
    
        if ( j%2 == 0 ) {

        
            cout << "Step = " << j+1 << ", time = " << t << endl ; 
        
            const Scalar* des[] = {&uu[j+1], &uu[j+1], &uu[j+1]} ;            
            double tab_theta0[] = {0.5*M_PI, 0.5*M_PI, 0.5*M_PI} ; 
            double tab_phi0[] = {0., 0.5*M_PI, M_PI} ; 
            
            des_profile_mult(des, 3, 0., 4., tab_theta0, tab_phi0, 1., false, 
            "U", 
            "U in the equatorial plane for phi=0, pi/2, pi", 0) ;
        
            double tab_theta1[] = {0, 0.25*M_PI, 0.5*M_PI} ; 
            double tab_phi1[] = {0., 0., 0.} ; 
 
            des_profile_mult(des, 3, 0., 4., tab_theta1, tab_phi1, 1., false, 
            "U", "U in the merional plane phi=0 for theta=0, pi/4, pi/2", 1) ;
        
        }
    }
    
    des_evol(uu_c, "U\\dc\\u", "Evolution of central value of U") ; 
    arrete() ; 
    
    
//  uu.visu_section('z', z0, -2., 2., -1.5, 1.5, "Potential", "uu") ;

    // Construction of a flat metric
    // -----------------------------

    Metric_flat mets(map, map.get_bvect_spher()) ; // spherical representation
    Metric_flat metc(map, map.get_bvect_cart()) ;  // Cartesian representation

    Vector duu = uu[j_max-1].derive_cov(metc) ; 
    
    // des_coupe_vect_z(duu, 0., -2., 0.5, 2, "Gradient of potential") ; 

    //duu.visu_arrows(-1., 1., -1., 1., -1., 1., "Gradient of potential", 
    //                 "gradient") ; 

    par.clean_all() ;

    return EXIT_SUCCESS ; 
}



