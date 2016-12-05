/*
 *  Tests of class Time_slice
 *
 *    (see file time_slice.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Eric Gourgoulhon, Jose Luis Jaramillo & Jerome Novak
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
 * Revision 1.8  2016/12/05 16:18:27  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2014/10/13 08:54:00  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:12:52  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2004/04/05 11:54:43  j_novak
 * First operational (but not tested!) version of checks of Eintein equation.
 *
 * Revision 1.4  2004/03/30 14:02:05  j_novak
 * Test of the class Tslide_dirac_max (preliminary version).
 *
 * Revision 1.3  2004/03/29 12:01:12  e_gourgoulhon
 * Added test for new class Time_slice_conf.
 *
 * Revision 1.2  2004/03/26 13:33:03  j_novak
 * New methods for accessing/updating members (nn(), beta(), gam_uu(), k_uu(), ...)
 *
 * Revision 1.1  2004/03/24 14:58:19  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>

// Lorene headers
#include "nbr_spx.h"
#include "tensor.h"
#include "metric.h"
#include "time_slice.h"
#include "utilitaires.h"
#include "graphique.h"

using namespace Lorene ;

int main() {

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
    
    // Construction of a time slice Sigma_t
    // ------------------------------------
    
    Time_slice sigma(map, map.get_bvect_spher()) ; 

    sigma.set_scheme_order(0) ; //stationary space-time

//     cout << sigma << endl ; 
//     arrete() ; 

//     cout << sigma.nn() ;
//     arrete() ;
//     cout << sigma.beta() ;
//     arrete() ;
//     cout << sigma.gam_uu() ;
//     arrete() ;
//     cout << sigma.gam_dd() ;
//     arrete() ;
//     cout << sigma.k_dd() ;
//     arrete() ;
//     cout << sigma.k_uu() << endl ;
    
//     cout << "nn : " << sigma.nn() << endl ; 
//     cout << "beta : " << sigma.beta() << endl ; 
//     cout << "gam_uu : " << sigma.gam_uu() << endl ; 
//     cout << "gam : " << sigma.gam() << endl ; 
//     cout << "gam_dd : " << sigma.gam_dd() << endl ; 
    
    
    // Construction of a time slice with conformal decomposition
    // ---------------------------------------------------------
    
    Time_slice_conf sigma_c(sigma.nn(), sigma.beta(), sigma.gam_uu(),
                sigma.k_uu(), map.flat_met_spher()) ; 
    
    //cout << sigma_c << endl ; 

    const Coord& x = map.x ; 
    const Coord& y = map.y ; 
    const Coord& r = map.r ; 

    Scalar khi_init(map) ; 
    khi_init = exp( - r*r ) * x*y ;
    khi_init.std_spectral_base() ; 
    khi_init.set_outer_boundary(nz-1, 0.) ;
    
    Scalar mu_init(map) ; 
    mu_init = 1. / (1.+r*r*r*r*r*r) ; 
    mu_init.std_spectral_base() ; 
    mu_init.mult_r() ; 
    mu_init.mult_r() ; 
    mu_init.mult_r() ; 
    mu_init.mult_cost() ; 
    
    Sym_tensor_tt hijtt(map, map.get_bvect_spher() , map.flat_met_spher()) ;
    hijtt.set_khi_mu( khi_init, mu_init ) ;

    Scalar trace_nulle(map) ;
    trace_nulle = 0 ;
    Sym_tensor_trans hij(map, map.get_bvect_spher() , map.flat_met_spher()) ;
    hij.set_tt_trace(hijtt, trace_nulle) ;
    
    Tslice_dirac_max dirac_slice(sigma_c.nn(), sigma_c.beta() , 
				 map.flat_met_spher(), sigma_c.psi(), 
				 hij, sigma_c.aa()) ;

    dirac_slice.set_scheme_order(0) ; //stationary space-time
    cout << dirac_slice ;

    dirac_slice.check_hamiltonian_constraint() ;

    dirac_slice.check_momentum_constraint() ;

    dirac_slice.check_dynamical_equations() ;
    

    return EXIT_SUCCESS ; 
}
