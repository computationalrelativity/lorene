/*
 *  Main code for time evolution of a wave packet within maximal slicing
 *  + Dirac gauge.
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

char wave_evol_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/04/07 07:59:22  e_gourgoulhon
 * Added check of constraints at the end.
 *
 * Revision 1.1  2004/04/05 21:26:25  e_gourgoulhon
 * First version (not ready yet !).
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Lorene headers
#include "time_slice.h"
#include "metric.h"
#include "evolution.h"
#include "param.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"

int main() {

    //======================================================================
    //      Parameters of the computation
    //======================================================================

    double pdt = 0.01 ; 
    int jmax = 1000 ; 
    int jstop = jmax ; 

    double ampli_h_init = 1. ;     // 0 = flat space
        

    //======================================================================
    //      Construction and initialization of the various objects
    //======================================================================

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 4 ; 	// Number of domains
    int nr = 9 ; 	// Number of collocation points in r in each domain
    int nt = 5 ; 	// Number of collocation points in theta in each domain
    int np = 8 ; 	// Number of collocation points in phi in each domain
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = SYM ; // no symmetry in phi
    bool compact = true ; // external domain is compactified
  
    // Multi-domain grid construction:
    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
    cout << mgrid << endl ; 

  
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------

    // radial boundaries of each domain:
    double r_limits[] = {0., 1., 2., 4., __infinity} ; 
    assert( nz == 4 ) ;  // since the above array described only 3 domains
  
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    cout << map << endl ;  
    
    // Flat metric f
    // -------------

    const Metric_flat& ff = map.flat_met_spher() ; 
    
    // Triad orthonormal with respect to the flat metric f
    // ----------------------------------------------------

    const Base_vect_spher& otriad = map.get_bvect_spher() ;
    
    
    // Construction of a time slice with maximal slicing and Dirac gauge
    // -----------------------------------------------------------------

    Tslice_dirac_max sigmat(map, otriad, ff) ;  

    // Set up of tensor h
    // ------------------
    
    Sym_tensor_trans hh_init(map, otriad, ff) ;  // hh is a transverse tensor
                                            // with respect to the flat metric
                                            // thanks to Dirac gauge
    const Coord& x = map.x ; 
    const Coord& y = map.y ; 
    // const Coord& z = map.z ; 
    const Coord& r = map.r ; 
    //const Coord& cost = map.cost ; 
    //const Coord& sint = map.sint ; 
    //const Coord& cosp = map.cosp ; 
    // const Coord& sinp = map.sinp ; 
    
    Scalar khi_init(map) ; 

    khi_init = ampli_h_init * exp( - r*r ) * x*y ;
    
    //khi_init = ampli_h_init * (3*cost*cost-1) / 
    //   ( (r*r + 1./(r*r)) * sqrt(1.+r*r) ) ; 
    khi_init.std_spectral_base() ; 

    khi_init.smooth_decay(2, 1) ; 
    
    //##
    // des_meridian(khi_init, 0., 3., "khi_init before", 1) ; 
    // arrete() ; 
    // khi_init.smooth_decay(3, 4) ; 
    // khi_init.spectral_display("khi_init") ;   
    //  des_meridian(khi_init, 0., 3., "khi_init after", 2) ; 
    //  arrete() ; 
    //## 
                 
    Scalar mu_init(map) ; 
    mu_init = 0. * ampli_h_init / (1+r*r*r*r*r*r) ; 
    mu_init.std_spectral_base() ; 
    mu_init.mult_r() ; 
    mu_init.mult_r() ; 
    mu_init.mult_r() ; 
    mu_init.mult_cost() ; 
    
    //##
    // des_meridian(mu_init, 0., 3., "mu_init before", 1) ; 
    // arrete() ; 
    // mu_init.smooth_decay(3, 4) ; 
    // mu_init.spectral_display("mu_init") ;   
    // des_meridian(mu_init, 0., 3., "mu_init after", 2) ; 
    // arrete() ; 
    //## 
                  

    Sym_tensor_tt htt_init(map, otriad, ff) ;  // htt is the TT part of hh
        
    htt_init.set_khi_mu(khi_init, mu_init) ; 
    
    hh_init = htt_init ; 
        
    // des_meridian(hh_init, 0., 5., "hh_init") ; 
    maxabs( hh_init.divergence(ff), "Divergence of hh_init") ; 
    maxabs( hh_init.trace(ff), "Trace of hh_init") ; 

    // arrete() ; 

    // Resolution of the initial data equations within 
    // the conformal thin sandwich framework
    // ------------------------------------------------

    // u^{ij} = d/dt h^{ij}
    Sym_tensor_trans uu_init(map, otriad, ff) ;  
    //uu_init.set_etat_zero() ; 
    uu_init = - 0.5 * hh_init ;
    uu_init.inc_dzpuis(2) ;  
    
    // tr K = K
    Scalar tmp(map) ; 
    tmp.set_etat_zero() ; 
    
    sigmat.initial_data_cts(hh_init, uu_init, tmp, tmp) ;
    
    cout << "sigmat : " << sigmat << endl ;  
    
    // Check of constraints:
    sigmat.check_hamiltonian_constraint() ;    
    sigmat.check_momentum_constraint() ; 

    return EXIT_SUCCESS ; 
}

