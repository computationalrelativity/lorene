/*
 *  Simple wave equation.
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon
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

// Lorene headers
#include "nbr_spx.h"
#include "tensor.h"
#include "metric.h"
#include "cmp.h"
#include "graphique.h"

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
        
    // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = map.r ;        // r field 
    const Coord& th = map.tet ;     // theta field 
    const Coord& phi = map.phi ;    // phi field 
    const Coord& x = map.x ;        // x field
    const Coord& y = map.y ;        // y field
    const Coord& z = map.z ;        // z field
    
    // Setup of a scalar field (source of the Poisson equation)
    // --------------------------------------------------------

    Scalar source(map) ;  // construction of an object of Lorene class Scalar
    
    source = exp( - r*r ) * (1 + x + x*y) ; 
    
    source.annule_domain(nz-1) ; // The source is set to zero in the last
                                 // domain 
    
    source.std_spectral_base() ; // sets the bases for the spectral expansions
                                 // to the standard ones for a scalar field
    //Scalar tmp(map) ; 
    //tmp = 2 ; 
    //source.set_domain(nz-1) = tmp.domain(nz-1) ; 
    //source.set_dzpuis(4) ; 
                        

    cout << source << endl ; 
    
    source.spectral_display() ;     // prints the spectral expansions
    
    // 2-D visualization via PGPLOT
    // ----------------------------

    des_coupe_z( Cmp(source), 0., 2, "Source") ; 
    

    // Output for 3-D visualization 
    // ----------------------------
    // Definition of the section plane :
    Tbl visuplane(3,3) ;        
    visuplane.set_etat_qcq() ;  // Memory allocation for the Tbl
     
    visuplane.set(0,0) = 0. ;   // Origin in the plane
    visuplane.set(0,1) = 0. ;   //  (absolute Cartesian coordinates)
    visuplane.set(0,2) = 0. ;   //
    
    visuplane.set(1,0) = 1. ;   // u-coordinate unit vector
    visuplane.set(1,1) = 0. ;   //  (absolute Cartesian components)
    visuplane.set(1,2) = 0. ;
    
    visuplane.set(2,0) = 0. ;   // v-coordinate unit vector
    visuplane.set(2,1) = 1. ;   //  (absolute Cartesian components)
    visuplane.set(2,2) = 0. ;
    
    source.visu_section(visuplane, -2., 2., -1., 2.) ;
    
    // Resolution of a Poisson equation 
    // --------------------------------
    
    Scalar pot = source.poisson() ; 
    
    cout << "Solution of the Poisson equation : " << endl ; 
    cout << pot << endl ; 
    
    pot.spectral_display() ;     // prints the spectral expansions 
                                     
    des_coupe_z( Cmp(pot), 0., 2, "Source") ; 

    // Construction of a flat metric
    // -----------------------------

    Metric_flat mets(map, map.get_bvect_spher()) ; // spherical representation
    Metric_flat metc(map, map.get_bvect_cart()) ;  // Cartesian representation

    Vector dpot = pot.derive_cov(metc) ; 

    return EXIT_SUCCESS ; 
}
