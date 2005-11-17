/*
 *  Test of arithmetic Scalar & Mtbl and Scalar & Coord
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
 * Revision 1.2  2005/11/17 16:32:35  e_gourgoulhon
 * Added test of ETATZERO and ETATUN cases.
 *
 * Revision 1.1  2005/11/17 15:31:53  e_gourgoulhon
 * Code for testing Scalar / Mtbl (or Coord) new arithmetics.
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// Lorene headers
#include "headcpp.h"    // standard input/output C++ headers (iostream, fstream)
#include "tensor.h"     // classes Tensor, etc...     
#include "nbr_spx.h"    // defines infinity as an ordinary number: __infinity
#include "graphique.h"  // for graphical outputs
#include "utilitaires.h"    // utilities 

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
    const Coord& z = map.z ;      
    
    // Some scalar field to be used as a conformal factor
    // --------------------------------------------------
    
    Scalar aa(map) ; 
    
    aa = 1 + y*y ;
    
    aa.std_spectral_base() ;  // Standar polynomial bases will be used 
                                // to perform the spectral expansions

    cout << endl << "aa : " << endl ; 
    aa.spectral_display() ; 
    
    Scalar bb = aa * x   ; 
    
    bb.std_spectral_base() ;
    //cout << endl << "bb : " << endl ; 
    //bb.spectral_display() ; 
    
    Scalar bb0(map) ; 

    bb0 = (1+y*y) * x ; 
    cout << "test Scalar * Coord : " << max(abs(bb-bb0)) << endl ;   
    
    bb = x * aa ; 
    bb0 = x * (1+y*y) ; 
    cout << "test Coord * Scalar : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa / (1+x*x) ; 
    bb0 = (1+y*y) / (1+x*x) ; 
    cout << "test Scalar / Mtbl : " << max(abs(bb-bb0)) << endl ;   
    
    bb = x / aa ; 
    bb0 = x / (1+y*y) ; 
    cout << "test Coord / Scalar : " << max(abs(bb-bb0)) << endl ;   
    
    bb = x + aa ; 
    bb0 = x + (1+y*y) ; 
    cout << "test Coord + Scalar : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa + x ; 
    bb0 = x + (1+y*y) ; 
    cout << "test Scalar + Coord : " << max(abs(bb-bb0)) << endl ;   
    
    bb = x - aa ; 
    bb0 = x - (1+y*y) ; 
    cout << "test Coord - Scalar : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa - x ; 
    bb0 = (1+y*y) - x ; 
    cout << "test Scalar - Coord : " << max(abs(bb-bb0)) << endl ;   
    
    
    arrete() ; 
    
    // Special case ETATUN

    cout << endl << "Special case ETATUN : " << endl ; 
    aa.set_etat_one() ; 
    bb = aa * x ; 
    bb0 = x ; 
    cout << "test Scalar(1) * Coord : " << max(abs(bb-bb0)) << endl ;   
    
    bb = x * aa ; 
    bb0 = x ; 
    cout << "test Coord * Scalar(1) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa / (1+x*x) ; 
    bb0 = 1 / (1+x*x) ; 
    cout << "test Scalar(1) / Mtbl : " << max(abs(bb-bb0)) << endl ;   
    
    bb = x / aa ; 
    bb0 = x  ; 
    cout << "test Coord / Scalar(1) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = x + aa ; 
    bb0 = x + 1 ; 
    cout << "test Coord + Scalar(1) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa + x ; 
    bb0 = x + 1 ; 
    cout << "test Scalar(1) + Coord : " << max(abs(bb-bb0)) << endl ;   
    
    bb = x - aa ; 
    bb0 = x - 1 ; 
    cout << "test Coord - Scalar(1) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa - x ; 
    bb0 = 1 - x ; 
    cout << "test Scalar(1) - Coord : " << max(abs(bb-bb0)) << endl ;   
    

    // Special case ETATZERO

    cout << endl << "Special case ETATZERO : " << endl ; 
    aa.set_etat_zero() ; 
    bb = aa * x ; 
    bb0 = 0 ; 
    cout << "test Scalar(0) * Coord : " << max(abs(bb-bb0)) << endl ;   
    
    bb = x * aa ; 
    bb0 = 0 ; 
    cout << "test Coord * Scalar(0) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa / (1+x*x) ; 
    bb0 = 0 ; 
    cout << "test Scalar(0) / Mtbl : " << max(abs(bb-bb0)) << endl ;   
    
    bb = x + aa ; 
    bb0 = x  ; 
    cout << "test Coord + Scalar(0) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa + x ; 
    bb0 = x  ; 
    cout << "test Scalar(0) + Coord : " << max(abs(bb-bb0)) << endl ;   
    
    bb = x - aa ; 
    bb0 = x  ; 
    cout << "test Coord - Scalar(0) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa - x ; 
    bb0 =  - x ; 
    cout << "test Scalar(0) - Coord : " << max(abs(bb-bb0)) << endl ;   
    
    arrete() ; 
    
    Mtbl mc(mgrid) ; 
    mc.set_etat_zero() ; 

    aa = 1 + y*y ;

    bb = aa * mc ; 
    bb0 = 0 ; 
    cout << "test Scalar * Mtbl(0) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = mc * aa ; 
    bb0 = 0 ; 
    cout << "test Mtbl(0) * Scalar : " << max(abs(bb-bb0)) << endl ;   
    
    bb = mc / aa ; 
    bb0 = 0 ; 
    cout << "test Mtbl(0) / Scalar : " << max(abs(bb-bb0)) << endl ;   
    
    bb = mc + aa ; 
    bb0 = (1+y*y) ; 
    cout << "test Mtbl(0) + Scalar : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa + mc ; 
    bb0 =  (1+y*y) ; 
    cout << "test Scalar + Mtbl(0) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = mc - aa ; 
    bb0 = - (1+y*y) ; 
    cout << "test Mtbl(0) - Scalar : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa - mc ; 
    bb0 = (1+y*y) ; 
    cout << "test Scalar - Mtbl(0) : " << max(abs(bb-bb0)) << endl ;   

    arrete() ;
    
    // Special case ETATUN

    cout << endl << "Special case Scalar = ETATUN  and Mtbl = ETATZERO: " << endl ; 
    aa.set_etat_one() ; 

    bb = aa * mc ; 
    bb0 = 0 ; 
    cout << "test Scalar(1) * Mtbl(0) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = mc * aa ; 
    bb0 = 0 ; 
    cout << "test Mtbl(0) * Scalar(1) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = mc / aa ; 
    bb0 = 0 ; 
    cout << "test Mtbl(0) / Scalar(1) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = mc + aa ; 
    bb0 = 1 ; 
    cout << "test Mtbl(0) + Scalar(1) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa + mc ; 
    bb0 =  1 ; 
    cout << "test Scalar(1) + Mtbl(0) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = mc - aa ; 
    bb0 = - 1 ; 
    cout << "test Mtbl(0) - Scalar(1) : " << max(abs(bb-bb0)) << endl ;   
    
    bb = aa - mc ; 
    bb0 = 1 ; 
    cout << "test Scalar(1) - Mtbl(0) : " << max(abs(bb-bb0)) << endl ;   



    return EXIT_SUCCESS ; 

}
