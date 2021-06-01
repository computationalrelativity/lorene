/*
 * Test program for the radial primitive computation with a Map_af mapping
 * 
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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
 * Revision 1.5  2021/06/01 14:21:48  j_novak
 * Updated the arguments of primr()
 *
 * Revision 1.4  2016/12/05 16:18:28  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:54:01  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:12:53  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/06/14 15:29:12  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

// headers C++
#include "headcpp.h"

// headers C
#include <cstdlib>
#include <cmath>

// headers Lorene
#include "tensor.h"
#include "nbr_spx.h"
#include "graphique.h"
#include "utilitaires.h"

//******************************************************************************

using namespace Lorene ;

int main(){
    
    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 4 ; 	// Number of domains
    int nzm1 = nz-1 ; 	
    int nr = 17 ; 	// Number of collocation points in r in each domain
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
    double r_limits[] = {0., 1.5, 2., 3., __infinity} ; 
    assert( nz == 4 ) ;  // since the above array described only 3 domains
  
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    cout << map << endl ;  
    
    // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = map.r ;        // r field 
    // const Coord& th = map.tet ;     // theta field 
    // const Coord& phi = map.phi ;    // phi field 
    const Coord& x = map.x ;        // x field
    const Coord& y = map.y ;        // y field
    const Coord& z = map.z ;        // z field
    
    
    //-----------------------------------------------------------------------
    //		Construction of a scalar field
    //-----------------------------------------------------------------------

    Scalar tmp(map) ; 
    tmp = r*r ;
    tmp.annule_domain(nzm1) ;  
    tmp.std_spectral_base() ; 
    
    double rced = map.val_r(nz-2, 1., 0., 0.) ; 
    double ray_des = 2. * rced ; // outermost radius for plots

    Scalar prim_ana(map) ;  
    
    prim_ana = r*( 1- x*x + z*z + x*y) ; 
    // prim_ana = r*r*r ; 
    Mtbl aux(mgrid) ; 
    aux = pow(rced,4) / r ; 
    prim_ana.set_domain(nzm1) = aux(nzm1) ; 
            
    prim_ana.set_spectral_base( tmp.dsdr().get_spectral_va().get_base() ) ; 

    des_meridian(prim_ana, 0., ray_des, "prim_ana", 9) ; 

    Scalar uu = prim_ana.dsdr() ; 

    uu.spectral_display("uu") ;     // prints the spectral expansions
    
    Scalar uu_des = uu ; 
    uu_des.dec_dzpuis(2) ; 
    des_meridian(uu_des, 0., ray_des, "U", 10) ; 
    
    Scalar prim = uu.primr(false) ; 
    prim.spectral_display("Primitive") ;     // prints the spectral expansions
    
    des_meridian(prim - prim_ana, 0., ray_des, "Primitive", 11) ; 
    
    maxabs(prim - prim_ana, "Absolute error prim - prim_ana") ; 
    
    maxabs(prim.dsdr() - uu, "Absolute error d/dr prim - uu") ; 
    
    arrete() ; 
    
    return EXIT_SUCCESS ; 
    
}

