/*
 *  Main code for test of template class Evolution
 *
 *    (see file evolution.h for documentation).
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

char test_evol_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2004/03/26 08:24:03  e_gourgoulhon
 * Takes into account the new setting of class Evolution.
 *
 * Revision 1.4  2004/03/06 21:13:52  e_gourgoulhon
 * Added test of time derivation.
 *
 * Revision 1.3  2004/02/16 10:30:04  e_gourgoulhon
 * Added #include <math.h>
 *
 * Revision 1.2  2004/02/15 22:01:36  e_gourgoulhon
 * New version to take into account the split of Evolution
 * into Evolution_full and Evolution_std.
 *
 * Revision 1.1  2004/02/13 15:54:03  e_gourgoulhon
 * First version.
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

// Lorene headers
#include "evolution.h"
#include "tensor.h"
#include "utilitaires.h"

int main() {


    //------------------------------------------------
    //      Test with a double
    //------------------------------------------------

    Evolution_full<double> aa(1.) ; 
    Evolution_std<double> bb(1., 3) ; 
    Evolution_std<double> bb1(1., 3) ; 
    
    cout << "aa[0] : " << aa[0] << endl ; 
    cout << "bb[0] : " << bb[0] << endl ; 
    
    // Test of time derivative: 
    
    double dt = 0.01 ; 
    double t = dt ; 
    double diffa = 0 ; 
    double diffb = 0 ; 
    double diffb1 = 0 ; 
    for (int j=1; j<400; j++) {
        
        cout << "j = " << j << "  t = " << t << " :" << endl ; 
        double f = cos( t ) ; 
        double df = - sin( t ) ; 
        aa.update(f, j, t) ; 
        bb.update(f, j, t) ; 
        bb1.update(f, j, t) ; 
        
        if (j>1) {
            double dfa = fabs( aa.time_derive(j) - df ); 
            double dfb = fabs( bb.time_derive(j) - df ); 
            double dfb1 = fabs( bb1.time_derive(j,1) - df ); 
            if (dfa>diffa) diffa = dfa ; 
            if (dfb>diffb) diffb = dfb ; 
            if (dfb1>diffb1) diffb1 = dfb1 ; 
            cout << "Error time derivative aa, bb, bb1 : " << dfa << " " << 
                dfb << " " << dfb1 << endl ; 
        }
        
        t += dt ; 
    }
    cout << endl << "Max error on time derivative aa : " << diffa << endl ; 
    cout << "Max error on time derivative bb : " << diffb << endl ; 
    cout << "Max error on time derivative bb1 : " << diffb1 << endl ; 
    arrete() ; 

    Evolution_std<double> cc(2., 3) ; 

    for (int j=1; j<20; j++) {
    
        double t_j = double(j) / double(10) ;
        cc.update(sqrt(double(j)), j, t_j) ; 
        
        cout << "time : " << cc.get_time(j) << " :  "  ; 
        if (j==1) cout << cc[0] << "  " << cc[1] << endl ; 
        if (j>=2) cout << cc[j-2] << "  " << cc[j-1] << "  " << cc[j] << endl ; 
        

    }
    cout << endl ; 

    arrete() ; 
    

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 2 ; 	// Number of domains
    int nr = 7 ; 	// Number of collocation points in r in each domain
    int nt = 5 ; 	// Number of collocation points in theta in each domain
    int np = 8 ; 	// Number of collocation points in phi in each domain
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = NONSYM ; // no symmetry in phi
    bool compact = false ; // external domain is not compactified
  
    // Multi-domain grid construction:
    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------

    // radial boundaries of each domain:
    double r_limits[] = {0., 0.5, 1.} ; 
    assert( nz == 2 ) ;  // since the above array described only 2 domains
  
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = map.r ;        // r field 
    const Coord& x = map.x ;        // x field
    const Coord& y = map.y ;        // y field
    
    // Setup of a scalar field (source of the Poisson equation)
    // --------------------------------------------------------

    Scalar source(map) ;  // construction of an object of Lorene class Scalar
    
    source = 2* exp( - r*r ) * (1 + x + x*y) ; 
        
    source.std_spectral_base() ; // sets the bases for the spectral expansions
                                 // to the standard ones for a scalar field

    
    Evolution_std<Scalar> evol(source, 3) ; 
    
    cout << evol[0] << endl ; 

    for (int j=0; j<20; j++) {
    
        double t_j = double(j) / double(10) ;
        Scalar tmp(map) ; 
        tmp =  sqrt(t_j) ; 
        tmp.std_spectral_base() ; 
        evol.update(tmp, j, t_j) ; 
        
        evol[j].spectral_display() ; 

    }
    cout << endl ; 


    return EXIT_SUCCESS ; 

}
