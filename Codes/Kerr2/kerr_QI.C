/*
 *  Kerr metric in Quasi-Isotropic coordinates
 */

/*
 *   Copyright (c) 2013 Claire Some, Eric Gourgoulhon
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

char kerrQI_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2013/04/03 12:11:19  e_gourgoulhon
 * Added member kk to Compobj; suppressed tkij
 *
 * Revision 1.1  2013/04/02 23:18:30  e_gourgoulhon
 * New code kerr_QI
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>
#include <cmath>

// Lorene headers
#include "compobj.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "proto.h"
#include "graphique.h"


int main() {

    // Parameters of the computation
    // -----------------------------
    
    ifstream fpar("par_kerr_QI.d") ;
    if ( !fpar.good() ) {
        cerr << "Problem in opening the file par_kerr_QI.d ! " << endl ;
        abort() ;
    }
    
    double mass ; // M
    fpar >> mass ; fpar.ignore(1000,'\n') ;

    double a_ov_m ; // Kerr parameter a/M
    fpar >> a_ov_m ; fpar.ignore(1000,'\n') ;
    
    int graphic_out ; // flag for graphical outputs
    fpar >> graphic_out ; fpar.ignore(1000,'\n') ;
   
    int nr ; // Number of collocation points in r in each domain
    fpar >> nr; fpar.ignore(1000,'\n') ;

    int nt ; // Number of collocation points in theta in each domain
    fpar >> nt; fpar.ignore(1000,'\n') ;

    int np ; // Number of collocation points in phi in each domain
    fpar >> np; fpar.ignore(1000,'\n') ;

    int nz ; // Number of domains
    fpar >> nz ; fpar.ignore(1000,'\n') ;
    int nzm1 = nz - 1 ; // Index of outermost domain

    fpar.ignore(1000,'\n') ; // skip title
    double* r_limits = new double[nz+1];  // radial boundaries of each domain in units of M      
    for (int l=0; l<nz; l++) {
        fpar >> r_limits[l]; 
    }
    r_limits[nz] = __infinity ;
    
    fpar.close();

    cout << "M = " << mass << ",  a/M = " << a_ov_m << endl ; 
    double r_hor = double(0.5)*mass*sqrt(double(1)-a_ov_m*a_ov_m) ; 
    cout << "Value of coordinate r at the event horizon : " << r_hor  << endl ; 
    cout << "r_limits : " ; 
    for (int l=0; l<nz+1; l++) {
      cout << r_limits[l] << "  " ; 
    }
    cout << endl ; 
    
    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = SYM ; // symmetry with respect to phi --> phi + pi
    bool compact = true ; // external domain is compactified

    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;

    cout << mgrid << endl ; 

  
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------
  
    Map_af map(mgrid, r_limits) ;

    // Construction of the Kerr_QI object:
    // ----------------------------------
    
    Kerr_QI bh(map, mass, a_ov_m) ; 
    
    bh.update_metric() ; 

    cout << bh << endl ; 
    
    // Drawings    
    if (graphic_out == 1) {
        des_meridian(bh.get_nn(), 0, 1.5*r_limits[nzm1], "N", 1) ; 

        des_meridian(bh.get_beta()(3), 0, 1.5*r_limits[nzm1], "beta\\u3\\d (ortho. comp.)", 3) ; 

        des_meridian(bh.get_nphi(), 0, 1.5*r_limits[nzm1], "Nphi", 3) ; 

        des_meridian(bh.get_gamma().cov()(1,1), 0, 1.5*r_limits[nzm1], "gamma_11", 4) ; 
        des_meridian(bh.get_gamma().cov()(2,2), 0, 1.5*r_limits[nzm1], "gamma_22", 5) ; 
        des_meridian(bh.get_gamma().cov()(3,3), 0, 1.5*r_limits[nzm1], "gamma_33", 6) ; 
    
        des_meridian(bh.get_kk()(1,3), 0, 1.5*r_limits[nzm1], "K_(r)(ph)", 7) ; 
        des_meridian(bh.get_kk()(2,3), 0, 1.5*r_limits[nzm1], "K_(th)(ph)", 8) ; 
    
        arrete() ; 
    }

    
    // Output file for GYOTO
    //----------------------
    
    bh.gyoto_data("gyoto_kerr_QI.d") ;
            
    return EXIT_SUCCESS ; 

}













