/*
 *  Code for reading alternative BH spacetime from a file
 */

/*
 *   Copyright (c) 2013 Odele Straub, Eric Gourgoulhon
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

char altBH_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2013/04/16 15:29:33  e_gourgoulhon
 * New code for reading Enrico's data
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
    
    ifstream fpar("par_altbh.d") ;
    if ( !fpar.good() ) {
        cerr << "Problem in opening the file par_altbh.d ! " << endl ;
        abort() ;
    }
    
    char file_name[256] ; 
    fpar.getline(file_name, 256) ;
    cout << "File to be read: " << file_name << endl ; 
    double a_spin ; 
    fpar >> a_spin ; fpar.ignore(1000,'\n') ;
    int graphic_out ; // flag for graphical outputs
    fpar >> graphic_out ; fpar.ignore(1000,'\n') ; 
    
    fpar.close();

    int nz = 4 ; 
    int nzm1 = nz-1 ; 
    double* r_limits = new double[nz+1];  // radial boundaries of each domain in units of M      
    r_limits[0] = 0 ; 
    for (int l=1; l<nz; l++) {
       r_limits[l] = pow(2,l); 
    }
    r_limits[nz] = __infinity ;

    cout << "r_limits : " ; 
    for (int l=0; l<nz+1; l++) {
      cout << r_limits[l] << "  " ; 
    }
    cout << endl ; 
    
    int nr = 33 ; 
    int nt = 5 ; 
    int np = 1 ; 
    
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

    // Construction of the AltBH_QI object:
    // ----------------------------------
    
    AltBH_QI bh(map, file_name, a_spin) ; 
    
    bh.update_metric() ; 

    cout.precision(15) ; 
    cout << bh << endl ; 

    // ISCO from Eq. (21) of Bardeen, Press & Teukolsky, ApJ 178, 347 (1972):
    cout << "Numerical value of r_ISCO (prograde orbits): " << 
        bh.r_isco(0) << endl ;  
        
    // Drawings    
    if (graphic_out == 1) {
        double r_max = 2*map.val_r(nzm1,-1.,0.,0.) ; 
        des_meridian(bh.get_nn(), 0, r_max, "N", 1) ; 

        des_meridian(bh.get_nphi(), 0, r_max, "Nphi", 3) ; 
        des_meridian(bh.get_nphi().dsdr(), 0, r_max, "dNphi/dr", 4) ; 

        des_meridian(bh.get_gamma().cov()(1,1), 0, r_max, "gamma_11", 5) ; 
        des_meridian(bh.get_gamma().cov()(3,3), 0, r_max, "gamma_33", 6) ; 
    
        des_meridian(bh.get_kk()(1,3), 0, r_max, "K_(r)(ph)", 7) ; 
        des_meridian(bh.get_kk()(2,3), 0, r_max, "K_(th)(ph)", 8) ; 
    
        arrete() ; 
    }

    
    // Output file for GYOTO
    //----------------------
    
    bh.gyoto_data("gyoto_altBH_QI.d") ;
            
    return EXIT_SUCCESS ; 

}













