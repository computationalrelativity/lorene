/*
 *  Code for reading Higgs monopole spacetime from a file
 */

/*
 *   Copyright (c) 2014 Marie Leroy,  Eric Gourgoulhon
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

char higgs_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/01/31 15:36:06  e_gourgoulhon
 * Added drawings
 *
 * Revision 1.1  2014/01/29 16:34:36  e_gourgoulhon
 * New main code for reading Higgs monopoles
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
#include "graphique.h"

int main() {

    // Parameters of the computation
    // -----------------------------
    
    ifstream fpar("par_higgs.d") ;
    if ( !fpar.good() ) {
        cerr << "Problem in opening the file par_higgs.d ! " << endl ;
        abort() ;
    }
    
    char file_name[256] ; 
    fpar.getline(file_name, 256) ;
    cout << "File to be read: " << file_name << endl ; 

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
    double* r_limits = new double[nz+1];  // inner boundaries of each domain in units of M      
    for (int l=0; l<nz+1; l++) {
        fpar >> r_limits[l]; 
    }
    //r_limits[nz] = __infinity ;
    
    fpar.close();
    
    cout << "r_limits : " ; 
    for (int l=0; l<nz+1; l++) {
        cout << r_limits[l] << "  " ; 
        }
    cout << endl ; 
    //arrete() ; 

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = SYM ; // symmetry with respect to phi --> phi + pi
    bool compact = false ; // external domain is compactified

    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;

    cout << mgrid << endl ; 
  
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------
  
    Map_af map(mgrid, r_limits) ;

    Mtbl r_grid(map.r) ; 
    ofstream grid_file("grid_file.d") ; 
    grid_file.precision(16) ; 
    for (int l = 0; l<nz; l++) {
        for (int i = 0; i<nr; i++) {
            grid_file << r_grid(l, 0, 0, i) << endl ; 
        }
    }
    grid_file.close() ; 
    

    // Construction of the AltBH_QI object:
    // ----------------------------------
    
    HiggsMonopole hmonop(map, file_name) ; 
    
    cout << hmonop << endl ;

// Drawings    
    if (graphic_out == 1) {
        double r_max = map.val_r(nzm1, 1.,0.,0.) ; 

        des_meridian(hmonop.get_nn(), 0, r_max, "N", 1) ; 
        des_meridian(hmonop.get_grr() , 0, r_max, "g_rr", 2) ; 
        des_meridian(hmonop.get_higgs() , 0, r_max, "h", 3) ; 
       des_meridian(hmonop.get_press() , 0, r_max, "P", 4) ; 
   
    arrete() ; 
    }



    return EXIT_SUCCESS ; 

}
