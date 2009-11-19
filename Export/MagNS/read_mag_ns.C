/*
 * Reads magnetized neutron star initial data, exporting Lorene structures
 * onto standard C arrays (double[]) on a Cartesian grid.
 */

/*
 *   Copyright (c) 2002  Eric Gourgoulhon
 *   Copyright (c) 2009 Jerome Novak
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

char read_mag_ns_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2009/11/19 16:15:21  j_novak
 * Export class for magnetized neutron stars.
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdlib.h>

// Definition of Mag_NS class
#include "mag_ns.h"

int main() {

    // Reads Cartesian grid parameters
    // -------------------------------
    ifstream fparam("read_mag_ns.par") ;
    char comment[120] ;
    char datafile[120] ;
    int nx, ny, nz ;
    double x_min, x_max, y_min, y_max, z_min, z_max ;
    fparam.getline(comment, 120) ;
    fparam.getline(comment, 120) ;
    fparam.getline(datafile, 120) ;
    fparam.getline(comment, 120) ;
    fparam.getline(comment, 120) ;
    fparam.getline(comment, 120) ;
    fparam >> nx ; fparam.getline(comment, 120) ;
    fparam >> ny ; fparam.getline(comment, 120) ; 
    fparam >> nz ; fparam.getline(comment, 120) ; 
    fparam >> x_min ; fparam >> x_max ; fparam.getline(comment, 120) ; 
    fparam >> y_min ; fparam >> y_max ; fparam.getline(comment, 120) ; 
    fparam >> z_min ; fparam >> z_max ; fparam.getline(comment, 120) ; 
    fparam.close() ; 
    
    cout << "File containing the initial data on the spectral grid:"
	 << endl ; 
    cout << datafile << endl << endl ;
    cout << "Cartesian grid : " << endl ; 
    cout << "---------------- " << endl ; 
    cout << "   Number of points in the x direction : " << nx << endl ; 
    cout << "   Number of points in the y direction : " << ny << endl ; 
    cout << "   Number of points in the z direction : " << nz << endl ; 
    cout << "   x_min, x_max : " << x_min << " , " << x_max << endl ; 
    cout << "   y_min, y_max : " << y_min << " , " << y_max << endl ; 
    cout << "   z_min, z_max : " << z_min << " , " << z_max << endl ; 
    
    // Construction of the Cartesian grid
    // ----------------------------------

    int nbp = nx*ny*nz ; 
    double* const xi = new double[nbp] ; 
    double* const yi = new double[nbp] ; 
    double* const zi = new double[nbp] ; 
    
    double dx = (x_max - x_min) / double(nx - 1) ;     
    double dy = (y_max - y_min) / double(ny - 1) ;     
    double dz = (z_max - z_min) / double(nz - 1) ;   
    
    double* pxi = xi ;
    double* pyi = yi ;
    double* pzi = zi ;

    for (int k=0; k<nz; k++) {

	double z = z_min + dz * k ;

	for (int j=0; j<ny; j++) {

	    double y = y_min + dy * j ;

	    for (int i=0; i<nx; i++) {

		*pxi = x_min + dx * i ;
		*pyi = y ;
		*pzi = z ;

		pxi++ ;
		pyi++ ;
		pzi++ ;
	    }
	}
    }

    // Read of the initial data file and computation of the
    //  fields at the Cartesian grid points
    // ----------------------------------------------------

    Mag_NS magstar(nbp, xi, yi, zi, datafile) ;

    cout << endl << magstar << endl ;


     // Save in a binary file
    // ---------------------
    FILE* finib = fopen("inib.d", "w") ;
    magstar.save_bin(finib) ;
    fclose( finib ) ;

    // Save in a formatted file
    // ------------------------
    ofstream finif("inif.d") ;
    magstar.save_form(finif) ;
    finif.close() ;

    // Create a new object from a binary file
    //---------------------------------------
    // finib = fopen("inib.d", "r") ;
    // Mag_NS magstar2(finib) ;
    // fclose( finib ) ;

    // cout << endl << "Magstar read in the binary file : " << magstar2 << endl ;

    // finib = fopen("inib2.d", "w") ;
    // magstar2.save_bin(finib) ;
    // fclose( finib ) ;

    // Create a new object from a formatted file
    //---------------------------------------
    // ifstream finif3("inif.d") ;
    // Mag_NS magstar3(finif3) ;
    // finif3.close() ;

    // cout << endl << "Binary read in the formatted file : " << magstar3 << endl ;

    // finif.open("inif3.d") ;
    // magstar3.save_form(finif) ;
    // finif.close() ;

    // Clean exit
    // ----------

    delete [] xi ;
    delete [] yi ;
    delete [] zi ;

    return EXIT_SUCCESS ;
}
