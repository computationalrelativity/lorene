/*
 * Reads binary neutron star initial data, exporting Lorene structures
 * onto standard C arrays (double[]) on a Cartesian grid.
 */

/*
 *   Copyright (c) 2002  Eric Gourgoulhon
 *   Copyright (c) 2002  Keisuke Taniguchi
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

char read_bin_ns_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/06 15:09:47  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2003/01/09 11:08:00  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.2  2002/01/15 09:35:33  e_gourgoulhon
 * Added README file
 * Suppressed outputs/inputs in read_bin_ns.C
 *
 * Revision 1.1  2002/01/11 17:16:14  e_gourgoulhon
 * Main code and documentation for exporting binary neutron stars
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>

// Definition of Bin_NS class
#include "bin_ns.h"

int main() {

    // Reads Cartesian grid parameters
    // -------------------------------
    ifstream fparam("read_bin_ns.par") ;
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

    Bin_NS binary(nbp, xi, yi, zi, datafile) ;

    cout << endl << binary << endl ;


     // Save in a binary file
    // ---------------------
    FILE* finib = fopen("inib.d", "w") ;
    binary.save_bin(finib) ;
    fclose( finib ) ;

    // Save in a formatted file
    // ------------------------
    ofstream finif("inif.d") ;
    binary.save_form(finif) ;
    finif.close() ;

    // Create a new object from a binary file
    //---------------------------------------
    // finib = fopen("inib.d", "r") ;
    // Bin_NS binary2(finib) ;
    // fclose( finib ) ;

    // cout << endl << "Binary read in the binary file : " << binary2 << endl ;

    // finib = fopen("inib2.d", "w") ;
    // binary2.save_bin(finib) ;
    // fclose( finib ) ;

    // Create a new object from a formatted file
    //---------------------------------------
    // ifstream finif3("inif.d") ;
    // Bin_NS binary3(finif3) ;
    // finif3.close() ;

    // cout << endl << "Binary read in the formatted file : " << binary3 << endl ;

    // finif.open("inif3.d") ;
    // binary3.save_form(finif) ;
    // finif.close() ;

    // Clean exit
    // ----------

    delete [] xi ;
    delete [] yi ;
    delete [] zi ;

    return EXIT_SUCCESS ;
}
