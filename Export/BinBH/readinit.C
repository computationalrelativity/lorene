/*
 * Reads binary black hole initial data, exporting Lorene structures
 * onto standard C arrays (double[]) on a Cartesian grid.
 */

/*
 *   Copyright (c) 2001-2002  Eric Gourgoulhon
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
 * Revision 1.9  2016/12/05 16:18:30  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.8  2014/10/13 08:54:04  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.7  2014/10/06 15:09:46  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2003/01/09 11:07:59  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.5  2002/03/20 08:24:56  e_gourgoulhon
 * Added the derivatives of Psi.
 *
 * Revision 1.4  2002/02/06 14:54:44  e_gourgoulhon
 * Update of bibliographical references
 *
 * Revision 1.3  2002/01/09 16:32:57  e_gourgoulhon
 * Added the parameter fill in readinit.par
 * Suppressed the multiple outputs in files ini* via method sauve
 * Output directly written in readinit.C
 *
 * Revision 1.2  2001/12/18 09:08:14  e_gourgoulhon
 * Adds the filling of the holes interiors
 *
 * Revision 1.1  2001/12/14 08:59:18  e_gourgoulhon
 * Exportation of Lorene Bhole_binaire object to a Cartesian grid
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>

// Definition of Bin_BH class
#include "bin_bh.h"

using namespace Lorene ;

int main() {

    // Reads Cartesian grid parameters
    // -------------------------------
    ifstream fparam("readinit.par") ; 
    char comment[120] ; 
    char datafile[120] ;
    int nx, ny, nz, fill ;
    double x_min, x_max, y_min, y_max, z_min, z_max ;
    fparam.getline(comment, 120) ;
    fparam.getline(comment, 120) ;
    fparam.getline(datafile, 120) ;
    fparam.getline(comment, 120) ;
    fparam.getline(comment, 120) ; 
    fparam.getline(comment, 120) ;
    fparam >> fill ; fparam.getline(comment, 120) ;
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
    
    double dx = (nx == 1) ? 0 : (x_max - x_min) / double(nx - 1) ;
    double dy = (ny == 1) ? 0 : (y_max - y_min) / double(ny - 1) ;
    double dz = (nz == 1) ? 0 : (z_max - z_min) / double(nz - 1) ;

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

    Bin_BH binary(nbp, xi, yi, zi, fill, datafile) ;

    cout << binary << endl ;

    ofstream fich("readout.d") ;

    fich.precision(13) ;
    fich.setf(ios::scientific, ios::floatfield) ;

    fich << binary.dist << " dist [a]" << endl ;
    fich << "nx, ny, nz : " << nx << "  " << ny << "  " << nz << endl ;
    fich << "x_min, dx : " << x_min << "  " << dx << endl ;
    fich << "y_min, dy : " << y_min << "  " << dy << endl ;
    fich << "z_min, dz : " << z_min << "  " << dz << endl ;

    int index = 0 ;

    for (int k=0; k<nz; k++) {

	for (int j=0; j<ny; j++) {

	    for (int i=0; i<nx; i++) {

                fich << i
                  << " " <<  binary.g_xx[index] << " " <<  binary.g_xy[index] << " " <<  binary.g_xz[index]
                  << " " <<  binary.g_yy[index] << " " <<  binary.g_yz[index] << " " <<  binary.g_zz[index]
                  << " " <<  binary.k_xx[index] << " " <<  binary.k_xy[index] << " " <<  binary.k_xz[index]
                  << " " <<  binary.k_yy[index] << " " <<  binary.k_yz[index] << " " <<  binary.k_zz[index]
                  << " " <<  binary.beta_x[index] << " " <<  binary.beta_y[index]
                  << " " <<  binary.beta_z[index]  << " " <<  binary.nnn[index] << endl ;

                index++ ;

                        }
                }

        }


    fich.close() ;



    // Save in a binary file
    // ---------------------
    FILE* finib = fopen("binbh.d", "w") ;
    binary.save_bin(finib) ;
    fclose( finib ) ;

    // Save in a formatted file
    // ------------------------
    // ofstream finif("inif.d") ;
    // binary.save_form(finif) ;
    // finif.close() ;


    // Clean exit
    // ----------

    delete [] xi ;
    delete [] yi ;
    delete [] zi ;

    return EXIT_SUCCESS ;

}
