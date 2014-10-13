/*
 *  Reads from a file a binary black hole configuration on a Cartesian grid
 *  (created by readinit), extract one field and save it in a formatted
 *  file for producing plots with Grace.
 *
 */

/*
 *   Copyright (c) 2002  Eric Gourgoulhon
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

char visu_1d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:54:04  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:09:46  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/01/09 11:07:59  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1  2002/03/20 08:24:56  e_gourgoulhon
 * Added the derivatives of Psi.
 *
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>
#include <cmath>

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

    cout << "File containing the initial data on the Cartesian grid:"
	 << endl ;
    cout << "bin_bh.d" << endl << endl ;
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

    double* const xx = new double[nx] ;
    double* const yy = new double[ny] ;
    double* const zz = new double[nz] ;

    double dx = (nx == 1) ? 0 : (x_max - x_min) / double(nx - 1) ;
    double dy = (ny == 1) ? 0 : (y_max - y_min) / double(ny - 1) ;
    double dz = (nz == 1) ? 0 : (z_max - z_min) / double(nz - 1) ;

    for (int i=0; i<nx; i++) {
	xx[i] = x_min + dx * i ;
    }

    for (int i=0; i<ny; i++) {
	yy[i] = y_min + dy * i ;
    }

    for (int i=0; i<nz; i++) {
	zz[i] = z_min + dz * i ;
    }

    // Read configuration from file
    FILE* fichbin = fopen("binbh.d", "r") ;
    Bin_BH binary(fichbin) ;
    fclose(fichbin) ;

    cout << endl << "Configuration read in file: " << endl ;
    cout << binary << endl ;

    cout << "Type of profiles : " << endl ;
    cout << "----------------   " << endl ;
    cout << "  1 : along x, for all y, at fixed value of z" << endl ;
    cout << "  2 : along x, for all z, at fixed value of y" << endl ;
    cout << "  3 : along y, for all x, at fixed value of z" << endl ;
    cout << "  4 : along y, for all z, at fixed value of x" << endl ;
    cout << "  5 : along z, for all x, at fixed value of y" << endl ;
    cout << "  6 : along z, for all y, at fixed value of x" << endl ;
    cout << "Your choice ?" << endl ;
    int type_plot ;
    cin >> type_plot ;
    int indplane ;
    switch (type_plot) {
    	case 1 : case 3 : {
    		cout << "z index (defining the plot plane) ?" << endl ;
    		cin >> indplane ;
    		if ( (indplane<0) || (indplane>=nz) ) {
    			cout << "visu_1d: z index out of range !" << endl ;
    			abort() ;
    		}
    		break ;
	}

    	case 2 : case 5 : {
    		cout << "y index (defining the plot plane) ?" << endl ;
    		cin >> indplane ;
    		if ( (indplane<0) || (indplane>=ny) ) {
    			cout << "visu_1d: y index out of range !" << endl ;
    			abort() ;
    		}
    		break ;
	}

    	case 4 : case 6 : {
    		cout << "x index (defining the plot plane) ?" << endl ;
    		cin >> indplane ;
    		if ( (indplane<0) || (indplane>=nx) ) {
    			cout << "visu_1d: x index out of range !" << endl ;
    			abort() ;
    		}
    		break ;
	}
	
	default : {
		cout << "visu_1d: Invalid type of profiles !" << endl ;
		abort() ;
	}

    }


    ofstream fich("plot.d") ;

    fich.precision(8) ;
    fich.setf(ios::scientific, ios::floatfield) ;

    fich << "# " << binary.dist << " dist [a]" << endl ;
    fich << "# nx, ny, nz : " << nx << "  " << ny << "  " << nz << endl ;
    fich << "# x_min, dx : " << x_min << "  " << dx << endl ;
    fich << "# y_min, dy : " << y_min << "  " << dy << endl ;
    fich << "# z_min, dz : " << z_min << "  " << dz << endl ;

    const double* uu = binary.nnn ;

    switch (type_plot) {

    	case 1 : {	

    		for (int i=0; i<nx; i++) {
    			fich << xx[i] ;
    			for (int j=0; j<ny; j++) {
    				int index = indplane*nx*ny + j*nx + i ;
    				fich << "  " << uu[index] ;
    			}
    			fich << endl ;
    		}
    		break ;
    	}

    	case 2 : {	

    		for (int i=0; i<nx; i++) {
    			fich << xx[i] ;
    			for (int j=0; j<nz; j++) {
    				int index = j*nx*ny + indplane*nx + i ;
    				fich << "  " << uu[index] ;
    			}
    			fich << endl ;
    		}
    		break ;
    	}

    	case 3 : {	

    		for (int i=0; i<ny; i++) {
    			fich << yy[i] ;
    			for (int j=0; j<nx; j++) {
    				int index = indplane*nx*ny + i*nx + j ;
    				fich << "  " << uu[index] ;
    			}
    			fich << endl ;
    		}
    		break ;
    	}

    	case 4 : {	

    		for (int i=0; i<ny; i++) {
    			fich << yy[i] ;
    			for (int j=0; j<nz; j++) {
    				int index = j*nx*ny + i*nx + indplane ;
    				fich << "  " << uu[index] ;
    			}
    			fich << endl ;
    		}
    		break ;
    	}

    	case 5 : {	

    		for (int i=0; i<nz; i++) {
    			fich << zz[i] ;
    			for (int j=0; j<nx; j++) {
    				int index = i*nx*ny + indplane*nx + j ;
    				fich << "  " << uu[index] ;
    			}
    			fich << endl ;
    		}
    		break ;
    	}

    	case 6 : {	

    		for (int i=0; i<nz; i++) {
    			fich << zz[i] ;
    			for (int j=0; j<ny; j++) {
    				int index = i*nx*ny + j*nx + indplane ;
    				fich << "  " << uu[index] ;
    			}
    			fich << endl ;
    		}
    		break ;
    	}

	default : {
		cout << "visu_1d: Invalid type of profiles !" << endl ;
		abort() ;
	}


    }
    fich.close() ;


    // Clean exit
    // ----------

    delete [] xx ;
    delete [] yy ;
    delete [] zz ;

    return EXIT_SUCCESS ;

}

