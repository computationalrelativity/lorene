/*
 *  Methods of class Grille_cart
 *
 */

/*
 *   Copyright (c) 2002 Nicolas Chamel
 *   Copyright (c) 2002 Eric Gourgoulhon
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

char grille_cart_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/03/07 15:37:05  n_chamel
 * First version
 *
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include <iostream.h>

// C headers
#include <stdlib.h>
#include <math.h>
#include <assert.h>


// Lorene headers
#include "grille_cart.h"
#include "type_parite.h"
#include "utilitaires.h"

                                                                //-------------------------//
                                                                //              Constructors                      //
                                                               //-------------------------//


// Standard constructor
// ---------------
Grille_cart::Grille_cart(int nx_i, int ny_i, int nz_i, int type_x_i, int type_y_i, int type_z_i)  :
        nx(nx_i),
        ny(ny_i),
        nz(nz_i),
        type_x(type_x_i),
        type_y(type_y_i),
        type_z(type_z_i)
        {

        x = new double [nx] ;
        y = new double [ny] ;
        z = new double [nz] ;

        // Computation of x[i]
        // --------------

        switch (type_x) {

                case UNIFORM : {

                        if (nx == 1) {
                                x[0] = 0.  ;
                        }
                        else {
                                double dx = 2. / double(nx-1) ;
                                for (int i=0; i<nx; i++) {
                                        x[i] = -1. + i * dx  ;
                                }
                        }
                        break ;
                }

                case FIN : {
                        assert( nx > 1) ;
                        double xx =  M_PI/double(nx-1) ;
                        for (int i=0; i<nx; i++) {
                                x[i] = -cos(xx*i)  ;
                        }
                        break ;
                }

               case RARE : {
                        assert( nx > 1) ;
                        double xx =  M_PI/double(2*(nx-1)) ;
                        for (int i=0; i<nx; i++) {
                                x[i] = sin(xx*i)  ;
                        }
                        break ;
                }

                default : {
                        cout << "  Grille_cart::Grille_cart : unknown sampling type in x !" << endl ;
                        abort() ;
                }
        }

       // Computation of y[i]
        // --------------

        switch (type_y) {

                case UNIFORM : {

                        if (ny == 1) {
                                y[0] = 0.  ;
                        }
                        else {
                                double dy = 2. / double(ny-1) ;
                                for (int i=0; i<ny; i++) {
                                        y[i] = -1. + i * dy  ;
                                }
                        }
                        break ;
                }

                case FIN : {
                        assert( ny > 1) ;
                        double yy =  M_PI/double(ny-1) ;
                        for (int i=0; i<ny; i++) {
                                y[i] = -cos(yy*i)  ;
                        }
                        break ;
                }

               case RARE : {
                        assert( ny > 1) ;
                        double yy =  M_PI/double(2*(ny-1)) ;
                        for (int i=0; i<ny; i++) {
                                y[i] = sin(yy*i)  ;
                        }
                        break ;
                }

                default : {
                        cout << "  Grille_cart::Grille_cart : unknown sampling type in y !" << endl ;
                        abort() ;
                }
        }

       // Computation of z[i]
        // --------------

        switch (type_z) {

                case UNIFORM : {

                        if (nz == 1) {
                                z[0] = 0.  ;
                        }
                        else {
                                double dz = 2. / double(nz-1) ;
                                for (int i=0; i<nz; i++) {
                                        z[i] = -1. + i * dz  ;
                                }
                        }
                        break ;
                }

                case FIN : {
                        assert( nz > 1) ;
                        double zz =  M_PI/double(nz-1) ;
                        for (int i=0; i<nz; i++) {
                                z[i] = -cos(zz*i)  ;
                        }
                        break ;
                }

               case RARE : {
                        assert( nz > 1) ;
                        double zz =  M_PI/double(2*(nz-1)) ;
                        for (int i=0; i<nz; i++) {
                                z[i] = sin(zz*i)  ;
                        }
                        break ;
                }

                default : {
                        cout << "  Grille_cart::Grille_cart : unknown sampling type in z !" << endl ;
                        abort() ;
                }
        }
}

// Copy constructor
// ------------
Grille_cart::Grille_cart(const Grille_cart& grid)  :
        nx(grid.nx),
        ny(grid.ny),
        nz(grid.nz),
        type_x(grid.type_x),
        type_y(grid.type_y),
        type_z(grid.type_z)
{

        x = new double [nx]  ;
        y = new double [ny]  ;
        z = new double [nz]  ;

        for (int i=0; i<nx; i++) {
                x[i] = grid.x[i] ;
        }

        for (int i=0; i<ny; i++) {
                y[i] = grid.y[i] ;
        }

        for (int i=0; i<nz; i++) {
                z[i] = grid.z[i] ;
        }

}


// Constructor from a file
// ----------------
Grille_cart::Grille_cart(FILE* fich)  {

                fread_be(&nx, sizeof(int), 1, fich) ;
                fread_be(&ny, sizeof(int), 1, fich) ;
                fread_be(&nz, sizeof(int), 1, fich) ;
                fread_be(&type_x, sizeof(int), 1, fich) ;
                fread_be(&type_y, sizeof(int), 1, fich) ;
                fread_be(&type_z, sizeof(int), 1, fich) ;


}

                                                         //--------------------//
                                                        //             Destructor                //
                                                       //--------------------//

Grille_cart::~Grille_cart(){

        delete [] x ;
        delete [] y ;
        delete [] z ;

}
                                                        //--------------------//
                                                        //              Outputs                   //
                                                       //--------------------//

// Write in a binary file
void Grille_cart::sauve(FILE* fich) const  {

                fwrite_be(&nx, sizeof(int), 1, fich) ;
                fwrite_be(&ny, sizeof(int), 1, fich) ;
                fwrite_be(&nz, sizeof(int), 1, fich) ;
                fwrite_be(&type_x, sizeof(int), 1, fich) ;
                fwrite_be(&type_y, sizeof(int), 1, fich) ;
                fwrite_be(&type_z, sizeof(int), 1, fich) ;

}


// Formatted output
        ostream& operator<<(ostream& ost, const Grille_cart& grid)  {

                ost << "Cartesian grid:" << endl ;
                ost << " Number of points and type of sampling in x: " << grid.get_nx() << "         "
                        << grid.get_type_x() << endl ;
                ost << " Number of points and type of sampling in y: " << grid.get_ny() << "         "
                        << grid.get_type_y() << endl ;
                ost << " Number of points and type of sampling in z: " << grid.get_nz() << "         "
                        << grid.get_type_z() << endl ;

                return ost ;
}






