/*
 *  Definition of Lorene class Grille_cart
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

#ifndef __Grille_cart_H_ 
#define __Grille_cart_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/03/07 15:41:12  n_chamel
 * New class for dealing with Cartesian grids
 * Added the sampling type UNIFORM in type_parite.h
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdio.h>



/**
 * Class for Cartesian cubic grid.
 *
 * @version #$Id$#
 */
class Grille_cart {

    // Data :
    // -----
    protected:

        int nx ; /// Number of grid points in the x direction
        int ny ; /// Number of grid points in the y direction
        int nz ; /// Number of grid points in the z direction
        int type_x ; /// Type of sampling in the x direction
        int type_y ; /// Type of sampling in the y direction
        int type_z ; /// Type of sampling in the z direction
        double* x ; /// Array (size {\tt nx}) of the values of x at the grid points
        double* y ; /// Array (size {\tt ny}) of the values of y at the grid points
        double* z ; /// Array (size {\tt nz}) of the values of z at the grid points

    // Constructors - Destructor
    // -------------------
    public:
                /**   Standard constructor
                 *    @param    nx_i    Number of grid points in the x direction
                 *    @param    ny_i    Number of grid points in the y direction
                 *    @param    nz_i    Number of grid points in the z direction
                 *    @param    type_x_i    Type of sampling in the x direction \\
                 *                      {\tt UNIFORM:}  uniform sampling        \\
                 *                      {\tt FIN:}  dense sampling near the boundaries (for Chebyshev expansions)        \\
                 *                      {\tt RARE:}  dense sampling near the right hand side boundary only  (for Chebyshev expansions)
                 *    @param    type_y_i    Type of sampling in the y direction (see {\tt tupe\_x\_i})
                 *    @param    type_z_i    Type of sampling in the z direction  (see {\tt tupe\_x\_i})
                 *
                 */
	Grille_cart(int nx_i, int ny_i, int nz_i, int type_x_i, int type_y_i, int type_z_i) ;

	Grille_cart(const Grille_cart& ) ;		/// Copy constructor

	/// Constructor from a file (see {\tt sauve(FILE* )})
	Grille_cart(FILE* ) ;

	virtual ~Grille_cart() ;			/// Destructor


    // Accessors
    // ---------
    public:
                int get_nx() const {return nx;} ;       /// Returns the number of grid points in the x direction
                int get_ny() const {return ny;} ;       /// Returns the number of grid points in the y direction
                int get_nz() const {return nz;} ;       /// Returns the number of grid points in the z direction
                int get_type_x() const {return type_x;} ;       /// Returns the type of sampling in the x direction
                int get_type_y() const {return type_y;} ;       /// Returns the type of sampling in the y direction
                int get_type_z() const {return type_z;} ;       /// Returns the type of sampling in the z direction
                double get_x(int i)  const {return x[i];} ;       /// Returns the {\tt i}th grid point in the x direction
                double get_y(int i)  const {return y[i];} ;       /// Returns the {\tt i}th grid point in the y direction
                double get_z(int i)  const {return z[i];} ;       /// Returns the {\tt i}th grid point in the z direction

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    /// Save in a file

	/// Display
	friend ostream& operator<<(ostream& , const Grille_cart& ) ;

};

#endif
