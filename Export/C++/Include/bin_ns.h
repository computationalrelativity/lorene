/*
 *  Definition of class Bin_NS (binary neutron star exportation)
 *
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

#ifndef __BIN_NS_H_
#define __BIN_NS_H_

/*
 * $Id$
 * $Log$
 * Revision 1.6  2014/10/06 15:13:25  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2010/07/14 16:47:30  e_gourgoulhon
 * Corrected error in the documentation for K_xx, K_xy, etc...:
 * the components are the covariant ones, not the contravariant ones.
 *
 * Revision 1.4  2004/10/20 15:01:37  e_gourgoulhon
 * Corrected error in the comments on the shift vector:
 * corotating coordinates -> non rotating coordinates.
 *
 * Revision 1.3  2003/10/24 15:49:03  e_gourgoulhon
 * Updated documentation.
 *
 * Revision 1.2  2003/01/09 11:08:00  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1  2002/01/11 17:03:02  e_gourgoulhon
 * Exportation of binary neutron stars configuration to a Cartesian grid
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdio>

#include <iostream>
#include <fstream>

using namespace std ;

/**
 * Binary neutron star configuration on a Cartesian grid.
 *
 * A binary black hole system is constructed on a Cartesian grid from
 * data stored in a file resulting from a computation by Taniguchi
 * and Gourgoulhon.
 *
 * Importation of Lorene data is performed by means of the constructor
 * {\tt Bin\_NS::Bin\_NS(int, const double*, const double*, const double*, const char*)}.
 * This constructor takes general arrays for the location of the Cartesian coordinates
 * $(x, y, z)$, i.e. it does not assume that the grid is a uniform one. Note also
 * that these arrays are 1-D, as well as all the metric fields,
 * in order to be use with any ordering of the 3-D storage.
 *
 *  This class is very simple, with all data members being public.
 *  A typical example of use is the following one
 *
 *  \begin{verbatim}
 *	    // Define the Cartesian grid by means of the arrays xg, yg, zg:
 *	    for (int i=0; i<nb_points; i++) {
 *           xg[i] = ...
 *           yg[i] = ...
 *           zg[i] = ...
 *	    }
 *
 *	    // Read the file containing the spectral data and evaluate
 *	    //  all the fields on the Cartesian grid :
 *
 *	    Bin_NS binary_system(nb_points, xg, yg, zg, datafile) ;
 *
 *	    // Extract what you need :
 *
 *	    double* gamma_xx = binary_system.g_xx ; // metric coefficient g_xx
 *
 *	    double* shift_x = binary_system.beta_x ; // x comp. of shift vector
 *
 *	    ...
 *
 *	    // Save everything in an ASCII file :
 *
 *	    ofstream file_ini("ini.d") ;
 *	    binary_system.save_form(file_ini) ;
 *	    file_ini.close() ;
 *
 *  \end{verbatim}
 *
 * @version #$Id$#
 */

class Bin_NS {

    // Data :
    // -----
    public:
        /// Eos name star 1
        char eos_name1[100] ;

        /// Adiabatic index of EOS 1 if it is polytropic (0 otherwise)
        double gamma_poly1 ;

        /**
	 *  Polytropic constant of EOS 1 if it is polytropic (0 otherwise)
	 *  [unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$]
	 */
        double kappa_poly1 ;

        /// Eos name star 2
        char eos_name2[100] ;

        /// Adiabatic index of EOS 2 if it is polytropic (0 otherwise)
        double gamma_poly2 ;

        /**
	 *  Polytropic constant of EOS 2 if it is polytropic (0 otherwise)
	 *  [unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$]
	 */
        double kappa_poly2 ;

	/// Orbital angular velocity [unit: rad/s]
	double omega ;

	/** Distance between the centers (maxiumum density) of the two neutron
	 *  stars [unit: km]
	 */
	double dist ;

	/** Distance between the center of masses of two neutron stars
	 *  [unit: km]
	 */
	double dist_mass ;

	/** Baryon mass of star 1 (less massive star)
	 *  [unit: $M_\odot$]
	 */
	double mass1_b ;

	/** Baryon mass of star 2 (massive star)
	 *  [unit: $M_\odot$]
	 */
	double mass2_b ;

	/** ADM mass of the binary system
	 *  [unit: $M_\odot$]
	 */
	double mass_adm ;

	/** Total angular momentum of the binary system
	 *  [unit: $GM_\odot^2/c$]
	 */
	double angu_mom ;

	/** Coordinate radius of star 1 (less massive star)
	 *  parallel to the x axis toward the companion star
	 *  [unit: km]
	 */
	double rad1_x_comp ;

	/** Coordinate radius of star 1 (less massive star)
	 *  parallel to the y axis [unit: km].
	 */
	double rad1_y ;

	/** Coordinate radius of star 1 (less massive star)
	 *  parallel to the z axis [unit: km].
	 */
	double rad1_z ;

	/** Coordinate radius of star 1 (less massive star)
	 *  parallel to the x axis opposite to the companion star
	 *  [unit: km].
	 */
	double rad1_x_opp ;

	/** Coordinate radius of star 2 (massive star)
	 *  parallel to the x axis toward the companion star
	 *  [unit: km].
	 */
	double rad2_x_comp ;

	/** Coordinate radius of star 2 (massive star)
	 *  parallel to the y axis [unit: km].
	 */
	double rad2_y ;

	/** Coordinate radius of star 2 (massive star)
	 *  parallel to the z axis [unit: km].
	 */
	double rad2_z ;

	/** Coordinate radius of star 2 (massive star)
	 *  parallel to the x axis opposite to the companion star
	 *  [unit: km].
	 */
	double rad2_x_opp ;


	/// Total number of grid points
	int np ;

	/// 1-D array storing the values of coordinate x of the {\tt np} grid points [unit: km]
	double* xx ;

	/// 1-D array storing the values of coordinate y of the {\tt np} grid points [unit: km]
	double* yy ;

	/// 1-D array storing the values of coordinate z of the {\tt np} grid points [unit: km]
	double* zz ;

	/// Lapse function $N$ at the {\tt np} grid points (1-D array)
	double* nnn ;

	/// Component $\beta^x$ of the shift vector of non rotating coordinates [unit: $c$]
	double* beta_x ;
	
	/// Component $\beta^y$ of the shift vector of non rotating coordinates [unit: $c$]
	double* beta_y ; 
	
	/// Component $\beta^z$ of the shift vector of non rotating coordinates [unit: $c$]
	double* beta_z ; 
	
	/// Metric coefficient $\gamma_{xx}$ at the grid points (1-D array)
	double* g_xx ; 

	/// Metric coefficient $\gamma_{xy}$ at the grid points (1-D array)
	double* g_xy ;

	/// Metric coefficient $\gamma_{xz}$ at the grid points (1-D array)
	double* g_xz ; 

	/// Metric coefficient $\gamma_{yy}$ at the grid points (1-D array)
	double* g_yy ; 

	/// Metric coefficient $\gamma_{yz}$ at the grid points (1-D array)
	double* g_yz ; 

	/// Metric coefficient $\gamma_{zz}$ at the grid points (1-D array)
	double* g_zz ; 

	/// Component $K_{xx}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
	double* k_xx ;

	/// Component $K_{xy}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
	double* k_xy ;

	/// Component $K_{xz}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
	double* k_xz ;

	/// Component $K_{yy}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
	double* k_yy ;

	/// Component $K_{yz}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
	double* k_yz ;

	/// Component $K_{zz}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
	double* k_zz ;

	// Hydro components
	//------------------
	/** Baryon density in the fluid frame at the {\tt np} grid points (1-D array)
         * [unit: ${\rm kg \, m}^{-3}$]
         */
        double* nbar ;

	/// Specific internal energy at the  {\tt np} grid points (1-D array) [unit: $c^2$]
	double* ener_spec ;

	/** Component $U^x$ of the fluid 3-velocity with respect to the Eulerian
         * observer, at the {\tt np} grid points (1-D array) [unit: $c$]
         */
        double* u_euler_x ;

	/** Component $U^y$ of the fluid 3-velocity with respect to the Eulerian
         * observer, at the {\tt np} grid points (1-D array) [unit: $c$]
         */
        double* u_euler_y ;

	/** Component $U^z$ of the fluid 3-velocity with respect to the Eulerian
         * observer, at the {\tt np} grid points (1-D array) [unit: $c$]
         */
        double* u_euler_z ;


    // Constructors - Destructor
    // -------------------------
    public:
	/** Constructor from Lorene spectral data.
	 *
	 * This constructor takes general arrays {\tt xi, yi, zi}
	 * for the location of the Cartesian coordinates
	 * $(x, y, z)$, i.e. it does not assume that the grid is a uniform one.
	 * These arrays are 1-D to deal with any ordering of a 3-D storage.
	 *
	 *  @param nbpoints [input] Total number of grid points
	 *  @param xi [input] 1-D array (size {\tt nbpoints}) storing the
	 *		values of coordinate x of the grid points [unit: km]
	 *  @param yi [input] 1-D array (size {\tt nbpoints}) storing the
	 *		values of coordinate y of the grid points [unit: km]
	 *  @param zi [input] 1-D array (size {\tt nbpoints}) storing the
	 *		values of coordinate z of the grid points [unit: km]
	 *  @param filename [input] Name of the (binary) file containing the result
	 *		of a computation by means of the multi-domain 
	 *		spectral method.
	 */
	Bin_NS(int nbpoints, const double* xi, const double* yi,
	       const double* zi, const char* filename) ;
	

	/** Constructor from a binary file 
	 *   (previously created by {\tt save\_bin})
	 */
	Bin_NS(FILE* ) ;

	/** Constructor from a formatted file
	 *   (previously created by {\tt save\_form})
	 */
	Bin_NS(ifstream& ) ;

	/// Destructor
	~Bin_NS() ;


    // Memory management
    // -----------------
    private:

	/// Allocate the memory for the arrays g\_ij, k\_ij, etc...
	void alloc_memory() ;

    // Outputs
    // -------
    public:
	/** Save in a binary file.
	 *  This file can be subsenquently read by the evolution code,
	 *  or by the constructor {\tt Bin\_NS::Bin\_NS(FILE* )}.
	 */
	void save_bin(FILE* ) const ;

	/** Save in a formatted file.
	 *  This file can be subsenquently read by the evolution code,
	 *  or by the constructor {\tt Bin\_NS::Bin\_NS(ifstream\& )}.
	 */
	void save_form(ofstream& ) const ;

	/// Display
	friend ostream& operator<<(ostream& , const Bin_NS& ) ;

};

#endif
