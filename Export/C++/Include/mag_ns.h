/*
 *  Definition of class Mag_NS (magnetized neutron star exportation)
 *
 */

/*
 *   Copyright (c) 2002  Eric Gourgoulhon
 *   Copyright (c) 2009  Jerome Novak
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

#ifndef __MAG_NS_H_
#define __MAG_NS_H_

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:54:05  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:25  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2009/11/19 16:15:21  j_novak
 * Export class for magnetized neutron stars.
 *
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

namespace Lorene {
/**
 * Magnetized neutron star configuration on a Cartesian grid.
 *
 * A magnetized neutron star is constructed on a Cartesian grid from
 * data stored in a file resulting from a computation obtained from the 
 * code magstar following the Bocquet et al. (1995) paper.
 *
 * Importation of Lorene data is performed by means of the constructor
 * {\tt Mag\_NS::Mag\_NS(int, const double*, const double*, const double*, const char*)}.
 * This constructor takes general arrays for the location of the Cartesian coordinates
 * $(x, y, z)$, i.e. it does not assume that the grid is a uniform one. Note also
 * that these arrays are 1-D, as well as all the metric fields,
 * in order to be use with any ordering of the 3-D storage. For the definitions of units,
 * see the file Lorene/C++/Include/unites.h
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
 *	    Mag_NS star_mag(nb_points, xg, yg, zg, datafile) ;
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

class Mag_NS {

    // Data :
    // -----
    public:
        /// Eos name star
        char eos_name[100] ;

        /// Adiabatic index of EOS if it is polytropic (0 otherwise)
        double gamma_poly ;

        /**
	 *  Polytropic constant of EOS if it is polytropic (0 otherwise)
	 *  [unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$]
	 */
        double kappa_poly ;

	/// Rotation frequency [unit: rad/s] (only rigid rotation)
	double omega ;

	/// Central density [unit: kg m${}^{-3}$]
	double rho_c ;

	/// Central specific internal energy [unit: c${}^2$]
	double eps_c ;

	/// Baryon mass [unit: solar mass]
	double mass_b ;

	/// Gravitational mass [unit: solar mass]
	double mass_g ;

	/// Coordinate equatorial radius [unit: km]
	double r_eq ;

	/// Coordinate polar radius [unit: km]
	double r_p ;

	/// Angular momentum [unit: $G M_{\textrm{sol}}^2/c$]
	double angu_mom ;

	/// Ratio T/W
	double T_over_W ;

	/// Magnetic momentum [unit: A m${}^2$]
	double magn_mom ;

	/// Magnetic field at the pole [unit: $10^9$ T]
	double b_z_pole ;

	/// Magnetic field at the equator [unit: $10^9$ T]
	double b_z_eq ;

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

	/// Component $K^{xx}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
	double* k_xx ;

	/// Component $K^{xy}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
	double* k_xy ;

	/// Component $K^{xz}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
	double* k_xz ;

	/// Component $K^{yy}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
	double* k_yy ;

	/// Component $K^{yz}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
	double* k_yz ;

	/// Component $K^{zz}$ of the extrinsic curvature at the grid points (1-D array) [unit: c/km]
	double* k_zz ;

	// Magneto-Hydro components
	//-------------------------
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

	/** Component $B^x$ of the magnetic field, at the {\tt np} grid points (1-D array) 
	 * [unit: $10^9$ T]
         */
        double* bb_x ;

	/** Component $B^y$ of the magnetic field, at the {\tt np} grid points (1-D array) 
	 * [unit: $10^9$ T]
         */
        double* bb_y ;

	/** Component $B^z$ of the magnetic field, at the {\tt np} grid points (1-D array) 
	 * [unit: $10^9$ T]
         */
        double* bb_z ;

	/** Component $j^t$ of the 4-current, at the {\tt np} grid points (1-D array)
	 * [unit: A m${}^{-2$}]
         */
        double* jj_t ;

	/** Component $j^x$ of the 4-current, at the {\tt np} grid points (1-D array) 
	 * [unit: A m${}^{-2$}]
         */
        double* jj_x ;

	/** Component $j^y$ of the 4-current, at the {\tt np} grid points (1-D array)
	 * [unit: A m${}^{-2$}]
         */
        double* jj_y ;

	/** Component $j^z$ of the 4-current, at the {\tt np} grid points (1-D array)
	 * [unit: A m${}^{-2$}]
         */
        double* jj_z ;

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
	Mag_NS(int nbpoints, const double* xi, const double* yi,
	       const double* zi, const char* filename) ;
	

	/** Constructor from a binary file 
	 *   (previously created by {\tt save\_bin})
	 */
	Mag_NS(FILE* ) ;

	/** Constructor from a formatted file
	 *   (previously created by {\tt save\_form})
	 */
	Mag_NS(ifstream& ) ;

	/// Destructor
	~Mag_NS() ;


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
	friend ostream& operator<<(ostream& , const Mag_NS& ) ;

};

}
#endif
