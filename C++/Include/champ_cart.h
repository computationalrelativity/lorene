/*
 *  Definition of Lorene class Champ_cart
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

#ifndef __Champ_cart_H_
#define __Champ_cart_H_

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/03/28 10:47:29  n_chamel
 * New class Champ_cart for fields on Cartesian grids.
 *
 *
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "map_cart.h"
#include "tbl.h"

/**   Scalar field on a Cartesian grid.
 *
 *
 * @version #$Id$#
 */
class Champ_cart{

    // Data :
    // -----
    protected:
       /// Pointer on the mapping (class {\tt Map\_cart}) on which {\tt this} is defined
        const Map_cart& map ;

        /**   Logical state of the field.
          *       {\tt etat = ETATNONDEF} : undefined state (typically uninitialized)    \\
          *       {\tt etat = ETATZERO} : null field                                                           \\
          *       {\tt etat = ETATQCQ} : ordinary state
          */
        int etat ;

        /** Bases on which the spectral expansions are performed.
          *     {\tt base[0]} : basis functions for the x expansion      \\
          *     {\tt base[1]} : basis functions for the y expansion      \\
          *     {\tt base[2]} : basis functions for the z expansion      \\
           * (see Lorene class {\tt  Base\_val} for the naming convention  of the various types of bases).
          */
        mutable int base[3] ;

        /// Pointer on the {\tt Tbl} containing the values of the field at the grid points
        Tbl*  val ;

        ///  Pointer on the {\tt Tbl} containing the spectral coefficients of the field
        Tbl* cf ;

    // Derived data :
    // ------------
    // protected:
    //	mutable Champ_cart* p_dx ;   /// Dervative d/dx

    // Constructors - Destructor
    // ------------------
    public:
                Champ_cart(const Map_cart& ) ;			/// Standard constructor
                Champ_cart(const Champ_cart& ) ;		/// Copy constructor

	/// Constructor from a file (see {\tt sauve(FILE* )})
                Champ_cart(const Map_cart&, FILE* ) ;

	virtual ~Champ_cart() ;			/// Destructor


    // Memory management
    // -----------------
    protected:
	/// Logical destructor
	virtual void del_all() ;

	/// Deletes all the derived quantities
	// virtual void del_deriv() const ;

	/// Sets to {\tt 0x0} all the pointers on derived quantities
	// virtual void set_der_0x0() const ;


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another {\tt Champ\_cart}
	void operator=(const Champ_cart&) ;

    // Accessors
    // ---------
    public:
                const Map_cart& get_map() const  {return map;} ;     /// Returns the mapping on which the field is defined

        /**   Returns the logical state of the field.
          *       {\tt etat = ETATNONDEF} : undefined state (typically uninitialized)    \\
          *       {\tt etat = ETATZERO} : null field                                                           \\
          *       {\tt etat = ETATQCQ} : ordinary state
          */
                int get_etat() const {return etat;};

        /** Returns the bases on which the spectral expansions are performed.
          *     {\tt base[0]} : basis functions for the x expansion      \\
          *     {\tt base[1]} : basis functions for the y expansion      \\
          *     {\tt base[2]} : basis functions for the z expansion      \\
           * (see Lorene class {\tt  Base\_val} for the naming convention  of the various types of bases).
          */
        const int* get_base() const {return base;} ;

        /// Returns the {\tt Tbl} containing the values of the field at the grid points
        const Tbl& get_val() const {return *val;}  ;

        ///  Returns the {\tt Tbl} containing the spectral coefficients of the field
        const Tbl& get_cf() const {return *cf;}  ;


    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    /// Save in a file

	/// Display
	friend ostream& operator<<(ostream& , const Champ_cart& ) ;

        // Computing functions
        // ---------------
        public:
        /**   Compute the spectral coefficients  from the values at the grid points, according to the
          *    specified bases
          *     @param base_x:  basis functions for the x expansion
          *     @param base_y:  basis functions for the y expansion
          *     @param base_z:  basis functions for the z expansion
            * (see Lorene class {\tt  Base\_val} for the naming convention  of the various types of bases).
          */
               void coef(int base_x, int base_y, int base_z) const ;

         /**   Compute the spectral coefficients  from the values at the grid points.
           *    Assumes that the spectral bases  have been already defined.
           */
              void coef() const ;

         ///   Compute the values at the grid points from the spectral coefficients
              void coef_i() const ;

};

#endif
