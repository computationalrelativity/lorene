/*
 *  Definition of Lorene class Map_cart
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

#ifndef __Map_cart_H_
#define __Map_cart_H_

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/03/15 13:16:23  n_chamel
 * Mapping between grid coordinates and physical coordinates
 *
 *
 *
 *
 * $Header$
 *
 */

// C++ headers

// Lorene headers
#include "grille_cart.h"

                                                                   //-------------------------------//
                                                                   //                   basis class Map_cart                //
                                                                   //-------------------------------//

/**
 * Mapping between Cartesian grid and physical space.
 * Defines the functions
 *   \begin{equation}
 *      x\mapsto X, \quad x\in [-1,1] \ \mbox{or}\ x\in[0,1]      \nonumber
 *   \end{equation}
 *   \begin{equation}
 *      y\mapsto Y, \quad y\in [-1,1] \ \mbox{or}\ y\in[0,1]           \nonumber
  *   \end{equation}
 *   \begin{equation}
 *      z\mapsto Z, \quad z\in [-1,1] \ \mbox{or}\ z\in[0,1]          \nonumber
 *   \end{equation}
 * where $(x,y,z)$ (resp. $(X,Y,Z)$) are the grid (resp. physical) coordinates.
 *
 * The class {\tt Map\_cart} is an abstract one: it cannot be instanciated.
 * Specific implementation of coordinate mappings will be performed by derived
 * classes of {\tt Map\_cart}.
  *
 * @version #$Id$#
 */
class Map_cart{

    // Data :
    // -----
    protected:
        /// Pointer on the grid (class {\tt Grille\_cart}) on which {\tt this} is defined
        const Grille_cart*  grid ;

    // Constructors - Destructor
    // -------------------------
    protected:
	Map_cart(const Grille_cart& ) ;	/// Standard constructor
	Map_cart(const Map_cart& ) ;		/// Copy constructor

              virtual ~Map_cart() ;			/// Destructor

    // Accessors
    // ---------
    public:
                const Grille_cart& get_grid() const  {return *grid;} ;     /// Returns the Cartesian on which the mapping is defined

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const  = 0 ;	    /// Save in a file

	/// Display
	friend ostream& operator<<(ostream& , const Map_cart& ) ;

    private:
	virtual ostream& operator>>(ostream &) const = 0 ;    /// Operator >>


};

                                                                  //-------------------------------//
                                                                   //                   class Map_cart_aff                    //
                                                                   //-------------------------------//

/**
 * Affine mapping between Cartesian grid and physical space:
 *   \begin{equation}
 *      x\mapsto \alpha_1 x + \beta_1, \quad x\in [-1,1] \ \mbox{or}\ x\in[0,1]      \nonumber
 *  \end{equation}
 *  \begin{equation}
 *    y\mapsto \alpha_2 y + \beta_2, \quad  y\in [-1,1] \ \mbox{or}\ y\in[0,1]    \nonumber
  *  \end{equation}
 *  \begin{equation}
*    z\mapsto \alpha_3 z + \beta_3, \quad z\in [-1,1] \ \mbox{or}\ z\in[0,1]      \nonumber
 *  \end{equation}
  * where $(x,y,z)$  are the grid coordinates.
 *
 * @version #$Id$#
 */
class Map_cart_aff : public Map_cart {

    // Data :
    // -----
    protected:
        double alpha1 ;         /// Coefficient  $\alpha_1$
        double beta1 ;         /// Coefficient  $\beta_1$
        double alpha2 ;         /// Coefficient  $\alpha_2$
        double beta2 ;         /// Coefficient  $\beta_2$
        double alpha3 ;         /// Coefficient  $\alpha_3$
        double beta3 ;         /// Coefficient  $\beta_3$

    // Constructors - Destructor
    // -------------------------
    public:
        /**    Standard constructor.
          *     @param grid_i Cartesian grid on which the mapping is defined
          *     @param  alpha1_i   coefficient  $\alpha_1$  of the mapping
          *     @param  beta1_i   coefficient  $\beta_1$  of the mapping
          *     @param  alpha2_i   coefficient  $\alpha_2$  of the mapping
          *     @param  beta2_i   coefficient  $\beta_2$  of the mapping
          *     @param  alpha3_i   coefficient  $\alpha_3$  of the mapping
          *     @param  beta3_i   coefficient  $\beta_3$  of the mapping
          *
          */
	Map_cart_aff(const Grille_cart&  grid_i,   double alpha1_i, double beta1_i,
                                                                                double alpha2_i, double beta2_i,
                                                                                double alpha3_i, double beta3_i) ;

	Map_cart_aff(const Map_cart_aff& ) ;		/// Copy constructor
	Map_cart_aff(const Grille_cart& , FILE* ) ;		/// Constructor  from file

              virtual ~Map_cart_aff() ;			/// Destructor

    // Accessors
    // ---------
    public:
                double get_alpha1() const  {return alpha1;} ;   /// Returns the coefficient  $\alpha_1$
                double get_beta1() const  {return beta1;} ;   /// Returns the coefficient  $\beta_1$
                double get_alpha2() const  {return alpha2;} ;   /// Returns the coefficient  $\alpha_2$
                double get_beta2() const  {return beta2;} ;   /// Returns the coefficient  $\beta_2$
                double get_alpha3() const  {return alpha3;} ;   /// Returns the coefficient  $\alpha_3$
                double get_beta3() const  {return beta3;} ;   /// Returns the coefficient  $\beta_3$

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const  ;	    /// Save in a file

    private:
	virtual ostream& operator>>(ostream &) const  ;    /// Operator >>


};



#endif
