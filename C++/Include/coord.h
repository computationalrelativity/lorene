/*
 *  Definition of Lorene class Coord
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


#ifndef __COORD_H_ 
#define __COORD_H_ 


/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.1  1999/10/15  09:15:56  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.0  1999/02/15  10:41:51  hyc
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/15  09:59:50  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1999/01/15  09:10:39  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Fichier includes
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include "mtbl.h"

class Map ;

/**
 * Active physical coordinates and mapping derivatives.
 * 
 * @version #$Id$#
 *
 */

class Coord {
    
	// Data :
	// -----
    public:
	const Map* mp ;			 /// Mapping on which the {\tt Coord} is defined
	Mtbl* (*met_fait)(const Map* ) ; /// Function to compute the coordinate
	mutable Mtbl* c ;		 /// The coordinate values at each grid point

	// Constructors, destructor : 
	// ------------------------
    public:
	Coord() ;			    /// Default constructor
	
	/**
	 * Constructor from a mapping and a method.
	 * @param mp [input] Mapping on which the Coord is defined
	 * @param construct [input] Method to construct the {\tt Coord}, i.e. to
	 *		    initialize the {\tt Mtbl} which contains the value of
	 *		    the coordinate or mapping derivative represented by 
	 *                  the {\tt Coord} 
	 */
	Coord(const Map* mp, Mtbl* (*construct)(const Map*) ) ;
	

    private:
	/** Copy constructor (private and not implemented to make {\tt Coord}
	 * a non-copyable class)
	 */ 
	Coord(const Coord & ) ;		    
	
    public: 
	~Coord() ;			    /// Destructor

	// Various methods :
	// ---------------	

	/** Assignement operator (private and not implemented to make 
	 *   {\tt Coord} a non-copyable class)
	 */
    private: 
	void operator=(const Coord& ) ;
	 	
    public: 
	/**
	 * Semi-constructor from a mapping and a method.
	 * This function is intended to complete the construction started by
	 * the default constructor.
	 * @param mp [input] Mapping on which the Coord is defined
	 * @param construct [input] Method to construct the {\tt Coord}, i.e. to
	 *		    initialize the {\tt Mtbl} which contains the value of
	 *		    the coordinate or mapping derivative represented by 
	 *                  the {\tt Coord} 
	 */
	void set(const Map* mp, Mtbl* (*construct)(const Map*) ) ;	    

	/**
	 * Computes, at each point of the grid, the value of the coordinate or 
	 * mapping derivative represented by the {\tt Coord}.
	 * The result is stored in the {\tt Mtbl} member {\tt *c}. 
	 */
	void fait() const ;	

	/// Logical destructor (deletes the {\tt Mtbl} member {\tt *c}).
 	void del_t() const ; 

	friend ostream& operator<<(ostream& , const Coord& ) ;	/// Display 

};

// Prototypage de l'arithmetique
/**
 * @name {\tt Coord} Arithmetics
 */
    //@{
Mtbl operator+(const Coord&) ;			/// + {\tt Coord}
Mtbl operator-(const Coord&) ;			/// - {\tt Coord}

Mtbl operator+(const Coord&, const Coord&) ;	/// {\tt Coord} + {\tt Coord}
Mtbl operator-(const Coord&, const Coord&) ;	/// {\tt Coord} - {\tt Coord} 
Mtbl operator*(const Coord&, const Coord&) ;	/// {\tt Coord} * {\tt Coord}

Mtbl operator+(const Coord&, const Mtbl&) ;	/// {\tt Coord} + {\tt Mtbl}
Mtbl operator-(const Coord&, const Mtbl&) ;	/// {\tt Coord} - {\tt Mtbl}
Mtbl operator*(const Coord&, const Mtbl&) ;	/// {\tt Coord} * {\tt Mtbl}

Mtbl operator+(const Mtbl&, const Coord&) ;	/// {\tt Mtbl} + {\tt Coord}
Mtbl operator-(const Mtbl&, const Coord&) ;	/// {\tt Mtbl} - {\tt Coord}
Mtbl operator*(const Mtbl&, const Coord&) ;	/// {\tt Mtbl} * {\tt Coord}
    //@}

#endif
