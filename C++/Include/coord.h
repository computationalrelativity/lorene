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
 * Revision 1.7  2014/10/13 08:52:33  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:09:39  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2004/03/22 13:12:40  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.4  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.3  2002/10/16 14:36:28  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/06/17 14:05:16  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
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
#include <cstdlib>
#include <cstdio>
#include "mtbl.h"

namespace Lorene {
class Map ;

/**
 * Active physical coordinates and mapping derivatives.
 * \ingroup (map)
 *
 */
class Coord {
    
	// Data :
	// -----
    public:
	const Map* mp ;			 ///< Mapping on which the \c Coord  is defined
	Mtbl* (*met_fait)(const Map* ) ; ///< Function to compute the coordinate
	mutable Mtbl* c ;		 ///< The coordinate values at each grid point

	// Constructors, destructor : 
	// ------------------------
    public:
	Coord() ;			    ///< Default constructor
	
	/**
	 * Constructor from a mapping and a method.
	 * @param mp [input] Mapping on which the Coord is defined
	 * @param construct [input] Method to construct the \c Coord , i.e. to
	 *		    initialize the \c Mtbl  which contains the value of
	 *		    the coordinate or mapping derivative represented by 
	 *                  the \c Coord  
	 */
	Coord(const Map* mp, Mtbl* (*construct)(const Map*) ) ;
	

    private:
	/** Copy constructor (private and not implemented to make \c Coord 
	 * a non-copyable class)
	 */ 
	Coord(const Coord & ) ;		    
	
    public: 
	~Coord() ;			    ///< Destructor

	// Various methods :
	// ---------------	

	/** Assignement operator (private and not implemented to make 
	 *   \c Coord  a non-copyable class)
	 */
    private: 
	void operator=(const Coord& ) ;
	 	
    public: 
	/**
	 * Semi-constructor from a mapping and a method.
	 * This function is intended to complete the construction started by
	 * the default constructor.
	 * @param mp [input] Mapping on which the Coord is defined
	 * @param construct [input] Method to construct the \c Coord , i.e. to
	 *		    initialize the \c Mtbl  which contains the value of
	 *		    the coordinate or mapping derivative represented by 
	 *                  the \c Coord  
	 */
	void set(const Map* mp, Mtbl* (*construct)(const Map*) ) ;	    

	/**
	 * Computes, at each point of the grid, the value of the coordinate or 
	 * mapping derivative represented by the \c Coord .
	 * The result is stored in the \c Mtbl  member \c *c . 
	 */
	void fait() const ;	

	/// Logical destructor (deletes the \c Mtbl  member \c *c ).
 	void del_t() const ; 

	friend ostream& operator<<(ostream& , const Coord& ) ;	///< Display 

 } ;
ostream& operator<<(ostream& , const Coord& ) ;	

// Prototypage de l'arithmetique
/**
 * \defgroup coord_ari Coord  Arithmetics.
 * \ingroup (map)
 * @{
 */

Mtbl operator+(const Coord&) ;			///< + \c Coord 
Mtbl operator-(const Coord&) ;			///< \c - \c Coord 

Mtbl operator+(const Coord& a, const Coord& b) ;	///< \c Coord  + \c Coord 
Mtbl operator-(const Coord& a, const Coord& b) ;	///< \c Coord  - \c Coord  
Mtbl operator*(const Coord& a, const Coord& b) ;	///< \c Coord  * \c Coord 

Mtbl operator+(const Coord& a, const Mtbl& b) ;	///< \c Coord  + \c Mtbl 
Mtbl operator-(const Coord& a, const Mtbl& b) ;	///< \c Coord  - \c Mtbl 
Mtbl operator*(const Coord& a, const Mtbl& b) ;	///< \c Coord  * \c Mtbl 

Mtbl operator+(const Mtbl& a, const Coord& b) ;	///< \c Mtbl  + \c Coord 
Mtbl operator-(const Mtbl& a, const Coord& b) ;	///< \c Mtbl  - \c Coord 
Mtbl operator*(const Mtbl& a, const Coord& b) ;	///< \c Mtbl  * \c Coord 
/** @} */

}
#endif
