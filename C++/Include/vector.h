/*
 *  Definition of Lorene class Vector
 *
 */

/*
 *   Copyright (c) 2003  Eric Gourgoulhon & Jerome Novak
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

#ifndef __VECTOR_H_ 
#define __VECTOR_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/09/29 12:52:56  j_novak
 * Methods for changing the triad are implemented.
 *
 * Revision 1.1  2003/09/26 08:07:32  j_novak
 * New class vector
 *
 *
 * $Header$
 *
 */


/**
 * Tensor field of valence 1.
 * 
 * @version #$Id$#
 */
class Vector: public Tensor {

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param tipe  the type {\tt COV} for a covariant vector (1-form) 
	 *		and {\tt CON} for a contravariant one
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *		    the vector components are defined 
	 */
	Vector(const Map& map, int tipe, const Base_vect& triad_i) ;

	/** Standard constructor with the triad passed as a pointer.
	 * 
	 * @param map   the mapping 
	 * @param tipe  the type {\tt COV} for a covariant vector (1-form) 
	 *		and {\tt CON} for a contravariant one
	 * @param triad_i  pointer on the vectorial basis (triad) 
	 * with respect to which the vector components are defined 
	 */
	Vector(const Map& map, int tipe, const Base_vect* triad_i) ;

	Vector(const Vector& ) ;       /// Copy constructor

	/** Constructor from a file (see {\tt Tensor::sauve(FILE* )}).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param fich  file which has been created by 
	 *			    the function {\tt sauve(FILE* )}.
	 */
	Vector (const Map& map, const Base_vect& triad_i, FILE* fich) ;

	virtual ~Vector() ;			/// Destructor
 

    // Mutators / assignment
    // ---------------------
    public:
	/** Sets a new vectorial basis (triad) of decomposition and modifies
	 *  the components accordingly. 
	 */
	virtual void change_triad(const Base_vect& ) ; 

	/// Assignment from a Tensor
	virtual void operator=(const Tensor&) ;	
	
    // Accessors
    // ---------
    public:
	Scalar& set(int ) ; ///Read/write access to a component

	const Scalar& operator()(int ) const; ///Readonly access to a component

	/**
	 * Returns the position in the {\tt Scalar} array {\tt cmp} of a 
	 * component given by its index.  
	 *
	 * @return position in the {\tt Scalar} array {\tt cmp}  
	 * corresponding to the index given in {\tt idx}. {\tt idx}
	 * must be a 1-D {\tt Itbl} of size 1, the element of which 
	 * must be one of the spatial indices 1, 2 or 3. 
	 */
	virtual int position(const Itbl& idx) const {
	  assert (idx.get_ndim() == 1) ;
	  assert (idx.get_dim(0) == 1) ;
	  assert ((idx(0) >= 1) && (idx(0) <= 3)) ;

	  return (idx(0) - 1) ;
    
	} ;

	/**
	 * Returns the index of a component given by its position in the 
	 * {\tt Scalar} array {\tt cmp}. 
	 *
	 * @return the index is stored in an 1-D array ({\tt Itbl}) of
	 *         size 1 giving its value for the component located at 
	 *         the position {\tt place} in the {\tt Scalar} array 
	 *         {\tt cmp}. The element of this {\tt Itbl} 
	 *	   corresponds to a spatial index 1, 2 or 3. 
	 */
	virtual Itbl indices(int place) const {
	  assert((place>=0) && (place<3)) ;
	  
	  Itbl res(1) ;
	  res = place + 1;
	  return res ;

	};
	
 
};

#endif
