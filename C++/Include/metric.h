 /*
 *  Definition of Lorene class Metric
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Philippe Grandclement (for previous class Metrique)
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

#ifndef __METRIC_H_ 
#define __METRIC_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.9  2014/10/13 08:52:35  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2004/11/18 12:23:42  jl_jaramillo
 * radial vector normal to a spherical slicing and pointing towards
 * spatial infinity
 *
 * Revision 1.7  2004/03/22 13:12:42  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.6  2003/12/30 23:04:58  e_gourgoulhon
 * Important reorganization of class Metric:
 *  -- suppression of virtual methods fait_* : the actual computations
 *     are now performed via the virtual methods con(), cov(), connect(),
 *     ricci(), ricci_scal(), determinant()
 *  -- the member p_connect is now treated as an ordinary derived data
 *     member
 *  -- the construction of the associated connection (member p_connect)
 *     is performed thanks to the new methods Map::flat_met_spher() and
 *     Map::flat_met_cart().
 *
 * Revision 1.5  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.4  2003/10/06 15:30:32  j_novak
 * Defined methods for flat metric.
 *
 * Revision 1.3  2003/10/06 13:58:45  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.2  2003/10/03 11:21:45  j_novak
 * More methods for the class Metric
 *
 * Revision 1.1  2003/10/02 15:45:48  j_novak
 * New class Metric
 *
 *
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "connection.h"

#define N_TENSOR_DEPEND 200

namespace Lorene {
/**
 * Metric for tensor calculation. \ingroup (tensor)
 * 
 */
class Metric {

    // Data : 
    // -----
    protected:
	const Map* const mp ;	///< Reference mapping.

	/**
	 * Pointer on the contravariant representation.
	 */
	mutable Sym_tensor* p_met_cov ;

	/**
	 * Pointer on the covariant representation.
	 */
	mutable Sym_tensor* p_met_con ;


    // Derived data : 
    // ------------
    protected:
    
	mutable Connection* p_connect ; ///< Connection associated with the metric

	/**
	 * Pointer on the Ricci scalar.
         * Remark: the Ricci tensor is stored in the connection (member
         *  \c p_connect->p_ricci ).
	 */
	mutable Scalar* p_ricci_scal ;

	/**
	 * Pointer to the radial vector normal to a spherical slicing and pointing 
	 * toward spatial infinity
	 */
	mutable Vector* p_radial_vect ;


	
	/**
	 * Pointer on the determinant.
	 */
	mutable Scalar* p_determinant ;
	
	/**
	 * Pointer on the dependancies, that means the array contains pointers
	 * on all the \c Tensor  whom derivative members have been calculated
	 * using \c *this .
	 */
	mutable const Tensor* tensor_depend[N_TENSOR_DEPEND] ;
	
    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor from a \c Sym_tensor .
	 *
	 *  The symmetric tensor can be either the covariant or
	 *  the contravariant representation of the metric.
	 */
	explicit Metric(const Sym_tensor& tens) ;  
         
	Metric(const Metric& met) ;		///< Copy constructor

	/// Constructor from a file (see \c sauve(FILE*) )
	Metric(const Map&, FILE* ) ;    		

    protected:
	/// Simplified constructor used by derived classes.
	explicit Metric(const Map& mpi) ;

    public:
	virtual ~Metric() ;			///< Destructor
 

    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	void del_deriv() const ; 
	
	/// Sets to \c 0x0  all the pointers on derived quantities
	void set_der_0x0() const ; 

	/**
	 * Deletes all the derivative members of the \c Tensor  contained in
	 * \c tensor_depend . Those quantities had been previously 
	 * calculated using \c *this .
	 */
	void del_tensor_depend() const ;
		
	///Sets all elements of \c tensor_depend  to 0x0.
	void set_tensor_depend_0x0() const ;
		

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Metric
	void operator=(const Metric& met) ;	

	/**
	 * Assignment from a \c Sym_tensor .
	 * The allocated representation depends on the type of the
	 * input tensor indices.
	 * All the other members are deleted.
	 */
	virtual void operator=(const Sym_tensor& tens) ;
	
    // Accessors
    // ---------
    public:
	/// Returns the mapping
	const Map& get_mp() const {return *mp ; } ;

	/// Read-only access to the covariant representation
	virtual const Sym_tensor& cov() const ;

	/// Read-only access to the contravariant representation
	virtual const Sym_tensor& con() const ;

	/// Returns the connection
	virtual const Connection& connect() const ;

	/** Returns the Ricci tensor (given by the \c Connection  
         *  \c p_connect )
         */
	const Sym_tensor& ricci() const ;
	
	/// Returns the Ricci scalar
	virtual const Scalar& ricci_scal() const ;

	/** Returns the  radial vector normal to a spherical slicing and pointing 
	 * toward spatial infinity
	 */

	virtual const Vector& radial_vect() const ;


	/**Returns the determinant.
	 * 
	 * This determinant is stored as a \c Scalar  although it
	 * a scalar density. To be a real scalar it must be divided
	 * by e.g. the determinant of a flat metric.
	 */
	virtual const Scalar& determinant() const ;





    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Metric& ) ;	

    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    


	friend class Tensor ;

};

/**
 * Flat metric for tensor calculation.\ingroup (tensor)
 * 
 */
class Metric_flat: public Metric {

    // Data : 
    // -----
    protected:
	
  /** Vectorial basis (triad) with respect to which the components of
   * the flat metric are defined.
   */
  const Base_vect* triad ; 

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor.
	 *
	 *  Standard constructor from a mapping and a triad.
	 */
	Metric_flat(const Map&, const Base_vect& ) ;   

	Metric_flat(const Metric_flat& ) ;		///< Copy constructor

	/// Constructor from a file (see \c sauve(FILE*) )
	Metric_flat(const Map&, FILE* ) ;    		

   public:
	virtual ~Metric_flat() ;			///< Destructor
 

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Metric_flat
	void operator=(const Metric_flat&) ;	

	/**
	 * Assignment from a \c Sym_tensor .
	 * In principle, this method should not be used for a \c Metric_flat .
	 */
	virtual void operator=(const Sym_tensor& tens) ;


    // Accessors
    // ---------
    public:
	/** Returns the vectorial basis (triad) on which the metric
	 *  is defined.  
	 */
	const Base_vect* get_triad() const {return triad;} ; 
    
	/// Read-only access to the covariant representation
	virtual const Sym_tensor& cov() const ;

	/// Read-only access to the contravariant representation
	virtual const Sym_tensor& con() const ;

	/// Returns the connection
	virtual const Connection& connect() const ;

	/// Returns the Ricci scalar
	virtual const Scalar& ricci_scal() const ;

	/**Returns the determinant.
	 * 
	 * This determinant is stored as a \c Scalar  although it
	 * a scalar density. To be a real scalar it must be divided
	 * by e.g. the determinant of a flat metric.
	 */
	virtual const Scalar& determinant() const ;


    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    


};



}
#endif
