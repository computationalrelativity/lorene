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

/**
 * Metric for tensor calculation.
 * 
 * @version #$Id$#
 */
class Metric {

    // Data : 
    // -----
    protected:
	const Map* const mp ;	/// Reference mapping.

	mutable Connection* p_connect ; ///Connection associated with the metric

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
	/**
	 * Pointer on the Ricci scalar.
	 */
	mutable Scalar* p_ricci_scal ;
	
	/**
	 * Pointer on the determinant.
	 */
	mutable Scalar* p_determinant ;
	
	/**
	 * Pointer on the dependancies, that means the array contains pointers
	 * on all the {\tt Tensor} whom derivative members have been calculated
	 * using {\tt *this}.
	 */
	mutable const Tensor* tensor_depend[N_TENSOR_DEPEND] ;
	
    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor from a {\tt Sym\_tensor}.
	 *
	 *  The symmetric tensor can be either the covariant or
	 *  the contravariant representation of the metric.
	 */
	Metric(const Sym_tensor& ) ;   
	Metric(const Metric& ) ;		/// Copy constructor

	/// Constructor from a file (see {\tt sauve(FILE* )})
	Metric(const Map&, FILE* ) ;    		

    protected:
	///Simplified constructor used by derived classes.
	Metric(const Map&) ;

    public:
	virtual ~Metric() ;			/// Destructor
 

    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to {\tt 0x0} all the pointers on derived quantities
	virtual void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Metric
	void operator=(const Metric&) ;	

	/**
	 * Assignment from a {\tt Sym\_tensor}.
	 * The allocated representation depends on the type of the
	 * input tensor indices.
	 * All the other members are deleted.
	 */
	virtual void operator= (const Sym_tensor& sym_t) ;
	
     protected:
	/**
	 * Deletes all the derivative members of the {\tt Tensor} contained in
	 * {\tt tensor_depend}. Those quantities had been previously 
	 * calculated using {\tt *this}.
	 */
	void del_tensor_depend() const ;
		
	///Sets all elements of {\tt tensor_depend} to 0x0.
	void set_tensor_depend_0x0() const ;
		
    // Accessors
    // ---------
    public:
	///Read-only access to the covariant representation
	const Sym_tensor& cov() const {
	  if (p_met_cov == 0x0) fait_cov() ;
	  return *p_met_cov ; } ;

	///Read-only access to the contravariant representation
	const Sym_tensor& con() const {
	  if (p_met_con == 0x0) fait_con() ;
	  return *p_met_con ; };

	///Returns the mapping
	const Map& get_mp() const {return *mp ; } ;

	///Returns the connection
	const Connection& get_connect() const ;

	///Returns the Ricci tensor (given by the {\tt Connection})
	const Sym_tensor& ricci() const ;
	
	///Returns the Ricci scalar
	const Scalar& ricci_scal() const ;

	/**Returns the determinant.
	 * 
	 * This determinant is stored as a {\tt Scalar} although it
	 * a scalar density. To be a real scalar it must be divided
	 * by e.g. the determinant of a flat metric.
	 */
	const Scalar& determinant() const ;


    // Computation of derived members
    // ------------------------------
    protected:
	/**
	 * Calculates, if needed, the covariant representation.
	 * The result is in {\tt *p\_met\_cov}.
	 */
	virtual void fait_cov() const ;

	/**
	 * Calculates, if needed, the contravariant representation.
	 * The result is in {\tt *p\_met\_con}.
	 */
	virtual void fait_con() const ;

	/**
	 * Calculates, when needed, the connection from this metric.
	 * The result is in {\tt *p_connect}.
	 */
	virtual void fait_connection() const ;

	/**
	 * Calculates, if needed, the Ricci-scalar.
	 * The result is in {\tt *p\_ricci\_scal}.
	 */
	virtual void fait_ricci_scal() const ;

	/**
	 * Calculates, if needed, the determinant.
	 * The result is in {\tt *p\_determinant}.
	 */
	virtual void fait_determinant() const ;


    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    /// Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Metric& ) ;	

    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    


	friend class Tensor ;

};
/**
 * Flat metric for tensor calculation.
 * 
 */
class Metric_flat: public Metric {

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor.
	 *
	 *  Standard constructor from a mapping and a triad.
	 */
	Metric_flat(const Map&, const Base_vect& ) ;   
	Metric_flat(const Metric_flat& ) ;		/// Copy constructor

	/// Constructor from a file (see {\tt sauve(FILE* )})
	Metric_flat(const Map&, FILE* ) ;    		

   public:
	virtual ~Metric_flat() ;			/// Destructor
 

    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to {\tt 0x0} all the pointers on derived quantities
	virtual void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Metric_flat
	void operator=(const Metric_flat&) ;	

	/**
	 * Assignment from a {\tt Sym\_tensor}.
	 * In principle, this method should not be used for a {\tt Metric\_flat}.
	 */
	virtual void operator= (const Sym_tensor& sym_t) ;
	

    // Computation of derived members
    // ------------------------------
    protected:
	/**
	 * Calculates, if needed, the contravariant representation.
	 * The result is in {\tt *p\_met\_con}.
	 */
	virtual void fait_con() const ;

	/**
	 * Calculates, if needed, the covariant representation.
	 * The result is in {\tt *p\_met\_cov}.
	 */
	virtual void fait_cov() const ;

	/**
	 * Calculates, when needed, the connection from this metric.
	 * The result is in {\tt *p_connect}.
	 */
	virtual void fait_connection() const ;

	/**
	 * Calculates, if needed, the Ricci-scalar.
	 * The result is in {\tt *p\_ricci\_scal}.
	 */
	virtual void fait_ricci_scal() const ;

	/**
	 * Calculates, if needed, the determinant.
	 * The result is in {\tt *p\_determinant}.
	 */
	virtual void fait_determinant() const ;


    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    /// Save in a file
    
    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    


};



#endif
