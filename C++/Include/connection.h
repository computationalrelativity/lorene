/*
 *  Definition of Lorene class Connection
 *
 */

/*
 *   Copyright (c) 2003	Eric Gourgoulhon & Jerome Novak
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

#ifndef __CONNECTION_H_ 
#define __CONNECTION_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.8  2003/10/06 13:58:45  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.7  2003/10/06 06:52:26  e_gourgoulhon
 * Corrected documentation.
 *
 * Revision 1.6  2003/10/05 21:04:25  e_gourgoulhon
 * Improved comments
 *
 * Revision 1.5  2003/10/03 14:07:23  e_gourgoulhon
 * Added derived class Connection_fcart.
 *
 * Revision 1.4  2003/10/02 21:31:11  e_gourgoulhon
 * Added methods fait_delta and update
 * flat_conn is now a modifiable pointer.
 *
 * Revision 1.3  2003/10/02 15:44:23  j_novak
 * The destructor is now public...
 *
 * Revision 1.2  2003/10/01 15:41:49  e_gourgoulhon
 * Added mapping
 *
 * Revision 1.1  2003/09/29 21:14:10  e_gourgoulhon
 * First version --- not ready yet.
 *
 *
 *
 * $Header$
 *
 */


// Lorene headers
#include "tensor.h"

class Metric ; 
class Connection_flat ; 

				//--------------------------//
				//     class Connection     // 
				//--------------------------//

/**
 * Class Connection.
 *
 * The class deals only with torsion-free connections. 
 * 
 * Note that we use the MTW convention for the indices of the connection 
 * coefficients with respect to a given triad $(e_i)$:
 * \begin{equation}
 *  \Gamma^i_{\ jk} := \langle e^i, \nabla_{e_k} \, e_j \rangle 
 * \end{equation} 
 *  
 * @version #$Id$#
 */
class Connection {

    // Data : 
    // -----
    protected:

	const Map* const mp ;	/// Reference mapping

	/** Triad $(e_i)$ with respect to which the connection coefficients 
	 * are defined.
	 */
	const Base_vect* const triad ; 
	
	/** Tensor $\Delta^i_{\ jk}$ which defines
	 *  the connection with respect to the flat one: $\Delta^i_{\ jk}$ 
	 * is the difference between the connection coefficients 
	 *  $\Gamma^i_{\ jk}$ and
	 * the connection coefficients ${\bar \Gamma}^i_{\ jk}$ of the
	 * flat connection. The connection coefficients with respect to the 
	 * triad $(e_i)$ are defined
	 * according to the MTW convention:
	 * \begin{equation}
	 *  \Gamma^i_{\ jk} := \langle e^i, \nabla_{e_k} \, e_j \rangle
	 * \end{equation} 
	 *
	 */
	Tensor_delta delta ; 

	/** Indicates whether the connection is associated with a metric
	 *  (in which case the Ricci tensor is symmetric, i.e. the
	 *	  actual type of {\tt p\_ricci} is a {\tt Sym\_tensor})
	 */
	bool assoc_metric ;


	private:

	/** Flat connection with respect to which $\Delta^i_{\ jk}$ 
	 *   (member {\tt delta}) is defined. 
	 *
	 */
	const Connection_flat* flat_conn ;


    // Derived data : 
    // ------------
    protected:

	/// Pointer of the Ricci tensor associated with the connection 
	mutable Tensor* p_ricci ;   

    // Constructors - Destructor
    // -------------------------
    public:
	
	/** Standard constructor ab initio.
	 *
	 * @param delta_i tensor $\Delta^i_{\ jk}$ which defines
	 *  the connection with respect to the flat one: $\Delta^i_{\ jk}$ 
	 * is the difference between the connection coefficients 
	 *  $\Gamma^i_{\ jk}$ and
	 * the connection coefficients ${\bar \Gamma}^i_{\ jk}$ of the
	 * flat connection. The connection coefficients with respect to the 
	 * triad $(e_i)$ are defined according to the MTW convention:
	 * \begin{equation}
	 *  \Gamma^i_{\ jk} := \langle e^i, \nabla_{e_k} \, e_j \rangle
	 * \end{equation} 
	 *
	 */
	explicit Connection(const Tensor_delta& delta_i) ;		
	
	/** Standard constructor for a connection associated with a metric. 
	 *
	 * @param met  Metric to which the connection will be associated
	 *
	 */
	explicit Connection(const Metric& met) ;		
	
	Connection(const Connection& ) ;		/// Copy constructor
	
	protected:

	/// Constructor for derived classes
	Connection(const Map&, const Base_vect& ) ; 		

        public:
	virtual ~Connection() ;			/// Destructor
 

    // Memory management
    // -----------------
    protected:	    

	/// Deletes all the derived quantities
	void del_deriv() const ; 
	
	/// Sets to {\tt 0x0} all the pointers on derived quantities
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:

	/// Assignment to another {\tt Connection}
	void operator=(const Connection&) ;	
	
    // Accessors
    // ---------
    public:

	const Map& get_mp() const {return *mp; } ;  /// Returns the mapping


	/** Returns the tensor $\Delta^i_{\ jk}$ which defines
	 *  the connection with respect to the flat one: $\Delta^i_{\ jk}$ 
	 * is the difference between the connection coefficients 
	 *  $\Gamma^i_{\ jk}$ and
	 * the connection coefficients ${\bar \Gamma}^i_{\ jk}$ of the
	 * flat connection. The connection coefficients with respect to the 
	 * triad $(e_i)$ are defined according to the MTW convention:
	 * \begin{equation}
	 *  \Gamma^i_{\ jk} := \langle e^i, \nabla_{e_k} \, e_j \rangle
	 * \end{equation} 
	 *
	 * @return {\tt delta}(i,j,k) = $\Delta^i_{\ jk}$
	 */
	const Tensor_delta& get_delta() const {return delta; } ; 

	// Computational methods
	// ---------------------
	
	public: 

	/// Covariant derivative of a tensor (with respect to the current connection)
	virtual Tensor derive_cov(const Tensor&) const ; 

	/** Pointer on the covariant derivative of a tensor 
	 * (with respect to the current connection). This method allocates memory
	 * that must be deallocated by the user afterwards.
	 */
	virtual Tensor* p_derive_cov(const Tensor&) const ; 

	/// Returns the Ricci tensor associated with the current connection
	const Tensor& ricci() const ; 
	
	/** Update the connection when it is defined ab initio.
	 *
	 * @param delta_i tensor $\Delta^i_{\ jk}$ which defines
	 *  the connection with respect to the flat one: $\Delta^i_{\ jk}$ 
	 * is the difference between the connection coefficients 
	 *  $\Gamma^i_{\ jk}$ and
	 * the connection coefficients ${\bar \Gamma}^i_{\ jk}$ of the
	 * flat connection. 
	 */
	void update(const Tensor_delta& delta_i) ;		
	
	/** Update the connection when it is associated with a metric. 
	 *
	 * @param met  Metric to which the connection is associated
	 *
	 */
	void update(const Metric& met) ;
	
			
	protected:

	/// Computes the Ricci tensor when necessary
	virtual void compute_ricci() const ; 

	private:
	/** Computes the difference $\Delta^i_{\ jk}$ between the
	 *  connection coefficients and that a the flat connection
	 *  in the case where the current connection is associated
	 *  with a metric
	 */
	void fait_delta(const Metric& ) ;  

};


				//-------------------------------//
				//     class Connection_flat     // 
				//-------------------------------//

/**
 * Class Connection\_flat.
 *
 * Abstract class for connections associated with a flat metric. 
 * 
 * @version #$Id$#
 */
class Connection_flat : public Connection {

  // Constructors - Destructor
  // -------------------------
 protected:
  
  /// Contructor from a triad, has to be defined in the derived classes
  Connection_flat(const Map&, const Base_vect&) ; 

 public:

	Connection_flat(const Connection_flat & ) ;		/// Copy constructor

  virtual ~Connection_flat() ; /// destructor


  // Mutators / assignment
  // ---------------------
 public:

  /// Assignment to another {\tt Connection\_flat}
  void operator=(const Connection_flat&) ;	
  

  // Computational methods
  // ---------------------
  
 public: 

  /// Covariant derivative of a tensor (with respect to the current connection)
  virtual Tensor derive_cov(const Tensor&) const = 0 ; 

  /** Pointer on the covariant derivative of a tensor 
   * (with respect to the current connection). This method allocates memory
   * that must be deallocated by the user afterwards.
   */
  virtual Tensor* p_derive_cov(const Tensor&) const = 0 ; 

 protected:

  /// Computes the Ricci tensor when necessary
  virtual void compute_ricci() const ; 
  
};


				//-------------------------------//
				//     class Connection_fspher   // 
				//-------------------------------//

/**
 * Class Connection\_fspher.
 *
 * Class for connections associated with a flat metric and given onto
 * an orthonormal spherical triad. 
 * 
 * @version #$Id$#
 */
class Connection_fspher : public Connection_flat {

  // Constructors - Destructor
  // -------------------------
  
 public:

  /// Contructor from a spherical flat-metric-orthonormal basis
	Connection_fspher(const Map&, const Base_vect_spher&) ; 

	Connection_fspher(const Connection_fspher& ) ;		/// Copy constructor
	
 public:

  virtual ~Connection_fspher() ; ///destructor


  // Mutators / assignment
  // ---------------------
 public:

  /// Assignment to another {\tt Connection\_fspher}
  void operator=(const Connection_fspher&) ;	
  

  // Computational methods
  // ---------------------
  
 public: 
  /// Covariant derivative of a tensor (with respect to the current connection)
  virtual Tensor derive_cov(const Tensor&) const ; 

  /** Pointer on the covariant derivative of a tensor 
   * (with respect to the current connection). This method allocates memory
   * that must be deallocated by the user afterwards.
   */
  virtual Tensor* p_derive_cov(const Tensor&) const ; 

  
};



				//-------------------------------//
				//     class Connection_fcart    // 
				//-------------------------------//

/**
 * Class Connection\_fcart.
 *
 * Class for connections associated with a flat metric and given onto
 * an orthonormal Cartesian triad. 
 * 
 * @version #$Id$#
 */
class Connection_fcart : public Connection_flat {

  // Constructors - Destructor
  // -------------------------
  
 public:

  /// Contructor from a Cartesian flat-metric-orthonormal basis
	Connection_fcart(const Map&, const Base_vect_cart&) ; 

	Connection_fcart(const Connection_fcart& ) ;		/// Copy constructor
	
 public:

  virtual ~Connection_fcart() ; ///destructor


  // Mutators / assignment
  // ---------------------
 public:

  /// Assignment to another {\tt Connection\_fcart}
  void operator=(const Connection_fcart&) ;	
  

  // Computational methods
  // ---------------------
  
 public: 
  /// Covariant derivative of a tensor (with respect to the current connection)
  virtual Tensor derive_cov(const Tensor&) const ; 

  /** Pointer on the covariant derivative of a tensor 
   * (with respect to the current connection). This method allocates memory
   * that must be deallocated by the user afterwards.
   */
  virtual Tensor* p_derive_cov(const Tensor&) const  ; 

  
};




#endif
