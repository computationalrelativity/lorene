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
 * Revision 1.18  2003/11/03 22:31:02  e_gourgoulhon
 * Class Vector_divfree: parameters of methods set_vr_eta_mu and set_vr_mu
 * are now all references.
 *
 * Revision 1.17  2003/11/03 14:02:17  e_gourgoulhon
 * Class Vector_divfree: the members p_eta and p_mu are no longer declared
 * "const".
 *
 * Revision 1.16  2003/10/20 15:12:17  j_novak
 * New method Vector::poisson
 *
 * Revision 1.15  2003/10/20 14:44:49  e_gourgoulhon
 * Class Vector_divfree: added method poisson().
 *
 * Revision 1.14  2003/10/20 09:32:10  j_novak
 * Members p_potential and p_div_free of the Helmholtz decomposition
 * + the method decompose_div(Metric).
 *
 * Revision 1.13  2003/10/17 16:36:53  e_gourgoulhon
 * Class Vector_divfree: Added new methods set_vr_eta_mu and set_vr_mu.
 *
 * Revision 1.12  2003/10/16 21:36:01  e_gourgoulhon
 * Added method Vector_divfree::update_vtvp().
 *
 * Revision 1.11  2003/10/16 15:25:00  e_gourgoulhon
 * Changes in documentation.
 *
 * Revision 1.10  2003/10/16 14:21:33  j_novak
 * The calculation of the divergence of a Tensor is now possible.
 *
 * Revision 1.9  2003/10/13 20:49:26  e_gourgoulhon
 * Corrected typo in the comments.
 *
 * Revision 1.8  2003/10/13 13:52:39  j_novak
 * Better managment of derived quantities.
 *
 * Revision 1.7  2003/10/12 20:33:15  e_gourgoulhon
 * Added new derived class Vector_divfree (preliminary).
 *
 * Revision 1.6  2003/10/06 13:58:45  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.5  2003/10/05 21:07:27  e_gourgoulhon
 * Method std_spectral_base() is now virtual.
 *
 * Revision 1.4  2003/10/03 14:08:03  e_gourgoulhon
 * Added constructor from Tensor.
 *
 * Revision 1.3  2003/09/29 13:48:17  j_novak
 * New class Delta.
 *
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

class Vector_divfree ;

/**
 * Tensor field of valence 1.
 * 
 * @version #$Id$#
 */
class Vector: public Tensor {

    // Derived data : 
    // ------------
    protected:
        /** The potential $\phi$giving the gradient part in the Helmholtz 
	 * decomposition of any 3D vector $\vec{V}: \quad \vec{V} = 
	 * \vec{\nabla} \phi + \vec{\nabla} \wedge \vec{\psi}$.
	 * Only in the case of contravariant vectors.
	 */
        mutable Scalar* p_potential[N_MET_MAX] ;

	/** The divergence-free vector $\vec{W} =  \vec{\nabla} \wedge 
	 * \vec{\psi}$ of the Helmholtz decomposition of any 3D vector 
	 *$\vec{V}: \quad \vec{V} = \vec{\nabla} \phi + \vec{\nabla} 
	 *\wedge \vec{\psi}$. Only in the case of contravariant vectors.
	 */
	mutable Vector_divfree* p_div_free[N_MET_MAX] ;

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

	/** Constructor from a {\tt Tensor}.
	 *  The {\tt Tensor} must be of valence one.
	 */
	Vector(const Tensor& ) ;	

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
	Vector(const Map& map, const Base_vect& triad_i, FILE* fich) ;

	virtual ~Vector() ;			/// Destructor

 
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const;	/// Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 

	/**
	 * Logical destructor of the derivatives depending on the i-th
	 * element of {\tt met\_depend} in the class {\tt vector}.
	 */	
	virtual void del_derive_met(int) const ;

	/**
	 * Sets all the i-th components of {\tt met\_depend} in the 
	 * class {\tt Vector} ({\tt p\_potential}, etc...) to 0x0.
	 */
	void set_der_met_0x0(int) const ;


    // Mutators / assignment
    // ---------------------
    public:
	/** Sets a new vectorial basis (triad) of decomposition and modifies
	 *  the components accordingly. 
	 */
	virtual void change_triad(const Base_vect& ) ; 

	/// Assignment from a Vector
	virtual void operator=(const Vector&) ;	
	
	/// Assignment from a Tensor
	virtual void operator=(const Tensor&) ;	

	/**Makes the Helmholtz decomposition (see documentation of
	 * {\tt p_potential}) of {\tt this} with respect to a given
	 * {\tt Metric}, only in the case of contravariant vectors.
	 */
	void decompose_div(const Metric&) const ;

	/**Returns the potential in the Helmholtz decomposition.
	 *
	 * It first makes the Helmholtz decomposition (see documentation of
	 * {\tt p_potential}) of {\tt this} with respect to a given
	 * {\tt Metric} and then returns $\phi$. Only in the case
	 * of contravariant vectors.
	 */
	const Scalar& potential(const Metric& ) const ;
	
	/**Returns the div-free vector in the Helmholtz decomposition.
	 *
	 * It first makes the Helmholtz decomposition (see documentation of
	 * {\tt p_potential}) of {\tt this} with respect to a given
	 * {\tt Metric} and then returns $\vec{W}$. Only in the case
	 * of contravariant vectors.	
	 */
	const Vector_divfree& div_free(const Metric& ) const;
	
    // Accessors
    // ---------
    public:
	Scalar& set(int ) ; /// Read/write access to a component

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
	
	/**
	 * Sets the standard spectal bases of decomposition for each component.
	 *
	 */
	virtual void std_spectral_base() ; 

    // Differential operators/ PDE solvers
    // -----------------------------------
    public:
	/**The divergence of {\tt this} with respect to a {\tt Metric}.
	 * The {\tt Vector} is assumed to be contravariant.
	 */
	const Scalar& divergence(const Metric&) const ; 

	/**Solves the vector Poisson equation with {\tt *this} as a source.
	 * 
	 * The equatiopn solved is $\Delta N^i +\lambda \nabla^i 
	 * \nabla_k N^k = S^i$.
	 * {\tt *this} must be given with {\tt dzpuis} = 4.
	 * It uses the Helmholtz decomposition (see documentation of
	 * {\tt p_potential}), with a flat metric, deduced from the triad.
	 *
	 * @param lambda [input] $\lambda$.
	 *
	 * @return the solution $N^i$.
	 */
	Vector poisson(const double lambda) const ;
     
	/**Solves the vector Poisson equation with {\tt *this} as a source
	 * and parameters controlling the solution.
	 * 
	 * The equatiopn solved is $\Delta N^i +\lambda \nabla^i 
	 * \nabla_k N^k = S^i$.
	 * {\tt *this} must be given with {\tt dzpuis} = 4.
	 * It uses the Helmholtz decomposition (see documentation of
	 * {\tt p_potential}), with a flat metric, deduced from the triad.
	 *
	 * @param lambda [input] $\lambda$.
	 *   @param par [input/output] possible parameters
	 *   @param uu [input/output] solution {\it u} with the 
	 *              boundary condition {\it u}=0 at spatial infinity. 
	 */
	void poisson(const double lambda, Param& par, Scalar& uu ) const ;
     
 
};


/**
 * Divergence-free vectors.
 *
 * This class is designed to store divergence-free vectors,
 * with the component expressed in a orthonormal spherical basis
 * $(e_r,e_\theta,e_\varphi)$.
 *
 * 
 * @version #$Id$#
 */
class Vector_divfree: public Vector {

    // Data : 
    // -----
    protected:
	/// Metric with respect to which the divergence is defined
	const Metric* const met_div ; 
	
	/** Field $\eta$ such that the angular components $(V^r,V^\theta)$
	 * of the vector are written:
	 * \begin{equation}
	 *	V^\theta =   {\partial \eta \over \partial\theta}
	 *		- {1\over\sin\theta} {\partial \mu \over \partial\varphi} 
	 * \end{equation} 
	 * \begin{equation}
	 *	V^\varphi =  {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} 
	 * \end{equation} 
	 */
	mutable Scalar* p_eta ;
	
	/** Field $\mu$ such that the angular components $(V^r,V^\theta)$
	 * of the vector are written:
	 * \begin{equation}
	 *	V^\theta =  {\partial \eta \over \partial\theta}
	 *		- {1\over\sin\theta} {\partial \mu \over \partial\varphi} 
	 * \end{equation} 
	 * \begin{equation}
	 *	V^\varphi =  {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} 
	 * \end{equation} 
	 */
	mutable Scalar* p_mu ;
	
    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *		    the vector components are defined 
	 * @param met the metric with respect to which the divergence is defined
	 */
	Vector_divfree(const Map& map, const Base_vect& triad_i, 
		const Metric& met) ;

	Vector_divfree(const Vector_divfree& ) ;       /// Copy constructor

	/** Constructor from a file (see {\tt Tensor::sauve(FILE* )}).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param met the metric with respect to which the divergence is defined
	 * @param fich  file which has been created by 
	 *			    the function {\tt sauve(FILE* )}.
	 */
	Vector_divfree(const Map& map, const Base_vect& triad_i, 
		const Metric& met, FILE* fich) ;

	virtual ~Vector_divfree() ;			/// Destructor

 
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const;	/// Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------

	public:
	/// Assignment from another {\tt Vector\_divfree}
	void operator=(const Vector_divfree&) ;	
	
	/// Assignment from a {\tt Vector}
	virtual void operator=(const Vector&) ;	
	
	/// Assignment from a {\tt Tensor}
	virtual void operator=(const Tensor&) ;	
	
	/** Sets the angular potentials $\eta$ and $\mu$ (see members
	 *  {\tt p\_eta} and {\tt p\_mu}), as well as the $V^r$ component
	 *  of the vector. 
	 *  The components $V^\theta$ and $V^\varphi$ are updated consistently
	 *  by a call to the method {\tt update\_vtvp()}.
	 *
	 *	@param vr_i [input] component $V^r$ of the vector
	 *	@param eta_i [input] angular potential $\eta$
	 *	@param mu_i [input] angular potential $\mu$
	 *
	 */
	void set_vr_eta_mu(const Scalar& vr_i, const Scalar& eta_i,
		const Scalar& mu_i) ; 

	/** Sets the angular potentials $\mu$ (see member
	 *  {\tt p\_mu}), and the $V^r$ component
	 *  of the vector. The potential $\eta$ is then deduced from
	 *  $V^r$ by the divergence-free condition. 
	 *  The components $V^\theta$ and $V^\varphi$ are updated consistently
	 *  by a call to the method {\tt update\_vtvp()}.
	 *
	 *	@param vr_i [input] component $V^r$ of the vector
	 *	@param mu_i [input] angular potential $\mu$
	 *
	 */
	void set_vr_mu(const Scalar& vr_i, const Scalar& mu_i) ; 
	

	// Computational methods
	// ---------------------
	public:
	/** Gives the field $\eta$ such that the angular components $(V^r,V^\theta)$
	 * of the vector are written:
	 * \begin{equation}
	 *	V^\theta =  {\partial \eta \over \partial\theta}
	 *		- {1\over\sin\theta} {\partial \mu \over \partial\varphi} 
	 * \end{equation} 
	 * \begin{equation}
	 *	V^\varphi =  {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} 
	 * \end{equation} 
	 */
	const Scalar& eta() const ;
	
	/** Gives the field $\mu$ such that the angular components $(V^r,V^\theta)$
	 * of the vector are written:
	 * \begin{equation}
	 *	V^\theta =  {\partial \eta \over \partial\theta}
	 *		- {1\over\sin\theta} {\partial \mu \over \partial\varphi}
	 * \end{equation} 
	 * \begin{equation}
	 *	V^\varphi =  {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} 
	 * \end{equation} 
	 */
	const Scalar& mu() const ;
	
	/** Computes the components $V^\theta$ and $V^\varphi$ from the
	 *  potential $\eta$ and  $\mu$, according to:
	 * \begin{equation}
	 *	V^\theta =  {\partial \eta \over \partial\theta}
	 *		- {1\over\sin\theta} {\partial \mu \over \partial\varphi}
	 * \end{equation} 
	 * \begin{equation}
	 *	V^\varphi =  {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} 
	 * \end{equation} 
	 */
	void update_vtvp() ;
	
	/** Computes the solution of a vectorial Poisson equation
	 *  with {\tt *this} $= \vec{V}$ as a source:
	 * \begin{equation}
	 *    \Delta \vec{W} = \vec{V}
	 * \end{equation} 
	 * 
	 * @return solution $\vec{W}$ of the above equation with the boundary
	 *	condition $\vec{W}=0$ at spatial infinity.
	 */
	Vector_divfree poisson() const ; 
	
	
};





#endif
