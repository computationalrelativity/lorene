/*
 *  Definition of Lorene class Compobj, Compobj_QI
 *
 */

/*
 *   Copyright (c) 2012 Claire Some, Eric Gourgoulhon
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

#ifndef __COMPOBJ_H_ 
#define __COMPOBJ_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.3  2012/11/16 16:13:12  c_some
 * Added new class Compobj_QI
 *
 * Revision 1.2  2012/11/15 20:50:41  e_gourgoulhon
 * Corrected the documentation
 *
 * Revision 1.1  2012/11/15 16:20:51  c_some
 * New class Compobj
 *
 *
 * $Header$
 *
 */


// Headers Lorene
#include "tensor.h"
#include "metric.h"


			//---------------------------//
			//    base class Compobj     //
			//---------------------------//

/**
 * Base class for stationary compact objects (***under development***). 
 * \ingroup(compactobjects)
 *
 * A \c Compobj describes a single compact object (star or black hole), in a stationary state.
 * 
 * The spacetime metric is written according to the 3+1 formalism :
 * \f[
 *   ds^2 = - N^2  dt^2 + \gamma_{ij} ( dx^i + \beta^i dt )
 *               (dx^j + \beta^j dt )
 * \f]
 * where \f$\gamma_{ij}\f$ is the 3-metric, described by a Lorene object of class \c Metric. 
 * 
 * The total energy-momentum tensor is orthogonally split with respect to the Eulerian observer as follows:
 * \f[
 *	T_{\alpha\beta} = E n_\alpha n_\beta + P_\alpha n_\beta + n_\alpha P_\beta + S_{\alpha\beta}
 * \f]
 */
class Compobj {

    // Data : 
    // -----
    protected:
	/// Mapping describing the coordinate system (r,theta,phi) 
	Map& mp ;  

	/// Lapse function \e N .
	Scalar nn ; 
	
	/// Shift vector \f$\beta^i\f$
	Vector beta ;
	
 	/// 3-metric  \f$\gamma_{ij}\f$
	Metric gamma ;

	/// Total energy density \e E in the Eulerian frame 
	Scalar ener_euler ; 

	/// Total 3-momentum density \f$P^i\f$ in the Eulerian frame 
	Vector mom_euler ; 

	/// Stress tensor \f$S_{ij}\f$  with respect to the Eulerian observer
	Sym_tensor stress_euler ;


   // Derived data : 
    // ------------
    protected:
	mutable double* p_mass_g ;	///< Gravitational mass (ADM mass) 

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the object is defined
	 * 
	 */
	Compobj(Map& map_i) ;

	Compobj(const Compobj& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE* )). 
	 * 
	 * @param mp_i Mapping on which the object is defined
	 * @param fich	input file (must have been created by the function
	 *	\c sauve)
	 */
	Compobj(Map& map_i, FILE* ) ;    		

	virtual ~Compobj() ;			///< Destructor
 

    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Compobj
	void operator=(const Compobj&) ;	
	
	/// Read/write of the mapping
	Map& set_mp() {return mp; } ; 


    // Accessors
    // ---------
    public:
	/// Returns the mapping
	const Map& get_mp() const {return mp; } ; 

	/// Returns the lapse function \e N .
	const Scalar& get_nn() const {return nn;} ;

	/// Returns the shift vector \f$\beta^i\f$.
	const Vector& get_beta() const {return beta;} ;
	
 	/// Returns the 3-metric \f$\gamma_{ij}\f$.
	const Metric& get_gamma() const {return gamma;} ;

	/// Returns the total energy density \e E in the Eulerian frame 
	const Scalar& get_ener_euler() const {return ener_euler;}  ; 

	/// Returns the total 3-momentum density \f$P^i\f$ in the Eulerian frame 
	const Vector& get_mom_euler() const {return mom_euler;} ; 

	/// Returns the stress tensor \f$S_{ij}\f$  with respect to the Eulerian observer
	const Sym_tensor& get_stress_euler() const {return stress_euler;} ;

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Compobj& ) ;	

    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

    // Global quantities
    // -----------------
    public:
	/// Gravitational mass (ADM mass)
    	virtual double mass_g() const ;
};


			//---------------------------//
			//    base class Compobj_QI  //
			//---------------------------//

/**
 * Base class for axisymmetric stationary compact objects in Quasi-Isotropic coordinates (***under development***). 
 * \ingroup(compactobjects)
 *
 * The metric is expressed in Quasi-Isotropic (QI) coordinates :
 * \f[
 *   ds^2 = - N^2 dt^2 + A^2 (dr^2 + r^2 d\theta^2)
 *		       + B^2 r^2 \sin^2\theta (d\varphi - N^\varphi dt)^2
 * \f]
 *
 * 
 */
class Compobj_QI : public Compobj {

    // Data : 
    // -----
    protected:

	/// Square of the metric factor \e A 
	Scalar a_car ; 

	/// Metric factor \e B 
	Scalar bbb ; 

	/// Square of the metric factor \e B 
	Scalar b_car ; 

	/// Metric coefficient \f$N^\varphi\f$
	Scalar nphi ; 

	/** Tensor \f${\tilde K_{ij}}\f$ related to the extrinsic curvature
	 *  tensor by \f${\tilde K_{ij}} = B^{-2} K_{ij}\f$.
	 *  \c tkij  contains the Cartesian components of
	 *  \f${\tilde K_{ij}}\f$. 
	 */
	Sym_tensor tkij ; 

	/** Scalar \f$A^2 K_{ij} K^{ij}\f$.
	 *  For axisymmetric stars, this quantity is related to the 
	 *  derivatives of \f$N^\varphi\f$ by
	 * \f[
	 *	A^2 K_{ij} K^{ij} = {B^2 \over 2 N^2} \, r^2\sin^2\theta \,  
	 *    \left[ \left( {\partial N^\varphi \over \partial r} \right) ^2
	 *	    + {1\over r^2} \left( {\partial N^\varphi \over 
	 *		    \partial \theta} \right) ^2 \right] \ . 
	 * \f]
	 * In particular it is related to the quantities \f$k_1\f$ and \f$k_2\f$
	 * introduced by Eqs.~(3.7) and (3.8) of 
	 * Bonazzola et al. \a Astron. \a Astrophys. \b 278 , 421 (1993)
	 * by 
	 * \f[
	 *	A^2 K_{ij} K^{ij} = 2 A^2 (k_1^2 + k_2^2) \ . 
	 * \f]
	 */
	 Scalar ak_car ; 


   // Derived data : 
    // ------------
    protected:
	mutable double* p_angu_mom ;	///< Angular momentum 
	mutable double* p_r_isco ;	///< Circumferential radius of the ISCO
	mutable double* p_f_isco ;	///< Orbital frequency of the ISCO
	/// Specific energy of a particle on the ISCO 
	mutable double* p_espec_isco ;	
	/// Specific angular momentum of a particle on the ISCO
	mutable double* p_lspec_isco ;	

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the object is defined
	 * 
	 */
	Compobj_QI(Map& map_i) ;

	Compobj_QI(const Compobj_QI& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE* )). 
	 * 
	 * @param mp_i Mapping on which the object is defined
	 * @param fich	input file (must have been created by the function
	 *	\c sauve)
	 */
	Compobj_QI(Map& map_i, FILE* ) ;    		

	virtual ~Compobj_QI() ;			///< Destructor
 

    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Compobj_QI
	void operator=(const Compobj_QI&) ;	
	

    // Accessors
    // ---------
    public:

	/// Returns the metric factor \e B 
	const Scalar& get_bbb() const {return bbb;} ; 

	/// Returns the square of the metric factor \e A 
	const Scalar& get_a_car() const {return a_car;} ; 

	/// Returns the square of the metric factor \e B 
	const Scalar& get_b_car() const {return b_car;} ; 

	/// Returns the metric coefficient \f$N^\varphi\f$
	const Scalar& get_nphi() const {return nphi;} ; 

	/** Returns the tensor \f${\tilde K_{ij}}\f$ related to the extrinsic 
	 *  curvature tensor by \f${\tilde K_{ij}} = B^{-2} K_{ij}\f$.
	 *  \c tkij  contains the Cartesian components of
	 *  \f${\tilde K_{ij}}\f$. 
	 */
	const Sym_tensor& get_tkij() const {return tkij;} ; 

	/** Returns the scalar \f$A^2 K_{ij} K^{ij}\f$.
	 *  For axisymmetric stars, this quantity is related to the 
	 *  derivatives of \f$N^\varphi\f$ by
	 * \f[
	 *	A^2 K_{ij} K^{ij} = {B^2 \over 2 N^2} \, r^2\sin^2\theta \,  
	 *    \left[ \left( {\partial N^\varphi \over \partial r} \right) ^2
	 *	    + {1\over r^2} \left( {\partial N^\varphi \over 
	 *		    \partial \theta} \right) ^2 \right] \ . 
	 * \f]
	 * In particular it is related to the quantities \f$k_1\f$ and \f$k_2\f$
	 * introduced by Eqs. (3.7) and (3.8) of 
	 * Bonazzola et al. \a Astron. \a Astrophys. \b 278 , 421 (1993)
	 * by 
	 * \f[
	 *	A^2 K_{ij} K^{ij} = 2 A^2 (k_1^2 + k_2^2) \ . 
	 * \f]
	 */
	 const Scalar& get_ak_car() const {return ak_car;} ; 

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

    // Global quantities
    // -----------------
    public:
	virtual double angu_mom() const ;	///< Angular momentum 

	/** Circumferential radius of the innermost stable circular orbit (ISCO).	
	 *
	 *  @param lmin index of the domain from which the ISCO is searched outwards ;
	 *  @param ost output stream to give details of the computation;
	 *		if set to 0x0 [default value], no details will be
	 *		given.
	 *  
	 */
	virtual double r_isco(int lmin, ostream* ost = 0x0) const ;	
 	
 	/// Orbital frequency at the innermost stable circular orbit (ISCO).	
 	virtual double f_isco(int lmin) const ;	

	/// Energy of a particle on the ISCO 
 	virtual double espec_isco(int lmin) const ;	
	
	/// Angular momentum of a particle on the ISCO
 	virtual double lspec_isco(int lmin) const ;	

    // Computational routines
    // ----------------------

	/** Updates the 3-metric \f$\gamma_{ij}\f$ from \e A and \e B 
	 *  and the shift vector \$f\beta^i\f$  from \f$N^\phi\f$. 
	 * 
	 */
	virtual void update_metric() ; 
		
	/** Computes \c tkij  and \c ak_car  from 
	 *  \c beta , \c nn  and \c b_car .
	 */
	virtual void extrinsic_curvature() ;
	
};

#endif
