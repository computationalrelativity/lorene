/*
 *  Definition of Lorene class Compobj
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
 * \]
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

	/// Total energy density \f$E\f$ in the Eulerian frame 
	Scalar ener_euler ; 

	/// Total 3-momentum density \f$P^\alpha\f$ in the Eulerian frame 
	Vector mom_euler ; 

	/// Stress tensor \f$S_{\alpha\beta}\f$  with respect to the Eulerian observer
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
	 * @param mp_i Mapping on which the star will be defined
	 * @param fich	input file (must have been created by the function
	 *	\c sauve)
	 */
	/// Constructor from a file (see \c sauve(FILE*) )
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

	/// Returns the total energy density \f$E\f$ in the Eulerian frame 
	const Scalar& get_ener_euler() const {return ener_euler;}  ; 

	/// Returns the total 3-momentum density \f$P^\alpha\f$ in the Eulerian frame 
	const Vector& get_mom_euler() const {return mom_euler;} ; 

	/// Returns the stress tensor \f$S_{\alpha\beta}\f$  with respect to the Eulerian observer
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
	/// Gravitational mass
    	virtual double mass_g() const ;
};

#endif
