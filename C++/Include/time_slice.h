/*
 *  Definition of Lorene class Time_slice
 *
 */

/*
 *   Copyright (c) 2004 Eric Gourgoulhon, Jose Luis Jaramillo & Jerome Novak
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

#ifndef __TIME_SLICE_H_ 
#define __TIME_SLICE_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/03/24 14:56:18  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

class Sym_tensor ; 
class Vector ; 
class Scalar ; 
class Metric ; 
class Base_vect ; 
class Map ; 

#include "headcpp.h"

#include "evolution.h"

/**
 * Spacelike time slice of a 3+1 spacetime (*** under development ***)
 * \ingroup(3p1)
 * 
 */
class Time_slice {

    // Data : 
    // -----
    protected:
        /// Number of stored time slices
        int depth ; 
        
        /// Time step index of the latest slice
        int jtime ; 
        
        /// Time label of each slice
        Evolution_std<double> the_time ;
        
        /** Values of the covariant components of the induced metric 
         * \f$ \gamma_{ij} \f$
         */
        Evolution_std<Sym_tensor> gamma_dd ; 
        
        /** Values of the contravariant components of the induced metric 
         * \f$ \gamma^{ij} \f$
         */        
        Evolution_std<Sym_tensor> gamma_uu ; 

        /** Values of the covariant components of the extrinsic curvature
         * tensor \f$ K_{ij} \f$
         */
        Evolution_std<Sym_tensor> kk_dd ; 

        /** Values of the contravariant components of the extrinsic curvature
         * tensor \f$ K^{ij} \f$
         */
        Evolution_std<Sym_tensor> kk_uu ; 

        /// Values of the lapse function \e N 
        Evolution_std<Scalar> lapse ; 
        
        /// Values of the shift vector \f$ \beta^i \f$
        Evolution_std<Vector> shift ; 
        

    // Derived data : 
    // ------------
    protected:
    /// Pointer on the induced metric 
	mutable Metric* p_gamma ;   

    // Constructors - Destructor
    // -------------------------
    public:
    
    /** General constructor (Hamiltonian-like). 
     *
     *  @param lapse_in lapse function \e N
     *  @param shift_in shift vector
     *  @param gamma_in induced metric (covariant or contravariant components) 
     *  @param kk_in extrinsic curvature (covariant or contravariant components)
     *  @param depth_in  number of stored time slices
     */
    Time_slice(const Scalar& lapse_in, const Vector& shift_in,
               const Sym_tensor& gamma_in, const Sym_tensor kk_in,
               int depth_in = 3) ; 
    
    /** General constructor (Lagrangian-like). 
     *
     *  @param lapse_in lapse function \e N
     *  @param shift_in shift vector
     *  @param gamma_in induced metric (covariant or contravariant components) 
     *          at various time steps
     */
    Time_slice(const Scalar& lapse_in, const Vector& shift_in,
               const Evolution_std<Sym_tensor>& gamma_in) ; 
    
    /** Constructor as standard time slice of flat spacetime (Minkowski). 
     *
     *  @param mp Mapping on which the various Lorene fields will be constructed
     *  @param triad vector basis with respect to which the various tensor
     *      components will be defined
     *  @param depth_in  number of stored time slices
     */
    Time_slice(const Map& mp, const Base_vect& triad, int depth_in = 3) ; 
    
    
	Time_slice(const Time_slice& ) ;		///< Copy constructor

	virtual ~Time_slice() ;			///< Destructor
 

    // Memory management
    // -----------------
    protected:
	    
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
	virtual void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Time_slice
	void operator=(const Time_slice&) ;	
	
    // Accessors
    // ---------
    public:

    // Outputs
    // -------
    public:

	/// Display
	friend ostream& operator<<(ostream& , const Time_slice& ) ;	


};

ostream& operator<<(ostream& , const Time_slice& ) ;	

#endif
