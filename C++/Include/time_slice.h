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
 * Revision 1.3  2004/03/28 21:27:57  e_gourgoulhon
 * Class Time_slice: - renamed the Evolution_std with suffix "_evol".
 *                   - added protected constructor for derived classes
 * Added class Time_slice_conf.
 *
 * Revision 1.2  2004/03/26 13:33:02  j_novak
 * New methods for accessing/updating members (nn(), beta(), gam_uu(), k_uu(), ...)
 *
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
class Metric_flat ; 
class Base_vect ; 
class Map ; 

#include "headcpp.h"

#include "evolution.h"

                    //---------------------------//
                    //      class Time_slice     //
                    //---------------------------//

/**
 * Spacelike time slice of a 3+1 spacetime (*** under development ***)
 * \ingroup (evol)
 * 
 */
class Time_slice {

    // Data : 
    // -----
    protected:
        /// Number of stored time slices
        int depth ; 
        
        /** Order of the finite-differences scheme for 
	 * the computation of time derivatives.
	 *
	 * This order is not constant and can be adjusted \e via
	 * \c set_scheme_order() .
	 */
        int scheme_order ; 
        
        /// Time step index of the latest slice
        int jtime ; 
        
        /// Time label of each slice
	Evolution_std<double> the_time ;
        
        /** Values at successive time steps of the covariant components of 
         * the induced metric \f$ \gamma_{ij} \f$
         */
	mutable Evolution_std<Sym_tensor> gam_dd_evol ; 
        
        /** Values at successive time steps of the contravariant components 
         * of the induced metric \f$ \gamma^{ij} \f$
         */        
	mutable Evolution_std<Sym_tensor> gam_uu_evol ; 

        /** Values at successive time steps of the covariant components 
         * of the extrinsic curvature tensor \f$ K_{ij} \f$
         */
	mutable Evolution_std<Sym_tensor> k_dd_evol ; 

        /** Values at successive time steps of the contravariant components 
         * of the extrinsic curvature tensor \f$ K^{ij} \f$
         */
	mutable Evolution_std<Sym_tensor> k_uu_evol ; 

        /// Values at successive time steps of the lapse function \e N 
	mutable Evolution_std<Scalar> n_evol ; 
        
        /// Values at successive time steps of the shift vector \f$ \beta^i \f$
	mutable Evolution_std<Vector> beta_evol ; 
        

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
     *  @param depth_in  number of stored time slices; this parameter is used
     *                   to set the \c scheme_order member with \c scheme_order
     *                   = \c depth_in - 1. \c scheme_order can be changed 
     *                   afterwards by the method \c set_scheme_order(int).
     */
    Time_slice(const Scalar& lapse_in, const Vector& shift_in,
               const Sym_tensor& gamma_in, const Sym_tensor kk_in,
               int depth_in = 3) ; 
    
    /** General constructor (Lagrangian-like). 
     *
     *  @param lapse_in lapse function \e N
     *  @param shift_in shift vector
     *  @param gamma_in induced metric (covariant or contravariant components) 
     *          at various time steps; note that the \c scheme_order member 
     *          is set to \c gamma_in.size - 1. \c scheme_order can be changed 
     *                   afterwards by the method \c set_scheme_order(int).
     */
    Time_slice(const Scalar& lapse_in, const Vector& shift_in,
               const Evolution_std<Sym_tensor>& gamma_in) ; 
    
    /** Constructor as standard time slice of flat spacetime (Minkowski). 
     *
     *  @param mp Mapping on which the various Lorene fields will be constructed
     *  @param triad vector basis with respect to which the various tensor
     *      components will be defined
     *  @param depth_in  number of stored time slices; this parameter is used
     *                   to set the \c scheme_order member with \c scheme_order
     *                   = \c depth_in - 1. \c scheme_order can be changed 
     *                   afterwards by the method \c set_scheme_order(int).
     */
    Time_slice(const Map& mp, const Base_vect& triad, int depth_in = 3) ; 
    
    
    Time_slice(const Time_slice& ) ;		///< Copy constructor
    
    protected:
    /** Special constructor for derived classes.
     *
     */
    Time_slice(int depth_in) ; 

    public:
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

	/// Sets the order of the finite-differences scheme.
	void set_scheme_order(int ord) { 
	  assert (0<= ord < 4) ;
	  scheme_order = ord ; } ; 
	
    // Accessors
    // ---------
    public:

	/// Gets the order of the finite-differences scheme.
	int get_scheme_order() const { return scheme_order ; } ;
	
	/// Lapse function \e N at the current time step (\c jtime )
	virtual const Scalar& nn() const ;
	
	/// shift vector \f$ \beta^i \f$ at the current time step (\c jtime )
	virtual const Vector& beta() const ;
	
	/// Induced metric \f$ \mathbf{\gamma} \f$ at the current time step (\c jtime )
	const Metric& gam() const ;
	
	/** Induced metric (covariant components \f$ \gamma_{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& gam_dd() const ;
	
	/** Induced metric (contravariant components \f$ \gamma^{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& gam_uu() const ;
	
	/** Extrinsic curvature tensor (covariant components \f$ K_{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& k_dd() const ;
	
	/** Extrinsic curvature tensor (contravariant components \f$ K^{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& k_uu() const ;
	
	
    // Outputs
    // -------
    public:
	
    /// Display
    friend ostream& operator<<(ostream& , const Time_slice& ) ;	


};

ostream& operator<<(ostream& , const Time_slice& ) ;	



                    //---------------------------//
                    //   class Time_slice_conf   //
                    //---------------------------//

/**
 * Spacelike time slice of a 3+1 spacetime with conformal decomposition
 * (*** under development ***)
 * \ingroup (evol)
 * 
 */
class Time_slice_conf : public Time_slice {

    // Data : 
    // -----
    protected: 
    
        /** Pointer on the flat metric \f$ f_{ij} \f$ with respect to
         * which the conformal decomposition is performed
         */
        const Metric_flat& ff ;  

        /** Values at successive time steps of the conformal factor 
         * \f$ \Psi \f$ relating the
         * physical metric \f$ \gamma_{ij} \f$ to the conformal one:
         * \f$ \gamma_{ij} = \Psi^4 \tilde\gamma_{ij} \f$.
         * \f$ \Psi \f$ is defined by
         *  \f[ \Psi := \left( \frac{\det\gamma_{ij}}{\det f_{ij}} 
         *      \right) ^{1/12} \r] 
         */
	mutable Evolution_std<Scalar> psi_evol ; 
        
        /** Values at successive time steps of the factor 
         * \f$ Q := \Psi^2 N \f$.
         */
	mutable Evolution_std<Scalar> qq_evol ; 
        
        
        /** Values at successive time steps of the components \f$ h^{ij} \f$ 
         * of the deviation 
         * of the conformal metric \f$ \tilde\gamma^{ij} \f$ from 
         * the flat metric \f$ f^{ij} \f$: 
         * \f$\tilde\gamma^{ij} = f^{ij} + h^{ij} \f$.
         */        
	mutable Evolution_std<Sym_tensor> hh_evol ; 

        /** Values at successive time steps of the components \f$ A^{ij} \f$
         * of the conformal representation of the traceless part
         * of the extrinsic curvature:
         * \f$ A^{ij} = \Psi^4 \left( K^{ij} - {1\over 3} K \gamma^{ij} 
         *  \right) \f$.
         */        
	mutable Evolution_std<Sym_tensor> aa_evol ; 

        /** Values at successive time steps of the trace \e K of the 
         *  extrinsic curvature
         */        
	mutable Evolution_std<Scalar> trk_evol ; 

        
    // Derived data : 
    // ------------
    protected:
        /// Pointer on the conformal metric \f$ \tilde\gamma_{ij} \f$
	mutable Metric* p_tgamma ; 
        
        /// Pointer on the factor \f$ \Psi^4 \f$
	mutable Scalar* p_psi4 ; 
        

    // Constructors - Destructor
    // -------------------------
    public:
    
    /** Constructor from conformal decomposition.
     *
     *  @param lapse_in lapse function \e N
     *  @param shift_in shift vector
     *  @param ff_in reference flat metric with respect to which the
     *           conformal decomposition is performed
     *  @param psi_in conformal factor \f$\Psi\f$ relating the
     *       physical metric \f$ \gamma_{ij} \f$ to the conformal one:
     *      \f$ \gamma_{ij} = \Psi^4 \tilde\gamma_{ij} \f$
     *  @param hh_in deviation \f$ h^{ij} \f$
     *      of the conformal metric \f$ \tilde\gamma^{ij} \f$ from 
     *      the flat metric \f$ f^{ij} \f$: 
     *      \f$\tilde\gamma^{ij} = f^{ij} + h^{ij} \f$.
     *  @param aa_in conformal representation \f$ A^{ij} \f$
     *      of the traceless part of the extrinsic curvature:
     *      \f$ A^{ij} = \Psi^4 \left( K^{ij} - {1\over 3} K \gamma^{ij}
     *  @param trk_in trace \e K of the extrinsic curvature 
     *  @param depth_in  number of stored time slices; this parameter is used
     *                   to set the \c scheme_order member with \c scheme_order
     *                   = \c depth_in - 1. \c scheme_order can be changed 
     *                   afterwards by the method \c set_scheme_order(int).
     */
    Time_slice_conf(const Scalar& lapse_in, const Vector& shift_in,
            const Metric_flat& ff_in, const Scalar& psi_in, 
            const Sym_tensor& hh_in, const Sym_tensor aa_in, 
            const Scalar& trk_in, int depth_in = 3) ; 
    
    
    /** Constructor from physical metric.
     *  The conformal decomposition is performed by the constructor. 
     *
     *  @param lapse_in lapse function \e N
     *  @param shift_in shift vector
     *  @param gamma_in induced metric (covariant or contravariant components) 
     *  @param kk_in extrinsic curvature (covariant or contravariant components)
     *  @param ff_in reference flat metric with respect to which the
     *           conformal decomposition is performed
     *  @param depth_in  number of stored time slices; this parameter is used
     *                   to set the \c scheme_order member with \c scheme_order
     *                   = \c depth_in - 1. \c scheme_order can be changed 
     *                   afterwards by the method \c set_scheme_order(int).
     */
    Time_slice_conf(const Scalar& lapse_in, const Vector& shift_in,
               const Sym_tensor& gamma_in, const Sym_tensor kk_in,
               const Metric_flat& ff_in, int depth_in = 3) ; 
               
               
    Time_slice_conf(const Time_slice_conf& ) ;	///< Copy constructor

    virtual ~Time_slice_conf() ;			///< Destructor
 

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
	/// Assignment to another \c Time_slice_conf
	void operator=(const Time_slice_conf&) ;

	/// Assignment to a \c Time_slice
	void operator=(const Time_slice&) ;
	
    // Accessors
    // ---------
    public:

        // Virtual functions from base class Time_slice:
        // ---------------------------------------------

	/// Lapse function \e N at the current time step (\c jtime )
	virtual const Scalar& nn() const ;
	
	/** Induced metric (covariant components \f$ \gamma_{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& gam_dd() const ;
	
	/** Induced metric (contravariant components \f$ \gamma^{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& gam_uu() const ;
	
	/** Extrinsic curvature tensor (covariant components \f$ K_{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& k_dd() const ;
	
	/** Extrinsic curvature tensor (contravariant components \f$ K^{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& k_uu() const ;
	
        // Virtual functions from this class:
        // ----------------------------------

        /** Conformal factor \f$ \Psi \f$ relating the
         * physical metric \f$ \gamma_{ij} \f$ to the conformal one:
         * \f$ \gamma_{ij} = \Psi^4 \tilde\gamma_{ij} \f$. 
         * \f$ \Psi \f$ is defined by
         *  \f[ \Psi := \left( \frac{\det\gamma_{ij}}{\det f_{ij}} 
         *      \right) ^{1/12} \r] 
         * Returns the value at the current time step (\c jtime ).
         */
	virtual const Scalar& psi() const ; 
        
        /// Factor \f$ \Psi^4 \f$ at the current time step (\c jtime ).
	const Scalar& psi4() const ; 
        
        /** Factor \f$ Q := \Psi^2 N \f$ at the current time step (\c jtime ).
         */
	virtual const Scalar& qq() const ; 
        
        /** Conformal metric 
         * \f$ \tilde\gamma_{ij} = \Psi^{-4} \gamma_{ij} \f$
         * Returns the value at the current time step (\c jtime ).
         */        
	const Metric& tgam() const ; 

        /** Deviation \f$ h^{ij} \f$ 
         * of the conformal metric \f$ \tilde\gamma^{ij} \f$ from 
         * the flat metric \f$ f^{ij} \f$: 
         * \f$\tilde\gamma^{ij} = f^{ij} + h^{ij} \f$.
         * Returns the value at the current time step (\c jtime ).
         */        
	virtual const Sym_tensor& hh() const ; 

        /** Conformal representation \f$ A^{ij} \f$ of the traceless part
         * of the extrinsic curvature:
         * \f$ A^{ij} = \Psi^4 \left( K^{ij} - {1\over 3} K \gamma^{ij} 
         *  \right) \f$.
         * Returns the value at the current time step (\c jtime ).
         */        
	virtual const Sym_tensor& aa() const ; 

        /** Trace \e K of the extrinsic curvature 
         *  at the current time step (\c jtime )
         */        
	virtual const Scalar& trk() const ; 

        
    // Outputs
    // -------
    public:
	

};
#endif
