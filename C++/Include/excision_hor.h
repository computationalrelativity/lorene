/*
 *  Definition of Lorene class Excision_hor, friend class of Spheroid.
 *
 */

/*
 *   Copyright (c) 2009  Jose-Luis Jaramillo & Nicolas Vasset
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

#ifndef __EXCISIONHOR_H_ 
#define __EXCISIONHOR_H_ 

/*
 * $Header$
 *
 */
#include "metric.h"
#include "spheroid.h"

/**
 * Surface where boundary conditions for quantities in the bulk will be calculated
 * It relies on geometrical properties of the associated Spheroid\ingroup (star)
 * (*** WARNING! under development***)
 * 
 */
class Excision_hor {


    // Data : 
    // -----
 protected:
  /// The associated Spheroid object
  Spheroid sph ;
  
  ///The value of the conformal factor on the 3-slice
  Scalar conf_fact ; 
  
  /** The lapse defined on the 3 slice*/
  Scalar lapse ;
  
  /** The Shift 3-vector on the slice **/
  Vector shift ;
  
  /** The 3-d metric on the slice*/
  Metric gamij ; 
  
  /** The 3-d extrinsic curvature on the slice */
  Sym_tensor Kij ;

  /** The time step for evolution in parabolic drivers */
  double delta_t;

  /** The internal number of timesteps for one iteration. */
  double no_of_steps;

  /** Value of the impulsion-energy tensor on the spheroid */
  Sym_tensor Tij; 
   

    // Derived data : 
    // ------------
    protected:
  mutable Scalar* p_get_BC_conf_fact ; ///< Source of Neumann BC on \f$ \psi \f$, derived from the vanishing expansion
  mutable Scalar* p_get_BC_bmN ; ///< Source of Dirichlet BC for (b-N).
  mutable Scalar* p_get_BC_bpN ; ///< Arbitrary source of Dirichlet BC for (b+N), case 0:  based on a parabolic driver towards a constant value, case 1:from a component of projected Einstein Equations. .
  mutable Vector* p_get_BC_shift ; ///< Source of Dirichlet BC for the shift, issued from BC on bpN and a gauge condition on the tangential shift.


    // Constructors - Destructor
    // -------------------------
    public:

	/** Constructor of an excision surface embedded in a 3-slice (\c Time_slice ) of 3+1 
	 * formalism. 
	 * This is done from the \c Time_slice data.
	 * @param h_in : the location of the surface \e r = h_in (WARNING:must be 
	                    defined on a mono-domain angular grid)
	 * @param gij : the 3-metric on the 3-slice
	 * @param Kij : the extrinsic curvature of the 3-slice 
	 *                 (covariant representation)
	 * @param timestep : time interval associated with the parabolic-driven boundary conditions.
	 * @param int_nos : Number of iterations to be done during timestep.
         * @param Tij : Value of the impulsion-energy tensor on the spheroid 
	 */
	Excision_hor(const Scalar& h_in, const Metric& gij, const Sym_tensor& Kij2, const Scalar& ppsi, const Scalar& nn, const Vector& beta,const Sym_tensor& Tij2, double timestep, int int_nos = 1) ;
	Excision_hor(const Excision_hor& ) ;		///< Copy constructor

	/// Constructor from a file (see \c sauve(FILE*) )
	Excision_hor(FILE* ) ;    		

	virtual ~Excision_hor() ;			///< Destructor
 

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
	/// Assignment to another Excision_hor
	void operator=(const Excision_hor&) ;	
	

    // Accessors
    // ---------
    public:
	/// Returns the spheroid
	const Spheroid& get_sph() const {return sph; }; 
	
	/// Returns the conformal factor associated with the surface
	const Scalar& get_conf_fact() const {return conf_fact; } ;

	/// Returns the lapse function
	const Scalar& get_lapse() const {return lapse ; } ;

	/// Returns the shift vector field.
	const Vector& get_shift() const {return shift ; } ;

	/// Returns the symmetric tensor \f$ gamij \f$
	const Metric& get_gamij() const {return gamij ; } ;

	/// returns the 3-d extrinsic curvature \f$ K_{ij}\f$
	const Sym_tensor& get_Kij() const {return Kij ; } ;

	/// Returns the timestep used for evolution.
	double get_delta_t() const {return delta_t ;};

	/// Returns the internal number of timesteps for one iteration.
	double get_no_of_steps() const {return no_of_steps ;};
 
        /// Returns the value of the impulsion-energy tensor 
        const Sym_tensor& get_Tij() const{return Tij; } ;

	/// Sets the value of the conformal factor
	Scalar& set_conf_fact() {del_deriv() ; return conf_fact ; } ;

	/// Sets the lapse function
	Scalar& set_lapse() {del_deriv() ; return lapse ; } ;

	/// Sets the shift vector field
	Vector& set_shift() {del_deriv() ; return shift ; } ;

	/// Sets the 3d metric of the TimeSlice
	Metric& set_gamij() {del_deriv() ; return gamij ; } ;

	/// Sets the extrinsic curvature
	Sym_tensor& set_Kij() {del_deriv() ; return Kij ; } ;
 
       /// Sets the value of the impulsion-energy tensor 
        Sym_tensor& set_Tij() {del_deriv() ; return Tij; } ;
	

	double& set_delta_t() {del_deriv() ; return delta_t ; } ;

	double& set_no_of_steps() {del_deriv() ; return no_of_steps ; } ;

	
    // Computational functions
    // -----------------------
    public:

	///  Source of Neumann BC on \f$ \psi \f$, derived from the vanishing expansion.
	const Scalar& get_BC_conf_fact() const ;
        /// Source of Dirichlet BC for (b-N): case 0: based on an entropy prescription, case 1: from a component of projected Einstein Equations.
	const Scalar& get_BC_bmN(int choice_bmN, double value = 1.) const ;
	///Case 0: Arbitrary source of Dirichlet BC for (b+N), based on a parabolic driver towards a constant value.
	///Case 1: Source of Dirichlet BC for (b+N), from a component of projected Einstein Equations.
	const Scalar& get_BC_bpN(int choice_bpN, double c_bpn_lap = 1., double c_bpn_fin = 1., Scalar *bpN_fin = 0x0) const ;
	///  Source of Dirichlet BC for the shift, issued from BC on bpN and a gauge condition on the tangential shift (based on a parabolic driver).
	const Vector& get_BC_shift(double c_V_lap) const ;
  
    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Spheroid& ) ;	

};

ostream& operator<<(ostream& , const Spheroid& ) ;


#endif
