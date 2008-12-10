/*
 *  Definition of Lorene class Excision_surf, friend class of Spheroid.
 *
 */

/*
 *   Copyright (c) 2008  Jose-Luis Jaramillo & Nicolas Vasset
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

#ifndef __EXCISIONSURF_H_ 
#define __EXCISIONSURF_H_ 

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
class Excision_surf {


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

  

    // Derived data : 
    // ------------
    protected:
	mutable Scalar* p_get_BC_conf_fact_1 ; ///< Source of Neumann boundary condition on \f$ \psi \f$, 
	mutable Scalar* p_get_BC_lapse_1 ;   ///< Source of Dirichlet boundary condition of \f$ N \f$
	mutable Vector* p_get_BC_shift_1 ; ///< Source of Dirichlet BC for the shift vector  \f$ \beta^{i} \f$
	mutable Scalar* p_get_BC_Npsi_1 ; ///<  Source of Neumann boundary condition on \f$ \psi \f$.
	mutable Scalar* p_get_BC_conf_fact_2 ; ///< Source of Neumann boundary condition on \f$ \psi \f$, 
	mutable Scalar* p_get_BC_conf_fact_3 ; ///< Source of Neumann boundary condition on \f$ \psi \f$, 
	mutable Scalar* p_get_BC_conf_fact_4 ; ///< Source of Birichlet boundary condition on \f$ \psi \f$, 
	mutable Scalar* p_get_BC_lapse_2 ;   ///< Source of Dirichlet boundary condition of \f$ N \f$
	mutable Vector* p_get_BC_shift_2 ; ///< Source of Dirichlet BC for the shift vector  \f$ \beta^{i} \f$
	mutable Scalar* p_get_BC_Npsi_2 ; ///<  Source of Neumann boundary condition on \f$ \psi \f$.

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
	 */
	Excision_surf(const Scalar& h_in, const Metric& gij, const Sym_tensor& Kij2, const Scalar& ppsi, const Scalar& nn, const Vector& beta, double timestep, int int_nos) ;
	Excision_surf(const Excision_surf& ) ;		///< Copy constructor

	/// Constructor from a file (see \c sauve(FILE*) )
	Excision_surf(FILE* ) ;    		

	virtual ~Excision_surf() ;			///< Destructor
 

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
	/// Assignment to another Excision_surf
	void operator=(const Excision_surf&) ;	
	

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
	const double get_delta_t() const {return delta_t ;};

	/// Returns the internal number of timesteps for one iteration.
	const double get_no_of_steps() const {return no_of_steps ;};

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

	double set_delta_t() {del_deriv() ; return delta_t ; } ;

	double set_no_of_steps() {del_deriv() ; return no_of_steps ; } ;

	
    // Computational functions
    // -----------------------
    public:


	/// Source for a Neumann BC on the conformal factor,based on the vanishing of the expansion
	const Scalar& get_BC_conf_fact_1() const ;
// Source for an arbitrary Dirichlet BC on the lapse	
	const Scalar& get_BC_lapse_1(double value) const ;

// Source for a global Dirichlet BC on the shift, imposing a conformal Killing symmetry on \f$ \varphi \f$.
	const Vector& get_BC_shift_1(double Omega) const ;
// Source for a Dirichlet arbitrary BC on (N*Psi1)
	const Scalar& get_BC_Npsi_1(double value) const ;
	/// Source for the Dirichlet BC on the conformal factor, based on a parabolic driver for the conformal factor
	const Scalar& get_BC_conf_fact_2(double c_psi_lap, double c_psi_fin, Scalar& expa_fin) const ;
/// Source for the Neumann BC on the conformal factor, based on a parabolic driver for the expansion
	const Scalar& get_BC_conf_fact_3(double c_theta_lap, double c_theta_fin, Scalar& expa_fin) const ;
/// Source for the Dirchlet BC on the conformal factor, based on the consistency condition derived from the trace
	const Scalar& get_BC_conf_fact_4() const ;
/// Source for Dirichlet BC on the lapse, based on a parabolic driver towards arbitrary constant value
	const Scalar& get_BC_lapse_2(double lapse_fin, double c_lapse_lap, double c_lapse_fi) const ;
/// Source for a Dirichlet BC on the shift, based on a Parabolic driver; no assumptions are made except a global conformal Killing symmetry.
	const Vector& get_BC_shift_2(double c_bb_lap, double c_bb_fin, double c_V_lap) const ;

/// Source for the Dirichlet BC on (N*Psi1), based on a parabolic driver.
	const Scalar& get_BC_Npsi_2(double value, double c_npsi_lap, double c_npsi_fin) const ;        

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Spheroid& ) ;	

};

ostream& operator<<(ostream& , const Spheroid& ) ;


#endif
