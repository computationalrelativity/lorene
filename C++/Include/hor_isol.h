/*
 *  Definition of Lorene class Hor_isol
 *
 */

/*
 *   Copyright (c) 2004 Jose Luis Jaramillo
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

#ifndef __Hor_isol_H_ 
#define __Hor_isol_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/09/09 16:49:40  f_limousin
 * First version
 *
 *
 *
 * $Header$
 *
 */

class Sym_tensor_trans ; 
class Sym_tensor ; 
class Vector ; 
class Scalar ; 
class Metric ; 
class Metric_flat ; 
class Base_vect ; 
class Map ; 
class Tbl ;
class Time_slice ;
class Time_slice_conf ;
class Tslice_dirac_max ;

#include "proto.h"
#include "headcpp.h"
#include "cmp.h"
#include "evolution.h"

                    
                    //----------------------------//
                    //       class Hor_isol       //
                    //----------------------------//

/**
 * Spacelike time-slice of an Isolated Horizon in a 3+1 spacetime with conformal decomposition
 * in the maximal slicing and Dirac gauge (*** under development ***)
 * \ingroup (evol)
 * 
 */
class Hor_isol : public Tslice_dirac_max {

  // Data : 
  // -----
 protected: 
  /** The time derivative \f$\dot{\gamma}_{ij} \f$ of the physical metric \f$.
   */
  mutable Evolution_std<Sym_tensor> gam_point_evol ;
  
  /** The time derivative \f$ \dot{\tilde{\gamma}}_{ij} \f$ of the conformal metric \f$.
   */
  mutable Evolution_std<Sym_tensor> gamt_point_evol ;

  //  /** Radius of the horizon 
  //  */
  //  Evolution_std<double> radius_evol ;  // NO SE PUEDE DEFINIR eVOLUTION DE UN DOUBLE: ?

  // Constructors - Destructor
  // -------------------------
 public:
  ///Standard 
  Hor_isol(const Scalar& lapse_in, const Vector& shift_in,
	   const Metric_flat& ff_in, const Scalar& psi_in, 
	   const Sym_tensor_trans& hh_in, const Sym_tensor aa_in, 
	   int depth_in = 3) ;	
  
  Hor_isol(const Hor_isol& ) ;   ///< Copy constructor
  
  virtual ~Hor_isol() ;			///< Destructor
  

  // Mutators / assignment
  // ---------------------
 public:
  /// Assignment to another Hor_isol
  void operator=(const Hor_isol&) ;	
	
  // Accessors
  // ---------
 public:
  // Virtual functions from base class Time_slice_conf:
  // -------------------------------------------------
  
  /** Deviation \f$ h^{ij} \f$ 
   * of the conformal metric \f$ \tilde\gamma^{ij} \f$ from 
   * the flat metric \f$ f^{ij} \f$: 
   * \f$\tilde\gamma^{ij} = f^{ij} + h^{ij} \f$.
   * Returns the value at the current time step (\c jtime ).
   */        
  virtual const Sym_tensor& hh() const ; 
  
  /** Trace \e K of the extrinsic curvature 
   *  at the current time step (\c jtime ).
   * It is null in the present case (maximal slicing)
   */        
  virtual const Scalar& trk() const ; 
  
  /** Vector \f$ H^i = {\cal D}_j \tilde\gamma^{ij} \f$ 
   * which vanishes in Dirac gauge.
   * It is null in the present case...
   */
  virtual const Vector& hdirac() const ; 
  
  // Virtual functions from this Tslice_dir_max:
  // ------------------------------------------
  
  /** Returns the \f$\chi \f$ potential of \f$ \bar{h}^{ij} \f$.
   *
   * It is given by \f$ \chi = r^2 \bar{h}^{rr}\f$.
   */
  virtual const Scalar& khi() const ; 
  
  /** Returns the \f$\mu \f$ potential of \f$ \bar{h}^{ij} \f$.
   *
   * See the documentation of \c Sym_tensor_tt for details.
   */
  virtual const Scalar& mu() const ;
  
  /** Returns the trace, with respect to the flat metric 
   * \c ff , of \f$ h^{ij} \f$.
   */
  virtual const Scalar& trh() const ;




  // Virtual functions of this class:
  // -------------------------------

  // Physical parameters
  //--------------------
 public:
 
   
  /** Vector radial normal */
  Vector radial_vect_hor() ;

  /** Vector beta for boundary conditions in cartesian  */
  Vector beta_bound_cart() ;


  /** Element of area of the horizon */
  Scalar darea_hor()  ;
  
  
  /** Radius of the horizon */
  double radius_hor()  ;

  
  /** Angular momentum (modulo)  */
  double ang_mom_hor()  ;
  
  /**   Mass      */
  double mass_hor()  ;
  
  /** Surface gravity   */
  double kappa_hor() ;
  
  /** Orbital velocity    */
  double omega_hor()  ;
  




  //Computational methods
  //---------------------
 public:
  void init_data(const Sym_tensor& uu, const Scalar& trk_in, 
		 const Scalar& trk_point, double precis = 1.e-12,
		 const Scalar* ener_dens=0x0, const Vector* mom_dens=0x0, 
		 const Scalar* trace_stress=0x0 ) ; 
        

  void init_data_schwar(const Sym_tensor& uu, const Scalar& trk_in, 
		 const Scalar& trk_point, double precis = 1.e-12,
		 const Scalar* ener_dens=0x0, const Vector* mom_dens=0x0, 
		 const Scalar* trace_stress=0x0 ) ; 
        
  void init_data_rot(const Sym_tensor& uu, const Scalar& trk_in, 
		 const Scalar& trk_point, double precis = 1.e-12,
		 const Scalar* ener_dens=0x0, const Vector* mom_dens=0x0, 
		 const Scalar* trace_stress=0x0 ) ; 
        




  //Sources
  //-------
  
  // Source Psi
  Scalar source_psi_hor(const Scalar* ener_dens=0x0, const Vector* mom_dens=0x0, 
		 const Scalar* trace_stress=0x0) ;

  // Source NN
  Scalar source_nn_hor( const Scalar& trk_in, const Scalar* ener_dens=0x0, const Vector* mom_dens=0x0, 
		 const Scalar* trace_stress=0x0) ;

  // Source beta
  Vector source_beta_hor(const Scalar* ener_dens=0x0, const Vector* mom_dens=0x0, 
		 const Scalar* trace_stress=0x0) ;
  



  
  // BOUNDARY CONDITIONS
  //--------------------
  
  /// Dirichlet boundary condition for Psi   
  Valeur boundary_psi_Dir() ;

  /// Neumann boundary condition for Psi   
  Valeur boundary_psi_Neu() ;


  /// Dirichlet boundary condition on nn using the extrinsic curvature
  ///< (No time evolution taken into account! Make this)
  ///<--------------------------------------------------------------------------
  Valeur boundary_nn_Dir_kk() ;

  /// Neumann boundary condition on nn using the extrinsic curvature
  /// (No time evolution taken into account! Make this)
  ///<--------------------------------------------------------------------------
  Valeur boundary_nn_Neu_kk() ;	


  /// Component r of boundary value of beta
  Valeur boundary_beta_r() ;
  
  /// Component theta of boundary value of beta
  Valeur boundary_beta_theta() ;
  
  /// Component phi of boundary value of beta
  Valeur boundary_beta_phi() ;
  
  /// Component x of boundary value of beta
  Valeur boundary_beta_x() ;
  
  /// Component theta of boundary value of beta
  Valeur boundary_beta_y() ;
  
  /// Component phi of boundary value of beta
  Valeur boundary_beta_z() ;




  /*
  /// boundary condition for Q
  
  Valeur boundary_qq() ;

  /// Component r of boundary value of beta
  Valeur boundary_beta_r() ;
  
  /// Component theta of boundary value of beta
  Valeur boundary_beta_theta() ;
  
  /// Component phi of boundary value of beta
  Valeur boundary_beta_phi() ;
  
  /// Component x of boundary value of beta
  Valeur boundary_beta_x() ;
  
  /// Component theta of boundary value of beta
  Valeur boundary_beta_y() ;
  
  /// Component phi of boundary value of beta
  Valeur boundary_beta_z() ;

	/// Dirichlet boundary condition on nn using the extrinsic curvature
	///< (No time evolution taken into account! Make this)
	///<--------------------------------------------------------------------------
	Valeur boundary_nn_Dir_kk() ;

	/// Neumann boundary condition on nn using the extrinsic curvature
	/// (No time evolution taken into account! Make this)
	///<--------------------------------------------------------------------------
	Valeur boundary_dnn_Neu_kk() ;	




  */
 



    // Outputs
    // -------
    protected:
	/// Operator >> (virtual function called by the operator<<). 
	virtual ostream& operator>>(ostream& ) const ;	


  
};
#endif
