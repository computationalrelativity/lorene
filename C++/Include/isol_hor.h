/*
 *  Definition of Lorene class Isol_Hor
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

#ifndef __Isol_hor_H_ 
#define __Isol_hor_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.9  2004/11/02 17:42:33  f_limousin
 * New method sauve(...) to save in a binary file.
 *
 * Revision 1.8  2004/11/02 16:15:12  f_limousin
 * Add new argument ang_vel in function init_dat(...).
 *
 * Revision 1.7  2004/10/29 15:46:14  jl_jaramillo
 * Remove 2 members, add ADM angular momentum and change name
 * of functions.
 *
 * Revision 1.6  2004/10/01 16:51:16  f_limousin
 * Pure Dirichlet boundary condition added
 *
 * Revision 1.5  2004/09/28 16:03:58  f_limousin
 * Add parameter niter in the parameter file par_hor.d. It appears in
 * argument of the function init_data_schwarz(...).
 *
 * Revision 1.4  2004/09/17 13:35:25  f_limousin
 * Introduction of relaxation in init_data_schwarz
 *
 * Revision 1.3  2004/09/16 08:36:57  f_limousin
 * New boundary conditions for lapse and psi.
 *
 * Revision 1.2  2004/09/09 17:04:27  jl_jaramillo
 * Elimination of _ih
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
//class Tslice_dirac_max ;

#include "proto.h"
#include "headcpp.h"
#include "cmp.h"
#include "evolution.h"

                    
                    //----------------------------//
                    //       class Isol_Hor       //
                    //----------------------------//

/**
 * Spacelike time-slice of an Isolated Horizon in a 3+1 spacetime with conformal decomposition
 * No gauge choice imposed  (*** under development ***)
 * \ingroup (evol)
 * 
 */
class Isol_hor : public Time_slice_conf {

  // Data : 
  // -----
 protected: 


  // Constructors - Destructor
  // -------------------------
 public:
  ///Standard 
  Isol_hor(const Scalar& lapse_in, const Vector& shift_in,
	   const Sym_tensor& gamma_in, const Sym_tensor kk_in, 
	   const Metric_flat& ff_in, int depth_in = 3) ;	
  
  Isol_hor(const Isol_hor& ) ;   /// Copy constructor

  Isol_hor (const Map& mp, const Base_vect& triad, 
	    const Metric_flat& ff_in, FILE* fich, 
	    bool partial_read, int depth_in) ;   ///  Constructor from a 
                                                  ///  binary file
  
  virtual ~Isol_hor() ;			/// Destructor
  

  // Mutators / assignment
  // ---------------------
 public:
  /// Assignment to another Hor_isol
  void operator=(const Isol_hor&) ;	
	

  // Physical parameters
  //--------------------
 public:
 
   
  /** Vector radial normal */
  Vector radial_vect_hor() ;

  /** Vector radial normal tilde */
  Vector tradial_vect_hor() ;

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
  
  /** ADM angular Momentum    */
  double ang_mom_adm() ;



  //Computational methods
  //---------------------
 public:
  
  void init_data(const Sym_tensor& uu, const Scalar& trk_in, 
		 const Scalar& trk_point, double precis = 1.e-12,
		 double relax = 1., int niter = 100, double ang_vel = 0.,
		 const Scalar* ener_dens=0x0, const Vector* mom_dens=0x0, 
		 const Scalar* trace_stress=0x0 ) ; 
        


  //Sources
  //-------
  
  // Source Psi
  Scalar source_psi(const Scalar* ener_dens=0x0, const Vector* mom_dens=0x0, 
		 const Scalar* trace_stress=0x0) ;

  // Source NN
  Scalar source_nn( const Scalar& trk_in, const Scalar* ener_dens=0x0, const Vector* mom_dens=0x0, 
		 const Scalar* trace_stress=0x0) ;

  // Source beta
  Vector source_beta(const Scalar* ener_dens=0x0, const Vector* mom_dens=0x0, 
		 const Scalar* trace_stress=0x0) ;
  

  // BOUNDARY CONDITIONS
  //--------------------
  
  /// Dirichlet boundary condition for Psi (evolution)
  Valeur boundary_psi_Dir_evol() ;

  /// Neumann boundary condition for Psi (evolution)
  Valeur boundary_psi_Neu_evol() ;

  /// Dirichlet boundary condition for Psi (spatial)
  Valeur boundary_psi_Dir_spat() ;

  /// Neumann boundary condition for Psi (spatial)  
  Valeur boundary_psi_Neu_spat() ;

  /// Dirichlet boundary condition on nn using the extrinsic curvature
  /// ----------------------------------------------------------------
  Valeur boundary_nn_Dir_kk() ;

  /// Neumann boundary condition on nn using the extrinsic curvature
  /// --------------------------------------------------------------
  Valeur boundary_nn_Neu_kk() ;	

  /// Dirichlet boundary condition on nn (eff)
  /// dnn + a nn = 0 
  /// ----------------------------------------------------------------
  Valeur boundary_nn_Dir_eff(double aa) ;

  /// Neumann boundary condition on nn (eff)
  /// dnn + a nn = 0 
  ///---------------------------------------------------------------
  Valeur boundary_nn_Neu_eff(double aa) ;	

  /// Dirichlet boundary condition nn = aa
  /// ----------------------------------------------------------------
  Valeur boundary_nn_Dir(double aa) ;


  /// Component r of boundary value of beta
  Valeur boundary_beta_r() ;
  
  /// Component theta of boundary value of beta
  Valeur boundary_beta_theta() ;
  
  /// Component phi of boundary value of beta
  Valeur boundary_beta_phi() ;
  
  /// Component x of boundary value of beta
  Valeur boundary_beta_x(double) ;
  
  /// Component theta of boundary value of beta
  Valeur boundary_beta_y(double) ;
  
  /// Component phi of boundary value of beta
  Valeur boundary_beta_z(double) ;

  /// Vector beta for boundary conditions in cartesian  
  Vector beta_bound_cart(double) ;


  // Outputs
  // -------
 protected:
  /// Operator >> (virtual function called by the operator<<). 
  virtual ostream& operator>>(ostream& ) const ;	
  

  public :
  /** Total or partial saves in a binary file.
   *  
   *  @param fich binary file 
   *  @param partial_save indicates whether the whole object must be
   *      saved.
   */
  virtual void sauve(FILE* fich, bool partial_save) const ; 
    
  
};
#endif
