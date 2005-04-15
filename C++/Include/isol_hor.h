/*
 *  Definition of Lorene class Isol_Hor
 *
 */

/*
 *   Copyright (c) 2004-2005 Jose Luis Jaramillo
 *                      Francois Limousin
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
 * Revision 1.33  2005/04/15 11:20:39  jl_jaramillo
 * Function adapt_hor(double c_min, double c_max)  for adapting a given surface
 * to the excised surface
 *
 * Revision 1.32  2005/04/08 12:15:38  f_limousin
 * Function set_psi(). And dependance in phi.
 *
 * Revision 1.31  2005/04/03 19:48:38  f_limousin
 * Implementation of set_psi(psi_in).
 *
 * Revision 1.30  2005/04/02 15:50:08  f_limousin
 * New data member nz (number of zones). Delete ww.
 *
 * Revision 1.29  2005/03/31 09:48:04  f_limousin
 * New functions compute_ww(..) and aa_kerr_ww() and new data ww.
 *
 * Revision 1.28  2005/03/28 19:45:41  f_limousin
 * Implement Isol_hor::aa_kerr_perturb(...) and new member aa_quad_evol.
 *
 * Revision 1.27  2005/03/24 16:50:40  f_limousin
 * Add parameters solve_shift and solve_psi in par_isol.d and in function
 * init_dat(...). Implement Isolhor::kerr_perturb().
 *
 * Revision 1.26  2005/03/10 16:57:01  f_limousin
 * Improve the convergence of the code coal_bh.
 *
 * Revision 1.25  2005/03/10 10:19:42  f_limousin
 * Add the regularisation of the shift in the case of a single black hole
 * and lapse zero on the horizon.
 *
 * Revision 1.24  2005/03/09 10:28:37  f_limousin
 * Delete functions init_data_b_neumann(...) and init_data_berlin(...)
 * --> New parameter solve_lapse in the function init_data(...).
 * New function update_aa().
 *
 * Revision 1.23  2005/03/06 16:59:33  f_limousin
 * New function Isol_hor::aa() (the one belonging to the class
 * Time_slice_conf need to compute the time derivative of hh and thus
 * cannot work in the class Isol_hor).
 *
 * Revision 1.22  2005/03/04 09:39:31  f_limousin
 * Implement the constructor from a file, operator>>, operator<<
 * and function sauve in the class Bin_hor.
 *
 * Revision 1.21  2005/03/03 10:25:16  f_limousin
 * In the class Isol_hor :
 *    - Add the boost in x and z-direction (members boost_x and boost_z,
 *    and functions get_boost_x(), set_boost_x(double))
 *    - Add function area_hor()
 *    - Put the boundary conditions for the lapse, psi and beta in
 *    the parameter file.
 * In the class bin_hor :
 *    -  Introduce function to compute global quantities as ADM mass,
 *    Komar mass and angular momentum.
 *
 * Revision 1.20  2005/02/24 17:22:53  f_limousin
 * Suppression of the function beta_bound_cart().
 * The boundary conditions for psi, N and beta are now some parameters
 * in par_init.D and par_coal.d.
 *
 * Revision 1.19  2005/02/07 10:30:09  f_limousin
 * Add the regularisation in the case N=0 on the horizon.
 *
 * Revision 1.18  2004/12/31 15:33:37  f_limousin
 * Change the constructor from a file and the standard constructor.
 *
 * Revision 1.17  2004/12/29 16:30:00  f_limousin
 * Improve comments for doxygen
 *
 * Revision 1.16  2004/12/29 16:10:25  f_limousin
 * Add the new class Bin_hor.
 *
 * Revision 1.14  2004/11/24 19:32:05  jl_jaramillo
 * Method for initial data with Berlin boundary conditions
 *
 * Revision 1.13  2004/11/18 10:53:03  jl_jaramillo
 * Declarations for Berlin boundary conditions
 *
 * Revision 1.12  2004/11/05 11:01:13  f_limousin
 * Delete arguments ener_dens, mom_dens and trace_stress in all functions
 * source_nn, source_psi, source_beta, init_data. Delete also
 * argument partial_save in function save.
 *
 * Revision 1.11  2004/11/05 10:11:23  f_limousin
 * The member Metric met_gamt replace Sym_tensor gamt.
 *
 * Revision 1.10  2004/11/03 17:15:46  f_limousin
 * Change the standart constructor. Add 4 memebers : trK, trK_point,
 * gamt and gamt_point.
 * Add also a constructor from a file.
 *
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

#include "time_slice.h"
#include "proto.h"
#include "headcpp.h"
#include "cmp.h"
#include "evolution.h"

                    
//----------------------------//
//       class Isol_Hor       //
//----------------------------//

/**
 * Spacelike time-slice of an Isolated Horizon in a 3+1 spacetime with conformal decomposition.
 * No gauge choice imposed. 
 * \ingroup (evol)
 * 
 */
class Isol_hor : public Time_slice_conf {

  // Data : 
  // -----
 protected: 
  /// Affine mapping.
  Map_af& mp ;  

  /// Number of zones.
  int nz ;

  /// Radius of the horizon in LORENE's units.
  double radius ; 

  /// Angular velocity in LORENE's units.
  double omega ; 
  
  /// Boost velocity in x-direction
  double boost_x ;
 
 /// Boost velocity in z-direction
  double boost_z ;

  /// Intensity of the correction on the shift vector. 
  double regul ; 

  /// Values at successive time steps of the lapse function \f$ N_{auto} \f$.
  mutable Evolution_std<Scalar> n_auto_evol ; 

  /// Values at successive time steps of the lapse function \f$ N_{comp} \f$.
  mutable Evolution_std<Scalar> n_comp_evol ; 

  /// Values at successive time steps of the conformal factor \f$ \Psi_{auto} \f$.
  mutable Evolution_std<Scalar> psi_auto_evol ; 

  /// Values at successive time steps of the lapse function \f$ \Psi_{comp} \f$.
  mutable Evolution_std<Scalar> psi_comp_evol ; 

  /// Values at successive time steps of the covariant derivative 
  /// of the lapse with respect to the flat metric \f$ \overline\nabla_i N \f$.
  mutable Evolution_std<Vector> dn_evol ; 

  /// Values at successive time steps of the covariant derivative
  /// of the conformal factor \f$ \overline\nabla_i \Psi \f$.
  mutable Evolution_std<Vector> dpsi_evol ;   

  /// Values at successive time steps of the shift function \f$ \beta^i_{auto} \f$.
  mutable Evolution_std<Vector> beta_auto_evol ; 

  /// Values at successive time steps of the shift function \f$ \beta^i_{comp} \f$.
  mutable Evolution_std<Vector> beta_comp_evol ; 

  /** Values at successive time steps of the components \f$ A^{ij}_{auto} \f$
   * of the conformal representation of the traceless part
   * of the extrinsic curvature:
   */        
  mutable Evolution_std<Sym_tensor> aa_auto_evol ; 

  /** Values at successive time steps of the components \f$ A^{ij}_{comp} \f$
   * of the conformal representation of the traceless part
   * of the extrinsic curvature:
   */        
  mutable Evolution_std<Sym_tensor> aa_comp_evol ; 
  
  /** Values at successive time steps of the components \f$ A^{ij}*2N \f$
   */        
  mutable Evolution_std<Sym_tensor> aa_nn ; 

  /// Values at successive time steps of the components \f$ A^{ij}A_{ij} \f$
  mutable Evolution_std<Scalar> aa_quad_evol ;

  /// 3 metric tilde
  Metric met_gamt ;
  
  /// Time derivative of the 3-metric tilde
  Sym_tensor gamt_point ;
  
  /// Trace of the extrinsic curvature
  Scalar trK ;
    
  /// Time derivative of the trace of the extrinsic curvature 
  Scalar trK_point ;

  /**
   * Function used to construct \f$ A^{ij}_{auto} \f$
   * from the total \f$A^{ij}\f$. Only used for a binary system.
   * 
   * Mainly this \c Scalar  is 1 around the hole and 0 around the companion
   * and the sum of \c decouple for the hole and his companion is 1 
   * everywhere.
   */
  Scalar decouple ;

  // Constructors - Destructor
  // -------------------------
 public:

  /** Standard constructor
   *  @param mpi affine mapping
   *  @param depth_in  number of stored time slices; this parameter is used
   *                   to set the \c scheme_order member with \c scheme_order
   *                   = \c depth_in - 1. \c scheme_order can be changed 
   *                   afterwards by the method \c set_scheme_order(int).
   */
  Isol_hor(Map_af& mpi, int depth_in = 3) ;

  /** Constructor from conformal decomposition
   *
   *  @param mpi affine mapping
   *  @param lapse_in lapse function \e N
   *  @param psi_in conformal factor \f$\Psi\f$ relating the
   *       physical metric \f$ \gamma_{ij} \f$ to the conformal one:
   *      \f$ \gamma_{ij} = \Psi^4 \tilde\gamma_{ij} \f$
   *  @param shift_in shift vector
   *  @param aa_in conformal representation \f$ A^{ij} \f$
   *      of the traceless part of the extrinsic curvature:
   *   \f$ A^{ij} = \Psi^4 \left( K^{ij} - \frac{1}{3} K \gamma^{ij} \right) \f$
   *  @param gamt 3-metric tilde 
   *  @param gamt_point time derivative of the 3-metric tilde    
   *  @param trK trace \e K of the extrinsic curvature 
   *  @param trK_point time derivative of the trace \e K of the extrinsic curvature 
   *  @param ff_in reference flat metric with respect to which the
   *           conformal decomposition is performed
   *  @param depth_in  number of stored time slices; this parameter is used
   *                   to set the \c scheme_order member with \c scheme_order
   *                   = \c depth_in - 1. \c scheme_order can be changed 
   *                   afterwards by the method \c set_scheme_order(int).
   */
  Isol_hor(Map_af& mpi, const Scalar& lapse_in, const Scalar& psi_in,
	   const Vector& shift_in, const Sym_tensor& aa_in, 
	   const Metric& gamt, const Sym_tensor& gamt_point, 
	   const Scalar& trK, const Scalar& trK_point, 
	   const Metric_flat& ff_in, int depth_in = 3) ;	
  
/// Copy constructor
  Isol_hor(const Isol_hor& ) ;  

  /**  Constructor from a binary file
   *  @param mpi affine mapping
   *  @param fich file containing the saved \c isol_hor
   *  @param partial_read indicates whether the full object must be read in 
   *                      file or whether the final construction is devoted 
   *                      to a constructor of a derived class
   *  @param depth_in  number of stored time slices; this parameter is used
   *                   to set the \c scheme_order member with \c scheme_order
   *                   = \c depth_in - 1. \c scheme_order can be changed 
   *                   afterwards by the method \c set_scheme_order(int).
   */
  Isol_hor (Map_af& mp, FILE* fich, 
	    bool partial_read, int depth_in = 3) ;  
  
/// Destructor
  virtual ~Isol_hor() ;			
  

  // Mutators / assignment
  // ---------------------
 public:
  /// Assignment to another Isol_hor
  void operator=(const Isol_hor&) ;	
	

 public:
  /// Returns the mapping (readonly).
  const Map_af& get_mp() const {return mp;} ; 
  
  /// Read/write of the mapping
  Map_af& set_mp() {return mp; } ;

  /**
   * Returns the radius of the horizon.
   */
  double get_radius() const {return radius;} ;
  
  /**
   * Sets the radius of the horizon to \c rad .
   */
  void set_radius(double rad) {radius = rad ;} ;
	
  /**
   * Returns the angular velocity.
   */
  double get_omega() const {return omega ;} ;
  /**
   * Sets the angular velocity to \c ome .
   */
  void set_omega(double ome) {omega = ome ;} ;

  /**
   * Returns the boost velocity in x-direction.
   */
  double get_boost_x() const {return boost_x ;} ;
  /**
   * Sets the boost velocity in x-direction to \c bo .
   */
  void set_boost_x(double bo) {boost_x = bo ;} ;

  /**
   * Returns the boost velocity in z-direction.
   */
  double get_boost_z() const {return boost_z ;} ;
  /**
   * Sets the boost velocity in z-direction to \c bo .
   */
  void set_boost_z(double bo) {boost_z = bo ;} ;


  
  // Accessors
  // ---------
 public:

  /// Lapse function \f$ N_{auto} \f$ at the current time step \c jtime 
  virtual const Scalar& n_auto() const ;

  /// Lapse function \f$ N_{comp} \f$ at the current time step \c jtime 
  virtual const Scalar& n_comp() const ;

  /// Conformal factor \f$ \Psi_{auto} \f$ at the current time step \c jtime 
  virtual const Scalar& psi_auto() const ;

  /// Conformal factor \f$ \Psi_{comp} \f$ at the current time step \c jtime 
  virtual const Scalar& psi_comp() const ;

  /// Covariant derivative of the lapse function \f$ \overline\nabla_i N \f$ at the 
  /// current time step \c jtime 
  virtual const Vector& dnn() const ;

  /// Covariant derivative with respect to the flat metric
  /// of the conformal factor \f$ \overline\nabla_i \Psi \f$ at the 
  /// current time step \c jtime 
  virtual const Vector& dpsi() const ;

  /// Shift function \f$ \beta^i_{auto} \f$ at the current time step \c jtime 
  virtual const Vector& beta_auto() const ;

  /// Shift function \f$ \beta^i_{comp} \f$ at the current time step \c jtime 
  virtual const Vector& beta_comp() const ;

  /** Conformal representation \f$ A^{ij}_{auto} \f$ of the traceless part
   * of the extrinsic curvature:
   * Returns the value at the current time step \c jtime.
   */        
  virtual const Sym_tensor& aa_auto() const ; 

  /** Conformal representation \f$ A^{ij}_{comp} \f$ of the traceless part
   * of the extrinsic curvature:
   * Returns the value at the current time step \c jtime.
   */        
  virtual const Sym_tensor& aa_comp() const ; 

  /** Conformal representation \f$ A^{ij}A_{ij} \f$.
   * Returns the value at the current time step \c jtime.
   */        
  virtual const Scalar& aa_quad() const ; 

  /** Conformal metric 
   * \f$ \tilde\gamma_{ij} = \Psi^{-4} \gamma_{ij} \f$
   * Returns the value at the current time step (\c jtime ).
   */        
  virtual const Metric& tgam() const {return met_gamt ;}

  /**
   * Returns the function used to construct \c tkij_auto  from \c tkij_tot .
   */
  const Scalar get_decouple() const {return decouple ;}


 public:
  /**
   * Imports the part of \e N  due to the companion hole \c comp . The 
   * total \e N  is then calculated.
   * 
   * It also imports the covariant derivative of \e N  and construct 
   * the total \f$\nabla_i N\f$.
   */
  void n_comp (const Isol_hor& comp) ;
  
  /**
   * Imports the part of \f$\Psi\f$ due to the companion hole \c comp . The 
   * total \f$\Psi\f$ is then calculated.
   * 
   * It also imports the covariant derivative of \f$\Psi\f$ and construct 
   * the total \f$\nabla_i \Psi\f$.
   */
  void psi_comp (const Isol_hor& comp) ;

  /**
   * Imports the part of \f$\beta^i\f$ due to the companion hole \c comp. The 
   * total \f$\beta^i\f$ is then calculated.
   * 
   */
  void beta_comp (const Isol_hor& comp) ;

  /**
   * Computes the viriel error, that is the difference between the ADM 
   * and the Komar 
   * masses,  calculated by the asymptotic behaviours of 
   * respectively \f$\Psi\f$ and \e N .
   * 
   * \b WARNING  this should only be used for an isolated black hole.
   */
  double viriel_seul () const ;

  /**
   * Sets the values of the fields to :
   * \li \c n_auto  \f$= \frac{1}{2}-2\frac{a}{r}\f$
   * \li \c n_comp   \f$= \frac{1}{2}\f$
   * \li \c psi_auto  \f$= \frac{1}{2}+\frac{a}{r}\f$
   * \li \c psi_comp   \f$= \frac{1}{2}\f$
   * 
   * \e a  being the radius of the hole, the other fields being set to zero.
   */
  void init_bhole () ;
  
  /**
   * Sets the 3-metric tilde to the flat metric and gamt_point,
   * trK and trK_point to zero.
   */

  void init_met_trK() ;

  /**
   * Initiates for a single black hole.
   * 
   * \b WARNING  It supposes that the boost is zero and should only be 
   * used for an isolated black hole..
   */
  void init_bhole_seul () ;

  /** Sets the conformal factor \f$ \Psi \f$ relating the
   * physical metric \f$ \gamma_{ij} \f$ to the conformal one:
   * \f$ \gamma_{ij} = \Psi^4 \tilde\gamma_{ij} \f$. 
   * \f$ \Psi \f$ is defined by
   *  \f[ \Psi := \left( \frac{\det\gamma_{ij}}{\det f_{ij}} 
   *      \right) ^{1/12} \f] 
   * Sets the value at the current time step (\c jtime ) and
   * delete all quantities which depend on \f$ \Psi \f$.
   */

  void set_psi(const Scalar& psi_in) ; 

  /// Sets the lapse
  void set_nn(const Scalar& nn_in) ; 
  
  /// Sets the conformal metric to gam_tilde.
  void set_gamt(const Metric& gam_tilde) ;

  // Physical parameters
  //--------------------
 public:
 
  
  /// Vector radial normal 
  const Vector radial_vect_hor() const ;

  /// Vector radial normal tilde 
  const Vector tradial_vect_hor() const ;

  /// Radial component of the shift with respect to the conformal metric
  const Scalar b_tilde() const ;

  /// Element of area of the horizon 
  const Scalar darea_hor() const ;

  /// Area of the horizon 
  double area_hor() const ;

  /// Radius of the horizon 
  double radius_hor() const ;
  
  /// Angular momentum (modulo)  
  double ang_mom_hor() const ;
  
  ///   Mass computed at the horizon     
  double mass_hor() const ;
  
  /// Surface gravity   
  double kappa_hor() const ;
  
  /// Orbital velocity    
  double omega_hor() const ;
  
  /// ADM angular Momentum    
  double ang_mom_adm() const ;

  /// Expansion of the outgoing null normal (\f$ \bf n + \bf s \f$)
  Scalar expansion() const ;


  //Computational methods
  //---------------------
 public:
  
  /* function to compute initial data for a single black hole
   *  @param bound_nn boundary condition for the lapse
   *  @param lim_nn value of the boundary condition for the lapse 
   *  @param bound_psi boundary condition for \f$ \Psi \f$
   *  @param bound_beta boundary condition for the shift
   *  @param solve_lapse do we solve the equation for the lapse ?
   *  @param precis precision for the convergence
   *  @param relax relaxation
   *  @param niter number of iterations
   */
  void init_data(int bound_nn, double lim_nn, int bound_psi, int bound_beta,
		 int solve_lapse, int solve_psi, int solve_shift, 
		 double precis = 1.e-12, double relax = 1., int niter = 100) ; 

  //Sources
  //-------
  
  /// Source for \f$ \Psi \f$
  const Scalar source_psi() const ;

  /// Source for \c N
  const Scalar source_nn() const ;

  /// Source for \f$ \beta \f$
  const Vector source_beta() const ;

  /// Source for \c b_tilde
  const Scalar source_b_tilde() const ;

  /// Source for \c vector_b
  const Vector source_vector_b() const ;
  
  
  // BOUNDARY CONDITIONS
  //--------------------
  
  /// Dirichlet boundary condition for \f$ \Psi \f$ (evolution)
  const Valeur boundary_psi_Dir_evol() const ;

  /// Neumann boundary condition for \f$ \Psi \f$  (evolution)
  const Valeur boundary_psi_Neu_evol() const ;

  /// Dirichlet boundary condition for  \f$ \Psi \f$ (spatial)
  const Valeur boundary_psi_Dir_spat() const ;

  /// Neumann boundary condition for  \f$ \Psi \f$ (spatial)  
  const Valeur boundary_psi_Neu_spat() const ;

  /// Neumann boundary condition for  \f$ \Psi \f$ (spatial)  
  const Valeur boundary_psi_app_hor() const ;

  /// Dirichlet boundary condition for \c N using the extrinsic curvature
  const Valeur boundary_nn_Dir_kk() const ;

  /// Neumann boundary condition for \c N using the extrinsic curvature
  const Valeur boundary_nn_Neu_kk() const ;	

  /// Dirichlet boundary condition for \c N (effectif)
  /// \f$ \partial_r N + a N = 0 \f$ 
  const Valeur boundary_nn_Dir_eff(double aa) const ;

  /// Neumann boundary condition on nn (effectif)
  ///  \f$ \partial_r N + a N = 0 \f$ 
  const Valeur boundary_nn_Neu_eff(double aa) const ;	

  /// Dirichlet boundary condition \f$ N = a \f$
  const Valeur boundary_nn_Dir(double aa) const ;

  /// Component r of boundary value of \f$ \beta \f$
  const Valeur boundary_beta_r() const ;
  
  /// Component theta of boundary value of \f$ \beta \f$
  const Valeur boundary_beta_theta() const ;
  
  /// Component phi of boundary value of \f$ \beta \f$
  const Valeur boundary_beta_phi() const ;
  
  /// Component x of boundary value of \f$ \beta \f$
  const Valeur boundary_beta_x(double om) const ;
  
  /// Component y of boundary value of \f$ \beta \f$
  const Valeur boundary_beta_y(double om) const ;
  
  /// Component z of boundary value of \f$ \beta \f$
  const Valeur boundary_beta_z() const ;

  /// Boundary value for a boost in x-direction
  const Valeur beta_boost_x() const ;  

  /// Boundary value for a boost in z-direction
  const Valeur beta_boost_z() const ;  

  ///  Vector \f$  V^i \f$ for boundary conditions in cartesian  
  const Vector vv_bound_cart(double om) const ;

  /// Component x of boundary value of \f$  V^i \f$
  const Valeur boundary_vv_x(double om) const ;

  /// Component y of boundary value of \f$  V^i \f$
  const Valeur boundary_vv_y(double om) const ;
  
  /// Component z of boundary value of \f$  V^i \f$
  const Valeur boundary_vv_z(double om) const ;

  /// Neumann boundary condition for \c b_tilde
  const Valeur boundary_b_tilde_Neu() const ;

  /// Dirichlet boundary condition for \c b_tilde
  const Valeur boundary_b_tilde_Dir() const ;
 
  /** Conformal representation \f$ A^{ij} \f$ of the traceless part
   * of the extrinsic curvature:
   * \f$ A^{ij} = \Psi^4 \left( K^{ij} - \frac{1}{3} K \gamma^{ij} \right) \f$.
   */        
  void update_aa() ; 
  
  /**
   * Corrects \c shift_auto  in such a way that the total \f$A^{ij}\f$ is 
   * equal to zero in the horizon,  which should ensure the regularity 
   * of \f$K^{ij}\f$.
   * 
   * \b WARNING :  this should only be used for a black hole in 
   * a binary system \c Bin_hor.
   * 
   * @param comp [input]: the part of \f$\beta^i\f$ generated by the companion 
   * hole.
   */
  double regularisation (const Vector& shift_auto, const Vector& shift_comp, 
			 double ang_vel) ;

  /**
   * Corrects the shift in the innermost shell, so that it remains \f$
   * {\mathcal{C}}^2\f$ and that \f$A^{ij}\f$ equals zero on the horizon.
   * 
   * return  the relative difference between the shift before
   * and after the regularisation.
   * 
   * \b WARNING  this should only be used for an isolated black hole.
   */
  double regularise_one() ;

  /// Initialisation of the metric tilde from equation (15) of
  /// Dain (2002). The determinant of this conformal metric is not one.
  void met_kerr_perturb() ;
  
  /* Perturbation of Kerr using  \f$ A^{ij}A_{ij} \f$ from 
   * equation (14) of Dain (2002).
   * @param mm mass of the Kerr black hole metric.
   * @param aa rotation parameter of the Kerr black hole metric.
   */
  void aa_kerr_ww(double mm, double aa) ;

  /*  Calculation of the outermost trapped surface and adaptation
   *  of all necessary quantities
   */
  void adapt_hor(double c_min, double c_max) ;
  

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
    
  friend class Bin_hor ; /// Binary systems 

};

/**
 * Binary black holes system. \ingroup (star) 
 * 
 * This class is intended for dealing with binary black holes configurations 
 * in the conformaly flat approximation.
 */


class Bin_hor {
    
    // data :
    private:
	// lThe two black holes
	Isol_hor hole1 ;	///< Black hole one
	Isol_hor hole2 ;	///< Black hole two
	
	///Array on the black holes
	Isol_hor* holes[2] ; 
	
	double omega ;	///< Angular velocity
	
    public:

	/** Standard constructor
	 *  @param mp1 affine mapping for the first black hole
	 *  @param mp2 affine mapping for the second black hole
	 *  @param depth_in  number of stored time slices; this parameter is 
	 *  used to set the \c scheme_order member with \c scheme_order
	 *  = \c depth_in - 1. \c scheme_order can be changed 
	 *  afterwards by the method \c set_scheme_order(int).
	 */
	Bin_hor(Map_af& mp1, Map_af& mp2, int depth_in) ; 

	Bin_hor(const Bin_hor& ) ;	///< Copy constructor

	/**  Constructor from a binary file
	 *  @param mp1 affine mapping of the first black hole
	 *  @param mp2 affine mapping of the second black hole
	 *  @param fich file containing the saved \c Bin_hor
	 *  @param partial_read indicates whether the full object must be 
	 *  read in file or whether the final construction is devoted 
	 *  to a constructor of a derived class
	 *  @param depth_in  number of stored time slices; this parameter 
	 *  is used to set the \c scheme_order member with \c scheme_order
	 *  = \c depth_in - 1. \c scheme_order can be changed 
	 *  afterwards by the method \c set_scheme_order(int).
	 */
	Bin_hor (Map_af& mp1, Map_af& mp2, FILE* fich, 
		 bool partial_read, int depth_in = 3) ;  
	
	virtual ~Bin_hor() ;  ///< Destructor
	
	// Outputs
	// -------
      private:
	/// Operator >> (function called by the operator<<). 
	ostream& operator>>(ostream& ) const ;	

      public :
	  /** Total or partial saves in a binary file.
	   *  
	   *  @param fich binary file 
	   *  @param partial_save indicates whether the whole object must be
	   *      saved.
	   */
	  void sauve(FILE* fich, bool partial_save) const ; 
      
      /// Display
      friend ostream& operator<<(ostream& , const Bin_hor& ) ;	

       public:
	
	void operator=(const Bin_hor&) ; ///< Affectation operator
	
	/**
	 * Read/write of a component of the system. \c i  must be equal to
	 * 1 or 2.
	 */
	Isol_hor& set(int i) 
	    { assert( (i==1) || (i==2) ); 
	      return *holes[i-1] ;} ; 
	/**
	 * Sets the orbital velocity to \c ome
	 */
	void set_omega(double ome) {omega = ome ; 
			     hole1.set_omega (ome) ;
			     hole2.set_omega (ome) ;} ;
	
    public:
	/**
	 * Read only of a component of the system. \c i must be equal to
	 * 1 or 2.
	 */	const Isol_hor& operator()(int i) const 
	    { assert( (i==1) || (i==2) ); 
	      return *holes[i-1] ;} ;
	      
	/// Returns the angular velocity 
	double get_omega() const {return omega; } ; 
    
	/**
	 * Initialisation of the system. Each hole is set close to a 
	 * Schwarzschild one 
	 * and the parts of the fields generated by 
	 * the companion are calculated.
	 * 
	 * The angular velocity is set to zero.
	 */
	void init_bin_hor() ;
	
	/**
	 * Computes the viriel error, that is the difference between 
	 * the ADM and the Komar 
	 * masses,  calculated by the asymptotic behaviours of 
	 * respectively \f$\Psi\f$ and \e N .
	 */
	double viriel() const ;

	/**
	 * Calculation of the extrinsic curvature tensor.
	 * \c aa_auto_evol and \c aa_auto_comp are also computed.
	 */
	void extrinsic_curvature () ;
	
	/**
	 * Calculates \c decouple which is used to obtain 
	 * \c tkij_auto  and \c tkij_comp
	 */
	void decouple () ;
    
    public:
	 /**
	  * Initialize the systeme to Misner Lindquist solution, 
	  * that is solving for \e N  and 
	  * \f$\Psi\f$ in the case \f$\Omega = 0\f$.
	  * @param precis [input] : precision for the convergence (on \e N ).
	  * @param relax [input] : relaxation parameter.
	  * @param bound_nn [input] : type of the boundary condition for 
	  * the lapse.
	  * @param lim_nn [input] : value (double) of the coefficient for
	  * the boundary condition. Only used for boundary_nn_Dir(double),
	  * boundary_nn_Neu_eff(double) and boundary_nn_Dir_eff(double). 
	  * @param bound_psi [input] : type of the boundary condition for 
	  * psi.
	  */
	 
	void set_statiques (double precis, double relax, int bound_nn,
			    double lim_nn, int bound_psi) ;
	 
	 /**
	  * Solves the equation for a particular angular velocity, 
	  * the systeme being 
	  * initialized to Misner-Lindquist solution.
	  * @param angu [input] : angular velocity used for 
	  * the boundary condition on 
	  * \f$\vec{\beta}\f$.
	  * @param relax [input] : relaxation parameter.
	  * @param nb_om [input] : number of intermediates 
	  * velocities to go from 0 to 
	  * \c omega , typically 10.
	  * @param nb_it [input] : number of iteration when omega is fixed 
	  * @param bound_nn [input] : type of the boundary condition for 
	  * the lapse.
	  * @param lim_nn [input] : value (double) of the coefficient for
	  * the boundary condition. Only used for boundary_nn_Dir(double),
	  * boundary_nn_Neu_eff(double) and boundary_nn_Dir_eff(double). 
	  * @param bound_psi [input] : type of the boundary condition for 
	  * psi.
	  * @param bound_nn [input] : type of the boundary condition for 
	  * the shift.
	  * @param step current step of the iteration 
	  * @param sortie [input] : flag for the output on files 
	  * (0 no output files).
	  * @returns : the virial error.
	  */
	double coal (double ang_vel, double relax, int nb_om,
		     int nb_it, int bound_nn, double lim_nn, 
		     int bound_psi, int bound_beta, 
		     ostream& fich_iteration, ostream& fich_correction,
		     ostream& fich_viriel, int step, const int sortie = 0) ;

	  
	  /**
	   * Solves the equation for the lapse :
	   * The fields are the total values except those with 
	   * subscript \f$_a\f$, which are 
	   * the fields generated by each holes (\e a  = 1, 2). 
	   * @param precis [input] : precision,  for the boudary conditions are
	   * obtained by iteration.
	   * @param relax [input] : relaxation parameter.
	   * @param bound_nn [input] : type of the boundary condition at the
	   * horizon.
	   * @param lim_nn [input] : value (double) of the coefficient for
	   * the boundary condition. Only used for boundary_nn_Dir(double),
	   * boundary_nn_Neu_eff(double) and boundary_nn_Dir_eff(double). 
	   */
	  void solve_lapse (double precis, double relax, int bound_nn,
			    double lim_nn) ;
	  
	  /**
	   * Solves the equation for the conformal factor :
	   * The fields are the total values excpet those with 
	   * subscript \f$_a\f$, which are 
	   * the fields generated by each holes (\e a  = 1, 2). 
	   * @param precis [input] : precision,  for the boudary 
	   * conditions are being 
	   * obtained by iteration.
	   * @param relax [input] : relaxation parameter.
	   * @param bound_psi [input] : type of the boundary condition at the
	   * horizon.
	   */
	  void solve_psi (double precis, double relax, int bound_psi) ;
	  
	  /**
	   * Solves the equation for the shift, using the 
	   * Oohara-Nakarmure scheme :
	   * The fields are the total values excpet those with 
	   * subscript \f$_a\f$, which are 
	   * the fields generated by each holes (\c a = 1, 2). 
	   * @param precis [input] : precision for the solver, the boundary 
	   * conditions and the inversion of the operator being 
	   * obtained by iteration.
	   * @param relax [input] : relaxation parameter.
	   * @param bound_beta [input] : type of the boundary condition at the
	   * horizon.
	   */
	  void solve_shift (double precis, double relax, int bound_beta) ;
	
	  /**
	   *  Calculates the ADM mass of the system.
	   */
	  double adm_mass() const ;
	  
	  /**
	   *  Calculates the Komar mass of the system using :
	   * \f$M = \frac{1}{4 \pi} \oint_{\infty} D^i N {\mathrm d} S_i\f$.
	   */
	  double komar_mass() const ;
	  
	  /**
	   * Calculates the angular momentum of the black hole.
	   */
	  double ang_mom_adm() ;

	  /**
	   * Calculation of the proper distance between the 
	   * two spheres of inversion, 
	   * along the x axis.
	   * @param nr [input] : number of points used for the calculation.
	   */
	  double proper_distance(const int nr = 65) const ;	
} ;



#endif
