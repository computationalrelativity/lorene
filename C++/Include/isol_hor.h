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
 * Spacelike time-slice of an Isolated Horizon in a 3+1 spacetime with conformal decomposition
 * No gauge choice imposed  (*** under development ***)
 * \ingroup (evol)
 * 
 */
class Isol_hor : public Time_slice_conf {

  // Data : 
  // -----
 protected: 

    Map_af& mp ;  /// Affine mapping.
    double radius ; /// Radius of the horizon in LORENE's units.
    double omega ; /// Angular velocity in LORENE's units.

  /// Values at successive time steps of the lapse function \f$ N_{auto} \f$.
  mutable Evolution_std<Scalar> n_auto_evol ; 

  /// Values at successive time steps of the lapse function \f$ N_{comp} \f$.
  mutable Evolution_std<Scalar> n_comp_evol ; 

  /// Values at successive time steps of the conformal factor \f$ \Psi_{auto} \f$.
  mutable Evolution_std<Scalar> psi_auto_evol ; 

  /// Values at successive time steps of the lapse function \f$ \Psi_{comp} \f$.
  mutable Evolution_std<Scalar> psi_comp_evol ; 

  /// Values at successive time steps of the covariant derivative
  /// of the lapse \f$ N \f$.
  mutable Evolution_std<Vector> dn_evol ; 

  /// Values at successive time steps of the covariant derivative
  /// of the conformal factor \f$ \Psi \f$.
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
  
  /// 3 metric tilde
  Metric met_gamt ;
  
  /// Time derivative of the 3-metric tilde
  Sym_tensor gamt_point ;
  
  /// Trace of the extrinsic curvature
  Scalar trK ;
    
  /// Time derivative of the trace of the extrinsic curvature 
  Scalar trK_point ;

  /**
   * Function used to construct the part of \f$K^{ij}\f$ generated by the hole 
   * from the total \f$K^{ij}\f$. Only used for a binary system.
   * 
   * Mainly this \c Scalar  is 1 around the hole and 0 around the companion
   * and the sum of \c decouple for the hole and his companion is 1 
   * everywhere.
   */
  Scalar decouple ;

 
  // Constructors - Destructor
  // -------------------------
 public:

// Standard constructor
  Isol_hor(Map_af& mpi, int depth_in) ;

// Constructor from conformal decomposition
  Isol_hor(Map_af& mpi, const Scalar& lapse_in, const Scalar& psi_in,
	   const Vector& shift_in, const Sym_tensor& aa_in, 
	   const Metric& gamt, const Sym_tensor& gamt_point, 
	   const Scalar& trK, const Scalar& trK_point, 
	   const Metric_flat& ff_in, int depth_in = 3) ;	
  
  Isol_hor(const Isol_hor& ) ;   /// Copy constructor

  Isol_hor (const Map& mp, const Base_vect& triad, 
	    const Metric_flat& ff_in, FILE* fich, 
	    bool partial_read, int depth_in) ;   ///  Constructor from a 
  /// binary file
  
  virtual ~Isol_hor() ;			/// Destructor
  

  // Mutators / assignment
  // ---------------------
 public:
  /// Assignment to another Hor_isol
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
  double get_omega() const {return omega;} ;
  /**
   * Sets the angular velocity to \c ome .
   */
  void set_omega(double ome) {omega = ome ;} ;


  // Accessors
  // ---------
 public:

  /// Lapse function \f$ N_{auto} \f$ at the current time step (\c jtime )
  virtual const Scalar& n_auto() const ;

  /// Lapse function \f$ N_{comp} \f$ at the current time step (\c jtime )
  virtual const Scalar& n_comp() const ;

  /// Conformal factor \f$ \Psi_{auto} \f$ at the current time step (\c jtime )
  virtual const Scalar& psi_auto() const ;

  /// Conformal factor \f$ \Psi_{comp} \f$ at the current time step (\c jtime )
  virtual const Scalar& psi_comp() const ;

  /// Covariant derivative of the lapse function \f$ N \f$ at the 
  /// current time step (\c jtime )
  virtual const Vector& dnn() const ;

  /// Covariant derivative of the conformal factor \f$ \Psi \f$ at the 
  /// current time step (\c jtime )
  virtual const Vector& dpsi() const ;

  /// Shift function \f$ \beta^i_{auto} \f$ at the current time step (\c jtime )
  virtual const Vector& beta_auto() const ;

  /// Shift function \f$ \beta^i_{comp} \f$ at the current time step (\c jtime )
  virtual const Vector& beta_comp() const ;

  /** Conformal representation \f$ A^{ij}_{auto} \f$ of the traceless part
   * of the extrinsic curvature:
   * Returns the value at the current time step (\c jtime ).
   */        
  virtual const Sym_tensor& aa_auto() const ; 

  /** Conformal representation \f$ A^{ij}_{comp} \f$ of the traceless part
   * of the extrinsic curvature:
   * Returns the value at the current time step (\c jtime ).
   */        
  virtual const Sym_tensor& aa_comp() const ; 

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
   * the total \f$\nabla N\f$.
   */
  void n_comp (const Isol_hor& comp) ;
  
  /**
   * Imports the part of \f$\Psi\f$ due to the companion hole \c comp . The 
   * total \f$\Psi\f$ is then calculated.
   * 
   * It also imports the covariant derivative of \f$\Psi\f$ and construct 
   * the total \f$\nabla \Psi\f$.
   */
  void psi_comp (const Isol_hor& comp) ;

  /**
   * Imports the part of \f$\beta^i\f$ due to the companion hole \c comp. The 
   * total \f$\beta^i\f$ is then calculated.
   * 
   */
  void beta_comp (const Isol_hor& comp) ;

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
   * Initiates for a single the black hole.
   * 
   * \b WARNING  It supposes that the boost is zero and should only be 
   * used for an isolated black hole..
   */
  void init_bhole_seul () ;

  
  // Physical parameters
  //--------------------
 public:
 
  
  /// Vector radial normal 
  Vector radial_vect_hor() ;

  /// Vector radial normal tilde 
  Vector tradial_vect_hor() ;

  /// Radial component of the shift respcto to the conformal metric
  Scalar b_tilde() ;

  /// Element of area of the horizon 
  Scalar darea_hor()  ;

  /// Radius of the horizon 
  double radius_hor()  ;
  
  /// Angular momentum (modulo)  
  double ang_mom_hor()  ;
  
  ///   Mass      
  double mass_hor()  ;
  
  /// Surface gravity   
  double kappa_hor() ;
  
  /// Orbital velocity    
  double omega_hor()  ;
  
  /// ADM angular Momentum    
  double ang_mom_adm() ;



  //Computational methods
  //---------------------
 public:
  
  /// function to compute initial data for a single black hole
  void init_data(double precis = 1.e-12,
		 double relax = 1., int niter = 100, double ang_vel = 0.) ; 

  /// function to compute initial data for a single black hole using
  /// Berlin boundary condition.
  void init_data_b_neumann(double precis = 1.e-12,
			   double relax = 1., int niter = 100, double ang_vel = 0.) ; 


  /// function to compute initial data for a single black hole using
  /// Berlin boundary condition.
  void init_data_berlin(double precis = 1.e-12,
			double relax = 1., int niter = 100, double ang_vel = 0.) ; 


  //Sources
  //-------
  
  /// Source Psi
  Scalar source_psi() ;

  /// Source N
  Scalar source_nn() ;

  /// Source beta
  Vector source_beta() ;

  /// Source b_tilde
  Scalar source_b_tilde() ;

  /// Source vector_b
  Vector source_vector_b() ;
  
  
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
  Valeur boundary_nn_Dir_kk() ;

  /// Neumann boundary condition on nn using the extrinsic curvature
  Valeur boundary_nn_Neu_kk() ;	

  /// Dirichlet boundary condition on nn (eff)
  /// dnn + a nn = 0 
  Valeur boundary_nn_Dir_eff(double aa) ;

  /// Neumann boundary condition on nn (eff)
  /// dnn + a nn = 0 
  Valeur boundary_nn_Neu_eff(double aa) ;	

  /// Dirichlet boundary condition nn = aa
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
  
  ///  Vector V^i for boundary conditions in cartesian  
  Vector vv_bound_cart(double velang) ;

  /// Component x of boundary value of v^i
  Valeur boundary_vv_x(double velang) ;

  /// Component y of boundary value of v^i
  Valeur boundary_vv_y(double velang) ;
  
  /// Component z of boundary value of v^i
  Valeur boundary_vv_z(double velang) ;

  /// Neumann boundary condition for b_tilde
  Valeur boundary_b_tilde_Neu() ;

  /// Dirichlet boundary condition for b_tilde
  Valeur boundary_b_tilde_Dir() ;


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
    void sauve(FILE* fich) const ; 
    
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
	// les deux trous noirs.
	Isol_hor hole1 ;	///< Black hole one
	Isol_hor hole2 ;	///< Black hole two
	
	// Tableau sur les deux trous.
	Isol_hor* holes[2] ; ///< Array on the black holes
	
	double omega ;	///< Angular velocity
	
    public:
	// constructeurs & destructeur
	/// Standard constructor.
	Bin_hor(Map_af& mp1, Map_af& mp2, int depth_in) ; 
	Bin_hor(const Bin_hor& ) ;	///< Copy
	~Bin_hor() ;  ///< Destructor
	
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
	 * Calculation af the extrinsic curvature tensor.
	 * 
	 * The regularisation of the shifts must have been done. All the 
	 * contributions of \f$A^{ij}\f$ are then calculated 
	 * and the total tensor 
	 * must be zero on both horizons. The computation is done 
	 * to avoid every singularity 
	 * close to the horizons (division by \e N =0) for it is 
	 * done in the coefficient space 
	 * in the two regions surrounding the holes.
	 */
	void extrinsic_curvature () ;
	
	/**
	 * Calculates {tt decouple} which is used to obtain 
	 * \c tkij_auto  by the formula : 
	 * \c tkij_auto  = \c decouple  * \c tkij_tot .
	 * (see the membre {tt Cmp decouple  for more precisions 
	 * about its value).
	 * 
	 */
	void decouple () ;
    
    public:
	 /**
	  * Initialize the systeme to Misner Lindquist solution, 
	  * that is solving for \e N  and 
	  * \f$Psi\f$ in the case \f$\Omega = 0\f$.
	  * @param precis [input] : precision for the convergence (on \e N ).
	  * @param relax [input] : relaxation parameter.
	  */
	 
	 void set_statiques (double precis, double relax) ;
	 
	 /**
	  * Solves the equation for a particular angular velocity, 
	  * the systeme being 
	  * initialized to Misner-Lindquist solution.
	  * @param angu [input] : angular velocity used for 
	  * the boundary condition on 
	  * \f$\vec{\beta}\f$.
	  * @param precis [input] : precision for the convergence 
	  * (on \f$\beta\f$).
	  * @param relax [input] : relaxation parameter.
	  * @param nbre_ome [input] : number of intermediates 
	  * velocities to go from 0 to 
	  * \c omega , typically 10.
	  * @param sortie [input] : flag for the output on files 
	  * (0 no output files).
	  * @returns : the virial error.
	  */
	  double coal (double angu, double precis, double relax,
	   	    double nbre_ome, const int sortie = 0) ;

	  
	  /**
	   * Solves the equation for the lapse~:
	   * The fields are the total values except those with 
	   * subscript \f$_a\f$, which are 
	   * the fields generated by each holes (\e a  = 1, 2). 
	   * The boundary conditions are 
	   * such that \e N =0 on both horizons.
	   * The companion contributions are then obtained.
	   * @param precis [input] : precision,  for the boudary conditions are
	   * obtained by iteration.
	   * @param relax [input] : relaxation parameter.
	   */
	  void solve_lapse (double precis, double relax) ;
	  
	  /**
	   * Solves the equation for the conformal factor~:
	   * The fields are the total values excpet those with 
	   * subscript \f$_a\f$, which are 
	   * the fields generated by each holes (\e a  = 1, 2). 
	   * @param precis [input] : precision,  for the boudary 
	   * conditions are being 
	   * obtained by iteration.
	   * @param relax [input] : relaxation parameter.
	   */
	  void solve_psi (double precis, double relax) ;
	  
	  /**
	   * Solves the equation for the shift, using the 
	   * Oohara-Nakarmure scheme~:
	   * The fields are the total values excpet those with 
	   * subscript \f$_a\f$, which are 
	   * the fields generated by each holes (\c a = 1, 2). 
	   * The companion contributions are then 
	   * obtained.
	   * @param precis [input] : precision for the solver, the boundary 
	   * conditions and the inversion of the operator being 
	   * obtained by iteration.
	   * @param relax [input] : relaxation parameter.
	   */
	  void solve_shift (double precis, double relax) ;
	  /// ADM angular Momentum    
	  double ang_mom_adm() ;
	  
	  /// Calculates the ADM mass of the system
	  double adm_mass() const ;
	  
	  /// Calculates the Komar mass of the system using :
	  double komar_mass() const ;
	
	  /**
	   * Calculation of the proper distance between the 
	   * two spheres of inversion, 
	   * along the x axis.
	   * @param nr [input] : number of points used for the calculation.
	   */
	  double proper_distance(const int nr = 65) const ;	
} ;



#endif
