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

  /// Radius of the horizon in LORENE's units.
  double radius ; 

  /// Angular velocity in LORENE's units.
  double omega ; 

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
  double get_omega() const {return omega;} ;
  /**
   * Sets the angular velocity to \c ome .
   */
  void set_omega(double ome) {omega = ome ;} ;


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

  
  // Physical parameters
  //--------------------
 public:
 
  
  /// Vector radial normal 
  Vector radial_vect_hor() ;

  /// Vector radial normal tilde 
  Vector tradial_vect_hor() ;

  /// Radial component of the shift with respect to the conformal metric
  Scalar b_tilde() ;

  /// Element of area of the horizon 
  Scalar darea_hor()  ;

  /// Radius of the horizon 
  double radius_hor()  ;
  
  /// Angular momentum (modulo)  
  double ang_mom_hor()  ;
  
  ///   Mass computed at the horizon     
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
  
  /// Source for \f$ \Psi \f$
  Scalar source_psi() ;

  /// Source for \c N
  Scalar source_nn() ;

  /// Source for \f$ \beta \f$
  Vector source_beta() ;

  /// Source for \c b_tilde
  Scalar source_b_tilde() ;

  /// Source for \c vector_b
  Vector source_vector_b() ;
  
  
  // BOUNDARY CONDITIONS
  //--------------------
  
  /// Dirichlet boundary condition for \f$ \Psi \f$ (evolution)
  Valeur boundary_psi_Dir_evol() ;

  /// Neumann boundary condition for \f$ \Psi \f$  (evolution)
  Valeur boundary_psi_Neu_evol() ;

  /// Dirichlet boundary condition for  \f$ \Psi \f$ (spatial)
  Valeur boundary_psi_Dir_spat() ;

  /// Neumann boundary condition for  \f$ \Psi \f$ (spatial)  
  Valeur boundary_psi_Neu_spat() ;

  /// Dirichlet boundary condition for \c N using the extrinsic curvature
  Valeur boundary_nn_Dir_kk() ;

  /// Neumann boundary condition for \c N using the extrinsic curvature
  Valeur boundary_nn_Neu_kk() ;	

  /// Dirichlet boundary condition for \c N (effectif)
  /// \f$ \partial_r N + a N = 0 \f$ 
  Valeur boundary_nn_Dir_eff(double aa) ;

  /// Neumann boundary condition on nn (effectif)
  ///  \f$ \partial_r N + a N = 0 \f$ 
  Valeur boundary_nn_Neu_eff(double aa) ;	

  /// Dirichlet boundary condition \f$ N = a \f$
  Valeur boundary_nn_Dir(double aa) ;

  /// Component r of boundary value of \f$ \beta \f$
  Valeur boundary_beta_r() ;
  
  /// Component theta of boundary value of \f$ \beta \f$
  Valeur boundary_beta_theta() ;
  
  /// Component phi of boundary value of \f$ \beta \f$
  Valeur boundary_beta_phi() ;
  
  /// Component x of boundary value of \f$ \beta \f$
  Valeur boundary_beta_x(double) ;
  
  /// Component theta of boundary value of \f$ \beta \f$
  Valeur boundary_beta_y(double) ;
  
  /// Component phi of boundary value of \f$ \beta \f$
  Valeur boundary_beta_z(double) ;

  /// Vector \f$ \beta^i \f$ for boundary conditions in cartesian  
  Vector beta_bound_cart(double) ;
  
  ///  Vector \f$  V^i \f$ for boundary conditions in cartesian  
  Vector vv_bound_cart(double velang) ;

  /// Component x of boundary value of \f$  V^i \f$
  Valeur boundary_vv_x(double velang) ;

  /// Component y of boundary value of \f$  V^i \f$
  Valeur boundary_vv_y(double velang) ;
  
  /// Component z of boundary value of \f$  V^i \f$
  Valeur boundary_vv_z(double velang) ;

  /// Neumann boundary condition for \c b_tilde
  Valeur boundary_b_tilde_Neu() ;

  /// Dirichlet boundary condition for \c b_tilde
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
	  * @param nb_om [input] : number of intermediates 
	  * velocities to go from 0 to 
	  * \c omega , typically 10.
	  * @param sortie [input] : flag for the output on files 
	  * (0 no output files).
	  * @returns : the virial error.
	  */
	  double coal (double ang_vel, double precis, double relax,
	   	    double nb_om, const int sortie = 0) ;

	  
	  /**
	   * Solves the equation for the lapse :
	   * The fields are the total values except those with 
	   * subscript \f$_a\f$, which are 
	   * the fields generated by each holes (\e a  = 1, 2). 
	   * @param precis [input] : precision,  for the boudary conditions are
	   * obtained by iteration.
	   * @param relax [input] : relaxation parameter.
	   */
	  void solve_lapse (double precis, double relax) ;
	  
	  /**
	   * Solves the equation for the conformal factor :
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
	   * Oohara-Nakarmure scheme :
	   * The fields are the total values excpet those with 
	   * subscript \f$_a\f$, which are 
	   * the fields generated by each holes (\c a = 1, 2). 
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
	  
	  /// Calculates the Komar mass of the system
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
