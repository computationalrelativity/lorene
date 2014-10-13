/*
 *  Definition of Lorene class Excised slice
 *				 
 */

/*
 *   Copyright (c) 2010 Nicolas Vasset
 */


#ifndef __EXCISEDSLICE_H_ 
#define __EXCISEDSLICE_H_ 


// Headers Lorene
#include "spheroid.h"

//---------------------------//
//    Class Excised_slice    //
//---------------------------//

namespace Lorene {
/** Class to compute single black hole spacetime excised slices.
 *  It takes as arguments:
 * - A Mapping \c mp
 * - The type of physical surface used to perform the excision:
 *                   - 0. General
 *                   - 1. Isolated horizon-like
 *                   - 2. Dynamical horizon-like
 *                   - 3. Inner trapped surface
 * - The type of metric fields that live in the slice:
 *                   - 1. FCF quantities are used (Bonazzola et al., 2004)
 *                   - 2. FCF-New scheme quantities are used (Cordero et al., 2009)
 * Members are obviously the metric fields living in the slice, and related variables.
 * The aim of this class is to coordinate Resolution of Einstein equations with the right 
 * set of boundary conditions, using classes like Spheroid() and Excision_surf();
 * Disclaimer: the class Isol_hole() is redundant with this class under a set of parameters;
 * therefore it has to go at some point.
 */

class Excised_slice {

  // Data : 
  // -----
 protected:
   const Map& mp ;	    ///< Mapping associated with the slice
  
  // Metric data
  // -----------------
  
  /// Chosen horizon type
   int type_hor;

   ///Chose field set type
   int field_set;

  /// Lapse function 
  Scalar lapse ; 
	
  // Conformal factor
  Scalar conf_fact;

  /// Shift vector
  Vector shift ;
	
  /** Deviation tensor( non-flat part of the conformal 3-metric on
   * the slice; see Bonazzola et al. (2004)).
   */
  Sym_tensor hij ;

  /**  Rescaled tracefree extrinsic curvature tensor: rescaled same way as Cordero et al.( 2009).
  *    In the New FCF scheme, this member will represent the TT part of the extrinsic curvature only. 
   */	
  Sym_tensor hatA ; 

  /** Longitudinal part of the extrinsic curvature. Set to zero if we are in the original FCF scheme.
   */

  Vector Xx ;
 
      
  // Derived data : 
  // ------------
 protected:

    
  /// Computation of the spheroid associated with the black hole horizon

  mutable Spheroid* p_hor ;
  ///   Computation of the ADM mass of the BH spacetime.
  mutable double* p_adm_mass ; 

  /** Computation of the Komar angular momentum w.r.t. assumed
   *rotational symmetry
   */
  mutable double* p_komar_angmom ;

  /**Computation of the Virial residual, as difference at infinity
   *between the ADM mass and the Komar integral associated to the mass.
   */
  mutable double* p_virial_residue ;
	

	


  // Constructors - Destructor
  // -------------------------
 public:
    
  /** Standard constructor. 
   * 
   * @param mp_i Mapping on which the black hole slice will be defined
   * @param hor_type Type of horizon used for excision:
   *                    - 0. General
   *                    - 1. Isolated horizon-like
   *                    - 2. Dynamical horizon-like
   *                    - 3. Inner trapped surface
   * @param scheme_type Type of scheme use for metric field computation
   *                    - 1. Original FCF scheme (Bonazzola et al., 2004)
   *                    - 2. New FCF scheme (Cordero et al., 2009)
   */


  Excised_slice(const Map& mp_i, int hor_type, int scheme_type) ;			
	
	
  Excised_slice(const Excised_slice& ) ;    ///< Copy constructor

  /** Constructor from a file (see \c sauve(FILE* )). 
   * 
   * @param mp_i Mapping on which the black hole slice will be defined
   * @param hor_type Type of horizon used for excision:
   *                    - 0. General
   *                    - 1. Isolated horizon-like
   *                    - 2. Dynamical horizon-like
   *                    - 3. Inner trapped surface
   * @param scheme_type Type of scheme use for metric field computation
   *                    - 1. Original FCF scheme (Bonazzola et al., 2004)
   *                    - 2. New FCF scheme (Cordero et al., 2009)
   * @param fich	input file (must have been created by the function
   *	\c sauve)
   */

  Excised_slice(const Map& mp_i, int hor_type, int scheme_type, FILE* fich) ;    		

  virtual ~Excised_slice() ;			///< Destructor

	
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
  /// Assignment to another \c Excised_slice
  void operator=(const Excised_slice&) ;	

  /** If hor_type=1, computes a quasi-stationary 3-slice from the chosen parameters.   
   *  @param precis [input] threshold in the relative difference between 
   *	the radial shift for two consecutive steps
   *	to stop the iterative procedure (default value: 1.e-11)
   * @param Omega_i: Rotation parameter for the isolated horizon slice. The principal Killing vector associated
   *  to this quantity will be by default the $e_{\varphi}$ vector.
   * @param NorKappa_i: boolean to tell whether boundary condition is put
   *  on the surface gravity (true) or the lapse (false)
   * @param NoK_i: The corresponding Dirichlet value for this quantity.
   * @param isCF: Tells whether the calculation is using the conformal flatness approximation (true) or not (false).
   * @param relax relaxation coefficient for all metric fields (0: no relaxation)
   * @param mer_max: Maximum number of iterations allowed in the main loop 
   * @param mer_max: Maximum number of iterations allowed in the possible secondary loop for conformal geometry computation.
   * @param isvoid: indicates whether or not the spacetime is vacuum.
   * This is useless for now since matter terms are not handled in the class.
   */

  void compute_stat_metric(double precis, double Omega_i, bool NorKappa_i, Scalar NoK_i, bool isCF_i=true, double relax=0.2, int mer_max = 2000, int mer_max2=200, bool isvoid = true) ; 

  /** If hor_type=1, computes the rhs of hyperbolic equation for conformal metric
   *assuming stationarity; 
   *WARNING; up to now, we are only able to handle void spacetimes.
   */

  void secmembre_kerr(Sym_tensor& source_hh);
 
 
  // Accessors
  // ---------
 public:
  /// Returns the mapping
  const Map& get_mp() const {return mp; } ; 

  ///Returns the type of horizon that performs the excision
  int get_type_hor() const {return type_hor; };
  
  ///Returns the field set chosen for the data

  int get_field_set() const {return field_set;};

  /// Returns the lapse function \e N
  const Scalar& get_lapse() const {return lapse;} ;
	
  /// Returns the conformal factor 
  const Scalar& get_conf_fact() const {return conf_fact;};

  /// Returns the shift vector \f$\beta^i\f$.
  const Vector& get_shift() const {return shift;} ;

  /// Returns the deviation tensor \f$ h^{ij} \f$.
  const Sym_tensor& get_hij() const {return hij;} ;

  /// Returns the rescaled tracefree extrinsic curvature \f$\hat{A}^{ij}\f$ (or its TT part if scheme_type=2)
  const Sym_tensor& get_hatA() const {return hatA;} ;
  
  /// Return the longitudinal part of Einstein Equations if the scheme is appropriate.
  const Vector& get_Xx() const {
    if (field_set==2) return Xx; 
    else if (field_set==1) cout <<"Error: the scheme used is the original FCF; this variable is irrelevant" << endl;
    else cout <<"error in the scheme definition; please check the class consistency" << endl;}

  ///Sets the horizon type
  int set_type_hor() {del_deriv() ; return type_hor ; } ;
  
  /// Sets the lapse
  Scalar& set_lapse() {del_deriv() ; return lapse ; } ;
  
  ///Sets the conformal factor
  Scalar& set_conf_fact() {del_deriv() ; return conf_fact ; } ;
  
  ///Sets the shift vector
  Vector& set_shift() {del_deriv() ; return shift ; } ;
  
  ///Sets the deviation tensor.
  Sym_tensor& set_hij() {del_deriv() ; return hij ; } ;
  
  /// Sets the rescaled tracefree extrinsic curvature \f$\hat{A}^{ij}\f$ (or its TT part if scheme_type=2)
  Sym_tensor& set_hatA() {del_deriv() ; return hatA ; } ;

  ///Sets the  longitudinal part of Einstein Equations if the scheme is apropriate
  Vector& set_Xx() {del_deriv() ; return Xx ; } ; // Include an assert here!!!!


  // Outputs
  // -------
 public:
  virtual void sauve(FILE* ) const ;	    ///< Save in a file

  void Einstein_errors();                     ///< Prints out errors in Einstein equations for the data obtained.

 
  // Global quantities
  // -----------------
 public:

      
  /** Spheroid at the excised boundary associated with the black hole
   * MOTS on the slice. Set by default at the position \f$ r=1 \f$.
   */

  Spheroid hor() ;
 
  ///   Computation of the ADM mass of the BH spacetime.
  double adm_mass() ; 

  /** Computation of the Komar angular momentum w.r.t. assumed
   *rotational symmetry
   */
  double komar_angmom() ;

  /**Computation of the Virial residual, as difference at infinity
   *between the ADM mass and the Komar integral associated to the mass.
   */
  double virial_residue() ;	

};

}
#endif
