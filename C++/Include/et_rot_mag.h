/*
 *  Definition of Lorene class Et_rot_mag
 *
 */

/*
 *   Copyright (c) 2002 Emmanuel Marcq
 *   Copyright (c) 2002,2013 Jerome Novak
 *   Copyright (c) 2013 Deberati Chatterjee
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


#ifndef __ET_ROT_MAG_H_ 
#define __ET_ROT_MAG_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.22  2014/05/13 10:06:12  j_novak
 * Change of magnetic units, to make the Lorene unit system coherent. Magnetic field is now expressed in Lorene units. Improvement on the comments on units.
 *
 * Revision 1.21  2014/04/29 13:46:06  j_novak
 * Addition of switches 'use_B_in_eos' and 'include_magnetisation' to control the model.
 *
 * Revision 1.20  2014/04/28 12:48:12  j_novak
 * Minor modifications.
 *
 * Revision 1.19  2013/12/13 16:36:51  j_novak
 * Addition and computation of magnetisation terms in the Einstein equations.
 *
 * Revision 1.18  2013/11/25 13:52:11  j_novak
 * New class Et_magnetisation to include magnetization terms in the stress energy tensor.
 *
 * Revision 1.17  2012/08/12 17:48:36  p_cerda
 * Magnetstar: New classes for magnetstar. Allowing for non-equatorial symmetry in Etoile et al. Adding B_phi in Et_rot_mag.
 *
 * Revision 1.16  2011/10/06 14:55:36  j_novak
 * equation_of_state() is now virtual to be able to call to the magnetized
 * Eos_mag.
 *
 * Revision 1.15  2005/06/02 11:35:27  j_novak
 * Added members for sving to a file and reading from it.
 *
 * Revision 1.14  2004/03/22 13:12:41  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.13  2002/10/11 11:47:35  j_novak
 * Et_rot_mag::MHD_comput is now virtual.
 * Use of standard constructor for Tenseur mtmp in Et_rot_mag::equilibrium_mag
 *
 * Revision 1.12  2002/10/09 07:54:29  j_novak
 * Et_rot_bifluid and Et_rot_mag inheritate virtually from Etoile_rot
 *
 * Revision 1.11  2002/08/02 15:07:41  j_novak
 * Member function determinant has been added to the class Metrique.
 * A better handling of spectral bases is now implemented for the class Tenseur.
 *
 * Revision 1.10  2002/06/05 15:15:59  j_novak
 * The case of non-adapted mapping is treated.
 * parmag.d and parrot.d have been merged.
 *
 * Revision 1.9  2002/06/03 13:23:16  j_novak
 * The case when the mapping is not adapted is now treated
 *
 * Revision 1.8  2002/06/03 13:00:45  e_marcq
 *
 * conduc parameter read in parmag.d
 *
 * Revision 1.7  2002/05/30 16:06:30  j_novak
 * added the right et_rot_mag.h
 *
 * Revision 1.6  2002/05/20 08:27:59  j_novak
 * *** empty log message ***
 *
 * Revision 1.5  2002/05/17 15:08:01  e_marcq
 *
 * Rotation progressive plug-in, units corrected, Q and a_j new member data
 *
 * Revision 1.4  2002/05/15 09:53:59  j_novak
 * First operational version
 *
 * Revision 1.3  2002/05/14 13:38:36  e_marcq
 *
 *
 * Unit update, new outputs
 *
 * Revision 1.1  2002/05/10 09:26:51  j_novak
 * Added new class Et_rot_mag for magnetized rotating neutron stars (under development)
 *
 *
 * $Header$
 *
 */

// Headers Lorene

#include "etoile.h"
#include "tensor.h"

// Local prototype (for determining the surface)
Cmp prolonge_c1(const Cmp& uu, const int nzet) ;

/**
 * Class for magnetized (isolator or perfect conductor), 
 * rigidly rotating stars. \ingroup (star)
 *
 * This is a child class of \c Etoile_rot , with the same metric
 * and overloaded member functions. Triaxial pertubrations are not 
 * operational.
 *
 */
class Et_rot_mag : public Etoile_rot {
  
  // Data : 
  // -----
 protected:

  ///t-component of the elecctromagnetic potential 1-form, divided by \f$\mu_0\f$.
  Cmp A_t ; 
  /**
   * \f$\varphi\f$-component of the electromagnetic potential 1-form
   * divided by \f$\mu_0\f$.
   */
  Cmp A_phi; 

  Cmp B_phi; ///< \f$\varphi\f$-component of the magnetic field
  Cmp j_t; ///< t-component of the current 4-vector
  Cmp j_phi; ///< \f$\varphi\f$-component of the current 4-vector

  Tenseur E_em; ///< electromagnetic energy density in the Eulerian frame

  /**
   * \f$\varphi\f$ component of the electromagnetic momentum density 3-vector,
   * as measured in the Eulerian frame.
   */
  Tenseur Jp_em; 

  ///rr component of the electromagnetic stress 3-tensor, as measured in the Eulerian frame. (not used and set to 0, should be supressed)
  Tenseur Srr_em;

  ///\f$\varphi \varphi\f$ component of the electromagnetic stress 3-tensor, as measured in the Eulerian frame. 
  Tenseur Spp_em; 

  /**
   * In the case of a perfect conductor, the requated baryonic charge.
   * For an isolator, the charge/baryon.
   */
  double Q ;
  double a_j ; ///<Amplitude of the curent/charge function
  int conduc ; ///<Flag: conduc=0->isolator, 1->perfect conductor

  // Constructors - Destructor
  // -------------------------
 public:

  /// Standard constructor
  Et_rot_mag(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i, 
	     const int cond); 


  Et_rot_mag(const Et_rot_mag& ) ;       ///< Copy constructor


  /** Constructor from a file (see \c sauve(FILE*) ). 
   * 
   * @param mp_i Mapping on which the star will be defined
   * @param eos_i Equation of state of the stellar matter
   * @param fich	input file (must have been created by the function
   *	\c sauve )
   * @param withbphi flag to create classes with toroidal field
   */
  Et_rot_mag(Map& mp_i, const Eos& eos_i, FILE* fich, int withbphi=0) ;    	
 		

  virtual ~Et_rot_mag() ;			///< Destructor
 

  // Memory management
  // -----------------
 protected:

  /// Deletes all the derived quantities
  virtual void del_deriv() const ; 
	
  /// Sets to \c 0x0  all the pointers on derived quantities
  virtual void set_der_0x0() const ; 

  /** Sets to \c ETATNONDEF  (undefined state) the hydrodynamical 
   *  quantities relative to the Eulerian observer.
   */
  virtual void del_hydro_euler() ; 
	

  // Mutators / assignment
  // ---------------------
 public:

  /// Assignment to another Et_rot_mag
  void operator=(const Et_rot_mag&) ;	

  /* /\** Computes the proper baryon and energy density, as well as */
  /*  *  pressure from the enthalpy. */
  /*  *\/ */
  /* virtual void equation_of_state() ;  */
	
  // Accessors
  // ---------
 public:
  /// Tells if the star is made of conducting or isolating material
  bool is_conduct() const {return (conduc==1) ;} ;
  /**
   * Returns the t component of the electromagnetic potential, 
   * divided by \f$\mu_0\f$.
   */
  const Cmp& get_At() const {return A_t ; } ; 
  /**
   * Returns the \f$\varphi\f$ component of the electromagnetic potential
   * divided by \f$\mu_0\f$.
   */
  const Cmp& get_Aphi() const {return A_phi ;} ;

  ///Returns the \f$\varphi\f$ component of the magnetic field
  const Cmp& get_Bphi() const {return B_phi ;} ;
  ///Returns the t component of the current 4-vector
  const Cmp& get_jt() const {return j_t ; } ;
  ///Returns the \f$\varphi\f$ component of the current 4-vector
  const Cmp& get_jphi() const {return j_phi ;} ;
  ///Returns the electromagnetic energy density in the Eulerian frame
  const Tenseur& get_Eem() const {return E_em ; } ;

  /** Returns the \f$\varphi\f$-component of the electromagnetic momentum
   * density 3-vector, as measured in the Eulerian frame.
   */
  const Tenseur& get_Jpem() const {return Jp_em ;} ;

  /** Returns the rr-component of the electromagnetic stress 3-tensor, 
   * as measured in the Eulerian frame. (not used and always equal to 0, 
   * should be supressed)
   */
  const Tenseur& get_Srrem() const {return Srr_em ; } ;

  /** Returns the \f$\varphi \varphi\f$ component of the electromagnetic 
   * stress 3-tensor, as measured in the Eulerian frame. 
   */
  const Tenseur& get_Sppem() const {return Spp_em ;} ;

  /**
   * Returns the requested electric charge in the case of a perfect conductor
   * and the charge/baryon for an isolator.
   */
  double get_Q() const {return Q ;} ;
  ///Returns the amplitude of the current/charge function
  double get_a_j() const {return a_j ;} ;

  // Outputs
  // -------
 public:
  virtual void sauve(FILE* ) const ;	    ///< Save in a file
   
  /// Operator >> (virtual function called by the operator <<). 
  virtual ostream& operator>>(ostream& ) const ;    

  // Global quantities
  // -----------------
 public:
  /// Computes the electric field spherical components in Lorene's units
  Tenseur Elec() const ; 
  /// Computes the magnetic field spherical components in Lorene's units
  Tenseur Magn() const ; 
  /// Computes the electromagnetic part of the stress-energy tensor
  virtual void MHD_comput() ; 
  virtual double mass_g() const ;	    ///< Gravitational mass
  virtual double angu_mom() const ;  ///< Angular momentum
  virtual double grv2() const ;	///< Error on the virial identity GRV2
  virtual double tsw() const ; ///< Ratio T/W
  double MagMom() const ; ///< Magnetic Momentum \f$\cal M\f$ in SI units
  /// Computed charge deduced from the asymptotic behaviour of At [SI units].
  double Q_comput() const; 

  /** Computed charge from the integration of charge density over the star
   * (i.e. without surface charge) [SI units].
   */
  double Q_int() const; 

  /// Gyromagnetic ratio \f$\sigma = \frac{2{\cal M}M}{QJ}\f$.
  double GyroMag() const ; 

  /** Error on the virial identity GRV3.
   *  The error is computed as the integral defined
   *  by Eq. (43) of [Gourgoulhon and Bonazzola, 
   *  \a Class. \a Quantum \a Grav. \b 11 , 443 (1994)] divided by
   *  the integral of the matter terms.
   * 
   *  @param ost output stream to give details of the computation;
   *		if set to 0x0 [default value], no details will be
   *		given.
   *   
   */
  virtual double grv3(ostream* ost = 0x0) const ;	

  /** Quadrupole moment.
   *  The quadrupole moment \e Q  is defined according to Eq. (7) of
   *  [Salgado, Bonazzola, Gourgoulhon and Haensel, \a Astron. \a Astrophys.
   *   \b 291 , 155 (1994)]. At the Newtonian limit it is related to
   *  the component \f${\bar I}_{zz}\f$ of the MTW (1973) reduced quadrupole 
   *  moment \f${\bar I}_{ij}\f$ by: \f$Q = -3/2 {\bar I}_{zz}\f$. 
   *  Note that \e Q  is the negative of the quadrupole moment defined 
   *  by Laarakkers and Poisson, \a Astrophys. \a J. \b 512 , 282 (1999).
   */
  virtual double mom_quad() const ;	

  // Computational routines
  // ----------------------
 public: 
  /** Computes the electromagnetic quantities solving the Maxwell
   *  equations (6) and (7) of [Bocquet, Bonazzola, Gourgoulhon and
   *  Novak, \a Astron. \a Astrophys. \b 301 , 757 (1995)]. In the case 
   *  of a perfect conductor, le electromagnetic potential may have
   *  a discontinuous derivative across star's surface.
   *
   *  @param conduc [input] flag: 0 for an isolator, 1 for a perfect 
   *                        conductor
   *  @param adapt_flag [input] flag: if 0 the mapping is NOT adapted
   *                             to star's surface
   *  @param f_j [input] current or charge coupling function 
   *                      (see Bocquet et al. 1995).
   *  @param par_poisson_At [input] parameters for controlling the 
   *                                  solution of the Poisson equation
   *                                  for At potential (see file
   *                                  et_rot_mag_equil.C)
   *  @param par_poisson_Avect [input] parameters for controlling the 
   *                                  solution of vector Poisson equation
   *                                  for magnetic potential (see file
   *                                  et_rot_mag_equil.C)
   */
  void magnet_comput(const int adapt_flag,
			     Cmp (*f_j)(const Cmp& x, const double),
		      Param& par_poisson_At, Param& par_poisson_Avect) ;


  /** Computes the electromagnetic quantities solving the Maxwell
   *  equations (6) and (7) of [Bocquet, Bonazzola, Gourgoulhon and
   *  Novak, \a Astron. \a Astrophys. \b 301 , 757 (1995)]. In the case 
   *  of a perfect conductor, le electromagnetic potential may have
   *  a discontinuous derivative across star's surface.
   *
   *  @param adapt_flag [input] flag: if 0 the mapping is NOT adapted
   *                             to star's surface
   *  @param initial_j [input] flag: initial current for the iteration:
   *                                 0= no current, 1=dipolar-like current
   *                                 , 2= quadrupolar-like current
   *  @param a_j0 [input] amplitude of the non-force free current
   *  @param f_j [input] current coupling function (non-FF part)
   *                      (see Bocquet et al. 1995).
   *  @param b_j0 [input] amplitude of the force free current
   *  @param g_j [input] current coupling function (FF-part)
   *  @param N_j [input] current coupling function (FF-part)
   *  @param par_poisson_At [input] parameters for controlling the 
   *                                  solution of the Poisson equation
   *                                  for At potential (see file
   *                                  et_rot_mag_equil.C)
   *  @param par_poisson_Avect [input] parameters for controlling the 
   *                                  solution of vector Poisson equation
   *                                  for magnetic potential (see file
   *                                  et_rot_mag_equil.C)
   */
  virtual void magnet_comput_plus(const int adapt_flag, const int initial_j,
      const Tbl an_j,
      Cmp (*f_j)(const Cmp& x, const Tbl),
      const Tbl bn_j,
      Cmp (*g_j)(const Cmp& x, const Tbl),
      Cmp (*N_j)(const Cmp& x, const Tbl),
      Param& par_poisson_At, Param& par_poisson_Avect) ;
	
  /** Computes an equilibrium configuration.
   *  
   *  @param ent_c  [input] Central enthalpy 
   *  @param omega0  [input] Requested angular velocity 
   *			     (if \c fact_omega=1. )
   *  @param fact_omega [input] 1.01 = search for the Keplerian frequency,
   *			      1. = otherwise.
   *  @param nzadapt  [input] Number of (inner) domains where the mapping 
   *			    adaptation to an iso-enthalpy surface
   *			    should be performed
   *  @param ent_limit [input] 1-D \c Tbl  of dimension \c nzet  which
   *				defines the enthalpy at the outer boundary
   *				of each domain
   *  @param icontrol [input] Set of integer parameters (stored as a
   *			    1-D \c Itbl  of size 8) to control the 
   *			    iteration: 
   *	\li \c icontrol(0) = mer_max  : maximum number of steps 
   *	\li \c icontrol(1) = mer_rot  : step at which the rotation is 
   *				      switched on 
   *	\li \c icontrol(2) = mer_change_omega  : step at which the rotation
   *			  velocity is changed to reach the final one  
   *	\li \c icontrol(3) = mer_fix_omega  :  step at which the final
   *			    rotation velocity must have been reached  
   *	\li \c icontrol(4) = mer_mass  : the absolute value of 
   *			    \c mer_mass  is the step from which the 
   *			    baryon mass is forced to converge, 
   *			    by varying the central enthalpy 
   *			    (\c mer_mass > 0 ) or the angular 
   *			    velocity (\c mer_mass < 0 ) 
   *	\li \c icontrol(5) = mermax_poisson  : maximum number of steps in 
   *				\c Map_et::poisson  
   *	\li \c icontrol(6) = mer_triax  : step at which the 3-D 
   *				perturbation is switched on 
   *	\li \c icontrol(7) = delta_mer_kep  : number of steps
   *			    after \c mer_fix_omega  when \c omega 
   *			    starts to be increased by \c fact_omega 
   *			    to search for the Keplerian velocity
   *    \li \c icontrol(8) = mer_mag  : step at which the electromagnetic
   *                        part is switched on 
   *    \li \c icontrol(9) = mer_change_mag  : step at which the amplitude
   *                        of the current/charge coupling function is changed
   *                        to reach a_j0 or Q
   *	\li \c icontrol(10) = mer_fix_mag  :  step at which the final
   *			    current/charge amplitude a_j0 or Q must have 
   *                        been reached  
   *    \li \c icontrol(11) = conduc  : flag 0 -> isolator material, 
   *                        1 -> perfect conductor 
   * 	 
   *  @param control [input] Set of parameters (stored as a 
   *			    1-D \c Tbl  of size 7) to control the 
   *			    iteration: 
   *	\li \c control(0) = precis  : threshold on the enthalpy relative 
   *				change for ending the computation 
   *	\li \c control(1) = omega_ini  : initial angular velocity, 
   *			    switched on only if \c mer_rot < 0 , 
   *			    otherwise 0 is used  
   *	\li \c control(2) = relax  : relaxation factor in the main 
   *				   iteration  
   *	\li \c control(3) = relax_poisson  : relaxation factor in 
   *				   \c Map_et::poisson  
   *	\li \c control(4) = thres_adapt  :  threshold on dH/dr for 
   *			    freezing the adaptation of the mapping 
   *	\li \c control(5) = ampli_triax  :  relative amplitude of 
   *			    the 3-D perturbation 
   *	\li \c control(6) = precis_adapt  : precision for 
   *			    \c Map_et::adapt 
   *	\li \c control(7) = Q_ini  : initial charge (total for the perfect 
   *                      conductor, per baryon for an isolator)
   *	\li \c control(8) = a_j_ini  : initial amplitude for the coupling
   *                      function
   * 
   *
   *  @param mbar_wanted [input] Requested baryon mass (effective only 
   *				if \c mer_mass>mer_max )
   *  @param aexp_mass [input] Exponent for the increase factor of the 
   *			      central enthalpy to converge to the 
   *			      requested baryon mass
   *  @param diff [output]   1-D \c Tbl  of size 1 for the storage of 
   *			    some error indicators : 
   *	    \li \c diff(0)  : Relative change in the enthalpy field
   *			      between two successive steps 
   *  @param Q0 [input] Requested electric charge for the case of a
   *                    perfect conductor. Charge per baryon for the case
   *                    of an isolator.
   *  @param a_j0 [input] Amplitude for the current/charge coupling function
   *
   *  @param f_j [input] current or charge coupling function 
   *                      (see Bocquet et al. 1995).
   * 
   *  @param M_j [input] primitive (null for zero) of current/charge 
   *                      coupling function (see Bocquet et al. 1995)
   *                      used for the first integral of stationary motion.
   */
  virtual void equilibrium_mag(double ent_c, double omega0, double fact_omega, 
		       int nzadapt, const Tbl& ent_limit, const Itbl& icontrol, 
		       const Tbl& control, double mbar_wanted, double aexp_mass, 
		       Tbl& diff, const double Q0, const double a_j0, 
		       Cmp (*f_j)(const Cmp& x, const double), 
		       Cmp (*M_j)(const Cmp& x,const double));


  /** Computes an equilibrium configuration.
   *  
   *  @param ent_c  [input] Central enthalpy 
   *  @param omega0  [input] Requested angular velocity 
   *			     (if \c fact_omega=1. )
   *  @param fact_omega [input] 1.01 = search for the Keplerian frequency,
   *			      1. = otherwise.
   *  @param nzadapt  [input] Number of (inner) domains where the mapping 
   *			    adaptation to an iso-enthalpy surface
   *			    should be performed
   *  @param ent_limit [input] 1-D \c Tbl  of dimension \c nzet  which
   *				defines the enthalpy at the outer boundary
   *				of each domain
   *  @param icontrol [input] Set of integer parameters (stored as a
   *			    1-D \c Itbl  of size 8) to control the 
   *			    iteration: 
   *	\li \c icontrol(0) = mer_max  : maximum number of steps 
   *	\li \c icontrol(1) = mer_rot  : step at which the rotation is 
   *				      switched on 
   *	\li \c icontrol(2) = mer_change_omega  : step at which the rotation
   *			  velocity is changed to reach the final one  
   *	\li \c icontrol(3) = mer_fix_omega  :  step at which the final
   *			    rotation velocity must have been reached  
   *	\li \c icontrol(4) = mer_mass  : the absolute value of 
   *			    \c mer_mass  is the step from which the 
   *			    baryon mass is forced to converge, 
   *			    by varying the central enthalpy 
   *			    (\c mer_mass > 0 ) or the angular 
   *			    velocity (\c mer_mass < 0 ) 
   *	\li \c icontrol(5) = mermax_poisson  : maximum number of steps in 
   *				\c Map_et::poisson  
   *	\li \c icontrol(6) = mer_triax  : step at which the 3-D 
   *				perturbation is switched on 
   *	\li \c icontrol(7) = delta_mer_kep  : number of steps
   *			    after \c mer_fix_omega  when \c omega 
   *			    starts to be increased by \c fact_omega 
   *			    to search for the Keplerian velocity
   *    \li \c icontrol(8) = mer_mag  : step at which the electromagnetic
   *                        part is switched on 
   *    \li \c icontrol(9) = mer_change_mag  : step at which the amplitude
   *                        of the current/charge coupling function is changed
   *                        to reach a_j0 or Q
   *	\li \c icontrol(10) = mer_fix_mag  :  step at which the final
   *			    current/charge amplitude a_j0 or Q must have 
   *                        been reached  
   *    \li \c icontrol(11) = conduc  : flag 0 -> isolator material, 
   *                        1 -> perfect conductor 
   * 	 
   *  @param control [input] Set of parameters (stored as a 
   *			    1-D \c Tbl  of size 7) to control the 
   *			    iteration: 
   *	\li \c control(0) = precis  : threshold on the enthalpy relative 
   *				change for ending the computation 
   *	\li \c control(1) = omega_ini  : initial angular velocity, 
   *			    switched on only if \c mer_rot < 0 , 
   *			    otherwise 0 is used  
   *	\li \c control(2) = relax  : relaxation factor in the main 
   *				   iteration  
   *	\li \c control(3) = relax_poisson  : relaxation factor in 
   *				   \c Map_et::poisson  
   *	\li \c control(4) = thres_adapt  :  threshold on dH/dr for 
   *			    freezing the adaptation of the mapping 
   *	\li \c control(5) = ampli_triax  :  relative amplitude of 
   *			    the 3-D perturbation 
   *	\li \c control(6) = precis_adapt  : precision for 
   *			    \c Map_et::adapt 
   *	\li \c control(7) = Q_ini  : initial charge (total for the perfect 
   *                      conductor, per baryon for an isolator)
   *	\li \c control(8) = a_j_ini  : initial amplitude for the coupling
   *                      function
   * 
   *
   *  @param mbar_wanted [input] Requested baryon mass (effective only 
   *				if \c mer_mass>mer_max )
   *  @param aexp_mass [input] Exponent for the increase factor of the 
   *			      central enthalpy to converge to the 
   *			      requested baryon mass
   *  @param diff [output]   1-D \c Tbl  of size 1 for the storage of 
   *			    some error indicators : 
   *	    \li \c diff(0)  : Relative change in the enthalpy field
   *			      between two successive steps 
   *  @param Q0 [input] Requested electric charge for the case of a
   *                    perfect conductor. Charge per baryon for the case
   *                    of an isolator.
   *  @param a_j0 [input] Amplitude for the current/charge coupling function
   *
   *  @param f_j [input] current or charge coupling function 
   *                      (see Bocquet et al. 1995).
   * 
   *  @param M_j [input] primitive (null for zero) of current/charge 
   *                      coupling function (see Bocquet et al. 1995)
   *                      used for the first integral of stationary motion.
   */
      void equilibrium_mag_plus( const Itbl& icontrol, const Tbl& control, 
		       Tbl& diff,
                       const int initial_j,
                       const Tbl an_j, 
		       Cmp (*f_j)(const Cmp& x, const Tbl), 
		       Cmp (*M_j)(const Cmp& x,const Tbl),
                       const Tbl bn_j, 
		       Cmp (*g_j)(const Cmp& x, const Tbl), 
		       Cmp (*N_j)(const Cmp& x,const Tbl),
                       const double relax_mag);
};

class Et_magnetisation : public Et_rot_mag {

  // Data : 
  // -----
 protected:

  ///Flag : true if the value of the magnetic field is used in the Eos.
  bool use_B_in_eos ; 

  ///Flag : true if magnetisation terms are included in the equations.
  bool include_magnetisation ;

  Scalar xmag ; ///< The magnetisation scalar.

  Scalar E_I; ///< Interaction (magnetisation) energy density.

  ///Interaction momentum density 3-vector.
  Vector J_I; 

  ///Interaction stress 3-tensor.
  Sym_tensor Sij_I;

  // Constructors - Destructor
  // -------------------------
 public:

  /// Standard constructor
  Et_magnetisation(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i,
		   bool include_mag=true, bool use_B = true); 


  Et_magnetisation(const Et_magnetisation& ) ;       ///< Copy constructor


  /** Constructor from a file (see \c sauve(FILE*) ). 
   * 
   * @param mp_i Mapping on which the star will be defined
   * @param eos_i Equation of state of the stellar matter
   * @param fich	input file (must have been created by the function
   *	\c sauve )
   */
  Et_magnetisation(Map& mp_i, const Eos& eos_i, FILE* fich) ;    	

  virtual ~Et_magnetisation() ;			///< Destructor
 
  // Mutators / assignment
  // ---------------------
 public:

  /// Assignment to another Et_rot_mag
  void operator=(const Et_magnetisation&) ;	

  /** Computes the proper baryon and energy density, as well as
   *  pressure from the enthalpy.
   */
  virtual void equation_of_state() ; 

  // Accessors
  // ---------
 public:
  ///Public accessor to the \c use_B_in_eos flag.
  bool B_in_eos() const {return use_B_in_eos;};

  ///Public accessor to the \c include_magnetisation flag.
  bool use_magnetisation() const {return include_magnetisation;} ;

  ///Accessor to the magnetisation scalar field.
  const Scalar& get_magnetisation() const {return xmag;} ;

  ///Accessor to the interaction energy density.
  const Scalar& get_E_I() const {return E_I;} ;

  ///Accessor to the interaction momentum vector.
  const Vector& get_J_I() const {return J_I;} ;

  ///Accessor to the interaction stress tensor.
  const Sym_tensor& get_Sij_I() const {return Sij_I;} ; 
  
  // Outputs
  // -------
 public:
  virtual void sauve(FILE* ) const ;	    ///< Save in a file
   
  /// Operator >> (virtual function called by the operator <<). 
  virtual ostream& operator>>(ostream& ) const ;    

  // Computational routines
  // ----------------------
 public: 
  /** Computes the electromagnetic quantities solving the Maxwell
   *  equations (6) and (7) of [Bocquet, Bonazzola, Gourgoulhon and
   *  Novak, \a Astron. \a Astrophys. \b 301 , 757 (1995)]. In the case 
   *  of a perfect conductor, le electromagnetic potential may have
   *  a discontinuous derivative across star's surface.
   *
   *  @param conduc [input] flag: 0 for an isolator, 1 for a perfect 
   *                        conductor
   *  @param adapt_flag [input] flag: if 0 the mapping is NOT adapted
   *                             to star's surface
   *  @param f_j [input] current or charge coupling function 
   *                      (see Bocquet et al. 1995).
   *  @param par_poisson_At [input] parameters for controlling the 
   *                                  solution of the Poisson equation
   *                                  for At potential (see file
   *                                  et_rot_mag_equil.C)
   *  @param par_poisson_Avect [input] parameters for controlling the 
   *                                  solution of vector Poisson equation
   *                                  for magnetic potential (see file
   *                                  et_rot_mag_equil.C)
   */
  virtual void magnet_comput(const int adapt_flag,
			     Cmp (*f_j)(const Cmp& x, const double),
		      Param& par_poisson_At, Param& par_poisson_Avect) ;

  /// Computes the electromagnetic part of the stress-energy tensor
  virtual void MHD_comput() ; 

  /** Computes an equilibrium configuration.
   *  
   *  @param ent_c  [input] Central enthalpy 
   *  @param omega0  [input] Requested angular velocity 
   *			     (if \c fact_omega=1. )
   *  @param fact_omega [input] 1.01 = search for the Keplerian frequency,
   *			      1. = otherwise.
   *  @param nzadapt  [input] Number of (inner) domains where the mapping 
   *			    adaptation to an iso-enthalpy surface
   *			    should be performed
   *  @param ent_limit [input] 1-D \c Tbl  of dimension \c nzet  which
   *				defines the enthalpy at the outer boundary
   *				of each domain
   *  @param icontrol [input] Set of integer parameters (stored as a
   *			    1-D \c Itbl  of size 8) to control the 
   *			    iteration: 
   *	\li \c icontrol(0) = mer_max  : maximum number of steps 
   *	\li \c icontrol(1) = mer_rot  : step at which the rotation is 
   *				      switched on 
   *	\li \c icontrol(2) = mer_change_omega  : step at which the rotation
   *			  velocity is changed to reach the final one  
   *	\li \c icontrol(3) = mer_fix_omega  :  step at which the final
   *			    rotation velocity must have been reached  
   *	\li \c icontrol(4) = mer_mass  : the absolute value of 
   *			    \c mer_mass  is the step from which the 
   *			    baryon mass is forced to converge, 
   *			    by varying the central enthalpy 
   *			    (\c mer_mass > 0 ) or the angular 
   *			    velocity (\c mer_mass < 0 ) 
   *	\li \c icontrol(5) = mermax_poisson  : maximum number of steps in 
   *				\c Map_et::poisson  
   *	\li \c icontrol(6) = mer_triax  : step at which the 3-D 
   *				perturbation is switched on 
   *	\li \c icontrol(7) = delta_mer_kep  : number of steps
   *			    after \c mer_fix_omega  when \c omega 
   *			    starts to be increased by \c fact_omega 
   *			    to search for the Keplerian velocity
   *    \li \c icontrol(8) = mer_mag  : step at which the electromagnetic
   *                        part is switched on 
   *    \li \c icontrol(9) = mer_change_mag  : step at which the amplitude
   *                        of the current/charge coupling function is changed
   *                        to reach a_j0 or Q
   *	\li \c icontrol(10) = mer_fix_mag  :  step at which the final
   *			    current/charge amplitude a_j0 or Q must have 
   *                        been reached  
   *    \li \c icontrol(11) = conduc  : flag 0 -> isolator material, 
   *                        1 -> perfect conductor 
   * 	 
   *  @param control [input] Set of parameters (stored as a 
   *			    1-D \c Tbl  of size 7) to control the 
   *			    iteration: 
   *	\li \c control(0) = precis  : threshold on the enthalpy relative 
   *				change for ending the computation 
   *	\li \c control(1) = omega_ini  : initial angular velocity, 
   *			    switched on only if \c mer_rot < 0 , 
   *			    otherwise 0 is used  
   *	\li \c control(2) = relax  : relaxation factor in the main 
   *				   iteration  
   *	\li \c control(3) = relax_poisson  : relaxation factor in 
   *				   \c Map_et::poisson  
   *	\li \c control(4) = thres_adapt  :  threshold on dH/dr for 
   *			    freezing the adaptation of the mapping 
   *	\li \c control(5) = ampli_triax  :  relative amplitude of 
   *			    the 3-D perturbation 
   *	\li \c control(6) = precis_adapt  : precision for 
   *			    \c Map_et::adapt 
   *	\li \c control(7) = Q_ini  : initial charge (total for the perfect 
   *                      conductor, per baryon for an isolator)
   *	\li \c control(8) = a_j_ini  : initial amplitude for the coupling
   *                      function
   * 
   *
   *  @param mbar_wanted [input] Requested baryon mass (effective only 
   *				if \c mer_mass>mer_max )
   *  @param aexp_mass [input] Exponent for the increase factor of the 
   *			      central enthalpy to converge to the 
   *			      requested baryon mass
   *  @param diff [output]   1-D \c Tbl  of size 1 for the storage of 
   *			    some error indicators : 
   *	    \li \c diff(0)  : Relative change in the enthalpy field
   *			      between two successive steps 
   *  @param Q0 [input] Requested electric charge for the case of a
   *                    perfect conductor. Charge per baryon for the case
   *                    of an isolator.
   *  @param a_j0 [input] Amplitude for the current/charge coupling function
   *
   *  @param f_j [input] current or charge coupling function 
   *                      (see Bocquet et al. 1995).
   * 
   *  @param M_j [input] primitive (null for zero) of current/charge 
   *                      coupling function (see Bocquet et al. 1995)
   *                      used for the first integral of stationary motion.
   */
  virtual void equilibrium_mag(double ent_c, double omega0, double fact_omega, 
		       int nzadapt, const Tbl& ent_limit, const Itbl& icontrol, 
		       const Tbl& control, double mbar_wanted, double aexp_mass, 
		       Tbl& diff, const double Q0, const double a_j0, 
		       Cmp (*f_j)(const Cmp& x, const double), 
		       Cmp (*M_j)(const Cmp& x,const double));

};

#endif

