/*
 *  Definition of Lorene class Et_rot_mag
 *
 */

/*
 *   Copyright (c) 2002 Emmanuel Marcq
 *   Copyright (c) 2002 Jerome Novak
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

// Local prototype (for determining the surface)
Cmp prolonge_c1(const Cmp& uu, const int nzet) ;

/**
 * Class for magnetized (isolator or perfect conductor), 
 * rigidly rotating stars.
 *
 * This is a child class of {\tt Etoile\_rot}, with the same metric
 * and overloaded member functions. Triaxial pertubrations are not 
 * operational.
 *
 * @version #$Id$#
 */
class Et_rot_mag : virtual public Etoile_rot {
  
  // Data : 
  // -----
 protected:

  Cmp A_t ; /// t-component of the elecctromagnetic potential 1-form
  Cmp A_phi; /// $\varphi$-component of the electromagnetic potential 1-form
  Cmp j_t; /// t-component of the current 4-vector
  Cmp j_phi; /// $\varphi$-component of the current 4-vector

  Tenseur E_em; /// electromagnetic energy density in the Eulerian frame

  ///$\varphi$ component of the electromagnetic momentum density 3-vector (as measured in the Eulerian frame).
  Tenseur Jp_em; 

  ///rr component of the electromagnetic stress 3-tensor, as measured in the Eulerian frame. (not used and set to 0, should be supressed)
  Tenseur Srr_em;

  ///$\varphi \varphi$ component of the electromagnetic stress 3-tensor, as measured in the Eulerian frame. 
  Tenseur Spp_em; 

  /**
   * In the case of a perfect conductor, the requated baryonic charge.\\
   * For an isolator, the charge/baryon.
   */
  double Q ;
  double a_j ; ///Amplitude of the curent/charge function
  int conduc ; ///Flag: conduc=0->isolator, 1->perfect conductor

  // Constructors - Destructor
  // -------------------------
 public:

  /// Standard constructor
  Et_rot_mag(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i, 
	     const int cond); 


  Et_rot_mag(const Et_rot_mag& ) ;       /// Copy constructor

  virtual ~Et_rot_mag() ;			/// Destructor
 

  // Memory management
  // -----------------
 protected:

  /// Deletes all the derived quantities
  virtual void del_deriv() const ; 
	
  /// Sets to {\tt 0x0} all the pointers on derived quantities
  virtual void set_der_0x0() const ; 

  /** Sets to {\tt ETATNONDEF} (undefined state) the hydrodynamical 
   *  quantities relative to the Eulerian observer.
   */
  virtual void del_hydro_euler() ; 
	

  // Mutators / assignment
  // ---------------------
 public:

  /// Assignment to another Et_rot_mag
  void operator=(const Et_rot_mag&) ;	

  // Accessors
  // ---------
 public:
  /// Tells if the star is made of conducting or isolating material
  bool is_conduct() const {return (conduc==1) ;} ;
  ///Returns the t component of the electromagnetic potential
  const Cmp& get_At() const {return A_t ; } ; 
  ///Returns the $\varphi$ component of the electromagnetic potential
  const Cmp& get_Aphi() const {return A_phi ;} ;
  ///Returns the t component of the current 4-vector
  const Cmp& get_jt() const {return j_t ; } ;
  ///Returns the $\varphi$ component of the current 4-vector
  const Cmp& get_jphi() const {return j_phi ;} ;
  ///Returns the electromagnetic energy density in the Eulerian frame
  const Tenseur& get_Eem() const {return E_em ; } ;

  /** Returns the $\varphi$-component of the electromagnetic momentum
   * density 3-vector, as measured in the Eulerian frame.
   */
  const Tenseur& get_Jpem() const {return Jp_em ;} ;

  /** Returns the rr-component of the electromagnetic stress 3-tensor, 
   * as measured in the Eulerian frame. (not used and always equal to 0, 
   * should be supressed)
   */
  const Tenseur& get_Srrem() const {return Srr_em ; } ;

  /** Returns the $\varphi \varphi$ component of the electromagnetic 
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
    
  /// Operator >> (virtual function called by the operator <<). 
  virtual ostream& operator>>(ostream& ) const ;    

  // Global quantities
  // -----------------
 public:

  Tenseur Elec() const ; /// Computes the electric field spherical components
  Tenseur Magn() const ; /// Computes the magnetic field spherical components
  /// Computes the electromagnetic part of the stress-energy tensor
  void MHD_comput() ; 
  virtual double mass_g() const ;	    /// Gravitational mass
  virtual double angu_mom() const ;  /// Angular momentum
  virtual double grv2() const ;	/// Error on the virial identity GRV2
  virtual double tsw() const ; /// Ratio T/W
  double MagMom() const ; /// Magnetic Momentum $\cal M$
  /// Computed charge deduced from the asymptotic behaviour of At.
  double Q_comput() const; 

  /** Computed charge from the integration of charge density over the star
   * (i.e. without surface charge)
   */
  double Q_int() const; 

  /// Gyromagnetic ratio $\sigma = \frac{2{\cal M}M}{QJ}$.
  double GyroMag() const ; 

  /** Error on the virial identity GRV3.
   *  The error is computed as the integral defined
   *  by Eq. (43) of [Gourgoulhon and Bonazzola, 
   *  Class. Quantum Grav. {\bf 11}, 443 (1994)] divided by
   *  the integral of the matter terms.
   * 
   *  @param ost output stream to give details of the computation;
   *		if set to 0x0 [default value], no details will be
   *		given.
   *   
   */
  virtual double grv3(ostream* ost = 0x0) const ;	

  /** Quadrupole moment.
   *  The quadrupole moment {\it Q} is defined according to Eq. (7) of
   *  [Salgado, Bonazzola, Gourgoulhon and Haensel, Astron. Astrophys.
   *   {\bf 291}, 155 (1994)]. At the Newtonian limit it is related to
   *  the component ${\bar I}_{zz}$ of the MTW (1973) reduced quadrupole 
   *  moment ${\bar I}_{ij}$ by: $Q = -3/2 {\bar I}_{zz}$. 
   *  Note that {\it Q} is the negative of the quadrupole moment defined 
   *  by Laarakkers and Poisson, Astrophys. J. {\bf 512}, 282 (1999).
   */
  virtual double mom_quad() const ;	

  // Computational routines
  // ----------------------
 public: 
  /** Computes the electromagnetic quantities solving the Maxwell
   *  equations (6) and (7) of [Bocquet, Bonazzola, Gourgoulhon and
   *  Novak, Astron. Astrophys. {\bf 301}, 757 (1995)]. In the case 
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
   *                                  et\_rot\_mag\_equil.C)
   *  @param par_poisson_Avect [input] parameters for controlling the 
   *                                  solution of vector Poisson equation
   *                                  for magnetic potential (see file
   *                                  et\_rot\_mag\_equil.C)
   */
  virtual void magnet_comput(const int adapt_flag,
			     Cmp (*f_j)(const Cmp& x, const double),
		      Param& par_poisson_At, Param& par_poisson_Avect) ;
	
  /** Computes an equilibrium configuration.
   *  
   *  @param ent_c  [input] Central enthalpy 
   *  @param omega0  [input] Requested angular velocity 
   *			     (if {\tt fact\_omega=1.})
   *  @param fact_omega [input] 1.01 = search for the Keplerian frequency,
   *			      1. = otherwise.
   *  @param nzadapt  [input] Number of (inner) domains where the mapping 
   *			    adaptation to an iso-enthalpy surface
   *			    should be performed
   *  @param ent_limit [input] 1-D {\tt Tbl} of dimension {\tt nzet} which
   *				defines the enthalpy at the outer boundary
   *				of each domain
   *  @param icontrol [input] Set of integer parameters (stored as a
   *			    1-D {\tt Itbl} of size 8) to control the 
   *			    iteration: \\
   *	{\tt icontrol(0) = mer\_max} : maximum number of steps \\
   *	{\tt icontrol(1) = mer\_rot} : step at which the rotation is 
   *				      switched on \\
   *	{\tt icontrol(2) = mer\_change\_omega} : step at which the rotation
   *			  velocity is changed to reach the final one  \\
   *	{\tt icontrol(3) = mer\_fix\_omega} :  step at which the final
   *			    rotation velocity must have been reached  \\
   *	{\tt icontrol(4) = mer\_mass} : the absolute value of 
   *			    {\tt mer\_mass} is the step from which the 
   *			    baryon mass is forced to converge, 
   *			    by varying the central enthalpy 
   *			    ({\tt mer\_mass > 0}) or the angular 
   *			    velocity ({\tt mer\_mass < 0}) \\
   *	{\tt icontrol(5) = mermax\_poisson} : maximum number of steps in 
   *				{\tt Map\_et::poisson} \\
   *	{\tt icontrol(6) = mer\_triax} : step at which the 3-D 
   *				perturbation is switched on \\
   *	{\tt icontrol(7) = delta\_mer\_kep} : number of steps
   *			    after {\tt mer\_fix\_omega} when {\tt omega}
   *			    starts to be increased by {\tt fact\_omega}
   *			    to search for the Keplerian velocity\\
   *    {\tt icontrol(8) = mer\_mag} : step at which the electromagnetic
   *                        part is switched on \\
   *    {\tt icontrol(9) = mer\_change\_mag} : step at which the amplitude
   *                        of the current/charge coupling function is changed
   *                        to reach a_j0 or Q\\
   *	{\tt icontrol(10) = mer\_fix\_mag} :  step at which the final
   *			    current/charge amplitude a_j0 or Q must have 
   *                        been reached  \\
   *    {\tt icontrol(11) = conduc} : flag 0 -> isolator material, 
   *                        1 -> perfect conductor \\
   * 	 
   *  @param control [input] Set of parameters (stored as a 
   *			    1-D {\tt Tbl} of size 7) to control the 
   *			    iteration: \\
   *	{\tt control(0) = precis} : threshold on the enthalpy relative 
   *				change for ending the computation \\
   *	{\tt control(1) = omega\_ini} : initial angular velocity, 
   *			    switched on only if {\tt mer\_rot < 0}, 
   *			    otherwise 0 is used \\ 
   *	{\tt control(2) = relax} : relaxation factor in the main 
   *				   iteration \\ 
   *	{\tt control(3) = relax\_poisson} : relaxation factor in 
   *				   {\tt Map\_et::poisson}\\ 
   *	{\tt control(4) = thres\_adapt} :  threshold on dH/dr for 
   *			    freezing the adaptation of the mapping \\
   *	{\tt control(5) = ampli\_triax} :  relative amplitude of 
   *			    the 3-D perturbation \\
   *	{\tt control(6) = precis\_adapt} : precision for 
   *			    {\tt Map\_et::adapt}
   *	{\tt control(7) = Q\_ini} : initial charge (total for the perfect 
   *                      conductor, per baryon for an isolator)\\
   *	{\tt control(8) = a\_j\_ini} : initial amplitude for the coupling
   *                      function\\
   * 
   *
   *  @param mbar_wanted [input] Requested baryon mass (effective only 
   *				if {\tt mer\_mass > mer\_max})
   *  @param aexp_mass [input] Exponent for the increase factor of the 
   *			      central enthalpy to converge to the 
   *			      requested baryon mass
   *  @param diff [output]   1-D {\tt Tbl} of size 1 for the storage of 
   *			    some error indicators : \\
   *	    {\tt diff(0)} : Relative change in the enthalpy field
   *			      between two successive steps \\
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
  void equilibrium_mag(double ent_c, double omega0, double fact_omega, 
		       int nzadapt, const Tbl& ent_limit, const Itbl& icontrol, 
		       const Tbl& control, double mbar_wanted, double aexp_mass, 
		       Tbl& diff, const double Q0, const double a_j0, 
		       Cmp (*f_j)(const Cmp& x, const double), 
		       Cmp (*M_j)(const Cmp& x,const double));

};

#endif

