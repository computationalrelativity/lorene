/*
 *  Definition of Lorene class Et_rot_bifluid
 *
 */

/*
 *   Copyright (c) 2001 Jerome Novak
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


#ifndef __ET_ROT_BIFLUID_H_ 
#define __ET_ROT_BIFLUID_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2002/10/09 07:54:29  j_novak
 * Et_rot_bifluid and Et_rot_mag inheritate virtually from Etoile_rot
 *
 * Revision 1.6  2002/09/13 09:17:33  j_novak
 * Modif. commentaires
 *
 * Revision 1.5  2002/04/05 09:09:36  j_novak
 * The inversion of the EOS for 2-fluids polytrope has been modified.
 * Some errors in the determination of the surface were corrected.
 *
 * Revision 1.4  2002/01/16 15:03:28  j_novak
 * *** empty log message ***
 *
 * Revision 1.3  2002/01/08 14:43:53  j_novak
 * better determination of surfaces for 2-fluid stars
 *
 * Revision 1.2  2002/01/03 15:30:27  j_novak
 * Some comments modified.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2001/10/03  09:49:06  novak
 * *** empty log message ***
 *
 * Revision 1.2  2001/08/28 14:14:10  novak
 * overrided l_surf function
 *
 * Revision 1.1  2001/06/22 15:38:52  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "eos_bifluid.h"
#include "etoile.h"

// Local prototype (for determining the surface)
Cmp prolonge_c1(const Cmp& uu, const int nzet) ;

/**
 * Class for two-fluid rotating relativistic stars.
 * 
 * This is a child class of {\tt Etoile\_rot}, with the same metric
 * and overloaded member functions. Still, there are two densities 
 * fields (and 2 log-enthalpies, see class {\tt Eos\_bifluid}), 
 * as well as two velocity fields.
 * Fluid 1 is supposed to correspond to neutrons, whereas fluid 2
 * may correspond to protons.
 * The quantity {\tt u_euler} in this class {\bf \em is not} the fluid
 * 3-velocity with respect to the eulerian observer; but 
 * $\sqrt{g_{\varphi\varphi}}J^{\varphi}$ expressed in a cartesian 
 * triad.
 * The quantity {\tt uuu} is still the norm of the fluid (no.1) 
 * 3-velocity, seen by the eulerian observer.
 *
 * @version #$Id$#
 */
class Et_rot_bifluid : virtual public Etoile_rot {
  
  // Data : 
  // -----
 protected:
  const Eos_bifluid& eos ; /// Equation of state for two-fluids model
  
  double omega2 ; /// Rotation angular velocity for fluid 2 ({\tt [f\_unit]})
  
  // Fluid quantities with respect to the fluid frame
  // ------------------------------------------------

  /// Log-enthalpy for the second fluid
  Tenseur ent2 ;
  
  Tenseur nbar2 ; /// Baryon density in the fluid frame, for fluid 2
  
  /**
   * In the relativistic case: relative Lorentz factor between both fluids;
   * in the newtonian case: square of the difference of velocities.
   * It is noted $\Delta^2$ in Prix et al. (see also {\tt Eos\_bifluid})
   */
  Tenseur xxx2 ; 
  
  // Fluid quantities with respect to the Eulerian frame
  // ---------------------------------------------------
  /// Lorentz factor between the fluid 2 and Eulerian observers   
  Tenseur gam_euler2 ;
  
  /// Norm of the (fluid no.2) 3-velocity with respect to the eulerian observer
  Tenseur uuu2 ;

  // Derived data : 
  // ------------
 protected:
  /// Coordinate radius at $\phi=0$, $\theta=\pi/2$. 
  mutable double* p_ray_eq2 ; 
  
  /// Coordinate radius at $\phi=\pi/2$, $\theta=\pi/2$. 
  mutable double* p_ray_eq2_pis2 ;
  
  /// Coordinate radius at $\phi=\pi$, $\theta=\pi/2$. 
  mutable double* p_ray_eq2_pi ;
  
  /// Coordinate radius at $\theta=0$. 
  mutable double* p_ray_pole2 ;
	
  /** Description of the surface of fluid 2: 2-D {\tt Itbl} containing the 
   *	values of the domain index {\it l} on the surface at the 
   *	collocation points in $(\theta', \phi')$
   */
  mutable Itbl* p_l_surf2 ; 
	
  /** Description of the surface of fluid 2: 2-D {\tt Tbl} containing the 
   *	values of the radial coordinate $\xi$ on the surface at the 
   *	collocation points in $(\theta', \phi')$
   */
  mutable Tbl* p_xi_surf2 ; 
	
  mutable double* p_r_circ2 ;	/// Circumferential radius of fluid no.2
  mutable double* p_aplat2 ;	/// Flatening r\_pole/r\_eq of fluid no.2

  // Constructors - Destructor
  // -------------------------
 public:

  Et_rot_bifluid(Map& mp_i, int nzet_i, bool relat, 
		 const Eos_bifluid& eos_i) ;     /// Standard constructor

  Et_rot_bifluid(const Et_rot_bifluid& ) ;       /// Copy constructor

  /** Constructor from a file (see {\tt sauve(FILE* )}) Works only for 
   *  relativistic stars.
   *  This has to be improved....
   */ 
  Et_rot_bifluid(Map& mp_i, const Eos_bifluid& eos_i, FILE* fich) ;    

  virtual ~Et_rot_bifluid() ;			/// Destructor
 

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

  /// Assignment to another Et_rot_bifluid
  void operator=(const Et_rot_bifluid&) ;	

  /// Sets both enthalpy profiles
  void set_enthalpies(const Cmp&, const Cmp&) ;

  /** Computes a spherical static configuration. 
   * 
   *  @param ent_c [input] central value of the enthalpy 1
   *  @param ent_c2 [input] central value of the enthalpy 2
   *  @param precis [input] threshold in the relative difference between 
   *	the enthalpy fields of two consecutive steps
   *	to stop the iterative procedure (default value: 1.e-14)
   */
  void equilibrium_spher_bi(double ent_c, double ent_c2, 
			 double precis = 1.e-14) ;
	
  /** Computes a spherical static configuration. 
   *  The sources for Poisson equations are regularized
   *  by extracting analytical diverging parts.
   * 
   *  @param ent_c [input] central value of the enthalpy 1
   *  @param ent_c2 [input] central value of the enthalpy 2
   *  @param precis [input] threshold in the relative difference between 
   *	the enthalpy fields of two consecutive steps
   *	to stop the iterative procedure (default value: 1.e-14)
   */
  void equil_spher_regular(double ent_c, double ent_c2, 
			   double precis = 1.e-14) ;

  // Accessors
  // ---------
 public:

  /// Returns the equation of state
  const Eos_bifluid& get_eos() const {return eos; } ;

  /// Returns the enthalpy field for fluid 2 
  const Tenseur& get_ent2() const {return ent2 ; } ;

  /// Returns the proper baryon density for fluid 2
  const Tenseur& get_nbar2() const {return nbar2 ; } ;
	
  /// Returns the coupling factor between both fluids $n_1n_2\Gamma_\Delta$
  const Tenseur& get_xxx2() const {return xxx2 ; } ;

  /// Returns the Lorentz factor between the fluid 2 and Eulerian observers
  const Tenseur& get_gam_euler2() const {return gam_euler2 ; } ;

  /// Returns the rotation angular velocity of fluid 2({\tt [f\_unit]}) 
  double get_omega2() const {return omega2 ; } ;

  /// Returns the norm of the fluid 2 3-velocity with respect to the eulerian frame
  const Tenseur& get_uuu2() const {return uuu2 ; } ;

  // Outputs
  // -------
 public:
  virtual void sauve(FILE *) const ;	    /// Save in a file
    
  /// Operator >> (virtual function called by the operator <<). 
  virtual ostream& operator>>(ostream& ) const ;    

  /// Printing of some informations, excluding all global quantities
  virtual void partial_display(ostream& ) const ;    

  // Global quantities
  // -----------------
 public:
	
  /** Description of the surface of fluid 1: returns a 2-D {\tt Itbl} 
   *	containing the 
   *	values of the domain index {\it l} on the surface at the 
   *	collocation points in $(\theta', \phi')$.
   *	The stellar surface is defined as the location where
   *	the enthalpy 2 (member {\tt ent2}) vanishes.
   */
  virtual const Itbl& l_surf() const ; 
	
  /** Description of the surface of fluid 2: returns a 2-D {\tt Itbl} 
   *	containing the 
   *	values of the domain index {\it l} on the surface at the 
   *	collocation points in $(\theta', \phi')$.
   *	The stellar surface is defined as the location where
   *	the enthalpy 2 (member {\tt ent2}) vanishes.
   */
  const Itbl& l_surf2() const ; 
	
  /** Description of the surface of fluid 2: returns a 2-D {\tt Tbl} 
   *	containing the values of the radial coordinate $\xi$ 
   *	on the surface at the 
   *	collocation points in $(\theta', \phi')$. 
   *	The stellar surface is defined as the location where
   *	the enthalpy 2 (member {\tt ent2}) vanishes.
   */
  const Tbl& xi_surf2() const ; 

  /// Coordinate radius for fluid 2 at $\phi=0$, $\theta=\pi/2$ [r\_unit].
  double ray_eq2() const ; 
	
  /// Coordinate radius for fluid 2 at $\phi=\pi/2$, $\theta=\pi/2$ [r\_unit].
  double ray_eq2_pis2() const ; 
	
  /// Coordinate radius for fluid 2 at $\phi=\pi$, $\theta=\pi/2$ [r\_unit].
  double ray_eq2_pi() const ; 
	
  /// Coordinate radius for fluid 2 at $\theta=0$ [r\_unit]. 
  double ray_pole2() const ; 
    
  virtual double mass_b() const ;	    /// Baryon mass
  virtual double mass_g() const ;	    /// Gravitational mass
  virtual double angu_mom() const ;  /// Angular momentum
  virtual double grv2() const ;	/// Error on the virial identity GRV2

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

  virtual double r_circ2() const ;     /// Circumferential radius for fluid 2
  virtual double aplat2() const ;      /// Flatening r\_pole/r\_eq for fluid 2

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
  /** Computes the hydrodynamical quantities relative to the Eulerian
   *  observer from those in the fluid frame.
   *
   *  The calculation is performed starting from the quantities
   *  {\tt ent}, {\tt ent2}, {\tt ener}, {\tt press}, and {\tt a\_car},  
   *  which are supposed to be up to date.  
   *  From these,  the following fields are updated:
   *  {\tt gam\_euler}, {\tt gam\_euler2}, {\tt u\_euler}, 
   *  {\tt ener\_euler}, {\tt s\_euler}. 
   * 
   */
  virtual void hydro_euler() ; 
	
  /** Computes the proper baryon and energy densities, as well as
   *  pressure from the enthalpies and both velocities.
   */
  void equation_of_state() ; 
	
  /** Computes an equilibrium configuration.
   *  
   *  @param ent_c  [input] Central enthalpy for fluid 1 
   *  @param ent_c2 [input] Central enthalpy for fluid 2
   *  @param omega0  [input] Requested angular velocity for fluid 1
   *  @param omega20  [input] Requested angular velocity for fluid 2
   *  @param ent_limit [input] 1-D {\tt Tbl} of dimension {\tt nzet} which
   *				defines the enthalpy for fluid 1 at the 
   *                            outer boundary of each domain
   *  @param ent2_limit [input] 1-D {\tt Tbl} of dimension {\tt nzet} which
   *				defines the enthalpy for fluid 2 at the 
   *                            outer boundary of each domain
   *  @param icontrol [input] Set of integer parameters (stored as a
   *			    1-D {\tt Itbl} of size 5) to control the 
   *			    iteration: \\
   *	{\tt icontrol(0) = mer\_max} : maximum number of steps \\
   *	{\tt icontrol(1) = mer\_rot} : step at which the rotation is 
   *				      switched on \\
   *	{\tt icontrol(2) = mer\_change\_omega} : step at which the rotation
   *			  velocity is changed to reach the final one  \\
   *	{\tt icontrol(3) = mer\_fix\_omega} :  step at which the final
   *			    rotation velocity must have been reached  \\
   *	{\tt icontrol(4) = mermax\_poisson} : maximum number of steps in 
   *				{\tt Map\_et::poisson} \\
   *  @param control [input] Set of parameters (stored as a 
   *			    1-D {\tt Tbl} of size 5) to control the 
   *			    iteration: \\
   *	{\tt control(0) = precis} : threshold on the enthalpy relative 
   *				change for ending the computation \\
   *	{\tt control(1) = omega\_ini} : initial angular velocity, 
   *			    switched on only if {\tt mer\_rot < 0}, 
   *			    otherwise 0 is used \\ 
   *	{\tt control(2) = omega2\_ini} : initial angular velocity, 
   *			    switched on only if {\tt mer\_rot < 0}, 
   *			    otherwise 0 is used \\ 
   *	{\tt control(3) = relax} : relaxation factor in the main 
   *				   iteration \\ 
   *	{\tt control(4) = relax\_poisson} : relaxation factor in 
   *				   {\tt Map\_et::poisson}\\ 
   *  @param diff [output]   1-D {\tt Tbl} of size 8 for the storage of 
   *			    some error indicators : \\
   *	    {\tt diff(0)} : Relative change in the enthalpy field 1
   *			      between two successive steps \\
   *	    {\tt diff(1)} : Relative change in the enthalpy field 2
   *			      between two successive steps \\
   *	    {\tt diff(2)} : Relative error in the resolution of the
   *			    Poisson equation for {\tt nuf} \\  
   *	    {\tt diff(3)} : Relative error in the resolution of the
   *			    Poisson equation for {\tt nuq} \\  
   *	    {\tt diff(4)} : Relative error in the resolution of the
   *			    Poisson equation for {\tt dzeta} \\  
   *	    {\tt diff(5)} : Relative error in the resolution of the
   *			    Poisson equation for {\tt tggg} \\  
   *	    {\tt diff(6)} : Relative error in the resolution of the
   *			    equation for {\tt shift} (x comp.) \\  
   *	    {\tt diff(7)} : Relative error in the resolution of the
   *			    equation for {\tt shift} (y comp.) \\  
   */
  void equilibrium_bi(double ent_c, double ent_c2, double omega0, 
		   double omega20, const Tbl& ent_limit, 
		   const Tbl& ent2_limit, const Itbl& icontrol, 
		   const Tbl& control, Tbl& diff) ;
	



};

#endif
