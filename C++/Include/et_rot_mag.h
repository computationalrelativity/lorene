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
 * Class for magnetized (perfect conductor), rigidly rotating stars.
 *
 * This is a child class of {\tt Etoile\_rot}, with the same metric
 * and overloaded member functions. Ecrire la suite bientot.
 *
 * @version #$Id$#
 */
class Et_rot_mag : public Etoile_rot {
  
  // Data : 
  // -----
 protected:

  Cmp A_t ; /// t-component of the elecctromagnetic potential 1-form
  Cmp A_phi; /// $\varphi$-component of the electromagnetic potential 1-form
  Cmp j_t;
  Cmp j_phi;

  Tenseur E_em;
  Tenseur Jp_em;
  Tenseur Srr_em; // Stt_em = - Srr_em...
  Tenseur Spp_em; 

  double Q ;
  double a_j ;

  
  // Constructors - Destructor
  // -------------------------
 public:

  /// Standard destructor

  Et_rot_mag(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i); 


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

  const Cmp& get_At() const {return A_t ; } ;
  const Cmp& get_Aphi() const {return A_phi ;} ;
  const Cmp& get_jt() const {return j_t ; } ;
  const Cmp& get_jphi() const {return j_phi ;} ;
  const Tenseur& get_Eem() const {return E_em ; } ;
  const Tenseur& get_Jpem() const {return Jp_em ;} ;
  const Tenseur& get_Srrem() const {return Srr_em ; } ;
  const Tenseur& get_Sppem() const {return Spp_em ;} ;
  double get_Q() const {return Q ;} ;
  double get_a_j() const {return a_j ;} ;

  // Outputs
  // -------
 public:
    
  /// Operator >> (virtual function called by the operator <<). 
  virtual ostream& operator>>(ostream& ) const ;    

  // Global quantities
  // -----------------
 public:

  /** Description of the surface of fluid 1: returns a 2-D {\tt Itbl} 
   *	containing the 
   *	values of the domain index $l$ on the surface at the 
   *	collocation points in $(\theta', \phi')$.
   *	The stellar surface is defined as the location where
   *	the enthalpy 2 (member {\tt ent2}) vanishes.
   */

  Tenseur Elec() const ; // computes E
  Tenseur Magn() const ; // computes B
  void MHD_comput() ; // computes T_{mu,nu}^{EM}
  virtual double mass_g() const ;	    /// Gravitational mass
  virtual double angu_mom() const ;  /// Angular momentum
  virtual double grv2() const ;	/// Error on the virial identity GRV2
  virtual double tsw() const ; // Ratio T/W
  double MagMom() const ; /// Magnetic Momentum
  double Q_comput() const; /// Computed charge
  double Q_int() const; /// Interior charge
  double GyroMag() const ; /// Gyromagnetic ratio

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
   *  The quadrupole moment $Q$ is defined according to Eq. (7) of
   *  [Salgado, Bonazzola, Gourgoulhon and Haensel, Astron. Astrophys.
   *   {\bf 291}, 155 (1994)]. At the Newtonian limit it is related to
   *  the component ${\bar I}_{zz}$ of the MTW (1973) reduced quadrupole 
   *  moment ${\bar I}_{ij}$ by: $Q = -3/2 {\bar I}_{zz}$. 
   *  Note that $Q$ is the negative of the quadrupole moment defined 
   *  by Laarakkers and Poisson, Astrophys. J. {\bf 512}, 282 (1999).
   */
  virtual double mom_quad() const ;	


  // Rajouter les nouvelles fonctions membres ici, comme :
  // Moment magnetique, champs E et B, pressions mag au pole et equateur



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

  virtual void magnet_comput(const int conduc, const int adapt_flag,
			     Cmp (*f_j)(const Cmp& x, const double a_j),
		      Param& par_poisson_At, Param& par_poisson_Avect) ;
	
  void equilibrium_mag(double ent_c, double omega0, double fact_omega, 
		 int nzadapt, const Tbl& ent_limit, const Itbl& icontrol, 
		 const Tbl& control, double mbar_wanted, double aexp_mass, 
		 Tbl& diff, const double Q0, const double a_j0, 
		       Cmp (*f_j)(const Cmp& x, const double a_j), 
		       Cmp (*M_j)(const Cmp& x,const double a_j));

};

#endif

