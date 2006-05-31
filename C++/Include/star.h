/*
 *  Definition of Lorene classes Star and Star_bin
 *				 
 */

/*
 *   Copyright (c) 2004 Francois Limousin
 *
 *   Copyright (c) 2000-2001 Eric Gourgoulhon (for preceding class Etoile)
 *   Copyright (c) 2000-2001 Keisuke Taniguchi (for preceding class Etoile)
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


#ifndef __STAR_H_ 
#define __STAR_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.22  2006/05/31 09:25:47  f_limousin
 * Modif. of the size of the different domains
 *
 * Revision 1.21  2006/04/11 14:26:12  f_limousin
 * New version of the code : improvement of the computation of some
 * critical sources, estimation of the dirac gauge, helical symmetry...
 *
 * Revision 1.20  2005/09/15 15:56:28  e_gourgoulhon
 * Made the documentation compliant with Doxygen.
 *
 * Revision 1.19  2005/09/13 19:38:32  f_limousin
 * Reintroduction of the resolution of the equations in cartesian coordinates.
 *
 * Revision 1.18  2005/08/13 16:14:11  m_saijo
 * Corrected the document of the total shift vector
 *
 * Revision 1.17  2005/04/08 12:36:45  f_limousin
 * Just to avoid warnings...
 *
 * Revision 1.16  2005/02/24 16:09:29  f_limousin
 * Change the name of some variables (for instance dcov_logn --> dlogn).
 * Add also member dlnq but delete dlnpsi_auto and dlogn_auto.
 *
 * Revision 1.15  2005/02/17 17:28:18  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.14  2005/02/11 18:11:16  f_limousin
 * Introduction of a member Map_af in derived class Star_bin.
 *
 * Revision 1.13  2004/11/11 16:29:48  j_novak
 * set_der_0x0 is no longer virtual (to be coherent with Tensor/Scalar classes).
 *
 * Revision 1.12  2004/11/10 16:31:53  j_novak
 * Star is no longer an abstract class (mass_b and mass_g are no longer
 * pure virtual). Modified comments to be readable by doxygen.
 *
 * Revision 1.11  2004/07/21 11:48:30  f_limousin
 * Remove function sprod.
 *
 * Revision 1.10  2004/06/22 12:47:01  f_limousin
 * Change qq, qq_auto and qq_comp to beta, beta_auto and beta_comp.
 *
 * Revision 1.9  2004/05/25 14:48:57  f_limousin
 * Add a parameter for the function equilibrium.
 *
 * Revision 1.8  2004/03/23 09:53:50  f_limousin
 * Minor changes
 *
 * Revision 1.7  2004/02/27 09:41:52  f_limousin
 * Scalars ssjm1_logn, ssjm1_qq ... for all metric coefficients have been
 * in class Star_bin for the resolution of Poisson equations.
 * The class Star is now abstract : the computational routines mass_b()
 * and mass_g() = 0.
 *
 * Revision 1.6  2004/01/22 10:06:33  f_limousin
 * Add methods set_logn_comp() and set_shift_auto().
 *
 * Revision 1.5  2004/01/20 15:26:00  f_limousin
 * New class star and star_bin.
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "tensor.h"
#include "metric.h"

class Eos ;

			//---------------------------//
			//    base class Star        //
			//---------------------------//

/**
 * Base class for stars. *** UNDER DEVELOPMENT ***  \ingroup (star)
 * 
 * A \c Star is constructed upon (i) a mapping 
 * (derived class of \c Map), the center of which defines the center of the 
 * star, and (ii) an equation of state (derived class of \c Eos).  
 * It contains tensor fields (class \c Tensor) which describle the
 * hydrodynamical quantities as well as the gravitational field (spacetime
 * metric). 
 * 
 * According to the 3+1 formalism, the spacetime metric is written
 * \f[
 *   ds^2 = - N^2  dt^2 + \gamma_{ij} ( dx^i + \beta^i dt )
 *               (dx^j + \beta^j dt )
 * \f]
 * where \f$\gamma_{ij}\f$ is the 3-metric.
 * 
 * The 3+1 formalism introduces two kinds of priviledged observers: the
 * fluid comoving observer and the Eulerian observer, whose 4-velocity
 * is the unit future directed normal to the \e t = const hypersurfaces. 
 * The hydrodynamical quantities measured by the fluid observer correspond
 * to the members \c ent, \c nbar, \c ener, and \c press. 
 * The hydrodynamical quantities measured by the Eulerian observer correspond
 * to the members \c ener_euler,  \c s_euler, \c gam_euler, and 
 * \c u_euler.  
 * 
 * @version #$Id$#
 */
class Star {

    // Data : 
    // -----
    protected:
	Map& mp ;	    ///< Mapping associated with the star

	/// Number of domains of \c *mp occupied by the star
	int nzet ;
	
	const Eos& eos ;   ///< Equation of state of the stellar matter
    
	// Fluid quantities with respect to the fluid frame
	// ------------------------------------------------

	Scalar ent ;	   ///< Log-enthalpy 
 	
	Scalar nbar ; 	   ///< Baryon density in the fluid frame
	Scalar ener ;	   ///< Total energy density in the fluid frame
	Scalar press ;	   ///< Fluid pressure

	// Fluid quantities with respect to the Eulerian frame
	// ---------------------------------------------------
	Scalar ener_euler ; ///< Total energy density in the Eulerian frame

	/// Trace of the stress scalar in the Eulerian frame
	Scalar s_euler ;   

	/// Lorentz factor between the fluid and Eulerian observers 
	Scalar gam_euler ; 
	
	/// Fluid 3-velocity with respect to the Eulerian observer
	Vector u_euler ; 
 
	/** Spatial part of the stress-energy tensor with respect
	 * to the Eulerian observer. 
	 */
	Tensor stress_euler ;
	
	// Metric potentials
	// -----------------
	
	/** Total of the logarithm of the lapse \e N .
	 *   In the Newtonian case, 
	 *   this is the Newtonian gravitational potential
	 *   (in units of \f$c^2\f$). 
	 */
	Scalar logn ;

	/// Total lapse function 
	Scalar nn ; 
	
	/// Total shift vector
	Vector beta ;
	
	/// Total \f$\ln Q = \ln(\psi^2  N)\f$
    Scalar lnq ;
	
	/// Total 3-metric 
	Metric gamma ;


    // Derived data : 
    // ------------
    protected:
	/// Coordinate radius at \f$\phi=0\f$, \f$\theta=\pi/2\f$. 
	mutable double* p_ray_eq ; 

	/// Coordinate radius at \f$\phi=\pi/2\f$, \f$\theta=\pi/2\f$. 
	mutable double* p_ray_eq_pis2 ;

	/// Coordinate radius at \f$\phi=\pi\f$, \f$\theta=\pi/2\f$. 
	mutable double* p_ray_eq_pi ;
	
	/// Coordinate radius at \f$\theta=0\f$. 
	mutable double* p_ray_pole ;
	
	/** Description of the stellar surface: 2-D \c Itbl containing the 
	 *	values of the domain index \e l on the surface at the 
	 *	collocation points in \f$(\theta', \phi')\f$
	 */
	mutable Itbl* p_l_surf ; 
	
	/** Description of the stellar surface: 2-D \c Tbl containing the 
	 *	values of the radial coordinate \f$\xi\f$ on the surface at the 
	 *	collocation points in \f$(\theta', \phi')\f$
	 */
	mutable Tbl* p_xi_surf ; 
	
	mutable double* p_mass_b ;	///< Baryon mass
	mutable double* p_mass_g ;	///< Gravitational mass


    // Constructors - Destructor
    // -------------------------
    public:
    
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param nzet_i Number of domains occupied by the star
	 * @param eos_i Equation of state of the stellar matter
	 * 
	 */
	Star(Map& mp_i, int nzet_i, const Eos& eos_i) ;			
	
	
	Star(const Star& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE* )). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param fich	input file (must have been created by the function
	 *	\c sauve)
	 */
	Star(Map& mp_i, const Eos& eos_i, FILE* fich) ;    		

	virtual ~Star() ;			///< Destructor

	
    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
	virtual void set_der_0x0() const ; 

	/** Sets to \c ETATNONDEF (undefined state) the hydrodynamical 
	 *  quantities relative to the Eulerian observer.
	 */
	virtual void del_hydro_euler() ; 
	

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Star
	void operator=(const Star&) ;	
	
	/// Read/write of the mapping
	Map& set_mp() {return mp; } ; 

	/// Assignment of the enthalpy field.
	void set_enthalpy(const Scalar& ) ; 
	
	/** Computes the proper baryon and energy density, as well as
	 *  pressure from the enthalpy.
	 */
	void equation_of_state() ; 
	
	/** Computes the hydrodynamical quantities relative to the Eulerian
	 *  observer from those in the fluid frame (\c nbar, \c ener
	 *  and \c press).
	 */
	 virtual void hydro_euler() ; 
	
	/** Computes a spherical static configuration. 
	 * 
	 *  @param ent_c [input] central value of the enthalpy
	 *  @param precis [input] threshold in the relative difference between 
	 *	the enthalpy fields of two consecutive steps
	 *	to stop the iterative procedure (default value: 1.e-14)
	 */
	void equilibrium_spher(double ent_c, double precis = 1.e-14) ; 

    // Accessors
    // ---------
    public:
	/// Returns the mapping
	const Map& get_mp() const {return mp; } ; 

	/// Returns the number of domains occupied by the star
	int get_nzet() const {return nzet; } ; 

	/// Returns the equation of state
	const Eos& get_eos() const {return eos; } ; 

	/// Returns the enthalpy field 
	const Scalar& get_ent() const {return ent;} ;

	/// Returns the proper baryon density
	const Scalar& get_nbar() const {return nbar;} ;

	/// Returns the proper total energy density
	const Scalar& get_ener() const {return ener;} ;

	/// Returns the fluid pressure
	const Scalar& get_press() const {return press;} ;

	/// Returns the total energy density with respect to the Eulerian observer
	const Scalar& get_ener_euler() const {return ener_euler;} ;

	/// Returns the trace of the stress tensor in the Eulerian frame
	const Scalar& get_s_euler() const {return s_euler;} ;

	/// Returns the Lorentz factor between the fluid and Eulerian observers
	const Scalar& get_gam_euler() const {return gam_euler;} ;

	/// Returns the fluid 3-velocity with respect to the Eulerian observer
	const Vector& get_u_euler() const {return u_euler;} ;

	 /** Returns the spatial part of the stress-energy tensor 
	  *  with respect to the Eulerian observer
	  */
	const Tensor& get_stress_euler() const {return stress_euler;} ;

	/** Returns the logarithm of the total lapse \e N.
	 *   In the Newtonian case, this is the Newtonian
	 *   gravitational potential (in units of \f$c^2\f$). 
	 */
	const Scalar& get_logn() const {return logn;} ;

	/// Returns the total lapse function \e N
	const Scalar& get_nn() const {return nn;} ;

	/// Returns the total shift vector \f$\beta^i\f$.
	const Vector& get_beta() const {return beta;} ;

	/// Returns the total scalar field \f$\ln Q\f$.
	const Scalar& get_lnq() const {return lnq;} ;

	/// Returns the total metric \f$\gamma\f$.
	const Metric& get_gamma() const {return gamma;} ;

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	    ///< Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Star& ) ;	

    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

    // Global quantities
    // -----------------
    public:
	/// Coordinate radius at \f$\phi=0\f$, \f$\theta=\pi/2\f$ [r_unit].
	double ray_eq() const ; 
	
	/// Coordinate radius at \f$\phi=\pi/2\f$, \f$\theta=\pi/2\f$ [r_unit].
	double ray_eq_pis2() const ; 
	
	/// Coordinate radius at \f$\phi=\pi\f$, \f$\theta=\pi/2\f$ [r_unit].
	double ray_eq_pi() const ; 
	
	/// Coordinate radius at \f$\theta=0\f$ [r_unit]. 
	double ray_pole() const ; 
    
	/** Description of the stellar surface: returns a 2-D \c Itbl 
	 *	containing the 
	 *	values of the domain index \e l on the surface at the 
	 *	collocation points in \f$(\theta', \phi')\f$.
	 *	The stellar surface is defined as the location where
	 *	the enthalpy (member \c ent) vanishes.
	 */
	virtual const Itbl& l_surf() const ; 
	
	/** Description of the stellar surface: returns a 2-D \c Tbl 
	 *	containing the values of the radial coordinate \f$\xi\f$ 
	 *	on the surface at the 
	 *	collocation points in \f$(\theta', \phi')\f$. 
	 *	The stellar surface is defined as the location where
	 *	the enthalpy (member \c ent) vanishes.
	 */
	const Tbl& xi_surf() const ; 

	/// Baryon mass
    	virtual double mass_b() const = 0  ;
	
	/// Gravitational mass
    	virtual double mass_g() const = 0  ;
	
};
ostream& operator<<(ostream& , const Star& ) ;	


			//---------------------------//
			//    class Star_bin       //
			//---------------------------//

/**
 * Class for stars in binary system. *** UNDER DEEVELOPMENT *** \ingroup (star)
 *
 * A \c Star_bin can be construted in two states, represented by
 * the \c bool member \c irrotational: (i) irrotational
 * (i.e. the fluid motion is irrotational) or (ii) rigidly corotating 
 * with respect to the orbital motion (synchronized binary). 
 *
 * @version #$Id$#
 */
class Star_bin : public Star {

    // Data : 
    // -----
    protected:
	/** \c true for an irrotational star, \c false for a
	 *  corotating one
	 */
	bool irrotational ; 
	
	/** Scalar potential \f$\Psi_0\f$ of the non-translational part of
	 *  fluid 4-velocity (in the irrotational case)
	 */
	Scalar psi0 ; 

	/** Gradient of \f$\Psi\f$ (in the irrotational case)
	 *  (Spherical components with respect to the mapping of the star)
	 */
	Vector d_psi ; 
	
	/** Spatial projection of the fluid 3-velocity with respect to  
	 *  the co-orbiting observer. 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	Vector wit_w ; 
	
	/** Logarithm of the Lorentz factor between the fluid and 
	 *  the co-orbiting observer.
	 */
	Scalar loggam ; 

	/** 3-vector shift, divided by \e N, of the rotating coordinates,
	 *  \f$\beta^i/N\f$. 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	Vector bsn ; 
	
	/// Centrifugal potential
	Scalar pot_centri ; 	


	/** Part of the lapse logarithm (gravitational potential at the
	 *  Newtonian limit) generated principaly by the star. 
	 */
	Scalar logn_auto ; 

	/** Part of the lapse logarithm (gravitational potential at the
	 *  Newtonian limit) generated principaly by the companion star. 
	 */
	Scalar logn_comp ; 

	/// Covariant derivative of the total logarithm of the lapse. 
	Vector dcov_logn ;

	/// Contravariant derivative of the total logarithm of the lapse. 
	Vector dcon_logn ;
	
	/** Scalar field \f$ Q = \psi^2 N \f$ generated principaly by the
	 *  star.
	 */
	Scalar lnq_auto ;

	/** Scalar field \f$ Q = \psi^2 N \f$ generated principaly by the
	 *  companion star.
	 */
	Scalar lnq_comp ;

	
	/// Conformal factor \f$\psi^4\f$
	Scalar psi4 ;

	/// Covariant derivative of the logarithm of the conformal factor
	Vector dcov_phi ;
	/// Contravariant derivative of the logarithm of the conformal factor
	Vector dcon_phi ;

	/** Flat metric defined on the mapping (Spherical components 
	 * with respect to the mapping of the star) .
	 */
	Metric_flat flat ;

	/// Conformal metric \f$\tilde \gamma_{ij}\f$
	Metric gtilde ;

	/** Part of the shift vector generated principaly by the star
	 *  (Spherical components with respect to the mapping of the star)
	 */
	Vector beta_auto ; 

	/** Part of the shift vector generated principaly by the star
	 *  (Spherical components with respect to the mapping of the star)
	 */
	Vector beta_comp ; 


	/** Total deviation of the inverse conformal metric 
	 *  \f$\tilde \gamma^{ij}\f$ from the inverse flat metric. 
	 */
	Sym_tensor hij ; 


	/** Deviation of the inverse conformal metric 
	 *  \f$\tilde \gamma^{ij}\f$ from the inverse flat metric generated
	 *  principaly by the star. 
	 */
	Sym_tensor hij_auto ;

	/** Deviation of the inverse conformal metric 
	 *  \f$\tilde \gamma^{ij}\f$ from the inverse flat metric generated
	 *  principaly by the companion star. 
	 */
	Sym_tensor hij_comp ;

	/** Part of the extrinsic curvature tensor \f$\tilde K^{ij}\f$
	 *  generated by \c beta_auto. 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	Sym_tensor tkij_auto ;
	
	/** Part of the extrinsic curvature tensor \f$\tilde K^{ij}\f$
	 *  generated by \c beta_comp. 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	Sym_tensor tkij_comp ;
	
	/** Part of the scalar \f$K_{ij} K^{ij}\f$
	 *  generated by \c beta_auto, i.e. 
	 *    \f$K_{ij}^{\rm auto} K^{ij}_{\rm auto}\f$ 
	 */
	Scalar kcar_auto ;
	
	/** Part of the scalar \f$K_{ij} K^{ij}\f$
	 *  generated by \c beta_auto and \c beta_comp, i.e. 
	 *    \f$K_{ij}^{\rm auto} K^{ij}_{\rm comp}\f$ 
	 */
	Scalar kcar_comp ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c logn_auto.
	 */
	Scalar ssjm1_logn ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c lnq_auto.
	 */
	Scalar ssjm1_lnq ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c khi. (second scalar equation 
	 *  for the resolution of the vectorial poisson equation for the shift)
	 */
	Scalar ssjm1_khi ;

	Vector ssjm1_wbeta ;
	
	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c h00_auto.
	 */
	Scalar ssjm1_h11 ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c h10_auto.
	 */
	Scalar ssjm1_h21 ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c h20_auto.
	 */
	Scalar ssjm1_h31 ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c h11_auto.
	 */
	Scalar ssjm1_h22 ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c h21_auto.
	 */
	Scalar ssjm1_h32 ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c h22_auto.
	 */
	Scalar ssjm1_h33 ;

	/**
	  * Function used to construct the part \f$lnq_auto\f$ generated by 
	  * the star  from the total \f$lnq\f$.
	  * Mainly this \c Scalar is 1 around the star and 0 around 
	  * the companion
	  * and the sum of \c decouple for the star and his companion is 1 
	  * everywhere.
	  */
	 Scalar decouple ;

	 /** \c true  if the 3-metric is conformally flat, \c false
	  *  for a more general metric. 
	  */
	 bool conf_flat ;
	
    // Derived data : 
    // ------------
    protected:
	/// Absolute coordinate X of the barycenter of the baryon density
	mutable double* p_xa_barycenter ; 
	 

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param nzet_i Number of domains occupied by the star
	 * @param eos_i Equation of state of the stellar matter
	 * @param irrot should be \c true for an irrotational star, 
	 *		    \c false for a corotating one
	 * @param conf_flat should be \c true for a conformally flat metric
	 *                  \c false for a general one
	 */
	Star_bin(Map& mp_i, int nzet_i,  const Eos& eos_i,
		   bool irrot, bool conf_flat) ;			
	
	
	Star_bin(const Star_bin& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE* )). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param fich	input file (must have been created by the function
	 *	\c sauve)
	 */
	Star_bin(Map& mp_i, const Eos& eos_i, FILE* fich) ;
    		
	virtual ~Star_bin() ;			///< Destructor


    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
	virtual void set_der_0x0() const ; 

	/** Sets to \c ETATNONDEF (undefined state) the hydrodynamical 
	 *  quantities relative to the Eulerian observer.
	 */
	virtual void del_hydro_euler() ; 
	

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Star_bin
	void operator=(const Star_bin& ) ;	
	
	/// Read/write the centrifugal potential
	Scalar& set_pot_centri() ;
	
        /** Read/write of the logarithm of the lapse generated 
	 *  principaly by the companion.
	 */ 
	Scalar& set_logn_comp() ;

	/// Assignment of a new logn_auto
	void set_logn_auto(const Scalar& logn_auto_new) {logn_auto = logn_auto_new ;
	return ;}

	/// Assignment of a new lnq_auto
	void set_lnq_auto(const Scalar& lnq_auto_new) {lnq_auto = lnq_auto_new ;
	return ;}

	/// Read/write of \f$beta_auto\f$
	Vector& set_beta_auto() ;	

	/// Read/write of \f$beta\f$
	Vector& set_beta() ;	
	
	/// Write if conformally flat 
	void set_conf_flat(bool confflat) {conf_flat = confflat ; return ;}

    // Accessors
    // ---------
    public:
	/** Returns \c true for an irrotational motion, \c false for 
	 *  a corotating one. 
	 */
	bool is_irrotational() const {return irrotational; } ; 

	/// Returns the non-translational part of the velocity potential
	const Scalar& get_psi0() const {return psi0;} ;

	/** Returns the covariant derivative of the velocity potential 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	const Vector& get_d_psi() const {return d_psi;} ;

	/** Returns the spatial projection of the fluid 3-velocity with 
	 *  respect to the co-orbiting observer. 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	const Vector& get_wit_w() const {return wit_w;} ;

	/** Returns the logarithm of the Lorentz factor between the fluid and 
	 *  the co-orbiting observer.
	 */
	const Scalar& get_loggam() const {return loggam;} ;

	/** Returns the shift vector, divided by \e N, of the rotating 
	 *   coordinates, \f$\beta^i/N\f$. 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	const Vector& get_bsn() const {return bsn;} ;

	/// Returns the centrifugal potential
	const Scalar& get_pot_centri() const {return pot_centri;} ;

	/** Returns the part of the lapse logarithm (gravitational potential 
	 *  at the Newtonian limit) generated principaly by the star. 
	 */
	const Scalar& get_logn_auto() const {return logn_auto;} ;

	/** Returns the part of the lapse logarithm (gravitational potential 
	 *  at the Newtonian limit) generated principaly by the companion star.
 	 */
	const Scalar& get_logn_comp() const {return logn_comp;} ;

	/** Returns the part of the shift vector \f$\beta^i\f$ generated 
	 *  principaly by the star.
	 *  (Spherical components with respect to the mapping of the star)
	 */
	const Vector& get_beta_auto() const {return beta_auto;} ;

	/** Returns the part of the shift vector \f$\beta^i\f$ generated 
	 *  principaly by the star.
	 *  (Spherical components with respect to the mapping of the star)
	 */
	const Vector& get_beta_comp() const {return beta_comp;} ;

	/** Returns the part of the vector field \f$Q\f$ generated principaly 
	 *   by the star. 
	 */
	const Scalar& get_lnq_auto() const {return lnq_auto;} ;

	/** Returns the part of the vector field \f$Q\f$ generated principaly 
	 *   by the companion star. 
	 */
	const Scalar& get_lnq_comp() const {return lnq_comp;} ;

	/** Returns the covariant derivative of \f$logn\f$. 
	 */
	const Vector& get_dcov_logn() const {return dcov_logn;} ;

	/** Returns the contravariant derivative of \f$logn\f$.
	 */
	const Vector& get_dcon_logn() const {return dcon_logn;} ;

	/** Returns the covariant derivative of \f$\Phi\f$
     *  (logarithm of the conformal factor). 
	 */
	const Vector& get_dcov_phi() const {return dcov_phi;} ;

	/** Returns the contravariant derivative of \f$\Phi\f$
     *  (logarithm of the conformal factor). 
	 */
	const Vector& get_dcon_phi() const {return dcon_phi;} ;

	/// Return the conformal factor \f$\psi^4\f$
	const Scalar& get_psi4() const {return psi4;} ;

	/** Return the flat metric defined on the mapping (Spherical
	 *  components  with respect to the mapping of the star)
	 */
	const Metric& get_flat() const {return flat;} ;	

	/// Return the conformal 3-metric \f$\tilde \gamma\f$
	const Metric& get_gtilde() const {return gtilde;} ;

	
	/** Return the total deviation of the inverse conformal metric 
	 *  \f$\tilde \gamma^{ij}\f$ from the inverse flat metric. 
	 */
	const Sym_tensor& get_hij() const {return hij;} ;

	/** Return the deviation of the inverse conformal metric 
	 *  \f$\tilde \gamma^{ij}\f$ from the inverse flat metric principaly
	 *  generated by the star.
	 */
	const Sym_tensor& get_hij_auto() const {return hij_auto;} ;

	/** Return the deviation of the inverse conformal metric 
	 *  \f$\tilde \gamma^{ij}\f$ from the inverse flat metric generated
	 *  principaly by the companion star.
	 */
	const Sym_tensor& get_hij_comp() const {return hij_comp;} ;


	/** Returns the part of the extrinsic curvature tensor 
	 *  \f$\tilde K^{ij}\f$ generated by \c beta_auto. 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	const Sym_tensor& get_tkij_auto() const {return tkij_auto;} ;

	/** Returns the part of the extrinsic curvature tensor 
	 *  \f$\tilde K^{ij}\f$ generated by \c beta_comp. 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	const Sym_tensor& get_tkij_comp() const {return tkij_comp;} ;

	/** Returns the part of 
	 *  \f$\tilde K^{ij} \tilde K_{ij}\f$ generated by \c beta_auto. 
	 */
	const Scalar& get_kcar_auto() const {return kcar_auto;} ;

	/** Returns the part of 
	 *  \f$\tilde K^{ij} \tilde K_{ij}\f$ generated by \c beta_comp. 
	 */
	const Scalar& get_kcar_comp() const {return kcar_comp;} ;

	/**
	 * Returns the function used to construct \c beta_auto 
	 from \c beta.
	*/
	const Scalar get_decouple() const {return decouple ;}

	bool is_conf_flat() const {return conf_flat; } ; 

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	    ///< Save in a file
    
    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

    // Global quantities
    // -----------------
    public:
	/// Baryon mass
    	virtual double mass_b() const ;
	
	/// Gravitational mass
    	virtual double mass_g() const ;
	
	/// Absolute coordinate X of the barycenter of the baryon density, 
    	virtual double xa_barycenter() const ;


    // Computational routines
    // ----------------------
    public: 

	/** Computes the hydrodynamical quantities relative to the Eulerian
	 *  observer from those in the fluid frame, as well as 
	 *  \c wit_w and \c loggam.  
	 *
	 *  The calculation is performed starting from the quantities
	 *  \c ent, \c ener, \c press, \c a_car and \c bsn,  
	 *  which are supposed to be up to date.  
	 *  From these,  the following fields are updated:
	 *  \c gam_euler, \c u_euler, \c ener_euler, 
	 *  \c s_euler, \c stress_euler,  
	 *  \c wit_w and \c loggam. 
	 * 
	 */
	virtual void hydro_euler() ; 
	
	/** Computes metric coefficients from known potentials,
	 * when the companion is another star.
	 *
	 *  The calculation is performed starting from the quantities
	 *  \c logn_auto,  \c lnq_auto, \c beta_auto, 
	 *  \c hij_auto, \c comp.logn_auto,  \c comp.lnq_auto,
	 *  \c comp.beta_auto, \c comp.hij_auto
	 *  which are supposed to be up to date.
	 *  From these,  the following fields are updated:
	 *  \c logn_comp, \c lnq_comp, \c beta_comp,
	 *  \c hij_comp, \c nn,  \c psi4,  \c beta,
	 *  
	 *  @param comp companion star.
	 *  @param omega  angular velocity with respect to an asymptotically 
	 *		  inertial observer
	 */
	void update_metric(const Star_bin& comp, double omega) ;
	
	/** Same as \c update_metric(const Star_bin\& ) but with
	 *  relaxation.
	 *
	 *  @param comp companion star.
	 *  @param star_prev previous value of the star. 
	 *  @param relax relaxation parameter. 
	 *  @param omega  angular velocity with respect to an asymptotically 
	 *		  inertial observer
	 */
	void update_metric(const Star_bin& comp, const Star_bin& star_prev, 
			   double relax, double omega) ; 
	
	/** Computes the derivative of metric functions related to the
	 *  companion star.
	 *  @param omega  angular velocity with respect to an asymptotically 
	 *		  inertial observer
	 */
	 void update_metric_der_comp(const Star_bin& comp, double omega) ;

	/** Computes the quantities \c bsn and \c pot_centri.
	 * 
	 *  The calculation is performed starting from the quantities
	 *  \c nn, \c beta,  \c Q,  
	 *  which are supposed to be up to date.  
	 * 
	 *  @param omega  angular velocity with respect to an asymptotically 
	 *		  inertial observer
	 *  @param x_axe  absolute X coordinate of the rotation axis
	 */
	void kinematics(double omega, double x_axe) ; 
	
	/** Computes the gradient of the total velocity potential \f$\psi\f$. 
	 * 
	 */
	void fait_d_psi() ; 
	
	/** Computes \c tkij_auto and \c akcar_auto from 
	 *  \c beta_auto, \c nn and \c Q.
	 *  @param omega  angular velocity with respect to an asymptotically 
	 *		  inertial observer
	 */
	void extrinsic_curvature(double omega) ; 
		
	
	/** Computes an equilibrium configuration.
	 * 
	 *  @param ent_c  [input] Central enthalpy
	 *  @param mermax [input] Maximum number of steps 
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    poisson scalar
	 *  @param relax_poisson [input]  Relaxation factor in poisson scalar
	 *  @param mermax_potvit [input]  Maximum number of steps in 
	 *				  Map_radial::poisson_compact
	 *  @param relax_potvit [input]   Relaxation factor in 
	 *				  Map_radial::poisson_compact
	 *  @param thres_adapt  [input]   Threshold on dH/dr for the adaptation 
	 *				  of the mapping
	 *  @param diff [output]   1-D \c Tbl for the storage of some
	 *			    error indicators 
	 */
	void equilibrium(double ent_c, int mermax, int mermax_potvit, 
			 int mermax_poisson, double relax_poisson, 
			 double relax_potvit, double thres_adapt, Tbl& diff,
			 double om) ;


	/** Computes the non-translational part of the velocity scalar potential
	 *  \f$\psi0\f$ by solving the continuity equation.
	 *  
	 *  @param mermax  [input] Maximum number of steps in the iteration
	 *  @param precis  [input] Required precision: the iteration will
	 *			   be stopped when the relative difference
	 *			   on \f$\psi0\f$ between two successive steps
	 *			   is lower than \c precis.
	 *  @param relax   [input] Relaxation factor.  
	 *
	 *  @return Relative error of the resolution obtained by comparing
	 *	    the operator acting on the solution with the source.
	 */
	 double velocity_potential(int mermax, double precis, double relax) ;

	/** Performs a relaxation on \c ent, \c logn_auto,
	 *  \c lnq_auto, \c beta_auto  and \c hij_auto. 
	 * 
	 *  @param star_prev   [input] star at the previous step.
	 *  @param relax_ent   [input] Relaxation factor for \c ent 
	 *  @param relax_met   [input] Relaxation factor for \c logn_auto,
	 *			       \c lnq_auto, \c beta_auto, 
	 *			       only if \c (mer \% fmer_met == 0).
	 *  @param mer	       [input] Step number
	 *  @param fmer_met    [input] Step interval between metric updates
	 */
	void relaxation(const Star_bin& star_prev, double relax_ent, 
			double relax_met, int mer, int fmer_met) ;
	

	/// Test if the gauge conditions we impose are well satisfied
	void test_K_Hi() const ;

	/// Test of the helical symmetry
	void helical(double omega) const ;

	friend class Binary ;

};

#endif
