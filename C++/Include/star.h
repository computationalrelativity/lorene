/*
 *  Definition of Lorene classes Star
 *				 Star_bin
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
 * Base class for stars. *** UNDER DEVELOPMENT ***
 * 
 * A {\tt Star} is constructed upon (i) a mapping 
 * (derived class of {\tt Map}), the center of which defines the center of the 
 * star, and (ii) an equation of state (derived class of {\tt Eos}).  
 * It contains tensor fields (class {\tt Tensor}) which describle the
 * hydrodynamical quantities as well as the gravitational field (spacetime
 * metric). 
 * 
 * According to the 3+1 formalism, the spacetime metric is written
 * \begin{equation}
 *   ds^2 = - N^2  dt^2 + \gamma_{ij} ( dx^i + \beta^i dt )
 *               (dx^j + \beta^j dt )
 * \end{equation}
 * where $\gamma_{ij}$ is the 3-metric.
 * 
 * The 3+1 formalism introduces two kinds of priviledged observers: the
 * fluid comoving observer and the Eulerian observer, whose 4-velocity
 * is the unit future directed normal to the {\it t} = const hypersurfaces. 
 * The hydrodynamical quantities measured by the fluid observer correspond
 * to the members {\tt ent}, {\tt nbar}, {\tt ener}, and {\tt press}. 
 * The hydrodynamical quantities measured by the Eulerian observer correspond
 * to the members {\tt ener\_euler}, {\tt s\_euler}, {\tt gam\_euler}, and 
 * {\tt u\_euler}.  
 * 
 * @version #$Id$#
 */
class Star {

    // Data : 
    // -----
    protected:
	Map& mp ;	    /// Mapping associated with the star

	/// Number of domains of {\tt *mp} occupied by the star
	int nzet ;
	
	const Eos& eos ;   /// Equation of state of the stellar matter
    
	// Fluid quantities with respect to the fluid frame
	// ------------------------------------------------

	Scalar ent ;	   /// Log-enthalpy 
 	
	Scalar nbar ; 	   /// Baryon density in the fluid frame
	Scalar ener ;	   /// Total energy density in the fluid frame
	Scalar press ;	   /// Fluid pressure

	// Fluid quantities with respect to the Eulerian frame
	// ---------------------------------------------------
	Scalar ener_euler ; /// Total energy density in the Eulerian frame

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
	
	/** Total of the logarithm of the part of the lapse {\it N} 
	 *   generated principaly by the star. In the Newtonian case, 
	 *   this is the Newtonian gravitational potential
	 *   (in units of $c^2$). 
	 */
	Scalar logn ;

	/// Total lapse function 
	Scalar nnn ; 
	
	/// Total shift vector
	Vector shift ;
	
	/* Total of the part of $q = \psi^2 * N$
	 *   generated principaly by the star. 
	 */
	Scalar qq ;
	
	/// Total 3-metric 
	Metric gamma ;


    // Derived data : 
    // ------------
    protected:
	/// Coordinate radius at $\phi=0$, $\theta=\pi/2$. 
	mutable double* p_ray_eq ; 

	/// Coordinate radius at $\phi=\pi/2$, $\theta=\pi/2$. 
	mutable double* p_ray_eq_pis2 ;

	/// Coordinate radius at $\phi=\pi$, $\theta=\pi/2$. 
	mutable double* p_ray_eq_pi ;
	
	/// Coordinate radius at $\theta=0$. 
	mutable double* p_ray_pole ;
	
	/** Description of the stellar surface: 2-D {\tt Itbl} containing the 
	 *	values of the domain index {\it l} on the surface at the 
	 *	collocation points in $(\theta', \phi')$
	 */
	mutable Itbl* p_l_surf ; 
	
	/** Description of the stellar surface: 2-D {\tt Tbl} containing the 
	 *	values of the radial coordinate $\xi$ on the surface at the 
	 *	collocation points in $(\theta', \phi')$
	 */
	mutable Tbl* p_xi_surf ; 
	
	mutable double* p_mass_b ;	/// Baryon mass
	mutable double* p_mass_g ;	/// Gravitational mass


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
	
	
	Star(const Star& ) ;		/// Copy constructor

	/** Constructor from a file (see {\tt sauve(FILE* )}). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param fich	input file (must have been created by the function
	 *	{\tt sauve})
	 */
	Star(Map& mp_i, const Eos& eos_i, FILE* fich) ;    		

	virtual ~Star() ;			/// Destructor

	
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
	/// Assignment to another {\tt Star}
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
	 *  observer from those in the fluid frame ({\tt nbar}, {\tt ener}
	 *  and {\tt press}).
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

	/** Returns the logarithm of the total lapse {\it N}.
	 *   In the Newtonian case, this is the Newtonian
	 *   gravitational potential (in units of $c^2$). 
	 */
	const Scalar& get_logn() const {return logn;} ;

	/// Returns the total lapse function {\it N}
	const Scalar& get_nnn() const {return nnn;} ;

	/// Returns the total shift vector $N^i$.
	const Vector& get_shift() const {return shift;} ;

	/// Returns the total scalar field Q.
	const Scalar& get_qq() const {return qq;} ;

	/// Returns the total metric $\gamma$.
	const Metric& get_gamma() const {return gamma;} ;

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	    /// Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Star& ) ;	

    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

    // Global quantities
    // -----------------
    public:
	/// Coordinate radius at $\phi=0$, $\theta=\pi/2$ [r\_unit].
	double ray_eq() const ; 
	
	/// Coordinate radius at $\phi=\pi/2$, $\theta=\pi/2$ [r\_unit].
	double ray_eq_pis2() const ; 
	
	/// Coordinate radius at $\phi=\pi$, $\theta=\pi/2$ [r\_unit].
	double ray_eq_pi() const ; 
	
	/// Coordinate radius at $\theta=0$ [r\_unit]. 
	double ray_pole() const ; 
    
	/** Description of the stellar surface: returns a 2-D {\tt Itbl} 
	 *	containing the 
	 *	values of the domain index {\it l} on the surface at the 
	 *	collocation points in $(\theta', \phi')$.
	 *	The stellar surface is defined as the location where
	 *	the enthalpy (member {\tt ent}) vanishes.
	 */
	virtual const Itbl& l_surf() const ; 
	
	/** Description of the stellar surface: returns a 2-D {\tt Tbl} 
	 *	containing the values of the radial coordinate $\xi$ 
	 *	on the surface at the 
	 *	collocation points in $(\theta', \phi')$. 
	 *	The stellar surface is defined as the location where
	 *	the enthalpy (member {\tt ent}) vanishes.
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
 * Class for stars in binary system. *** UNDER DEEVELOPMENT ***
 *
 * A {\tt Star\_bin} can be construted in two states, represented by
 * the {\tt bool} member {\tt irrotational}: (i) irrotational
 * (i.e. the fluid motion is irrotational) or (ii) rigidly corotating 
 * with respect to the orbital motion (synchronized binary). 
 *
 * @version #$Id$#
 */
class Star_bin : public Star {

    // Data : 
    // -----
    protected:
	/** {\tt true} for an irrotational star, {\tt false} for a
	 *  corotating one
	 */
	bool irrotational ; 
	
	/** Scalar potential $\Psi_0$ of the non-translational part of
	 *  fluid 4-velocity (in the irrotational case)
	 */
	Scalar psi0 ; 

	/** Gradient of $\Psi$ (in the irrotational case)
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

	/** 3-vector shift, divided by {\it N}, of the rotating coordinates,
	 *  $\beta^i/N$. 
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
	
	/** Scalar field $ Q = \psi^2 N $ generated principaly by the
	 *  star.
	 */
	Scalar qq_auto ;

	/** Scalar field $ Q = \psi^2 N $ generated principaly by the
	 *  companion star.
	 */
	Scalar qq_comp ;

	
	/// Conformal factor $\psi^4$
	Scalar psi4 ;

	/// Covariant derivative of the logarithm of the conformal factor
	Vector dcov_lnpsi ;
	/// Contravariant derivative of the logarithm of the conformal factor
	Vector dcon_lnpsi ;

	/** Flat metric defined on the mapping (Spherical components 
	 * with respect to the mapping of the star }.
	 */
	Metric_flat flat ;

	/// Conformal metric $\tilde \gamma_{ij}$
	Metric gtilde ;

	/** Part of the shift vector generated principaly by the star
	 *  (Spherical components with respect to the mapping of the star)
	 */
	Vector shift_auto ; 

	/** Part of the shift vector generated principaly by the star
	 *  (Spherical components with respect to the mapping of the star)
	 */
	Vector shift_comp ; 


	/** Total deviation of the inverse conformal metric 
	 *  $\tilde \gamma^{ij}$ from the inverse flat metric. 
	 */
	Sym_tensor hij ; 


	/** Deviation of the inverse conformal metric 
	 *  $\tilde \gamma^{ij}$ from the inverse flat metric generated
	 *  principaly by the star. 
	 */
	Sym_tensor hij_auto ;

	/** Deviation of the inverse conformal metric 
	 *  $\tilde \gamma^{ij}$ from the inverse flat metric generated
	 *  principaly by the companion star. 
	 */
	Sym_tensor hij_comp ;

	/** Part of the extrinsic curvature tensor $\tilde K^{ij}$
	 *  generated by {\tt shift\_auto}. 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	Sym_tensor tkij_auto ;
	
	/** Part of the extrinsic curvature tensor $\tilde K^{ij}$
	 *  generated by {\tt shift\_comp}. 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	Sym_tensor tkij_comp ;
	
	/** Part of the scalar $K_{ij} K^{ij}$
	 *  generated by {\tt shift\_auto}, i.e. 
	 *    $K_{ij}^{\rm auto} K^{ij}_{\rm auto}$ 
	 */
	Scalar kcar_auto ;
	
	/** Part of the scalar $K_{ij} K^{ij}$
	 *  generated by {\tt shift\_auto} and {\tt shift\_comp}, i.e. 
	 *    $K_{ij}^{\rm auto} K^{ij}_{\rm comp}$ 
	 */
	Scalar kcar_comp ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt logn\_auto}.
	 */
	Scalar ssjm1_logn ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt qq\_auto}.
	 */
	Scalar ssjm1_qq ;


	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt phi}. (first scalar equation 
	 *  for the resolution of the vectorial poisson equation for the shift)
	 */
	Scalar ssjm1_phi ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt khi}. (second scalar equation 
	 *  for the resolution of the vectorial poisson equation for the shift)
	 */
	Scalar ssjm1_khi ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt mu}. (third scalar equation 
	 *  for the resolution of the vectorial poisson equation for the shift)
	 */
	Scalar ssjm1_mu ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt h00\_auto}.
	 */
	Scalar ssjm1_h11 ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt h10\_auto}.
	 */
	Scalar ssjm1_h21 ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt h20\_auto}.
	 */
	Scalar ssjm1_h31 ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt h11\_auto}.
	 */
	Scalar ssjm1_h22 ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt h21\_auto}.
	 */
	Scalar ssjm1_h32 ;

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt h22\_auto}.
	 */
	Scalar ssjm1_h33 ;

	/**
	  * Function used to construct the part $qq_auto$ generated by 
	  * the star  from the total $qq$.
	  * Mainly this {\tt Scalar} is 1 around the star and 0 around 
	  * the companion
	  * and the sum of {tt decouple} for the star and his companion is 1 
	  * everywhere.
	  */
	 Scalar decouple ;

	 /** {\tt true } if the 3-metric is conformally flat, {\tt false}
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
	 * @param irrot should be {\tt true} for an irrotational star, 
	 *		    {\tt false} for a corotating one
	 * @param conf_flat should be {\tt true} for a conformally flat metric
	 *                  {\tt false} for a general one
	 */
	Star_bin(Map& mp_i, int nzet_i,  const Eos& eos_i,
		   bool irrot, bool conf_flat) ;			
	
	
	Star_bin(const Star_bin& ) ;		/// Copy constructor

	/** Constructor from a file (see {\tt sauve(FILE* )}). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param fich	input file (must have been created by the function
	 *	{\tt sauve})
	 */
	Star_bin(Map& mp_i, const Eos& eos_i, FILE* fich) ;
    		
	virtual ~Star_bin() ;			/// Destructor


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
	/// Assignment to another {\tt Star\_bin}
	void operator=(const Star_bin& ) ;	
	
	/// Read/write the centrifugal potential
	Scalar& set_pot_centri() ;
	
        /** Read/write of the logarithm of the lapse generated 
	 *  principaly by the companion.
	 */ 
	Scalar& set_logn_comp() ;

	/// Read/write of $shift\_auto$
	Vector& set_shift_auto() ;	

	/// Read/write of $shift$
	Vector& set_shift() ;	
	
    // Accessors
    // ---------
    public:
	/** Returns {\tt true} for an irrotational motion, {\tt false} for 
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

	/** Returns the shift vector, divided by {\it N}, of the rotating 
	 *   coordinates, $\beta^i/N$. 
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

	/** Returns the part of the shift vector $\beta^i$ generated 
	 *  principaly by the star.
	 *  (Spherical components with respect to the mapping of the star)
	 */
	const Vector& get_shift_auto() const {return shift_auto;} ;

	/** Returns the part of the shift vector $\beta^i$ generated 
	 *  principaly by the star.
	 *  (Spherical components with respect to the mapping of the star)
	 */
	const Vector& get_shift_comp() const {return shift_comp;} ;

	/** Returns the part of the vector field $Q$ generated principaly 
	 *   by the star. 
	 */
	const Scalar& get_qq_auto() const {return qq_auto;} ;

	/** Returns the part of the vector field $Q$ generated principaly 
	 *   by the companion star. 
	 */
	const Scalar& get_qq_comp() const {return qq_comp;} ;

	/// Return the conformal factor $\psi^4$
	const Scalar& get_psi4() const {return psi4;} ;

	/** Return the flat metric defined on the mapping (Spherical
	 *  components  with respect to the mapping of the star)
	 */
	const Metric& get_flat() const {return flat;} ;	

	/// Return the conformal 3-metric $\tilde \gamma$
	const Metric& get_gtilde() const {return gtilde;} ;

	
	/** Return the total deviation of the inverse conformal metric 
	 *  $\tilde \gamma^{ij}$ from the inverse flat metric. 
	 */
	const Sym_tensor& get_hij() const {return hij;} ;

	/** Return the deviation of the inverse conformal metric 
	 *  $\tilde \gamma^{ij}$ from the inverse flat metric principaly
	 *  generated by the star.
	 */
	const Sym_tensor& get_hij_auto() const {return hij_auto;} ;

	/** Return the deviation of the inverse conformal metric 
	 *  $\tilde \gamma^{ij}$ from the inverse flat metric generated
	 *  principaly by the companion star.
	 */
	const Sym_tensor& get_hij_comp() const {return hij_comp;} ;


	/** Returns the part of the extrinsic curvature tensor 
	 *  $\tilde K^{ij}$ generated by {\tt shift\_auto}. 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	const Sym_tensor& get_tkij_auto() const {return tkij_auto;} ;

	/** Returns the part of the extrinsic curvature tensor 
	 *  $\tilde K^{ij}$ generated by {\tt shift\_comp}. 
	 *  (Spherical components with respect to the mapping of the star)
	 */
	const Sym_tensor& get_tkij_comp() const {return tkij_comp;} ;

	/**
	 * Returns the function used to construct {\tt qq\_auto} 
	 from {qq}.
	*/
	const Scalar get_decouple() const {return decouple ;}

	bool is_conf_flat() const {return conf_flat; } ; 

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	    /// Save in a file
    
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
	/** Performs the scalar product of two tensors by contracting
	 *  the last index of {\tt t1} with the first index of {\tt t2}.
	 */
	Tensor sprod(const Tensor& t1, const Tensor& t2) const ; 

	/** Computes the hydrodynamical quantities relative to the Eulerian
	 *  observer from those in the fluid frame, as well as 
	 *  {\tt wit\_w} and {\tt loggam}.  
	 *
	 *  The calculation is performed starting from the quantities
	 *  {\tt ent}, {\tt ener}, {\tt press}, {\tt a\_car} and {\tt bsn},  
	 *  which are supposed to be up to date.  
	 *  From these,  the following fields are updated:
	 *  {\tt gam\_euler}, {\tt u\_euler}, {\tt ener\_euler}, 
	 *  {\tt s\_euler}, {\tt stress\_euler},  
	 *  {\tt wit\_w} and {\tt loggam}. 
	 * 
	 */
	virtual void hydro_euler() ; 
	
	/** Computes metric coefficients from known potentials,
	 * when the companion is another star.
	 *
	 *  The calculation is performed starting from the quantities
	 *  {\tt logn\_auto},  {\tt qq\_auto}, {\tt shift\_auto}, 
	 *  {\tt hij\_auto}, {\tt comp.logn\_auto},  {\tt comp.qq\_auto},
	 *  {\tt comp.shift\_auto}, {\tt comp.hij\_auto}
	 *  which are supposed to be up to date.
	 *  From these,  the following fields are updated:
	 *  {\tt logn\_comp}, {\tt qq\_comp}, {\tt shift\_comp},
	 *  {\tt hij\_comp}, {\tt nnn},  {\tt psi4},  {\tt shift},
	 *  
	 *  @param comp companion star.
	 *
	 */
	void update_metric(const Star_bin& comp) ;
	
	/** Function used to update {\tt qq\_auto} and {\tt logn\_auto}
	 *  at the beginning of coal. 
	 */
	void update_metric_init1(const Star_bin& comp) ;

	/** Function used to update {\tt qq\_auto} and {\tt logn\_auto}
	 *  at the beginning of coal. 
	 */
	void update_metric_init2(const Star_bin& comp) ;

	/** Function used to update {\tt qq\_auto} and {\tt logn\_auto}
	 *  at the beginning of coal. 
	 */
	void update_metric_init3() ;

	/** Same as {\tt update\_metric(const Star\_bin\& )} but with
	 *  relaxation.
	 *
	 *  @param comp companion star.
	 *  @param star_prev previous value of the star. 
	 *  @param relax relaxation parameter. 
	 * 
	 */
	void update_metric(const Star_bin& comp, const Star_bin& star_prev, 
			   double relax) ; 
	
	/** Computes the derivative of metric functions related to the
	 *  companion star.
	 */
	 void update_metric_der_comp(const Star_bin& comp) ;

	/** Computes the quantities {\tt bsn} and {\tt pot\_centri}.
	 * 
	 *  The calculation is performed starting from the quantities
	 *  {\tt nnn}, {\tt shift},  {\tt Q},  
	 *  which are supposed to be up to date.  
	 * 
	 *  @param omega  angular velocity with respect to an asymptotically 
	 *		  inertial observer
	 *  @param x_axe  absolute X coordinate of the rotation axis
	 */
	void kinematics(double omega, double x_axe) ; 
	
	/** Computes the gradient of the total velocity potential $\psi$. 
	 * 
	 */
	void fait_d_psi() ; 
	
	/** Computes {\tt tkij\_auto} and {\tt akcar\_auto} from 
	 *  {\tt shift\_auto}, {\tt nnn} and {\tt Q}.
	 */
	void extrinsic_curvature() ; 
		
	
	/** Computes an equilibrium configuration.
	 * 
	 *  @param ent_c  [input] Central enthalpy
	 *  @param mermax [input] Maximum number of steps 
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    poisson scalar
	 *  @param relax_poisson [input]  Relaxation factor in poisson scalar
	 *  @param mermax_potvit [input]  Maximum number of steps in 
	 *				  Map\_radial::poisson\_compact
	 *  @param relax_potvit [input]   Relaxation factor in 
	 *				  Map\_radial::poisson\_compact
	 *  @param thres_adapt  [input]   Threshold on dH/dr for the adaptation 
	 *				  of the mapping
	 *  @param fact [input]    1-D {\tt Tbl} for the input
	 *                          of some factors : \\
	 *          {\tt fact(0)} : A resizing factor for the first shell
	 *  @param diff [output]   1-D {\tt Tbl} for the storage of some
	 *			    error indicators : \\
	 */
	void equilibrium(double ent_c, int mermax, int mermax_potvit, 
			 int mermax_poisson, double relax_poisson, 
			 double relax_potvit, double thres_adapt, 
			 const Tbl& fact, Tbl& diff) ;


	/** Computes the non-translational part of the velocity scalar potential
	 *  $\psi0$ by solving the continuity equation.
	 *  
	 *  @param mermax  [input] Maximum number of steps in the iteration
	 *  @param precis  [input] Required precision: the iteration will
	 *			   be stopped when the relative difference
	 *			   on $\psi0$ between two successive steps
	 *			   is lower than {\tt precis}.
	 *  @param relax   [input] Relaxation factor.  
	 *
	 *  @return Relative error of the resolution obtained by comparing
	 *	    the operator acting on the solution with the source.
	 */
	 double velocity_potential(int mermax, double precis, double relax) ;

	/** Performs a relaxation on {\tt ent}, {\tt logn\_auto},
	 *  {\tt qq\_auto}, {\tt shift\_auto}  and {\tt hij\_auto}. 
	 * 
	 *  @param star_prev   [input] star at the previous step.
	 *  @param relax_ent   [input] Relaxation factor for {\tt ent} 
	 *  @param relax_met   [input] Relaxation factor for {\tt logn\_auto},
	 *			       {\tt qq\_auto}, {\tt shift\_auto}, 
	 *			       only if {\tt (mer \% fmer\_met == 0)}.
	 *  @param mer	       [input] Step number
	 *  @param fmer_met    [input] Step interval between metric updates
	 */
	void relaxation(const Star_bin& star_prev, double relax_ent, 
			double relax_met, int mer, int fmer_met) ;
	

	friend class Binary ;

};

#endif
