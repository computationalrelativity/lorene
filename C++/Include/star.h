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
 * Base class for stars.
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
 *   ds^2 = - (N^2 - N_i N^i) dt^2 - 2 N_i \,  dt\,  dx^i
 *	    + A^2  \,  {\tilde h}_{ij} \,  dx^i dx^j
 * \end{equation}
 * where ${\tilde h}_{ij}$ is a 3-metric.
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

	/** Number of domains of {\tt *mp} occupied by the star
	 * 
	 */
	int nzet ;
	
	const Eos& eos ;   /// Equation of state of the stellar matter
    
	// Fluid quantities with respect to the fluid frame
	// ------------------------------------------------
	/// Log-enthalpy 
	Scalar ent ;	  
	
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
	
	/* Total of the part of $q = psi^2 * N$
	 *   generated principaly by the star. 
	 */
	Scalar qq ;
	
	/// Total 3-metric gamma
	Metric gamij ;


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
	 * @param relat should be {\tt true} for a relativistic
	 *			star,  {\tt false} for a Newtonian one
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

	/** Returns the logarithm of the part of the lapse {\it N} generated 
	 *   principaly by the star.
	 *   In the Newtonian case, this is the Newtonian
	 *   gravitational potential (in units of $c^2$). 
	 */
	const Scalar& get_logn() const {return logn;} ;

	/// Returns the total lapse function {\it N}
	const Scalar& get_nnn() const {return nnn;} ;

	/// Returns the total shift vector $N^i$
	const Vector& get_shift() const {return shift;} ;

	const Scalar& get_qq() const {return qq;} ;
	const Metric& get_gamij() const {return gamij;} ;

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
    	virtual double mass_b() const ;
	
	/// Gravitational mass
    	virtual double mass_g() const ;
	
};
ostream& operator<<(ostream& , const Star& ) ;	


			//---------------------------//
			//    class Star_bin       //
			//---------------------------//

/**
 * Class for stars in binary system.
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
	 *  (Spherical components with respect to {\tt map\_triad})
	 */
	Vector d_psi ; 
	
	/** Spatial projection of the fluid 3-velocity with respect to  
	 *  the co-orbiting observer. 
	 *  (Spherical components with respect to {\tt map\_triad})
	 */
	Vector wit_w ; 
	
	/** Logarithm of the Lorentz factor between the fluid and 
	 *  the co-orbiting observer.
	 */
	Scalar loggam ; 

	/** 3-vector shift, divided by {\it N}, of the rotating coordinates,
	 *  $B^i/N$. 
	 *  (Spherical components with respect to {\tt map\_triad})
	 */
	Vector bsn ; 
	
	/// Centrifugal potential
	Scalar pot_centri ; 	

	Tensor stress ;

	/** Part of the lapse logarithm (gravitational potential at the
	 *  Newtonian limit) generated principaly by the star. 
	 */
	Scalar logn_auto ; 

	/** Part of the lapse logarithm (gravitational potential at the
	 *  Newtonian limit) generated principaly by the companion star. 
	 */
	Scalar logn_comp ; 

	Vector dcov_logn ;
	Vector dcon_logn ;
	
	/** Part of the shift vector $N^i$ generated principaly by the star.
	 *  (Spherical components with respect to {\tt map\_triad)
	 */
	Vector shift_auto ; 
	
	/** Part of the shift vector $N^i$ generated principaly by the 
	 *  companion star. 
	 *  (Spherical components with respect to {\tt map\_triad})
	 */
	Vector shift_comp ; 


	Scalar qq_auto ;
	Scalar qq_comp ;
	Scalar psi4 ;
	Vector dcov_lnpsi ;
	Vector dcon_lnpsi ;
	Metric_flat flat ;
	Metric gtilde ;
	Sym_tensor hij ; 
	Sym_tensor hij_auto ;
	Sym_tensor hij_comp ;

	/** Part of the extrinsic curvature tensor 
	 *  $\tilde K^{ij} = A^2 K^{ij}$
	 *  generated by {\tt shift\_auto}. 
	 *  (Spherical components with respect to {\tt map\_triad})
	 */
	Sym_tensor tkij_auto ;
	
	/** Part of the extrinsic curvature tensor 
	 *  $\tilde K^{ij} = A^2 K^{ij}$
	 *  generated by {\tt shift\_comp}. 
	 *  (Spherical components with respect to {\tt map\_triad})
	 */
	Sym_tensor tkij_comp ;
	
	/** Part of the scalar $A^2 K_{ij} K^{ij}$
	 *  generated by {\tt shift\_auto}, i.e. 
	 *    $A^2 K_{ij}^{\rm auto} K^{ij}_{\rm auto}$ 
	 */
	Scalar kcar_auto ;
	
	/** Part of the scalar $A^2 K_{ij} K^{ij}$
	 *  generated by {\tt shift\_auto} and {\tt shift\_comp}, i.e. 
	 *    $A^2 K_{ij}^{\rm auto} K^{ij}_{\rm comp}$ 
	 */
	Scalar kcar_comp ;

	 /**
	  * Function used to construct the part of $K^{ij}$ generated by 
	  * the star  from the total $K^{ij}$. Only used for a binary system
	  * where the other member is a black hole.
	  * 
	  * Mainly this {\tt Scalar} is 1 around the hole and 0 around the companion
	  * and the sum of {tt decouple} for the hole and his companion is 1 
	  * everywhere.
	  */
	 Scalar decouple ;

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
	 * @param conf_flat should be {\tt true} for a conformally flat star
	 *                  {\tt false} for a general one
	 * @param flat a flat metric 
	 * @param source the source for gtilde 
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
	
	/** Read/write the part of the lapse logarithm (gravitational potential 
	 *  at the Newtonian limit) generated principaly by the companion star. 
	 */
	Scalar& set_logn_comp() ;

	/// Read/write the centrifugal potential
	Scalar& set_pot_centri() ;
	
    // Accessors
    // ---------
    public:
	/** Returns {\tt true} for an irrotational motion, {\tt false} for 
	 *  a corotating one. 
	 */
	bool is_irrotational() const {return irrotational; } ; 

	/// Returns the non-translational part of the velocity potential
	const Scalar& get_psi0() const {return psi0;} ;

	/** Returns the gradient of the velocity potential 
	 *  (Spherical components with respect to {\tt ref\_triad})
	 */
	const Vector& get_d_psi() const {return d_psi;} ;

	/** Returns the spatial projection of the fluid 3-velocity with 
	 *  respect to the co-orbiting observer. 
	 *  (Spherical components with respect to {\tt ref\_triad})
	 */
	const Vector& get_wit_w() const {return wit_w;} ;

	/** Returns the logarithm of the Lorentz factor between the fluid and 
	 *  the co-orbiting observer.
	 */
	const Scalar& get_loggam() const {return loggam;} ;

	/** Returns the shift vector, divided by {\it N}, of the rotating 
	 *   coordinates, $B^i/N$. 
	 *  (Spherical components with respect to {\tt ref\_triad})
	 */
	const Vector& get_bsn() const {return bsn;} ;

	/// Returns the centrifugal potential
	const Scalar& get_pot_centri() const {return pot_centri;} ;

	const Tensor& get_stress() const {return stress;} ;

	/** Returns the part of the lapse logarithm (gravitational potential 
	 *  at the Newtonian limit) generated principaly by the star. 
	 */
	const Scalar& get_logn_auto() const {return logn_auto;} ;

	/** Returns the part of the lapse logarithm (gravitational potential 
	 *  at the Newtonian limit) generated principaly by the companion star.
 	 */
	const Scalar& get_logn_comp() const {return logn_comp;} ;

	/** Returns the part of the shift vector $N^i$ generated principaly 
	 * by the star.
	 *  (Spherical components with respect to {\tt ref\_triad})
	 */
	const Vector& get_shift_auto() const {return shift_auto;} ;

	/** Returns the part of the shift vector $N^i$ generated principaly 
	 *   by the companion star. 
	 *  (Spherical components with respect to {\tt ref\_triad})
	 */
	const Vector& get_shift_comp() const {return shift_comp;} ;

	const Scalar& get_qq_auto() const {return qq_auto;} ;
	const Scalar& get_qq_comp() const {return qq_comp;} ;
	const Scalar& get_psi4() const {return psi4;} ;
	const Metric& get_flat() const {return flat;} ;	
	const Metric& get_gtilde() const {return gtilde;} ;
	const Sym_tensor& get_hij() const {return hij;} ;
	const Sym_tensor& get_hij_auto() const {return hij_auto;} ;
	const Sym_tensor& get_hij_comp() const {return hij_comp;} ;


	/** Returns the part of the extrinsic curvature tensor 
	 *  $\tilde K^{ij} = A^2 K^{ij}$ generated by {\tt shift\_auto}. 
	 *  (Spherical components with respect to {\tt map\_triad})
	 */
	const Sym_tensor& get_tkij_auto() const {return tkij_auto;} ;

	/** Returns the part of the extrinsic curvature tensor 
	 *  $\tilde K^{ij} = A^2 K^{ij}$ generated by {\tt shift\_comp}. 
	 *  (Spherical components with respect to {\tt map\_triad})
	 */
	const Sym_tensor& get_tkij_comp() const {return tkij_comp;} ;


	/**
	 * Returns the function used to construct {\tt tkij\_auto} 
	 from {\tt tkij\_tot}.
	*/
	const Scalar get_decouple() const {return decouple ;}

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
	
	/** Absolute coordinate X of the barycenter of the baryon density, 
	 *  defined according to the formula
	 *  \begin{equation}
	 *    X_G := \int A^3 \Gamma_{\rm n} \,  n \,  X \, d^3x \ ,  
	 *  \end{equation}
	 *  where $\Gamma_{\rm n}$ is the Lorentz factor between the fluid 
	 *  and Eulerian observers.
	 */
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
	 *  {\tt gam\_euler}, {\tt u\_euler}, {\tt ener\_euler}, {\tt s\_euler}, 
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
	 *  {\tt comp.shift\_auto}, {\TT COMP.hij\_auto}
	 *  which are supposed to be up to date.
	 *  From these,  the following fields are updated:
	 *  {\tt logn\_comp}, {\tt beta\_comp}, {\tt shift\_comp},
	 *  {\tt nnn},  {\tt a\_car},  {\tt shift},
	 *  {\tt d\_logn\_auto}, {\tt d\_beta\_auto}, {\tt tkij\_auto},
	 *  {\tt akcar\_auto}.
	 *
	 *  @param comp companion star.
	 *
	 */
	void update_metric(const Star_bin& comp) ;
	void update_metric_init(const Star_bin& comp) ;

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
	 *  {\tt nnn}, {\tt shift},  {\tt qq\_car},  
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
	 *  {\tt shift\_auto}, {\tt nnn} and {\tt a\_car}.
	 */
	void extrinsic_curvature() ; 
		
	
	/** Computes an equilibrium configuration.
	 * 
	 *  The values of {\tt logn\_comp}, {\tt beta\_comp}, {\tt pot\_centri}
	 *  are held fixed during the iteration. 
	 *  
	 *  @param ent_c  [input] Central enthalpy
	 *  @param mermax [input] Maximum number of steps 
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    Map\_et::poisson
	 *  @param relax_poisson [input]  Relaxation factor in Map\_et::poisson
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
	 *	    {\tt diff(0)} : Relative change in the enthalpy field
	 *			      between two successive steps \\
	 *	    {\tt diff(1)} : Relative error returned by the routine
	 *				{\tt Star\_bin::velocity\_potential}  
	 *	    {\tt diff(2)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt logn\_auto} \\  
	 *	    {\tt diff(3)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt beta\_auto} \\  
	 *	    {\tt diff(4)} : Relative error in the resolution of the
	 *			    equation for {\tt shift\_auto} (x comp.) \\  
	 *	    {\tt diff(5)} : Relative error in the resolution of the
	 *			    equation for {\tt shift\_auto} (y comp.) \\  
	 *	    {\tt diff(6)} : Relative error in the resolution of the
	 *			    equation for {\tt shift\_auto} (z comp.)   
	 */
	void equilibrium(double ent_c, int mermax, int mermax_potvit, 
			 double relax_potvit, double thres_adapt, 
			 const Tbl& fact, Tbl& diff) ;


//	double velocity_potential(int mermax, double precis, double relax) ;

	/** Performs a relaxation on {\tt ent}, {\tt logn\_auto},
	 *  {\tt beta\_auto} and {\tt shift\_auto}. 
	 * 
	 *  @param star_prev   [input] star at the previous step.
	 *  @param relax_ent   [input] Relaxation factor for {\tt ent} 
	 *  @param relax_met   [input] Relaxation factor for {\tt logn\_auto},
	 *			       {\tt beta\_auto}, {\tt shift\_auto}, 
	 *			       only if {\tt (mer \% fmer\_met == 0)}.
	 *  @param mer	       [input] Step number
	 *  @param fmer_met    [input] Step interval between metric updates
	 */
	void relaxation(const Star_bin& star_prev, double relax_ent, 
			double relax_met, int mer, int fmer_met) ;
	

	friend class Binary ;

};

#endif
