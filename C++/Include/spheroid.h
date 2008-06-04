/*
 *  Definition of Lorene class Spheroid
 *
 */

/*
 *   Copyright (c) 2006  Nicolas Vasset, Jerome Novak & Jose-Luis Jaramillo
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

#ifndef __SPHEROID_H_ 
#define __SPHEROID_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.9  2008/06/04 12:31:56  n_vasset
 * new functions multipole_mass and multipole_angu.
 *
 * Revision 1.8  2008/06/03 14:32:23  n_vasset
 * dzpuis corrections (provisory). New function mass implemented.
 *
 * Revision 1.7  2006/09/07 08:39:43  j_novak
 * Minor changes.
 *
 * Revision 1.6  2006/07/03 10:12:51  n_vasset
 * adding of flag isphere
 *
 * Revision 1.5  2006/06/13 08:07:55  n_vasset
 *   Addition of 2 members
 *
 * Revision 1.4  2006/06/08 09:50:49  n_vasset
 * corrected version for coordinates transformation
 *
 * Revision 1.3  2006/06/07 14:32:24  n_vasset
 * modified constructor for 3d projector (Sym_tensor to Tensor)
 *
 * Revision 1.2  2006/05/26 13:25:42  j_novak
 * Removed a 'const'
 *
 * Revision 1.1  2006/05/26 13:20:42  j_novak
 * New class spheroid.
 *
 *
 *
 * $Header$
 *
 */
#include "metric.h"

/**
 * Spheroidal 2-surfaces embedded in a time-slice of the 3+1 formalism.\ingroup (star)
 * 
 */
class Spheroid {


    // Data : 
    // -----
    protected:
        ///The location of the 2-surface as \e r = \c h_surf\f$(\theta, \varphi)\f$
        Scalar h_surf ; 

        /** The jacobian of the adaptation of coordinates
	    (contravariant/covariant representation)*/

	Tensor jac2d ;

 	/** The 3-d projector on the 2-surface (contravariant-covariant form).
	 **/
    	Tensor proj ;

	/** The 3-d covariant degenerated 2-metric on the surface*/
 
	Sym_tensor qq ; 
   
        /** The adapted normal vector field to spheroid in the 3-slice */
 
        Vector ss ;

	Metric qab ;  /// Induced metric on the 2-surface \f$ q_{ab} \f$
        /** Extrinsic curvature of the 2-surface in the 3-slice.
	 * \f$ H_{ab} \f$ (covariant representation)
	 */
	Sym_tensor hh ; 

	Scalar trk ; ///< Trace \e K of the extrinsic curvature of the 3-slice

	/** Normal-tangent component of the extrinsic curvature of the 3-slice.
	 *  \f$ L_a \f$ (covariant representation)
	 */
	Vector ll ;  

	/** Tangent components of the extrinsic curvature of the 3-slice.
	 * \f$ J_{ab} \f$ (covariant representation)
	 */
	Sym_tensor jj ; 

        /** Normalization function for the outgoing null vector \b l. 
	 */ 
	Scalar fff;

	/** Normalization function for the ingoing null vector \b k.
	 */

	Scalar ggg;

	/** Flag to know whether the horizon is geometrically round or distorted*/
 
        bool issphere ;  


    // Derived data : 
    // ------------
    protected:
	mutable Scalar* p_sqrt_q ; ///< Surface element \f$\sqrt{\det q_{ab}} \f$
	mutable double* p_area ;   ///< The area of the 2-surface
	mutable double* p_angu_mom ; ///< The angular momentum
	mutable double* p_mass ; ///< Mass defined from angular momentum
	mutable double* p_multipole_mass ;///< Mass multipole for the spheroid.
	mutable double* p_multipole_angu ;///< Angular momentum multipole for the spheroid.
	mutable Scalar* p_theta_plus ; ///< Null outgoing expansion
	mutable Scalar* p_theta_minus ; ///< Null ingoing expansion
	mutable Sym_tensor* p_shear ; ///< The shear tensor

        mutable Sym_tensor* p_ricci ; ///< The ricci tensor for the spheroid
        mutable Tensor* p_delta ; 



    // Constructors - Destructor
    // -------------------------
    public:
        /** Standard constructor.
	 * The input mapping must be defined on a mono-domain angular grid
	 * (see \c Mg3d::get_angu_mono_domain() for details).
	 */
	Spheroid(const Map_af& map, double radius) ;  
	/** Constructor of a spheroid embedded in a 3-slice (\c Time_slice ) of 3+1 
	 * formalism. 
	 * This is done from the \c Time_slice data.
	 * @param h_in : the location of the surface \e r = h_in (WARNING:must be 
	                    defined on a mono-domain angular grid)
	 * @param gamij : the 3-metric on the 3-slice
	 * @param Kij : the extrinsic curvature of the 3-slice 
	 *                 (covariant representation)
	 */
	Spheroid(const Scalar& h_in, const Metric& gamij, const Sym_tensor& Kij) ;
	Spheroid(const Spheroid& ) ;		///< Copy constructor

	/// Constructor from a file (see \c sauve(FILE*) )
	Spheroid(FILE* ) ;    		

	virtual ~Spheroid() ;			///< Destructor
 

    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Spheroid
	void operator=(const Spheroid&) ;	
	
    // Accessors
    // ---------
    public:
	/// Returns the field \c h_surf
	const Scalar& get_hsurf() const {return h_surf ; } ;

	/// Returns the metric \f$ q_{ab} \f$
	const Metric& get_qab() const {return qab ; } ;

	/// Returns the symmetric tensor \f$ H_{ab} \f$
	const Sym_tensor& get_hh() const {return hh ; } ;

	/// returns the 3-d degenerate 2-metric \f$ Q_{ab}\f$
	const Sym_tensor& get_qq() const {return qq ; } ;

	/// returns the 3-d projector on 2-surface \f$ Pi \f$
	const Tensor& get_proj() const {return proj ; } ;

	/// returns the 2-d jacobian of coordinate transformation \f$ J \f$
	const Tensor& get_jac2d() const {return jac2d ; } ;

	/// Returns the trace \c K on the 2-surface
	const Scalar& get_trk() const {return trk; } ;

	/// Returns the vector \f$ L_a \f$
	const Vector& get_ll() const {return ll;} ;

	/// Returns the vector \f$ S_a \f$
	const Vector& get_ss() const {return ss;} ;

	///Returns the symmetric tensor \f$ J_{ab} \f$
	const Sym_tensor& get_jj() const {return jj;} ;

       	///Returns the normalization scalar \c F
	const Scalar& get_fff() const {return fff;} ;
    
	///Returns the normalization scalar \c G
	const Scalar& get_ggg() const {return ggg;} ;

        /// Returns the flag saying whether or not the horizon is geometrically round
	bool get_issphere() const{return issphere;}; 

	/// Sets the field \c h_surf
	Scalar& set_hsurf() {del_deriv() ; return h_surf ; } ;

	/// Sets the modified metric (non degenerated) \f$ q_{ab} \f$
	Metric& set_qab() {del_deriv() ; return qab ; } ;

	/// Sets the degenerated metric \f$ Q_{ab} \f$
	Sym_tensor& set_qq() {del_deriv() ; return qq ; } ;

	/// Sets the projector \f$ Pi  \f$
	Tensor& set_proj() {del_deriv() ; return proj ; } ;

	/// Sets the symmetric tensor \f$ H_{ab} \f$
	Sym_tensor& set_hh() {del_deriv() ; return hh ; } ;

	/// Sets the trace \c K on the 2-surface
	Scalar& set_trk() {del_deriv() ; return trk; } ;

	/// Sets the vector \f$ L_a \f$
	Vector& set_ll() {del_deriv() ; return ll;} ;

	/// Sets the vector \f$ s_a \f$
	Vector& set_ss() {del_deriv() ; return ss;} ;

	///Sets the symmetric tensor \f$ J_{ab} \f$
	Sym_tensor& set_jj() {del_deriv() ; return jj;} ;

       	///Sets the normalization factor \c F
	Scalar& set_fff() {del_deriv() ; return fff;} ;

      	///Sets the normalization factor \c G
	Scalar& set_ggg() {del_deriv() ; return ggg;} ;

        /// Sets the boolean linked to geometrical shape of the horizon
        bool set_issphere() {del_deriv() ; return issphere; };  
      
        /// Updates from the 3-slice data
        void update_from_tslice(const Metric& gamij, const Sym_tensor& Kij) ;
	
    // Computational functions
    // -----------------------
    public:

	/// Computes the normal vector field to the 2-surface
//	Vector compute_normal(const Metric& gamij) ;

	/// Computes the square root of the determinant of \f$ q_{ab} \f$.
	const Scalar& sqrt_q() const ;

	/** Computes the area of the 2-surface.
	 * This is defined as \f[
	 * {\cal A} = \int h^2 \sqrt{\det q_{ab}} \sin \theta {\rm d}\theta
	 * {\rm d}\varphi \f] */
	double area() const ;

	/** Computes the angular momentum with respect to a divergence-free vector field
	 * tangent to the 2-surface.
	 * This is defined as \f[
	 * {\cal J} = \int h^2 L_i \phi^i \sqrt{\det q_{ab}} \sin \theta {\rm d}\theta
	 * {\rm d}\varphi \f]
	 * @param  phi : the divergence-free vector field \f$ \phi \f$
	 */
	double angu_mom(const Vector& phi) const ;

	/** Computes the mass as defined from the calculus of angular momentum,
	 * done with respect to a divergence free tangent vector field \f$ phi \f$.
	 * Spheroid has to be a real sphere (flag issphere true), of constant radius \f[R_{s} \f].
	 * defined as \f[ M = \frac{1}{2 R_{s}} \sqrt{R_{s}^{4} + 4{\cal J}^{2}} \f]
         */

	double mass(const Vector& phi) const;

	/** Computes the mass multipole of a given order for the spheroid,
	 * assumed to be spherical.
	 */

	double multipole_mass(const int order, const Vector& phi) const;


	/** Computes the angular multipole of a given order for the spheroid,
	 * assumed to be spherical. \f$ phi \f$ is a divergence free tangent vector field.
	 * order has to be strictly higher than zero (no topological defects here...)
	 */

	double multipole_angu(const int order, const Vector& phi) const;

	/// Computes the outgoing null expansion \f$ \theta_+ \f$.
        const Scalar& theta_plus() const ;

	/// Computes the ingoing null expansion \f$ \theta_- \f$.
	const Scalar& theta_minus() const ;

	/// Computes the shear of the 2-surface \f$ \sigma_{ab} \f$.
	const Sym_tensor& shear () const ;

        ///Computes the Ricci tensor of the 2-surface
 
	const Sym_tensor& ricci () const; 

	/// Computes the round covariant derivative on the spheroid 
 
	Tensor derive_cov2dflat(const Tensor& uu) const ; 

	///Computes the delta coefficients for covariant derivative
  
	const Tensor& delta() const; 

	/// Computes the total covariant derivative on the spheroid 
 
	Tensor derive_cov2d(const Tensor& uu) const ; 

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Spheroid& ) ;	

};

ostream& operator<<(ostream& , const Spheroid& ) ;

 //------------------
 //  Class App_hor
 //------------------

/**
 * Class describing an apparent horizon. It is derived from the class \c Spheroid,
 * with the property that the outgoing null expansion \f$ \theta_+ = 0\f$.
 */
class App_hor : public Spheroid {

    // Constructors - Destructor
    // -------------------------
    public:
	App_hor(const Mg3d& grid_angu, double radius) ;  ///< Standard constructor
	/** Constructor of an apparent horizon embedded in a 3-slice (\c Time_slice ) 
	 * of 3+1 formalism. 
	 * This is done from the \c Time_slice data.
	 * @param h_in : the location of the surface \e r = h_in
	 * @param gamij : the 3-metric on the 3-slice
	 * @param Kij : the extrinsic curvature of the 3-slice 
	 *                 (covariant representation)
	 */
	App_hor(const Scalar& h_in, const Metric& gamij, const Sym_tensor& Kij) ;
	App_hor(const App_hor& ) ;		///< Copy constructor

	/// Constructor from a file (see \c sauve(FILE*) )
	App_hor(FILE* ) ;    		
       
	virtual ~App_hor() ;			///< Destructor
 
    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another App_hor
	void operator=(const App_hor&) ;

	bool check_expansion(double thres = 1.e-7) const ;

        /// Lie derivative of shear tensor with respect to the evolution vector field.
	const Sym_tensor& lie_derive_shear(const Scalar& bb, const Scalar& lapse) ;

        /** Lie derivative of the null outgoing expansion rate with respect to the
	 * evolution vector field.
	 */
	const Sym_tensor& lie_derive_theta_plus(const Scalar& bb, const Scalar& lapse);

        /** Lie derivative of the null ingoing expansion rate with respect to the
	 * evolution vector field.
	 */
	const Sym_tensor& lie_derive_theta_minus(const Scalar& bb, const Scalar& lapse);

        /// Lie derivative of 2-metric with respect to the evolution vector field.
	const Sym_tensor& lie_derive_q_ab(const Scalar& bb, const Scalar& lapse) ;

        /** non-affinity (or surface gravity) with respect to the outgoing null vector
	 * field
	 */
	const Scalar& l_non_affinity(const Scalar& bb, const Scalar& lapse) ;
        
} ;
	

#endif
