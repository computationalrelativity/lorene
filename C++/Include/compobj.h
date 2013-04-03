/*
 *  Definition of Lorene class Compobj, Compobj_QI, Star_QI, Kerr_QI
 *
 */

/*
 *   Copyright (c) 2012, 2013 Claire Some, Eric Gourgoulhon
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

#ifndef __COMPOBJ_H_ 
#define __COMPOBJ_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.9  2013/04/03 12:08:57  e_gourgoulhon
 * Added member kk to Compobj; suppressed tkij
 *
 * Revision 1.8  2013/04/02 23:17:18  e_gourgoulhon
 * New class Kerr_QI
 *
 * Revision 1.7  2012/12/03 15:26:14  c_some
 * Added data member m2
 *
 * Revision 1.6  2012/11/22 16:02:18  c_some
 * *** empty log message ***
 *
 * Revision 1.5  2012/11/21 14:52:13  c_some
 * Documentation corrected
 *
 * Revision 1.4  2012/11/20 16:21:16  c_some
 * Added new class Star_QI
 *
 * Revision 1.3  2012/11/16 16:13:12  c_some
 * Added new class Compobj_QI
 *
 * Revision 1.2  2012/11/15 20:50:41  e_gourgoulhon
 * Corrected the documentation
 *
 * Revision 1.1  2012/11/15 16:20:51  c_some
 * New class Compobj
 *
 *
 * $Header$
 *
 */


// Headers Lorene
#include "tensor.h"
#include "metric.h"


			//---------------------------//
			//    base class Compobj     //
			//---------------------------//

/**
 * Base class for stationary compact objects (***under development***). 
 * \ingroup(compactobjects)
 *
 * A \c Compobj describes a single compact object (star or black hole), in a stationary state.
 * 
 * The spacetime metric is written according to the 3+1 formalism :
 * \f[
 *   ds^2 = - N^2  dt^2 + \gamma_{ij} ( dx^i + \beta^i dt )
 *               (dx^j + \beta^j dt )
 * \f]
 * where \f$\gamma_{ij}\f$ is the 3-metric, described by a Lorene object of class \c Metric. 
 * 
 * The total energy-momentum tensor is orthogonally split with respect to the Eulerian observer as follows:
 * \f[
 *	T_{\alpha\beta} = E n_\alpha n_\beta + P_\alpha n_\beta + n_\alpha P_\beta + S_{\alpha\beta}
 * \f]
 */
class Compobj {

    // Data : 
    // -----
    protected:
	/// Mapping describing the coordinate system (r,theta,phi) 
	Map& mp ;  

	/// Lapse function \e N .
	Scalar nn ; 
	
	/// Shift vector \f$\beta^i\f$
	Vector beta ;
	
 	/// 3-metric  \f$\gamma_{ij}\f$
	Metric gamma ;

	/// Total energy density \e E in the Eulerian frame 
	Scalar ener_euler ; 

	/// Total 3-momentum density \f$P^i\f$ in the Eulerian frame 
	Vector mom_euler ; 

    /// Stress tensor \f$S_{ij}\f$  with respect to the Eulerian observer
    Sym_tensor stress_euler ;

    /// Extrinsic curvature tensor \f$K_{ij}\f$  
    Sym_tensor kk ;

   // Derived data : 
    // ------------
    protected:
	mutable double* p_adm_mass ;	///< ADM mass 

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the object is defined
	 * 
	 */
	Compobj(Map& map_i) ;

	Compobj(const Compobj& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE* )). 
	 * 
	 * @param mp_i Mapping on which the object is defined
	 * @param fich	input file (must have been created by the function
	 *	\c sauve)
	 */
	Compobj(Map& map_i, FILE* ) ;    		

	virtual ~Compobj() ;			///< Destructor
 

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
	/// Assignment to another Compobj
	void operator=(const Compobj&) ;	
	
	/// Read/write of the mapping
	Map& set_mp() {return mp; } ; 


    // Accessors
    // ---------
    public:
	/// Returns the mapping
	const Map& get_mp() const {return mp; } ; 

	/// Returns the lapse function \e N .
	const Scalar& get_nn() const {return nn;} ;

	/// Returns the shift vector \f$\beta^i\f$.
	const Vector& get_beta() const {return beta;} ;
	
 	/// Returns the 3-metric \f$\gamma_{ij}\f$.
	const Metric& get_gamma() const {return gamma;} ;

	/// Returns the total energy density \e E in the Eulerian frame 
	const Scalar& get_ener_euler() const {return ener_euler;}  ; 

	/// Returns the total 3-momentum density \f$P^i\f$ in the Eulerian frame 
	const Vector& get_mom_euler() const {return mom_euler;} ; 

    /// Returns the stress tensor \f$S_{ij}\f$  with respect to the Eulerian observer
    const Sym_tensor& get_stress_euler() const {return stress_euler;} ;

    /// Returns the extrinsic curvature tensor \f$K_{ij}\f$ 
    const Sym_tensor& get_kk() const {return kk;} ;

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
    void gyoto_data(const char* file_name) const ; ///< Save in a file for GYOTO
    
	/// Display
	friend ostream& operator<<(ostream& , const Compobj& ) ;	

    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

    // Computational methods
    // ---------------------
    public:
    /// Computation of the extrinsic curvature 
    virtual void extrinsic_curvature() ; 
    
    
    // Global quantities
    // -----------------
    public:
	/// ADM mass (computed as a surface integral at spatial infinity)
    virtual double adm_mass() const ;
};


			//---------------------------//
			//    base class Compobj_QI  //
			//---------------------------//

/**
 * Base class for axisymmetric stationary compact objects in Quasi-Isotropic coordinates (***under development***). 
 * \ingroup(compactobjects)
 *
 * The metric is expressed in Quasi-Isotropic (QI) coordinates :
 * \f[
 *   ds^2 = - N^2 dt^2 + A^2 (dr^2 + r^2 d\theta^2)
 *		       + B^2 r^2 \sin^2\theta (d\varphi - N^\varphi dt)^2
 * \f]
 *
 * 
 */
class Compobj_QI : public Compobj {

    // Data : 
    // -----
    protected:

	/// Square of the metric factor \e A 
	Scalar a_car ; 

	/// Metric factor \e B 
	Scalar bbb ; 

	/// Square of the metric factor \e B 
	Scalar b_car ; 

	/// Metric coefficient \f$N^\varphi\f$
	Scalar nphi ; 

	/** Scalar \f$A^2 K_{ij} K^{ij}\f$.
	 *  For axisymmetric stars, this quantity is related to the 
	 *  derivatives of \f$N^\varphi\f$ by
	 * \f[
	 *	A^2 K_{ij} K^{ij} = {B^2 \over 2 N^2} \, r^2\sin^2\theta \,  
	 *    \left[ \left( {\partial N^\varphi \over \partial r} \right) ^2
	 *	    + {1\over r^2} \left( {\partial N^\varphi \over 
	 *		    \partial \theta} \right) ^2 \right] \ . 
	 * \f]
	 * In particular it is related to the quantities \f$k_1\f$ and \f$k_2\f$
	 * introduced by Eqs.~(3.7) and (3.8) of 
	 * Bonazzola et al. \a Astron. \a Astrophys. \b 278 , 421 (1993)
	 * by 
	 * \f[
	 *	A^2 K_{ij} K^{ij} = 2 A^2 (k_1^2 + k_2^2) \ . 
	 * \f]
	 */
	 Scalar ak_car ; 


   // Derived data : 
    // ------------
    protected:
	mutable double* p_angu_mom ;	///< Angular momentum 
	mutable double* p_r_isco ;	///< Circumferential radius of the ISCO
	mutable double* p_f_isco ;	///< Orbital frequency of the ISCO
	/// Specific energy of a particle on the ISCO 
	mutable double* p_espec_isco ;	
	/// Specific angular momentum of a particle on the ISCO
	mutable double* p_lspec_isco ;	

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the object is defined
	 * 
	 */
	Compobj_QI(Map& map_i) ;

	Compobj_QI(const Compobj_QI& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE* )). 
	 * 
	 * @param mp_i Mapping on which the object is defined
	 * @param fich	input file (must have been created by the function
	 *	\c sauve)
	 */
	Compobj_QI(Map& map_i, FILE* ) ;    		

	virtual ~Compobj_QI() ;			///< Destructor
 

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
	/// Assignment to another Compobj_QI
	void operator=(const Compobj_QI&) ;	
	

    // Accessors
    // ---------
    public:

	/// Returns the metric factor \e B 
	const Scalar& get_bbb() const {return bbb;} ; 

	/// Returns the square of the metric factor \e A 
	const Scalar& get_a_car() const {return a_car;} ; 

	/// Returns the square of the metric factor \e B 
	const Scalar& get_b_car() const {return b_car;} ; 

	/// Returns the metric coefficient \f$N^\varphi\f$
	const Scalar& get_nphi() const {return nphi;} ; 


	/** Returns the scalar \f$A^2 K_{ij} K^{ij}\f$.
	 *  For axisymmetric stars, this quantity is related to the 
	 *  derivatives of \f$N^\varphi\f$ by
	 * \f[
	 *	A^2 K_{ij} K^{ij} = {B^2 \over 2 N^2} \, r^2\sin^2\theta \,  
	 *    \left[ \left( {\partial N^\varphi \over \partial r} \right) ^2
	 *	    + {1\over r^2} \left( {\partial N^\varphi \over 
	 *		    \partial \theta} \right) ^2 \right] \ . 
	 * \f]
	 * In particular it is related to the quantities \f$k_1\f$ and \f$k_2\f$
	 * introduced by Eqs. (3.7) and (3.8) of 
	 * Bonazzola et al. \a Astron. \a Astrophys. \b 278 , 421 (1993)
	 * by 
	 * \f[
	 *	A^2 K_{ij} K^{ij} = 2 A^2 (k_1^2 + k_2^2) \ . 
	 * \f]
	 */
	 const Scalar& get_ak_car() const {return ak_car;} ; 

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

    // Global quantities
    // -----------------
    public:
	virtual double angu_mom() const ;	///< Angular momentum 

	/** Circumferential radius of the innermost stable circular orbit (ISCO).	
	 *
	 *  @param lmin index of the domain from which the ISCO is searched outwards ;
	 *  @param ost output stream to give details of the computation;
	 *		if set to 0x0 [default value], no details will be
	 *		given.
	 *  
	 */
	virtual double r_isco(int lmin, ostream* ost = 0x0) const ;	
 	
 	/// Orbital frequency at the innermost stable circular orbit (ISCO).	
 	virtual double f_isco(int lmin) const ;	

	/// Energy of a particle on the ISCO 
 	virtual double espec_isco(int lmin) const ;	
	
	/// Angular momentum of a particle on the ISCO
 	virtual double lspec_isco(int lmin) const ;	

    // Computational routines
    // ----------------------

	/** Updates the 3-metric \f$\gamma_{ij}\f$ from \e A and \e B 
	 *  and the shift vector \f$\beta^i\f$  from \f$N^\phi\f$. 
	 * 
	 */
	virtual void update_metric() ; 
		
	/** Computes the extrinsic curvature  and \c ak_car  from 
	 *  \c nphi , \c nn  and \c b_car .
	 */
	virtual void extrinsic_curvature() ;
	
};


			//--------------------------//
			//   base class Star_QI     //
			//--------------------------//

/**
 *Base class for axisymmetric stationary compact stars in Quasi-Isotropic coordinates (***under development***). 
 * \ingroup(compactobjects)
 *
 * The time slice \f$t=\mathrm{const}\f$ has the topology of \f$R^3\f$ 
 * and the metric is expressed in Quasi-Isotropic (QI) coordinates :
 * \f[
 *   ds^2 = - N^2 dt^2 + A^2 (dr^2 + r^2 d\theta^2)
 *		       + B^2 r^2 \sin^2\theta (d\varphi - N^\varphi dt)^2
 * \f]
 *
 * 
 */
class Star_QI : public Compobj_QI {

    // Data : 
    // -----
    protected:

	/** Logarithm of the lapse \e N .
	 */
	Scalar logn ;

	/** Component \f$\tilde N^\varphi = N^\varphi r\sin\theta\f$ of the
	 *  shift vector
	 */
	Scalar tnphi ; 

	/** Part of the Metric potential \f$\nu = \ln N\f$ = \c logn 
	 *  generated by the matter terms
	 */
	Scalar nuf ;	

	/** Part of the Metric potential \f$\nu = \ln N\f$ = \c logn 
	 *  generated by the quadratic terms
	 */
	Scalar nuq ;	

	/// Metric potential \f$\zeta = \ln(AN)\f$ 
	Scalar dzeta ;	

	/// Metric potential \f$\tilde G = (NB-1) r\sin\theta\f$
	Scalar tggg ; 

	/** Vector \f$W^i\f$ used in the decomposition of \c shift ,
	 *  following Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 * NB: \c w_shift  contains the components of \f$W^i\f$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping \c mp . 
	 */
	Vector w_shift ; 
	
	/** Scalar \f$\chi\f$ used in the decomposition of \c shift ,
	 *  following Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 */
	Scalar khi_shift ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c nuf  by means of
	 *  \c Map_et::poisson .
	 */
	Scalar ssjm1_nuf ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c nuq  by means of
	 *  \c Map_et::poisson .
	 */
	Scalar ssjm1_nuq ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c dzeta .
	 */
	Scalar ssjm1_dzeta ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c tggg .
	 */
	Scalar ssjm1_tggg ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for the scalar \f$\chi\f$ by means of
	 *  \c Map_et::poisson . 
	 *  \f$\chi\f$ is an intermediate quantity for the resolution of the
	 *  elliptic equation for the shift vector \f$N^i\f$
	 */
	 Scalar ssjm1_khi ; 
	 
	/** Effective source at the previous step for the resolution of 
	 *  the vector Poisson equation for \f$W^i\f$.
	 *  \f$W^i\f$ is an intermediate quantity for the resolution of the
	 *  elliptic equation for the shift vector \f$N^i\f$
	 *  (Components with respect to the Cartesian triad associated with 
	 *   the mapping \c mp )
	 */
	 Vector ssjm1_wshift ; 
	 

    // Derived data : 
    // ------------
    protected:
	
	mutable double* p_grv2 ;	///< Error on the virial identity GRV2
	mutable double* p_grv3 ;	///< Error on the virial identity GRV3
	mutable double* p_mom_quad ;	///< Quadrupole moment	
	mutable double* p_mass_g ;	///< Gravitational mass (ADM mass as a volume integral)
	 

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star is contructed
     *
	 */
	Star_QI(Map& mp_i) ;			
	
	
	Star_QI(const Star_QI& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE*) ). 
	 * 
	 * @param mp_i Mapping on which the star is constructed
	 * @param fich	input file (must have been created by the function
	 *	\c Star_QI::sauve )
	 */
	Star_QI(Map& mp_i, FILE* fich) ;    		

	virtual ~Star_QI() ;			///< Destructor


    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0  all the pointers on derived quantities
	virtual void set_der_0x0() const ; 

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Star_QI 
	void operator=(const Star_QI& ) ;	
	
    // Accessors
    // ---------
    public:

	/** Returns the logarithm of the lapse \e N.
	 */
	const Scalar& get_logn() const {return logn;} ;


	/** Returns the component \f$\tilde N^\varphi = N^\varphi r\sin\theta\f$ 
	 *  of the shift vector
	 */
	const Scalar& get_tnphi() const {return tnphi;} ; 
	
	/** Returns the part of the Metric potential \f$\nu = \ln N\f$ = \c logn 
	 *  generated by the matter terms
	 */
	const Scalar& get_nuf() const {return nuf;} ;	

	/** Returns the Part of the Metric potential \f$\nu = \ln N\f$ = \c logn 
	 *  generated by the quadratic terms
	 */
	const Scalar& get_nuq() const {return nuq;} ;	

	/// Returns the Metric potential \f$\zeta = \ln(AN)\f$ 
	const Scalar& get_dzeta() const {return dzeta;} ;	

	/// Returns the Metric potential \f$\tilde G = (NB-1) r\sin\theta\f$
	const Scalar& get_tggg() const {return tggg;} ; 

	/** Returns the vector \f$W^i\f$ used in the decomposition of 
	 *  \c shift ,
	 *  following Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 * NB: \c w_shift  contains the components of \f$W^i\f$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping \c mp . 
	 */
	const Vector& get_w_shift() const {return w_shift;} ; 
	
	/** Returns the scalar \f$\chi\f$ used in the decomposition of 
	 *  \c shift  
	 *  following Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 * NB: \c w_shift  contains the components of \f$W^i\f$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping \c mp . 
	 */
	const Scalar& get_khi_shift() const {return khi_shift;} ; 


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
		
	virtual double mass_g() const ;	    ///< Gravitational mass
	virtual double angu_mom() const ;	///< Angular momentum 

	/** Error on the virial identity GRV2.
	 *  This indicator is only valid for relativistic computations.
	 */
	virtual double grv2() const ;	

	/** Error on the virial identity GRV3.
	 *  The error is computed as the integral defined
	 *  by Eq. (43) of [Gourgoulhon and Bonazzola, 
	 *  \a Class. \a Quantum \a Grav. \b 11, 443 (1994)] divided by
	 *  the integral of the matter terms.
	 * 
	 *  @param ost output stream to give details of the computation;
	 *		if set to 0x0 [default value], no details will be
	 *		given.
	 *   
	 */
	virtual double grv3(ostream* ost = 0x0) const ;	
    
	/** Quadrupole moment.
	 *  The quadrupole moment \e Q is defined according to Eq. (7) of
	 *  [Salgado, Bonazzola, Gourgoulhon and Haensel, \a Astron. \a Astrophys.
	 *   \b 291 , 155 (1994)]. At the Newtonian limit it is related to
	 *  the component \f${\bar I}_{zz}\f$ of the MTW (1973) reduced quadrupole 
	 *  moment \f${\bar I}_{ij}\f$ by: \f$Q = -3/2 {\bar I}_{zz}\f$. 
	 *  Note that \e Q is the negative of the quadrupole moment defined 
	 *  by Laarakkers and Poisson, \a Astrophys. \a J. \b 512 , 282 (1999).
	 */
	virtual double mom_quad() const ;	
	

    // Computational routines
    // ----------------------
    public: 
	
	/** Computes metric coefficients from known potentials. 
	 * 
	 *  The calculation is performed starting from the quantities
	 *  \c logn ,  \c dzeta , \c tggg  and \c shift , 
	 *  which are supposed to be up to date.  
	 *  From these,  the following fields are updated:
	 *  \c nnn , \c a_car ,  \c bbb  and \c b_car, as well as 
	 *  the 3-metric \c gamma. 
	 * 
	 */
	void update_metric() ; 
		
	/** Computes \c shift  from \c w_shift  and \c khi_shift 
	 *  according to Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 */
	void fait_shift() ; 
	
	/** Computes \c tnphi  and \c nphi  from the Cartesian 
	 *   components of the shift, stored in \c shift .
	 */
	void fait_nphi() ; 
			
	/** Computes the coefficient \f$\lambda\f$ which ensures that the
	 *	GRV2 virial identity is satisfied.
	 *  \f$\lambda\f$ is the coefficient by which one must multiply
	 *  the quadratic source term \f$\sigma_q\f$ of the 2-D Poisson equation
	 *	\f[
	 *		\Delta_2 u = \sigma_m + \sigma_q
	 *	\f]
	 *  in order that the total source does not contain any monopolar term,
	 *  i.e. in order that
	 *  \f[
	 *		\int_0^{2\pi} \int_0^{+\infty} \sigma(r, \theta)
	 *				\, r \, dr \, d\theta = 0	    \ ,
	 *  \f]
	 *  where \f$\sigma = \sigma_m + \sigma_q\f$.
	 *	\f$\lambda\f$ is computed according to the formula
	 *  \f[
	 *		\lambda = - { \int_0^{2\pi} \int_0^{+\infty} \sigma_m(r, \theta)
	 *				\, r \, dr \, d\theta	    \over
	 * 			\int_0^{2\pi} \int_0^{+\infty} \sigma_q(r, \theta)
	 *				\, r \, dr \, d\theta } \ .
	 *  \f]
	 *  Then, by construction, the new source
	 *	\f$\sigma' = \sigma_m + \lambda \sigma_q\f$ has a vanishing monopolar
	 *  term.
	 *
	 *	@param sou_m [input] matter source term \f$\sigma_m\f$
	 *	@param sou_q [input] quadratic source term \f$\sigma_q\f$
	 *  @return	value of \f$\lambda\f$
	 */
	static double lambda_grv2(const Scalar& sou_m, const Scalar& sou_q) ;
		
};


			//---------------------//
			//   class Kerr_QI     //
			//---------------------//

/**
 * Kerr spacetime in Quasi-Isotropic coordinates (***under development***). 
 * \ingroup(compactobjects)
 *
 * The metric is expressed in Quasi-Isotropic (QI) coordinates :
 * \f[
 *   ds^2 = - N^2 dt^2 + A^2 (dr^2 + r^2 d\theta^2)
 *		       + B^2 r^2 \sin^2\theta (d\varphi - N^\varphi dt)^2
 * \f]
 *
 * 
 */
class Kerr_QI : public Compobj_QI {

    // Data : 
    // -----
    protected:

	/** mass parameter \f$M\f$
	 */
	double mm ;

	/** angular momentum parameter \f$a\f$
	 */
	double aa ; 


    // Derived data : 
    // ------------
    protected:
		 

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star is contructed
	 * @param mass Black hole mass M
	 * @param a_over_m Black hole reduced angular momentum a/M (dimensionless)
     *
	 */
	Kerr_QI(Map& mp_i, double mass, double a_over_m) ;			
	
	
	Kerr_QI(const Kerr_QI& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE*) ). 
	 * 
	 * @param mp_i Mapping on which the star is constructed
	 * @param fich	input file (must have been created by the function
	 *	\c Kerr_QI::sauve )
	 */
	Kerr_QI(Map& mp_i, FILE* fich) ;    		

	virtual ~Kerr_QI() ;			///< Destructor

    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0  all the pointers on derived quantities
	virtual void set_der_0x0() const ; 

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Kerr_QI 
	void operator=(const Kerr_QI& ) ;	
	
    // Accessors
    // ---------
    public:

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
			

    // Computational routines
    // ----------------------
    public: 
	

};


#endif
