/*
 *  Definition of Lorene class Star_rot_Dirac
 *
 */

/*
 *   Copyright (c) 2005 Lap-Ming Lin & Jerome Novak
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

#ifndef __STAR_ROT_DIRAC_H_ 
#define __STAR_ROT_DIRAC_H_ 

/*
 *
 * $Header$
 *
 */


// Headers Lorene
#include "star.h"


/**
 * Class for relativistic rotating stars in Dirac gauge and maximal slicing. 
 * (*** Under development ***) \ingroup (star)
 * 
 */
class Star_rot_Dirac : public Star {

    // Data : 
    // -----
    protected:

       double omega ;  ///< Rotation angular velocity (\c [f_unit] )

       // Quantities related to the conformal factor and lapse
       //----------------------------------

       Scalar psi4 ;   ///< Conformal factor \f$\Psi^4\f$
       Scalar psi2 ;   ///< \f$\Psi^2\f$
       Scalar qqq ;    ///< \f$Q = \Psi^2 N\f$ 
       Scalar ln_psi ; ///< \f$ln(\Psi)\f$


       // Fluid quantities
       //-----------------------------------------

       /**
	* Momentum density 3-vector with respect to the 
	* Eulerian observer
	*/
       Vector j_euler ; 
       Scalar v2 ; ///< \f$\gamma_{ij}v^i v^j\f$

       // Metric stuff 
       //-------------------------------------

       /// flat metric \f$f_{ij}\f$ (spherical components)
       const Metric_flat& flat ; 

       Metric tgamma  ;  ///< \f$\tilde{\gamma}_{ij}\f$
       Sym_tensor aa ; ///< \f$A^{ij}\f$ 
       Sym_tensor taa ; ///< \f$\tilde{A}_{ij}\f$
       Scalar aa_quad ; ///< \f$\tilde{A}_{ij} A^{ij}\f$ 

       /**
	* \f$h^{ij}\f$ is defined by \f$\tilde{\gamma}^{ij}=f^{ij}+h^{ij}\f$.
	* We impose the Dirac gauge \f$D_j h^{ij}\f$ explicitly by 
	* defining \f$h^{ij}\f$ to be a symmetric transverse tensor. 
	*/
       Sym_tensor_trans hh ; 


    // Derived data : 
    // ------------
    protected:

       // More to come later.....
       //----------------------------

       mutable double* p_angu_mom ; ///< Angular momentum 
       mutable double* p_grv2 ; ///< Error on the virial identity GRV2
       mutable double* p_grv3 ; ///< Error on the virial identity GRV3


    // Constructors - Destructor
    // -------------------------
    public:

	/**
	 *
	 * Standard constructor.
	 *
	 * @param mp_i Mapping on which the star will be defined
	 * @param nzet_i Number of domains occupied by the star
         * @param eos_i Equation of state of the stellar matter
	 *
	 */
       Star_rot_Dirac(Map& mp_i, int nzet_i, const Eos& eos_i) ;  

       Star_rot_Dirac(const Star_rot_Dirac& ) ; ///< Copy constructor


        /** Constructor from a file (see \c sauve(FILE*) ).
         *
         * @param mp_i Mapping on which the star will be defined
         * @param eos_i Equation of state of the stellar matter
         * @param fich  input file (must have been created by the function
         *      \c sauve )
         */
       Star_rot_Dirac(Map& mp_i, const Eos& eos_i, FILE* fich) ;

	
       virtual ~Star_rot_Dirac() ;	 ///< Destructor
 

    // Memory management
    // -----------------
     protected:

	/// Deletes all the derived quantities
       virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
       void set_der_0x0() const ; 

	/** Sets to \c ETATNONDEF  (undefined state) the hydrodynamical
         *  quantities relative to the Eulerian observer.
         */
       virtual void del_hydro_euler() ;


    // Mutators / assignment
    // ---------------------
    public:
	
	/// Assignment to another \c Star_rot_Dirac
       void operator=(const Star_rot_Dirac& ) ;	
	
    // Accessors
    // ---------
    public:

       /** 
	* Returns the rotation angular velocity \f$\Omega\f$
	*/
       double get_omega() const {return omega;} ;


       /**
	* Returns the conformal factor \f$\Psi^4\f$
	*/
       const Scalar& get_psi4() const {return psi4;} ;

       /**
	* Returns \f$\Psi^2\f$
	*/
       const Scalar& get_psi2() const {return psi2;} ;

       /**
	* Returns \f$Q=\Psi^2 N\f$
	*/
       const Scalar& get_qqq() const {return qqq;} ;

       /**
	* Returns \f$ln(\Psi)\f$
	*/
       const Scalar& get_ln_psi() const {return ln_psi;} ;


       // Fluid stuff
       //------------------

       /**
	* Returns the momentum density 3-vector with respect to the 
	* Eulerian observer
	*/
       const Vector& get_j_euler() const {return j_euler;} ;

       /**
	* Reutrns \f$\gamma_{ij}v^i v^j\f$
	*/
       const Scalar& get_v2() const {return v2;} ;


       //Metric stuff
       //-------------------
       /**
	* Returns the conformal metric \f$\tilde{\gamma}_{ij}\f$
	*/
       const Metric get_tgamma() const {return tgamma;} ;

       // documentation comes later......

       /**
	* Returns \f$A^{ij}\f$
	*/
       const Sym_tensor get_aa() const {return aa;} ;

       /**
	* Returns \f$\tilde{A}_{ij}\f$
	*/
       const Sym_tensor get_taa() const {return taa;} ;
       
       /** 
	* Returns \f$\tilde{A}_{ij} A^{ij}\f$
	*/
       const Scalar get_aa_quad() const {return aa_quad;} ;
       
       /** 
	* Returns \f$h^{ij}\f$
	*/
       const Sym_tensor_trans get_hh() const {return hh;} ;




    // Outputs
    // -------
    public:

       virtual void sauve(FILE* ) const ;	    ///< Save in a file
	
    protected:

       virtual ostream& operator>>(ostream& ) const ;

   
    // Global quantities
    //-------------------------
    public:

       virtual double mass_b() const ; ///< Baryonic mass 
       virtual double mass_g() const ; ///< Gravitational mass
       virtual double angu_mom() const ; ///< Angular momentum
       virtual double grv2() const ;  ///< Error on the virial identity GRV2
       virtual double grv3() const ; ///< Error on the virial identity GRV3

       //    More to come later.....



    // Computational routines
    //--------------------------
    public:

       /**
	* Computes the hydrodynamical quantities relative to the Eulerian 
	* observer from those in the fluid frame.
	* 
	* More later......
	*/ 
       virtual void hydro_euler() ; 

       
       /**
	* Computes metric quantities from known potentials.
	* 
	* The calculation is performed starting from \c qqq, \c logn,
	* \c shift, \c hh, which are supposed to be up to date.
	* From these, the following fields are updated: \c nnn, 
	* \c psi4, \c psi2, \c ln_psi, \c tgamma, \c aa, \c taa, and 
	* \c aa_quad. 
	*/
       void update_metric() ;


       /**
	* Computes an equilibrium configuration 
	*/ 
       void equilibrium(double ent_c, double omega0, double fact_omega, 
			int nzadapt, const Tbl& ent_limit,
			const Itbl& icontrol, const Tbl& control,
			double mbar_wanted, double aexp_mass, 
			Tbl& diff)  ;

       /**
	* Solution of the two scalar Poisson equations for rotating 
	* stars in Dirac gauge
	*/
       void solve_logn_f(Scalar& ln_f_new) const ;

       /**
	* Solution of the two scalar Poisson equations for rotating 
	* stars in Dirac gauge
	*/
       void solve_logn_q(Scalar& ln_q_new) const ;

       /**
	* Solution of the two scalar Poisson equations for rotating 
	* stars in Dirac gauge
	*/
       void solve_qqq(Scalar& q_new) const ;

       /**
	* Solution of the shift equation for rotating stars 
	* in Dirac gauge
	*/
       void solve_shift(Vector& shift_new) const ;

       /**
	* Solution of the tensor Poisson equation for rotating stars 
	* in Dirac gauge
	*/
       void solve_hij(Sym_tensor_trans& hij_new) const ;
              
       /**
	* Computes the coefficient \f$\lambda\f$ which ensures that 
	* the GRV2 virial identity is satisfied. 
	* See Etoile_rot::lambda_grv2 for detail.
	*/
       static double lambda_grv2(const Scalar& sou_m, const Scalar& sou_q) ;

};

#endif
