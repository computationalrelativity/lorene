/*
 *  Definition of Lorene class Star_rot_dirac_diff
 *
 */

/*
 *   Copyright (c) 2005 Motoyuki Saijo
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


#ifndef __STAR_ROT_DIRAC_DIFF_H_ 
#define __STAR_ROT_DIRAC_DIFF_H_ 

/*
 *
 * $Header$
 *
 */

// Headers Lorene
#include "star_rot_dirac.h"

namespace Lorene {

/**
 * Class for relativistic differentially rotating stars in Dirac gauge and 
 * maximal slicing.  (*** Under development ***) \ingroup (star)
 * 
 *
 */

class Star_rot_Dirac_diff : public Star_rot_Dirac {

    // Data : 
    // -----
    protected:
	/** Function \f$F(\Omega)\f$ defining the rotation profile.
	 *  This function is linked to the components of the fluid 4-velocity
	 *  by 
	 *  \f[
	 *	F(\Omega) = u^t u_\varphi \ . 
	 *  \f]
	 *  The first argument of \c frot  must be \f$\Omega\f$; the second
	 *  argument contains the parameters, the first of which being
	 *  necessarily the central value of \f$\Omega\f$.   
	 */
	double (*frot)(double, const Tbl&) ;
	
	/** Primitive of the function \f$F(\Omega)\f$, which vanishes at the
	 *  stellar center. 
	 *  The first argument of \c primfrot  must be \f$\Omega\f$; the second
	 *  argument contains the parameters, the first of which being
	 *  necessarily the central value of \f$\Omega\f$.   
	 */
	double (*primfrot)(double, const Tbl&) ;
	
	/** Parameters of the function \f$F(\Omega)\f$.
	 * 
	 * To be used as the second argument of functions \c frot 
	 * and \c primfrot . 
	 * The parameter \c par_frot(0)  must always be the central angular
	 * velocity. 
	 */
	Tbl par_frot ; 
	
	/// Field \f$\Omega(r,\theta)\f$
	Scalar omega_field ; 
	
	double omega_min ;   ///< Minimum value of \f$\Omega\f$
	double omega_max ;   ///< Maximum value of \f$\Omega\f$
	
	/// Field \f$\int_{\Omega_{\rm c}}^\Omega F(\Omega') \, d\Omega' \f$
	Scalar prim_field ; 
	
    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param nzet_i Number of domains occupied by the star
	 * @param relat should be \c true  for a relativistic
	 *			star,  \c false  for a Newtonian one
	 * @param eos_i Equation of state of the stellar matter
	 * @param frot_i Function \f$F(\Omega)\f$ defining the rotation profile.
	 * @param primfrot_i Primitive of \f$F(\Omega)\f$ which vanishes at the
	 *   stellar center
	 * @param par_frot_i Parameters of functions \c frot_i 
	 *     and \c primfrot_i , 
	 *	    \c par_frot_i(0)  being the central value of
	 *	    \f$\Omega\f$. 
	 */
	Star_rot_Dirac_diff(Map& mp_i, int nzet_i, const Eos& eos_i, 
		    double (*frot_i)(double, const Tbl&), 
		    double (*primfrot_i)(double, const Tbl&), 
		    const Tbl& par_frot_i) ;			

	Star_rot_Dirac_diff(const Star_rot_Dirac_diff& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE*) ). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param fich	input file (must have been created by the function
	 *	\c sauve )
	 * @param frot_i Function \f$F(\Omega)\f$ defining the rotation profile.
	 * @param primfrot_i Primitive of \f$F(\Omega)\f$ which vanishes at the
	 *   stellar center
	 */
	Star_rot_Dirac_diff(Map& mp_i, const Eos& eos_i, FILE* fich,
		    double (*frot_i)(double, const Tbl&), 
		    double (*primfrot_i)(double, const Tbl&) ) ;			
	
	virtual ~Star_rot_Dirac_diff() ;			///< Destructor


    // Memory management
    // -----------------

    // Everything is inherited from Star_rot_dirac

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Star_rot_Dirac_diff 
	void operator=(const Star_rot_Dirac_diff& ) ;	
	
    // Accessors
    // ---------
    public:
	/// Returns the angular velocity field \f$\Omega\f$
	const Scalar& get_omega_field() const {return omega_field;} ; 

	/** Returns the central value of the rotation angular velocity 
	 *  (\c [f_unit] )
	 */ 
       virtual double get_omega_c() const ;


    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
/* 	/// Display in polytropic units */
/* 	virtual void display_poly(ostream& ) const ; */

	/// Display
	friend ostream& operator<<(ostream& , const Star& ) ;	


    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    


    // Computational routines
    // ----------------------
    public: 

	virtual double tsw() const ;		///< Ratio T/W

	/** Computes the hydrodynamical quantities relative to the Eulerian
	 *  observer from those in the fluid frame.
	 *
	 *  The calculation is performed starting from the quantities
	 *  \c omega_field , \c ent , \c ener , \c press , 
	 *  and \c a_car ,  
	 *  which are supposed to be up to date.  
	 *  From these,  the following fields are updated:
	 *  \c gam_euler , \c u_euler , \c ener_euler , \c s_euler . 
	 * 
	 */
	virtual void hydro_euler() ; 

	/** Computes \f$\Omega(r,\theta)\f$ (member \c omega_field ).
	 *  
	 *  The computation amounts to solving the equation
	 *  \f[
	 *	F(\Omega) - {\gamma_{\varphi \varphi} (\beta^\varphi + \Omega)
	 *	    \over N^2 - \gamma_{\varphi \varphi}(\beta^\varphi + \Omega)^2} = 0
	 *  \f]
	 *  for \f$\Omega\f$. 
	 *
	 *  @param omeg_min [input] Lower bound of the interval for
	 *			     searching omega
	 *  @param omeg_max [input] Higher bound of the interval for
	 *			     searching omega
	 *  @param precis [input] Required precision in the determination of
	 *			the zero by the secant method
	 *  @param nitermax [input] Maximum number of iterations in the secant 
	 *			    method to compute the zero. 
	 */
	void fait_omega_field(double omeg_min, double omeg_max,
			      double precis, int nitermax) ;
	
	/// Computes the member \c prim_field  from \c omega_field . 
	void fait_prim_field() ;
	
	/** Evaluates \f$F(\Omega)\f$, where \e F is the function
	 *  defining the rotation profile.
	 *  This function is linked to the components of the fluid 4-velocity
	 *  by
	 *  \f[
	 *	F(\Omega) = u^t u_\varphi \ .
	 *  \f]
	 *
	 *  \c funct_omega  calls \c frot  with the parameters
	 *  \c par_frot .
	 *
	 *  @param omeg [input] value of \f$\Omega\f$
	 *  @return value of \f$F(\Omega)\f$
	 *
	 */
	double funct_omega(double omeg) const ;  	

	/** Evaluates the primitive of \f$F(\Omega)\f$, where \e F is the function
	 *  defining the rotation profile.
	 *
	 *  @param omeg [input] value of \f$\Omega\f$
	 *  @return value of 
	 *	\f$\int_{\Omega_{\rm c}}^\Omega F(\Omega') \, d\Omega' \f$
	 *
	 */
	double prim_funct_omega(double omeg) const ;  	

	/** Computes an equilibrium configuration.
	 *  
	 *  @param ent_c  [input] Central enthalpy 
	 *  @param omega0  [input] Requested central angular velocity 
	 *			     (if \c fact_omega=1. )
	 *  @param fact_omega [input] 1.01 = search for the Keplerian frequency,
	 *			      1. = otherwise.
	 *  @param nzadapt  [input] Number of (inner) domains where the mapping 
	 *			    adaptation to an iso-enthalpy surface
	 *			    should be performed
	 *  @param ent_limit [input] 1-D \c Tbl of dimension \c nzet  which
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
	 * 	 
	 *  @param control [input] Set of parameters (stored as a 
	 *			    1-D \c Tbl  of size 7) to control the 
	 *			    iteration: 
	 *	\li \c control(0) = precis  : threshold on the enthalpy relative 
	 *				change for ending the computation 
	 *	\li \c control(1) = omega_ini  : initial central angular velocity, 
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
	 *
	 *  @param mbar_wanted [input] Requested baryon mass (effective only 
	 *				if \c mer_mass > mer_max )
	 *  @param aexp_mass [input] Exponent for the increase factor of the 
	 *			      central enthalpy to converge to the 
	 *			      requested baryon mass
	 *  @param diff [output]   1-D \c Tbl  of size 7 for the storage of 
	 *			    some error indicators : 
	 *	    \li \c diff(0)  : Relative change in the enthalpy field
	 *			      between two successive steps 
	 *	    \li \c diff(1)  : Relative error in the resolution of the
	 *			    Poisson equation for \c nuf    
	 *	    \li \c diff(2)  : Relative error in the resolution of the
	 *			    Poisson equation for \c nuq    
	 *	    \li \c diff(3)  : Relative error in the resolution of the
	 *			    Poisson equation for \c dzeta    
	 *	    \li \c diff(4)  : Relative error in the resolution of the
	 *			    Poisson equation for \c tggg    
	 *	    \li \c diff(5)  : Relative error in the resolution of the
	 *			    equation for \c shift  (x comp.)   
	 *	    \li \c diff(6)  : Relative error in the resolution of the
	 *			    equation for \c shift  (y comp.)   
	 */
	virtual void equilibrium(double ent_c, double omega0, double fact_omega, 
			 int nzadapt, const Tbl& ent_limit,
			 const Itbl& icontrol, const Tbl& control,
			 double mbar_wanted, double aexp_mass, 
			 Tbl& diff) ;
	
 };

}
#endif
