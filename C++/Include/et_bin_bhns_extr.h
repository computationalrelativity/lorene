/*
 *  Definition of Lorene class Et_bin_bhns_extr
 *
 */

/*
 *   Copyright (c) 2004-2005 Keisuke Taniguchi
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

#ifndef __ET_BIN_BHNS_EXTR_H_ 
#define __ET_BIN_BHNS_EXTR_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2005/02/28 23:04:56  k_taniguchi
 * Addition of two indicators for the backfround metric and the boundary
 * condition, and some codes for the conformally flat case
 *
 * Revision 1.3  2005/01/31 20:25:22  k_taniguchi
 * Change the argument of equil_bhns_extr_ylm.
 *
 * Revision 1.2  2004/12/29 16:26:02  k_taniguchi
 * Addition of a comutational method to calculate equilibrium figures
 * with a multipole falloff condition at the outer boundary.
 *
 * Revision 1.1  2004/11/30 20:36:31  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "etoile.h"

/**
 * Class for a neutron star in black hole - neutron star binary systems.
 * \ingroup(star)
 *
 * In this class, we assume that the mass ratio of NS to BH is extreme,
 * and the effect from a black hole is treated as the external,
 * background field of the Kerr-Schild metric or the conformally flat metric.
 * They are splitted by the boolian indicator.
 *
 */
class Et_bin_bhns_extr : public Etoile_bin {

    // Data
    // ----
    protected:

        /** Indicator of the background metric:
	 *  \c true  for the Kerr-Shild metric,
	 *  \c false  for the conformally flat one
	 */
        bool kerrschild ;

	/** Indicator of the boundary condition:
	 *  \c true  for the multipole falloff condition,
	 *  \c false  for the \f$1/r\f$ one
	 */
	bool multipole ;

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
	 * @param irrot should be \c true  for an irrotational star, 
	 *		    \c false  for a corotating one
	 * @param ref_triad_i  Reference triad ("absolute frame"), with
	 *	    respect to which the components of all the member 
	 *	    \c Tenseur 's are defined, except for \c w_shift 
	 *	    and \c ssjm1_wshift  whose components are defined
	 *	    with respect to the mapping \c mp  Cartesian triad. 
	 * @param kerrs should be \c true  for the Kerr-Schild background
	 *              metric, \c false  for the conformally flat one
	 * @param multi should be \c true  for the multipol falloff boundary
	 *              condition, \c false  for the \f$1/r\f$ one
	 */
	Et_bin_bhns_extr(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i,
			 bool irrot, const Base_vect& ref_triad_i,
			 bool kerrs, bool multi) ;

	Et_bin_bhns_extr(const Et_bin_bhns_extr& ) ;	///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE*) ). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param ref_triad_i  Reference triad ("absolute frame"), with
	 *	    respect to which the components of all the member 
	 *	    \c Tenseur 's are defined, except for \c w_shift 
	 *	    and \c ssjm1_wshift  whose components are defined
	 *	    with respect to the mapping \c mp  Cartesian triad. 
	 * @param fich	input file (must have been created by the function
	 *	\c sauve )
	 */
	Et_bin_bhns_extr(Map& mp_i, const Eos& eos_i,
			 const Base_vect& ref_triad_i, FILE* fich) ;

	virtual ~Et_bin_bhns_extr() ;			///< Destructor
 

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Et_bin_bhns_extr
	void operator=(const Et_bin_bhns_extr&) ;	

    // Accessors
    // ---------
    public:

	/** Returns \c true  for the Kerr-Schild background metric,
	 *  \c false  for the conformally flat one
	 */
	bool in_kerrschild() const {return kerrschild ;} ;

	/** Returns \c true  for the multipole falloff boundary condition,
	 *  \c false  for the \f$1/r\f$ one
	 */
	bool with_multipole() const {return multipole ;} ;

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
    // Computational routines
    // ----------------------
    public:
	/** Computes the hydrodynamical quantities relative to the Eulerian
	 *  observer from those in the fluid frame, as well as 
	 *  \c wit_w  and \c loggam in the Kerr-Schild background metric
	 *  or in the conformally flat one
	 *
	 *  The calculation is performed starting from the quantities
	 *  \c ent , \c ener , \c press , \c a_car  and \c bsn ,  
	 *  which are supposed to be up to date.  
	 *  From these,  the following fields are updated:
	 *  \c gam_euler , \c u_euler , \c ener_euler , \c s_euler , 
	 *  \c wit_w  and \c loggam . 
	 *
	 *  @param mass mass of the BH
	 *  @param sepa separation between NS and BH
	 *
	 */
	void hydro_euler_extr(const double& mass, const double& sepa) ;

	/** Computes metric coefficients from known potentials,
	 *  when the companion is a black hole with the Kerr-Schild metric
	 *  or with the conformally flat one
	 *
	 *  The calculation is performed starting from the quantities
	 *  \c logn_auto ,  \c beta_auto , \c shift_auto
	 *  which are supposed to be up to date.
	 *  Using the analytical form of the lapse of BH,
	 *  the conformal factor of BH, and the shift vector of BH,
	 *  the following fields are updated:
	 *  \c nnn ,  \c a_car ,  \c shift ,
	 *  \c d_logn_auto , \c d_beta_auto , \c tkij_auto ,
	 *  \c akcar_auto .
	 *
	 *  @param mass mass of the BH
	 *  @param sepa separation between NS and BH
	 *
	 */
	void update_metric_extr(const double& mass, const double& sepa) ;

	/** Computes the derivative of metric functions related to the
	 *  companion black hole with the Kerr-Schild metric
	 *  or with the conformally flat one
	 *
	 *  The calculation is performed starting from the quantities
	 *  \c comp.d_logn_auto ,  \c comp.d_beta_auto ,
	 *  \c comp.tkij_auto 
	 *  which are supposed to be up to date.
	 *  From these,  the following fields are updated:
	 *  \c d_logn_comp , \c d_beta_comp , \c tkij_comp ,
	 *  \c akcar_comp .
	 *
	 *  @param mass mass of the BH
	 *  @param sepa separation between NS and BH
	 *
	 */
	 void update_metric_der_comp_extr(const double& mass,
					  const double& sepa) ;

	 /** Computes the quantities \c bsn  and \c pot_centri
	 *  in the Kerr-Schild background metric
	 *  or in the conformally flat one
	 * 
	 *  The calculation is performed starting from the quantities
	 *  \c nnn , \c shift ,  \c a_car ,  
	 *  which are supposed to be up to date.  
	 * 
	 *  @param omega  angular velocity with respect to an asymptotically 
	 *		  inertial observer
	 *  @param mass mass of the BH
	 *  @param sepa separation between NS and BH
	 *
	 */
	void kinematics_extr(double omega, const double& mass,
			     const double& sepa) ;

	/** Computes \c tkij_auto  and \c akcar_auto  from
	 *  \c shift_auto , \c nnn  and \c a_car .
	 *  in the Kerr-Schild background metric
	 *  or in the conformally flat one
	 *
	 *  @param mass mass of the BH
	 *  @param sepa separation between NS and BH
	 *
	 */
	void extrinsic_curv_extr(const double& mass, const double& sepa) ;

	/** Computes an equilibrium configuration of a BH-NS binary system
	 *  in the Kerr-Schild background metric using the \f$1/r\f$
	 *  falloff boundary condition
	 * 
	 *  The values of \c logn_comp , \c beta_comp , \c pot_centri 
	 *  are held fixed during the iteration. 
	 *  
	 *  @param ent_c  [input] Central enthalpy
	 *  @param mass   [input] Mass of BH
	 *  @param sepa   [input] Orbital separation
	 *  @param mermax [input] Maximum number of steps 
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    Map_et::poisson
	 *  @param relax_poisson [input]  Relaxation factor in Map_et::poisson
	 *  @param mermax_potvit [input]  Maximum number of steps in 
	 *				  Map_radial::poisson_compact
	 *  @param relax_potvit [input]   Relaxation factor in 
	 *				  Map_radial::poisson_compact
	 *  @param np_filter [input]  Number of coefficients in phi which are
	 *                            deleted by filter
	 *  @param thres_adapt  [input]   Threshold on dH/dr for the adaptation 
	 *				  of the mapping
	 *  @param diff [output]   1-D \c Tbl  for the storage of some
	 *			    error indicators : 
	 *	    \li \c diff(0)  : Relative change in the enthalpy field
	 *			      between two successive steps 
	 *	    \li \c diff(1)  : Relative error returned by the routine
	 *				\c Etoile_bin::velocity_potential   
	 *	    \li \c diff(2)  : Relative error in the resolution of the
	 *			    Poisson equation for \c logn_auto    
	 *	    \li \c diff(3)  : Relative error in the resolution of the
	 *			    Poisson equation for \c beta_auto    
	 *	    \li \c diff(4)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (x comp.)   
	 *	    \li \c diff(5)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (y comp.)   
	 *	    \li \c diff(6)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (z comp.)   
	 */
	void equil_bhns_extr_ks(double ent_c, const double& mass,
				const double& sepa, int mermax,
				int mermax_poisson, 
				double relax_poisson, int mermax_potvit, 
				double relax_potvit, int np_filter,
				double thres_adapt, Tbl& diff) ;

	/** Computes an equilibrium configuration of a BH-NS binary system
	 *  in the conformally flat background metric using the \f$1/r\f$
	 *  falloff boundary condition
	 * 
	 *  The values of \c logn_comp , \c beta_comp , \c pot_centri 
	 *  are held fixed during the iteration. 
	 *  
	 *  @param ent_c  [input] Central enthalpy
	 *  @param mass   [input] Mass of BH
	 *  @param sepa   [input] Orbital separation
	 *  @param mermax [input] Maximum number of steps 
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    Map_et::poisson
	 *  @param relax_poisson [input]  Relaxation factor in Map_et::poisson
	 *  @param mermax_potvit [input]  Maximum number of steps in 
	 *				  Map_radial::poisson_compact
	 *  @param relax_potvit [input]   Relaxation factor in 
	 *				  Map_radial::poisson_compact
	 *  @param np_filter [input]  Number of coefficients in phi which are
	 *                            deleted by filter
	 *  @param thres_adapt  [input]   Threshold on dH/dr for the adaptation 
	 *				  of the mapping
	 *  @param diff [output]   1-D \c Tbl  for the storage of some
	 *			    error indicators : 
	 *	    \li \c diff(0)  : Relative change in the enthalpy field
	 *			      between two successive steps 
	 *	    \li \c diff(1)  : Relative error returned by the routine
	 *				\c Etoile_bin::velocity_potential   
	 *	    \li \c diff(2)  : Relative error in the resolution of the
	 *			    Poisson equation for \c logn_auto    
	 *	    \li \c diff(3)  : Relative error in the resolution of the
	 *			    Poisson equation for \c beta_auto    
	 *	    \li \c diff(4)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (x comp.)   
	 *	    \li \c diff(5)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (y comp.)   
	 *	    \li \c diff(6)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (z comp.)   
	 */
	void equil_bhns_extr_cf(double ent_c, const double& mass,
				const double& sepa, int mermax,
				int mermax_poisson, 
				double relax_poisson, int mermax_potvit, 
				double relax_potvit, int np_filter,
				double thres_adapt, Tbl& diff) ;

	/** Computes an equilibrium configuration of a BH-NS binary system
	 *  in the Kerr-Schild background metric using the multipole falloff
	 *  boundary condition
	 * 
	 *  The values of \c logn_comp , \c beta_comp , \c pot_centri 
	 *  are held fixed during the iteration. 
	 *  
	 *  @param ent_c  [input] Central enthalpy
	 *  @param mass   [input] Mass of BH
	 *  @param sepa   [input] Orbital separation
	 *  @param nu_int [input] Multipole moment for logn in the previous
	 *                        step
	 *  @param beta_int [input] Multipole moment for beta in the
	 *                          previous step
	 *  @param shift_int [input] Multipole moment for shift in the
	 *                           previous step
	 *  @param mermax [input] Maximum number of steps 
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    Map_et::poisson
	 *  @param relax_poisson [input]  Relaxation factor in Map_et::poisson
	 *  @param relax_ylm [input] Relaxation factor on the outer boundary
	                             condition
	 *  @param mermax_potvit [input]  Maximum number of steps in 
	 *				  Map_radial::poisson_compact
	 *  @param relax_potvit [input]   Relaxation factor in 
	 *				  Map_radial::poisson_compact
	 *  @param np_filter [input]  Number of coefficients in phi which are
	 *                            deleted by filter
	 *  @param thres_adapt  [input]   Threshold on dH/dr for the adaptation 
	 *				  of the mapping
	 *  @param diff [output]   1-D \c Tbl  for the storage of some
	 *			    error indicators : 
	 *	    \li \c diff(0)  : Relative change in the enthalpy field
	 *			      between two successive steps 
	 *	    \li \c diff(1)  : Relative error returned by the routine
	 *				\c Etoile_bin::velocity_potential   
	 *	    \li \c diff(2)  : Relative error in the resolution of the
	 *			    Poisson equation for \c logn_auto    
	 *	    \li \c diff(3)  : Relative error in the resolution of the
	 *			    Poisson equation for \c beta_auto    
	 *	    \li \c diff(4)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (x comp.)   
	 *	    \li \c diff(5)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (y comp.)   
	 *	    \li \c diff(6)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (z comp.)   
	 */
	void equil_bhns_extr_ylm_ks(double ent_c, const double& mass,
				    const double& sepa, double* nu_int,
				    double* beta_int, double* shift_int,
				    int mermax, int mermax_poisson, 
				    double relax_poisson, double relax_ylm,
				    int mermax_potvit, double relax_potvit,
				    int np_filter,
				    double thres_adapt, Tbl& diff) ;

	/** Computes an equilibrium configuration of a BH-NS binary system
	 *  in the conformally flat background metric using the multipole
	 *  falloff boundary condition
	 * 
	 *  The values of \c logn_comp , \c beta_comp , \c pot_centri 
	 *  are held fixed during the iteration. 
	 *  
	 *  @param ent_c  [input] Central enthalpy
	 *  @param mass   [input] Mass of BH
	 *  @param sepa   [input] Orbital separation
	 *  @param nu_int [input] Multipole moment for logn in the previous
	 *                        step
	 *  @param beta_int [input] Multipole moment for beta in the
	 *                          previous step
	 *  @param shift_int [input] Multipole moment for shift in the
	 *                           previous step
	 *  @param mermax [input] Maximum number of steps 
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    Map_et::poisson
	 *  @param relax_poisson [input]  Relaxation factor in Map_et::poisson
	 *  @param relax_ylm [input] Relaxation factor on the outer boundary
	                             condition
	 *  @param mermax_potvit [input]  Maximum number of steps in 
	 *				  Map_radial::poisson_compact
	 *  @param relax_potvit [input]   Relaxation factor in 
	 *				  Map_radial::poisson_compact
	 *  @param np_filter [input]  Number of coefficients in phi which are
	 *                            deleted by filter
	 *  @param thres_adapt  [input]   Threshold on dH/dr for the adaptation 
	 *				  of the mapping
	 *  @param diff [output]   1-D \c Tbl  for the storage of some
	 *			    error indicators : 
	 *	    \li \c diff(0)  : Relative change in the enthalpy field
	 *			      between two successive steps 
	 *	    \li \c diff(1)  : Relative error returned by the routine
	 *				\c Etoile_bin::velocity_potential   
	 *	    \li \c diff(2)  : Relative error in the resolution of the
	 *			    Poisson equation for \c logn_auto    
	 *	    \li \c diff(3)  : Relative error in the resolution of the
	 *			    Poisson equation for \c beta_auto    
	 *	    \li \c diff(4)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (x comp.)   
	 *	    \li \c diff(5)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (y comp.)   
	 *	    \li \c diff(6)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (z comp.)   
	 */
	void equil_bhns_extr_ylm_cf(double ent_c, const double& mass,
				    const double& sepa, double* nu_int,
				    double* beta_int, double* shift_int,
				    int mermax, int mermax_poisson, 
				    double relax_poisson, double relax_ylm,
				    int mermax_potvit, double relax_potvit,
				    int np_filter,
				    double thres_adapt, Tbl& diff) ;

	/** Tests the resolution of the Poisson equations
	 *  when the NS has no matter source.
	 *  The solution should be the same as the Kerr-Schild metric
	 * 
	 *  @param mass   [input] Mass of BH
	 *  @param sepa   [input] Orbital separation
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    Map_et::poisson
	 *  @param relax_poisson [input]  Relaxation factor in Map_et::poisson
	 *  @param mermax_potvit [input]  Maximum number of steps in 
	 *				  Map_radial::poisson_compact
	 *  @param relax_potvit [input]   Relaxation factor in 
	 *				  Map_radial::poisson_compact
	 *  @param diff [output]   1-D \c Tbl  for the storage of some
	 *			    error indicators : 
	 *	    \li \c diff(0)  : Relative error returned by the routine
	 *				\c Etoile_bin::velocity_potential   
	 *	    \li \c diff(1)  : Relative error in the resolution of the
	 *			    Poisson equation for \c logn_auto    
	 *	    \li \c diff(2)  : Relative error in the resolution of the
	 *			    Poisson equation for \c beta_auto    
	 *	    \li \c diff(3)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (x comp.)   
	 *	    \li \c diff(4)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (y comp.)   
	 *	    \li \c diff(5)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (z comp.)   
	 */
	void test_bhns_extr(const double& mass,
			    const double& sepa, int mermax_poisson,
			    double relax_poisson, int mermax_potvit,
			    double relax_potvit, Tbl& diff) ;

	/** Computes the non-translational part of the velocity scalar
	 *  potential \f$\psi0\f$ by solving the continuity equation
	 *  in the Kerr-Schild background metric
	 *  or in the conformally flat one
	 *
	 *  @param mass    [input] Mass of BH
	 *  @param sepa    [input] Orbital separation
	 *  @param mermax  [input] Maximum number of steps in the iteration
	 *  @param precis  [input] Required precision: the iteration will
	 *			   be stopped when the relative difference
	 *			   on \f$\psi0\f$ between two successive steps
	 *			   is lower than \c precis .
	 *  @param relax   [input] Relaxation factor.  
	 *
	 *  @return Relative error of the resolution obtained by comparing
	 *	    the operator acting on the solution with the source.
	 */
	double velocity_pot_extr(const double& mass, const double& sepa,
				 int mermax, double precis, double relax) ;

	/** Searches the position of the maximum enthalpy
	 *  @param xx [output] x-coordinate of the maximum enthalpy
	 *  @param yy [output] y-coordinate of the maximum enthalpy
	 */
	void ent_max_search(double& xx, double& yy) const ;

	/** Searches the position (phi) of the longest radius of NS
	 *  from the position of the maximum enthalpy
	 *  @param x_max [input] x-coordinate of the maximum enthalpy
	 *  @param y_max [input] y-coordinate of the maximum enthalpy
	 */
	double phi_longest_rad(double x_max, double y_max) const ;

	/// Constructs spherical harmonics
	void get_ylm(int nylm, Cmp** ylmvec) const ;

	/// Computes multipole moments
	void get_integrals(int nylm, double* intvec, Cmp** ylmvec,
			   Cmp source) const ;

	friend class Bin_bhns_extr ;

};

#endif
