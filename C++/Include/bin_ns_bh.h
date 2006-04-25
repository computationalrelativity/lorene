/*
 *  Definition of Lorene class Bin_ns_bh
 *
 */

/*
 *   Copyright (c) 2002  Philippe Grandclement, Keisuke Taniguchi,
 *              Eric Gourgoulhon
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

#ifndef __BIN_NS_BH_H_
#define __BIN_NS_BH_H_

/*
 * $Id$
 * $Log$
 * Revision 1.13  2006/04/25 07:21:54  p_grandclement
 * Various changes for the NS_BH project
 *
 * Revision 1.12  2005/12/01 12:59:08  p_grandclement
 * Files for bin_ns_bh project
 *
 * Revision 1.11  2005/11/30 11:09:03  p_grandclement
 * Changes for the Bin_ns_bh project
 *
 * Revision 1.10  2005/10/18 13:12:31  p_grandclement
 * update of the mixted binary codes
 *
 * Revision 1.9  2005/08/29 15:10:12  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.8  2004/06/09 06:23:50  k_taniguchi
 * Introduce analytical_omega() and analytical_shift().
 *
 * Revision 1.7  2003/11/25 07:25:44  k_taniguchi
 * Change the attribute of fait_decouple() from private to public.
 *
 * Revision 1.6  2003/10/24 16:56:30  k_taniguchi
 * Add the method for the calculation of the orbital angular velocity
 *
 * Revision 1.5  2003/10/24 12:45:22  k_taniguchi
 * Change the class for the star from Etoile_bin to Et_bin_nsbh
 *
 * Revision 1.4  2003/02/13 16:40:24  p_grandclement
 * Addition of various things for the Bin_ns_bh project, non of them being
 * completely tested
 *
 * Revision 1.3  2002/12/19 14:46:16  e_gourgoulhon
 *
 * Modified prototype of functions set_omega and set_x_axe
 *
 * Revision 1.2  2002/12/18 10:29:18  e_gourgoulhon
 *
 * Added set_omega() and set_x_axe()
 *
 * Revision 1.1  2002/12/17 13:10:49  e_gourgoulhon
 * Definition of class Bin_ns_bh
 *
 *
 *
 *
 * $Header$
 *
 */


// Lorene headers
#include "et_bin_nsbh.h"
#include "bhole.h"

/**
 * Neutron star - black hole binary system.
 *
 * The class {\tt Bin\_ns\_bh} is composed of an object of class
 * {\tt Et\_bin\_nsbh} and an object of class {\tt Bhole}.
 *
 * @version #$Id$#
 */
class Bin_ns_bh {

    // Data :
    // -----
    private:
	/** Cartesian triad of the absolute reference frame
	 *
	 */
	const Base_vect_cart ref_triad ;

	/// The neutron star
	Et_bin_nsbh star ;

	/// The black hole
	Bhole hole ;

	/** Angular velocity with respect to an asymptotically inertial
	 *  observer
	 */
	double omega ;

	/** Absolute X coordinate of the rotation axis
	 */
	double x_axe ;

    // Derived data :
    // ------------
    private:
	/// Total ADM mass of the system
	mutable double* p_mass_adm ;

	/// Total Komar mass of the system
	mutable double* p_mass_kom ;

	/// Total angular momentum of the system
	mutable Tbl* p_angu_mom ;

	/// Total energy of the system
	mutable double* p_total_ener ;

	/// Virial theorem error
	mutable double* p_virial ;

	/// Virial theorem error by E.Gourgoulhon and S.Bonazzola.
	mutable double* p_virial_gb ;

	/// Virial theorem error by J.L.Friedman, K.Uryu, and M.Shibata.
	mutable double* p_virial_fus ;

	/// Relative error on the Hamiltonian constraint
	mutable double* p_ham_constr ;

	/// Relative error on the momentum constraint
	mutable Tbl* p_mom_constr ;


    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor.
	 *
	 * @param mp_ns Mapping on which {\tt star} will be defined
	 * @param nzet Number of domains occupied by {\tt star}
	 * @param eos Equation of state of {\tt star}
	 * @param irrot_ns should be {\tt true} if {\tt star} is irrotational,
	 *		    {\tt false} if {\tt star} is corotating
	 * @param mp_bh Mapping on which {\tt bhole} will be defined
	 *
	 */
	Bin_ns_bh(Map& mp_ns, int nzet, const Eos& eos, bool irrot_ns,
	        Map_af& mp_bh) ;


	Bin_ns_bh(const Bin_ns_bh& ) ;		/// Copy constructor

	/** Constructor from a file (see {\tt sauve(FILE* )}).
	 *
	 * @param mp_ns Mapping on which {\tt star} will be defined
	 * @param eos Equation of state of {\tt star}
	 * @param mp_bh Mapping on which {\tt star} will be defined
	 * @param fich	input file (must have been created by the function
	 *	{\tt sauve})
	 */
	Bin_ns_bh(Map& mp_ns, const Eos& eos, Map_af& mp_bh, FILE* fich) ;

	virtual ~Bin_ns_bh() ;			/// Destructor


    // Memory management
    // -----------------
    private:
	/// Deletes all the derived quantities
	void del_deriv() const ;

	/// Sets to {\tt 0x0} all the pointers on derived quantities
	void set_der_0x0() const ;


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Bin_ns_bh
	void operator=(const Bin_ns_bh& ) ;

	/// Read/write of the neutron star
	Et_bin_nsbh& set_ns()
	    { del_deriv() ;
               return star ;} ;

	/// Read/write of the black hole
	Bhole& set_bh()
	    { del_deriv() ;
	       return hole ;} ;

	/// Sets the orbital angular velocity [{\tt f\_unit}]
	void set_omega(double ) ;

	/// Sets the absolute coordinate X of the rotation axis [{\tt r\_unit}]
	void set_x_axe(double ) ;

    // Accessors
    // ---------
    public:
	/// Returns a constant reference to the neutron star
	const Et_bin_nsbh& get_ns() const
	    { return star ;} ;

	/// Returns a constant reference to the black hole
	const Bhole& get_bh() const
	    { return hole ;} ;

	/// Returns the orbital velocity
	double get_omega() const
	    { return omega ;} ;

	/// Returns a constant reference to the black hole
	double get_x_axe() const
	    { return x_axe ;} ;

	/// Return the separation
	double separation() const  ;

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    /// Save in a file

	/// Display
	friend ostream& operator<<(ostream& , const Bin_ns_bh& ) ;

    private:
	/// Operator >> (function called by the operator <<).
	ostream& operator>>(ostream& ) const ;

    public:
	/** Function used to compute the {\tt decouple} functions for both
	 * the NS and the BH
	 **/
	void fait_decouple() ;

	/** Computation of the extrinsic curvature tensor for both
	 * {\tt star} and {\tt bhole}.
	 **/
	void fait_tkij() ;

    // Computational routines
    // ----------------------
    public:
	/** Computes the orbital angular velocity {\tt omega}
	 *
	 *  @param fact_omeg_min [input] : determines the lower bound of the 
	 *		interval {\tt [omega\_min, omega\_max]} in which 
	 *		{\tt omega} is searched by 
	 *		{\tt omega\_min = fact\_omeg\_min * omega}, 
	 *		where {\tt omega} is the previous value of the 
	 *		angular velocity 
	 *		(typical value : {\tt fact\_omeg\_min = 0.5})
	 *
	 *  @param fact_omeg_max [input] : determines the higher bound of the 
	 *		interval {\tt [omega\_min, omega\_max]} in which 
	 *		{\tt omega} is searched by 
	 *		{\tt omega\_max = fact\_omeg\_max * omega}, 
	 *		where {\tt omega} is the previous value of the 
	 *		angular velocity.
	 *		(typical value : {\tt fact\_omeg\_max = 1.5})
	 *
	 */
	void orbit_omega(double fact_omeg_min, double fact_omeg_max) ;

	/** Sets the orbital angular velocity of the neutron star
	 *   to some 2-PN analytical value
	 */
	void analytical_omega() ;

	/** Sets some analytical template for the shift vector (via the
	 *   members {\tt w\_shift} and {\tt khi\_shift} of the neutron star
	 */
	void analytical_shift() ;

	void init_auto () ;
	void affecte (const Bin_ns_bh&) ;
	void pseudo_misner (int&, int, double, double) ;
	double adm_systeme() const ;
	double adm_systeme_volume() const ;
	double komar_systeme() const ;
	double moment_systeme_inf() const ;
	double moment_systeme_hor() const ; 
	double smarr() const ;
	Tbl linear_momentum_systeme_inf() const ;
	double viriel() const ;
	void coal (double, double, int, int,  double, double, double, double, double, int, const int) ;
	double distance_propre_axe_bh (const int nr  = 65) const ;
	double distance_propre_axe_ns (const int nr  = 65) const ;
	
};
ostream& operator<<(ostream& , const Bin_ns_bh& ) ;

#endif
