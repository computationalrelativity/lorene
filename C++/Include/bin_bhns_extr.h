/*
 *  Definition of Lorene class Bin_bhns_extr
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

#ifndef __BIN_BHNS_EXTR_H_ 
#define __BIN_BHNS_EXTR_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.2  2005/02/28 23:02:51  k_taniguchi
 * Addition of some codes for the conformally flat case
 *
 * Revision 1.1  2004/11/30 20:36:00  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Lorene header
#include "et_bin_bhns_extr.h"
#include "etoile.h"

/**
 * Class for computing a Black hole - Neutron star binary system
 * with an extreme mass ratio.
 * \ingroup(star)
 * 
 */
class Bin_bhns_extr {

    // Data : 
    // -----
    private:
        /// Cartesian triad of the absolute reference frame
        const Base_vect_cart ref_triad ;

	/// Neutron star
	Et_bin_bhns_extr star ;

	/** Angular velocity with respect to an asymptotically inertial
	 *  observer
	 */
	double omega ;

	/// Absolute orbital separation between two centers of BH and NS
	double separ ;

	/// Gravitational mass of BH
	double mass_bh ;

    // Derived data :
    // ------------
    private :
        /** Absolute coordinate X of the barycenter of the baryon density
	 *  in the Kerr-Schild background metric or in the conformally flat
	 *  one
	 */
	mutable double* p_xa_barycenter_extr ;

        /** Absolute coordinate Y of the barycenter of the baryon density
	 *  in the Kerr-Schild background metric
	 */
	mutable double* p_ya_barycenter_extr ;

        /** Baryon mass of the neutron star in the Kerr-Schild background
	 *  metric or in the conformally flat one
	 */
        mutable double* p_mass_b_extr ;

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor
	 *
	 * @param mp Mapping on which the neutron star will be defined
	 * @param nzet Number of domains occupied by the neutron star
	 * @param eos Equation of state of the neutron star
	 * @param irrot should be {\tt true} if NS is irrotational,
	 *                       {\tt false} if NS is corotating
	 * @param relat should be {\tt true} for a relativistic configuration,
	 *                       {\tt false} for a Newtonian one
	 * @param kerrs should be {\tt true} for the Kerr-Schild background
	 *              metric, {\tt false} for the conformally flat one
	 * @param multi should be {\tt true} for the multipole falloff
	 *              boundary condition, {\tt true} for the \f$1/r\f$ one
	 *
	 */
	Bin_bhns_extr(Map& mp, int nzet, const Eos& eos,
		      bool irrot, bool relat, bool kerrs, bool multi) ;

	Bin_bhns_extr(const Bin_bhns_extr& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE*) )
	 *
	 * @param mp Mapping on which the neutron star will be defined
	 * @param eos Equation of state of the neutron star
	 *
	 */
	Bin_bhns_extr(Map& mp, const Eos& eos, FILE* fich) ;

	~Bin_bhns_extr() ;			///< Destructor
 

    // Memory management
    // -----------------
    private:
	/// Deletes all the derived quantities
	void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Bin_bhns_extr
	void operator=(const Bin_bhns_extr&) ;

	/// Read/write of the neutron star
	Et_bin_bhns_extr& set_ns()
	    { del_deriv() ;
	      return star ; } ;

	/// Sets the orbital angular velocity [{\tt f\_unit}]
	double& set_omega() { return omega ; } ;

	/// Sets the orbital separation [{\tt r\_unit}]
	double& set_separ() { return separ ; } ;

	/// Sets the gravitational mass of BH [{\tt m\_unit}]
	double& set_mass_bh() { return mass_bh ; } ;

    // Accessors
    // ---------
    public:
	/// Returns a reference to the neutron star
	const Et_bin_bhns_extr& get_ns() const
	     { return star ; } ;

	/// Returns the orbital angular velocity [{\tt f\_unit}]
	double get_omega() const { return omega ; } ;

	/** Returns the coordinate separation of the binary system
	 *  [{\tt r\_unit}]
	 */
	double get_separ() const { return separ ; } ;

	/// Returns the gravitational mass of BH [{\tt m\_unit}]
	double get_mass_bh() const { return mass_bh ; } ;

    // Outputs
    // -------
    public:
	void sauve(FILE* ) const ;	    ///< Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Bin_bhns_extr& ) ;	

	/// Display in polytropic units
	void display_poly(ostream& ) const ;

    private:
	/// Operator >> (function called by the operator <<)
	ostream& operator>>(ostream& ) const ;

    // Computational routines
    // ----------------------
    public:
	/** Absolute coordinate X of the barycenter of the baryon density
	 *  in the Kerr-Schild background metric or in the conformally flat
	 *  one
	 */
	double xa_barycenter_extr() const ;

	// Absolute coordinate Y of the barycenter of the baryon density
	/// in the Kerr-Schild background metric
	double ya_barycenter_extr() const ;

	/** Baryon mass of the neutron star in the Kerr-Schild background
	 *  metric or in the conformally flat
	 */
	double mass_b_extr() const ;

	/** Computes the orbital angular velocity {\tt omega}
	 *  in the Kerr-Schild background metric
	 *
	 * @param fact_omeg_min [input] : determines the lower bound of the
	 *             interval {\tt [omega\_min, omega\_max]} in which
	 *             {\tt omega} is searched by
	 *             {\tt omega\_min = fact\_omeg\_min * omega},
	 *             where {\tt omega} is the previous value of the
	 *             angular velocity
	 *             (typical value : {\tt fact\_omeg\_min = 0.5})
	 *
	 * @param fact_omeg_max [input] : determines the higher bound of the
	 *             interval {\tt [omega\_min, omega\_max]} in which
	 *             {\tt omega} is searched by
	 *             {\tt omega\_max = fact\_omeg\_max * omega},
	 *             where {\tt omega} is the previous value of the
	 *             angular velocity.
	 *             (typical value : {\tt fact\_omeg\_max = 1.5})
	 *
	 */
	void orbit_omega_ks(double fact_omeg_min, double fact_omeg_max) ;

	/** Computes the orbital angular velocity {\tt omega}
	 *  in the conformally flat background metric
	 *
	 * @param fact_omeg_min [input] : determines the lower bound of the
	 *             interval {\tt [omega\_min, omega\_max]} in which
	 *             {\tt omega} is searched by
	 *             {\tt omega\_min = fact\_omeg\_min * omega},
	 *             where {\tt omega} is the previous value of the
	 *             angular velocity
	 *             (typical value : {\tt fact\_omeg\_min = 0.5})
	 *
	 * @param fact_omeg_max [input] : determines the higher bound of the
	 *             interval {\tt [omega\_min, omega\_max]} in which
	 *             {\tt omega} is searched by
	 *             {\tt omega\_max = fact\_omeg\_max * omega},
	 *             where {\tt omega} is the previous value of the
	 *             angular velocity.
	 *             (typical value : {\tt fact\_omeg\_max = 1.5})
	 *
	 */
	void orbit_omega_cf(double fact_omeg_min, double fact_omeg_max) ;

	/** Sets the orbital angular velocity to some 2-PN analytical
	 *  value (Keplerian value in the Newtonian case)
	 */
	void analytical_omega() ;

	/** Sets some analytical template for the shift vector
	 *  (via the members {\tt w\_shift} and {\tt khi\_shift}
	 *   of {\tt Etoile\_bin})
	 */
	void analytical_shift() ;

};
ostream& operator<<(ostream& , const Bin_bhns_extr& ) ;

#endif
