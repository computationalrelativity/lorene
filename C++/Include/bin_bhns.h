/*
 *  Definition of Lorene class Bin_bhns
 *
 */

/*
 *   Copyright (c) 2005-2007 Keisuke Taniguchi
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

#ifndef __BIN_BHNS_H_ 
#define __BIN_BHNS_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.2  2008/05/15 18:50:06  k_taniguchi
 * Addition of new global quantities.
 *
 * Revision 1.1  2007/06/22 01:03:50  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// External classes which appear in the declaration of class Bin_bhns:
class Hole_bhns ;
class Star_bhns ;

// Lorene headers
#include "hole_bhns.h"
#include "star_bhns.h"


/**
 * Class for computing a black hole - neutron star binary system
 * with comparable mass
 * \ingroup(star)
 * 
 */
class Bin_bhns {

    // Data : 
    // -----
    protected:
        /// Cartesian triad of the absolute reference frame
        const Base_vect_cart ref_triad ;

        /// Black hole
        Hole_bhns hole ;

	/// Neutron star
	Star_bhns star ;

	/** Angular velocity with respect to an asymptotically inertial
	 *   observer
	 */
	double omega ;

	/// Absolute orbital separation between two centers of BH and NS
	double separ ;

	/// Absolute X coordinate of the rotation axis
	double x_rot ;

	/// Absolute Y coordinate of the rotation axis
	double y_rot ;

    // Derived data : 
    // ------------
    protected:
	/** Total ADM mass of the system calculated by the surface integral
	 *  at infinity
	 */
	mutable double* p_mass_adm_bhns_surf ;

	/** Total ADM mass of the system calculated by the volume integral
	 *  and the surface integral at the apparent horizon
	 */
	mutable double* p_mass_adm_bhns_vol ;

	/** Total Komar mass of the system calculated by the surface integral
	 *  at infinity
	 */
	mutable double* p_mass_kom_bhns_surf ;

	/** Total Komar mass of the system calculated by the volume integral
	 *  and the surface integral at the apparent horizon
	 */
	mutable double* p_mass_kom_bhns_vol ;

	/// Total linear momentum of the system
	mutable Tbl* p_line_mom_bhns ;

	/// Total angular momentum of the system
	mutable Tbl* p_angu_mom_bhns ;

	/** Virial theorem error calculated by the ADM mass and the Komar
	 *  mass of the surface integral at infinity
	 */
	mutable double* p_virial_bhns_surf ;

	/** Virial theorem error calculated by the ADM mass and the Komar
	 *  mass of the volume integral
	 */
	mutable double* p_virial_bhns_vol ;

	/// Absolute coordinate X of the barycenter of the baryon density
	mutable double* p_xa_barycenter ;

	/// Absolute coordinate Y of the barycenter of the baryon density
	mutable double* p_ya_barycenter ;

	/// Orbital angular velocity derived from another method
	mutable double* p_omega_two_points ;

	/// Relative error on the Hamiltonian constraint
	//	mutable double* p_ham_constr_bhns ;

	/// Relative error on the momentum constraint
	//	mutable Tbl* p_mom_constr_bhns ;


    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor
	 *
	 *  @param mp_bh Mapping on which the black hole will be defined
	 *  @param mp_ns Mapping on which the neutron star will be defined
	 *  @param nzet Number of domains occupied by the neutron star
	 *  @param eos Equation of state of the neutron star
	 *  @param irrot_ns should be {\tt true} if NS is irrotational,
	 *                           {\tt false} if NS is corotating
	 *  @param kerrschild should be {\tt true} if the background metric
	 *                    is Kerr-Schild, {\tt false} if the background
	 *                    metric is conformally flat
	 *  @param bc_lapse_nd should be {\tt true} if the BC type for lapse
	 *     is Neumann, {\tt false} if the BC type is Dirichlet
	 *  @param bc_lapse_fs should be {\tt true} if the BC is first type
	 *                              {\tt false} if the BC is second one
	 *  @param irrot_bh should be {\tt true} if BH is irrotational,
	 *                           {\tt false} if NS is corotating
	 *  @param mass_bh Black hole mass which appears in the background
	 *                 metric
	 */
	Bin_bhns(Map& mp_bh, Map& mp_ns, int nzet, const Eos& eos,
		 bool irrot_ns, bool kerrschild,
		 bool bc_lapse_nd, bool bc_lapse_fs, bool irrot_bh,
		 double mass_bh) ;

	Bin_bhns(const Bin_bhns& ) ;		///< Copy constructor

	/// Constructor from a file (see \c sauve(FILE*) )
	Bin_bhns(Map& mp_bh, Map& mp_ns, const Eos& eos, FILE* fich) ;

	virtual ~Bin_bhns() ;			///< Destructor
 

    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Bin_bhns
	void operator=(const Bin_bhns&) ;

	/// Read/write of the black hole
	Hole_bhns& set_bh()
	    { del_deriv() ;
	      return hole ; } ;

	/// Read/write of the neutron star
	Star_bhns& set_ns()
	    { del_deriv() ;
	      return star ; } ;

	/// Sets the orbital angular velocity [{\tt f\_unit}]
	double& set_omega() { return omega ; } ;

	/// Sets the orbital separation [{\tt r\_unit}]
	double& set_separ() { return separ ; } ;

	/// Sets the absolute coordinate X of the rotation axis [{\tt r\_unit}]
	double& set_x_rot() {return x_rot; } ;

	/// Sets the absolute coordinate Y of the rotation axis [{\tt r\_unit}]
	double& set_y_rot() {return y_rot; } ;

	
    // Accessors
    // ---------
    public:
	/// Returns a reference to the black hole
	const Hole_bhns& get_bh() const { return hole ; } ;

	/// Returns a reference to the neutron star
	const Star_bhns& get_ns() const { return star ; } ;

	/// Returns the orbital angular velocity [{\tt f\_unit}]
	double get_omega() const { return omega ; } ;

	/** Returns the coordinate separation of the binary system
	 *  [{\tt r\_unit}]
	 */
	double get_separ() const { return separ ; } ;

	/** Returns the absolute coordinate X of the rotation axis
	 *  [{\tt r\_unit}]
	 */
	double get_x_rot() const {return x_rot; } ;

	/** Returns the absolute coordinate Y of the rotation axis
	 *  [{\tt r\_unit}]
	 */
	double get_y_rot() const {return y_rot; } ;


    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Bin_bhns& ) ;	

	/// Display in polytropic units
	void display_poly(ostream& ) const ;

    private:
	/// Operator >> (function called by the operator <<)
	ostream& operator>>(ostream& ) const ;

    // Computational routines
    // ----------------------
    public:

	/// Total ADM mass
	double mass_adm_bhns_surf() const ;

	double mass_adm_bhns_vol() const ;

	/// Total Komar mass
	double mass_kom_bhns_surf() const ;

	double mass_kom_bhns_vol() const ;

	/** Total linear momentum.
	 *
	 *  @return 1-D {\tt Tbl} of size 3, according to \\
	 *   {\tt line\_mom()(0)} = $P^x$, \\
	 *   {\tt line\_mom()(1)} = $P^y$, \\
	 *   {\tt line\_mom()(2)} = $P^z$.
	 */
	const Tbl& line_mom_bhns() const ;

	/** Total angular momentum.
	 *
	 *  @return 1-D {\tt Tbl} of size 3, according to \\
	 *   {\tt angu\_mom()(0)} = $J^x$, \\
	 *   {\tt angu\_mom()(1)} = $J^y$, \\
	 *   {\tt angu\_mom()(2)} = $J^z$.
	 */
	const Tbl& angu_mom_bhns() const ;

	/** Estimates the relative error on the virial theorem
	 *  $|1 - M_{\rm Komar} / M_{\rm ADM}|$
	 */
	double virial_bhns_surf() const ;

	/** Estimates the relative error on the virial theorem
	 *  $|1 - M_{\rm Komar} / M_{\rm ADM}|$
	 */
	double virial_bhns_vol() const ;

	/// Absolute coordinate X of the barycenter of the baryon density
	double xa_barycenter() const ;

	/// Absolute coordinate Y of the barycenter of the baryon density
	double ya_barycenter() const ;

	/// Orbital angular velocity derived from another method
	double omega_two_points() const ;

	/** Estimates the relative error on the Hamiltonian constraint
	 */
	//	double ham_constr() const ;

	/** Estimates the relative error on the momentum constraint
	 */
	//	const Tbl& mom_constr() const ;

	/** Computes the orbital angular velocity {\tt omega}
	 *
	 *  @para fact_omeg_min [input] : determines the lower bound of the 
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
	 */
	void orbit_omega(double fact_omeg_min, double fact_omeg_max) ;

	/** Computes the position of the rotation axis X
	 *
	 *  @param rot_exp_x [input] : exponent of the factor which modifies
	 *                the position of the two stars from the rotation axis
	 *
	 */
	void rotation_axis_x(double rot_exp_x) ;

	/** Computes the position of the rotation axis Y
	 *
	 *  @param thres_rot [input] : threshold to stop moving to the Y dir.
	 *  @param rot_exp_y [input] : exponent of the factor which modifies
	 *                the Y position of the neutron star coordinate
	 *  @param fact [input] : factor to multiply to Y_NS
	 *
	 */
	void rotation_axis_y(double thres_rot, double rot_exp_y, double fact) ;

	/** Sets some analytical template for the initial shift vector
	 *
	 *  @param reduce_shift_bh [input] : factor to reduce the initial
	 *                                   ansatz
	 *  @param reduce_shift_ns [input] : factor to reduce the initial
	 *                                   ansatz
	 *
	 */
	void shift_analytic(double reduce_shift_bh, double reduce_shift_ns) ;


};
ostream& operator<<(ostream& , const Bin_bhns& ) ;

#endif
