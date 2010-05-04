/*
 *  Definition of Lorene class Binary_xcts
 *
 */

/*
 *   Copyright (c) 2010  Michal Bejger
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
 
#ifndef __BINARY_XCTS_H_ 
#define __BINARY_XCTS_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.1  2010/05/04 07:44:02  m_bejger
 * Initial version
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "star.h"

/**
 * Binary systems in eXtended Conformal Thin Sandwich formulation. 
 * *** UNDER DEVELOPMENT *** \ingroup (star)
 * 
 */

class Binary_xcts {

    // Data : 
    // -----
    protected:

	/// First star of the system
	Star_bin_xcts star1 ;
	 
	/// Second star of the system
	Star_bin_xcts star2 ;
	
	/** Array of the two stars (to perform loops on the stars):
	 *  \c et[0] contains the address of \c star1 and \c et[1]
	 *  that of \c star2.
	 */
	Star_bin_xcts* et[2] ; 
	
	/** Angular velocity with respect to an asymptotically inertial 
	 *  observer
	 */
	double omega ;
	
	/** Absolute X coordinate of the rotation axis
	 */
	double x_axe ;


    // Derived data : 
    // ------------
    protected:
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

	/// Relative error on the Hamiltonian constraint
	mutable double* p_ham_constr ;
	
	/// Relative error on the momentum constraint
	mutable Tbl* p_mom_constr ;

	

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp1 Mapping on which \c star1 will be defined
	 * @param nzet1 Number of domains occupied by \c star1
	 * @param eos1 Equation of state of \c star1
	 * @param irrot1 should be \c true if \c star1 is irrotational, 
	 *		    \c false if \c star1 is corotating
	 * @param mp2 Mapping on which \c star2 will be defined
	 * @param nzet2 Number of domains occupied by \c star2
	 * @param eos2 Equation of state of \c star2
	 * @param irrot2 should be \c true if \c star2 is irrotational, 
	 *		    \c false if \c star2 is corotating
	 */
	Binary_xcts(Map& mp1, int nzet1, const Eos& eos1, int irrot1, 
		        Map& mp2, int nzet2, const Eos& eos2, int irrot2) ;			


	Binary_xcts(const Binary_xcts& ) ;      	///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE* )). 
	 * 
	 * @param mp1 Mapping on which \c star1 will be defined
	 * @param eos1 Equation of state of \c star1
	 * @param mp2 Mapping on which \c star2 will be defined
	 * @param eos2 Equation of state of \c star2
	 * @param fich	input file (must have been created by the function
	 *	\c sauve)
	 */
	Binary_xcts(Map& mp1, const Eos& eos1, Map& mp2, const Eos& eos2, 
		   		FILE* fich) ;			

	~Binary_xcts() ;			///< Destructor
 


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
	/// Assignment to another  \c Binary_xcts
	void operator=(const Binary_xcts&) ;	
	
	/// Read/write of the star no. i
	Star_bin_xcts& set(int i) { 
			assert( (i==1) || (i==2) ); 
	      	del_deriv() ; 
	      	return *et[i-1] ;
	      	
	} ; 

	/// Sets the orbital angular velocity [\c f_unit] 
	double& set_omega() {return omega; } ; 

	/// Sets the absolute coordinate X of the rotation axis [\c r_unit]
	double& set_x_axe() {return x_axe; } ; 
	

    // Accessors
    // ---------
    public:
	/// Returns a reference to the star no. i
	const Star_bin_xcts& operator()(int i) const 
	    { assert( (i==1) || (i==2) ); 
	      return *et[i-1] ;} ; 

	/// Returns the orbital angular velocity [\c f_unit]
	double get_omega() const {return omega; } ; 
	
	/// Returns the absolute coordinate X of the rotation axis [\c r_unit]
	double get_x_axe() const {return x_axe; } ; 

	/** Returns the coordinate separation of the two stellar 
	 *	centers [\c r_unit]
	 */
	double separation() const ; 


    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    ///< Save in a file
    
	/// Display
	friend ostream& operator<<(ostream&, const Binary_xcts& ) ;	

	/// Display in polytropic units
	void display_poly(ostream& ) const ; 

	/** Write global quantities in a formatted file. 
	 * This file can be read by an external program. 
	 */
	void write_global(ostream& ) const  ;

      private:
	/// Operator >> (function called by the operator <<). 
	ostream& operator>>(ostream& ) const ;    


    // Computational routines
    // ----------------------
    public:

	/// Total ADM mass
    	double mass_adm() const ;

	/// Total ADM mass (computed by a volume integral)
    	double mass_adm_vol() const ;

	/// Total Komar mass
    	double mass_kom() const ;
	
	/// Total Komar mass (computed by a volume integral)
    	double mass_kom_vol() const ;

	/** Total angular momentum.
	 *
	 *  @return 1-D \c Tbl of size 3, according to
	 *   \li \c angu_mom()(0) = \f$J^r\f$, 
	 *   \li \c angu_mom()(1) = \f$J^t\f$, 
	 *   \li \c angu_mom()(2) = \f$J^p\f$. 
	 */
    	const Tbl& angu_mom() const ;	

	/** Total energy (excluding the rest mass energy).
	 * 
	 *  In the Newtonian case, it is defined as the sum of kinetic, 
	 *  internal, and gravitational potential energies. 
	 * 
	 *  In the relativistic case, it is defined as 
	 *  \f$M_{\rm ADM} - M_{\rm bar,1} - M_{\rm bar,2}\f$.
	 */
    	double total_ener() const ;	

	/** Estimates the relative error on the virial theorem
	 */
    	double virial() const ;	

	/** Estimates the relative error on the Hamiltonian constraint 
	 */
    	double ham_constr() const ;	

	/** Estimates the relative error on the momentum constraint 
	 */
    	const Tbl& mom_constr() const ;	

	/** Computes the orbital angular velocity \c omega and the 
	 *  position of the rotation axis \c x_axe. 
	 *
	 *  @param fact_omeg_min [input] : determines the lower bound of the 
	 *		interval \c [omega_min, omega_max] in which 
	 *		\c omega is searched by 
	 *		\c omega_min = fact_omeg_min * omega, 
	 *		where \c omega is the previous value of the 
	 *		angular velocity 
	 *		(typical value : \c fact_omeg_min = 0.5)
	 *
	 *  @param fact_omeg_max [input] : determines the higher bound of the 
	 *		interval \c [omega_min, omega_max] in which 
	 *		\c omega is searched by 
	 *		\c omega_max = fact_omeg_max * omega, 
	 *		where \c omega is the previous value of the 
	 *		angular velocity.
	 *		(typical value : \c fact_omeg_max = 1.5)
	 *
	 *  @param xgg1 [output] : x coordinate (relative to star 1 mapping)
	 *		    of the ``center of mass'' of star 1
	 *
	 *  @param xgg2 [output] : x coordinate (relative to star 2 mapping)
	 *		    of the ``center of mass'' of star 2
	 *
	 */
	void orbit(double fact_omeg_min, double fact_omeg_max, double& xgg1, 
		   double& xgg2) ;
	   
	/** Sets the orbital angular velocity to some 2-PN analytical
	 *  value (Keplerian value in the Newtonian case)
	 */
	void analytical_omega() ; 

	/** Sets some analytical template for the shift vector (via the
	 *   members \c w_shift and \c khi_shift of the two
	 *   \c Star_bin. 
	 */
	void analytical_shift() ; 

	/**
	 * Calculates \c decouple which is used to obtain 
	 * \c qq_auto by the formula : 
	 * \c qq_auto = \c decouple * \c qq.
	 * (see the membre \c Scalar \c decouple for more 
	 * precisions about its value).
	 * 
	 */
	void fait_decouple() ;
 
} ;

ostream& operator<<(ostream& , const Binary_xcts& ) ;	

#endif
