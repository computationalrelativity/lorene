/*
 *  Definition of Lorene class Binary
 *
 */

/*
 *   Copyright (c) 2004  Limousin Francois
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

#ifndef __BINARY_H_ 
#define __BINARY_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/01/20 15:25:39  f_limousin
 * Nex class binary
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "star.h"

/**
 * Binary systems. *** UNDER DEVELOPMENT ***
 * 
 * @version #$Id$#
 */

class Binary {

    // Data : 
    // -----
    protected:

	/// First star of the system
	Star_bin star1 ;
	 
	/// Second star of the system
	Star_bin star2 ;
	
	/** Array of the two stars (to perform loops on the stars):
	 *  {\tt et[0]} contains the address of {\tt star1} and {\tt et[1]}
	 *  that of {\tt star2}.
	 */
	Star_bin* et[2] ; 
	
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
	 * @param mp1 Mapping on which {\tt star1} will be defined
	 * @param nzet1 Number of domains occupied by {\tt star1}
	 * @param eos1 Equation of state of {\tt star1}
	 * @param irrot1 should be {\tt true} if {\tt star1} is irrotational, 
	 *		    {\tt false} if {\tt star1} is corotating
	 * @param mp2 Mapping on which {\tt star2} will be defined
	 * @param nzet2 Number of domains occupied by {\tt star2}
	 * @param eos2 Equation of state of {\tt star2}
	 * @param irrot2 should be {\tt true} if {\tt star2} is irrotational, 
	 *		    {\tt false} if {\tt star2} is corotating
	 * @param conf_flat should be {\tt true} for a 3-metric conformally
	 *            flat and {\tt false} for a more general one.
	 */
	Binary(Map& mp1, int nzet1, const Eos& eos1, int irrot1, 
		   Map& mp2, int nzet2, const Eos& eos2, int irrot2,
		   int conf_flat) ;			


	Binary(const Binary& ) ;      	/// Copy constructor

	/** Constructor from a file (see {\tt sauve(FILE* )}). 
	 * 
	 * @param mp1 Mapping on which {\tt star1} will be defined
	 * @param eos1 Equation of state of {\tt star1}
	 * @param mp2 Mapping on which {\tt star2} will be defined
	 * @param eos2 Equation of state of {\tt star2}
	 * @param fich	input file (must have been created by the function
	 *	{\tt sauve})
	 */
	Binary(Map& mp1, const Eos& eos1, Map& mp2, const Eos& eos2, 
		   FILE* fich) ;			

	~Binary() ;			/// Destructor
 


    // Memory management
    // -----------------
    protected:
	    
	/// Deletes all the derived quantities
	void del_deriv() const ; 
	
	/// Sets to {\tt 0x0} all the pointers on derived quantities
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another  {\tt Binary}
	void operator=(const Binary&) ;	
	
	/// Read/write of the star no. i
	Star_bin& set(int i) 
	    { assert( (i==1) || (i==2) ); 
	      del_deriv() ; 
	      return *et[i-1] ;} ; 

	/// Sets the orbital angular velocity [{\tt f\_unit}] 
	double& set_omega() {return omega; } ; 

	/// Sets the absolute coordinate X of the rotation axis [{\tt r\_unit}]
	double& set_x_axe() {return x_axe; } ; 
	

    // Accessors
    // ---------
    public:
	/// Returns a reference to the star no. i
	const Star_bin& operator()(int i) const 
	    { assert( (i==1) || (i==2) ); 
	      return *et[i-1] ;} ; 

	/// Returns the orbital angular velocity [{\tt f\_unit}]
	double get_omega() const {return omega; } ; 
	
	/// Returns the absolute coordinate X of the rotation axis [{\tt r\_unit}]
	double get_x_axe() const {return x_axe; } ; 

	/** Returns the coordinate separation of the two stellar 
	 *	centers [{\tt r\_unit}]
	 */
	double separation() const ; 


    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    /// Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Binary& ) ;	

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

	/// Total Komar mass
    	double mass_kom() const ;
	
	/** Total angular momentum.
	 *
	 *  @return 1-D {\tt Tbl} of size 3, according to \\
	 *   {\tt angu\_mom()(0)} = $J^r$, \\
	 *   {\tt angu\_mom()(1)} = $J^t$, \\
	 *   {\tt angu\_mom()(2)} = $J^p$. 
	 */
    	const Tbl& angu_mom() const ;	

	/** Total energy (excluding the rest mass energy).
	 * 
	 *  In the Newtonian case, it is defined as the sum of kinetic, 
	 *  internal, and gravitational potential energies. 
	 * 
	 *  In the relativistic case, it is defined as 
	 *  $M_{\rm ADM} - M_{\rm bar,1} - M_{\rm bar,2}$.
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

	/** Computes the orbital angular velocity {\tt omega} and the 
	 *  position of the rotation axis {\tt x\_axe}. 
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
	 *   members {\tt w\_shift} and {\tt khi\_shift} of the two
	 *   {\tt Etoile\_bin}. 
	 */
	void analytical_shift() ; 

	/**
	 * Calculates {tt decouple} which is used to obtain 
	 * {\tt qq\_auto} by the formula : 
	 * {\tt qq\_auto} = {\tt decouple} * {\tt qq}.
	 * (see the membre {tt Scalar decouple} for more 
	 * precisions about its value).
	 * 
	 */
	void fait_decouple () ;
 
};
ostream& operator<<(ostream& , const Binary& ) ;	

#endif
