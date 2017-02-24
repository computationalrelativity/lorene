/*
 *  Definition of Lorene class Binaire
 *
 */

/*
 *   Copyright (c) 2000-2003 Eric Gourgoulhon
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


#ifndef __BINAIRE_H_ 
#define __BINAIRE_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2017/02/24 15:34:59  j_novak
 * Removal of spurious comments
 *
 * Revision 1.6  2014/10/13 08:52:32  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2009/06/18 18:42:13  k_taniguchi
 * Defined a slightly modified code to determine
 * the orbital angular velocity.
 *
 * Revision 1.4  2003/09/15 15:09:47  e_gourgoulhon
 * Added the member function write_global.
 *
 * Revision 1.3  2002/06/17 14:05:16  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.2  2001/12/20 13:00:31  k_taniguchi
 * Addition of the Komar mass, the virial error by Gourgoulhon and Bonazzola, and the virial error by Friedman, Uryu, and Shibata.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.11  2001/02/28  14:27:32  keisuke
 * *** empty log message ***
 *
 * Revision 2.10  2000/07/07  14:10:04  eric
 * AJout de display_poly.
 *
 * Revision 2.9  2000/03/17  15:26:53  eric
 * Ajout de la fonction analytical_omega.
 *
 * Revision 2.8  2000/03/15  16:43:12  eric
 * Ajout de la fonction analytical_shift().
 *
 * Revision 2.7  2000/03/13  14:24:50  eric
 * Ajout des membres p_ham_constr et p_mom_constr ainsi que des
 * fonctions de calcul associees (verification des equations de contrainte).
 *
 * Revision 2.6  2000/02/18  14:52:05  eric
 * Ajout des membres p_virial et p_total_ener et des fonctions de calcul
 * associees.
 *
 * Revision 2.5  2000/02/12  17:36:39  eric
 * Ajout de la fonction separation().
 *
 * Revision 2.4  2000/02/12  17:09:02  eric
 * Ajout de la fonction orbit.
 *
 * Revision 2.3  2000/02/04  17:15:05  eric
 *  Ajout du membre ref_triad.
 *
 * Revision 2.2  2000/02/02  10:13:08  eric
 * Ajout des fonctions de lecture/ecriture omega, x_axe.
 *
 * Revision 2.1  2000/01/31  17:01:46  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/01/31  15:57:39  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

#include "etoile.h"

namespace Lorene {
/**
 * Binary systems.
 * 
 * @version #$Id$#
 */
class Binaire {

    // Data : 
    // -----
    private: 
	/** Cartesian triad of the absolute reference frame
	 *  
	 */
	const Base_vect_cart ref_triad ; 
	
	/// First star of the system
	Etoile_bin star1 ;
	 
	/// Second star of the system
	Etoile_bin star2 ;
	
	/** Array of the two stars (to perform loops on the stars):
	 *  {\tt et[0]} contains the address of {\tt star1} and {\tt et[1]}
	 *  that of {\tt star2}.
	 */
	Etoile_bin* et[2] ; 
	
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
	 * @param relat should be {\tt true} for a relativistic configuration, 
	 *		{\tt false} for a Newtonian one
	 * 
	 */
	Binaire(Map& mp1, int nzet1, const Eos& eos1, int irrot1, 
		Map& mp2, int nzet2, const Eos& eos2, int irrot2,
		int relat) ;			


	Binaire(const Binaire& ) ;		/// Copy constructor

	/** Constructor from a file (see {\tt sauve(FILE* )}). 
	 * 
	 * @param mp1 Mapping on which {\tt star1} will be defined
	 * @param eos1 Equation of state of {\tt star1}
	 * @param mp2 Mapping on which {\tt star2} will be defined
	 * @param eos2 Equation of state of {\tt star2}
	 * @param fich	input file (must have been created by the function
	 *	{\tt sauve})
	 */
	Binaire(Map& mp1, const Eos& eos1, Map& mp2, const Eos& eos2, 
		FILE* fich) ;			

	~Binaire() ;			/// Destructor
 

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
	/// Assignment to another {\tt Binaire}
	void operator=(const Binaire&) ;	
	
	/// Read/write of the star no. i
	Etoile_bin& set(int i) 
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
	const Etoile_bin& operator()(int i) const 
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
    
	// Display
	friend ostream& operator<<(ostream& , const Binaire& ) ;	
	
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
	 *   {\tt angu\_mom()(0)} = $J^X$, \\
	 *   {\tt angu\_mom()(1)} = $J^Y$, \\
	 *   {\tt angu\_mom()(2)} = $J^Z$. 
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
	 *  (for a relativistic one,
	 *   it returns $|1 - M_{\rm Komar} / M_{\rm ADM}|$)
	 */
    	double virial() const ;	

	/** Estimates the relative error on the virial theorem
	 *  calculated by E.Gourgoulhon and S.Bonazzola
	 *  (Class. Quantum Grav. 11, 443 (1994): Eq.(29))
	 *  normalized by $2 \pi G$.
	 */
    	double virial_gb() const ;	

	/** Estimates the relative error on the virial theorem
	 *  calculated by J.L.Friedman, K.Uryu, and M.Shibata
	 *  (PRD accepted, gr-qc/0108070)
	 *
	 *  The expression used in the LORENE is Eq.(5.7) in the paper
	 *  by M.Shibata and K.Uryu (PRD64, 104017 (2001))
	 */
    	double virial_fus() const ;	

	/** Estimates the relative error on the Hamiltonian constraint 
	 *  equation by comparing $\underline\Delta\ln A$ with
	 *  \begin{equation}
	 *    -4\pi A^2 E - {A^2\over 4} K_{ij} K^{ij} - {1\over 2} 
	 *	    {\overline\nabla}_i \ln A {\overline\nabla}^i \ln A
	 *  \end{equation} 
	 */
    	double ham_constr() const ;	

	/** Estimates the relative error on the momentum constraint 
	 *  equation by comparing ${\overline\nabla}_j K^{ij}$ with
	 *  \begin{equation}
	 *    8\pi J^i - 5 K^{ij} {\overline\nabla}_j \ln A
	 *  \end{equation} 
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
	 *  @param mass1 [input] : baryon rest mass of NS1
	 *
	 *  @param mass2 [input] : baryon rest mass of NS2
	 *
	 *  @param xgg1 [output] : x coordinate (relative to star 1 mapping)
	 *		    of the ``center of mass'' of star 1
	 *
	 *  @param xgg2 [output] : x coordinate (relative to star 2 mapping)
	 *		    of the ``center of mass'' of star 2
	 *
	 */
	void orbit_eqmass(double fact_omeg_min, double fact_omeg_max,
			  double mass1, double mass2,
			  double& xgg1, double& xgg2) ;

	/** Sets the orbital angular velocity to some 2-PN analytical
	 *  value (Keplerian value in the Newtonian case)
	 */
	void analytical_omega() ; 

	/** Sets some analytical template for the shift vector (via the
	 *   members {\tt w\_shift} and {\tt khi\_shift} of the two
	 *   {\tt Etoile\_bin}. 
	 */
	void analytical_shift() ; 

};
ostream& operator<<(ostream& , const Binaire& ) ;	

}
#endif
