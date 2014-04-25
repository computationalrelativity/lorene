/*
 *  Definition of Lorene classes Eos_bifluid
 *				 Eos_bf_poly
 *
 */

/*
 *   Copyright (c) 2001-2002 Jerome Novak
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


#ifndef __EOS_BIFLUID_H_ 
#define __EOS_BIFLUID_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.18  2014/04/25 10:43:50  j_novak
 * The member 'name' is of type string now. Correction of a few const-related issues.
 *
 * Revision 1.17  2004/09/01 09:49:46  r_prix
 * adapted to change in read_variable() for strings
 *
 * Revision 1.16  2004/03/22 13:12:41  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.15  2004/01/30 13:21:29  r_prix
 * add documentation about 'special' 2-fluid typeos=5: == type0 + slow-rot style inversion
 *
 * Revision 1.14  2003/12/17 23:12:30  r_prix
 * replaced use of C++ <string> by standard ANSI char* to be backwards compatible
 * with broken compilers like MIPSpro Compiler 7.2 on SGI Origin200. ;-)
 *
 * Revision 1.13  2003/12/05 15:08:38  r_prix
 * - use read_variable() to read eos_bifluid from file
 * - changed 'contructor from file' to take filename as an argument instead of ifstream
 * - changed 'name' member of Eos_bifluid to C++-type string (for convenience&safety)
 *
 * Revision 1.12  2003/12/04 14:13:32  r_prix
 * added method get_typeos {return typeos}; and fixed some comments.
 *
 * Revision 1.11  2003/11/18 18:25:15  r_prix
 * moved particle-masses m_1, m_2 of the two fluids into class eos_bifluid (from eos_bf_poly)
 *
 * Revision 1.10  2002/10/16 14:36:29  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.9  2002/09/13 09:17:31  j_novak
 * Modif. commentaires
 *
 * Revision 1.8  2002/06/17 14:05:16  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.7  2002/06/03 13:23:16  j_novak
 * The case when the mapping is not adapted is now treated
 *
 * Revision 1.6  2002/05/31 16:13:36  j_novak
 * better inversion for eos_bifluid
 *
 * Revision 1.5  2002/05/02 15:16:22  j_novak
 * Added functions for more general bi-fluid EOS
 *
 * Revision 1.4  2002/01/16 15:03:27  j_novak
 * *** empty log message ***
 *
 * Revision 1.3  2002/01/11 14:09:34  j_novak
 * Added newtonian version for 2-fluid stars
 *
 * Revision 1.2  2001/11/29 15:05:26  j_novak
 * The entrainment term in 2-fluid eos is modified
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.6  2001/08/31  15:47:35  novak
 * The flag tronc has been added to nbar_ent.. functions
 *
 * Revision 1.5  2001/08/27 12:21:11  novak
 * The delta2 Cmp argument put to const
 *
 * Revision 1.4  2001/08/27 09:50:15  novak
 * New formula for "polytrope"
 *
 * Revision 1.3  2001/06/22 15:36:11  novak
 * Modification de Eos_bifluid::trans2Eos
 *
 * Revision 1.2  2001/06/22 11:52:44  novak
 * *** empty log message ***
 *
 * Revision 1.1  2001/06/21 15:21:22  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Standard C++
#include "headcpp.h"
#include <string>

// Headers C
#include <cstdio>

// Lorene classes
#include "param.h"
class Tbl ;
class Param ;
class Cmp ;
class Eos ;
class Eos_poly ;

		    //------------------------------------//
		    //   base class Eos for two fluids	  //
		    //------------------------------------//
#define MAX_EOSNAME 100

/**
 * 2-fluids equation of state base class. 
 *
 * Fluid 1 is supposed to correspond to neutrons, whereas fluid 2
 * corresponds to e.g. protons. Neutron 4-velocity is \f$u^\alpha_{\rm n}\f$
 * and proton one is \f$u^\alpha_{\rm p}\f$
 *
 * Therefore, the EOS is defined by giving two log-enthalpies AND a 
 * relative velocity as inputs. The output are then: two baryonic
 * densities, the total energy density and pressure
 * The enthalpies \f$f_1\f$ and \f$f_2\f$ are obtained through the formula
 * \f[
 * {\rm d}{\cal E}=f^1{\rm d}n_1+f^2{\rm d}n_2+\alpha{\rm d}\Delta^2
 * \label{eeosbfdefent}
 * \f]
 * see Comer, Novak \& Prix. Log-enthalpies are then defined as
 * \f$H_1 = \ln\left( \frac{f^1}{m_1c^2} \right)\f$, where \f$m_1\f$ is the
 * mass of a particle of the first fluid. (same for \f$H_2\f$)
 * The relative velocity \f$\Delta^2\f$
 * is defined as in Comer, Novak \& Prix. It can be seen as
 * the neutron velocity seen in the frame of protons (or vice-versa): 
 * \f$\Gamma_\Delta = -g_{\alpha\beta} u^\alpha_{\rm n} u^\beta_{\rm p}\f$
 * \ingroup (eos)
 */
class Eos_bifluid {

    // Data :
    // -----

    protected: 
	string name;	    ///< EOS name

	/** Individual particle mass \f$m_1\f$  
	 *  [unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$]. 
	 */
	double m_1 ; 

	/** Individual particle mass \f$m_2\f$  
	 *  [unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$]. 
	 */
	double m_2 ; 


    // Constructors - Destructor
    // -------------------------
    protected:
	Eos_bifluid() ;			///< Standard constructor

	/// Standard constructor with name and the two masses per particle
	explicit Eos_bifluid(const char* name_i, double mass1, double mass2) ; 

	Eos_bifluid(const Eos_bifluid& ) ;	///< Copy constructor	

    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 * \c Eos_bifluid::eos_from_file(FILE*) . 
	 */
	Eos_bifluid(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Eos_bifluid::eos_from_file(const char*). 
	 *
	 *  The following fields have to be present in the config-file:\\
	 *  name: [string] name of the EOS
	 *  m_1, m_2: [double] baryon masses of the 2-fluids
	 *
	 */
	Eos_bifluid (const char *fname ) ; 
	
	
    public:
	virtual ~Eos_bifluid() ;			///< Destructor

    // Assignment
    // ----------
	/// Assignment to another \c Eos_bifluid 
	void operator=(const Eos_bifluid& ) ;


    // Name manipulation
    // -----------------
    public:
	/// Returns the EOS name
	string get_name() const {return name;} ; 

    // Miscellaneous
    // -------------
    public:

	/** Return the individual particule mass \f$m_1\f$  
	 *  
	 *  [unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$]. 
	 */
	double get_m1() const {return m_1 ;}; 

	/** Return the individual particule mass \f$m_2\f$  
	 *  
	 *  [unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$]. 
	 */
	double get_m2() const {return m_2 ;}; 

	/** Construction of an EOS from a binary file.
	 *  The file must have been created by the function 
	 *  \c sauve(FILE*) .
	 */
	static Eos_bifluid* eos_from_file (FILE* ) ; 
	
	/** Construction of an EOS from a formatted file.
	 * 
	 *  The following field has to be present:\\
	 *  ident: [int] identifying the type of 2-fluid EOS
	 *	1 = relativistic polytropic EOS (class \c Eos_bf_poly ). \\
	 *      2 = Newtonian polytropic EOS (class \c Eos_bf_poly_newt ).
	 */
	static Eos_bifluid* eos_from_file ( const char *fname ) ; 
	
	/// Comparison operator (egality)
	virtual bool operator==(const Eos_bifluid& ) const = 0 ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos_bifluid& ) const = 0 ; 
    
	/** Returns a number to identify the sub-classe of \c Eos_bifluid  
	 *  the object belongs to. 
	 */
	virtual int identify() const = 0 ; 

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	///< Save in a file

	/// Display
	friend ostream& operator<<(ostream& , const Eos_bifluid& ) ;	

    protected: 
	virtual ostream& operator>>(ostream &) const = 0 ;    ///< Operator >>


    // Computational functions
    // -----------------------
    public:
	/**  General computational method for \c Cmp 's, it computes
	 *   both baryon densities, energy and pressure profiles.
	 * 
	 *  @param ent1 [input] the first log-enthalpy field \f$H_1\f$.  
	 *  @param ent2 [input] the second log-enthalpy field \f$H_2\f$.
	 *  @param delta2 [input] the relative velocity field \f$\Delta^2 \f$
	 *  @param nbar1 [output] baryonic density of the first fluid
	 *  @param nbar2 [output] baryonic density of the second fluid
	 *  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *  @param ener [output] total energy density \f$\cal E\f$ 
	 *                             of both fluids together
	 *  @param press [output] pressure \e p  of both fluids together
	 *  @param nzet  [input] number of domains where \c resu  is to be
	 *	computed. 
	 *  @param l_min [input] index of the innermost domain is which 
	 *      \c resu  is to be computed [default value: 0]; 
	 *      \c resu  is computed only in domains whose indices are 
	 *      in \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero. 
	 */
	void calcule_tout(const Cmp& ent1, const Cmp& ent2, const Cmp& delta2,
			  Cmp& nbar1, Cmp& nbar2, Cmp& ener, Cmp& press,
			  int nzet, int l_min = 0) const ; 

 	/** Computes both baryon densities from the log-enthalpies 
	 *  (virtual function implemented in the derived classes). 
	 * 
	 *  @param ent1 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_1\f$ 
	 *  @param ent2 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_2\f$ 
	 *  @param delta2 [input,  unit: \f$c^2\f$] relative velocity \f$\Delta^2\f$ 
	 * 
	 *  @param nbar1 [output] baryonic density of the first fluid
	 *  @param nbar2 [output] baryonic density of the second fluid
	 *  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *  @return true if the 2-fluids model is correct, false otherwise.
	 */
    	virtual bool nbar_ent_p(const double ent1, const double ent2, 
				const double delta2, double& nbar1, 
				double& nbar2) const = 0 ;

	/** Computes baryon density out of the log-enthalpy asuming
	 *  that only fluid 1 is present (virtual function implemented 
	 *  in the derived classes).
	 *  @param ent1 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_1\f$ 
	 *  @return nbar1 baryonic density of the first fluid
	 */
	virtual double nbar_ent_p1(const double ent1) const = 0 ;

	/** Computes baryon density out of the log-enthalpy assuming
	 *  that only fluid 2 is present (virtual function implemented 
	 *  in the derived classes).
	 *  @param ent2 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_1\f$ 
	 *  @return nbar1 baryonic density of the first fluid
	 */
	virtual double nbar_ent_p2(const double ent2) const = 0 ;

	/** Computes both baryon density fields from the log-enthalpy fields
	 *  and the relative velocity.
	 * 
	 *  @param ent1 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_1\f$ 
	 *  @param ent2 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_2\f$ 
	 *  @param delta2 [input,  unit: \f$c^2\f$] relative velocity \f$\Delta^2\f$ 
	 *
	 *  @param nbar1 [output] baryonic density of the first fluid
	 *  @param nbar2 [output] baryonic density of the second fluid
	 *  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *  @param nzet  number of domains where the baryon density is to be
	 *	computed. 
	 *  @param l_min  index of the innermost domain is which the baryon 
	 *	density is
	 *	to be computed [default value: 0]; the baryon density is 
	 *	computed only in domains whose indices are in 
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero. 
	 * 
	 */
    	void nbar_ent(const Cmp& ent1, const Cmp& ent2, const Cmp& delta2, 
		      Cmp& nbar1, Cmp& nbar2, int nzet, int l_min = 0) 
	  const  ; 
    
 	/** Computes the total energy density from the baryonic densities
	 *  and the relative velocity.
	 *  (virtual function implemented in the derived classes). 
	 * 
	 *  @param nbar1 [input] baryonic density of the first fluid
	 *  @param nbar2 [input] baryonic density of the second fluid
	 *  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *  @param delta2 [input,  unit: \f$c^2\f$] relative velocity \f$\Delta^2\f$ 
	 * 
	 *  @return energy density \f$\cal E\f$ [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_nbar_p(const double nbar1, const double nbar2, 
				   const double delta2) const = 0 ; 
       
	/** Computes the total energy density from the log-enthalpy fields
	 *  and the relative velocity.
	 * 
	 *  @param ent1 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_1\f$ 
	 *  @param ent2 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_2\f$ 
	 *  @param delta2 [input,  unit: \f$c^2\f$] relative velocity \f$\Delta^2\f$ 
	 *  @param nzet  number of domains where the energy density is to be
	 *	computed. 
	 *  @param l_min  index of the innermost domain is which the energy 
	 *	density is
	 *	to be computed [default value: 0]; the energy density is 
	 *	computed only in domains whose indices are in 
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero. 
	 * 
	 *  @return  energy density field [unit: \f$\rho_{\rm nuc} c^2\f$], 
	 *      where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	Cmp ener_ent(const Cmp& ent1, const Cmp& ent2, const Cmp& delta2, 
		      int nzet, int l_min = 0) const ; 
    
 	/** Computes the pressure from the baryonic densities
	 *  and the relative velocity.
	 *  (virtual function implemented in the derived classes). 
	 * 
	 *  @param nbar1 [input] baryonic density of the first fluid
	 *  @param nbar2 [input] baryonic density of the second fluid
	 *  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *  @param delta2 [input,  unit: \f$c^2\f$] relative velocity \f$\Delta^2\f$ 
	 * 
	 *  @return pressure \e p  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_nbar_p(const double nbar1, const double nbar2, 
				    const double delta2) const = 0 ; 
       
	/** Computes the pressure from the log-enthalpy fields
	 *  and the relative velocity.
	 * 
	 *  @param ent1 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_1\f$ 
	 *  @param ent2 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_2\f$ 
	 *  @param delta2 [input,  unit: \f$c^2\f$] relative velocity \f$\Delta^2\f$ 
	 *  @param nzet  number of domains where the pressure is to be
	 *	computed. 
	 *  @param l_min  index of the innermost domain is which the pressure
	 *    is to be computed [default value: 0]; the pressure is computed 
	 *      only in domains whose indices are in 
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero. 
	 * 
	 *  @return pressure field [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 * 
	 */
    	Cmp press_ent(const Cmp& ent1, const Cmp& ent2, const Cmp& delta2, 
		      int nzet, int l_min = 0) const ; 

	/** Computes the derivative of the energy with respect to
	 * (baryonic density 1)\f$^2\f$.
	 *  (virtual function implemented in the derived classes). 
	 *
	 *  @param n1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *           relative Lorentz factor\f$\times\f$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative \f$K^{11}=2\frac{\partial{\cal{E}}}{\partial{n_1^2}}\f$ 
	 */
	virtual double get_K11(const double n1, const double n2, const
			       double x)  const = 0 ;

	/** Computes the derivative of the energy with respect to 
	 *  \f$x^2=n_1n_2\Gamma_\Delta\f$.
	 *  (virtual function implemented in the derived classes). 
	 *
	 *  @param n1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *           relative Lorentz factor\f$\times\f$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative \f$K^{12}=\frac{\partial {\cal E}}{\partial (n_1n_2\Gamma_\Delta)}\f$ 
	 */
	virtual double get_K12(const double n1, const double n2,const
			       double x) const = 0 ;

	/** Computes the derivative of the energy/(baryonic density 2)\f$^2\f$.
	 *  (virtual function implemented in the derived classes). 
	 *
	 *  @param n1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *           relative Lorentz factor\f$\times\f$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative \f$K^{22} = 2\frac{\partial {\cal E}}{\partial n_2^2}\f$ 
	 */
	virtual double get_K22(const double n1, const double n2, const
			       double x) const = 0 ;

	/** Computes the derivatives of the energy/(baryonic density 1)\f$^2\f$.
	 *
	 *  @param nbar1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density field of fluid 1 at which 
	 *           the derivatives are to be computed
	 *  @param nbar2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density field of fluid 2 at which 
	 *           the derivatives are to be computed
	 *  @param x2 [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *             relative velocity\f$\times\f$both densities at which 
	 *           the derivative is to be computed	
	 *  @param nzet  number of domains where the derivatives are to be
	 *	computed. 
	 *  @param l_min  index of the innermost domain is which the
	 *	derivatives are
	 *	to be computed [default value: 0]; the derivatives are
	 *	computed only in domains whose indices are in 
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero. 
	 * 
	 *  @return derivative \f$K^{11}\f$ field (see \c get_K11 ) 
	 */
    	Cmp get_Knn(const Cmp& nbar1, const Cmp& nbar2, const Cmp& x2,
		    int nzet, int l_min = 0) const  ; 

	/** Computes the derivatives of the energy/(baryonic density 2)\f$^2\f$.
	 *
	 *  @param nbar1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density field of fluid 1 at which 
	 *           the derivatives are to be computed
	 *  @param nbar2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density field of fluid 2 at which 
	 *           the derivatives are to be computed
	 *  @param x2 [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *             relative velocity\f$\times\f$both densities at which 
	 *           the derivative is to be computed	
	 *  @param nzet  number of domains where the derivatives are to be
	 *	computed. 
	 *  @param l_min  index of the innermost domain is which the
	 *	derivatives are
	 *	to be computed [default value: 0]; the derivatives are
	 *	computed only in domains whose indices are in 
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero. 
	 * 
	 *  @return derivative \f$K^{22}\f$ field (see \c get_K12 ) 
	 */
    	Cmp get_Kpp(const Cmp& nbar1, const Cmp& nbar2, const Cmp&
		    x2, int nzet, int l_min = 0) const ; 

	/** Computes the derivatives of the energy with respect to
	 *  \f$x^2=n_1n_2\Gamma_\Delta^2\f$.
	 *
	 *  @param nbar1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density field of fluid 1 at which 
	 *           the derivatives are to be computed
	 *  @param nbar2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density field of fluid 2 at which 
	 *           the derivatives are to be computed
	 *  @param x2 [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *             relative velocity\f$\times\f$both densities at which 
	 *           the derivative is to be computed	
	 *  @param nzet  number of domains where the derivatives are to be
	 *	computed. 
	 *  @param l_min  index of the innermost domain is which the
	 *	derivatives are
	 *	to be computed [default value: 0]; the derivatives are
	 *	computed only in domains whose indices are in 
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero. 
	 * 
	 *  @return derivative \f$K^{12}\f$ field (see \c get_K12 ) 
	 */
     	Cmp get_Knp(const Cmp& nbar1, const Cmp& nbar2, const Cmp& x2,
		    int nzet, int l_min = 0) const ; 

	/**  General computational method for \c Cmp 's (\f$K^{ij}\f$'s).
	 * 
	 *  @param nbar1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density field of fluid 1 at which 
	 *           the derivatives are to be computed.
	 *  @param nbar2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density field of fluid 2 at which 
	 *           the derivatives are to be computed
	 *  @param x2 [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *             relative velocity\f$\times\f$both densities at which 
	 *           the derivative is to be computed	
	 *  @param nzet  [input] number of domains where \c resu  is to be
	 *	computed. 
	 *  @param l_min [input] index of the innermost domain is which \c resu 
	 *	is to be computed [default value: 0]; \c resu  is 
	 *	computed only in domains whose indices are in 
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero. 
	 *  @param fait [input] pointer on the member function of class 
	 *	  \c Eos_bifluid  which performs the pointwise calculation.
	 *  @param resu [output] result of the computation. 
	 */
	void calcule(const Cmp& nbar1, const Cmp& nbar2, const Cmp&
		     x2, int nzet, int l_min, double
		     (Eos_bifluid::*fait)(double, double, double) const, 
		     Cmp& resu)
	  const ; 
 	
   // Conversion functions
    // ---------------------

	/** Makes a translation from \c Eos_bifluid  to \c Eos . 
	 *  (virtual function implemented in the derived classes). 
	 *
	 *  This is only useful for the construction of a 
	 *  \c Et_rot_bifluid 
	 *  star and ought not to be used in other situations.
	 *  @param relat [input] Relativistic EOS or not. 
	 */
	virtual Eos* trans2Eos() const = 0 ;
    
};
ostream& operator<<(ostream& , const Eos_bifluid& ) ;	


		    //------------------------------------//
		    //	      class Eos_bf_poly		  //
		    //------------------------------------//

/**
 * Analytic equation of state for two fluids (relativistic case).
 * 
 * This equation of state (EOS) corresponds to two types of relativistic
 * particles of rest mass is \f$m_1\f$ and \f$m_2\f$,  whose total energy density 
 * \f$\cal{E}\f$ is related to their numerical densities \f$n_1\f$, \f$n_2\f$ and 
 * relative velocity 
 * \f[
 * \Gamma_\Delta = -g_{\alpha\beta} u^\alpha_{\rm n} u^\beta_{\rm p}
 * \label{e:defgamamdelta}
 * \f]
 * (\f$u^\alpha_{\rm n}\f$ and \f$u^\alpha_{\rm p}\f$ being the 4-velocities of 
 * both fluids), by
 * \f[ \label{eeosbfpolye}
 *   {\cal{E}} = \frac{1}{2}\kappa_1 n_1^{\gamma_1} + m_1 c^2 \, n_1 
 *          \  + \frac{1}{2}\kappa_2 n_2^{\gamma_2} + m_2 c^2 \, n_2 
 *          \  + \kappa_3 n_1^{\gamma_3} n_2^{\gamma_4} 
 *          \  + \Delta^2 \beta n_1^{\gamma_5} n_2^{\gamma_6}\ .  
 * \f]
 * The relativistic (i.e. including rest mass energy) chemical potentials
 * are then
 * \f[ 
 * \mu_1 := {{\rm d}{\cal{E}} \over {\rm d}n_1} = \frac{1}{2}\gamma_1\kappa_1
 *         n_1^{\gamma_1-1} + m_1 c^2 + \gamma_3 \kappa_3 
 *         n_1^{\gamma_3-1} n_2^{\gamma_4} + \Delta^2 \gamma_5 \beta 
 *         n_1^{\gamma_5-1} n_2^{\gamma_6}\ , 
 * \f]
 * \f[
 * \mu_2 := {{\rm d}{\cal{E}} \over {\rm d}n_2} = \frac{1}{2}\gamma_2\kappa_2
 *         n_2^{\gamma_2-1} + m_2 c^2 + \gamma_4 \kappa_3 
 *         n_1^{\gamma_3} n_2^{\gamma_4-1} + \Delta^2 \gamma_6 \beta 
 *         n_1^{\gamma_5} n_2^{\gamma_6-1} \ .
 * \f]
 * The pressure is given by the (zero-temperature) First Law of Thermodynamics:
 * \f$p = \mu_1 n_1 + \mu_2 n_2 - {\cal E}\f$, so that
 * \f[ 
 *   p = \frac{1}{2} (\gamma_1 -1)\kappa_1 n_1^{\gamma_1} +
 *  \frac{1}{2}(\gamma_2-1)\kappa_2 n_2^{\gamma_2} + (\gamma_3 +\gamma_4
 *  -1)\kappa_3 n_1^{\gamma_3}n_2^{\gamma_4} + \Delta^2 \beta \left(
 *  (\gamma_5 + \gamma_6 - 1) n_1^{\gamma_5} n_2^{\gamma_6} \right) \ .  
 * \f]
 * The log-enthalpies are defined as the logarithm of the ratio of the enthalpy
 * per particle (see Eq.~\ref{eeosbfdefent}) by the particle rest mass energy :
 * \f[\label{eeosbfentanal} 
 *   H_i := c^2 \ln \left( \frac{f_i}{m_ic^2} \right)   \ .  
 * \f]
 * From this system, the particle densities are obtained in term of 
 * the log-enthalpies. (The system (\ref{eeosbfentanal}) is a linear one
 * if \f$\gamma_1 = \gamma_2 = 2\f$ and \f$\gamma_3 = \gamma_4 = \gamma_5 =
 * \gamma_6 = 1\f$).
 *
 * The energy density \f$\cal E\f$and pressure \e p  can then be obtained
 * as functions of baryonic densities. \ingroup (eos)
 *
 */
class Eos_bf_poly : public Eos_bifluid {

    // Data :
    // -----

    protected: 
	/// Adiabatic indexes \f$\gamma_1\f$, see Eq.~\ref{eeosbfpolye}
	double gam1 ;

	/// Adiabatic indexes \f$\gamma_2\f$, see  Eq.~\ref{eeosbfpolye}
	double gam2 ;

	/// Adiabatic indexes \f$\gamma_3\f$, see  Eq.~\ref{eeosbfpolye}
	double gam3 ;
	
	/// Adiabatic indexes \f$\gamma_4\f$, see Eq.~\ref{eeosbfpolye}
	double gam4 ;

	/// Adiabatic indexes \f$\gamma_5\f$, see  Eq.~\ref{eeosbfpolye}
	double gam5 ;

	/// Adiabatic indexes \f$\gamma_6\f$, see  Eq.~\ref{eeosbfpolye}
	double gam6 ;
	
	/** Pressure coefficient \f$\kappa_1\f$  , see Eq.~\ref{eeosbfpolye}
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$. 
	 */
	double kap1 ; 

	/** Pressure coefficient \f$\kappa_2\f$  , see Eq.~\ref{eeosbfpolye}
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$. 
	 */
	double kap2 ; 

	/** Pressure coefficient \f$\kappa_3\f$  , see Eq.~\ref{eeosbfpolye}
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$. 
	 */
	double kap3 ; 
	
	/** Coefficient \f$\beta\f$  , see Eq.~\ref{eeosbfpolye}
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$. 
	 */
	double beta ; 

	double gam1m1 ;	    ///< \f$\gamma_1-1\f$
	double gam2m1 ;	    ///< \f$\gamma_2-1\f$
	double gam34m1 ;    ///< \f$\gamma_3+\gamma_4-1\f$
	double gam56m1 ;    ///< \f$\gamma_5+\gamma_6-1\f$

 protected:
	/** The bi-fluid analytical EOS type:
	 * 
	 *  0 - \f$\gamma_1 = \gamma_2 = 2\f$ and 
	 *  \f$\gamma_3 = \gamma_4 = \gamma_5 = \gamma_6 = 1\f$. In this case, 
	 *  the EOS can be inverted analytically.
	 *
	 *  1 - \f$\gamma_3 = \gamma_4 = \gamma_5 = \gamma_6 = 1\f$, but
	 *  \f$\gamma_1 \not= 2\f$ or \f$\gamma_2 \not= 2\f$. 
	 *
	 *  2 - \f$\gamma_3 = \gamma_5 = 1\f$, but none of the previous cases.
	 *  
	 *  3 - \f$\gamma_4 = \gamma_6 = 1\f$, but none of the previous cases.
	 * 
	 *  4 - None of the previous cases (the most general)
	 *
	 *  5 - special case of comparison to slow-rotation approximation:
	 *      this is identical to typeos=0, but using a modified
	 *      EOS-inversion method, namely we don't switch to a 1-fluid
	 *      EOS in 1-fluid regions.
	 * 
	 **/
	int typeos ; 

	/** Parameters needed for some inversions of the EOS.
	 *  In particular, it is used for type 4 EOS:
	 *  contains the relaxation parameter needed in the iteration
	 */
	double relax ;

	double precis ; ///< contains the precision required in zerosec_b
	
	///contains the precision required in the relaxation nbar_ent_p
	double ecart ; 

	
    // Constructors - Destructor
    // -------------------------
    public:
    
	/** Standard constructor.
	 *
	 *  The adiabatic indexes \f$\gamma_1\f$ and \f$\gamma_2\f$ are set to 2.
	 *  All other adiabatic indexes \f$\gamma_i\f$, \f$i\in 3\dots 6\f$ are
	 *  set to 1.
	 *  The individual particle masses \f$m_1\f$ and \f$m_2\f$ are set to the 
	 *  mean baryon mass \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$. 
	 *  The inversion parameters are set to their default values
	 *  (see hereafter the consrtuctor with all parameters).
	 *  
	 *  @param kappa1  pressure coefficient \f$\kappa_1\f$  
	 *  @param kappa2  pressure coefficient \f$\kappa_2\f$  
	 *  @param kappa3  pressure coefficient \f$\kappa_3\f$  
	 *  @param beta coefficient in the entrainment term \f$\beta\f$  
	 *		(cf. Eq.~(\ref{eeosbfpolye}))
	 *		[unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *		\f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ 
	 */
	Eos_bf_poly(double kappa1, double kappa2, double kappa3, double beta) ;

	/** Standard constructor with all parameters. 
	 * 
	 *  @param gamma1  adiabatic index \f$\gamma_1\f$ 
	 *  @param gamma2  adiabatic index \f$\gamma_2\f$ 
	 *  @param gamma3  adiabatic index \f$\gamma_3\f$ 
	 *  @param gamma4  adiabatic index \f$\gamma_4\f$ 
	 *  @param gamma5  adiabatic index \f$\gamma_5\f$ 
	 *  @param gamma6  adiabatic index \f$\gamma_6\f$ 
	 *				(cf. Eq.~(\ref{eeosbfpolye}))
	 *  @param kappa1  pressure coefficient \f$\kappa_1\f$  
	 *  @param kappa2  pressure coefficient \f$\kappa_2\f$  
	 *  @param kappa3  pressure coefficient \f$\kappa_3\f$  
	 *  @param beta coefficient in the entrainment term \f$\beta\f$  
	 *		(cf. Eq.~(\ref{eeosbfpolye}))
	 *		[unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *		\f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ 
	 *  @param mass1  individual particule mass \f$m_1\f$ (neutrons)  
	 *  @param mass2  individual particule mass \f$m_2\f$ (protons)  
	 *		[unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$]
	 *  @param relax relaxation parameter (see \c par_inv)
	 *  @param precis precision parameter for zerosec_b 
	 *                (see \c par_inv)
	 *  @param relax precision parameter for relaxation 
	 *               procedure (see \c par_inv)
	 *		 
	 */
	Eos_bf_poly(double gamma1, double gamma2, double gamma3,
		    double gamma4, double gamma5, double gamma6,
		    double kappa1, double kappa2, double kappa3,
		    double beta, double mass1=1, double mass2=1, 
		    double relax=0.5, double precis = 1.e-9,
		    double ecart = 1.e-8) ;	

	Eos_bf_poly(const Eos_bf_poly& ) ;	///< Copy constructor	
	
    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 * \c Eos_bifluid::eos_from_file(FILE*) . 
	 */
	Eos_bf_poly(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Eos_bifluid::eos_from_file(const char*) . 
	 *
	 */
	Eos_bf_poly (const char *fname) ; 
	
	/// The construction functions from a file
	friend Eos_bifluid* Eos_bifluid::eos_from_file(FILE* ) ; 
	friend Eos_bifluid* Eos_bifluid::eos_from_file(const char *fname ) ; 

        public:
	virtual ~Eos_bf_poly() ;			///< Destructor

    // Assignment
    // ----------
	/// Assignment to another \c Eos_bf_poly 
	void operator=(const Eos_bf_poly& ) ;


    // Miscellaneous
    // -------------
	public : 
	/// Comparison operator (egality)
	virtual bool operator==(const Eos_bifluid& ) const ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos_bifluid& ) const ; 
    
	/** Returns a number to identify the sub-classe of \c Eos_bifluid  
	 *  the object belongs to. 
	 */
	virtual int identify() const ; 

	/// Returns the adiabatic index \f$\gamma_1\f$ 
	double get_gam1() const {return gam1 ;};

	/// Returns the adiabatic index \f$\gamma_2\f$ 
	double get_gam2() const {return gam2 ;};

	/// Returns the adiabatic index \f$\gamma_3\f$ 
	double get_gam3() const {return gam3 ;};
	
	/// Returns the adiabatic index \f$\gamma_4\f$ 
	double get_gam4() const {return gam4 ;};
	
	/// Returns the adiabatic index \f$\gamma_5\f$ 
	double get_gam5() const {return gam5 ;};
	
	/// Returns the adiabatic index \f$\gamma_6\f$ 
	double get_gam6() const {return gam6 ;};
	
	/** Returns the pressure coefficient \f$\kappa_1\f$  
	 *  [unit: \f$\rho_{\rm nuc} c^2 \f$], where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$.
	 */
	double get_kap1() const {return kap1 ;}; 

	/** Returns the pressure coefficient \f$\kappa_2\f$  
	 *  [unit: \f$\rho_{\rm nuc} c^2 \f$], where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$.
	 */
	double get_kap2() const {return kap2 ;}; 

	/** Returns the pressure coefficient \f$\kappa_3\f$  
	 *  [unit: \f$\rho_{\rm nuc} c^2 \f$], where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$.
	 */
	double get_kap3() const {return kap3 ;}; 

	/** Returns the coefficient \f$\beta\f$  
	 *  [unit: \f$\rho_{\rm nuc} c^2 \f$], where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$.
	 */
	double get_beta() const {return beta ;}; 

	// Returns (sub)type of bifluid-eos (member \c typeos})
	int get_typeos() const {return typeos;};

    protected:
	/** Computes the auxiliary quantities \c gam1m1 , \c gam2m1 
	 *  and \c gam3m1 
	 */
	void set_auxiliary() ; 
    
	/// Determines the type of the analytical EOS (see \c typeos )
	void determine_type() ;

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	///< Save in a file

    protected: 
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>


    // Computational functions
    // -----------------------

    public: 

 	/** Computes both baryon densities from the log-enthalpies
	 * 
	 *  @param ent1 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_1\f$ 
	 *  @param ent2 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_2\f$ 
	 *  @param delta2 [input,  unit: \f$c^2\f$] relative velocity \f$\Delta^2\f$ 
	 * 
	 *  @param nbar1 [output] baryonic density of the first fluid
	 *  @param nbar2 [output] baryonic density of the second fluid
	 *  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 * 
	 */
    	virtual bool nbar_ent_p(const double ent1, const double ent2, 
				const double delta2, double& nbar1, 
				double& nbar2) const ; 
       
	/** Computes baryon density out of the log-enthalpy asuming
	 *  that only fluid 1 is present.
	 *  @param ent1 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_1\f$ 
	 *  @return nbar1 baryonic density of the first fluid
	 */
	virtual double nbar_ent_p1(const double ent1) const  ;

	/** Computes baryon density out of the log-enthalpy assuming
	 *  that only fluid 2 is present.
	 *  @param ent2 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_1\f$ 
	 *  @return nbar1 baryonic density of the first fluid
	 */
	virtual double nbar_ent_p2(const double ent2) const  ;

 	/** Computes the total energy density from the baryonic densities
	 *  and the relative velocity. 
	 * 
	 *  @param nbar1 [input] baryonic density of the first fluid
	 *  @param nbar2 [input] baryonic density of the second fluid
	 *  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *  @param delta2 [input,  unit: \f$c^2\f$] relative velocity \f$\Delta^2\f$ 
	 * 
	 *  @return energy density \f$\cal E\f$ [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_nbar_p(const double nbar1, const double nbar2, 
				   const double delta2) const  ; 
       
 	/** Computes the pressure from the baryonic densities
	 *  and the relative velocity.
	 * 
	 *  @param nbar1 [input] baryonic density of the first fluid
	 *  @param nbar2 [input] baryonic density of the second fluid
	 *  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *  @param delta2 [input,  unit: \f$c^2\f$] relative velocity \f$\Delta^2\f$ 
	 * 
	 *  @return pressure \e p  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_nbar_p(const double nbar1, const double nbar2, 
				    const double delta2) const ; 
     // Conversion functions
     // ---------------------

	/** Makes a translation from \c Eos_bifluid  to \c Eos . 
	 *
	 *  This is only useful for the construction of a 
	 *  \c Et_rot_bifluid 
	 *  star and ought not to be used in other situations.
	 */
	virtual Eos* trans2Eos() const ;
       
	/** Computes the derivative of the energy with respect to
	 * (baryonic density 1)\f$^2\f$.
	 *
	 *  @param n1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *           relative Lorentz factor\f$\times\f$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative \f$K^{11}=2\frac{\partial{\cal{E}}}{\partial{n_1^2}}\f$ 
	 */
	virtual double get_K11(const double n1, const double n2, const
			       double delta2)  const  ;

	/** Computes the derivative of the energy with respect to 
	 *  \f$x^2=n_1n_2\Gamma_\Delta\f$.
	 *
	 *  @param n1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *           relative Lorentz factor\f$\times\f$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative \f$K^{12}=\frac{\partial {\cal E}}{\partial (n_1n_2\Gamma_\Delta)}\f$ 
	 */
	virtual double get_K12(const double n1, const double n2,const
			       double delta2) const ;

	/** Computes the derivative of the energy/(baryonic density 2)\f$^2\f$.
	 *
	 *  @param n1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *           relative Lorentz factor\f$\times\f$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative \f$K^{22} = 2\frac{\partial {\cal E}}{\partial n_2^2}\f$ 
	 */
	virtual double get_K22(const double n1, const double n2, const
			       double delta2) const ;

};

		    //------------------------------------//
		    //	      class Eos_bf_poly_newt	  //
		    //------------------------------------//

/**
 * Analytic equation of state for two fluids (Newtonian case).
 * 
 * This equation of state (EOS) corresponds to two types of non-relativistic
 * particles of rest mass is \f$m_1\f$ and \f$m_2\f$,  whose total energy density 
 * \f$\cal{E}\f$ is related to their numerical densities \f$n_1\f$, \f$n_2\f$ and 
 * relative velocity \f$\Delta^2\f$
 * \f[
 * \Delta = \left( \vec{v}_n - \vec{v}_p \right)^2
 * \label{e:defdeltan}
 * \f]
 * by
 * \f[ \label{eeosbfnewte}
 *   {\cal{E}} = \frac{1}{2}\kappa_1 n_1^{\gamma_1} 
 *          \  + \frac{1}{2}\kappa_2 n_2^{\gamma_2} 
 *          \  + \kappa_3 n_1^{\gamma_3} n_2^{\gamma_4} 
 *          \  + \Delta^2 \beta n_1^{\gamma_5} n_2^{\gamma_6}\ .  
 * \f]
 * The non-relativistic chemical potentials are then
 * \f[ 
 * \mu_1 := {{\rm d}{\cal{E}} \over {\rm d}n_1} = \frac{1}{2}\gamma_1\kappa_1
 *         n_1^{\gamma_1-1} + \gamma_3 \kappa_3 
 *         n_1^{\gamma_3-1} n_2^{\gamma_4} + \Delta^2 \gamma_5 \beta 
 *         n_1^{\gamma_5-1} n_2^{\gamma_6}\ , 
 * \f]
 * \f[
 * \mu_2 := {{\rm d}{\cal{E}} \over {\rm d}n_2} = \frac{1}{2}\gamma_2\kappa_2
 *         n_2^{\gamma_2-1} + \gamma_4 \kappa_3 
 *         n_1^{\gamma_3} n_2^{\gamma_4-1} + \Delta^2 \gamma_6 \beta 
 *         n_1^{\gamma_5} n_2^{\gamma_6-1} \ .
 * \f]
 * The pressure is given by the (zero-temperature) First Law of Thermodynamics:
 * \f$p = \mu_1 n_1 + \mu_2 n_2 - {\cal E}\f$, so that
 * \f[ 
 *   p = \frac{1}{2} (\gamma_1 -1)\kappa_1 n_1^{\gamma_1} +
 *  \frac{1}{2}(\gamma_2-1)\kappa_2 n_2^{\gamma_2} + (\gamma_3 +\gamma_4
 *  -1)\kappa_3 n_1^{\gamma_3}n_2^{\gamma_4} + \Delta^2 \beta \left(
 *  (\gamma_5 + \gamma_6 - 1) n_1^{\gamma_5} n_2^{\gamma_6} \right) \ .  
 * \f]
 * The
 * specific enthalpies are related to the chemical potentials by
 * \f[
 * h_i = \frac{\mu_i}{m_i}
 * \f]
 *
 * From this system, the particle densities are obtained in term of 
 * the enthalpies. (The system is a linear one
 * if \f$\gamma_1 = \gamma_2 = 2\f$ and \f$\gamma_3 = \gamma_4 = \gamma_5 =
 * \gamma_6 = 1\f$). \ingroup (eos)
 *
 * The energy density \f$\cal E\f$and pressure \e p  can then be obtained.
 *
 */
class Eos_bf_poly_newt : public Eos_bf_poly {

    // Data :
    // -----

    // no new data with respect to Eos_bf_poly	

    // Constructors - Destructor
    // -------------------------
    public:
    
	/** Standard constructor.
	 *
	 *  The adiabatic indexes \f$\gamma_1\f$ and \f$\gamma_2\f$ are set to 2.
	 *  All other adiabatic indexes \f$\gamma_i\f$, \f$i\in 3\dots 6\f$ are
	 *  set to 1.
	 *  The individual particle masses \f$m_1\f$ and \f$m_2\f$ are set to the 
	 *  mean baryon mass \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$. 
	 *  The inversion parameters are set to their default values
	 *  (see hereafter the consrtuctor with all parameters).
	 *  
	 *  @param kappa1  pressure coefficient \f$\kappa_1\f$  
	 *  @param kappa2  pressure coefficient \f$\kappa_2\f$  
	 *  @param kappa3  pressure coefficient \f$\kappa_3\f$  
	 *  @param beta coefficient in the entrainment term \f$\beta\f$  
	 *		(cf. Eq.~(\ref{eeosbfpolye}))
	 *		[unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *		\f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ 
	 */
	Eos_bf_poly_newt(double kappa1, double kappa2, double kappa3,
		    double beta) ;

	/** Standard constructor with all parameters. 
	 * 
	 *  @param gamma1  adiabatic index \f$\gamma_1\f$ 
	 *  @param gamma2  adiabatic index \f$\gamma_2\f$ 
	 *  @param gamma3  adiabatic index \f$\gamma_3\f$ 
	 *  @param gamma4  adiabatic index \f$\gamma_4\f$ 
	 *  @param gamma5  adiabatic index \f$\gamma_5\f$ 
	 *  @param gamma6  adiabatic index \f$\gamma_6\f$ 
	 *				(cf. Eq.~(\ref{eeosbfpolye}))
	 *  @param kappa1  pressure coefficient \f$\kappa_1\f$  
	 *  @param kappa2  pressure coefficient \f$\kappa_2\f$  
	 *  @param kappa3  pressure coefficient \f$\kappa_3\f$  
	 *  @param beta coefficient in the entrainment term \f$\beta\f$  
	 *		(cf. Eq.~(\ref{eeosbfpolye}))
	 *		[unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *		\f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ 
	 *  @param mass1  individual particule mass \f$m_1\f$ (neutrons)  
	 *  @param mass2  individual particule mass \f$m_2\f$ (protons)  
	 *  @param relax relaxation parameter (see \c par_inv)
	 *  @param precis precision parameter for zerosec_b 
	 *                (see \c par_inv)
	 *  @param relax precision parameter for relaxation 
	 *               procedure (see \c par_inv)
	 *		 
	 *		[unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$]
	 */
	Eos_bf_poly_newt(double gamma1, double gamma2, double gamma3,
			 double gamma4, double gamma5, double gamma6,
			 double kappa1, double kappa2, double kappa3,
			 double beta, double mass1, double mass2, 
			 double relax=0.5, double precis = 1.e-9,
			 double ecart = 1.e-8) ;	


	Eos_bf_poly_newt(const Eos_bf_poly_newt& ) ;	///< Copy constructor	
	
    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 * \c Eos_bifluid::eos_from_file(FILE*) . 
	 */
	Eos_bf_poly_newt(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Eos_bifluid::eos_from_file(const char* ) . 
	 */
	Eos_bf_poly_newt(const char *fname ) ; 
	
	/// The construction functions from a file
	friend Eos_bifluid* Eos_bifluid::eos_from_file(FILE* ) ; 
	friend Eos_bifluid* Eos_bifluid::eos_from_file(const char *fname) ; 

    public:
	virtual ~Eos_bf_poly_newt() ;			///< Destructor

    // Assignment
    // ----------
	/// Assignment to another \c Eos_bf_poly_newt 
	void operator=(const Eos_bf_poly_newt& ) ;


    // Miscellaneous
    // -------------

    public : 
	/// Comparison operator (egality)
	virtual bool operator==(const Eos_bifluid& ) const ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos_bifluid& ) const ; 
    
	/** Returns a number to identify the sub-classe of \c Eos_bifluid  
	 *  the object belongs to. 
	 */
	virtual int identify() const ; 

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	///< Save in a file

    protected: 
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>


    // Computational functions
    // -----------------------

    public: 

 	/** Computes both baryon densities from the log-enthalpies
	 * 
	 *  @param ent1 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_1\f$ 
	 *  @param ent2 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_2\f$ 
	 *  @param delta2 [input,  unit: \f$c^2\f$] relative velocity \f$\Delta^2\f$ 
	 * 
	 *  @param nbar1 [output] baryonic density of the first fluid
	 *  @param nbar2 [output] baryonic density of the second fluid
	 *  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 * 
	 */
    	virtual bool nbar_ent_p(const double ent1, const double ent2, 
				const double delta2, double& nbar1, 
				double& nbar2) const ; 
       
	/** Computes baryon density out of the log-enthalpy asuming
	 *  that only fluid 1 is present (virtual function implemented 
	 *  in the derived classes).
	 *  @param ent1 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_1\f$ 
	 *  @return nbar1 baryonic density of the first fluid
	 */
	virtual double nbar_ent_p1(const double ent1) const  ;

	/** Computes baryon density out of the log-enthalpy assuming
	 *  that only fluid 2 is present.
	 *  @param ent2 [input,  unit: \f$c^2\f$] log-enthalpy \f$H_1\f$ 
	 *  @return nbar1 baryonic density of the first fluid
	 */
	virtual double nbar_ent_p2(const double ent2) const  ;

 	/** Computes the total energy density from the baryonic densities
	 *  and the relative velocity. 
	 * 
	 *  @param nbar1 [input] baryonic density of the first fluid
	 *  @param nbar2 [input] baryonic density of the second fluid
	 *  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *  @param delta2 [input,  unit: \f$c^2\f$] relative velocity \f$\Delta^2\f$ 
	 * 
	 *  @return energy density \f$\cal E\f$ [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_nbar_p(const double nbar1, const double nbar2, 
				   const double delta2) const  ; 
       
 	/** Computes the pressure from the baryonic densities
	 *  and the relative velocity.
	 * 
	 *  @param nbar1 [input] baryonic density of the first fluid
	 *  @param nbar2 [input] baryonic density of the second fluid
	 *  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *  @param delta2 [input,  unit: \f$c^2\f$] relative velocity \f$\Delta^2\f$ 
	 * 
	 *  @return pressure \e p  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_nbar_p(const double nbar1, const double nbar2, 
				    const double delta2) const ; 
     // Conversion functions
     // ---------------------

	/** Makes a translation from \c Eos_bifluid  to \c Eos . 
	 *
	 *  This is only useful for the construction of a 
	 *  \c Et_rot_bifluid 
	 *  star and ought not to be used in other situations.
	 */
	virtual Eos* trans2Eos() const ;
       
	/** Computes the derivative of the energy with respect to
	 * (baryonic density 1)\f$^2\f$.
	 *
	 *  @param n1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *           relative Lorentz factor\f$\times\f$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative \f$K^{11}=2\frac{\partial{\cal{E}}}{\partial{n_1^2}}\f$ 
	 */
	virtual double get_K11(const double n1, const double n2, const
			       double delta2)  const  ;

	/** Computes the derivative of the energy with respect to 
	 *  \f$x^2=n_1n_2\Gamma_\Delta\f$.
	 *
	 *  @param n1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *           relative Lorentz factor\f$\times\f$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative \f$K^{12}=\frac{\partial {\cal E}}{\partial (n_1n_2\Gamma_\Delta)}\f$ 
	 */
	virtual double get_K12(const double n1, const double n2,const
			       double delta2) const ;

	/** Computes the derivative of the energy/(baryonic density 2)\f$^2\f$.
	 *
	 *  @param n1 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit \f$n_{\rm nuc}^2c^2\f$]
	 *           relative Lorentz factor\f$\times\f$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative \f$K^{22} = 2\frac{\partial {\cal E}}{\partial n_2^2}\f$ 
	 */
	virtual double get_K22(const double n1, const double n2, const
			       double delta2) const ;

};

		

#endif
