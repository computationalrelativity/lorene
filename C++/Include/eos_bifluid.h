/*
 *  Definition of Lorene classes Eos_bifluid
 *				 Eos_bf_poly
 *
 */

/*
 *   Copyright (c) 2001 Jerome Novak
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
class ostream ; 
class ifstream ; 

// Headers C
#include <stdio.h>

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

/**
 * 2-fluids equation of state base class.
 *
 * Fluid 1 is supposed to correspond to neutrons, whereas fluid 2
 * corresponds to e.g. protons. Neutron 4-velocity is $u^\alpha_{\rm n}$
 * and proton one is $u^\alpha_{\rm p}$
 *
 * Therefore, the EOS is defined by giving two log-enthalpies AND a 
 * relative velocity as inputs. The output are then: two baryonic
 * densities, the total energy density and pressure
 * The enthalpies $f_1$ and $f_2$ are obtained through the formula
 * \begin{equation}
 * {\rm d}{\cal E}=f^1{\rm d}n_1+f^2{\rm d}n_2+\alpha{\rm d}\Delta^2
 * \label{eeosbfdefent}
 * \end{equation}
 * see Comer, Novak \& Prix. Log-enthalpies are then defined as
 * $H_1 = \ln\left( \frac{f^1}{m_1c^2} \right)$, where $m_1$ is the
 * mass of a particle of the first fluid. (same for $H_2$)
 * The relative velocity $\Delta^2$
 * is defined as in Comer, Novak \& Prix. It can be seen as
 * the neutron velocity seen in the frame of protons (or vice-versa): 
 * $\Gamma_\Delta = -g_{\alpha\beta} u^\alpha_{\rm n} u^\beta_{\rm p}$
 *
 * @version #$Id$#
 */
class Eos_bifluid {

    // Data :
    // -----

    protected: 
	char name[100] ;	    /// EOS name

    // Constructors - Destructor
    // -------------------------
    protected:
	Eos_bifluid() ;			/// Standard constructor

	/// Standard constructor with name
	explicit Eos_bifluid(const char* name_i) ; 

	Eos_bifluid(const Eos_bifluid& ) ;	/// Copy constructor	

    protected:
	/** Constructor from a binary file (created by the function 
	 *  {\tt sauve(FILE* )}). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  {\tt Eos\_bifluid::eos\_from\_file(FILE* )}. 
	 */
	Eos_bifluid(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  {\tt Eos\_bifluid::eos\_from\_file(ifstream\&)}. 
	 */
	Eos_bifluid(ifstream& ) ; 
	
	
    public:
	virtual ~Eos_bifluid() ;			/// Destructor


    // Name manipulation
    // -----------------
    public:
	const char* get_name() const ;	/// Returns the EOS name

	/// Sets the EOS name
	void set_name(const char* name_i) ; 

    // Miscellaneous
    // -------------
    public:
	/** Construction of an EOS from a binary file.
	 *  The file must have been created by the function 
	 *  {\tt sauve(FILE* )}.
	 */
	static Eos_bifluid* eos_from_file(FILE* ) ; 
	
	/** Construction of an EOS from a formatted file.
	 * 
	 *  The first line of the file must start by the EOS number, according 
	 *  to the following conventions: \\
	 *	1 = analytic EOS (class {\tt Eos\_bf\_poly}). \\
	 *      Up to now (21.06.2001) only analytic EOS is implemented
	 */
	static Eos_bifluid* eos_from_file(ifstream& ) ; 
	
	/// Comparison operator (egality)
	virtual bool operator==(const Eos_bifluid& ) const = 0 ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos_bifluid& ) const = 0 ; 
    
	/** Returns a number to identify the sub-classe of {\tt Eos\_bifluid} 
	 *  the object belongs to. 
	 */
	virtual int identify() const = 0 ; 

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	/// Save in a file

	/// Display
	friend ostream& operator<<(ostream& , const Eos_bifluid& ) ;	

    protected: 
	virtual ostream& operator>>(ostream &) const = 0 ;    /// Operator >>


    // Computational functions
    // -----------------------
    public:
	/**  General computational method for {\tt Cmp}'s, it computes
	 *   both baryon densities, energy and pressure profiles.
	 * 
	 *  @param ent1 [input] the first log-enthalpy field $H_1$.  
	 *  @param ent2 [input] the second log-enthalpy field $H_2$.
	 *  @param delta2 [input] the relative velocity field $\Delta^2 $
	 *  @param nbar1 [output] baryonic density of the first fluid
	 *  @param nbar2 [output] baryonic density of the second fluid
	 *  [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *  @param ener [output] total energy density $\cal E$ 
	 *                             of both fluids together
	 *  @param press [output] pressure $p$ of both fluids together
	 *  @param nzet  [input] number of domains where {\tt resu} is to be
	 *	computed. 
	 *  @param l_min [input] index of the innermost domain is which 
	 *      {\tt resu} is to be computed [default value: 0]; 
	 *      {\tt resu} is computed only in domains whose indices are 
	 *      in {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero. 
	 */
	void calcule_tout(const Cmp& ent1, const Cmp& ent2, const Cmp& delta2,
			  Cmp& nbar1, Cmp& nbar2, Cmp& ener, Cmp& press,
			  int nzet, int l_min = 0) const ; 

 	/** Computes both baryon densities from the log-enthalpies 
	 *  (virtual function implemented in the derived classes). 
	 * 
	 *  @param ent1 [input,  unit: $c^2$] log-enthalpy $H_1$ 
	 *  @param ent2 [input,  unit: $c^2$] log-enthalpy $H_2$ 
	 *  @param delta2 [input,  unit: $c^2$] relative velocity $\Delta^2$ 
	 *  @param tronc [input] if true, the resulting densities cannot be
	 *      negative.
	 * 
	 *  @param nbar1 [output] baryonic density of the first fluid
	 *  @param nbar2 [output] baryonic density of the second fluid
	 *  [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 * 
	 */
    	virtual void nbar_ent_p(const double ent1, const double ent2, 
				const double delta2, double& nbar1, 
				double& nbar2, bool tronc = true) const = 0 ; 

	/** Computes both baryon density fields from the log-enthalpy fields
	 *  and the relative velocity.
	 * 
	 *  @param ent1 [input,  unit: $c^2$] log-enthalpy $H_1$ 
	 *  @param ent2 [input,  unit: $c^2$] log-enthalpy $H_2$ 
	 *  @param delta2 [input,  unit: $c^2$] relative velocity $\Delta^2$ 
	 *
	 *  @param nbar1 [output] baryonic density of the first fluid
	 *  @param nbar2 [output] baryonic density of the second fluid
	 *  [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *  @param nzet  number of domains where the baryon density is to be
	 *	computed. 
	 *  @param l_min  index of the innermost domain is which the baryon 
	 *	density is
	 *	to be computed [default value: 0]; the baryon density is 
	 *	computed only in domains whose indices are in 
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero. 
	 *  @param tronc  if true, the resulting densities cannot be
	 *      negative.
	 * 
	 */
    	void nbar_ent(const Cmp& ent1, const Cmp& ent2, const Cmp& delta2, 
		      Cmp& nbar1, Cmp& nbar2, int nzet, int l_min = 0,
		      bool tronc = true) const  ; 
    
 	/** Computes the total energy density from the baryonic densities
	 *  and the relative velocity.
	 *  (virtual function implemented in the derived classes). 
	 * 
	 *  @param nbar1 [input] baryonic density of the first fluid
	 *  @param nbar2 [input] baryonic density of the second fluid
	 *  [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *  @param delta2 [input,  unit: $c^2$] relative velocity $\Delta^2$ 
	 * 
	 *  @return energy density $\cal E$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double ener_nbar_p(const double nbar1, const double nbar2, 
				   const double delta2) const = 0 ; 
       
	/** Computes the total energy density from the log-enthalpy fields
	 *  and the relative velocity.
	 * 
	 *  @param ent1 [input,  unit: $c^2$] log-enthalpy $H_1$ 
	 *  @param ent2 [input,  unit: $c^2$] log-enthalpy $H_2$ 
	 *  @param delta2 [input,  unit: $c^2$] relative velocity $\Delta^2$ 
	 *  @param nzet  number of domains where the energy density is to be
	 *	computed. 
	 *  @param l_min  index of the innermost domain is which the energy 
	 *	density is
	 *	to be computed [default value: 0]; the energy density is 
	 *	computed only in domains whose indices are in 
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero. 
	 * 
	 *  @return  energy density field [unit: $\rho_{\rm nuc} c^2$], 
	 *      where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	Cmp ener_ent(const Cmp& ent1, const Cmp& ent2, const Cmp& delta2, 
		      int nzet, int l_min = 0) const ; 
    
 	/** Computes the pressure from the baryonic densities
	 *  and the relative velocity.
	 *  (virtual function implemented in the derived classes). 
	 * 
	 *  @param nbar1 [input] baryonic density of the first fluid
	 *  @param nbar2 [input] baryonic density of the second fluid
	 *  [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *  @param delta2 [input,  unit: $c^2$] relative velocity $\Delta^2$ 
	 * 
	 *  @return pressure $p$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double press_nbar_p(const double nbar1, const double nbar2, 
				    const double delta2) const = 0 ; 
       
	/** Computes the pressure from the log-enthalpy fields
	 *  and the relative velocity.
	 * 
	 *  @param ent1 [input,  unit: $c^2$] log-enthalpy $H_1$ 
	 *  @param ent2 [input,  unit: $c^2$] log-enthalpy $H_2$ 
	 *  @param delta2 [input,  unit: $c^2$] relative velocity $\Delta^2$ 
	 *  @param nzet  number of domains where the pressure is to be
	 *	computed. 
	 *  @param l_min  index of the innermost domain is which the pressure
	 *    is to be computed [default value: 0]; the pressure is computed 
	 *      only in domains whose indices are in 
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero. 
	 * 
	 *  @return pressure field [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 * 
	 */
    	Cmp press_ent(const Cmp& ent1, const Cmp& ent2, const Cmp& delta2, 
		      int nzet, int l_min = 0) const ; 

	/** Computes the derivative of the energy with respect to
	 * (baryonic density 1)$^2$.
	 *  (virtual function implemented in the derived classes). 
	 *
	 *  @param n1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit $n_{\rm nuc}^2c^2$]
	 *           relative Lorentz factor$\times$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative $K^{11}=2\frac{\partial{\cal{E}}}{\partial{n_1^2}}$ 
	 */
	virtual double get_K11(const double n1, const double n2, const
			       double x)  const = 0 ;

	/** Computes the derivative of the energy with respect to 
	 *  $x^2=n_1n_2\Gamma_\Delta$.
	 *  (virtual function implemented in the derived classes). 
	 *
	 *  @param n1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit $n_{\rm nuc}^2c^2$]
	 *           relative Lorentz factor$\times$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative $K^{12}=\frac{\partial {\cal E}}{\partial (n_1n_2\Gamma_\Delta)}$ 
	 */
	virtual double get_K12(const double n1, const double n2,const
			       double x) const = 0 ;

	/** Computes the derivative of the energy/(baryonic density 2)$^2$.
	 *  (virtual function implemented in the derived classes). 
	 *
	 *  @param n1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit $n_{\rm nuc}^2c^2$]
	 *           relative Lorentz factor$\times$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative $K^{22} = 2\frac{\partial {\cal E}}{\partial n_2^2}$ 
	 */
	virtual double get_K22(const double n1, const double n2, const
			       double x) const = 0 ;

	/** Computes the derivatives of the energy/(baryonic density 1)$^2$.
	 *
	 *  @param nbar1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density field of fluid 1 at which 
	 *           the derivatives are to be computed
	 *  @param nbar2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density field of fluid 2 at which 
	 *           the derivatives are to be computed
	 *  @param x2 [input, unit $n_{\rm nuc}^2c^2$]
	 *             relative velocity$\times$both densities at which 
	 *           the derivative is to be computed	
	 *  @param nzet  number of domains where the derivatives are to be
	 *	computed. 
	 *  @param l_min  index of the innermost domain is which the
	 *	derivatives are
	 *	to be computed [default value: 0]; the derivatives are
	 *	computed only in domains whose indices are in 
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero. 
	 * 
	 *  @return derivative $K^{11}$ field (see {\tt get\_K11}) 
	 */
    	Cmp get_Knn(const Cmp& nbar1, const Cmp& nbar2, const Cmp& x2,
		    int nzet, int l_min = 0) const  ; 

	/** Computes the derivatives of the energy/(baryonic density 2)$^2$.
	 *
	 *  @param nbar1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density field of fluid 1 at which 
	 *           the derivatives are to be computed
	 *  @param nbar2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density field of fluid 2 at which 
	 *           the derivatives are to be computed
	 *  @param x2 [input, unit $n_{\rm nuc}^2c^2$]
	 *             relative velocity$\times$both densities at which 
	 *           the derivative is to be computed	
	 *  @param nzet  number of domains where the derivatives are to be
	 *	computed. 
	 *  @param l_min  index of the innermost domain is which the
	 *	derivatives are
	 *	to be computed [default value: 0]; the derivatives are
	 *	computed only in domains whose indices are in 
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero. 
	 * 
	 *  @return derivative $K^{22}$ field (see {\tt get\_K12}) 
	 */
    	Cmp get_Kpp(const Cmp& nbar1, const Cmp& nbar2, const Cmp&
		    x2, int nzet, int l_min = 0) const ; 

	/** Computes the derivatives of the energy with respect to
	 *  $x^2=n_1n_2\Gamma_\Delta^2$.
	 *
	 *  @param nbar1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density field of fluid 1 at which 
	 *           the derivatives are to be computed
	 *  @param nbar2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density field of fluid 2 at which 
	 *           the derivatives are to be computed
	 *  @param x2 [input, unit $n_{\rm nuc}^2c^2$]
	 *             relative velocity$\times$both densities at which 
	 *           the derivative is to be computed	
	 *  @param nzet  number of domains where the derivatives are to be
	 *	computed. 
	 *  @param l_min  index of the innermost domain is which the
	 *	derivatives are
	 *	to be computed [default value: 0]; the derivatives are
	 *	computed only in domains whose indices are in 
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero. 
	 * 
	 *  @return derivative $K^{12}$ field (see {\tt get\_K12}) 
	 */
     	Cmp get_Knp(const Cmp& nbar1, const Cmp& nbar2, const Cmp& x2,
		    int nzet, int l_min = 0) const ; 

	/**  General computational method for {\tt Cmp}'s ($K^{ij}$'s).
	 * 
	 *  @param nbar1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density field of fluid 1 at which 
	 *           the derivatives are to be computed.
	 *  @param nbar2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density field of fluid 2 at which 
	 *           the derivatives are to be computed
	 *  @param x2 [input, unit $n_{\rm nuc}^2c^2$]
	 *             relative velocity$\times$both densities at which 
	 *           the derivative is to be computed	
	 *  @param nzet  [input] number of domains where {\tt resu} is to be
	 *	computed. 
	 *  @param l_min [input] index of the innermost domain is which {\tt resu}
	 *	is to be computed [default value: 0]; {\tt resu} is 
	 *	computed only in domains whose indices are in 
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero. 
	 *  @param fait [input] pointer on the member function of class 
	 *	  {\tt Eos\_bifluid} which performs the pointwise calculation.
	 *  @param resu [output] result of the computation. 
	 */
	void calcule(const Cmp& nbar1, const Cmp& nbar2, const Cmp&
		     x2, int nzet, int l_min, double
		     (Eos_bifluid::*fait)(double, double, double) const, 
		     Cmp& resu)
	  const ; 
 	
   // Conversion functions
    // ---------------------

	/** Makes a translation from {\tt Eos\_bifluid} to {\tt Eos}. 
	 *  (virtual function implemented in the derived classes). 
	 *
	 *  This is only useful for the construction of a 
	 *  {\tt Et\_rot\_bifluid}
	 *  star and ought not to be used in other situations.
	 *  @param relat [input] Relativistic EOS or not. 
	 */
	virtual Eos* trans2Eos() const = 0 ;
    
};


		    //------------------------------------//
		    //	      class Eos_bf_poly		  //
		    //------------------------------------//

/**
 * Analytic equation of state for two fluids (relativistic case).
 * 
 * This equation of state (EOS) corresponds to two types of relativistic
 * particles of rest mass is $m_1$ and $m_2$,  whose total energy density 
 * $\cal{E}$ is related to their numerical densities $n_1$, $n_2$ and 
 * relative velocity 
 * \begin{equation}
 * \Gamma_\Delta = -g_{\alpha\beta} u^\alpha_{\rm n} u^\beta_{\rm p}
 * \label{e:defgamamdelta}
 * \end{equation}
 * ($u^\alpha_{\rm n}$ and $u^\alpha_{\rm p}$ being the 4-velocities of 
 * both fluids), by
 * \begin{equation} \label{eeosbfpolye}
 *   {\cal{E}} = \frac{1}{2}\kappa_1 n_1^{\gamma_1} + m_1 c^2 \, n_1 
 *          \  + \frac{1}{2}\kappa_2 n_2^{\gamma_2} + m_2 c^2 \, n_2 
 *          \  + \kappa_3 n_1^{\gamma_3} n_2^{\gamma_4} 
 *          \  + \Delta^2 \beta n_1^{\gamma_5} n_2^{\gamma_6}\ .  
 * \end{equation}
 * The relativistic (i.e. including rest mass energy) chemical potentials
 * are then
 * \begin{equation} 
 * \mu_1 := {{\rm d}{\cal{E}} \over {\rm d}n_1} = \frac{1}{2}\gamma_1\kappa_1
 *         n_1^{\gamma_1-1} + m_1 c^2 + \gamma_3 \kappa_3 
 *         n_1^{\gamma_3-1} n_2^{\gamma_4} + \Delta^2 \gamma_5 \beta
 *         n_1^{\gamma_5-1} n_2^{\gamma_6}\ , 
 * \end{equation}
 * \begin{equation}
 * \mu_2 := {{\rm d}{\cal{E}} \over {\rm d}n_2} = \frac{1}{2}\gamma_2\kappa_2
 *         n_2^{\gamma_2-1} + m_2 c^2 + \gamma_4 \kappa_3 
 *         n_1^{\gamma_3} n_2^{\gamma_4-1} + \Delta^2 \gamma_6 \beta
 *         n_1^{\gamma_5} n_2^{\gamma_6-1} \ .
 * \end{equation}
 * The pressure is given by the (zero-temperature) First Law of Thermodynamics:
 * $p = \mu_1 n_1 + \mu_2 n_2 - {\cal E}$, so that
 * \begin{equation} 
 *   p = \frac{1}{2} (\gamma_1 -1)\kappa_1 n_1^{\gamma_1} +
 *  \frac{1}{2}(\gamma_2-1)\kappa_2 n_2^{\gamma_2} + (\gamma_3 +\gamma_4
 *  -1)\kappa_3 n_1^{\gamma_3}n_2^{\gamma_4} + \Delta^2 \beta \left(
 *  (\gamma_5 + \gamma_6 - 1) n_1^{\gamma_5} n_2^{\gamma_6} \right) \ .  
 * \end{equation}
 * The log-enthalpies are defined as the logarithm of the ratio of the enthalpy
 * per particle (see Eq.~\ref{eeosbfdefent}) by the particle rest mass energy :
 * \begin{equation}\label{eeosbfentanal} 
 *   H_i := c^2 \ln \left( \frac{f_i}{m_ic^2} \right)   \ .  
 * \end{equation}
 * From this system, the particle densities are obtained in term of 
 * the log-enthalpies. (The system (\ref{eeosbfentanal}) is a linear one
 * if $\gamma_1 = \gamma_2 = 2$ and $\gamma_3 = \gamma_4 = \gamma_5 =
 * \gamma_6 = 1$).
 *
 * The energy density $\cal E$and pressure $p$ can then be obtained
 * as functions of baryonic densities.
 *
 * @version #$Id$#
 */
class Eos_bf_poly : public Eos_bifluid {

    // Data :
    // -----

    protected: 
	/// Adiabatic indexes $\gamma_1$, see Eq.~\ref{eeosbfpolye}
	double gam1 ;

	/// Adiabatic indexes $\gamma_2$, see  Eq.~\ref{eeosbfpolye}
	double gam2 ;

	/// Adiabatic indexes $\gamma_3$, see  Eq.~\ref{eeosbfpolye}
	double gam3 ;
	
	/// Adiabatic indexes $\gamma_4$, see Eq.~\ref{eeosbfpolye}
	double gam4 ;

	/// Adiabatic indexes $\gamma_5$, see  Eq.~\ref{eeosbfpolye}
	double gam5 ;

	/// Adiabatic indexes $\gamma_6$, see  Eq.~\ref{eeosbfpolye}
	double gam6 ;
	
	/** Pressure coefficient $\kappa_1$  , see Eq.~\ref{eeosbfpolye}
	 *  [unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$], where
	 *  $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ and
	 *  $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$. 
	 */
	double kap1 ; 

	/** Pressure coefficient $\kappa_2$  , see Eq.~\ref{eeosbfpolye}
	 *  [unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$], where
	 *  $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ and
	 *  $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$. 
	 */
	double kap2 ; 

	/** Pressure coefficient $\kappa_3$  , see Eq.~\ref{eeosbfpolye}
	 *  [unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$], where
	 *  $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ and
	 *  $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$. 
	 */
	double kap3 ; 
	
	/** Coefficient $\beta$  , see Eq.~\ref{eeosbfpolye}
	 *  [unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$], where
	 *  $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ and
	 *  $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$. 
	 */
	double beta ; 

	/** Individual particule mass $m_1$  
	 *  [unit: $m_B = 1.66\ 10^{-27} \ {\rm kg}$]. 
	 */
	double m_1 ; 

	/** Individual particule mass $m_2$  
	 *  [unit: $m_B = 1.66\ 10^{-27} \ {\rm kg}$]. 
	 */
	double m_2 ; 

	double gam1m1 ;	    /// $\gamma_1-1$
	double gam2m1 ;	    /// $\gamma_2-1$
	double gam34m1 ;    /// $\gamma_3+\gamma_4-1$
	double gam56m1 ;    /// $\gamma_5+\gamma_6-1$

 protected:
	/** The bi-fluid analytical EOS type:
	 * 
	 *  0 - $\gamma_1 = \gamma_2 = 2$ and 
	 *  $\gamma_3 = \gamma_4 = \gamma_5 = \gamma_6 = 1$. In this case, 
	 *  the EOS can be inverted analytically.
	 *
	 *  1 - $\gamma_3 = \gamma_4 = \gamma_5 = \gamma_6 = 1$, but
	 *  $\gamma_1 \not= 2$ or $\gamma_2 \not= 2$. 
	 *
	 *  2 - $\gamma_3 = \gamma_5 = 1$, but none of the previous cases.
	 *  
	 *  3 - $\gamma_4 = \gamma_6 = 1$, but none of the previous cases.
	 * 
	 *  4 - None of the previous cases (the most general)
	 **/
	int typeos ; 

	/** Parameters needed for some inversions of the EOS.
	 *  In particular, it is used for type 4 EOS:
	 *  contains the relaxation parameter needed in the iteration
	 */
	double relax ;

	double precis ; /// contains the precision required in zerosec\_b
	
	///contains the precision required in the relaxation nbar\_ent\_p
	double ecart ; 

	
    // Constructors - Destructor
    // -------------------------
    public:
    
	/** Standard constructor.
	 *
	 *  The adiabatic indexes $\gamma_1$ and $\gamma_2$ are set to 2.
	 *  All other adiabatic indexes $\gamma_i$, $i\in 3\dots 6$ are
	 *  set to 1.
	 *  The individual particle masses $m_1$ and $m_2$ are set to the 
	 *  mean baryon mass $m_B = 1.66\ 10^{-27} \ {\rm kg}$. 
	 *  The inversion parameters are set to their default values
	 *  (see hereafter the consrtuctor with all parameters).
	 *  
	 *  @param kappa1  pressure coefficient $\kappa_1$  
	 *  @param kappa2  pressure coefficient $\kappa_2$  
	 *  @param kappa3  pressure coefficient $\kappa_3$  
	 *  @param beta coefficient in the entrainment term $\beta$  
	 *		(cf. Eq.~(\ref{eeosbfpolye}))
	 *		[unit: $\rho_{\rm nuc} c^2$], where
	 *		$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ 
	 */
	Eos_bf_poly(double kappa1, double kappa2, double kappa3,
		    double beta) ;

	/** Standard constructor with all parameters. 
	 * 
	 *  @param gamma1  adiabatic index $\gamma_1$ 
	 *  @param gamma2  adiabatic index $\gamma_2$ 
	 *  @param gamma3  adiabatic index $\gamma_3$ 
	 *  @param gamma4  adiabatic index $\gamma_4$ 
	 *  @param gamma5  adiabatic index $\gamma_5$ 
	 *  @param gamma6  adiabatic index $\gamma_6$ 
	 *				(cf. Eq.~(\ref{eeosbfpolye}))
	 *  @param kappa1  pressure coefficient $\kappa_1$  
	 *  @param kappa2  pressure coefficient $\kappa_2$  
	 *  @param kappa3  pressure coefficient $\kappa_3$  
	 *  @param beta coefficient in the entrainment term $\beta$  
	 *		(cf. Eq.~(\ref{eeosbfpolye}))
	 *		[unit: $\rho_{\rm nuc} c^2$], where
	 *		$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ 
	 *  @param mass1  individual particule mass $m_1$ (neutrons)  
	 *  @param mass2  individual particule mass $m_2$ (protons)  
	 *		[unit: $m_B = 1.66\ 10^{-27} \ {\rm kg}$]
	 *  @param relax relaxation parameter (see {\tt par\_inv})
	 *  @param precis precision parameter for zerosec\_b 
	 *                (see {\tt par\_inv})
	 *  @param relax precision parameter for relaxation 
	 *               procedure (see {\tt par\_inv})
	 *		 
	 */
	Eos_bf_poly(double gamma1, double gamma2, double gamma3,
		    double gamma4, double gamma5, double gamma6,
		    double kappa1, double kappa2, double kappa3,
		    double beta, double mass1, double mass2, 
		    double relax=0.5, double precis = 1.e-9,
		    double ecart = 1.e-8) ;	

	Eos_bf_poly(const Eos_bf_poly& ) ;	/// Copy constructor	
	
    protected:
	/** Constructor from a binary file (created by the function 
	 *  {\tt sauve(FILE* )}). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  {\tt Eos\_bifluid::eos\_from\_file(FILE* )}. 
	 */
	Eos_bf_poly(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  {\tt Eos\_bifluid::eos\_from\_file(ifstream\& )}. 
	 */
	Eos_bf_poly(ifstream& ) ; 
	
	/// The construction functions from a file
	friend Eos_bifluid* Eos_bifluid::eos_from_file(FILE* ) ; 
	friend Eos_bifluid* Eos_bifluid::eos_from_file(ifstream& ) ; 

        public:
	virtual ~Eos_bf_poly() ;			/// Destructor

    // Assignment
    // ----------
	/// Assignment to another {\tt Eos\_bf\_poly}
	void operator=(const Eos_bf_poly& ) ;


    // Miscellaneous
    // -------------
	public : 
	/// Comparison operator (egality)
	virtual bool operator==(const Eos_bifluid& ) const ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos_bifluid& ) const ; 
    
	/** Returns a number to identify the sub-classe of {\tt Eos\_bifluid} 
	 *  the object belongs to. 
	 */
	virtual int identify() const ; 

	/// Returns the adiabatic index $\gamma_1$ 
	double get_gam1() const {return gam1 ;};

	/// Returns the adiabatic index $\gamma_2$ 
	double get_gam2() const {return gam2 ;};

	/// Returns the adiabatic index $\gamma_3$ 
	double get_gam3() const {return gam3 ;};
	
	/// Returns the adiabatic index $\gamma_4$ 
	double get_gam4() const {return gam4 ;};
	
	/// Returns the adiabatic index $\gamma_5$ 
	double get_gam5() const {return gam5 ;};
	
	/// Returns the adiabatic index $\gamma_6$ 
	double get_gam6() const {return gam6 ;};
	
	/** Returns the pressure coefficient $\kappa_1$  
	 *  [unit: $\rho_{\rm nuc} c^2 $], where
	 *  $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$.
	 */
	double get_kap1() const {return kap1 ;}; 

	/** Returns the pressure coefficient $\kappa_2$  
	 *  [unit: $\rho_{\rm nuc} c^2 $], where
	 *  $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$.
	 */
	double get_kap2() const {return kap2 ;}; 

	/** Returns the pressure coefficient $\kappa_3$  
	 *  [unit: $\rho_{\rm nuc} c^2 $], where
	 *  $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$.
	 */
	double get_kap3() const {return kap3 ;}; 

	/** Returns the coefficient $\beta$  
	 *  [unit: $\rho_{\rm nuc} c^2 $], where
	 *  $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$.
	 */
	double get_beta() const {return beta ;}; 

	/** Return the individual particule mass $m_1$  
	 *  
	 *  [unit: $m_B = 1.66\ 10^{-27} \ {\rm kg}$]. 
	 */
	double get_m1() const {return m_1 ;}; 

	/** Return the individual particule mass $m_2$  
	 *  
	 *  [unit: $m_B = 1.66\ 10^{-27} \ {\rm kg}$]. 
	 */
	double get_m2() const {return m_2 ;}; 

    protected:
	/** Computes the auxiliary quantities {\tt gam1m1}, {\tt gam2m1}
	 *  and {\tt gam3m1}
	 */
	void set_auxiliary() ; 
    
	/// Determines the type of the analytical EOS (see {\tt typeos})
	void determine_type() ;

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	/// Save in a file

    protected: 
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


    // Computational functions
    // -----------------------

    public: 

 	/** Computes both baryon densities from the log-enthalpies
	 * 
	 *  @param ent1 [input,  unit: $c^2$] log-enthalpy $H_1$ 
	 *  @param ent2 [input,  unit: $c^2$] log-enthalpy $H_2$ 
	 *  @param delta2 [input,  unit: $c^2$] relative velocity $\Delta^2$ 
	 * 
	 *  @param nbar1 [output] baryonic density of the first fluid
	 *  @param nbar2 [output] baryonic density of the second fluid
	 *  [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 * 
	 */
    	virtual void nbar_ent_p(const double ent1, const double ent2, 
				const double delta2, double& nbar1, 
				double& nbar2, bool tronc = true) const ; 
       
 	/** Computes the total energy density from the baryonic densities
	 *  and the relative velocity. 
	 * 
	 *  @param nbar1 [input] baryonic density of the first fluid
	 *  @param nbar2 [input] baryonic density of the second fluid
	 *  [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *  @param delta2 [input,  unit: $c^2$] relative velocity $\Delta^2$ 
	 * 
	 *  @return energy density $\cal E$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double ener_nbar_p(const double nbar1, const double nbar2, 
				   const double delta2) const  ; 
       
 	/** Computes the pressure from the baryonic densities
	 *  and the relative velocity.
	 * 
	 *  @param nbar1 [input] baryonic density of the first fluid
	 *  @param nbar2 [input] baryonic density of the second fluid
	 *  [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *  @param delta2 [input,  unit: $c^2$] relative velocity $\Delta^2$ 
	 * 
	 *  @return pressure $p$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double press_nbar_p(const double nbar1, const double nbar2, 
				    const double delta2) const ; 
     // Conversion functions
     // ---------------------

	/** Makes a translation from {\tt Eos\_bifluid} to {\tt Eos}. 
	 *
	 *  This is only useful for the construction of a 
	 *  {\tt Et\_rot\_bifluid}
	 *  star and ought not to be used in other situations.
	 */
	virtual Eos* trans2Eos() const ;
       
	/** Computes the derivative of the energy with respect to
	 * (baryonic density 1)$^2$.
	 *
	 *  @param n1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit $n_{\rm nuc}^2c^2$]
	 *           relative Lorentz factor$\times$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative $K^{11}=2\frac{\partial{\cal{E}}}{\partial{n_1^2}}$ 
	 */
	virtual double get_K11(const double n1, const double n2, const
			       double delta2)  const  ;

	/** Computes the derivative of the energy with respect to 
	 *  $x^2=n_1n_2\Gamma_\Delta$.
	 *
	 *  @param n1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit $n_{\rm nuc}^2c^2$]
	 *           relative Lorentz factor$\times$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative $K^{12}=\frac{\partial {\cal E}}{\partial (n_1n_2\Gamma_\Delta)}$ 
	 */
	virtual double get_K12(const double n1, const double n2,const
			       double delta2) const ;

	/** Computes the derivative of the energy/(baryonic density 2)$^2$.
	 *
	 *  @param n1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit $n_{\rm nuc}^2c^2$]
	 *           relative Lorentz factor$\times$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative $K^{22} = 2\frac{\partial {\cal E}}{\partial n_2^2}$ 
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
 * particles of rest mass is $m_1$ and $m_2$,  whose total energy density 
 * $\cal{E}$ is related to their numerical densities $n_1$, $n_2$ and 
 * relative velocity $\Delta^2$
 * \begin{equation}
 * \Delta = \left( \vec{v}_n - \vec{v}_p \right)^2
 * \label{e:defdeltan}
 * \end{equation}
 * by
 * \begin{equation} \label{eeosbfnewte}
 *   {\cal{E}} = \frac{1}{2}\kappa_1 n_1^{\gamma_1} 
 *          \  + \frac{1}{2}\kappa_2 n_2^{\gamma_2} 
 *          \  + \kappa_3 n_1^{\gamma_3} n_2^{\gamma_4} 
 *          \  + \Delta^2 \beta n_1^{\gamma_5} n_2^{\gamma_6}\ .  
 * \end{equation}
 * The non-relativistic chemical potentials are then
 * \begin{equation} 
 * \mu_1 := {{\rm d}{\cal{E}} \over {\rm d}n_1} = \frac{1}{2}\gamma_1\kappa_1
 *         n_1^{\gamma_1-1} + \gamma_3 \kappa_3 
 *         n_1^{\gamma_3-1} n_2^{\gamma_4} + \Delta^2 \gamma_5 \beta
 *         n_1^{\gamma_5-1} n_2^{\gamma_6}\ , 
 * \end{equation}
 * \begin{equation}
 * \mu_2 := {{\rm d}{\cal{E}} \over {\rm d}n_2} = \frac{1}{2}\gamma_2\kappa_2
 *         n_2^{\gamma_2-1} + \gamma_4 \kappa_3 
 *         n_1^{\gamma_3} n_2^{\gamma_4-1} + \Delta^2 \gamma_6 \beta
 *         n_1^{\gamma_5} n_2^{\gamma_6-1} \ .
 * \end{equation}
 * The pressure is given by the (zero-temperature) First Law of Thermodynamics:
 * $p = \mu_1 n_1 + \mu_2 n_2 - {\cal E}$, so that
 * \begin{equation} 
 *   p = \frac{1}{2} (\gamma_1 -1)\kappa_1 n_1^{\gamma_1} +
 *  \frac{1}{2}(\gamma_2-1)\kappa_2 n_2^{\gamma_2} + (\gamma_3 +\gamma_4
 *  -1)\kappa_3 n_1^{\gamma_3}n_2^{\gamma_4} + \Delta^2 \beta \left(
 *  (\gamma_5 + \gamma_6 - 1) n_1^{\gamma_5} n_2^{\gamma_6} \right) \ .  
 * \end{equation}
 * The
 * specific enthalpies are related to the chemical potentials by
 * \begin{equation}
 * h_i = \frac{\mu_i}{m_i}
 * \end{equation}
 *
 * From this system, the particle densities are obtained in term of 
 * the enthalpies. (The system is a linear one
 * if $\gamma_1 = \gamma_2 = 2$ and $\gamma_3 = \gamma_4 = \gamma_5 =
 * \gamma_6 = 1$).
 *
 * The energy density $\cal E$and pressure $p$ can then be obtained.
 *
 * @version #$Id$
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
	 *  The adiabatic indexes $\gamma_1$ and $\gamma_2$ are set to 2.
	 *  All other adiabatic indexes $\gamma_i$, $i\in 3\dots 6$ are
	 *  set to 1.
	 *  The individual particle masses $m_1$ and $m_2$ are set to the 
	 *  mean baryon mass $m_B = 1.66\ 10^{-27} \ {\rm kg}$. 
	 *  The inversion parameters are set to their default values
	 *  (see hereafter the consrtuctor with all parameters).
	 *  
	 *  @param kappa1  pressure coefficient $\kappa_1$  
	 *  @param kappa2  pressure coefficient $\kappa_2$  
	 *  @param kappa3  pressure coefficient $\kappa_3$  
	 *  @param beta coefficient in the entrainment term $\beta$  
	 *		(cf. Eq.~(\ref{eeosbfpolye}))
	 *		[unit: $\rho_{\rm nuc} c^2$], where
	 *		$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ 
	 */
	Eos_bf_poly_newt(double kappa1, double kappa2, double kappa3,
		    double beta) ;

	/** Standard constructor with all parameters. 
	 * 
	 *  @param gamma1  adiabatic index $\gamma_1$ 
	 *  @param gamma2  adiabatic index $\gamma_2$ 
	 *  @param gamma3  adiabatic index $\gamma_3$ 
	 *  @param gamma4  adiabatic index $\gamma_4$ 
	 *  @param gamma5  adiabatic index $\gamma_5$ 
	 *  @param gamma6  adiabatic index $\gamma_6$ 
	 *				(cf. Eq.~(\ref{eeosbfpolye}))
	 *  @param kappa1  pressure coefficient $\kappa_1$  
	 *  @param kappa2  pressure coefficient $\kappa_2$  
	 *  @param kappa3  pressure coefficient $\kappa_3$  
	 *  @param beta coefficient in the entrainment term $\beta$  
	 *		(cf. Eq.~(\ref{eeosbfpolye}))
	 *		[unit: $\rho_{\rm nuc} c^2$], where
	 *		$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ 
	 *  @param mass1  individual particule mass $m_1$ (neutrons)  
	 *  @param mass2  individual particule mass $m_2$ (protons)  
	 *  @param relax relaxation parameter (see {\tt par\_inv})
	 *  @param precis precision parameter for zerosec\_b 
	 *                (see {\tt par\_inv})
	 *  @param relax precision parameter for relaxation 
	 *               procedure (see {\tt par\_inv})
	 *		 
	 *		[unit: $m_B = 1.66\ 10^{-27} \ {\rm kg}$]
	 */
	Eos_bf_poly_newt(double gamma1, double gamma2, double gamma3,
			 double gamma4, double gamma5, double gamma6,
			 double kappa1, double kappa2, double kappa3,
			 double beta, double mass1, double mass2, 
			 double relax=0.5, double precis = 1.e-9,
			 double ecart = 1.e-8) ;	


	Eos_bf_poly_newt(const Eos_bf_poly_newt& ) ;	/// Copy constructor	
	
    protected:
	/** Constructor from a binary file (created by the function 
	 *  {\tt sauve(FILE* )}). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  {\tt Eos\_bifluid::eos\_from\_file(FILE* )}. 
	 */
	Eos_bf_poly_newt(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  {\tt Eos\_bifluid::eos\_from\_file(ifstream\& )}. 
	 */
	Eos_bf_poly_newt(ifstream& ) ; 
	
	/// The construction functions from a file
	friend Eos_bifluid* Eos_bifluid::eos_from_file(FILE* ) ; 
	friend Eos_bifluid* Eos_bifluid::eos_from_file(ifstream& ) ; 

    public:
	virtual ~Eos_bf_poly_newt() ;			/// Destructor

    // Assignment
    // ----------
	/// Assignment to another {\tt Eos\_bf\_poly\_newt}
	void operator=(const Eos_bf_poly_newt& ) ;


    // Miscellaneous
    // -------------

    public : 
	/// Comparison operator (egality)
	virtual bool operator==(const Eos_bifluid& ) const ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos_bifluid& ) const ; 
    
	/** Returns a number to identify the sub-classe of {\tt Eos\_bifluid} 
	 *  the object belongs to. 
	 */
	virtual int identify() const ; 

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	/// Save in a file

    protected: 
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


    // Computational functions
    // -----------------------

    public: 

 	/** Computes both baryon densities from the log-enthalpies
	 * 
	 *  @param ent1 [input,  unit: $c^2$] log-enthalpy $H_1$ 
	 *  @param ent2 [input,  unit: $c^2$] log-enthalpy $H_2$ 
	 *  @param delta2 [input,  unit: $c^2$] relative velocity $\Delta^2$ 
	 * 
	 *  @param nbar1 [output] baryonic density of the first fluid
	 *  @param nbar2 [output] baryonic density of the second fluid
	 *  [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 * 
	 */
    	virtual void nbar_ent_p(const double ent1, const double ent2, 
				const double delta2, double& nbar1, 
				double& nbar2, bool tronc = true) const ; 
       
 	/** Computes the total energy density from the baryonic densities
	 *  and the relative velocity. 
	 * 
	 *  @param nbar1 [input] baryonic density of the first fluid
	 *  @param nbar2 [input] baryonic density of the second fluid
	 *  [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *  @param delta2 [input,  unit: $c^2$] relative velocity $\Delta^2$ 
	 * 
	 *  @return energy density $\cal E$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double ener_nbar_p(const double nbar1, const double nbar2, 
				   const double delta2) const  ; 
       
 	/** Computes the pressure from the baryonic densities
	 *  and the relative velocity.
	 * 
	 *  @param nbar1 [input] baryonic density of the first fluid
	 *  @param nbar2 [input] baryonic density of the second fluid
	 *  [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *  @param delta2 [input,  unit: $c^2$] relative velocity $\Delta^2$ 
	 * 
	 *  @return pressure $p$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double press_nbar_p(const double nbar1, const double nbar2, 
				    const double delta2) const ; 
     // Conversion functions
     // ---------------------

	/** Makes a translation from {\tt Eos\_bifluid} to {\tt Eos}. 
	 *
	 *  This is only useful for the construction of a 
	 *  {\tt Et\_rot\_bifluid}
	 *  star and ought not to be used in other situations.
	 */
	virtual Eos* trans2Eos() const ;
       
	/** Computes the derivative of the energy with respect to
	 * (baryonic density 1)$^2$.
	 *
	 *  @param n1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit $n_{\rm nuc}^2c^2$]
	 *           relative Lorentz factor$\times$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative $K^{11}=2\frac{\partial{\cal{E}}}{\partial{n_1^2}}$ 
	 */
	virtual double get_K11(const double n1, const double n2, const
			       double delta2)  const  ;

	/** Computes the derivative of the energy with respect to 
	 *  $x^2=n_1n_2\Gamma_\Delta$.
	 *
	 *  @param n1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit $n_{\rm nuc}^2c^2$]
	 *           relative Lorentz factor$\times$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative $K^{12}=\frac{\partial {\cal E}}{\partial (n_1n_2\Gamma_\Delta)}$ 
	 */
	virtual double get_K12(const double n1, const double n2,const
			       double delta2) const ;

	/** Computes the derivative of the energy/(baryonic density 2)$^2$.
	 *
	 *  @param n1 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 1 at which the derivative 
	 *           is to be computed
	 *  @param n2 [input, unit $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *           baryonic density of fluid 2 at which the derivative 
	 *           is to be computed
	 *  @param x [input, unit $n_{\rm nuc}^2c^2$]
	 *           relative Lorentz factor$\times$both densities at which 
	 *           the derivative is to be computed
	 *
	 *  @return derivative $K^{22} = 2\frac{\partial {\cal E}}{\partial n_2^2}$ 
	 */
	virtual double get_K22(const double n1, const double n2, const
			       double delta2) const ;

};

		

#endif
