/*
 *  Definition of Lorene classes Eos
 *				 Eos_poly
 *				 Eos_poly_newt
 *				 Eos_incomp
 *				 Eos_incomp_newt
 *				 Eos_strange
 *				 Eos_strange_c
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


#ifndef __EOS_H_ 
#define __EOS_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2002/06/17 14:05:16  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.4  2002/04/11 13:27:48  e_gourgoulhon
 * *** empty log message ***
 *
 * Revision 1.3  2002/04/09 15:19:03  e_gourgoulhon
 * Add EOS number 100 in the comments of eos_from_file
 *
 * Revision 1.2  2002/04/09 14:32:14  e_gourgoulhon
 * 1/ Added extra parameters in EOS computational functions (argument par)
 * 2/ New class MEos for multi-domain EOS
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.15  2001/09/12  15:54:17  eric
 * Modif documentation eos_from_file.
 *
 * Revision 2.14  2001/06/13  14:12:18  eric
 *  Modif commentaires (mise en conformite Doc++ 3.4.7)
 *
 * Revision 2.13  2001/02/07  09:33:42  eric
 * Suppression des fonctions derent_ent et derent_ent_p.
 * Ajout des fonctions donnant les derivees de l'EOS:
 * 	der_nbar_ent
 * 	der_ener_ent
 * 	der_press_ent
 *
 * Revision 2.12  2000/11/23  14:45:33  eric
 * Ajout de l'EOS Eos_strange_cr.
 *
 * Revision 2.11  2000/11/22  19:28:45  eric
 * Ajout de #include "eos_tabul.h" a la fin.
 *
 * Revision 2.10  2000/10/25  10:55:08  eric
 * Eos_strange: modif commentaires.
 *
 * Revision 2.9  2000/10/24  15:28:43  eric
 * Ajout de l'EOS matiere etrange (Eos_strange).
 *
 * Revision 2.8  2000/06/20  08:34:20  eric
 * Ajout des membres get_gam(), ... a Eos_ploy
 *
 * Revision 2.7  2000/02/14  14:43:22  eric
 * Modif commentaires.
 *
 * Revision 2.6  2000/02/14  14:32:46  eric
 * Ajout des constructeurs par lecture de fichier formate.
 *
 * Revision 2.5  2000/01/21  16:21:12  eric
 * Modif commentaires.
 *
 * Revision 2.4  2000/01/21  15:16:10  eric
 * Ajout de la fonction identify()
 * Ajout de la fonction de construction a partir d'un fichier
 *   static Eos* Eos::eos_from_file(FILE* ).
 * Ajout des operateurs de comparaison == et !=
 *
 * Revision 2.3  2000/01/18  16:10:57  eric
 * Ajout des EOS Eos_incomp et Eos_incomp_newt.
 *
 * Revision 2.2  2000/01/18  15:13:28  eric
 * Ajout de l'equation d'etat Eos_poly_newt.
 *
 * Revision 2.1  2000/01/18  13:46:50  eric
 * Premiere version operationnelle
 *
 * Revision 2.0  2000/01/18  10:46:08  eric
 * *** empty log message ***
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
class Tbl ;
class Cmp ;
class Param ;

		    //------------------------------------//
		    //		base class Eos		  //
		    //------------------------------------//

/**
 * Equation of state base class.
 * 
 *
 * @version #$Id$#
 */
class Eos {

    // Data :
    // -----

    protected: 
	char name[100] ;	    /// EOS name


    // Constructors - Destructor
    // -------------------------
    protected:
	Eos() ;			/// Standard constructor

	/// Standard constructor with name
	explicit Eos(const char* name_i) ; 

	Eos(const Eos& ) ;	/// Copy constructor	

    protected:
	/** Constructor from a binary file (created by the function 
	 *  {\tt sauve(FILE* )}). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  {\tt Eos::eos\_from\_file(FILE* )}. 
	 */
	Eos(FILE* ) ; 

	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  {\tt Eos::eos\_from\_file(ifstream\&)}. 
	 */
	Eos(ifstream& ) ; 
	
	
    public:
	virtual ~Eos() ;			/// Destructor


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
	 *  The file must have been created by the function {\tt sauve(FILE* )}.
	 */
	static Eos* eos_from_file(FILE* ) ; 
	
	/** Construction of an EOS from a formatted file.
	 * 
	 *  The fist line of the file must start by the EOS number, according 
	 *  to the following conventions: \\
	 *	1 = relativistic polytropic EOS (class {\tt Eos\_poly}). \\
	 *      2 = Newtonian polytropic EOS (class {\tt Eos\_poly\_newt}). \\
	 *	3 = Relativistic incompressible EOS (class {\tt Eos\_incomp}). \\
	 *	4 = Newtonian incompressible EOS 
	 *		    (class {\tt Eos\_incomp\_newt}). \\
	 *	5 = Strange matter (MIT Bag model) \\
	 *	6 = Strange matter (MIT Bag model) with crust \\
	 *     10 = SLy4 (Douchin \& Haensel 2001)  \\
	 *     11 = FPS (Friedman-Pandharipande + Skyrme) \\
	 *     12 = BPAL12 (Bombaci et al. 1995) \\
	 *     13 = AkmalPR (Akmal, Pandharipande \& Ravenhall 1998) \\
	 *     14 = BBB2 (Baldo, Bombaci \& Burgio 1997) \\
	 *     15 = BalbN1H1 (Balberg 2000) \\
	 *     100 = Multi-domain EOS (class {\tt MEos}) \\
	 *  The second line in the file must be the EOS name.
	 *  The following lines should contain the EOS parameters (one
	 *  parameter per line), in the same order than in the class declaration.
	 */
	static Eos* eos_from_file(ifstream& ) ; 
	
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const = 0 ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const = 0 ; 
    
	/** Returns a number to identify the sub-classe of {\tt Eos} the
	 *  object belongs to. 
	 */
	virtual int identify() const = 0 ; 

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	/// Save in a file

	/// Display
	friend ostream& operator<<(ostream& , const Eos& ) ;	

    protected: 
	virtual ostream& operator>>(ostream &) const = 0 ;    /// Operator >>


    // Computational functions
    // -----------------------
    protected:
	/**  General computational method for {\tt Cmp}'s
	 *
	 *   @param thermo [input] thermodynamical quantity (for instance the
	 *	    enthalpy field)from which the
	 *          thermodynamical quantity {\tt resu} is to be computed.
	 *  @param nzet  [input] number of domains where {\tt resu} is to be
	 *	computed.
	 *  @param l_min [input] index of the innermost domain is which {\tt resu}
	 *	is to be computed [default value: 0]; {\tt resu} is
	 *	computed only in domains whose indices are in
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero.
	 *  @param fait [input] pointer on the member function of class
	 *		{\tt Eos} which performs the pointwise calculation.
         * @param par possible extra parameters of the EOS
	 *  @param resu [output] result of the computation.
	 */
	void calcule(const Cmp& thermo, int nzet, int l_min,
		     double (Eos::*fait)(double, const Param*) const, const Param* par, Cmp& resu) const ;


    public:
 	/** Computes the baryon density from the log-enthalpy and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *    $H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) $,
	 *    where $e$ is the (total) energy density, $p$ the pressure,
	 *    $n$ the baryon density, and $m_B$ the baryon mass
        *  @param par possible extra parameters of the EOS
	 *
	 *  @return baryon density [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const = 0 ;

	/** Computes the baryon density field from the log-enthalpy field and
        * extra parameters
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *    $H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) $,
	 *    where $e$ is the (total) energy density, $p$ the pressure,
	 *    $n$ the baryon density, and $m_B$ the baryon mass
	 *  @param nzet  number of domains where the baryon density is to be
	 *	computed.
	 *  @param l_min  index of the innermost domain is which the baryon
	 *	density is
	 *	to be computed [default value: 0]; the baryon density is
	 *	computed only in domains whose indices are in
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero.
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return baryon density [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *
	 */
    	Cmp nbar_ent(const Cmp& ent, int nzet, int l_min = 0, const Param* par=0x0) const  ;

 	/** Computes the total energy density from the log-enthalpy and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *    $H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) $,
	 *    where $e$ is the (total) energy density, $p$ the pressure,
	 *    $n$ the baryon density, and $m_B$ the baryon mass
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return energy density $e$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const = 0 ;

	/** Computes the total energy density from the log-enthalpy and extra parameters.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *    $H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) $,
	 *    where $e$ is the (total) energy density, $p$ the pressure,
	 *    $n$ the baryon density, and $m_B$ the baryon mass
	 *  @param nzet  number of domains where the energy density is to be
	 *	computed.
	 *  @param l_min  index of the innermost domain is which the energy
	 *	density is
	 *	to be computed [default value: 0]; the energy density is
	 *	computed only in domains whose indices are in
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero.
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return energy density [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	Cmp ener_ent(const Cmp& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

	/** Computes the pressure from the log-enthalpy and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *    $H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) $,
	 *    where $e$ is the (total) energy density, $p$ the pressure,
	 *    $n$ the baryon density, and $m_B$ the baryon mass
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return pressure $p$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const = 0 ;


	/** Computes the pressure from the log-enthalpy and extra parameters
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *    $H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) $,
	 *    where $e$ is the (total) energy density, $p$ the pressure,
	 *    $n$ the baryon density, and $m_B$ the baryon mass
	 *  @param nzet  number of domains where the pressure is to be
	 *	computed.
	 *  @param l_min  index of the innermost domain is which the pressure is
	 *	to be computed [default value: 0]; the pressure is computed
	 *      only in domains whose indices are in
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero.
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return pressure [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 *
	 */
    	Cmp press_ent(const Cmp& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative $d\ln n/d\ln H$
	 *  from the log-enthalpy and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *    $H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) $,
	 *    where $e$ is the (total) energy density, $p$ the pressure,
	 *    $n$ the baryon density, and $m_B$ the baryon mass
         *  @param par possible extra parameters of the EOS
         *
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const = 0 ;

	/** Computes the logarithmic derivative $d\ln n/d\ln H$
	 *  from the log-enthalpy and extra parameters
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *    $H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) $,
	 *    where $e$ is the (total) energy density, $p$ the pressure,
	 *    $n$ the baryon density, and $m_B$ the baryon mass
	 *  @param nzet  number of domains where the derivative
	 *	dln(n)/dln(H) is to be computed.
	 *  @param l_min  index of the innermost domain is which the
	 *	   coefficient dln(n)/dln(H) is
	 *	to be computed [default value: 0]; the derivative
	 *	dln(n)/dln(H) is
	 *	computed only in domains whose indices are in
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero.
         *  @param par possible extra parameters of the EOS
	 *
	 *  @return dln(n)/dln(H)
	 *
	 */
    	Cmp der_nbar_ent(const Cmp& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative $d\ln e/d\ln H$
	 *  from the log-enthalpy with extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *    $H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) $,
	 *    where $e$ is the (total) energy density, $p$ the pressure,
	 *    $n$ the baryon density, and $m_B$ the baryon mass
         *  @param par possible extra parameters of the EOS
  	 *
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const = 0 ;

	/** Computes the logarithmic derivative $d\ln e/d\ln H$
	 *  from the log-enthalpy and extra parameters
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *    $H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) $,
	 *    where $e$ is the (total) energy density, $p$ the pressure,
	 *    $n$ the baryon density, and $m_B$ the baryon mass
	 *  @param nzet  number of domains where the derivative
	 *	dln(e)/dln(H) is to be computed.
	 *  @param l_min  index of the innermost domain is which the
	 *	   coefficient dln(n)/dln(H) is
	 *	to be computed [default value: 0]; the derivative
	 *	dln(e)/dln(H) is
	 *	computed only in domains whose indices are in
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero.
         *  @param par possible extra parameters of the EOS
	 *
	 *  @return dln(e)/dln(H)
	 *
	 */
    	Cmp der_ener_ent(const Cmp& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative $d\ln p/d\ln H$
	 *  from the log-enthalpy and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *    $H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) $,
	 *    where $e$ is the (total) energy density, $p$ the pressure,
	 *    $n$ the baryon density, and $m_B$ the baryon mass
        *  @param par possible extra parameters of the EOS
	 *
	 *  @return dln(p)/dln(H)
	 */
    	virtual double der_press_ent_p(double ent, const Param* par=0x0) const = 0 ;

	/** Computes the logarithmic derivative $d\ln p/d\ln H$
	 *  from the log-enthalpy and extra parameters
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *    $H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) $,
	 *    where $e$ is the (total) energy density, $p$ the pressure,
	 *    $n$ the baryon density, and $m_B$ the baryon mass
	 *  @param nzet  number of domains where the derivative
	 *	dln(p)/dln(H) is to be computed.
        *  @param par possible extra parameters of the EOS
	 *  @param l_min  index of the innermost domain is which the
	 *	   coefficient dln(n)/dln(H) is
	 *	to be computed [default value: 0]; the derivative
	 *	dln(p)/dln(H) is
	 *	computed only in domains whose indices are in
	 *      {\tt [l\_min, l\_min + nzet-1]}. In the other
	 *	domains, it is set to zero.
	 *
	 *  @return dln(p)/dln(H)
	 *
	 */
    	Cmp der_press_ent(const Cmp& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

};
ostream& operator<<(ostream& , const Eos& ) ;	


		    //------------------------------------//
		    //		class Eos_poly		  //
		    //------------------------------------//


/**
 * Polytropic equation of state (relativistic case).
 *
 * This equation of state (EOS) corresponds to identical relativistic
 * particles of rest mass is $m_0$,  whose total energy density $e$ is
 * related to their numerical density $n$ by
 * \begin{equation} \label{eeospolye}
 *   e(n) = {\kappa \over \gamma-1} n^\gamma + \mu_0 \, n \ ,
 * \end{equation}
 * where $\mu_0$ is the chemical potential at zero pressure.
 * The relativistic (i.e. including rest mass energy) chemical potential is
 * then
 * \begin{equation}  \label{eeospolymu}
 *   \mu(n) := {de\over dn} = {\kappa \gamma \over \gamma-1} n^{\gamma-1}
 *		+ \mu_0 \ .
 * \end{equation}
 * The pressure is given by the (zero-temperature) First Law of Thermodynamics:
 * $p = \mu n - e$, so that
 * \begin{equation} \label{eeospolyp}
 *   p(n) = \kappa n^\gamma  \ .
 * \end{equation}
 * The log-enthalpy is defined as the logarithm of the ratio of the enthalpy
 * par particle by the partical rest mass energy :
 * \begin{equation} \label{eeospolyh}
 *   H(n) := c^2 \ln \left( {e+p \over m_0 c^2\, n} \right)   \ .
 * \end{equation}
 * According to the (zero-temperature) First Law of Thermodynamics, the
 * log-enthalpy is related to the chemical potential by
 * \begin{equation}
 *   H = c^2 \ln \left( {\mu \over m_0 c^2} \right) \ .
 * \end{equation}
 * From this expression and relation (\ref{eeospolymu}), the expression
 * of the particle density in term of the log-enthalpy is
 * \begin{equation}
 *   n(H) = \left[ {\gamma-1\over \gamma} {m_0 c^2 \over \kappa}
                \left( \exp(H) - {\mu_0\over m_0 c^2} \right)
 *	    \right] ^{1/(\gamma-1)}  \ .
 * \end{equation}
 * The energy density and pressure as functions of $H$ can then be obtained
 * by inserting this relation into Eqs.~(\ref{eeospolye}) and
 * (\ref{eeospolyp}).
 *
 * @version #$Id$#
 */
class Eos_poly : public Eos {

    // Data :
    // -----

    protected:
	/// Adiabatic index $\gamma$ (cf. Eq.~(\ref{eeospolyp}))
	double gam ;

	/** Pressure coefficient $\kappa$  (cf. Eq.~(\ref{eeospolyp}))
	 *  [unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$], where
	 *  $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ and
	 *  $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$.
	 */
	double kap ; 

	/** Individual particule mass $m_0$  (cf. Eq.~(\ref{eeospolye}))
	 *  [unit: $m_B = 1.66\ 10^{-27} \ {\rm kg}$].
	 */
	double m_0 ;

        /** Relativistic chemical potential at zero pressure
	 *  [unit: $m_B c^2$, with $m_B = 1.66\ 10^{-27} \ {\rm kg}$].
         * (standard value: 1)
        */
        double mu_0 ;



	double gam1 ;	    /// $\gamma-1$
	double unsgam1 ;    /// $1/(\gamma-1)$
	double gam1sgamkap ; /// $(\gamma-1) / (\gamma \kappa) m_0$
        double rel_mu_0 ;       /// $\mu_0/m_0$
        double ent_0 ;          /// Enthalpy at zero pressure ($\ln (\mu_0/m_0)$)

    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor (sets both {\tt m\_0} and {\tt mu\_0} to 1).
	 *
	 *  The individual particle mass $m_0$ is set to the mean baryon
	 *  mass $m_B = 1.66\ 10^{-27} \ {\rm kg}$.
	 *
	 *  @param gamma  adiabatic index $\gamma$
	 *				(cf. Eq.~(\ref{eeospolyp}))
	 *  @param kappa  pressure coefficient $\kappa$  
	 *		(cf. Eq.~(\ref{eeospolyp}))
	 *		[unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$], where
	 *		$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ and
	 *		$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$
	 */
	Eos_poly(double gamma, double kappa) ;

	/** Standard constructor with individual particle mass
	*   (sets {\tt mu\_0} to 1).
	 *  @param gamma  adiabatic index $\gamma$ (cf. Eq.~(\ref{eeospolyp}))
	 *  @param kappa  pressure coefficient $\kappa$
	 *		(cf. Eq.~(\ref{eeospolyp}))
	 *		[unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$], where
	 *		$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ and
	 *		$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$
	 *  @param mass  individual particule mass $m_0$
	 *		 (cf. Eq.~(\ref{eeospolye}))
	 *		[unit: $m_B = 1.66\ 10^{-27} \ {\rm kg}$]
	 */
	Eos_poly(double gamma, double kappa, double mass) ;

	/** Standard constructor with individual particle mass and zero-pressure
         * chemical potential
	 *
	 *  @param gamma  adiabatic index $\gamma$ (cf. Eq.~(\ref{eeospolyp}))
	 *  @param kappa  pressure coefficient $\kappa$
	 *		(cf. Eq.~(\ref{eeospolyp}))
	 *		[unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$], where
	 *		$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ and
	 *		$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$
	 *  @param mass  individual particule mass $m_0$
	 *		 (cf. Eq.~(\ref{eeospolye}))
	 *		[unit: $m_B = 1.66\ 10^{-27} \ {\rm kg}$]
        *  @param mu_zero  Relativistic chemical potential at zero pressure
	 *  [unit: $m_B c^2$, with $m_B = 1.66\ 10^{-27} \ {\rm kg}$].
         * (standard value: 1)
        */
	Eos_poly(double gamma, double kappa, double mass, double mu_zero) ;

	Eos_poly(const Eos_poly& ) ;	/// Copy constructor
	
    protected:
	/** Constructor from a binary file (created by the function
	 *  {\tt sauve(FILE* )}).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  {\tt Eos::eos\_from\_file(FILE* )}.
	 */
	Eos_poly(FILE* ) ;

	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}. 
	 */
	Eos_poly(ifstream& ) ;

	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ; 

    public:
	virtual ~Eos_poly() ;			/// Destructor

    // Assignment
    // ----------
	/// Assignment to another {\tt Eos\_poly}
	void operator=(const Eos_poly& ) ;


    // Miscellaneous
    // -------------

    public :
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ;

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ;

	/** Returns a number to identify the sub-classe of {\tt Eos} the
	 *  object belongs to.
	 */
	virtual int identify() const ; 

	/// Returns the adiabatic index $\gamma$ (cf. Eq.~(\ref{eeospolyp}))
	double get_gam() const ;

	/** Returns the pressure coefficient $\kappa$  (cf. Eq.~(\ref{eeospolyp}))
	 *  [unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$], where
	 *  $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ and
	 *  $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$.
	 */
	double get_kap() const ;
	
	/** Return the individual particule mass $m_0$
	 *  (cf. Eq.~(\ref{eeospolye}))
	 *  [unit: $m_B = 1.66\ 10^{-27} \ {\rm kg}$].
	 */
	double get_m_0() const ;

	/** Return the relativistic chemical potential at zero pressure
	 *  [unit: $m_B c^2$, with $m_B = 1.66\ 10^{-27} \ {\rm kg}$].
	 */
	double get_mu_0() const ;

    protected:
	/** Computes the auxiliary quantities {\tt gam1}, {\tt unsgam1},
	 *  {\tt gam1sgamkap} from the values of {\tt gam} and {\tt kap}
	 */
	void set_auxiliary() ;


    // Outputs
    // -------

    public:
	virtual void sauve(FILE* ) const ;	/// Save in a file

    protected:
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


    // Computational functions
    // -----------------------

    public:
	/** Computes the baryon density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return baryon density $n$ [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the total energy density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return energy density $e$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the pressure from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return pressure $p$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative $d\ln n/d\ln H$
	 * from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative $d\ln e/d\ln H$
	 * from the log-enthalpy.
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative $d\ln p/d\ln H$
	 * from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return dln(p)/dln(H)
	 */
    	virtual double der_press_ent_p(double ent, const Param* par=0x0) const ;
       
};

		    //------------------------------------//
		    //		class Eos_poly_newt	  //
		    //------------------------------------//



/**
 * Polytropic equation of state (Newtonian case).
 *
 * This equation of state (EOS) corresponds to identical non relativistic
 * particles of rest mass is $m_0$,  whose internal energy density $\epsilon$ is
 * related to their numerical density $n$ by
 * \begin{equation} \label{eeospolynewte}
 *   \epsilon(n) = {\kappa \over \gamma-1} n^\gamma  \ .
 * \end{equation}
 * The (non-relativistic) chemical potential is
 * then
 * \begin{equation}  \label{eeospolynewtmu}
 *   \mu(n) := {d\epsilon\over dn} = {\kappa \gamma \over \gamma-1} n^{\gamma-1}
 *		  \ .
 * \end{equation}
 * The pressure is given by the (zero-temperature) First Law of Thermodynamics:
 * $p = \mu n - \epsilon$, so that
 * \begin{equation} \label{eeospolynewtp}
 *   p(n) = \kappa n^\gamma  \ .  
 * \end{equation}
 * The (non-relativistic) specific enthalpy is  :
 * \begin{equation} \label{eeospolynewth}
 *   h(n) := {\epsilon + p \over m_0 n}   \ .  
 * \end{equation}
 * According to the (zero-temperature) First Law of Thermodynamics, the
 * specific enthalpy is related to the chemical potential by
 * \begin{equation}
 *   h = {\mu \over m_0}  \ . 
 * \end{equation}
 * From this expression and relation (\ref{eeospolynewtmu}), the expression
 * of the particle density in term of the  specific enthalpy is 
 * \begin{equation}
 *   n(h) = \left[ {\gamma-1\over \gamma} {m_0 \over \kappa} h
 *	    \right] ^{1/(\gamma-1)}  \ .  
 * \end{equation}
 * The energy density and pressure as functions of $h$ can then be obtained
 * by inserting this relation into Eqs.~(\ref{eeospolynewte}) and 
 * (\ref{eeospolynewtp}). 
 *
 * @version #$Id$#
 */
class Eos_poly_newt : public Eos_poly {

    // Data :
    // -----

    // no new data with respect to Eos_poly	

    // Constructors - Destructor
    // -------------------------
    public:
    
	/** Standard constructor.
	 * 
	 *  The individual particle mass $m_0$ is set to the mean baryon
	 *  mass $m_B = 1.66\ 10^{-27} \ {\rm kg}$. 
	 *  
	 *  @param gamma  adiabatic index $\gamma$ 
	 *				(cf. Eq.~(\ref{eeospolynewtp}))
	 *  @param kappa  pressure coefficient $\kappa$  
	 *		(cf. Eq.~(\ref{eeospolynewtp}))
	 *		[unit: $\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma$], where
	 *		$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ and
	 *		$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$
	 */
	Eos_poly_newt(double gamma, double kappa) ;	

	Eos_poly_newt(const Eos_poly_newt& ) ;	/// Copy constructor	

    protected:
	/** Constructor from a binary file (created by the function 
	 *  {\tt sauve(FILE* )}). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  {\tt Eos::eos\_from\_file(FILE* )}. 
	 */
	Eos_poly_newt(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}. 
	 */
	Eos_poly_newt(ifstream& ) ; 
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ; 
	friend Eos* Eos::eos_from_file(ifstream& ) ; 


    public:
	virtual ~Eos_poly_newt() ;			/// Destructor

    // Assignment
    // ----------
	/// Assignment to another {\tt Eos\_poly\_newt}
	void operator=(const Eos_poly_newt& ) ;

    // Miscellaneous
    // -------------

    public : 
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ; 
    
	/** Returns a number to identify the sub-classe of {\tt Eos} the
	 *  object belongs to. 
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
	/** Computes the baryon density from the specific enthalpy.
	 * 
	 *  @param ent [input,  unit: $c^2$] specific enthalpy $h$ defined by
	 *				     Eq. (\ref{eeospolynewth})
	 *
	 *  @return baryon density $n$ [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 * 
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ; 
    
 	/** Computes the total energy density from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] specific enthalpy $h$ defined by
	 *				     Eq. (\ref{eeospolynewth})
	 *
	 *  @return energy density $e$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ; 
       
 	/** Computes the pressure from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] specific enthalpy $h$ defined by
	 *				     Eq. (\ref{eeospolynewth})
	 *
	 *  @return pressure $p$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln n/d\ln h$ 
	 * from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] specific enthalpy $h$ defined by
	 *				     Eq. (\ref{eeospolynewth})
	 *
	 *  @return dln(n)/dln(h)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln e/d\ln h$ 
	 * from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] specific enthalpy $h$ defined by
	 *				     Eq. (\ref{eeospolynewth})
	 *
	 *  @return dln(e)/dln(h)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln p/d\ln h$ 
	 * from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] specific enthalpy $h$ defined by
	 *				     Eq. (\ref{eeospolynewth})
	 *
	 *  @return dln(p)/dln(h)
	 */
    	virtual double der_press_ent_p(double ent, const Param* par=0x0) const ; 
       
};


		    //------------------------------------//
		    //		class Eos_incomp	  //
		    //------------------------------------//


/**
 * Equation of state of incompressible matter (relativistic case).
 * 
 * This equation of state (EOS) corresponds to a constant density 
 * matter $e = \rho_0 \,  c^2$. 
 * 
 *
 * @version #$Id$#
 */
class Eos_incomp : public Eos {

    // Data :
    // -----

    protected: 
	/// Constant density $\rho_0$ 
	double rho0 ;
	
	/** Log-enthalpy threshold for setting the energy density to 
	 *  a non zero value (should be negative).   
	 */
	double ent0 ; 	

    // Constructors - Destructor
    // -------------------------
    public:
    
	/** Standard constructor.
	 * 
	 *  The log-enthalpy threshold for non-zero density is set 
	 *  -1.e-6. 
	 *  
	 *  @param rho_c  constant density $\rho_0$ 
	 *		    [unit: $\rho_{\rm nuc} c^2$], where
	 *		    $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
	Eos_incomp(double rho_c) ;	

	/** Standard constructor with log-enthalpy threshold. 
	 * 
	 *  @param rho_c  constant density $\rho_0$
	 *		    [unit: $\rho_{\rm nuc} c^2$], where
	 *		    $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 *  @param ent_c  log-enthalpy threshold for non-zero density 
	 *		  [unit: $c^2$] :
	 *		  the energy density is set to {\tt rho\_c} for
	 *		  $H > {\tt ent\_c}$. 
	 */
	Eos_incomp(double rho_c, double ent_c) ;	

	Eos_incomp(const Eos_incomp& ) ;	/// Copy constructor	

    protected:
	/** Constructor from a binary file (created by the function 
	 *  {\tt sauve(FILE* )}). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  {\tt Eos::eos\_from\_file(FILE* )}. 
	 */
	Eos_incomp(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}. 
	 */
	Eos_incomp(ifstream& ) ; 
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ; 
	friend Eos* Eos::eos_from_file(ifstream& ) ; 

    public:
	virtual ~Eos_incomp() ;			/// Destructor

    // Assignment
    // ----------
	/// Assignment to another {\tt Eos\_incomp}
	void operator=(const Eos_incomp& ) ;
    

    // Miscellaneous
    // -------------

    public : 
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ; 
    
	/** Returns a number to identify the sub-classe of {\tt Eos} the
	 *  object belongs to. 
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
	/** Computes the baryon density from the log-enthalpy.
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return baryon density $n$ [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 * 
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ; 
    
 	/** Computes the total energy density from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return energy density $e$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ; 
       
 	/** Computes the pressure from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return pressure $p$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln n/d\ln H$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln e/d\ln H$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln p/d\ln H$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return dln(p)/dln(H)
	 */
    	virtual double der_press_ent_p(double ent, const Param* par=0x0) const ; 

};

		    //------------------------------------//
		    //	   class Eos_incomp_newt	  //
		    //------------------------------------//


/**
 * Equation of state of incompressible matter (Newtonian case).
 * 
 * This equation of state (EOS) corresponds to a constant density 
 * matter $\rho = \rho_0$. 
 * 
 *
 * @version #$Id$#
 */
class Eos_incomp_newt : public Eos_incomp {

    // Data :
    // -----

    // no new data with respect to Eos_incomp
    
    // Constructors - Destructor
    // -------------------------
    public:
    
	/** Standard constructor.
	 * 
	 *  The log-enthalpy threshold for non-zero density is set 
	 *  -1.e-6. 
	 *  
	 *  @param rho_c  constant density $\rho_0$ 
	 *		    [unit: $\rho_{\rm nuc} c^2$], where
	 *		    $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
	Eos_incomp_newt(double rho_c) ;	

	/** Standard constructor with specific enthalpy threshold. 
	 * 
	 *  @param rho_c  constant density $\rho_0$
	 *		    [unit: $\rho_{\rm nuc} c^2$], where
	 *		    $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 *  @param ent_c  specific enthalpy threshold for non-zero density 
	 *		  [unit: $c^2$] :
	 *		  the energy density is set to {\tt rho\_c} for
	 *		  $h > {\tt ent\_c}$. 
	 */
	Eos_incomp_newt(double rho_c, double ent_c) ;	

	Eos_incomp_newt(const Eos_incomp_newt& ) ;	/// Copy constructor	

    protected:
	/** Constructor from a binary file (created by the function 
	 *  {\tt sauve(FILE* )}). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  {\tt Eos::eos\_from\_file(FILE* )}. 
	 */
	Eos_incomp_newt(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}. 
	 */
	Eos_incomp_newt(ifstream& ) ; 
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ; 
	friend Eos* Eos::eos_from_file(ifstream& ) ; 

    public:
	virtual ~Eos_incomp_newt() ;			/// Destructor

    // Assignment
    // ----------
	/// Assignment to another {\tt Eos\_incomp\_newt}
	void operator=(const Eos_incomp_newt& ) ;
    

    // Miscellaneous
    // -------------

    public : 
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ; 
    
	/** Returns a number to identify the sub-classe of {\tt Eos} the
	 *  object belongs to. 
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
	/** Computes the baryon density from the specific enthalpy.
	 * 
	 *  @param ent [input,  unit: $c^2$] specific enthalpy $h$
	 *
	 *  @return baryon density $n$ [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 * 
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ; 
    
 	/** Computes the total energy density from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] specific enthalpy $h$ 
	 *
	 *  @return energy density $e$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ; 
       
 	/** Computes the pressure from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] specific enthalpy $h$ 
	 *
	 *  @return pressure $p$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln n/d\ln h$ 
	 * from the specific enthalpy.
	 * 
	 *  @param ent [input,  unit: $c^2$] specific enthalpy $h$ 
	 *
	 *  @return dln(n)/dln(h)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln e/d\ln h$ 
	 * from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] specific enthalpy $h$ 
	 *
	 *  @return dln(e)/dln(h)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln p/d\ln h$ 
	 * from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] specific enthalpy $h$ 
	 *
	 *  @return dln(p)/dln(h)
	 */
    	virtual double der_press_ent_p(double ent, const Param* par=0x0) const ; 
       
};


		    //------------------------------------//
		    //		class Eos_strange	  //
		    //------------------------------------//


/**
 * Strange matter EOS (MIT Bag model).
 * 
 * This equation of state (EOS) corresponds to u,d,s degenerated symetric
 * matter in the MIT bag model, according to approximate formula
 * given in Zdunik, Astron. Astrophys. {\bf 359}, 311 (2000). 
 *
 * @version #$Id$#
 */
class Eos_strange : public Eos {

    // Data :
    // -----

    protected: 
	/** Baryon density at zero pressure divided by $B_{60}^{3/4}$.
	 *   [unit: ${\rm fm}^{-3}$]
	 */
	double n0_b60 ; 
	
	/// Bag constant [unit: $60\ {\rm MeV\, fm}^{-3}$]
	double b60 ;
	
	/** Log-enthalpy threshold for setting the energy density to 
	 *  a non zero value (should be negative).   
	 */
	double ent0 ; 	
	
	/** Fitting parameter $\epsilon_{\rm fit}$ related to the
	 *  square of sound velocity by $c_s^2 = 1/3(1+\epsilon_{\rm fit})$.
	 *  [cf. Zdunik, Astron. Astrophys. {\bf 359}, 311 (2000)]
	 */
	double eps_fit ; 
	
	/** Energy density at zero pressure divided by $B_{60}$.
	 *  [unit: ${\rm MeV\,  fm^{-3}}$]
	 */
	double rho0_b60 ; 
	
	/** Baryon density at zero pressure.
	 *   [unit: $0.1{\rm \  fm}^{-3}$ (Lorene's unit)]
	 */
	double n0 ; 
	
	/** Energy density at zero pressure. 
	 *  [unit: $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$ 
	 *   (Lorene's unit)]
	 */
	double rho0 ; 

	/** $B_{60}^{3/4}$
	 * 
	 */
	double b34 ;
	
	/**  Factor $(4+\epsilon_{\rm fit})/(1+\epsilon_{\rm fit})$
	 * 
	 */
	double fach ; 

    // Constructors - Destructor
    // -------------------------
    public:
    
	/** Standard constructor.
	 * 
	 *  @param n0_b60_i Baryon density at zero pressure divided 
	 *		    by $B_{60}^{3/4}$ [unit: ${\rm fm}^{-3}$]
	 *  @param b60_i Bag constant [unit: $60\ {\rm MeV\, fm}^{-3}$]
	 *  @param ent0_i Log-enthalpy threshold for setting the energy 
	 *		    density to a non zero value (should be negative)
	 *  @param eps_fit_i Fitting parameter $\epsilon_{\rm fit}$ related 
	 *		     to the square of sound velocity by 
	 *		     $c_s^2 = 1/3(1+\epsilon_{\rm fit})$
	 *		[cf. Zdunik, Astron. Astrophys. {\bf 359}, 311 (2000)]
	 *  @param rho0_b60_i Energy density at zero pressure divided by 
	 *		     $B_{60}$ [unit: ${\rm MeV\,  fm^{-3}}$]
	 *
	 */
	Eos_strange(double n0_b60_i, double b60_i, double ent0_i, 
		    double eps_fit_i, double rho0_b60_i) ;	


	Eos_strange(const Eos_strange& ) ;	/// Copy constructor	
	
    protected:
	/** Constructor from a binary file (created by the function 
	 *  {\tt sauve(FILE* )}). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  {\tt Eos::eos\_from\_file(FILE* )}. 
	 */
	Eos_strange(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}. 
	 */
	Eos_strange(ifstream& ) ; 
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ; 
	friend Eos* Eos::eos_from_file(ifstream& ) ; 

    public:
	virtual ~Eos_strange() ;			/// Destructor

    // Assignment
    // ----------
	/// Assignment to another {\tt Eos\_strange}
	void operator=(const Eos_strange& ) ;


    // Miscellaneous
    // -------------

    public : 
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ; 
    
	/** Returns a number to identify the sub-classe of {\tt Eos} the
	 *  object belongs to. 
	 */
	virtual int identify() const ; 

	/** Returns the baryon density at zero pressure divided by $B_{60}^{3/4}$
	 *   [unit: ${\rm fm}^{-3}$]
	 */
	double get_n0_b60() const {return n0_b60;} ; 
	
	/// Returns the bag constant [unit: $60\ {\rm MeV\, fm}^{-3}$]
	double get_b60() const {return b60;} ;
	
	/** Returns the log-enthalpy threshold for setting the energy density to 
	 *  a non zero value (should be negative).   
	 */
	double get_ent0() const {return ent0;} ; 	
	
	/** Returns the fitting parameter $\epsilon_{\rm fit}$ related to the
	 *  square of sound velocity by $c_s^2 = 1/3(1+\epsilon_{\rm fit})$.
	 *  [cf. Zdunik, Astron. Astrophys. {\bf 359}, 311 (2000)]
	 */
	double get_eps_fit() const {return eps_fit;} ; 
	
	/** Returns the energy density at zero pressure divided by $B_{60}$.
	 *  [unit: $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$]
	 */
	double get_rho0_b60() const {return rho0_b60;} ; 
	
    protected:
	/** Computes the auxiliary quantities {\tt n0}, {\tt rh0},
	 *  {\tt b34} and {\tt fach} from the values of the other 
	 *  parameters
	 */
	void set_auxiliary() ; 
    

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	/// Save in a file

    protected: 
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


    // Computational functions
    // -----------------------

    public: 
	/** Computes the baryon density from the log-enthalpy.
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return baryon density $n$ [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 * 
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ; 
    
 	/** Computes the total energy density from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return energy density $e$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ; 
       
 	/** Computes the pressure from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return pressure $p$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln n/d\ln H$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln e/d\ln H$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln p/d\ln H$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return dln(p)/dln(H)
	 */
    	virtual double der_press_ent_p(double ent, const Param* par=0x0) const ; 

};


		    //------------------------------------//
		    //		class Eos_strange_cr	  //
		    //------------------------------------//


/**
 * Strange matter EOS (MIT Bag model) with crust.
 *
 * For liquid core, this equation of state (EOS) corresponds to u,d,s
 * degenerated symetric
 * matter in the MIT bag model, according to approximate formula
 * given in Zdunik, Astron. Astrophys. {\bf 359}, 311 (2000).
 * The EOS for crust is a polytropic approximation of the BPS
 * model up to neutron drip point.
 *
 * @version #$Id$#
 */
class Eos_strange_cr : public Eos {

    // Data :
    // -----

    protected:
	/** Baryon density at zero pressure divided by $B_{60}^{3/4}$.
	 *   [unit: ${\rm fm}^{-3}$]
	 */
	double n0_b60 ;
	
	/// Bag constant [unit: $60\ {\rm MeV\, fm}^{-3}$]
	double b60 ;
	
	/** Log-enthalpy threshold for setting the energy density to
	 *  a non zero value (should be negative).
	 */
	double ent0 ; 	
	
	/** Fitting parameter $\epsilon_{\rm fit}$ related to the
	 *  square of sound velocity by $c_s^2 = 1/3(1+\epsilon_{\rm fit})$.
	 *  [cf. Zdunik, Astron. Astrophys. {\bf 359}, 311 (2000)]
	 */
	double eps_fit ;
	
	/** Energy density at zero pressure divided by $B_{60}$.
	 *  [unit: ${\rm MeV\,  fm^{-3}}$]
	 */
	double rho0_b60 ;
	
	/** Log-enthalpy at neutron drip point, defining the
	 *  boundary between crust and core.
	 *
	 */
	 double ent_nd ;
	
	/** Energy density at neutron drip point, defining the
	 *  boundary between crust and core
	 *  [unit: ${\rm MeV\,  fm^{-3}}$]
	 *
	 */
	 double rho_nd ;
	
	/** Adiabatic index for the crust model
	 *
	 */
	double gam ;

	// Derived data:
	// -------------

	/** Baryon density at zero pressure.
	 *   [unit: $0.1{\rm \  fm}^{-3}$ (Lorene's unit)]
	 */
	double n0 ;
	
	/** Energy density at zero pressure.
	 *  [unit: $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 *   (Lorene's unit)]
	 */
	double rho0 ;

	/** $B_{60}^{3/4}$
	 *
	 */
	double b34 ;
	
	/**  Factor $(4+\epsilon_{\rm fit})/(1+\epsilon_{\rm fit})$
	 *
	 */
	double fach ;
	
	/** Energy density at neutron drip point, defining the
	 *  boundary between crust and core
	 *  [unit: {\tt rho\_unit}]
	 *
	 */
	 double rho_nd_nucl ;
	
	/**  Ratio of pressure to energy density at neutron drip
 	 *   point
	 */
	double x_nd ;

	/// Rescaled number density at neutron drip point
	double ncr_nd ;

	/// Enthalpy shift in quark phase
	double delent ;

	/// $1/(\gamma-1)$
	double unsgam1 ;

	/// $ (\gamma - 1 -x_{\rm nd}) / \gamma / x_{\rm nd}$
	double gam1sx ;

	
    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor.
	 *
	 *  @param n0_b60_i Baryon density at zero pressure divided
	 *		    by $B_{60}^{3/4}$ [unit: ${\rm fm}^{-3}$]
	 *  @param b60_i Bag constant [unit: $60\ {\rm MeV\, fm}^{-3}$]
	 *  @param ent0_i Log-enthalpy threshold for setting the energy
	 *		    density to a non zero value (should be negative)
	 *  @param eps_fit_i Fitting parameter $\epsilon_{\rm fit}$ related
	 *		     to the square of sound velocity by
	 *		     $c_s^2 = 1/3(1+\epsilon_{\rm fit})$
	 *		[cf. Zdunik, Astron. Astrophys. {\bf 359}, 311 (2000)]
	 *  @param rho0_b60_i Energy density at zero pressure divided by
	 *		     $B_{60}$ [unit: ${\rm MeV\,  fm^{-3}}$]
	 *  @param ent_nd_i Log-enthalpy at neutron drip point,
	 *		  defining the boundary between crust and core	
	 *  @param rho_nd_i Energy density at neutron drip point,
	 *		  defining the boundary between crust and core	
	 *                [unit: ${\rm MeV\,  fm^{-3}}$]
	 *  @param gam_i Adiabatic index for the crust model
	 *
	 */
	Eos_strange_cr(double n0_b60_i, double b60_i, double ent0_i,
		    double eps_fit_i, double rho0_b60_i,
		    double ent_nd_i, double rho_nd_i,
		    double gam_i) ;	


	Eos_strange_cr(const Eos_strange_cr& ) ; /// Copy constructor	
	
    protected:
	/** Constructor from a binary file (created by the function
	 *  {\tt sauve(FILE* )}).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 *  {\tt Eos::eos\_from\_file(FILE* )}.
	 */
	Eos_strange_cr(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}.
	 */
	Eos_strange_cr(ifstream& ) ;
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_strange_cr() ;		/// Destructor

    // Assignment
    // ----------
	/// Assignment to another {\tt Eos\_strange}
	void operator=(const Eos_strange_cr& ) ;


    // Miscellaneous
    // -------------

    public :
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ;

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ;

	/** Returns a number to identify the sub-classe of {\tt Eos} the
	 *  object belongs to.
	 */
	virtual int identify() const ;

	/** Returns the baryon density at zero pressure divided by $B_{60}^{3/4}$
	 *   [unit: ${\rm fm}^{-3}$]
	 */
	double get_n0_b60() const {return n0_b60;} ;
	
	/// Returns the bag constant [unit: $60\ {\rm MeV\, fm}^{-3}$]
	double get_b60() const {return b60;} ;
	
	/** Returns the log-enthalpy threshold for setting the energy density to
	 *  a non zero value (should be negative).
	 */
	double get_ent0() const {return ent0;} ; 	
	
	/** Returns the fitting parameter $\epsilon_{\rm fit}$ related to the
	 *  square of sound velocity by $c_s^2 = 1/3(1+\epsilon_{\rm fit})$.
	 *  [cf. Zdunik, Astron. Astrophys. {\bf 359}, 311 (2000)]
	 */
	double get_eps_fit() const {return eps_fit;} ;
	
	/** Returns the energy density at zero pressure divided by $B_{60}$.
	 *  [unit: $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$]
	 */
	double get_rho0_b60() const {return rho0_b60;} ;
	
	/** Returns the log-enthalpy at neutron drip point,
	 *  defining the boundary between crust and core.
	 */
	double get_ent_nd() const {return ent_nd;} ;
	
	/** Returns the energy density at neutron drip point,
	 *  defining the boundary between crust and core
	 *  [unit: ${\rm MeV\,  fm^{-3}}$].
	 *
	 */
	double get_rho_nd() const {return rho_nd;} ;
	
	/** Returns the adiabatic index for the crust model
	 */
	double get_gam() const {return gam;} ;
	


    protected:
	/** Computes the auxiliary quantities {\tt n0}, {\tt rh0},
	 *  {\tt b34} and {\tt fach} from the values of the other
	 *  parameters
	 */
	void set_auxiliary() ;


    // Outputs
    // -------

    public:
	virtual void sauve(FILE* ) const ;	/// Save in a file

    protected:
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


    // Computational functions
    // -----------------------

    public:
	/** Computes the baryon density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$
	 *
	 *  @return baryon density $n$ [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the total energy density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$
	 *
	 *  @return energy density $e$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the pressure from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$
	 *
	 *  @return pressure $p$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative $d\ln n/d\ln H$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln e/d\ln H$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln p/d\ln H$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ 
	 *
	 *  @return dln(p)/dln(H)
	 */
    	virtual double der_press_ent_p(double ent, const Param* par=0x0) const ;

};


                //------------------------------//
                //  EOS with domain dependency //
              //------------------------------//

/**
 * EOS with domain dependency.
 *
 *
 * @version #$Id$#
 */
class MEos : public Eos {

    // Data :
    // -----

    protected:
    /// Array (upon the domains) containing the various EOS
    const Eos** mono_eos ;

    /// Number of domains
    int ndom ;

    /// Indicates wether the EOS has been constructed from a file
    bool constructed_from_file ;

    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor.
         *  @param ndom_i number of domains
         * @param mono_eos_i array (size {\tt ndom\_i}) of pointers on the various EOS
	 */
	MEos(int ndom_i, const Eos** mono_eos_i) ;

        /// Constructor for 2 domains
        MEos(const Eos& eos1, const Eos& eos2) ;

        /// Constructor for 3 domains
        MEos(const Eos& eos1, const Eos& eos2, const Eos& eos3) ;

         /// Constructor for 4 domains
       MEos(const Eos& eos1, const Eos& eos2, const Eos& eos3, const Eos& eos4) ;

	MEos(const MEos& ) ;	/// Copy constructor

    protected:
	/** Constructor from a binary file (created by the function
	 *  {\tt sauve(FILE* )}).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 *  {\tt Eos::eos\_from\_file(FILE* )}.
	 */
	MEos(FILE* ) ;

	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}.
	 */
	MEos(ifstream& ) ;

	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~MEos() ;			/// Destructor

    // Assignment
    // ----------
	/// Assignment to another {\tt MEos}
	void operator=(const MEos& ) ;


    // Miscellaneous
    // -------------

    public :
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ;

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ;

	/** Returns a number to identify the sub-classe of {\tt Eos} the
	 *  object belongs to.
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
	/** Computes the baryon density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return baryon density $n$ [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the total energy density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return energy density $e$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the pressure from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return pressure $p$ [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative $d\ln n/d\ln H$
	 * from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative $d\ln e/d\ln H$
	 * from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative $d\ln p/d\ln H$
	 * from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy $H$ defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return dln(p)/dln(H)
	 */
    	virtual double der_press_ent_p(double ent, const Param* par=0x0) const ;

};


		//------------------//
		//  Remaining EOS   //
		//------------------//
		
#include "eos_tabul.h"
		

#endif
