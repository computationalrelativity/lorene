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
 * Revision 1.13  2004/09/26 18:51:47  k_taniguchi
 * Introduction of new classes Eos_fitting, Eos_fit_SLy4, and Eos_fit_FPS
 *
 * Revision 1.12  2004/05/07 08:08:29  k_taniguchi
 * Add the case of Eos_multi_poly.C
 *
 * Revision 1.11  2004/03/22 16:10:20  j_novak
 * Excluding some files
 *
 * Revision 1.10  2004/03/22 13:12:40  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.9  2004/01/14 15:52:26  f_limousin
 * Added methods calcule, nbar_ent, der_nbar_ent and der_ener_ent for Scalar.
 *
 * Revision 1.8  2003/12/08 15:48:12  m_bejger
 * GlendNH3 (Glendenning 1985, case 3) added
 *
 * Revision 1.7  2002/10/16 14:36:28  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.6  2002/09/13 09:17:31  j_novak
 * Modif. commentaires
 *
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
#include "headcpp.h"

// Headers C
#include <stdio.h>

// Lorene classes
class Tbl ;
class Cmp ;
class Scalar ;
class Param ;

		    //------------------------------------//
		    //		base class Eos		  //
		    //------------------------------------//

/**
 * Equation of state base class. \ingroup (eos)
 * 
 *
 */
class Eos {

    // Data :
    // -----

    protected: 
	char name[100] ;	    ///< EOS name


    // Constructors - Destructor
    // -------------------------
    protected:
	Eos() ;			///< Standard constructor

	/// Standard constructor with name
	explicit Eos(const char* name_i) ; 

	Eos(const Eos& ) ;	///< Copy constructor	

    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Eos::eos_from_file(FILE*) . 
	 */
	Eos(FILE* ) ; 

	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Eos::eos_from_file(ifstream&) . 
	 */
	Eos(ifstream& ) ; 
	
	
    public:
	virtual ~Eos() ;			///< Destructor


    // Name manipulation
    // -----------------
    public:
	const char* get_name() const ;	///< Returns the EOS name

	/// Sets the EOS name
	void set_name(const char* name_i) ; 
	
    // Miscellaneous
    // -------------
    public:
	/** Construction of an EOS from a binary file.
	 *  The file must have been created by the function \c sauve(FILE*) .
	 */
	static Eos* eos_from_file(FILE* ) ; 
	
	/** Construction of an EOS from a formatted file.
	 * 
	 *  The fist line of the file must start by the EOS number, according 
	 *  to the following conventions: \\
	 *	1 = relativistic polytropic EOS (class \c Eos_poly ). \\
	 *      2 = Newtonian polytropic EOS (class \c Eos_poly_newt ). \\
	 *	3 = Relativistic incompressible EOS (class \c Eos_incomp ). \\
	 *	4 = Newtonian incompressible EOS 
	 *		    (class \c Eos_incomp_newt ). \\
	 *	5 = Strange matter (MIT Bag model) \\
	 *	6 = Strange matter (MIT Bag model) with crust \\
	 *     10 = SLy4 (Douchin \& Haensel 2001)  \\
	 *     11 = FPS (Friedman-Pandharipande + Skyrme) \\
	 *     12 = BPAL12 (Bombaci et al. 1995) \\
	 *     13 = AkmalPR (Akmal, Pandharipande \& Ravenhall 1998) \\
	 *     14 = BBB2 (Baldo, Bombaci \& Burgio 1997) \\
	 *     15 = BalbN1H1 (Balberg 2000) \\
         *     16 = GlendNH3 (Glendenning 1985, case 3) \\ 
	 *     100 = Multi-domain EOS (class \c MEos ) \\
	 *     110 = Multi-polytropic EOS (class \c Eos_multi_poly ) \\
	 *     120 = Fitted SLy4 (Shibata 2004) \\
	 *     121 = Fitted FPS (Shibata 2004) \\
	 *  The second line in the file must be the EOS name.
	 *  The following lines should contain the EOS parameters (one
	 *  parameter per line), in the same order than in the class declaration.
	 */
	static Eos* eos_from_file(ifstream& ) ; 
	
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const = 0 ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const = 0 ; 
    
	/** Returns a number to identify the sub-classe of \c Eos the
	 *  object belongs to. 
	 */
	virtual int identify() const = 0 ; 

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	///< Save in a file

	/// Display
	friend ostream& operator<<(ostream& , const Eos& ) ;	

    protected: 
	virtual ostream& operator>>(ostream &) const = 0 ;    ///< Operator >>


    // Computational functions
    // -----------------------
    protected:
	/**  General computational method for \c Cmp 's
	 *
	 *   @param thermo [input] thermodynamical quantity (for instance the
	 *	    enthalpy field)from which the
	 *          thermodynamical quantity \c resu  is to be computed.
	 *  @param nzet  [input] number of domains where \c resu  is to be
	 *	computed.
	 *  @param l_min [input] index of the innermost domain is which \c resu 
	 *	is to be computed [default value: 0]; \c resu  is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
	 *  @param fait [input] pointer on the member function of class
	 *		\c Eos which performs the pointwise calculation.
         * @param par possible extra parameters of the EOS
	 *  @param resu [output] result of the computation.
	 */
	void calcule(const Cmp& thermo, int nzet, int l_min,
		     double (Eos::*fait)(double, const Param*) const, const Param* par, Cmp& resu) const ;

	/**  General computational method for \c Scalar 's
	 *
	 *   @param thermo [input] thermodynamical quantity (for instance the
	 *	    enthalpy field)from which the
	 *          thermodynamical quantity \c resu  is to be computed.
	 *  @param nzet  [input] number of domains where \c resu  is to be
	 *	computed.
	 *  @param l_min [input] index of the innermost domain is which \c resu 
	 *	is to be computed [default value: 0]; \c resu  is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
	 *  @param fait [input] pointer on the member function of class
	 *		\c Eos which performs the pointwise calculation.
         * @param par possible extra parameters of the EOS
	 *  @param resu [output] result of the computation.
	 */
	

	void calcule(const Scalar& thermo, int nzet, int l_min,
		     double (Eos::*fait)(double, const Param*) const, const Param* par, Scalar& resu) const ;

    public:
 	/** Computes the baryon density from the log-enthalpy and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
        *  @param par possible extra parameters of the EOS
	 *
	 *  @return baryon density [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const = 0 ;

	/** Computes the baryon density field from the log-enthalpy field and
        * extra parameters
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
	 *  @param nzet  number of domains where the baryon density is to be
	 *	computed.
	 *  @param l_min  index of the innermost domain is which the baryon
	 *	density is
	 *	to be computed [default value: 0]; the baryon density is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return baryon density [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *
	 */
    	Cmp nbar_ent(const Cmp& ent, int nzet, int l_min = 0, const Param* par=0x0) const  ;

	/** Computes the baryon density field from the log-enthalpy field and
        * extra parameters
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
	 *  @param nzet  number of domains where the baryon density is to be
	 *	computed.
	 *  @param l_min  index of the innermost domain is which the baryon
	 *	density is
	 *	to be computed [default value: 0]; the baryon density is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return baryon density [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *
	 */

    	Scalar nbar_ent(const Scalar& ent, int nzet, int l_min = 0, const Param* par=0x0) const  ;

	/** Computes the total energy density from the log-enthalpy and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const = 0 ;

	/** Computes the total energy density from the log-enthalpy and extra parameters.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
	 *  @param nzet  number of domains where the energy density is to be
	 *	computed.
	 *  @param l_min  index of the innermost domain is which the energy
	 *	density is
	 *	to be computed [default value: 0]; the energy density is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return energy density [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	Cmp ener_ent(const Cmp& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;
 
 	/** Computes the total energy density from the log-enthalpy and extra parameters.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
	 *  @param nzet  number of domains where the energy density is to be
	 *	computed.
	 *  @param l_min  index of the innermost domain is which the energy
	 *	density is
	 *	to be computed [default value: 0]; the energy density is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return energy density [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */

  	Scalar ener_ent(const Scalar& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

	/** Computes the pressure from the log-enthalpy and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const = 0 ;


	/** Computes the pressure from the log-enthalpy and extra parameters
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
	 *  @param nzet  number of domains where the pressure is to be
	 *	computed.
	 *  @param l_min  index of the innermost domain is which the pressure is
	 *	to be computed [default value: 0]; the pressure is computed
	 *      only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return pressure [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 *
	 */
    	Cmp press_ent(const Cmp& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

	/** Computes the pressure from the log-enthalpy and extra parameters
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
	 *  @param nzet  number of domains where the pressure is to be
	 *	computed.
	 *  @param l_min  index of the innermost domain is which the pressure is
	 *	to be computed [default value: 0]; the pressure is computed
	 *      only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
          *  @param par possible extra parameters of the EOS
	 *
	 *  @return pressure [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 *
	 */
   
    	Scalar press_ent(const Scalar& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln n/d\ln H\f$
	 *  from the log-enthalpy and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
         *  @param par possible extra parameters of the EOS
         *
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const = 0 ;

	/** Computes the logarithmic derivative \f$d\ln n/d\ln H\f$
	 *  from the log-enthalpy and extra parameters
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
	 *  @param nzet  number of domains where the derivative
	 *	dln(n)/dln(H) is to be computed.
	 *  @param l_min  index of the innermost domain is which the
	 *	   coefficient dln(n)/dln(H) is
	 *	to be computed [default value: 0]; the derivative
	 *	dln(n)/dln(H) is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
         *  @param par possible extra parameters of the EOS
	 *
	 *  @return dln(n)/dln(H)
	 *
	 */
    	Cmp der_nbar_ent(const Cmp& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln n/d\ln H\f$
	 *  from the log-enthalpy and extra parameters
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
	 *  @param nzet  number of domains where the derivative
	 *	dln(n)/dln(H) is to be computed.
	 *  @param l_min  index of the innermost domain is which the
	 *	   coefficient dln(n)/dln(H) is
	 *	to be computed [default value: 0]; the derivative
	 *	dln(n)/dln(H) is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
         *  @param par possible extra parameters of the EOS
	 *
	 *  @return dln(n)/dln(H)
	 *
	 */
    
   	Scalar der_nbar_ent(const Scalar& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln e/d\ln H\f$
	 *  from the log-enthalpy with extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
         *  @param par possible extra parameters of the EOS
  	 *
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const = 0 ;

	/** Computes the logarithmic derivative \f$d\ln e/d\ln H\f$
	 *  from the log-enthalpy and extra parameters
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
	 *  @param nzet  number of domains where the derivative
	 *	dln(e)/dln(H) is to be computed.
	 *  @param l_min  index of the innermost domain is which the
	 *	   coefficient dln(n)/dln(H) is
	 *	to be computed [default value: 0]; the derivative
	 *	dln(e)/dln(H) is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
         *  @param par possible extra parameters of the EOS
	 *
	 *  @return dln(e)/dln(H)
	 *
	 */
    	Cmp der_ener_ent(const Cmp& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;
  
	/** Computes the logarithmic derivative \f$d\ln e/d\ln H\f$
	 *  from the log-enthalpy and extra parameters
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
	 *  @param nzet  number of domains where the derivative
	 *	dln(e)/dln(H) is to be computed.
	 *  @param l_min  index of the innermost domain is which the
	 *	   coefficient dln(n)/dln(H) is
	 *	to be computed [default value: 0]; the derivative
	 *	dln(e)/dln(H) is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
         *  @param par possible extra parameters of the EOS
	 *
	 *  @return dln(e)/dln(H)
	 *
	 */
  	
	Scalar der_ener_ent(const Scalar& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln p/d\ln H\f$
	 *  from the log-enthalpy and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
        *  @param par possible extra parameters of the EOS
	 *
	 *  @return dln(p)/dln(H)
	 */
    	virtual double der_press_ent_p(double ent, const Param* par=0x0) const = 0 ;

	/** Computes the logarithmic derivative \f$d\ln p/d\ln H\f$
	 *  from the log-enthalpy and extra parameters
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
	 *  @param nzet  number of domains where the derivative
	 *	dln(p)/dln(H) is to be computed.
        *  @param par possible extra parameters of the EOS
	 *  @param l_min  index of the innermost domain is which the
	 *	   coefficient dln(n)/dln(H) is
	 *	to be computed [default value: 0]; the derivative
	 *	dln(p)/dln(H) is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
	 *
	 *  @return dln(p)/dln(H)
	 *
	 */
    	Cmp der_press_ent(const Cmp& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

 	/** Computes the logarithmic derivative \f$d\ln p/d\ln H\f$
	 *  from the log-enthalpy and extra parameters
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
	 *  @param nzet  number of domains where the derivative
	 *	dln(p)/dln(H) is to be computed.
        *  @param par possible extra parameters of the EOS
	 *  @param l_min  index of the innermost domain is which the
	 *	   coefficient dln(n)/dln(H) is
	 *	to be computed [default value: 0]; the derivative
	 *	dln(p)/dln(H) is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
	 *
	 *  @return dln(p)/dln(H)
	 *
	 */
 
  	Scalar der_press_ent(const Scalar& ent, int nzet, int l_min = 0, const Param* par=0x0) const ;

};
ostream& operator<<(ostream& , const Eos& ) ;	


		    //------------------------------------//
		    //		class Eos_poly		  //
		    //------------------------------------//


/**
 * Polytropic equation of state (relativistic case).
 *
 * This equation of state (EOS) corresponds to identical relativistic
 * particles of rest mass is \f$m_0\f$,  whose total energy density \e e  is
 * related to their numerical density \e n  by
 * \f[ \label{eeospolye}
 *   e(n) = {\kappa \over \gamma-1} n^\gamma + \mu_0 \, n \ ,
 * \f]
 * where \f$\mu_0\f$ is the chemical potential at zero pressure.
 * The relativistic (i.e. including rest mass energy) chemical potential is
 * then
 * \f[  \label{eeospolymu}
 *   \mu(n) := {de\over dn} = {\kappa \gamma \over \gamma-1} n^{\gamma-1}
 *		+ \mu_0 \ .
 * \f]
 * The pressure is given by the (zero-temperature) First Law of Thermodynamics:
 * \f$p = \mu n - e\f$, so that
 * \f[ \label{eeospolyp}
 *   p(n) = \kappa n^\gamma  \ .
 * \f]
 * The log-enthalpy is defined as the logarithm of the ratio of the enthalpy
 * par particle by the partical rest mass energy :
 * \f[ \label{eeospolyh}
 *   H(n) := c^2 \ln \left( {e+p \over m_0 c^2\, n} \right)   \ .
 * \f]
 * According to the (zero-temperature) First Law of Thermodynamics, the
 * log-enthalpy is related to the chemical potential by
 * \f[
 *   H = c^2 \ln \left( {\mu \over m_0 c^2} \right) \ .
 * \f]
 * From this expression and relation (\ref{eeospolymu}), the expression
 * of the particle density in term of the log-enthalpy is
 * \f[
 *   n(H) = \left[ {\gamma-1\over \gamma} {m_0 c^2 \over \kappa}
                \left( \exp(H) - {\mu_0\over m_0 c^2} \right)
 *	    \right] ^{1/(\gamma-1)}  \ .
 * \f]
 * The energy density and pressure as functions of \e H  can then be obtained
 * by inserting this relation into Eqs.~(\ref{eeospolye}) and
 * (\ref{eeospolyp}).
 *\ingroup (eos)
 */
class Eos_poly : public Eos {

    // Data :
    // -----

    protected:
	/// Adiabatic index \f$\gamma\f$ (cf. Eq.~(\ref{eeospolyp}))
	double gam ;

	/** Pressure coefficient \f$\kappa\f$  (cf. Eq.~(\ref{eeospolyp}))
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$.
	 */
	double kap ; 

	/** Individual particule mass \f$m_0\f$  (cf. Eq.~(\ref{eeospolye}))
	 *  [unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
	 */
	double m_0 ;

        /** Relativistic chemical potential at zero pressure
	 *  [unit: \f$m_B c^2\f$, with \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
         * (standard value: 1)
        */
        double mu_0 ;



	double gam1 ;	    ///< \f$\gamma-1\f$
	double unsgam1 ;    ///< \f$1/(\gamma-1)\f$
	double gam1sgamkap ; ///< \f$(\gamma-1) / (\gamma \kappa) m_0\f$
        double rel_mu_0 ;       ///< \f$\mu_0/m_0\f$
        double ent_0 ;          ///< Enthalpy at zero pressure (\f$\ln (\mu_0/m_0)\f$)

    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor (sets both \c m_0 and \c mu_0  to 1).
	 *
	 *  The individual particle mass \f$m_0\f$ is set to the mean baryon
	 *  mass \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$.
	 *
	 *  @param gamma  adiabatic index \f$\gamma\f$
	 *				(cf. Eq.~(\ref{eeospolyp}))
	 *  @param kappa  pressure coefficient \f$\kappa\f$  
	 *		(cf. Eq.~(\ref{eeospolyp}))
	 *		[unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *		\f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *		\f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$
	 */
	Eos_poly(double gamma, double kappa) ;

	/** Standard constructor with individual particle mass
	*   (sets \c mu_0  to 1).
	 *  @param gamma  adiabatic index \f$\gamma\f$ (cf. Eq.~(\ref{eeospolyp}))
	 *  @param kappa  pressure coefficient \f$\kappa\f$
	 *		(cf. Eq.~(\ref{eeospolyp}))
	 *		[unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *		\f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *		\f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$
	 *  @param mass  individual particule mass \f$m_0\f$
	 *		 (cf. Eq.~(\ref{eeospolye}))
	 *		[unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$]
	 */
	Eos_poly(double gamma, double kappa, double mass) ;

	/** Standard constructor with individual particle mass and zero-pressure
         * chemical potential
	 *
	 *  @param gamma  adiabatic index \f$\gamma\f$ (cf. Eq.~(\ref{eeospolyp}))
	 *  @param kappa  pressure coefficient \f$\kappa\f$
	 *		(cf. Eq.~(\ref{eeospolyp}))
	 *		[unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *		\f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *		\f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$
	 *  @param mass  individual particule mass \f$m_0\f$
	 *		 (cf. Eq.~(\ref{eeospolye}))
	 *		[unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$]
        *  @param mu_zero  Relativistic chemical potential at zero pressure
	 *  [unit: \f$m_B c^2\f$, with \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
         * (standard value: 1)
        */
	Eos_poly(double gamma, double kappa, double mass, double mu_zero) ;

	Eos_poly(const Eos_poly& ) ;	///< Copy constructor
	
    protected:
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Eos::eos_from_file(FILE*) .
	 */
	Eos_poly(FILE* ) ;

	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c Eos::eos_from_file(ifstream&) . 
	 */
	Eos_poly(ifstream& ) ;

	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ; 

    public:
	virtual ~Eos_poly() ;			///< Destructor

    // Assignment
    // ----------
	/// Assignment to another \c Eos_poly 
	void operator=(const Eos_poly& ) ;


    // Miscellaneous
    // -------------

    public :
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ;

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ;

	/** Returns a number to identify the sub-classe of \c Eos the
	 *  object belongs to.
	 */
	virtual int identify() const ; 

	/// Returns the adiabatic index \f$\gamma\f$ (cf. Eq.~(\ref{eeospolyp}))
	double get_gam() const ;

	/** Returns the pressure coefficient \f$\kappa\f$  (cf. Eq.~(\ref{eeospolyp}))
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$.
	 */
	double get_kap() const ;
	
	/** Return the individual particule mass \f$m_0\f$
	 *  (cf. Eq.~(\ref{eeospolye}))
	 *  [unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
	 */
	double get_m_0() const ;

	/** Return the relativistic chemical potential at zero pressure
	 *  [unit: \f$m_B c^2\f$, with \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
	 */
	double get_mu_0() const ;

    protected:
	/** Computes the auxiliary quantities \c gam1 , \c unsgam1 ,
	 *  \c gam1sgamkap  from the values of \c gam  and \c kap 
	 */
	void set_auxiliary() ;


    // Outputs
    // -------

    public:
	virtual void sauve(FILE* ) const ;	///< Save in a file

    protected:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>


    // Computational functions
    // -----------------------

    public:
	/** Computes the baryon density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return baryon density \e n  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the total energy density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the pressure from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln n/d\ln H\f$
	 * from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln e/d\ln H\f$
	 * from the log-enthalpy.
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln p/d\ln H\f$
	 * from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
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
 * Polytropic equation of state (Newtonian case). \ingroup (eos)
 *
 * This equation of state (EOS) corresponds to identical non relativistic
 * particles of rest mass is \f$m_0\f$,  whose internal energy density \f$\epsilon\f$ is
 * related to their numerical density \e n  by
 * \f[ \label{eeospolynewte}
 *   \epsilon(n) = {\kappa \over \gamma-1} n^\gamma  \ .
 * \f]
 * The (non-relativistic) chemical potential is
 * then
 * \f[  \label{eeospolynewtmu}
 *   \mu(n) := {d\epsilon\over dn} = {\kappa \gamma \over \gamma-1} n^{\gamma-1}
 *		  \ .
 * \f]
 * The pressure is given by the (zero-temperature) First Law of Thermodynamics:
 * \f$p = \mu n - \epsilon\f$, so that
 * \f[ \label{eeospolynewtp}
 *   p(n) = \kappa n^\gamma  \ .  
 * \f]
 * The (non-relativistic) specific enthalpy is  :
 * \f[ \label{eeospolynewth}
 *   h(n) := {\epsilon + p \over m_0 n}   \ .  
 * \f]
 * According to the (zero-temperature) First Law of Thermodynamics, the
 * specific enthalpy is related to the chemical potential by
 * \f[
 *   h = {\mu \over m_0}  \ . 
 * \f]
 * From this expression and relation (\ref{eeospolynewtmu}), the expression
 * of the particle density in term of the  specific enthalpy is 
 * \f[
 *   n(h) = \left[ {\gamma-1\over \gamma} {m_0 \over \kappa} h
 *	    \right] ^{1/(\gamma-1)}  \ .  
 * \f]
 * The energy density and pressure as functions of \e H  can then be obtained
 * by inserting this relation into Eqs.~(\ref{eeospolynewte}) and 
 * (\ref{eeospolynewtp}). 
 *
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
	 *  The individual particle mass \f$m_0\f$ is set to the mean baryon
	 *  mass \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$. 
	 *  
	 *  @param gamma  adiabatic index \f$\gamma\f$ 
	 *				(cf. Eq.~(\ref{eeospolynewtp}))
	 *  @param kappa  pressure coefficient \f$\kappa\f$  
	 *		(cf. Eq.~(\ref{eeospolynewtp}))
	 *		[unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *		\f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *		\f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$
	 */
	Eos_poly_newt(double gamma, double kappa) ;	

	Eos_poly_newt(const Eos_poly_newt& ) ;	///< Copy constructor	

    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Eos::eos_from_file(FILE*) . 
	 */
	Eos_poly_newt(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Eos::eos_from_file(ifstream&) . 
	 */
	Eos_poly_newt(ifstream& ) ; 
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ; 
	friend Eos* Eos::eos_from_file(ifstream& ) ; 


    public:
	virtual ~Eos_poly_newt() ;			///< Destructor

    // Assignment
    // ----------
	/// Assignment to another \c Eos_poly_newt 
	void operator=(const Eos_poly_newt& ) ;

    // Miscellaneous
    // -------------

    public : 
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ; 
    
	/** Returns a number to identify the sub-classe of \c Eos the
	 *  object belongs to. 
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
	/** Computes the baryon density from the specific enthalpy.
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] specific enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolynewth})
	 *
	 *  @return baryon density \e n  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 * 
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ; 
    
 	/** Computes the total energy density from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] specific enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolynewth})
	 *
	 *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ; 
       
 	/** Computes the pressure from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] specific enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolynewth})
	 *
	 *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln n/d\ln h\f$ 
	 * from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] specific enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolynewth})
	 *
	 *  @return dln(n)/dln(h)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln e/d\ln h\f$ 
	 * from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] specific enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolynewth})
	 *
	 *  @return dln(e)/dln(h)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln p/d\ln h\f$ 
	 * from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] specific enthalpy \e H  defined by
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
 * matter \f$e = \rho_0 \,  c^2\f$. 
 * \ingroup (eos)
 *
 */
class Eos_incomp : public Eos {

    // Data :
    // -----

    protected: 
	/// Constant density \f$\rho_0\f$ 
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
	 *  @param rho_c  constant density \f$\rho_0\f$ 
	 *		    [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *		    \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
	Eos_incomp(double rho_c) ;	

	/** Standard constructor with log-enthalpy threshold. 
	 * 
	 *  @param rho_c  constant density \f$\rho_0\f$
	 *		    [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *		    \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 *  @param ent_c  log-enthalpy threshold for non-zero density 
	 *		  [unit: \f$c^2\f$] :
	 *		  the energy density is set to \c rho_c} for
	 *		  \f$H > {\tt ent_c}\f$. 
	 */
	Eos_incomp(double rho_c, double ent_c) ;	

	Eos_incomp(const Eos_incomp& ) ;	///< Copy constructor	

    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Eos::eos_from_file(FILE*) . 
	 */
	Eos_incomp(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Eos::eos_from_file(ifstream&) . 
	 */
	Eos_incomp(ifstream& ) ; 
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ; 
	friend Eos* Eos::eos_from_file(ifstream& ) ; 

    public:
	virtual ~Eos_incomp() ;			///< Destructor

    // Assignment
    // ----------
	/// Assignment to another \c Eos_incomp 
	void operator=(const Eos_incomp& ) ;
    

    // Miscellaneous
    // -------------

    public : 
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ; 
    
	/** Returns a number to identify the sub-classe of \c Eos the
	 *  object belongs to. 
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
	/** Computes the baryon density from the log-enthalpy.
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return baryon density \e n  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 * 
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ; 
    
 	/** Computes the total energy density from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ; 
       
 	/** Computes the pressure from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln n/d\ln H\f$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln e/d\ln H\f$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln p/d\ln H\f$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
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
 * matter \f$\rho = \rho_0\f$. 
 * \ingroup (eos)
 *
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
	 *  @param rho_c  constant density \f$\rho_0\f$ 
	 *		    [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *		    \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
	Eos_incomp_newt(double rho_c) ;	

	/** Standard constructor with specific enthalpy threshold. 
	 * 
	 *  @param rho_c  constant density \f$\rho_0\f$
	 *		    [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *		    \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 *  @param ent_c  specific enthalpy threshold for non-zero density 
	 *		  [unit: \f$c^2\f$] :
	 *		  the energy density is set to \c rho_c  for
	 *		  \f$h > {\tt ent_c}\f$. 
	 */
	Eos_incomp_newt(double rho_c, double ent_c) ;	

	Eos_incomp_newt(const Eos_incomp_newt& ) ;	///< Copy constructor	

    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Eos::eos_from_file(FILE*) . 
	 */
	Eos_incomp_newt(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Eos::eos_from_file(ifstream&) . 
	 */
	Eos_incomp_newt(ifstream& ) ; 
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ; 
	friend Eos* Eos::eos_from_file(ifstream& ) ; 

    public:
	virtual ~Eos_incomp_newt() ;			///< Destructor

    // Assignment
    // ----------
	/// Assignment to another \c Eos_incomp_newt 
	void operator=(const Eos_incomp_newt& ) ;
    

    // Miscellaneous
    // -------------

    public : 
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ; 
    
	/** Returns a number to identify the sub-classe of \c Eos the
	 *  object belongs to. 
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
	/** Computes the baryon density from the specific enthalpy.
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] specific enthalpy \e H 
	 *
	 *  @return baryon density \e n  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 * 
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ; 
    
 	/** Computes the total energy density from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] specific enthalpy \e H  
	 *
	 *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ; 
       
 	/** Computes the pressure from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] specific enthalpy \e H  
	 *
	 *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln n/d\ln h\f$ 
	 * from the specific enthalpy.
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] specific enthalpy \e H  
	 *
	 *  @return dln(n)/dln(h)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln e/d\ln h\f$ 
	 * from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] specific enthalpy \e H  
	 *
	 *  @return dln(e)/dln(h)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln p/d\ln h\f$ 
	 * from the specific enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] specific enthalpy \e H  
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
 * given in Zdunik, Astron. Astrophys. \b 359 , 311 (2000). 
 *\ingroup (eos)
 */
class Eos_strange : public Eos {

    // Data :
    // -----

    protected: 
	/** Baryon density at zero pressure divided by \f$B_{60}^{3/4}\f$.
	 *   [unit: \f${\rm fm}^{-3}\f$]
	 */
	double n0_b60 ; 
	
	/// Bag constant [unit: \f$60\ {\rm MeV\, fm}^{-3}\f$]
	double b60 ;
	
	/** Log-enthalpy threshold for setting the energy density to 
	 *  a non zero value (should be negative).   
	 */
	double ent0 ; 	
	
	/** Fitting parameter \f$\epsilon_{\rm fit}\f$ related to the
	 *  square of sound velocity by \f$c_s^2 = 1/3(1+\epsilon_{\rm fit})\f$.
	 *  [cf. Zdunik, Astron. Astrophys. \b 359 , 311 (2000)]
	 */
	double eps_fit ; 
	
	/** Energy density at zero pressure divided by \f$B_{60}\f$.
	 *  [unit: \f${\rm MeV\,  fm^{-3}}\f$]
	 */
	double rho0_b60 ; 
	
	/** Baryon density at zero pressure.
	 *   [unit: \f$0.1{\rm \  fm}^{-3}\f$ (Lorene's unit)]
	 */
	double n0 ; 
	
	/** Energy density at zero pressure. 
	 *  [unit: \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ 
	 *   (Lorene's unit)]
	 */
	double rho0 ; 

	/** \f$B_{60}^{3/4}\f$
	 * 
	 */
	double b34 ;
	
	/**  Factor \f$(4+\epsilon_{\rm fit})/(1+\epsilon_{\rm fit})\f$
	 * 
	 */
	double fach ; 

    // Constructors - Destructor
    // -------------------------
    public:
    
	/** Standard constructor.
	 * 
	 *  @param n0_b60_i Baryon density at zero pressure divided 
	 *		    by \f$B_{60}^{3/4}\f$ [unit: \f${\rm fm}^{-3}\f$]
	 *  @param b60_i Bag constant [unit: \f$60\ {\rm MeV\, fm}^{-3}\f$]
	 *  @param ent0_i Log-enthalpy threshold for setting the energy 
	 *		    density to a non zero value (should be negative)
	 *  @param eps_fit_i Fitting parameter \f$\epsilon_{\rm fit}\f$ related 
	 *		     to the square of sound velocity by 
	 *		     \f$c_s^2 = 1/3(1+\epsilon_{\rm fit})\f$
	 *		[cf. Zdunik, Astron. Astrophys. \b 359 , 311 (2000)]
	 *  @param rho0_b60_i Energy density at zero pressure divided by 
	 *		     \f$B_{60}\f$ [unit: \f${\rm MeV\,  fm^{-3}}\f$]
	 *
	 */
	Eos_strange(double n0_b60_i, double b60_i, double ent0_i, 
		    double eps_fit_i, double rho0_b60_i) ;	


	Eos_strange(const Eos_strange& ) ;	///< Copy constructor	
	
    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Eos::eos_from_file(FILE*) . 
	 */
	Eos_strange(FILE* ) ; 
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Eos::eos_from_file(ifstream&) . 
	 */
	Eos_strange(ifstream& ) ; 
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ; 
	friend Eos* Eos::eos_from_file(ifstream& ) ; 

    public:
	virtual ~Eos_strange() ;			///< Destructor

    // Assignment
    // ----------
	/// Assignment to another \c Eos_strange 
	void operator=(const Eos_strange& ) ;


    // Miscellaneous
    // -------------

    public : 
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ; 
    
	/** Returns a number to identify the sub-classe of \c Eos the
	 *  object belongs to. 
	 */
	virtual int identify() const ; 

	/** Returns the baryon density at zero pressure divided by \f$B_{60}^{3/4}\f$
	 *   [unit: \f${\rm fm}^{-3}\f$]
	 */
	double get_n0_b60() const {return n0_b60;} ; 
	
	/// Returns the bag constant [unit: \f$60\ {\rm MeV\, fm}^{-3}\f$]
	double get_b60() const {return b60;} ;
	
	/** Returns the log-enthalpy threshold for setting the energy density to 
	 *  a non zero value (should be negative).   
	 */
	double get_ent0() const {return ent0;} ; 	
	
	/** Returns the fitting parameter \f$\epsilon_{\rm fit}\f$ related to the
	 *  square of sound velocity by \f$c_s^2 = 1/3(1+\epsilon_{\rm fit})\f$.
	 *  [cf. Zdunik, Astron. Astrophys. \b 359 , 311 (2000)]
	 */
	double get_eps_fit() const {return eps_fit;} ; 
	
	/** Returns the energy density at zero pressure divided by \f$B_{60}\f$.
	 *  [unit: \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$]
	 */
	double get_rho0_b60() const {return rho0_b60;} ; 
	
    protected:
	/** Computes the auxiliary quantities \c n0 , \c rh0 ,
	 *  \c b34  and \c fach  from the values of the other 
	 *  parameters
	 */
	void set_auxiliary() ; 
    

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	///< Save in a file

    protected: 
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>


    // Computational functions
    // -----------------------

    public: 
	/** Computes the baryon density from the log-enthalpy.
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return baryon density \e n  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 * 
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ; 
    
 	/** Computes the total energy density from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ; 
       
 	/** Computes the pressure from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln n/d\ln H\f$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln e/d\ln H\f$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln p/d\ln H\f$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
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
 * given in Zdunik, Astron. Astrophys. \b 359 , 311 (2000).
 * The EOS for crust is a polytropic approximation of the BPS
 * model up to neutron drip point.
 *\ingroup (eos)
 */
class Eos_strange_cr : public Eos {

    // Data :
    // -----

    protected:
	/** Baryon density at zero pressure divided by \f$B_{60}^{3/4}\f$.
	 *   [unit: \f${\rm fm}^{-3}\f$]
	 */
	double n0_b60 ;
	
	/// Bag constant [unit: \f$60\ {\rm MeV\, fm}^{-3}\f$]
	double b60 ;
	
	/** Log-enthalpy threshold for setting the energy density to
	 *  a non zero value (should be negative).
	 */
	double ent0 ; 	
	
	/** Fitting parameter \f$\epsilon_{\rm fit}\f$ related to the
	 *  square of sound velocity by \f$c_s^2 = 1/3(1+\epsilon_{\rm fit})\f$.
	 *  [cf. Zdunik, Astron. Astrophys. \b 359 , 311 (2000)]
	 */
	double eps_fit ;
	
	/** Energy density at zero pressure divided by \f$B_{60}\f$.
	 *  [unit: \f${\rm MeV\,  fm^{-3}}\f$]
	 */
	double rho0_b60 ;
	
	/** Log-enthalpy at neutron drip point, defining the
	 *  boundary between crust and core.
	 *
	 */
	 double ent_nd ;
	
	/** Energy density at neutron drip point, defining the
	 *  boundary between crust and core
	 *  [unit: \f${\rm MeV\,  fm^{-3}}\f$]
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
	 *   [unit: \f$0.1{\rm \  fm}^{-3}\f$ (Lorene's unit)]
	 */
	double n0 ;
	
	/** Energy density at zero pressure.
	 *  [unit: \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 *   (Lorene's unit)]
	 */
	double rho0 ;

	/** \f$B_{60}^{3/4}\f$
	 *
	 */
	double b34 ;
	
	/**  Factor \f$(4+\epsilon_{\rm fit})/(1+\epsilon_{\rm fit})\f$
	 *
	 */
	double fach ;
	
	/** Energy density at neutron drip point, defining the
	 *  boundary between crust and core
	 *  [unit: \c rho_unit ]
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

	/// \f$1/(\gamma-1)\f$
	double unsgam1 ;

	/// \f$ (\gamma - 1 -x_{\rm nd}) / \gamma / x_{\rm nd}\f$
	double gam1sx ;

	
    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor.
	 *
	 *  @param n0_b60_i Baryon density at zero pressure divided
	 *		    by \f$B_{60}^{3/4}\f$ [unit: \f${\rm fm}^{-3}\f$]
	 *  @param b60_i Bag constant [unit: \f$60\ {\rm MeV\, fm}^{-3}\f$]
	 *  @param ent0_i Log-enthalpy threshold for setting the energy
	 *		    density to a non zero value (should be negative)
	 *  @param eps_fit_i Fitting parameter \f$\epsilon_{\rm fit}\f$ related
	 *		     to the square of sound velocity by
	 *		     \f$c_s^2 = 1/3(1+\epsilon_{\rm fit})\f$
	 *		[cf. Zdunik, Astron. Astrophys. \b 359 , 311 (2000)]
	 *  @param rho0_b60_i Energy density at zero pressure divided by
	 *		     \f$B_{60}\f$ [unit: \f${\rm MeV\,  fm^{-3}}\f$]
	 *  @param ent_nd_i Log-enthalpy at neutron drip point,
	 *		  defining the boundary between crust and core	
	 *  @param rho_nd_i Energy density at neutron drip point,
	 *		  defining the boundary between crust and core	
	 *                [unit: \f${\rm MeV\,  fm^{-3}}\f$]
	 *  @param gam_i Adiabatic index for the crust model
	 *
	 */
	Eos_strange_cr(double n0_b60_i, double b60_i, double ent0_i,
		    double eps_fit_i, double rho0_b60_i,
		    double ent_nd_i, double rho_nd_i,
		    double gam_i) ;	


	Eos_strange_cr(const Eos_strange_cr& ) ; ///< Copy constructor	
	
    protected:
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 *  \c Eos::eos_from_file(FILE*) .
	 */
	Eos_strange_cr(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c Eos::eos_from_file(ifstream&) .
	 */
	Eos_strange_cr(ifstream& ) ;
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_strange_cr() ;		///< Destructor

    // Assignment
    // ----------
	/// Assignment to another \c Eos_strange 
	void operator=(const Eos_strange_cr& ) ;


    // Miscellaneous
    // -------------

    public :
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ;

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ;

	/** Returns a number to identify the sub-classe of \c Eos the
	 *  object belongs to.
	 */
	virtual int identify() const ;

	/** Returns the baryon density at zero pressure divided by \f$B_{60}^{3/4}\f$
	 *   [unit: \f${\rm fm}^{-3}\f$]
	 */
	double get_n0_b60() const {return n0_b60;} ;
	
	/// Returns the bag constant [unit: \f$60\ {\rm MeV\, fm}^{-3}\f$]
	double get_b60() const {return b60;} ;
	
	/** Returns the log-enthalpy threshold for setting the energy density to
	 *  a non zero value (should be negative).
	 */
	double get_ent0() const {return ent0;} ; 	
	
	/** Returns the fitting parameter \f$\epsilon_{\rm fit}\f$ related to the
	 *  square of sound velocity by \f$c_s^2 = 1/3(1+\epsilon_{\rm fit})\f$.
	 *  [cf. Zdunik, Astron. Astrophys. \b 359 , 311 (2000)]
	 */
	double get_eps_fit() const {return eps_fit;} ;
	
	/** Returns the energy density at zero pressure divided by \f$B_{60}\f$.
	 *  [unit: \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$]
	 */
	double get_rho0_b60() const {return rho0_b60;} ;
	
	/** Returns the log-enthalpy at neutron drip point,
	 *  defining the boundary between crust and core.
	 */
	double get_ent_nd() const {return ent_nd;} ;
	
	/** Returns the energy density at neutron drip point,
	 *  defining the boundary between crust and core
	 *  [unit: \f${\rm MeV\,  fm^{-3}}\f$].
	 *
	 */
	double get_rho_nd() const {return rho_nd;} ;
	
	/** Returns the adiabatic index for the crust model
	 */
	double get_gam() const {return gam;} ;
	


    protected:
	/** Computes the auxiliary quantities \c n0 , \c rh0 ,
	 *  \c b34  and \c fach  from the values of the other
	 *  parameters
	 */
	void set_auxiliary() ;


    // Outputs
    // -------

    public:
	virtual void sauve(FILE* ) const ;	///< Save in a file

    protected:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>


    // Computational functions
    // -----------------------

    public:
	/** Computes the baryon density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H 
	 *
	 *  @return baryon density \e n  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the total energy density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H 
	 *
	 *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the pressure from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H 
	 *
	 *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln n/d\ln H\f$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln e/d\ln H\f$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative \f$d\ln p/d\ln H\f$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
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
         * @param mono_eos_i array (size \c ndom_i ) of pointers on the various EOS
	 */
	MEos(int ndom_i, const Eos** mono_eos_i) ;

        /// Constructor for 2 domains
        MEos(const Eos& eos1, const Eos& eos2) ;

        /// Constructor for 3 domains
        MEos(const Eos& eos1, const Eos& eos2, const Eos& eos3) ;

         /// Constructor for 4 domains
       MEos(const Eos& eos1, const Eos& eos2, const Eos& eos3, const Eos& eos4) ;

	MEos(const MEos& ) ;	///< Copy constructor

    protected:
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 *  \c Eos::eos_from_file(FILE*) .
	 */
	MEos(FILE* ) ;

	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c Eos::eos_from_file(ifstream&) .
	 */
	MEos(ifstream& ) ;

	///< The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~MEos() ;			///< Destructor

    // Assignment
    // ----------
	/// Assignment to another \c MEos 
	void operator=(const MEos& ) ;


    // Miscellaneous
    // -------------

    public :
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ;

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ;

	/** Returns a number to identify the sub-classe of \c Eos the
	 *  object belongs to.
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
	/** Computes the baryon density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return baryon density \e n  [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the total energy density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the pressure from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln n/d\ln H\f$
	 * from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln e/d\ln H\f$
	 * from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
	 *				     Eq. (\ref{eeospolyh})
	 *
        *  @param par possible extra parameters of the EOS
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln p/d\ln H\f$
	 * from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
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
