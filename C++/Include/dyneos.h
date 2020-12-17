/*
 *  Definition of Lorene classes Dyn_eos
 *				 Dyn_eos_poly
 *                               Dyn_eos_tab
 *                               Dyn_eos_cons
 *
 */

/*
 *   Copyright (c) 2019 Jerome Novak
 *             (c) 2000 Eric Gourgoulhon for Eos classes
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


#ifndef __DYNEOS_H_ 
#define __DYNEOS_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.2  2020/12/17 17:00:27  j_novak
 * Output of sound speed squared, instead of sound speed.
 *
 * Revision 1.1  2019/12/06 14:30:50  j_novak
 * New classes Dyn_eos... for cold Eos's with baryon density as input.
 *
 *
 *
 * $Header$
 *
 */

// Standard C++
#include "headcpp.h"

// Lorene classes
namespace Lorene {
  class Tbl ;
  class Param ;
  class Scalar ;

		    //------------------------------------//
		    //	     base class Dyn_eos		  //
		    //------------------------------------//

/**
 * Equation of state for use in dynamical code base class. \ingroup(eos)
 * 
 * This class describes equation of states with baryon density as parameter
 * as opposed to the classes of the group \c Eos where it is the 
 * log-enthalpy which is used.
 *
 */
class Dyn_eos {

    // Data :
    // -----

    protected: 
	string name ;	    ///< EOS name


    // Constructors - Destructor
    // -------------------------
    protected:
	Dyn_eos() ;			///< Standard constructor

	/// Standard constructor with name
	explicit Dyn_eos(const string&) ; 

	Dyn_eos(const Dyn_eos& ) ;	///< Copy constructor	

    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Dyn_eos::eos_from_file(FILE*) . 
	 */
	Dyn_eos(FILE* ) ; 

	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Dyn_eos::eos_from_file(ifstream&) . 
	 */
	Dyn_eos(ifstream& ) ; 
	
	
    public:
	virtual ~Dyn_eos() ;			///< Destructor


    // Name manipulation
    // -----------------
    public:
	const string& get_name() const ;	///< Returns the EOS name

	/// Sets the EOS name
	void set_name(const string&) ; 
	
    // Miscellaneous
    // -------------
    public:
	/** Construction of an EOS from a binary file.
	 *  The file must have been created by the function \c sauve(FILE*) .
	 */
	static Dyn_eos* eos_from_file(FILE* ) ; 
	
	/** Construction of an EOS from a formatted file.
	 * 
	 *  The fist line of the file must start by the EOS number, according 
	 *  to the following conventions (same as fo the classes \c Eos ):
	 *  - 1 = relativistic polytropic EOS (class \c Dyn_eos_poly ). 
	 *  - 2 = Newtonian polytropic EOS (class \c Dyn_eos_poly_newt ). 
	 *  - 17 = Tabulated EOS (class \c Dyn_eos_tab ).
	 *  - 20 = Consistent EOS from table (class \c Dyn_eos_cons ).
	 *
	 *  The second line in the file should contain a name given by the user to the EOS.
	 *  The following lines should contain the EOS parameters (one
	 *  parameter per line), in the same order than in the class declaration.
	 */
	static Dyn_eos* eos_from_file(ifstream& ) ; 
	
	/// Comparison operator (egality)
	virtual bool operator==(const Dyn_eos& ) const = 0 ; 

	/// Comparison operator (difference)
	virtual bool operator!=(const Dyn_eos& ) const = 0 ; 
    
	/** Returns a number to identify the sub-classe of \c Dyn_eos the
	 *  object belongs to. 
	 */
	virtual int identify() const = 0 ; 

    // Outputs
    // -------

    public: 
	virtual void sauve(FILE* ) const ;	///< Save in a file

	/// Display
	friend ostream& operator<<(ostream& , const Dyn_eos& ) ;	

    protected: 
	virtual ostream& operator>>(ostream &) const = 0 ;    ///< Operator >>


    // Computational functions
    // -----------------------
    protected:
	/**  General computational method for \c Scalar 's
	 *
	 *   @param thermo [input] thermodynamical quantity (for instance the
	 *	    density field) from which the
	 *          thermodynamical quantity \c resu  is to be computed.
	 *  @param nzet  [input] number of domains where \c resu  is to be
	 *	computed.
	 *  @param l_min [input] index of the innermost domain is which \c resu 
	 *	is to be computed [default value: 0]; \c resu  is
	 *	computed only in domains whose indices are in
	 *      \c [l_min,l_min+nzet-1] . In the other
	 *	domains, it is set to zero.
	 *  @param fait [input] pointer on the member function of class
	 *		\c Dyn_eos which performs the pointwise calculation.
         * @param par possible extra parameters of the EOS
	 *  @param resu [output] result of the computation.
	 */
	void calcule(const Scalar& thermo, int nzet, int l_min,
		     double (Dyn_eos::*fait)(double, const Param*) const,
		     Param* par, Scalar& resu) const ;

    public:
 	/** Computes the log-enthalpy from the baryon density and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *         baryon density
         *  @param par possible extra parameters of the EOS
	 *
	 *  @return ent log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass.
	 *
	 */
    	virtual double ent_nbar_p(double nbar, const Param* par=0x0) const = 0 ;

	/** Computes the log-enthalpy field from the baryon density field and 
	 *  extra parameters
	 *
	 *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *         baryon density
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
	 *  @return ent log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass.
	 *
	 */
    	Scalar ent_nbar(const Scalar& nbar, int nzet, int l_min = 0, Param* par=0x0) const  ;

	/** Computes the total energy density from the baryon density and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *         baryon density 
         *  @param par possible extra parameters of the EOS
	 *
	 *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_nbar_p(double nbar, const Param* par=0x0) const = 0 ;

	/** Computes the total energy density from the baryon density and extra parameters.
	 *
	 *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *         baryon density 
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
  	Scalar ener_nbar(const Scalar& nbar, int nzet, int l_min = 0, Param* par=0x0) const ;

	/** Computes the pressure from the baryon density and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *         baryon density 
	 *  @param par possible extra parameters of the EOS
	 *
	 *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_nbar_p(double nbar, const Param* par=0x0) const = 0 ;


	/** Computes the pressure from the baryon density and extra parameters
	 *
	 *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *         baryon density 
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
    	Scalar press_nbar(const Scalar& nbar, int nzet, int l_min = 0, Param* par=0x0) const ;

	/** Computes the sound speed squared \f$ c_s^2 = c^2 \frac{dp}{de}\f$
	 *  from the baryon density with extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *         baryon density 
         *  @param par possible extra parameters of the EOS
  	 *
	 *  @return \f$c_s^2 \f$ [unit: \e c^2]
	 */
    	virtual double csound_square_nbar_p(double nbar, const Param* par=0x0) const = 0 ;

	/** Computes the sound speed squared \f$ c_s^2 = c^2 \frac{dp}{de}\f$
	 *  from the baryon density with extra parameters
	 *
	 *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *         baryon density 
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
	 *  @return \f$c_s^2 \f$ [unit: \e c^2]
	 *
	 */
	Scalar csound_square_nbar(const Scalar& nbar, int nzet, int l_min = 0, Param* par=0x0) const ;


};
ostream& operator<<(ostream& , const Dyn_eos& ) ;	


		    //------------------------------------//
		    //	       class Dyn_eos_poly         //
		    //------------------------------------//


/**
 * Polytropic equation of state (relativistic case) for use in dynamical code.
 *
 * This equation of state (EOS) corresponds to identical relativistic
 * particles of rest mass is \f$m_0\f$,  whose total energy density \e e  is
 * related to their numerical density \e n  by
 * \f[ 
 *   e(n) = {\kappa \over \gamma-1} n^\gamma + \mu_0 \, n \ , \qquad \qquad (1)
 * \f]
 * where \f$\mu_0\f$ is the chemical potential at zero pressure.
 * The relativistic (i.e. including rest mass energy) chemical potential is
 * then
 * \f[  
 *   \mu(n) := {de\over dn} = {\kappa \gamma \over \gamma-1} n^{\gamma-1}
 *		+ \mu_0 \ .\qquad \qquad (2)
 * \f]
 * The pressure is given by the (zero-temperature) First Law of Thermodynamics:
 * \f$p = \mu n - e\f$, so that
 * \f[ 
 *   p(n) = \kappa n^\gamma  \ . \qquad \qquad (3)
 * \f]
 * The log-enthalpy is defined as the logarithm of the ratio of the enthalpy
 * par particle by the partical rest mass energy :
 * \f[ 
 *   H(n) := c^2 \ln \left( {e+p \over m_0 c^2\, n} \right)   \ . \qquad \qquad (4)
 * \f]
 * According to the (zero-temperature) First Law of Thermodynamics, the
 * log-enthalpy is related to the chemical potential by
 * \f[
 *   H = c^2 \ln \left( {\mu \over m_0 c^2} \right) \ .  \qquad \qquad (5)
 * \f]
 *
 *\ingroup (eos)
 */
class Dyn_eos_poly : public Dyn_eos {

    // Data :
    // -----

    protected:
	/// Adiabatic index \f$\gamma\f$ (cf. Eq. (3))
	double gam ;

	/** Pressure coefficient \f$\kappa\f$  (cf. Eq. (3))
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$.
	 */
	double kap ; 

	/** Individual particule mass \f$m_0\f$  (cf. Eq. (1))
	 *  [unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
	 */
	double m_0 ;

        /** Relativistic chemical potential at zero pressure
	 *  [unit: \f$m_B c^2\f$, with \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
         * (standard value: 1)
        */
        double mu_0 ;



	double gam1 ;	    ///< \f$\gamma-1\f$
	double kapsgam1 ;    ///< \f$\kappa/(\gamma-1)\f$
	double gamkapsgam1 ; ///< \f$(\gamma \kappa) / [(\gamma - 1)*m_0]\f$
        double rel_mu_0 ;       ///< \f$\mu_0/m_0\f$

    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor (sets both \c m_0 and \c mu_0  to 1).
	 *
	 *  The individual particle mass \f$m_0\f$ is set to the mean baryon
	 *  mass \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$.
	 *
	 *  @param gamma  adiabatic index \f$\gamma\f$
	 *				(cf. Eq. (3))
	 *  @param kappa  pressure coefficient \f$\kappa\f$  
	 *		(cf. Eq. (3))
	 *		[unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *		\f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *		\f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$
	 */
	Dyn_eos_poly(double gamma, double kappa) ;

	/** Standard constructor with individual particle mass
	*   (sets \c mu_0  to 1).
	 *  @param gamma  adiabatic index \f$\gamma\f$ (cf. Eq. (3))
	 *  @param kappa  pressure coefficient \f$\kappa\f$
	 *		(cf. Eq. (3))
	 *		[unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *		\f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *		\f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$
	 *  @param mass  individual particule mass \f$m_0\f$
	 *		 (cf. Eq. (1)
	 *		[unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$]
	 */
	Dyn_eos_poly(double gamma, double kappa, double mass) ;

	/** Standard constructor with individual particle mass and zero-pressure
         * chemical potential
	 *
	 *  @param gamma  adiabatic index \f$\gamma\f$ (cf. Eq. (3))
	 *  @param kappa  pressure coefficient \f$\kappa\f$
	 *		(cf. Eq. (3))
	 *		[unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *		\f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *		\f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$
	 *  @param mass  individual particule mass \f$m_0\f$
	 *		 (cf. Eq. (1))
	 *		[unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$]
        *  @param mu_zero  Relativistic chemical potential at zero pressure
	 *  [unit: \f$m_B c^2\f$, with \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
         * (standard value: 1)
        */
	Dyn_eos_poly(double gamma, double kappa, double mass, double mu_zero) ;

	Dyn_eos_poly(const Dyn_eos_poly& ) ;	///< Copy constructor
	
    protected:
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Dyn_eos::eos_from_file(FILE*) .
	 */
	Dyn_eos_poly(FILE* ) ;

	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c Dyn_eos::eos_from_file(ifstream&) . 
	 */
	Dyn_eos_poly(ifstream& ) ;

	/// The construction functions from a file
	friend Dyn_eos* Dyn_eos::eos_from_file(FILE* ) ;
	friend Dyn_eos* Dyn_eos::eos_from_file(ifstream& ) ; 

    public:
	virtual ~Dyn_eos_poly() ;			///< Destructor

    // Assignment
    // ----------
	/// Assignment to another \c Dyn_eos_poly 
	void operator=(const Dyn_eos_poly& ) ;


    // Miscellaneous
    // -------------

    public :
	/// Comparison operator (egality)
	virtual bool operator==(const Dyn_eos& ) const ;

	/// Comparison operator (difference)
	virtual bool operator!=(const Dyn_eos& ) const ;

	/** Returns a number to identify the sub-classe of \c Dyn_eos the
	 *  object belongs to.
	 */
	virtual int identify() const ; 

	/// Returns the adiabatic index \f$\gamma\f$ (cf. Eq. (3))
	double get_gam() const ;

	/** Returns the pressure coefficient \f$\kappa\f$  (cf. Eq. (3))
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$.
	 */
	double get_kap() const ;
	
	/** Return the individual particule mass \f$m_0\f$
	 *  (cf. Eq. (1))
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
 	/** Computes the log-enthalpy from the baryon density and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *         baryon density
         *  @param par possible extra parameters of the EOS
	 *
	 *  @return ent log-enthalpy \e H  defined by
	 *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
	 *    where \e e  is the (total) energy density, \e p the pressure,
	 *    \e n  the baryon density, and \f$m_B\f$ the baryon mass.
	 *
	 */
    	virtual double ent_nbar_p(double nbar, const Param* par=0x0) const ;

	/** Computes the total energy density from the baryon density and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *         baryon density 
         *  @param par possible extra parameters of the EOS
	 *
	 *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_nbar_p(double nbar, const Param* par=0x0) const ;

	/** Computes the pressure from the baryon density and extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *         baryon density 
	 *  @param par possible extra parameters of the EOS
	 *
	 *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_nbar_p(double nbar, const Param* par=0x0) const ;

	/** Computes the sound speed squared \f$ c_s^2 = c^2 \frac{dp}{de}\f$
	 *  from the baryon density with extra parameters.
	 *
	 *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *         baryon density 
         *  @param par possible extra parameters of the EOS
  	 *
	 *  @return \f$c_s^2 \f$ [unit: \e c^2]
	 */
    	virtual double csound_square_nbar_p(double nbar, const Param* par=0x0) const  ;
       
};

 		    //------------------------------------//
		    //        class Dyn_eos_tab		  //
		    //------------------------------------//


/**
 * Class for tabulated equations of state for use in dynamical code. \ingroup (eos)
 *
 * Data are to be stored in a formatted file either in the standard LORENE format,
 * or taken from the <a href="http://compose.obspm.fr">CompOSE</a> database.
 * Standard format : the first five lines are comments, preceded by hashes. 
 * Then is given the number of data lines. After, three other lines of comments, 
 * data are given in four columns. 
 * The first one is just a number, not used by LORENE. Second is the baryon number 
 * density in \f${\rm fm}^{-3}\f$, then is the total energy density (including rest mass) 
 * and pressure, both in cgs units. Here is an example from the file \c eos_fps.d:
 * \verbatim #  Date: Tue, 21 Nov 2000 17:24:31 +0100
#  From: xxx
#  FPS e.o.s.:   BPS below n.drip, then FPS. Supplied by N.Stergioulas
#  (June 1998) crust bottom at n=0.0957 , rho=1.60E14 1.E50,1.60E14
#
129  <-- Number of lines
#
#            n_B [fm^{-3}]       rho [g/cm^3]        p [dyn/cm^2]
#  
   6     0.2722041559E-13    0.4518590000E+02    0.1701030000E+15
   7     0.1279071657E-12    0.2123260000E+03    0.5817150000E+16  \endverbatim
 * CompOSE data: use the files XXX.nb and XXX.thermo taken from the database.
 * When built with \c Dyn_eos::eos_from_file(), the file must 
 * be composed of the following lines:
 * \verbatim 17	Type of the EOS 
1	0: standard format	1: CompOSE format 
Tabulated EoS
/full/path/to/the/eos/table/name_of_the_table \endverbatim
 * On the second line '0' means that the table has the standard LORENE format
 * for tabulated EoSs.
 * Note that tables must be ordered from lowest densities up and that
 * units used in the table file \b are \b not those of LORENE.
 * The interpolation through the tables is a cubic Hermite interpolation, which is
 * thermodynamically consistent, i.e. preserves the Gibbs-Duhem relation. 
 * It is defined in [Nozawa, Stergioulas, Gourgoulhon \& Eriguchi,
 * \a Astron. \a Astrophys. Suppl. Ser.  \b 132 , 431 (1998)],
 * and derives from a general technique presented in
 * [Swesty, \a J. \a Comp. \a Phys.  \b 127 , 118 (1996)].
 */
 class Dyn_eos_tab : public Dyn_eos {
   
   // Data :
   // -----
   
 protected:
   /// Name of the file containing the tabulated data
   string tablename ;
   
   string authors ; ///<Authors - reference for the table
   
   bool compose_format ; ///< Are(is) the table(s) in CompOSE format?
   
   /// Lower boundary of the baryon density interval
   double nbmin ;
   
   /// Upper boundary of the baryon density interval
   double nbmax ;
   
   /// Table of \f$\log n_b\f$
   Tbl* lognb ;
   
   /// Table of \f$\log e\f$
   Tbl* loge ;
   
   /// Table of \f$d\log e/d\log n_b\f$
   Tbl* dlesdlnb ;
   
   /// Table of \f$c_s^2 = c^2 \frac{d p}{d e}\f$
   Tbl* c_sound2 ;

   // Constructors - Destructor
   // -------------------------
 public:
   
   /** Standard constructor.
    *
    * @param name_i Name of the EoS
    * @param table_name Name of the file containing the table
    *                   (case compose = false), or the prefix
    *                   i.e. without the .nb or .thermo of the
    *                   two files (case compose = true)
    * @param compose Are the file(s) in CompOSE format?
    */
   Dyn_eos_tab(const string& name_i, const string& table_name,
	       bool compose = true) ;	
   
 protected:
   /// Default constructor to be called by derived classes.
   Dyn_eos_tab() ;
   
   /** Constructor from a binary file (created by the function
    *  \c sauve(FILE*) ).
    *  This constructor is protected because any EOS construction
    *  from a binary file must be done via the function
    * \c Dyn_eos::eos_from_file(FILE*) .
    */
   Dyn_eos_tab(FILE* ) ;
   
   /** Constructor from a formatted file.
    *  This constructor is protected because any EOS construction
    *  from a formatted file must be done via the function
    *  \c  Dyn_eos::eos_from_file(ifstream\& ) .
    *
    *   @param ist input file stream containing a name as first line
    *		and the full filename (including the path) containing 
    *              the EOS file as second line
    */
   Dyn_eos_tab(ifstream& ist) ;
	
 private:
   /** Copy constructor (private to make \c  Dyn_eos_tab
    *  a non-copiable class)
    */	
   Dyn_eos_tab(const Dyn_eos_tab& ) ;
   
   
   /// The construction functions from a file
   friend Dyn_eos* Dyn_eos::eos_from_file(FILE* ) ;
   friend Dyn_eos* Dyn_eos::eos_from_file(ifstream& ) ;

 public:
   virtual ~Dyn_eos_tab() ;			///< Destructor
   
   // Miscellaneous
   // -------------
   
 protected: 	
   /** Reads the files .nb and .thermo containing the table in
    * CompOSE format and initializes
    * the arrays \c  lognb , \c  loge  and \c  dlesdlnb .
    */
   virtual void read_table_compose() ;
   
   /** Reads the file containing the table in LORENE format and initializes
    * the arrays \c  lognb , \c  loge  and \c  dlesdlnb .
    */
   virtual void read_table_lorene() ;

 public :
   /// Comparison operator (egality)
   virtual bool operator==(const Dyn_eos& ) const ;
   
   /// Comparison operator (difference)
   virtual bool operator!=(const Dyn_eos& ) const ;
   
   /** Returns a number to identify the sub-classe of \c Dyn_eos  the
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
   /** Computes the log-enthalpy from the baryon density and extra parameters
    *  (virtual function implemented in the derived classes).
    *
    *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
    *         baryon density
    *  @param par possible extra parameters of the EOS
    *
    *  @return ent log-enthalpy \e H  defined by
    *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
    *    where \e e  is the (total) energy density, \e p the pressure,
    *    \e n  the baryon density, and \f$m_B\f$ the baryon mass.
    *
    */
   virtual double ent_nbar_p(double nbar, const Param* par=0x0) const ;
   
   /** Computes the total energy density from the baryon density and extra parameters
    *  (virtual function implemented in the derived classes).
    *
    *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
    *         baryon density 
    *  @param par possible extra parameters of the EOS
    *
    *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
    *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
    */
   virtual double ener_nbar_p(double nbar, const Param* par=0x0) const ;
   
   /** Computes the pressure from the baryon density and extra parameters
    *  (virtual function implemented in the derived classes).
    *
    *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
    *         baryon density 
    *  @param par possible extra parameters of the EOS
    *
    *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
    *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
    */
   virtual double press_nbar_p(double nbar, const Param* par=0x0) const ;
   
   /** Computes the sound speed squared \f$ c_s^2 = c^2 \frac{dp}{de}\f$
    *  from the baryon density with extra parameters.
    *
    *  @param nbar [input, unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
    *         baryon density 
    *  @param par possible extra parameters of the EOS
    *
    *  @return \f$c_s^2 \f$ [unit: \e c^2]
    */
   virtual double csound_square_nbar_p(double nbar, const Param* par=0x0) const  ;
       
 };

		    //------------------------------------//
		    //	     class Dyn_eos_cons      	  //
		    //------------------------------------//


 /**
  * Equation of state for the <a href="http://compose.obspm.fr">CompOSE</a> 
  * database with a consistent computation of the baryon density.
  * (derived from \c Dyn_eos_tab ). 
  *
  * General tabulated EOS, reading a table passed as an argument to the 
  * constructor. The baryon density \f$n_B\f$ is computed to ensure the relation 
  * \f$ de/dn = h \f$, thus eventually modifying the table. See documentation of
  * the class \c Dyn_eos_tab for details of the table format.
  */
 class Dyn_eos_cons : public Dyn_eos_tab {

   // Constructors - Destructor
   // -------------------------
 public:

   /** Standard constructor.
    *
    * @param name_i Name of the EoS
    * @param table_name Name of the file containing the table
    *                   (case compose = false), or the prefix
    *                   i.e. without the .nb or .thermo of the
    *                   two files (case compose = true)
    * @param compose Are the file(s) in CompOSE format?
    */
   Dyn_eos_cons(const string& name_i, const string& table_name,
	       bool compose = true) ;
  
 protected:
   /** Constructor from a binary file (created by the function
    *  \c sauve(FILE*) ).
    *  This constructor is protected because any EOS construction
    *  from a binary file must be done via the function
    * \c Dyn_eos::eos_from_file(FILE*) .
    */
   Dyn_eos_cons(FILE* ) ;
   
   /** Constructor from a formatted file.
    *  This constructor is protected because any EOS construction
    *  from a formatted file must be done via the function
    *  \c  Dyn_eos::eos_from_file(ifstream\& ) .
    *
    *   @param ist input file stream containing a name as first line
    *		and the full filename (including the path) containing 
    *              the EOS file as second line
    */
   Dyn_eos_cons(ifstream& ist) ;
	
 private:
   /** Copy constructor (private to make \c  Dyn_eos_cons
    *  a non-copiable class)
    */	
   Dyn_eos_cons(const Dyn_eos_cons& ) ;
   
   /// The construction functions from a file
   friend Dyn_eos* Dyn_eos::eos_from_file(FILE* ) ;
   friend Dyn_eos* Dyn_eos::eos_from_file(ifstream& ) ;
  
 public:
  virtual ~Dyn_eos_cons() ;			///< Destructor
  
  // Miscellaneous
  // -------------
  
 protected: 	
   /** Reads the files .nb and .thermo containing the table in
    * CompOSE format and initializes
    * the arrays \c  lognb , \c  loge  and \c  dlesdlnb .
    */
   virtual void read_table_compose() ;
   
   /** Reads the file containing the table in LORENE format and initializes
    * the arrays \c  lognb , \c  loge  and \c  dlesdlnb .
    */
   virtual void read_table_lorene() ;

 public :
  /// Comparison operator (egality)
  virtual bool operator==(const Dyn_eos& ) const ;
  
  /// Comparison operator (difference)
  virtual bool operator!=(const Dyn_eos& ) const ;
  
  /** Returns a number to identify the sub-classe of \c  Dyn_eos  the
   *  object belongs to.
   */
  virtual int identify() const ;
  
  // Outputs
  // -------
  
 protected:
  virtual ostream& operator>>(ostream &) const ;    ///< Operator >>
  
 };

 

}

#endif
