/*
 *  Definition of Lorene class Eos_multi_poly
 *
 */

/*
 *   Copyright (c) 2004 Keisuke Taniguchi
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

#ifndef __EOS_MULTI_POLY_H_
#define __EOS_MULTI_POLY_H_

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/05/14 11:35:17  k_taniguchi
 * Minor changes in some comments.
 *
 * Revision 1.2  2004/05/07 13:04:01  j_novak
 * Forgotten #include<assert.h>
 *
 * Revision 1.1  2004/05/07 08:09:56  k_taniguchi
 * Initial revision
 *
 *
 *
 * $Header$
 *
 */

// Standard C++
#include "headcpp.h"

// Headers C
#include <stdio.h>
#include <assert.h>

// Lorene classes
#include "eos.h"
#include "param.h"
class Tbl ;
class Cmp ;
class Param ;
class Eos ;

		    //-------------------------------------------//
		    //   base class Eos for multiple polytrope   //
		    //-------------------------------------------//

/**
 * Base class for multiple polytropic equation of state.
 *
 * This equation of state mimics some realistic, tabulated EOS.
 * Then the polytropic indices and constants of the pressure are
 * different from those of the energy density.
 * \ingroup (eos)
 * 
 */
class Eos_multi_poly : public Eos {

    // Data : 
    // -----

    protected:
	/// Number of domains for the pressure
	int ndom_press ;

	/// Number of domains for the energy density
	int ndom_epsil ;

	/** Array (size \c ndom_press - 1) of the number density
	 *   at which the polytropic EOS
	 *   for the pressure changes its index and constant
	 */
	double* nb_press ;

	/** Array (size \c ndom_epsil - 1) of the number density
	 *   at which the polytropic EOS
	 *   for the energy density changes its index and constant
	 */
	double* nb_epsil ;

	/** Array (size \c ndom_press - 1) of the enthalpy
	 *   at the points where the polytropic EOS for the pressure
	 *   changes its index and constant
	 */
	double* ent_crit_p ;

	/** Array (size \c ndom_epsil - 1) of the enthalpy
	 *   at the points where the polytropic EOS for the energy density
	 *   changes its index and constant
	 */
	double* ent_crit_e ;

	/** Array (size: \c ndom_press) of adiabatic index \f$\gamma\f$
	 *   for the pressure
	 */
	double* gam_p ;

	/** Array (size: \c ndom_press) of adiabatic index \f$\gamma\f$
	 *   for the energy density
	 */
	double* gam_e ;

	/** Array (size: \c ndom_press) of pressure coefficient \f$\kappa\f$
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$]
	 *  for the pressure, where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$.
	 */
	double* kap_p ;

	/** Array (size: \c ndom_press) of pressure coefficient \f$\kappa\f$
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$]
	 *  for the energy density, where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$.
	 */
	double* kap_e ;

	/** Individual particule mass \f$m_0\f$
	 *  [unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
	 */
	double m_0 ;

	/** Relativistic chemical potential at zero pressure
	 *  [unit: \f$m_B c^2\f$, with \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
         * (standard value: 1)
        */
        double mu_0 ;


    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor (sets both \c m_0 and \c mu_0 to 1).
	 *
	 *  The individual particle mass \f$m_0\f$ is set to the mean baryon
	 *  mass \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$.
	 *
	 *  @param ndom_p  number of polytropes for the pressure
	 *  @param gamma_p  array of adiabatic index \f$\gamma\f$
	 *         for the pressure
	 *  @param kappa_p  array of pressure coefficient \f$\kappa\f$  
	 *	   [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$]
	 *         for the pressure,
	 *	   where \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 *	   and \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$
	 *  @param ndom_e  number of polytropes for the energy density
	 *  @param gamma_e  array of adiabatic index \f$\gamma\f$
	 *         for the energy density
	 *  @param kappa_e  array of pressure coefficient \f$\kappa\f$  
	 *	   [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$]
	 *         for the energy density,
	 *	   where \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 *	   and \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$
	 */
	Eos_multi_poly(int ndom_p, double* gamma_p, double* kappa_p,
		       int ndom_e, double* gamma_e, double* kappa_e) ;

	Eos_multi_poly(const Eos_multi_poly& ) ;     ///< Copy constructor

    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Eos::eos_from_file(FILE*) . 
	 */
	Eos_multi_poly(FILE* ) ;    		

	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Eos::eos_from_file(ifstream&) . 
	 */
	Eos_multi_poly(ifstream& ) ;

	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_multi_poly() ;			///< Destructor
 

    // Assignment
    // ----------
    public:
	/// Assignment to another \c Eos_multi_poly
	void operator=(const Eos_multi_poly&) ;	


    // Miscellaneous
    // -------------
    public:
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ;

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ;

	/** Returns a number to identify the sub-classe of \c Eos
	 *  the object belongs to.
	 */
	virtual int identify() const ;

	/// Returns the number of polytropes \c ndom_press
	int get_ndom_press() const { return ndom_press ; } ;

	/// Returns the number of polytropes \c ndom_epsil
	int get_ndom_epsil() const { return ndom_epsil ; } ;

	/// Returns the adiabatic index \f$\gamma\f$ for the pressure
	double get_gam_press(int n) const { 
	    assert(n>=0 && n<ndom_press) ;
	    return gam_p[n] ;
	} ;

	/// Returns the adiabatic index \f$\gamma\f$ for the energy density
	double get_gam_epsil(int n) const { 
	    assert(n>=0 && n<ndom_epsil) ;
	    return gam_e[n] ;
	} ;

	/** Returns the pressure coefficient \f$\kappa\f$
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$]
	 *  for the pressure, where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$.
	 */
	double get_kap_press(int n) const { 
	    assert(n>=0 && n<ndom_press) ;
	    return kap_p[n] ;
	} ;

	/** Returns the pressure coefficient \f$\kappa\f$
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$]
	 *  for the energy density, where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$.
	 */
	double get_kap_epsil(int n) const { 
	    assert(n>=0 && n<ndom_epsil) ;
	    return kap_e[n] ;
	} ;

    protected:
	/** Computes the auxiliary quantities \c gam1 , \c unsgam1 ,
	 *  \c gam1sgamkap from the values of \c gam and \c kap 
	 */
	void set_auxiliary() ;


    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
    protected: 
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>


    // Computational functions
    // -----------------------
    public:
	/** Computes the baryon density from the log-enthalpy.
	 *
	 *  @param ent [input, unit: \f$c^2\f$] log-enthalpy \e H
	 *
	 *  @param par possible extra parameters of the EOS
	 *  @return baryon density \e n
	 *      [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
	 *
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the total energy density from the log-enthalpy.
	 *
	 *  @param ent [input, unit: \f$c^2\f$] log-enthalpy \e H
	 *
	 *  @param par possible extra parameters of the EOS
	 *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$],
	 *      where \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the pressure from the log-enthalpy.
	 *
	 *  @param ent [input, unit: \f$c^2\f$] log-enthalpy \e H
	 *
	 *  @param par possible extra parameters of the EOS
	 *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
	 *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln n/d\ln H\f$
	 *   from the log-enthalpy.
	 *
	 *  @param ent [input, unit: \f$c^2\f$] log-enthalpy \e H
	 *
	 *  @param par possible extra parameters of the EOS
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln e/d\ln H\f$
	 *   from the log-enthalpy.
	 * 
	 *  @param ent [input, unit: \f$c^2\f$] log-enthalpy \e H
	 *
	 *  @param par possible extra parameters of the EOS
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative \f$d\ln p/d\ln H\f$
	 *   from the log-enthalpy.
	 *
	 *  @param ent [input, unit: \f$c^2\f$] log-enthalpy \e H
	 *
	 *  @param par possible extra parameters of the EOS
	 *  @return dln(p)/dln(H)
	 */
    	virtual double der_press_ent_p(double ent, const Param* par=0x0) const ;

};

#endif
