/*
 *  Definition of Lorene class Eos_multi_poly
 *
 */

/*
 *   Copyright (c) 2009 Keisuke Taniguchi
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
 * Revision 1.5  2014/10/06 15:09:39  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2009/06/23 14:33:31  k_taniguchi
 * Completely revised.
 *
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
#include <cstdio>
#include <cassert>

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
 * Base class for a multiple polytropic equation of state.
 *
 * This equation of state mimics some realistic, tabulated EOSs.
 * \ingroup (eos)
 * 
 */
class Eos_multi_poly : public Eos {

    // Data : 
    // -----

    protected:
	/// Number of polytropic equations of state
	int npeos ;

	/// Array (size: \c npeos) of adiabatic index \f$\gamma\f$
	double* gamma ;

	/** Pressure coefficient for the crust
	 *   [unit: \f$({\rm g/cm^3})^{1-\gamma_0}\f$]
	 */
	double kappa0 ;

	/// Exponent of the pressure at the fiducial density \f$\rho_1\f$
	double logP1 ;

	/// Array (size: \c npeos - 1) of the exponent of fiducial densities
	double* logRho ;

	/** Array (size: \c npeos) of pressure coefficient \f$\kappa\f$
	 *   [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$],
	 *   where
	 *   \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *   \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$.
	 */
	double* kappa ;

	/** Array (size \c npeos - 1) of the number density
	 *   at which the polytropic EOS changes its index and constant
	 */
	double* nbCrit ;

	/** Array (size \c npeos - 1) of the critical enthalpy
	 *   at which the polytropic EOS changes its index and constant
	 */
	double* entCrit ;

	/** Array (size \c npeos - 1) of the percentage which detemines
	 *   the terminating enthalpy for lower density and the starting
	 *   enthalpy for higher density
	 */
	double* decInc ;

	/** Individual particule mass \f$m0\f$
	 *   [unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
	 */
	double m0 ;

	/** Array (size: \c npeos) of the relativistic chemical potential
	 *   at zero pressure
	 *   [unit: \f$m_B c^2\f$,
	 *          with \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
         * (The value for the EOS which covers the lowest density: 1)
        */
        double* mu0 ;


    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor (sets \c m0 to 1).
	 *
	 *  The individual particle mass \f$m0\f$ is set to the mean baryon
	 *  mass \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$.
	 *
	 *  @param npoly  number of polytropes
	 *  @param gamma_i  array of adiabatic index \f$\gamma\f$
	 *  @param kappa0_i  pressure coefficient for the crust
	 *  @param logP1_i exponent of the pressure at the fiducial density
	 *  @param logRho_i array of the exponent of fiducial densities
	 *  @param decInc_i array of percentage
	 */
	Eos_multi_poly(int npoly, double* gamma_i, double kappa0_i,
		       double logP1_i, double* logRho_i, double* decInc_i) ;

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

	/// Read/write kappa
	//	double& set_kappa(int n) ;


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

	/// Returns the number of polytropes \c npeos
	const int& get_npeos() const { return npeos ; } ;

	/// Returns the adiabatic index \f$\gamma\f$
	const double& get_gamma(int n) const { 
	    assert(n>=0 && n<npeos) ;
	    return gamma[n] ;
	} ;

	/// Returns the pressure coefficient for the crust
	const double& get_kappa0() const { return kappa0 ; } ;

	/// Returns the exponent of the pressure at the fiducial density
	const double& get_logP1() const { return logP1 ; } ;

	/// Returns the exponent of fiducial densities
	const double& get_logRho(int n) const {
	    assert(n>=0 && n<npeos-1) ;
	    return logRho[n] ;
	} ;

	/** Returns the pressure coefficient \f$\kappa\f$
	 *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$],
	 *  where
	 *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
	 *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$.
	 */
	const double& get_kappa(int n) const { 
	    assert(n>=0 && n<npeos) ;
	    return kappa[n] ;
	} ;

	/// Returns the critical number density
	const double& get_nbCrit(int n) const { 
	    assert(n>=0 && n<npeos-1) ;
	    return nbCrit[n] ;
	} ;

	/// Returns the critical enthalpy
	const double& get_entCrit(int n) const { 
	    assert(n>=0 && n<npeos-1) ;
	    return entCrit[n] ;
	} ;

    protected:
	/// Computes the auxiliary quantities
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

	/** Computes the logarithmic derivative \f$d\ln p/d\ln n\f$
	 *   from the log-enthalpy.
	 *
	 *  @param ent [input, unit: \f$c^2\f$] log-enthalpy \e H
	 *
	 *  @param par possible extra parameters of the EOS
	 *  @return dln(p)/dln(n)
	 */
    	virtual double der_press_nbar_p(double ent, const Param* par=0x0) const ;
};

#endif
