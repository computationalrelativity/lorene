/*
 *  Definition of Lorene class Eos_fitting
 *                             Eos_fit_SLy4
 *                             Eos_fit_FPS
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

#ifndef __EOS_FITTING_H_ 
#define __EOS_FITTING_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/09/26 18:50:00  k_taniguchi
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

// External classes which appear in the declaration of class Eos_fitting:
class Tbl ; 
class Cmp ;
class Param ;
class Eos ;

//-----------------------------------------------------------------------//
//     class Eos_fitting for the analytical fitting of realistic EOS     //
//-----------------------------------------------------------------------//

/**
 * Base class for the analytically fitted equation of state.
 * \ingroup(eos)
 * 
 */
class Eos_fitting : public Eos {

    // Data : 
    // -----
    protected:
        /// Name of the file containing the fitting data
        char dataname[160] ;

	/// Array of the coefficients of the fitting data
	double* pp ;

    // Constructors - Destructor
    // -------------------------
    protected:

        /** Standard constructor
	 *
	 * @param name Name of the fitted EOS
	 * @param data Name of the file containing the fitting data
	 * @param path Path to the directory containing the EOS file
	 */
	Eos_fitting(const char* name_i, const char* data, const char* path) ;

	//	Eos_fitting(const Eos_fitting& ) ;	///< Copy constructor

    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Eos::eos_from_file(FILE*) . 
	 */
	Eos_fitting(FILE* ) ;

	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Eos::eos_from_file(ifstream&) .
	 *
	 * @param ist input file stream containing a name as first line
	 * @param data Name of the file containing the fitting data
	 */
	Eos_fitting(ifstream& ist, const char* data) ;

	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_fitting() ;		///< Destructor
 
    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file

    // Miscellaneous
    // -------------
    protected:
	/** Reading coefficients of the fitting equation for the energy
	 *  density, the pressure, and the enthalpy
	 */
	void read_coef() ;

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


                    //--------------------------------------//
                    //          class Eos_fit_SLy4          //
                    //--------------------------------------//

/**
 * Fitted equation of state of SLy4. \ingroup(eos)
 *
 */
class Eos_fit_SLy4 : public Eos_fitting {

    // Constructors - Destructor
    // -------------------------
    public:

        /** Standard constructor
	 *
	 * @param path Path to the directory containing the EOS file
	 */
	Eos_fit_SLy4(const char* path) ;

    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Eos::eos_from_file(FILE*) . 
	 */
	Eos_fit_SLy4(FILE* ) ;

	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Eos::eos_from_file(ifstream&) .
	 *
	 * @param ist input file stream containing a name as first line
	 * @param data Name of the file containing the fitting data
	 */
	Eos_fit_SLy4(ifstream& ) ;

	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_fit_SLy4() ;		///< Destructor
 
    // Outputs
    // -------
    protected:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

    // Miscellaneous
    // -------------
    public :
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ;

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ;

	/** Returns a number to identify the sub-classe of \c  Eos  the
	 *  object belongs to.
	 */
	virtual int identify() const ;

};


                    //-------------------------------------//
                    //          class Eos_fit_FPS          //
                    //-------------------------------------//

/**
 * Fitted equation of state of FPS. \ingroup(eos)
 *
 */
class Eos_fit_FPS : public Eos_fitting {

    // Constructors - Destructor
    // -------------------------
    public:

        /** Standard constructor
	 *
	 * @param path Path to the directory containing the EOS file
	 */
	Eos_fit_FPS(const char* path) ;

    protected:
	/** Constructor from a binary file (created by the function 
	 *  \c sauve(FILE*) ). 
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Eos::eos_from_file(FILE*) . 
	 */
	Eos_fit_FPS(FILE* ) ;

	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function 
	 *  \c Eos::eos_from_file(ifstream&) .
	 *
	 * @param ist input file stream containing a name as first line
	 * @param data Name of the file containing the fitting data
	 */
	Eos_fit_FPS(ifstream& ) ;

	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_fit_FPS() ;		///< Destructor
 
    // Outputs
    // -------
    protected:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

    // Miscellaneous
    // -------------
    public :
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ;

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ;

	/** Returns a number to identify the sub-classe of \c  Eos  the
	 *  object belongs to.
	 */
	virtual int identify() const ;

};


#endif
