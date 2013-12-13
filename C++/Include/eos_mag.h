/*
 *  Definition of Lorene classes Eos_mag
 */

/*
 *   Copyright (c) 2011 Jerome Novak
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


#ifndef __EOS_MAG_H_
#define __EOS_MAG_H_

/*
 * $Id$
 * $Log$
 * Revision 1.3  2013/12/13 16:36:51  j_novak
 * Addition and computation of magnetisation terms in the Einstein equations.
 *
 * Revision 1.2  2011/10/04 16:05:18  j_novak
 * Update of Eos_mag class. Suppression of loge, re-definition of the derivatives
 * and use of interpol_herm_2d.
 *
 * Revision 1.1  2011/06/16 10:49:18  j_novak
 * New class Eos_mag for EOSs depending on density and magnetic field.
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
class Tbl ;
class Cmp ;

		
		    //------------------------------------//
		    //		class Eos_mag		  //
		    //------------------------------------//


/**
 * Class for a magnetized (tabulated) equation of state. \ingroup (eos)
 *
 * This EOS depends on two variables \f$(n, B)\f$: density (or log-enthalpy) 
 * and the magnetic field amplitude. The interpolation through the tables is
 * a cubic Hermite interpolation, which is thermodynamically consistent, i.e. 
 * preserves the Gibbs-Duhem relation. 
 *
 */
class Eos_mag : public Eos {

    // Data :
    // -----

    protected:
    	/// Name of the file containing the tabulated data
    	string tablename ;
    	
    	/// Lower boundary of the log-enthalpy interval
    	double hmin ;
    	
    	/// Upper boundary of the log-enthalpy interval
    	double hmax ;
    	
    	/// Upper boundary of the magnetic field interval
    	double Bmax ;
    	
    	/// Table of \f$\log p\f$
    	Tbl* logp ;
    	
    	/// Table of \f$\log h\f$
    	Tbl* logh ;
    	
    	/// Table of \f$B\f$
    	Tbl* Bfield ;
    	
    	/// Table of \f$\frac{d \log p}{d \log h}\f$
    	Tbl* dlpsdlh ;
    	
    	/// Table of \f$\frac{d \log p}{d \log B}\f$
    	Tbl* dlpsdB ;
    	
    	/// Table of \f$\frac{d^2 \log p}{d \log h d \log B}\f$
    	Tbl* d2lp ;
    	        
                
    // Constructors - Destructor
    // -------------------------
    protected:

	/** Standard constructor.
	 *
	 * @param name_i Name of the equation of state
	 * @param table Name of the file containing the EOS table
	 * @param path Path to the directory containing the EOS file
	 */
	Eos_mag(const char* name_i, const char* table, const char* path) ;	

	/** Standard constructor from the full filename.
	 *
	 * @param name_i Name of the equation of state
	 * @param table Full name of the file containing the EOS table
	 *              (including the absolute path).
	 */
	Eos_mag(const char* name_i, const char* file_name) ;

	Eos_mag(const Eos_mag& ) ;	///< Copy constructor	
	
    protected:
	
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 * \c Eos::eos_from_file(FILE*) .
	 */
	Eos_mag(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c  Eos::eos_from_file(ifstream\& ) .
	 *
	 *   @param ist input file stream containing a name as first line
	 *		and the path to the directory containing the EOS file
	 *		as second line
	 *   @param table Name of the file containing the EOS table
	 */
	Eos_mag(ifstream& ist, const char* table) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c  Eos::eos_from_file(ifstream\& ) .
	 *
	 *   @param ist input file stream containing a name as first line
	 *		and the full filename (including the path) containing 
	 *              the EOS file as second line
	 */
	Eos_mag(ifstream& ist) ;
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_mag() ;			///< Destructor


    // Miscellaneous
    // -------------

    protected: 	
    	/** Reads the file containing the table and initializes
    	 *  in the arrays \c  logh , \c  logp  and \c  dlpsdlh .
    	 */
    	void read_table() ;

     public :
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ;

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ;

	/** Returns a number to identify the sub-classe of \c  Eos  the
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
	/** Computes the magnetisation.
	 */
	double mag_ent_p(double ent, const Param* par=0x0) const ;

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
	 *  @return pressure \e p  [unit: \f$\rho_{\rm nuc} c^2\f$], where
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

        /** Computes the logarithmic derivative \f$d\ln p/d\ln n\f$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  
	 *
	 *  @return dln(p)/dln(n)
	 */
    	virtual double der_press_nbar_p(double ent, const Param* par=0x0) const ; 

};


#endif

