/*
 *  Definition of Lorene class Eos_compose_fit
 */

/*
 *   Copyright (c) 2022 Jerome Novak
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


#ifndef __EOS_COMPOSE_FIT_H_
#define __EOS_COMPOSE_FIT_H_

/*
 * $Id$
 * $Log$
 * Revision 1.4  2023/01/27 16:10:35  j_novak
 * A polytrope (class Eos_poly) is used for low and high enthalpies.
 *
 * Revision 1.3  2022/07/21 12:33:51  j_novak
 * Improved comments
 *
 * Revision 1.2  2022/07/20 12:59:43  j_novak
 * Added methods for saving into files and constructor from a file
 *
 * Revision 1.1  2022/04/15 13:39:24  j_novak
 * New class Eos_compose_fit to generate fitted EoSs from CompOSE tables.
 *
 *
 *
 *
 * $Header$
 *
 */

// Standard C++
#include "eos.h"
#include "scalar.h"

// Lorene classes
namespace Lorene {


		    //------------------------------------//
		    //	     class Eos_compose_fit    	  //
		    //------------------------------------//


/**
 * Equation of state for fitting the 1-parameter EoSs from 
 * the <a href="http://compose.obspm.fr">CompOSE</a> database. \ingroup(eos)
 *
 * A polynomial fit is done ont the adidabatic index for the highest densities. 
 * At the lowest ones, a simple polytrope is fitted. In between, a linear 
 * interpolation is done in terms of adiabatic index. When built with 
 * \c Eos::eos_from_file(), the file must be composed of the following lines:
 * \verbatim 7	Type of the EOS
EoS fitted from CompOSE data
/full/path/to/the/eos/table/
0.01   1.e-9  # limiting density values
6             # polynomial degree for the regression
33            # number of Chebyshev point sin each domain \endverbatim
 * On the second line the name is just for output. The path to 
 * the directory containing CompOSE EoS (third line), does not contain 
 * the name of the tables, it is assumed to be called \c eos.nb and 
 * \c eos.thermo (see CompOSE documentation).
 * The fourth line gives limit values in density for the fit: above the first 
 * value, the polynomial regression is done. Below the second value, 
 * a polytrope is fitted. 
 */
class Eos_compose_fit : public Eos {

  // Data
  //--------
  
 protected:
  /// Name of the file containing the tabulated data
  string tablename ;
    	
  /** Number of coeficients for polynomial regression. 
   *  Note: this is in general different from the number of spectral coefficients
   *  used for the representation of thermodynamical fields.
   */
  int n_coefs ;

  /// Lower bound in baryon density, below which the EoS is assumed to be a polytrope.
  double nb_min ;

  /** Middle bound in baryon density, above which which the EoS is determined 
   *  from the polynomial fit to the adiabatic index.
   */
  double nb_mid ;

  /// Higher bound on the density, above which the EoS is assumed to be a polytrope.
  double nb_max ;

  /// Values of enthalpy corresponding to nb_min and nb_max
  double hmin, hmax ;

  /// Pointer on a polytropic EoS for the description of low densities (nb<nb_min)
  const Eos_poly* p_eos_low ;

  /// Pointer on a polytropic EoS for the description of high densities (nb>nb_max)
  const Eos_poly* p_eos_high ;

  /// Multi-grid defining the number of domains and spectral coefficients
  const Mg3d* mg ;
  
  /// Mapping in \f$ x = \log H\f$
  const Map_af* mp ;
    	
  /// Table of \f$\log p\f$
  Scalar* log_p ;
    	
  /// Table of \f$\log e\f$
  Scalar* log_e ;

  /// Table of \f$\log n_b\f$
  Scalar* log_nb ;
    	
  /// Table of \f$\log c_s^2 = \log \left( c^2 \frac{d p}{d e} \right) \f$
  Scalar* log_cs2 ;

  
  // Constructors - Destructor
  // -------------------------
 public:
  
  /** Constructor from a parameter file.
   *
   * @param files_path Absolute name (including path) of the parameter file.
   */
  Eos_compose_fit(const string& param_file) ;	

	
 protected:
  /** Constructor from a binary file (created by the function
   *  \c sauve(FILE*) ).
   *  This constructor is protected because any EOS construction
   *  from a binary file must be done via the function
   * \c Eos::eos_from_file(FILE*) .
   */
  Eos_compose_fit(FILE* ) ;
  
  /** Constructor from a formatted file.
   *  This constructor is protected because any EOS construction
   *  from a formatted file must be done via the function
   *  \c  Eos::eos_from_file(ifstream\& ) .
   */
  Eos_compose_fit(ifstream&) ;
  
 private:	
  /** Copy constructor (private to make \c  Eos_compose_fit 
   *  a non-copiable class)
   */	
  Eos_compose_fit(const Eos_compose_fit& ) ;	

  /// The construction functions from a file
  friend Eos* Eos::eos_from_file(FILE* ) ;
  friend Eos* Eos::eos_from_file(ifstream& ) ;
  
 public:
  virtual ~Eos_compose_fit() ;			///< Destructor

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

protected:

  /// Reads the Compose data and makes the fit.
  void read_and_compute(ifstream&) ;

  /// Reads Compose data and stores the values of thermodynamic quantities. 
  void read_compose_data(int& nbp, Tbl*& logh, Tbl*& logp, Tbl*& loge,
			 Tbl*& lognb, Tbl*& gam1) ;

  /** From the read values, makes the fit on the adiabatic index and deduces 
   *  the other quantities from thermodynamic relations and definitions.
   */
  void adiabatic_index_fit(int i_min, int i_max, const Tbl& logh_read,
			   const Tbl& logp_read, const Tbl& loge_read,
			   const Tbl& lognb_read, const Tbl& gam1_read) ;
  
  // Outputs
  // -------
 public:

  virtual void sauve(FILE* ) const ;	///< Save in a file.

  /// Save into a table in Lorene format.
  void write_lorene_table(const string&, int nlines = 200) const ; 
  
 protected:
  
  virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

  // Accessors
  //-------------

  /// Returns the name (given in the parameter file, see the introduction of the class).
  const string& get_tablename() const {return tablename ;} ;

  double get_nbmin() const { return nb_min ;} ;  ///< Lower bound in density
  double get_nbmax() const { return nb_max ;} ;  ///< Higher bound in density
  
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
		
	/** Computes the sound speed squared \f$ c_s^2 = c^2 \frac{dp}{de}\f$
	 *  from the enthapy with extra parameters
	 *  (virtual function implemented in the derived classes).
	 *
	 *  @param ent [input, unit: \e c^2]
	 *         enthalpy 
	 *  @param par possible extra parameters of the EOS
	 *
	 *  @return \f$c_s^2 \f$ [unit: \e c^2]
	 */
	virtual double csound_square_ent_p(double, const Param* par=0x0) const ;
  
};



}
#endif

