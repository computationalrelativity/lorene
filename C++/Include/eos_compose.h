/*
 *  Definition of Lorene classes Eos_CompOSE
 *                               Eos_consistent
 */

/*
 *   Copyright (c) 2014-2015 Jerome Novak
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


#ifndef __EOS_COMPOSE_H_
#define __EOS_COMPOSE_H_

/*
 * $Id$
 * $Log$
 * Revision 1.3  2019/04/09 12:50:22  j_novak
 * Improved documentation
 *
 * Revision 1.2  2019/03/28 13:41:01  j_novak
 * Improved managed of saved EoS (functions sauve and constructor form FILE*)
 *
 * Revision 1.1  2015/08/04 14:41:28  j_novak
 * Back to previous version for Eos_CompOSE. Enthalpy-consistent EoS can be accessed using Eos_consistent class (derived from Eos_CompOSE).
 *
 *
 *
 * $Header$
 *
 */

// Standard C++
#include <string>
#include "headcpp.h"

// Headers C
#include <cstdio>

// Lorene classes
namespace Lorene {
class Tbl ;


		    //------------------------------------//
		    //	     class Eos_CompOSE        	  //
		    //------------------------------------//


/**
 * Equation of state for the <a href="http://compose.obspm.fr">CompOSE</a> 
 * database. 
 *
 * General tabulated EOS, reading a table passed as an argument to the 
 * constructor. When built with \c Eos::eos_from_file(), the file must 
 * be composed of the following lines:
 * \verbatim 17	Type of the EOS 
1	0: standard format	1: CompOSE format 
Tabulated EoS
/full/path/to/the/eos/table/name_of_the_table \endverbatim
 * On the second line '0' means that the table has the standard LORENE format
 * for tabulated EoSs (see class \c Eos_tabul for details). 
 * '1' means that the files from the CompOSE database
 * are used and that the 'name_of_the_table' should be without suffix:
 * e.g. \c great_eos would stand for files \c great_eos.nb and 
 * \c great_eos.thermo (see CompOSE documentation).
 */
class Eos_CompOSE : public Eos_tabul {

  // Data
  //--------
  
 protected:
  int format ; ///< 0 for standard (old) LORENE format, 1 for CompOSE format
  
  // Constructors - Destructor
  // -------------------------
 public:
  
  /** Constructor from CompOSE data.
   *
   * @param files_path Absolute name (including path), but without
   *                   extensions of the CompOSE data, e.g.
   *                   \c /home/foo/eos/my_eos standing for
   *                   files \c my_eos.nb and \c my_eos.thermo
   *                   as dowloaded from the <a href="http://compose.obspm.fr">
   *                   CompOSE server</a>.
   */
  Eos_CompOSE(const string& files_path) ;
  
  /** Standard constructor.
   *
   * @param file_name Absolute name (including path) containing 
   *                  the EOS file. This file should contain a header 
   *                  with the following first lines:
   * \verbatim# Comments ...
# Name of the authors / reference
# Comments ...
# Comments ...
# Comments
XXX  <-- Number of lines
#
#  index      n_B [fm^{-3}]  rho [g/cm^3]   p [dyn/cm^2]
# \endverbatim
   * where 'XXX' is the number of following lines of the EoS, each line containing 
   * an index (integer), baryon density, energy total density and pressure 
   * in the units given above.
   */
  Eos_CompOSE(const char* file_name) ;	

	
 protected:
  /** Constructor from a binary file (created by the function
   *  \c sauve(FILE*) ).
   *  This constructor is protected because any EOS construction
   *  from a binary file must be done via the function
   * \c Eos::eos_from_file(FILE*) .
   */
  Eos_CompOSE(FILE* ) ;
  
  /** Constructor from a formatted file.
   *  This constructor is protected because any EOS construction
   *  from a formatted file must be done via the function
   *  \c  Eos::eos_from_file(ifstream\& ) .
   */
  Eos_CompOSE(ifstream&) ;
  
 private:	
  /** Copy constructor (private to make \c  Eos_CompOSE 
   *  a non-copiable class)
   */	
  Eos_CompOSE(const Eos_CompOSE& ) ;	
  
	
  /// The construction functions from a file
  friend Eos* Eos::eos_from_file(FILE* ) ;
  friend Eos* Eos::eos_from_file(ifstream& ) ;
  
 public:
  virtual ~Eos_CompOSE() ;			///< Destructor
  
  // Miscellaneous
  // -------------
  
 protected: 	
  /** Reads the files containing the table and initializes
   *  in the arrays \c  logh , \c  logp  and \c  dlpsdlh (CompOSE format).
   */
  virtual void read_compose_data() ;

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
  
  
};

		    //------------------------------------//
		    //	     class Eos_consistent      	  //
		    //------------------------------------//


/**
 * Equation of state for the <a href="http://compose.obspm.fr">CompOSE</a> 
 * database with a consistent computation of the log-enthalpy 
 * (derived from \c Eos_CompOSE ). 
 *
 * General tabulated EOS, reading a table passed as an argument to the 
 * constructor. The log-enthalpy \f$h\f$ is computed to ensure the relation 
 * \f$ dp = (e+p) dh \f$, thus eventually modifying the table. When built 
 * with \c Eos::eos_from_file(), the file must be composed of the following 
 * lines: \verbatim 20	Type of the EOS 
1	0: standard format	1: CompOSE format 
Tabulated EoS
/full/path/to/the/eos/table/name_of_the_table \endverbatim
 * On the second line '0' means that the table has the standard LORENE format
 * for tabulated EoSs. '1' means that the files from the CompOSE database
 * are used and that the 'name_of_the_table' should be without suffix:
 * e.g. \c great_eos would stand for files \c great_eos.nb and 
 * \c great_eos.thermo (see CompOSE documentation).
 */
class Eos_consistent : public Eos_CompOSE {

    // Constructors - Destructor
    // -------------------------
    public:

  /** Constructor from CompOSE data.
   *
   * @param files_path Absolute name (including path), but without
   *                   extensions of the CompOSE data, e.g.
   *                   \c /home/foo/eos/my_eos standing for
   *                   files \c my_eos.nb and \c my_eos.thermo
   *                   as dowloaded from the <a href="http://compose.obspm.fr">
   *                   CompOSE server</a>.
   */
  Eos_consistent(const string& files_path) ;
  
  /** Standard constructor.
   *
   * @param file_name Absolute name (including path) containing 
   *                  the EOS file. This file should contain a header 
   *                  with the following first lines:
   * \verbatim# Comments ...
# Name of the authors / reference
# Comments ...
# Comments ...
# Comments
XXX  <-- Number of lines
#
#  index      n_B [fm^{-3}]  rho [g/cm^3]   p [dyn/cm^2]
# \endverbatim
   * where 'XXX' is the number of following lines of the EoS, each line containing 
   * an index (integer), baryon density, energy total density and pressure 
   * in the units given above.
   */
  Eos_consistent(const char* file_name) ;	

	
 protected:
  /** Constructor from a binary file (created by the function
   *  \c sauve(FILE*) ).
   *  This constructor is protected because any EOS construction
   *  from a binary file must be done via the function
   * \c Eos::eos_from_file(FILE*) .
   */
  Eos_consistent(FILE* ) ;
  
  /** Constructor from a formatted file.
   *  This constructor is protected because any EOS construction
   *  from a formatted file must be done via the function
   *  \c  Eos::eos_from_file(ifstream\& ) .
   */
  Eos_consistent(ifstream&) ;
  
 private:	
  /** Copy constructor (private to make \c  Eos_consistent 
   *  a non-copiable class)
   */	
  Eos_consistent(const Eos_consistent& ) ;	
  
	
  /// The construction functions from a file
  friend Eos* Eos::eos_from_file(FILE* ) ;
  friend Eos* Eos::eos_from_file(ifstream& ) ;
  
 public:
  virtual ~Eos_consistent() ;			///< Destructor
  
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
  /** Reads the file containing the table and initializes
   *  in the arrays \c  logh , \c  logp  and \c  dlpsdlh .
   */
  virtual void read_table() ;
  
  /** Reads the files containing the table and initializes
   *  in the arrays \c  logh , \c  logp  and \c  dlpsdlh (CompOSE format).
   */
  virtual void read_compose_data() ;

  // Outputs
  // -------
  
 protected:
  virtual ostream& operator>>(ostream &) const ;    ///< Operator >>
  
  // Computational functions
  // -----------------------
  
 public:
  /** Computes the baryon density from the log-enthalpy.
   *
   *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H 
   *
   *  @return baryon density \e n [unit:\f$n_{\rm nuc}:=0.1\ {\rm fm}^{-3}\f$]
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
  
};



}
#endif

