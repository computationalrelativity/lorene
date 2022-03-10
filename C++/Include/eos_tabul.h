/*
 *  Definition of Lorene classes Eos_tabul
 *			         Eos_SLy4
 *			         Eos_FPS
 *				 Eos_BPAL12
 *				 Eos_AkmalPR
 *				 Eos_BBB2
 *				 Eos_BalbN1H1
 *                               Eos_GlendNH3
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *             (c) 2014 Jerome Novak
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


#ifndef __EOS_TABUL_H_
#define __EOS_TABUL_H_

/*
 * $Id$
 * $Log$
 * Revision 1.21  2022/03/10 16:38:39  j_novak
 * log(cs^2) is tabulated instead of cs^2.
 *
 * Revision 1.20  2021/05/14 15:39:22  g_servignat
 * Added sound speed computation from enthalpy to Eos class and tabulated+polytropic derived classes
 *
 * Revision 1.19  2021/05/06 14:33:17  j_novak
 * New conversion function from Eos to Dyn_eos.
 *
 * Revision 1.18  2019/04/09 12:50:22  j_novak
 * Improved documentation
 *
 * Revision 1.17  2019/03/28 13:41:01  j_novak
 * Improved managed of saved EoS (functions sauve and constructor form FILE*)
 *
 * Revision 1.16  2015/08/04 14:41:28  j_novak
 * Back to previous version for Eos_CompOSE. Enthalpy-consistent EoS can be accessed using Eos_consistent class (derived from Eos_CompOSE).
 *
 * Revision 1.15  2015/01/27 14:22:38  j_novak
 * New methods in Eos_tabul to correct for EoS themro consistency (optional).
 *
 * Revision 1.14  2014/10/13 08:52:34  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.13  2014/07/01 09:26:20  j_novak
 * Improvement of comments
 *
 * Revision 1.12  2014/06/30 16:13:18  j_novak
 * New methods for reading directly from CompOSE files.
 *
 * Revision 1.11  2014/03/06 15:53:34  j_novak
 * Eos_compstar is now Eos_compOSE. Eos_tabul uses strings and contains informations about authors.
 *
 * Revision 1.10  2010/02/02 14:26:10  j_novak
 * *** empty log message ***
 *
 * Revision 1.9  2010/02/02 13:21:52  j_novak
 * New class Eos_Compstar.
 *
 * Revision 1.8  2004/03/22 13:12:41  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.7  2003/12/08 15:48:43  m_bejger
 * GlendNH3 EOS (Glendenning 1985, case 3) added
 *
 * Revision 1.6  2003/11/25 13:44:15  m_bejger
 * Declared some vectors for Eos_tabul::read_table()
 *
 * Revision 1.5  2003/11/21 16:19:09  m_bejger
 * Added new tables: lognb, dlpsdlnb
 *
 * Revision 1.4  2002/10/16 14:36:29  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.3  2002/09/13 09:17:31  j_novak
 * Modif. commentaires
 *
 * Revision 1.2  2002/04/09 14:32:15  e_gourgoulhon
 * 1/ Added extra parameters in EOS computational functions (argument par)
 * 2/ New class MEos for multi-domain EOS
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.5  2001/09/11  16:15:46  eric
 * Ajout des classes Eos_BBB2 et Eos_BalbN1H1
 *
 * Revision 2.4  2001/09/11  15:05:48  eric
 * Ajout de la classe Eos_AkmalPR
 *
 * Revision 2.3  2001/03/23  13:40:23  eric
 * Modifs commentaires.
 *
 * Revision 2.2  2001/02/07  09:45:28  eric
 * Suppression de la fonction derent_ent_p.
 * Ajout des fonctions donnant les derivees de l'EOS:
 *  	der_nbar_ent_p
 * 	der_ener_ent_p
 * 	der_press_ent_p
 *
 * Revision 2.1  2000/11/23  22:33:48  eric
 * Ajout de Eos_BPAL12.
 *
 * Revision 2.0  2000/11/22  19:29:18  eric
 * *** empty log message ***
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
		    //		class Eos_tabul		  //
		    //------------------------------------//


/**
 * Base class for tabulated equations of state. \ingroup (eos)
 *
 * EoS data are to be stored in a formatted file in the following format.
 * The first five lines are comments, preceded by hashes. Then is given the number of 
 * data lines. After, three other lines of comments, data are given in four columns. 
 * The first one is just a number, not used by Lorene. Second is the baryon number 
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
 * Note that units used in the table file \b are \b not those of Lorene.
 * The interpolation through the tables is
 * a cubic Hermite interpolation, which is
 * thermodynamically consistent, i.e. preserves the
 * Gibbs-Duhem relation. It is defined in
 * [Nozawa, Stergioulas, Gourgoulhon \& Eriguchi,
 * \a Astron. \a Astrophys. Suppl. Ser.  \b 132 , 431 (1998)],
 * and derives from a general technique presented in
 * [Swesty, \a J. \a Comp. \a Phys.  \b 127 , 118 (1996)].
 *
 *
 */
class Eos_tabul : public Eos {

    // Data :
    // -----

    protected:
    	/// Name of the file containing the tabulated data
    	string tablename ;
    	
	string authors ; ///<Authors - reference for the table

    	/// Lower boundary of the enthalpy interval
    	double hmin ;
    	
    	/// Upper boundary of the enthalpy interval
    	double hmax ;
    	
    	/// Table of \f$\log H\f$
    	Tbl* logh ;
    	
    	/// Table of \f$\log p\f$
    	Tbl* logp ;
    	
    	/// Table of \f$d\log P/d\log H\f$
    	Tbl* dlpsdlh ;

    	/// Table of \f$\log n_b\f$
    	Tbl* lognb ;
    	
        /// Table of \f$d\log P/d\log nb\f$
        Tbl* dlpsdlnb ;
		
	/// Table of \f$\log c_s^2 = \log \left( c^2 \frac{d p}{d e} \right) \f$
	Tbl* log_cs2 ;

        double* press ; 
        double* nb ; 
        double* ro ;
        
                
    // Constructors - Destructor
    // -------------------------
    protected:

	/** Standard constructor.
	 *
	 * @param name_i Name of the equation of state
	 * @param table Name of the file containing the EOS table
	 * @param path Path to the directory containing the EOS file
	 */
	Eos_tabul(const char* name_i, const char* table, const char* path) ;	

	/** Standard constructor from the full filename.
	 *
	 * @param name_i Name of the equation of state
	 * @param table Full name of the file containing the EOS table
	 *              (including the absolute path).
	 */
	Eos_tabul(const char* name_i, const char* file_name) ;

	Eos_tabul(const Eos_tabul& ) ;	///< Copy constructor	

	Eos_tabul(const char* name_i) ; ///< Standard constructor with a name
	
    protected:
	
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 * \c Eos::eos_from_file(FILE*) .
	 */
	Eos_tabul(FILE* ) ;
	
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
	Eos_tabul(ifstream& ist, const char* table) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c  Eos::eos_from_file(ifstream\& ) .
	 *
	 *   @param ist input file stream containing a name as first line
	 *		and the full filename (including the path) containing 
	 *              the EOS file as second line
	 *   @param table Name of the file containing the EOS table
	 */
	Eos_tabul(ifstream& ist) ;
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_tabul() ;			///< Destructor


    // Miscellaneous
    // -------------
    public:
  const string& get_tablename() const { return tablename ; } ;
  
  const Tbl& get_logh() const {return *logh ; } ;
  
    protected: 	
    	/** Reads the file containing the table and initializes
    	 *  in the arrays \c  logh , \c  logp  and \c  dlpsdlh .
    	 */
    	virtual void read_table() ;


    // Outputs
    // -------

    public:
	virtual void sauve(FILE* ) const ;	///< Save in a file


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
	/* mb test version for linear interpolation 

    	virtual double nbar_ent_p_mbtest(double ent, const Param* par=0x0) const ;
         */       

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
	/* mb test version for linear interpolation 
    	virtual double press_ent_p_mbtest(double ent, const Param* par=0x0) const ;
        */

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
	virtual double csound_square_ent_p(double, const Param*) const ;
};

		    //------------------------------------//
		    //		class Eos_SLy4 	  	  //
		    //------------------------------------//


/**
 * Equation of state SLy4 (Douchin \& Haensel 2001).\ingroup (eos)
 *
 * Interior: neutrons, protons, electrons and muons described by the Skyrme
 * Lyon 4 potential.
 * 
 * Inner crust: Douchin \& Haensel 2001
 *
 * Outer crust: Haensel \& Pichon,  \a Astron. \a Astrophys. \b 283 ,  
 *  313 (1994)
 *
 */
class Eos_SLy4 : public Eos_tabul {


    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor.
	 *
	 * @param path Path to the directory containing the EOS file
	 */
	Eos_SLy4(const char* path) ;	

	
    protected:
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 * \c Eos::eos_from_file(FILE*) .
	 */
	Eos_SLy4(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c  Eos::eos_from_file(ifstream\& ) .
	 */
	Eos_SLy4(ifstream& ) ;

    private:	
	/** Copy constructor (private to make \c  Eos_SLy4 
	 *  a non-copiable class)
	 */	
	Eos_SLy4(const Eos_SLy4& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_SLy4() ;			///< Destructor

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

    // Outputs
    // -------

    protected:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

};

		    //------------------------------------//
		    //		class Eos_FPS 	  	  //
		    //------------------------------------//


/**
 * Equation of state FPS (Friedman-Pandharipande + Skyrme).\ingroup (eos)
 * Authors: Lorenz, Ravenhall \& Pethick (unpublished).
 *
 *
 */
class Eos_FPS : public Eos_tabul {


    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor.
	 *
	 * @param path Path to the directory containing the EOS file
	 */
	Eos_FPS(const char* path) ;	

	
    protected:
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 * \c Eos::eos_from_file(FILE*) .
	 */
	Eos_FPS(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c  Eos::eos_from_file(ifstream\& ) .
	 */
	Eos_FPS(ifstream& ) ;

    private:	
	/** Copy constructor (private to make \c  Eos_FPS 
	 *  a non-copiable class)
	 */	
	Eos_FPS(const Eos_FPS& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_FPS() ;			///< Destructor

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

    // Outputs
    // -------

    protected:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

};

		    //------------------------------------//
		    //		class Eos_BPAL12 	  //
		    //------------------------------------//


/**
 * Equation of state BPAL12 (Bombaci et al 1995).\ingroup (eos)
 *
 *
 */
class Eos_BPAL12 : public Eos_tabul {


    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor.
	 *
	 * @param path Path to the directory containing the EOS file
	 */
	Eos_BPAL12(const char* path) ;	

	
    protected:
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 * \c Eos::eos_from_file(FILE*) .
	 */
	Eos_BPAL12(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c  Eos::eos_from_file(ifstream\& ) .
	 */
	Eos_BPAL12(ifstream& ) ;

    private:	
	/** Copy constructor (private to make \c  Eos_BPAL12 
	 *  a non-copiable class)
	 */	
	Eos_BPAL12(const Eos_BPAL12& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_BPAL12() ;			///< Destructor

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

    // Outputs
    // -------

    protected:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

};


		    //------------------------------------//
		    //		class Eos_AkmalPR 	  //
		    //------------------------------------//


/**
 * Equation of state AkmalPR (Akmal, Pandharipande \& Ravenhall 1998).
 *
 * Interior: neutrons, protons, electrons and muons described by the 
 *  A18+dv+UIX* model of Akmal, Pandharipande \& Ravenhall, 
 *  Phys. Rev. C \b 58 ,  1804 (1998) \ingroup (eos)
 * 
 * Inner crust: SLy4
 *
 * Outer crust: BPS + Haensel \& Pichon,  \a Astron. \a Astrophys. \b 283 ,  
 *  313 (1994)
 *
 */
class Eos_AkmalPR : public Eos_tabul {


    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor.
	 *
	 * @param path Path to the directory containing the EOS file
	 */
	Eos_AkmalPR(const char* path) ;	

	
    protected:
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 * \c Eos::eos_from_file(FILE*) .
	 */
	Eos_AkmalPR(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c  Eos::eos_from_file(ifstream\& ) .
	 */
	Eos_AkmalPR(ifstream& ) ;

    private:	
	/** Copy constructor (private to make \c  Eos_AkmalPR 
	 *  a non-copiable class)
	 */	
	Eos_AkmalPR(const Eos_AkmalPR& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_AkmalPR() ;			///< Destructor

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

    // Outputs
    // -------

    protected:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

};

		    //------------------------------------//
		    //		class Eos_BBB2 	          //
		    //------------------------------------//


/**
 * Equation of state BBB2 (Baldo, Bombaci \& Burgio 1997).
 *
 * Interior: BHF (Paris +TBF) model of Baldo, Bombaci \& Burgio, 
 *  \a Astron. \a Astrophys. \b 328 , 274 (1997)
 * 
 * Crust: SLy \ingroup (eos)
 *
 *
 */
class Eos_BBB2 : public Eos_tabul {


    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor.
	 *
	 * @param path Path to the directory containing the EOS file
	 */
	Eos_BBB2(const char* path) ;	

	
    protected:
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 * \c Eos::eos_from_file(FILE*) .
	 */
	Eos_BBB2(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c  Eos::eos_from_file(ifstream\& ) .
	 */
	Eos_BBB2(ifstream& ) ;

    private:	
	/** Copy constructor (private to make \c  Eos_BBB2 
	 *  a non-copiable class)
	 */	
	Eos_BBB2(const Eos_BBB2& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_BBB2() ;			///< Destructor

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

    // Outputs
    // -------

    protected:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

};


		    //------------------------------------//
		    //		class Eos_BalbN1H1 	  //
		    //------------------------------------//


/**
 * Equation of state BalbN1H1 (Balberg 2000). \ingroup (eos)
 *
 *
 */
class Eos_BalbN1H1 : public Eos_tabul {


    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor.
	 *
	 * @param path Path to the directory containing the EOS file
	 */
	Eos_BalbN1H1(const char* path) ;	

	
    protected:
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 * \c Eos::eos_from_file(FILE*) .
	 */
	Eos_BalbN1H1(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c  Eos::eos_from_file(ifstream\& ) .
	 */
	Eos_BalbN1H1(ifstream& ) ;

    private:	
	/** Copy constructor (private to make \c  Eos_BalbN1H1 
	 *  a non-copiable class)
	 */	
	Eos_BalbN1H1(const Eos_BalbN1H1& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_BalbN1H1() ;			///< Destructor

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

    // Outputs
    // -------

    protected:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

};



		    //------------------------------------//
		    //		class Eos_GlendNH3	  //
		    //------------------------------------//


/**
 * Equation of state GlendNH3 (Glendenning 1985, case 3 ).
 *
 *\ingroup (eos)
 * 
 */
class Eos_GlendNH3 : public Eos_tabul {


    // Constructors - Destructor
    // -------------------------
    public:

	/** Standard constructor.
	 *
	 * @param path Path to the directory containing the EOS file
	 */
	Eos_GlendNH3(const char* path) ;	

	
    protected:
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 * \c Eos::eos_from_file(FILE*) .
	 */
	Eos_GlendNH3(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  \c  Eos::eos_from_file(ifstream\& ) .
	 */
	Eos_GlendNH3(ifstream& ) ;

    private:	
	/** Copy constructor (private to make \c  Eos_GlendNH3 
	 *  a non-copiable class)
	 */	
	Eos_GlendNH3(const Eos_GlendNH3& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_GlendNH3() ;			///< Destructor

    // Miscellaneous
    // -------------

    public :
	/// Comparison operator (egality)
	virtual bool operator==(const Eos& ) const ;

	/// Comparison operator (difference)
	virtual bool operator!=(const Eos& ) const ;

	/** Returns a number to identify the sub-classe of \c Eos  the
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

