/*
 *  Definition of Lorene classes Eos_tabul
 *			         Eos_SLy4
 *			         Eos_FPS
 *				 Eos_BPAL12
 *				 Eos_AkmalPR
 *				 Eos_BBB2
 *				 Eos_BalbN1H1
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


#ifndef __EOS_TABUL_H_
#define __EOS_TABUL_H_

/*
 * $Id$
 * $Log$
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
class ostream ;
class ifstream ;

// Headers C
#include <stdio.h>

// Lorene classes
class Tbl ;
class Cmp ;

		
		    //------------------------------------//
		    //		class Eos_tabul		  //
		    //------------------------------------//


/**
 * Base class for tabulated equations of state.
 *
 * The interpolation through the tables is
 * a cubic Hermite interpolation, which is
 * thermodynamically consistent, i.e. preserves the
 * Gibbs-Duhem relation. It is defined in
 * [Nozawa, Stergioulas, Gourgoulhon \& Eriguchi,
 * {\sl Astron. Astrophys. Suppl. Ser.} {\bf 132}, 431 (1998)],
 * and derives from a general technique presented in
 * [Swesty, {\bf J. Comp. Phys.} {\bf 127}, 118 (1996)].
 *
 *
 * @version #$Id$#                                                              * @version #$Id$#
 */
class Eos_tabul : public Eos {

    // Data :
    // -----

    protected:
    	/// Name of the file containing the tabulated data
    	char tablename[160] ;
    	
    	/// Lower boundary of the enthalpy interval
    	double hmin ;
    	
    	/// Upper boundary of the enthalpy interval
    	double hmax ;
    	
    	/// Table of $\log H$
    	Tbl* logh ;
    	
    	/// Table of $\log p$
    	Tbl* logp ;
    	
    	/// Table of $d\log P/d\log H$
    	Tbl* dlpsdlh ;

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

	Eos_tabul(const Eos_tabul& ) ;	/// Copy constructor	
	
    protected:
	
	/** Constructor from a binary file (created by the function
	 *  {\tt sauve(FILE* )}).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 *  {\tt Eos::eos\_from\_file(FILE* )}.
	 */
	Eos_tabul(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}.
	 *
	 *   @param ist input file stream containing a name as first line
	 *		and the path to the directory containing the EOS file
	 *		as second line
	 *   @param table Name of the file containing the EOS table
	 */
	Eos_tabul(ifstream& ist, const char* table) ;
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_tabul() ;			/// Destructor


    // Miscellaneous
    // -------------

    protected: 	
    	/** Reads the file containing the table and initializes
    	 *  in the arrays {\tt logh}, {\tt logp} and {\tt dlpsdlh}.
    	 */
    	void read_table() ;


    // Outputs
    // -------

    public:
	virtual void sauve(FILE* ) const ;	/// Save in a file


    // Computational functions
    // -----------------------

    public:
	/** Computes the baryon density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy {\it H}
	 *
	 *  @return baryon density {\it n} [unit: $n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}$]
	 *
	 */
    	virtual double nbar_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the total energy density from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy {\it H}
	 *
	 *  @return energy density {\it e} [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double ener_ent_p(double ent, const Param* par=0x0) const ;

 	/** Computes the pressure from the log-enthalpy.
	 *
	 *  @param ent [input,  unit: $c^2$] log-enthalpy {\it H}
	 *
	 *  @return pressure {\it p} [unit: $\rho_{\rm nuc} c^2$], where
	 *      $\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3$
	 */
    	virtual double press_ent_p(double ent, const Param* par=0x0) const ;

	/** Computes the logarithmic derivative $d\ln n/d\ln H$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy {\it H} 
	 *
	 *  @return dln(n)/dln(H)
	 */
    	virtual double der_nbar_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln e/d\ln H$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy {\it H} 
	 *
	 *  @return dln(e)/dln(H)
	 */
    	virtual double der_ener_ent_p(double ent, const Param* par=0x0) const ; 
       
	/** Computes the logarithmic derivative $d\ln p/d\ln H$ 
	 * from the log-enthalpy. 
	 * 
	 *  @param ent [input,  unit: $c^2$] log-enthalpy {\it H} 
	 *
	 *  @return dln(p)/dln(H)
	 */
    	virtual double der_press_ent_p(double ent, const Param* par=0x0) const ; 

};

		    //------------------------------------//
		    //		class Eos_SLy4 	  	  //
		    //------------------------------------//


/**
 * Equation of state SLy4 (Douchin \& Haensel 2001).
 *
 * Interior: neutrons, protons, electrons and muons described by the Skyrme
 * Lyon 4 potential.
 * 
 * Inner crust: Douchin \& Haensel 2001
 *
 * Outer crust: Haensel \& Pichon,  Astron. Astrophys. {\bf 283},  
 *  313 (1994)
 *
 * @version #$Id$#                                                              * @version #$Id$#
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
	 *  {\tt sauve(FILE* )}).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 *  {\tt Eos::eos\_from\_file(FILE* )}.
	 */
	Eos_SLy4(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}.
	 */
	Eos_SLy4(ifstream& ) ;

    private:	
	/** Copy constructor (private to make {\tt Eos\_SLy4}
	 *  a non-copiable class)
	 */	
	Eos_SLy4(const Eos_SLy4& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_SLy4() ;			/// Destructor

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

    protected:
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


};

		    //------------------------------------//
		    //		class Eos_FPS 	  	  //
		    //------------------------------------//


/**
 * Equation of state FPS (Friedman-Pandharipande + Skyrme).
 * Authors: Lorenz, Ravenhall \& Pethick (unpublished).
 *
 *
 * @version #$Id$#                                                              * @version #$Id$#
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
	 *  {\tt sauve(FILE* )}).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 *  {\tt Eos::eos\_from\_file(FILE* )}.
	 */
	Eos_FPS(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}.
	 */
	Eos_FPS(ifstream& ) ;

    private:	
	/** Copy constructor (private to make {\tt Eos\_FPS}
	 *  a non-copiable class)
	 */	
	Eos_FPS(const Eos_FPS& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_FPS() ;			/// Destructor

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

    protected:
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


};

		    //------------------------------------//
		    //		class Eos_BPAL12 	  //
		    //------------------------------------//


/**
 * Equation of state BPAL12 (Bombaci et al 1995).
 *
 *
 * @version #$Id$#                                                              * @version #$Id$#
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
	 *  {\tt sauve(FILE* )}).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 *  {\tt Eos::eos\_from\_file(FILE* )}.
	 */
	Eos_BPAL12(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}.
	 */
	Eos_BPAL12(ifstream& ) ;

    private:	
	/** Copy constructor (private to make {\tt Eos\_BPAL12}
	 *  a non-copiable class)
	 */	
	Eos_BPAL12(const Eos_BPAL12& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_BPAL12() ;			/// Destructor

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

    protected:
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


};


		    //------------------------------------//
		    //		class Eos_AkmalPR 	  //
		    //------------------------------------//


/**
 * Equation of state AkmalPR (Akmal, Pandharipande \& Ravenhall 1998).
 *
 * Interior: neutrons, protons, electrons and muons described by the 
 *  A18+dv+UIX* model of Akmal, Pandharipande \& Ravenhall, 
 *  Phys. Rev. C {\bf 58},  1804 (1998)
 * 
 * Inner crust: SLy4
 *
 * Outer crust: BPS + Haensel \& Pichon,  Astron. Astrophys. {\bf 283},  
 *  313 (1994)
 *
 * @version #$Id$#                                                              * @version #$Id$#
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
	 *  {\tt sauve(FILE* )}).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 *  {\tt Eos::eos\_from\_file(FILE* )}.
	 */
	Eos_AkmalPR(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}.
	 */
	Eos_AkmalPR(ifstream& ) ;

    private:	
	/** Copy constructor (private to make {\tt Eos\_AkmalPR}
	 *  a non-copiable class)
	 */	
	Eos_AkmalPR(const Eos_AkmalPR& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_AkmalPR() ;			/// Destructor

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

    protected:
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


};

		    //------------------------------------//
		    //		class Eos_BBB2 	          //
		    //------------------------------------//


/**
 * Equation of state BBB2 (Baldo, Bombaci \& Burgio 1997).
 *
 * Interior: BHF (Paris +TBF) model of Baldo, Bombaci \& Burgio, 
 *  Astron. Astrophys. {\bf 328}, 274 (1997)
 * 
 * Crust: SLy
 *
 *
 * @version #$Id$#                                                              * @version #$Id$#
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
	 *  {\tt sauve(FILE* )}).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 *  {\tt Eos::eos\_from\_file(FILE* )}.
	 */
	Eos_BBB2(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}.
	 */
	Eos_BBB2(ifstream& ) ;

    private:	
	/** Copy constructor (private to make {\tt Eos\_BBB2}
	 *  a non-copiable class)
	 */	
	Eos_BBB2(const Eos_BBB2& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_BBB2() ;			/// Destructor

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

    protected:
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


};


		    //------------------------------------//
		    //		class Eos_BalbN1H1 	  //
		    //------------------------------------//


/**
 * Equation of state BalbN1H1 (Balberg 2000).
 *
 *
 * @version #$Id$#                                                              * @version #$Id$#
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
	 *  {\tt sauve(FILE* )}).
	 *  This constructor is protected because any EOS construction
	 *  from a binary file must be done via the function
	 *  {\tt Eos::eos\_from\_file(FILE* )}.
	 */
	Eos_BalbN1H1(FILE* ) ;
	
	/** Constructor from a formatted file.
	 *  This constructor is protected because any EOS construction
	 *  from a formatted file must be done via the function
	 *  {\tt Eos::eos\_from\_file(ifstream\& )}.
	 */
	Eos_BalbN1H1(ifstream& ) ;

    private:	
	/** Copy constructor (private to make {\tt Eos\_BalbN1H1}
	 *  a non-copiable class)
	 */	
	Eos_BalbN1H1(const Eos_BalbN1H1& ) ;	
	
	
	/// The construction functions from a file
	friend Eos* Eos::eos_from_file(FILE* ) ;
	friend Eos* Eos::eos_from_file(ifstream& ) ;

    public:
	virtual ~Eos_BalbN1H1() ;			/// Destructor

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

    protected:
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>


};



#endif

