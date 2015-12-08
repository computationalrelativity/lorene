/*
 *  Definition of Lorene class Hot_eos.
 *
 */

/*
 *   Copyright (c) 2015 Jerome Novak
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

#ifndef __HOTEOS_H_ 
#define __HOTEOS_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.3  2015/12/08 10:52:17  j_novak
 * New class Hoteos_tabul for tabulated temperature-dependent EoSs.
 *
 * Revision 1.2  2015/09/10 13:28:00  j_novak
 * New methods for the class Hot_Eos
 *
 * Revision 1.1  2015/03/17 14:19:59  j_novak
 * New class Hot_eos to deal with temperature-dependent EOSs.
 *
 *
 * $Header$
 *
 */

//C++ headers
#include "headcpp.h"

//C headers
#include<cstdio>
#include "tbl.h"

namespace Lorene{

class Scalar ;
class Param ;
class Eos ;  
		    //------------------------------------//
		    //		class Hot_eos		  //
		    //------------------------------------//

  /**
   * Base class for temperature-dependent equations of state (abstract class).
   * \ingroup(eos)
   * 
   */
  class Hot_eos {
    
    // Data : 
    // -----
  protected:
    string name ;      ///< EOS name
    
    // Constructors - Destructor
    // -------------------------
  protected:
    Hot_eos() ;			///< Standard constructor
    
    /// Standard constructor from a name (string)
    explicit Hot_eos(const string&) ;
    
    /// Standard constructor from a name (char*)
    explicit Hot_eos(const char*) ;
    
    Hot_eos(const Hot_eos& ) ;		///< Copy constructor
    
    /** Constructor from a binary file (created by the function 
     *  \c sauve(FILE*) ). 
     *  This constructor is protected because any hot EOS construction
     *  from a binary file must be done via the function 
     *  \c Hot_eos::hoteos_from_file(FILE*) . 
     */
    Hot_eos(FILE* ) ;    		
    
    /** Constructor from a formatted file.
     *  This constructor is protected because any hot EOS construction
     *  from a formatted file must be done via the function 
     *  \c Hot_eos::hoteos_from_file(ifstream&) . 
     */
    Hot_eos(ifstream& ) ; 
    
  public:
    virtual ~Hot_eos() ;			///< Destructor
 
    // Derived data : 
    // ------------
  protected:
    mutable Eos* p_cold_eos ;     ///< Corresponding cold Eos.
    
    /// Deletes all the derived quantities
    virtual void del_deriv() const ; 

    /// Sets to \c 0x0 all the pointers on derived quantities
    void set_der_0x0() const ; 
	
    // Name manipulation
    // -----------------
  public:
    /// Returns the hot EOS name
    const string& get_name() const {return name; };
    
    /// Sets the hot EOS name
    void set_name(const char* ) ; 
	
    // Miscellaneous
    // -------------
  public:
    /** Construction of an EOS from a binary file.
     *  The file must have been created by the function \c sauve(FILE*) .
     */
    static Hot_eos* hoteos_from_file(FILE* ) ; 
	
    /** Construction of a hot EOS from a formatted file.
     * 
     *  The fist line of the file must start by the EOS number, according 
     *  to the following conventions:
     *  - 1 = relativistic ideal gas (class \c Ideal_gas ). 
     *  - 2 = non-relativistic ideal gas (class \c Ideal_gas_norel ). 
     *
     *  The second line in the file should contain a name given by the user to the EOS.
     *  The following lines should contain the EOS parameters (one
     *  parameter per line), in the same order than in the class declaration.
     */
    static Hot_eos* hoteos_from_file(ifstream& ) ; 
	
    /// Comparison operator (egality)
    virtual bool operator==(const Hot_eos& ) const = 0 ; 

    /// Comparison operator (difference)
    virtual bool operator!=(const Hot_eos& ) const = 0 ; 
    
    /** Returns a number to identify the sub-classe of \c Hot_eos the
     *  object belongs to. 
     */
    virtual int identify() const = 0 ; 
    
    // Outputs
    // -------
    
  public: 
    virtual void sauve(FILE* ) const ;	///< Save in a file

    /// Display
    friend ostream& operator<<(ostream& , const Hot_eos& ) ;	

  protected: 
    virtual ostream& operator>>(ostream &) const = 0 ;    ///< Operator >>

  public:
    /// Returns the corresponding cold \c Eos.
    virtual const Eos& new_cold_Eos() const = 0 ; 


    // Computational functions
    // -----------------------
  protected:
    /**  General computational method for \c Scalar 's
     *
     *  @param thermo1 [input] first thermodynamical quantity (for instance the
     *	    enthalpy field) from which the thermodynamical quantity \c resu  
     *      is to be computed.
     *  @param thermo2 [input] second thermodynamical quantity (for instance the
     *	    entropy field) from which the thermodynamical quantity \c resu  
     *      is to be computed.
     *  @param nzet  [input] number of domains where \c resu  is to be
     *	    computed.
     *  @param l_min [input] index of the innermost domain is which \c resu 
     *	    is to be computed [default value: 0]; \c resu  is computed only in 
     *      domains whose indices are in \c [l_min,l_min+nzet-1] . In the other
     *	    domains, it is set to zero.
     *  @param fait [input] pointer on the member function of class
     *	    \c Hot_eos which performs the pointwise calculation.
     *  @param resu [output] result of the computation.
     */
    void calcule(const Scalar& thermo1, const Scalar& thermo2, int nzet, int l_min,
		 double (Hot_eos::*fait)(double, double) const, Scalar& resu) const ;

  public:
    /** Computes the baryon density from the log-enthalpy and entropy per baryon
     *  (virtual function implemented in the derived classes).
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} (to be modified) \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *
     *  @return baryon density [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
     *
     */
    virtual double nbar_Hs_p(double ent, double sb) const = 0 ;

    /** Computes the baryon density field from the log-enthalpy field and
     * entropy per baryon.
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *  @param nzet  number of domains where the baryon density is to be
     *	  computed.
     *  @param l_min  index of the innermost domain is which the baryon
     *	density is
     *	to be computed [default value: 0]; the baryon density is
     *	computed only in domains whose indices are in
     *      \c [l_min,l_min+nzet-1] . In the other
     *	domains, it is set to zero.
     *
     *  @return baryon density [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
     *
     */
    Scalar nbar_Hs(const Scalar& ent, const Scalar& sb, int nzet, int l_min = 0) const  ;

    /** Computes the total energy density from the log-enthalpy and entropy per baryon
     *  (virtual function implemented in the derived classes).
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *
     *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
     *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
     */
    virtual double ener_Hs_p(double ent, double sb) const = 0 ;
 
    /** Computes the total energy density from the log-enthalpy and entropy per baryon.
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *  @param nzet  number of domains where the energy density is to be
     *	computed.
     *  @param l_min  index of the innermost domain is which the energy
     *	density is
     *	to be computed [default value: 0]; the energy density is
     *	computed only in domains whose indices are in
     *      \c [l_min,l_min+nzet-1] . In the other
     *	domains, it is set to zero.
     *
     *  @return energy density [unit: \f$\rho_{\rm nuc} c^2\f$], where
     *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
     */
    Scalar ener_Hs(const Scalar& ent, const Scalar& sb, int nzet, int l_min = 0) const ;

    /** Computes the pressure from the log-enthalpy and entropy per baryon
     *  (virtual function implemented in the derived classes).
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *
     *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
     *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
     */
    virtual double press_Hs_p(double ent, double sb) const = 0 ;

    /** Computes the pressure from the log-enthalpy and entropy per baryon.
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *  @param nzet  number of domains where the pressure is to be
     *	computed.
     *  @param l_min  index of the innermost domain is which the pressure is
     *	  to be computed [default value: 0]; the pressure is computed
     *    only in domains whose indices are in \c [l_min,l_min+nzet-1] . 
     *    In the other domains, it is set to zero.
     *
     *  @return pressure [unit: \f$\rho_{\rm nuc} c^2\f$], where
     *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
     *
     */
    Scalar press_Hs(const Scalar& ent, const Scalar& sb, int nzet, int l_min = 0) const ;

    /** Computes the temperature from the log-enthalpy and entropy per baryon
     *  (virtual function implemented in the derived classes).
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} (to be modified) \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *
     *  @return temperature [unit: MeV]
     *
     */
    virtual double temp_Hs_p(double ent, double sb) const = 0 ;

    /** Computes the temperature field from the log-enthalpy field and
     * entropy per baryon.
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *  @param nzet  number of domains where the baryon density is to be
     *	  computed.
     *  @param l_min  index of the innermost domain is which the baryon
     *	density is
     *	to be computed [default value: 0]; the baryon density is
     *	computed only in domains whose indices are in
     *      \c [l_min,l_min+nzet-1] . In the other
     *	domains, it is set to zero.
     *
     *  @return temperature [unit: MeV]
     *
     */
    Scalar temp_Hs(const Scalar& ent, const Scalar& sb, int nzet, int l_min = 0) const  ;

  };
  ostream& operator<<(ostream& , const Hot_eos& ) ;	

                    //------------------------------------//
		    //		class Ideal_gas		  //
		    //------------------------------------//

  /**
   * Ideal-gas (temperature-dependent) equation of state, with mass-term 
   * in the energy density. 
   * 
   * \f[
   *    p(n, s_b) = \kappa n^\gamma e^{(\gamma-1)s_b}\ .\qquad (1)
   * \f]
   * and \f[
   *    e(n, s_b) = \frac{\kappa}{\gamma - 1} n^\gamma e^{(\gamma-1)s_b} + m_0\, n\ .
   *    \qquad (2) \f]
   * ### (to be written...)
   *
   *\ingroup (eos) 
   *
   */
  class Ideal_gas : public Hot_eos {
    
    // Data :
    //-------
    
  protected:
    /// Adiabatic index \f$\gamma\f$
    double gam ;
    
    /** Pressure coefficient \f$\kappa\f$  (cf. Eq. (1))
     *  [unit: \f$\rho_{\rm nuc} c^2 / n_{\rm nuc}^\gamma\f$], where
     *  \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$ and
     *  \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$.
     */
    double kap ; 
    
    /** Individual particule mass \f$m_0\f$  (cf. Eq. (2))
     *  [unit: \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$].
     */
    double m_0 ;
    
    double gam1 ;	    ///< \f$\gamma-1\f$
    double unsgam1 ;    ///< \f$1/(\gamma-1)\f$
    double gam1sgamkap ; ///< \f$(\gamma-1) / (\gamma \kappa) m_0\f$
    
    // Constructors - Destructor
    // -------------------------
  public:
    
    /** Standard constructor.
     *
     *  Unless specified, the individual particle mass \f$m_0\f$ is set 
     *  to the mean baryon mass \f$m_B = 1.66\ 10^{-27} \ {\rm kg}\f$.
     *
     *  @param gamma  adiabatic index \f$\gamma\f$
     *  @param kappa  pressure coefficient \f$\kappa\f$
     *  @param mass  individual particule mass \f$m_0\f$
     */
    Ideal_gas(double gamma, double kappa, double mass=1.) ;

    Ideal_gas(const Ideal_gas& ) ;	///< Copy constructor
	
    protected:
	/** Constructor from a binary file (created by the function
	 *  \c sauve(FILE*) ).
	 *  This constructor is protected because any hot EOS construction
	 *  from a binary file must be done via the function 
	 *  \c Hot_eos::eos_from_file(FILE*) .
	 */
    Ideal_gas(FILE* ) ;

    /** Constructor from a formatted file.
     *  This constructor is protected because any EOS construction
     *  from a formatted file must be done via the function
     *  \c Hot_eos::hoteos_from_file(ifstream&) . 
     */
    Ideal_gas(ifstream& ) ;

    /// The construction functions from a file
    friend Hot_eos* Hot_eos::hoteos_from_file(FILE* ) ;
    friend Hot_eos* Hot_eos::hoteos_from_file(ifstream& ) ; 

  public:
    virtual ~Ideal_gas() ;			///< Destructor

    // Assignment
    // ----------
    /// Assignment to another \c Ideal_gas 
    void operator=(const Ideal_gas& ) ;

    // Miscellaneous
    // -------------
    
  public :
    /// Comparison operator (egality)
    virtual bool operator==(const Hot_eos& ) const ;
    
    /// Comparison operator (difference)
    virtual bool operator!=(const Hot_eos& ) const ;
    
    /** Returns a number to identify the sub-classe of \c Hot_eos the
     *  object belongs to.
     */
    virtual int identify() const ; 
    
    /// Returns the adiabatic index \f$\gamma\f$ (cf. Eq. (1)).
    double get_gam() const ;

    /// Returns the pressure coefficient \f$\kappa\f$  (cf. Eq. (1)).
    double get_kap() const ;
	
    /** Return the individual particule mass \f$m_0\f$
     *  (cf. Eq. (1))
     */
    double get_m_0() const ;

    virtual const Eos& new_cold_Eos() const ;

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
    /** Computes the baryon density from the log-enthalpy and entropy per baryon
     *  (virtual function implemented in the derived classes).
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} (to be modified) \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *
     *  @return baryon density [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
     *
     */
    virtual double nbar_Hs_p(double ent, double sb) const ;

    /** Computes the total energy density from the log-enthalpy and entropy per baryon
     *  (virtual function implemented in the derived classes).
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *
     *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
     *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
     */
    virtual double ener_Hs_p(double ent, double sb) const ;
 
    /** Computes the pressure from the log-enthalpy and entropy per baryon
     *  (virtual function implemented in the derived classes).
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *
     *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
     *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
     */
    virtual double press_Hs_p(double ent, double sb) const ;

    /** Computes the temperature from the log-enthalpy and entropy per baryon
     *  (virtual function implemented in the derived classes).
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} (to be modified) \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *
     *  @return temperature [unit: MeV]
     *
     */
    virtual double temp_Hs_p(double ent, double sb) const ;

};

                    //------------------------------------//
		    //	      class Hoteos_tabul	  //
		    //------------------------------------//

  /**
   * Hot (temperature-dependent) tabulated equation of state, read from a file. 
   * 
   *
   *\ingroup (eos) 
   *
   */
  class Hoteos_tabul : public Hot_eos {
    
    // Data :
    //-------
    
  protected:
    /// Name of the file containing the tabulated data
    string tablename ;
    
    string authors ; ///<Authors - reference for the table

    /// Lower boundary of the enthalpy interval
    double hmin ;
    	
    /// Upper boundary of the enthalpy interval
    double hmax ;
    
    /// Lower boundary of the entropy interval
    double sbmin ;
    	
    /// Upper boundary of the entropy interval
    double sbmax ;
    
    /// Table of \f$H = \log ( e + P ) / n_B\f$
    Tbl* hhh ;
    
    /// Table of \f$s_B\f$, entropy per baryon (in units of Boltzmann constant).
    Tbl* s_B ;
    
    /// Table of pressure $P$
    Tbl* ppp ;
    
    /// Table of \f$\partial P/\partial H\f$
    Tbl* dpdh ;
    
    /// Table of \f$\partial P/\partial s_B\f$
    Tbl* dpds ;
    
    /// Table of \f$\partial^2 P/\partial s_B \partial H\f$
    Tbl* d2p ;
    
    // Constructors - Destructor
    // -------------------------
  public:
    
    /** Standard constructor from a filename.
     */
    Hoteos_tabul(const string& filename) ;

    Hoteos_tabul(const Hoteos_tabul& ) ;	///< Copy constructor
	
    protected:
    /** Constructor from a binary file (created by the function
     *  \c sauve(FILE*) ).
     *  This constructor is protected because any hot EOS construction
     *  from a binary file must be done via the function 
     *  \c Hot_eos::eos_from_file(FILE*) .
     */
    Hoteos_tabul(FILE* ) ;

    /** Constructor from a formatted file.
     *  This constructor is protected because any EOS construction
     *  from a formatted file must be done via the function
     *  \c Hot_eos::hoteos_from_file(ifstream&) . 
     */
    Hoteos_tabul(ifstream& ) ;

    /// The construction functions from a file
    friend Hot_eos* Hot_eos::hoteos_from_file(FILE* ) ;
    friend Hot_eos* Hot_eos::hoteos_from_file(ifstream& ) ; 

  public:
    virtual ~Hoteos_tabul() ;			///< Destructor

    /// Assignment to another \c Hoteos_tabul 
    void operator=(const Hoteos_tabul& ) ;

    // Miscellaneous
    // -------------

  protected: 	
    /** Reads the file containing the table and initializes
     *  in the arrays \c hhh , \c s_B, \c ppp, ...
     */
    void read_table() ;

    /// Sets all the arrays to the null pointer.
    void set_arrays_0x0() ;
    
  public :
    /// Comparison operator (egality)
    virtual bool operator==(const Hot_eos& ) const ;
    
    /// Comparison operator (difference)
    virtual bool operator!=(const Hot_eos& ) const ;
    
    /** Returns a number to identify the sub-classe of \c Hot_eos the
     *  object belongs to.
     */
    virtual int identify() const ; 
    
    virtual const Eos& new_cold_Eos() const ;

    // Outputs
    // -------

  public:
    virtual void sauve(FILE* ) const ;	///< Save in a file

  protected:
    virtual ostream& operator>>(ostream &) const ;    ///< Operator >>


    // Computational functions
    // -----------------------

  public:
    /** Computes the baryon density from the log-enthalpy and entropy per baryon
     *  (virtual function implemented in the derived classes).
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} (to be modified) \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *
     *  @return baryon density [unit: \f$n_{\rm nuc} := 0.1 \ {\rm fm}^{-3}\f$]
     *
     */
    virtual double nbar_Hs_p(double ent, double sb) const ;

    /** Computes the total energy density from the log-enthalpy and entropy per baryon
     *  (virtual function implemented in the derived classes).
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *
     *  @return energy density \e e  [unit: \f$\rho_{\rm nuc} c^2\f$], where
     *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
     */
    virtual double ener_Hs_p(double ent, double sb) const ;
 
    /** Computes the pressure from the log-enthalpy and entropy per baryon
     *  (virtual function implemented in the derived classes).
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *
     *  @return pressure \e p [unit: \f$\rho_{\rm nuc} c^2\f$], where
     *      \f$\rho_{\rm nuc} := 1.66\ 10^{17} \ {\rm kg/m}^3\f$
     */
    virtual double press_Hs_p(double ent, double sb) const ;

    /** Computes the temperature from the log-enthalpy and entropy per baryon
     *  (virtual function implemented in the derived classes).
     *
     *  @param ent [input,  unit: \f$c^2\f$] log-enthalpy \e H  defined by
     *    \f$H = c^2 \ln\left( {e+p \over m_B c^2 n} (to be modified) \right) \f$,
     *    where \e e  is the (total) energy density, \e p the pressure,
     *    \e n  the baryon density, and \f$m_B\f$ the baryon mass
     *  @param sb [input,  unit: \f$k_B\f$] entropy per baryon \f$s_b\f$
     *
     *  @return temperature [unit: MeV]
     *
     */
    virtual double temp_Hs_p(double ent, double sb) const ;

};

}
#endif
