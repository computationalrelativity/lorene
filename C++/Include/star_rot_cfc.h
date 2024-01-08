/*
 *  Definition of Lorene class Star_rot_CFC
 *
 */

/*
 *   Copyright (c) 2024 Jerome Novak
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

#ifndef __STAR_ROT_CFC_H_ 
#define __STAR_ROT_CFC_H_ 

// Headers Lorene
#include "star.h"


namespace Lorene {
  /**
   * Class for relativistic rotating stars in Conformal Flatness Condition 
   * and maximal slicing. \ingroup (star)
   * 
   */
  class Star_rot_CFC : public Star {
    
    // Data : 
    // -----
  protected:
    /** 
     * Relativistic flag.
     * \c 0 : Newtonian theory (not implemented)
     * \c 3 : CFC with \f$\hat{A}^{ij}_{TT} = 0\f$ (see Cordero-Carrion 
     *        et al. 2009).
     */
    int relat_type ;
  
    /**
     * Spectral exponential filtering order.
     * If 0, no filtering is done (see also \c Scalar::exponential_filter_r).
     * Filtering is performed only in shells containing matter, i.e. for
     * domain numbers \e l such that \f$0 < l <\f$ \c nzet .
     */
    int spectral_filter ;
    
    double omega ;  ///< Rotation angular velocity (\c [f_unit] )
    
    // Quantities related to the conformal factor
    //--------------------------------------------
    
    Scalar psi ; ///< Conformal factor \f$\Psi\f$

    
    // Fluid quantities
    //-----------------------------------------
    
    /**
     * Momentum density 3-vector with respect to the 
     * Eulerian observer
     */
    Vector j_euler ; 
    Scalar v2 ; ///< \f$\gamma_{ij}v^i v^j\f$
    
    // Metric stuff 
    //-------------------------------------
    
    Scalar psi4 ;   ///< Conformal factor \f$\Psi^4\f$
    Scalar psi2 ;   ///< \f$\Psi^2\f$
    
    /// flat metric \f$f_{ij}\f$ (spherical components)
    const Metric_flat& flat ; 
    
    Sym_tensor hatA ; ///< \f$\hat{A}^{ij}\f$ 
    Scalar hatA_quad ; ///< \f$f_{il}\, f_{jm}\, \hat{A}_{ij} \hat{A}^{lm}\f$ 

    // Derived data : 
    // ------------
  protected:
    
    // More to come later.....
    //----------------------------
    
    mutable double* p_angu_mom ; ///< Angular momentum. 
    mutable double* p_grv2 ; ///< Error on the virial identity GRV2.
    mutable double* p_grv3 ; ///< Error on the virial identity GRV3.
    mutable double* p_tsw ; ///< Ratio T/W.
    mutable double* p_r_circ ; ///<Circumferential equatorial radius.
    mutable double* p_rp_circ ; ///<Circumferential polar radius.
    

    // Constructors - Destructor
    // -------------------------
  public:
    
    /**
     *
     * Standard constructor.
     *
     * @param mp_i Mapping on which the star will be defined
     * @param nzet_i Number of domains occupied by the star
     * @param eos_i Equation of state of the stellar matter
     * @param relat_i Reltivity parameter (see \c relat_type )
     * @param filter order for spectral exponential filtering
     */
    Star_rot_CFC(Map& mp_i, int nzet_i, const Eos& eos_i, int relat_i=3,
		   int filter=0) ;  

    Star_rot_CFC(const Star_rot_CFC& ) ; ///< Copy constructor

    
    /** Constructor from a file (see \c sauve(FILE*) ).
     *
     * @param mp_i Mapping on which the star will be defined
     * @param eos_i Equation of state of the stellar matter
     * @param fich  input file (must have been created by the function
     *      \c sauve )
     */
    Star_rot_CFC(Map& mp_i, const Eos& eos_i, FILE* fich) ;
    
    
    virtual ~Star_rot_CFC() ;	 ///< Destructor
    
    
    // Memory management
    // -----------------
  protected:
    
    /// Deletes all the derived quantities
    virtual void del_deriv() const ; 
    
    /// Sets to \c 0x0 all the pointers on derived quantities
    void set_der_0x0() const ; 
    
    /** Sets to \c ETATNONDEF  (undefined state) the hydrodynamical
     *  quantities relative to the Eulerian observer.
     */
    virtual void del_hydro_euler() ;
    
    
    // Mutators / assignment
    // ---------------------
  public:
    
    /// Assignment to another \c Star_rot_CFC
    void operator=(const Star_rot_CFC& ) ;	
    
    // Accessors
    // ---------
  public:

    /// Returns the relativity parameter
    int get_relat() const {return relat_type;};

    /// Checks whether the star is computed using a relativistic theory
    bool is_relativistic() const {return (relat_type > 0) ;}; 
    
    /// Returns the filtering order
    int spectral_filter_order() const {return spectral_filter;};
    
    /** 
     * Returns the rotation angular velocity \f$\Omega\f$
     */
    double get_omega() const {return omega;} ;
    
    
    /**
     * Returns the conformal factor \f$\Psi^4\f$
     */
    const Scalar& get_psi4() const {return psi4;} ;
    
    /**
     * Returns \f$\Psi^2\f$
     */
    const Scalar& get_psi2() const {return psi2;} ;
    
    /**
     * Returns \f$\Psi\f$
     */
    const Scalar& get_psi() const {return psi;} ;
    
    // Fluid stuff
    //------------------
    
    /**
     * Returns the momentum density 3-vector with respect to the 
     * Eulerian observer
     */
    const Vector& get_j_euler() const {return j_euler;} ;
    
    /**
     * Returns \f$\gamma_{ij}v^i v^j\f$
     */
    const Scalar& get_v2() const {return v2;} ;
    
    
    //Metric stuff
    //-------------------
    /**
     * Returns \f$\hat{A}^{ij}\f$
     */
    const Sym_tensor get_hatA() const {return hatA;} ;
    
    /** 
     * Returns \f$\tilde{A}_{ij} A^{ij}\f$
     */
    const Scalar get_hatA_quad() const {return hatA_quad;} ;
    
    
    // Outputs
    // -------
  public:
    
    virtual void sauve(FILE* ) const ;	    ///< Save in a file
    
  protected:
    
    virtual ostream& operator>>(ostream& ) const ;
    
   
    // Global quantities
    //-------------------------
  public:
    
    virtual double mass_b() const ; ///< Baryonic mass 
    virtual double mass_g() const ; ///< Gravitational mass
    virtual double angu_mom() const ; ///< Angular momentum
    virtual double grv2() const ;  ///< Error on the virial identity GRV2
    virtual double grv3() const ; ///< Error on the virial identity GRV3
    virtual double tsw() const ; ///< Ratio T/W
    virtual double aplat() const ; ///< Flattening r_pole/r_eq
    virtual double r_circ() const ; ///< Circumferential equatorial radius. 
    virtual double rp_circ() const ; ///< Circumferential polar radius. 
    /**
     * Ellipticity \e e.
     * Defined as \f$ e = \sqrt{1 - \left( \frac{R^c_e}{R^c_p} \right)^2} \f$,
     * where \f$R^c_e\f$ and \f$R^c_p\f$ are, respectively, the equatorial 
     * and polar circumferential radius, given by \c r_circ() and \c rp_circ().
     */
    virtual double ellipt() const ;




    // Computational routines
    //--------------------------
  public:

    /**
     * Computes the hydrodynamical quantities relative to the Eulerian 
     * observer from those in the fluid frame.
     * 
     * More later......
     */ 
    virtual void hydro_euler() ; 

       
    /**
     * Computes metric quantities from known potentials.
     * 
     * The calculation is performed starting from \c psi, \c logn,
     * \c hatA, which are supposed to be up to date.
     * From these, the following fields are updated: \c nnn, 
     * \c psi4, \c psi2,\c shift, and \c hatA_quad. 
     */
    void update_metric() ;


    /**
     * Computes an equilibrium configuration 
     */ 
    void equilibrium(double ent_c, double omega0, double fact_omega, 
		     int nzadapt, const Tbl& ent_limit,
		     const Itbl& icontrol, const Tbl& control,
		     double mbar_wanted, double aexp_mass, 
		     Tbl& diff)  ;

    /**
     * Solution of the ``matter'' part of the Poisson equation for the lapse
     * for rotating stars in CFC.
     */
    void solve_logn_f(Scalar& ln_f_new) const ;

    /**
     * Solution of the quadratic part of the Poisson equation for the lapse
     * for rotating stars in CFC.
     */
    void solve_logn_q(Scalar& ln_q_new) const ;

    /**
     * Solution of the equations for the conformal factor for rotating 
     * stars in CFC
     */
    void solve_psi(Scalar& psi_new) ;

    /**
     * Solution of the shift equation for rotating stars 
     * in CFC
     */
    void solve_shift(Vector& shift_new) const ;

  };

}
#endif
