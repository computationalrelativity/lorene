/*
 *  Definition of Lorene class Time_slice
 *
 */

/*
 *   Copyright (c) 2004 Eric Gourgoulhon, Jose Luis Jaramillo & Jerome Novak
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

#ifndef __TIME_SLICE_H_ 
#define __TIME_SLICE_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.30  2010/10/20 07:58:09  j_novak
 * Better implementation of the explicit time-integration. Not fully-tested yet.
 *
 * Revision 1.29  2008/12/04 18:22:49  j_novak
 * Enhancement of the dzpuis treatment + various bug fixes.
 *
 * Revision 1.28  2008/12/02 15:02:21  j_novak
 * Implementation of the new constrained formalism, following Cordero et al. 2009
 * paper. The evolution eqs. are solved as a first-order system. Not tested yet!
 *
 * Revision 1.27  2008/08/19 06:41:59  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.26  2007/11/06 14:47:06  j_novak
 * New constructor from a rotating star in Dirac gauge (class Star_rot_Dirac).
 * Evolution can take into account matter terms.
 *
 * Revision 1.25  2007/04/25 15:20:59  j_novak
 * Corrected an error in the initialization of tildeB in
 * Tslice_dirac_max::initial_dat_cts. + New method for solve_hij_AB.
 *
 * Revision 1.24  2007/03/21 14:51:48  j_novak
 * Introduction of potentials A and tilde(B) of h^{ij} into Tslice_dirac_max.
 *
 * Revision 1.23  2005/03/28 19:44:00  f_limousin
 * Function tgam() is now virtual.
 *
 * Revision 1.22  2004/06/24 20:36:07  e_gourgoulhon
 * Class Time_slice_conf: added method check_psi_dot.
 *
 * Revision 1.21  2004/06/14 20:46:35  e_gourgoulhon
 * Added argument method_poisson to Tslice_dirac_max::solve_hij.
 *
 * Revision 1.20  2004/05/31 20:28:20  e_gourgoulhon
 * -- Class Time_slice : added inline functions get_latest_j() and
 *    get_time()
 * -- Class Tslice_dirac_max: method hh_det_one takes now a time step
 *    as argument, to compute h^{ij} from khi and mu at an arbitrary
 *    time step and not only the latest one.
 *
 * Revision 1.19  2004/05/27 15:22:28  e_gourgoulhon
 * Added functions save and sauve, as well as constructors from file.
 *
 * Revision 1.18  2004/05/20 20:30:37  e_gourgoulhon
 * Added arguments check_mod and save_mod to method Tsclice_dirac_max::evolve.
 *
 * Revision 1.17  2004/05/17 19:52:16  e_gourgoulhon
 * -- Method initial_data_cts: added arguments graph_device and
 *    method_poisson_vect.
 * -- Method Tslice_dirac_max::solve_beta : added argument method
 * -- Method Tslice_dirac_max::solve_hij : added argument graph_device
 * -- Method Tslice_dirac_max::evolve : added arguments
 *    method_poisson_vect, nopause and graph_device.
 *
 * Revision 1.16  2004/05/12 15:16:25  e_gourgoulhon
 * Added #include "metric.h" before #include "evolution.h".
 *
 * Revision 1.15  2004/05/09 20:56:09  e_gourgoulhon
 * Added member adm_mass_evol and corresponding virtual method adm_mass().
 *
 * Revision 1.14  2004/05/06 15:23:10  e_gourgoulhon
 * initial_data_cts is know a virtual function of class Time_slice_conf
 * and is implemented also for class Tslice_dirac_max.
 *
 * Revision 1.13  2004/05/03 14:46:11  e_gourgoulhon
 * Class Tslice_dirac_max: -- changed prototype of method solve_hij
 *                         -- added new method evolve
 *
 * Revision 1.12  2004/04/30 14:36:15  j_novak
 * Added the method Tslice_dirac_max::solve_hij(...)
 * NOT READY YET!!!
 *
 * Revision 1.11  2004/04/30 10:51:38  e_gourgoulhon
 * Class Tslice_dirac_max: added methods solve_n, solve_q and solve_beta
 * for resolution of the elliptic part of Einstein equations.
 *
 * Revision 1.10  2004/04/29 17:07:27  e_gourgoulhon
 * Added argument pdt to Time_slice_conf::initial_data_cts.
 *
 * Revision 1.9  2004/04/08 16:42:11  e_gourgoulhon
 * Many changes:
 * -- class Time_slice_conf: added methods set_*, changed argument list of
 *    method initial_data_cts.
 * -- class Tslice_dirac_max: added methods set_* and  hh_det_one().
 *
 * Revision 1.8  2004/04/05 21:21:51  e_gourgoulhon
 * class Time_slice_conf: added method initial_data_cts (computation of
 *  initial data from conformally thin sandwich method).
 * classes  Time_slice_conf and Tslice_dirac_max: added constructor as
 *  standard time slice of Minkowski spacetime.
 *
 * Revision 1.7  2004/04/01 16:09:01  j_novak
 * Trace of K_ij is now member of Time_slice (it was member of Time_slice_conf).
 * Added new methods for checking 3+1 Einstein equations (preliminary).
 *
 * Revision 1.6  2004/03/30 14:00:30  j_novak
 * New class Tslide_dirac_max (first version).
 *
 * Revision 1.5  2004/03/29 11:58:53  e_gourgoulhon
 * Many modif. to class Time_slice_conf.
 * Minor modif. to class Time_slice.
 *
 * Revision 1.4  2004/03/28 21:33:14  e_gourgoulhon
 * Constructor  Time_slice::Time_slice(int depth_in) declared "explicit".
 *
 * Revision 1.3  2004/03/28 21:27:57  e_gourgoulhon
 * Class Time_slice: - renamed the Evolution_std with suffix "_evol".
 *                   - added protected constructor for derived classes
 * Added class Time_slice_conf.
 *
 * Revision 1.2  2004/03/26 13:33:02  j_novak
 * New methods for accessing/updating members (nn(), beta(), gam_uu(), k_uu(), ...)
 *
 * Revision 1.1  2004/03/24 14:56:18  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

#include "star_rot_dirac.h"
#include "evolution.h"

                    //---------------------------//
                    //      class Time_slice     //
                    //---------------------------//

/**
 * Spacelike time slice of a 3+1 spacetime.
 * \ingroup (evol)
 * 
 */
class Time_slice {

    // Data : 
    // -----
    protected:
        /// Number of stored time slices
        int depth ; 
        
        /** Order of the finite-differences scheme for 
	 * the computation of time derivatives.
	 *
	 * This order is not constant and can be adjusted \e via
	 * \c set_scheme_order() .
	 */
        int scheme_order ; 
        
        /// Time step index of the latest slice
        int jtime ; 
        
        /// Time label of each slice
	Evolution_std<double> the_time ;
        
        /** Values at successive time steps of the covariant components of 
         * the induced metric \f$ \gamma_{ij} \f$
         */
	mutable Evolution_std<Sym_tensor> gam_dd_evol ; 
        
        /** Values at successive time steps of the contravariant components 
         * of the induced metric \f$ \gamma^{ij} \f$
         */        
	mutable Evolution_std<Sym_tensor> gam_uu_evol ; 

        /** Values at successive time steps of the covariant components 
         * of the extrinsic curvature tensor \f$ K_{ij} \f$
         */
	mutable Evolution_std<Sym_tensor> k_dd_evol ; 

        /** Values at successive time steps of the contravariant components 
         * of the extrinsic curvature tensor \f$ K^{ij} \f$
         */
	mutable Evolution_std<Sym_tensor> k_uu_evol ; 

        /// Values at successive time steps of the lapse function \e N 
	mutable Evolution_std<Scalar> n_evol ; 
        
        /// Values at successive time steps of the shift vector \f$ \beta^i \f$
	mutable Evolution_std<Vector> beta_evol ; 
        
        /** Values at successive time steps of the trace \e K of the 
         *  extrinsic curvature
         */        
	mutable Evolution_std<Scalar> trk_evol ; 

        /** ADM mass at each time step, since the creation of the slice.
         * At a given time step \c j, \c adm_mass_evol[j] is a 1-D \c Tbl
         *  of size the number \c nz of domains, containing the "ADM mass" 
         *  evaluated at the outer boundary of each domain. The true ADM
         *  mass is thus the last value, i.e. \c adm_mass_evol[j](nz-1). 
         *
         */
        mutable Evolution_full<Tbl> adm_mass_evol ; 

    // Derived data : 
    // ------------
    protected:
        /// Pointer on the induced metric at the current time step (\c jtime) 
	mutable Metric* p_gamma ;   

    // Constructors - Destructor
    // -------------------------
    public:
    
	/** General constructor (Hamiltonian-like). 
	 *
	 *  @param lapse_in lapse function \e N
	 *  @param shift_in shift vector
	 *  @param gamma_in induced metric (covariant or contravariant components) 
	 *  @param kk_in extrinsic curvature (covariant or contravariant components)
	 *  @param depth_in  number of stored time slices; this parameter is used
	 *                   to set the \c scheme_order member with \c scheme_order
	 *                   = \c depth_in - 1. \c scheme_order can be changed 
	 *                   afterwards by the method \c set_scheme_order(int).
	 */
	Time_slice(const Scalar& lapse_in, const Vector& shift_in,
		   const Sym_tensor& gamma_in, const Sym_tensor& kk_in,
		   int depth_in = 3) ; 
    
	/** General constructor (Lagrangian-like). 
	 *
	 *  @param lapse_in lapse function \e N
	 *  @param shift_in shift vector
	 *  @param gamma_in induced metric (covariant or contravariant components) 
	 *          at various time steps; note that the \c scheme_order member 
	 *          is set to \c gamma_in.size - 1. \c scheme_order can be changed 
	 *                   afterwards by the method \c set_scheme_order(int).
	 */
	Time_slice(const Scalar& lapse_in, const Vector& shift_in,
		   const Evolution_std<Sym_tensor>& gamma_in) ; 
    
	/** Constructor as standard time slice of flat spacetime (Minkowski). 
	 *
	 *  @param mp Mapping on which the various Lorene fields will be constructed
	 *  @param triad vector basis with respect to which the various tensor
	 *      components will be defined
	 *  @param depth_in  number of stored time slices; this parameter is used
	 *                   to set the \c scheme_order member with \c scheme_order
	 *                   = \c depth_in - 1. \c scheme_order can be changed 
	 *                   afterwards by the method \c set_scheme_order(int).
	 */
	Time_slice(const Map& mp, const Base_vect& triad, int depth_in = 3) ; 
	
	/** Constructor from binary file.
	 *
	 *  The binary file must have been created by method \c save.
	 *
	 *  @param mp Mapping on which the various Lorene fields will be constructed
	 *  @param triad vector basis with respect to which the various tensor
	 *      components will be defined
	 *  @param fich file containing the saved \c Time_slice
	 *  @param partial_read indicates whether the full object must
	 *      be read in file or whether the final construction is
	 *      devoted to a constructor of a derived class
	 *  @param depth_in  number of stored time slices; the given value must
	 *      coincide with that stored in the file.
	 */
	Time_slice(const Map& mp, const Base_vect& triad, FILE* fich, 
		   bool partial_read, int depth_in = 3) ; 
	
	Time_slice(const Time_slice& ) ;   ///< Copy constructor
	
    protected:
	/** Special constructor for derived classes.
	 *
	 */
	explicit Time_slice(int depth_in) ; 
	
    public:
	virtual ~Time_slice() ;			///< Destructor
 

    // Memory management
    // -----------------
    protected:
	    
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Time_slice
	void operator=(const Time_slice&) ;

	/// Sets the order of the finite-differences scheme.
	void set_scheme_order(int ord) { 
	    assert ((0<= ord)&&(ord < 4)) ;
	  scheme_order = ord ; } ; 
	
    // Accessors
    // ---------
    public:

	/// Gets the order of the finite-differences scheme.
	int get_scheme_order() const { return scheme_order ; } ;
        
        /// Gets the latest value of time step index
        int get_latest_j() const {return jtime; } ;
        
        /// Gets the time coordinate \e t at successive time steps
        const Evolution_std<double>& get_time() const {return the_time; } ; 
	
	/// Lapse function \e N at the current time step (\c jtime )
	virtual const Scalar& nn() const ;
	
	/// shift vector \f$ \beta^i \f$ at the current time step (\c jtime )
	virtual const Vector& beta() const ;
	
	/// Induced metric \f$ \mathbf{\gamma} \f$ at the current time step (\c jtime )
	const Metric& gam() const ;
	
	/** Induced metric (covariant components \f$ \gamma_{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& gam_dd() const ;
	
	/** Induced metric (contravariant components \f$ \gamma^{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& gam_uu() const ;
	
	/** Extrinsic curvature tensor (covariant components \f$ K_{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& k_dd() const ;
	
	/** Extrinsic curvature tensor (contravariant components \f$ K^{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& k_uu() const ;
	
        /** Trace \e K of the extrinsic curvature 
         *  at the current time step (\c jtime )
         */        
        virtual const Scalar& trk() const ; 
        
	
    // Computational functions
    // -----------------------
    public:
	/** 
	 *  Checks the level at which the hamiltonian constraint is verified.
	 * 
	 * \f[
	 * R + K^2 - K_{ij}K^{ij} = 16\pi E 
	 * \f]
	 * @param energy_density : a pointer on the energy density \e E
	 * measured by the Eulerian observer of 4-velocity 
	 * \f$\mbox{\boldmath{$n $}}\f$ ; if this
	 * is the null pointer, it is assumed that \e E = 0 (vacuum).
	 * @param ost : output stream for a formatted output of the result
	 * @return Tbl  of size the number of domains containing the 
	 * absolute ( if \e E = 0 ) or the relative (in presence of matter)
	 * error in max version.
	 */
	 Tbl check_hamiltonian_constraint(const Scalar* energy_density = 0x0,
					  ostream& ost = cout, bool verb=true) const ;
	
	/** 
	 *  Checks the level at which the momentum constraints are verified.
	 * 
	 * \f[
	 * D_j K_i^{\ j} - D_i K = 8 \pi J_i
	 * \f]
	 * @param momentum_density : a pointer on the momentum density 
	 * \f$ J_i \f$ measured by the Eulerian observer of 4-velocity 
	 * \f$\mbox{\boldmath{$n $}}\f$ ; if this is the null pointer, 
	 * it is assumed that\f$ J_i \f$ = 0 (vacuum).
	 * @param ost : output stream for a formatted output of the result
	 * @return Tbl 2D of size the number of domains times 3 
	 * (components)containing the absolute ( if \f$ J_i \f$ = 0 ) 
	 * or the relative (in presence of matter) error in max version.
	 */
	 Tbl check_momentum_constraint(const Vector* momentum_density = 0x0,
				       ostream& ost = cout, bool verb=true) const ;
	
	/** 
	 *  Checks the level at which the dynamical equations are verified.
	 * 
	 * \f[
	 * 	\frac{\partial K_{ij}}{\partial t} - 
	 * \pounds_{\mbox{\boldmath{$\beta $}}} K_{ij} =  - D_i D_j N 
	 *	+ N \left[ R_{ij} - 2 K_{ik} K^k_{\ j} + K K_{ij} 
	 *	+ 4\pi \left( (S-E)\gamma_{ij} - 2 S_{ij} \right) 
	 *	\right]
	 * \f]
	 * @param strain_tensor : a pointer on the strain_tensor 
	 * \f$ S_{ij} \f$ measured by the Eulerian observer of 4-velocity 
	 * \f$\mbox{\boldmath{$n $}}\f$ ; if this is the null pointer, 
	 * it is assumed that \f$ S_{ij} \f$ = 0 (vacuum).
	 * @param energy_density : a pointer on the energy density \e E
	 * (see \c check_hamiltonian_constraint)
	 * @param ost : output stream for a formatted output of the result
	 * @return Tbl 3D of size the number of domains times 3 times 3
	 * (corresponding to the rank-2 tensor, with the symmetry in the 
	 * components) containing the absolute ( if  \f$ J_i \f$ = 0 ) or 
	 * the relative (in presence of matter) error in max version.
	 */
	 Tbl check_dynamical_equations(const Sym_tensor* strain_tensor = 0x0,
				       const Scalar* energy_density = 0x0,
				       ostream& ost = cout, bool verb=true) const  ; 
	
        /** Returns the ADM mass (geometrical units) at the current step.
         * Moreover this method updates \c adm_mass_evol if
         * necessary. 
         */
        virtual double adm_mass() const ; 
        
    // Outputs
    // -------
    protected:
	/// Operator >> (virtual function called by the operator<<). 
	virtual ostream& operator>>(ostream& ) const ; 
	
    /// Display
	friend ostream& operator<<(ostream& , const Time_slice& ) ;	

    public:
	/** Saves in a binary file.
	 *  The saved data is sufficient to restore the whole time slice
	 *  via the constructor from file. 
	 *  @param rootname root for the file name; the current time step index
	 *      will be appended to it. 
	 */
	void save(const char* rootname) const ; 
    
    protected:
	/** Total or partial saves in a binary file.
	 *  This protected method is to be called either from public method 
	 *  \c save or from method \c sauve of a derived class. 
	 *  
	 *  @param fich binary file 
	 *  @param partial_save indicates whether the whole object must be
	 *      saved.
	 */
	virtual void sauve(FILE* fich, bool partial_save) const ; 

};

ostream& operator<<(ostream& , const Time_slice& ) ;	



                    //---------------------------//
                    //   class Time_slice_conf   //
                    //---------------------------//

/**
 * Spacelike time slice of a 3+1 spacetime with conformal decomposition.
 * \ingroup (evol)
 * 
 */
class Time_slice_conf : public Time_slice {

    // Data : 
    // -----
    protected: 
    
        /** Pointer on the flat metric \f$ f_{ij} \f$ with respect to
         * which the conformal decomposition is performed
         */
        const Metric_flat& ff ;  

        /** Values at successive time steps of the conformal factor 
         * \f$ \Psi \f$ relating the
         * physical metric \f$ \gamma_{ij} \f$ to the conformal one:
         * \f$ \gamma_{ij} = \Psi^4 \tilde\gamma_{ij} \f$.
         * \f$ \Psi \f$ is defined by
         *  \f[ \Psi := \left( \frac{\det\gamma_{ij}}{\det f_{ij}} 
         *      \right) ^{1/12} \f] 
         */
	mutable Evolution_std<Scalar> psi_evol ; 
        
        /** Values at successive time steps of the factor 
         * \f$ N\Psi \f$.
         */
	mutable Evolution_std<Scalar> npsi_evol ; 
        
        
        /** Values at successive time steps of the components \f$ h^{ij} \f$. 
         * It is the deviation of the conformal metric \f$ \tilde\gamma^{ij} \f$ from 
         * the flat metric \f$ f^{ij} \f$: 
         * \f$\tilde\gamma^{ij} = f^{ij} + h^{ij} \f$.
         */        
	mutable Evolution_std<Sym_tensor> hh_evol ; 

        /** Values at successive time steps of the components \f$ \hat{A}^{ij} \f$.
         * It is the conformal representation of the traceless part
         * of the extrinsic curvature:
         * \f$ \hat{A}^{ij} = \Psi^{10} \left( K^{ij} - \frac{1}{3} K \gamma^{ij} \right) \f$.
	 * One can uniquely (up to the boundary conditions) define the decomposition:
	 * \f$ \hat{A}^{ij} = {\cal D}^i X^j + {\cal D}^j X^i - \frac{2}{3} 
	 * {\cal D}_k X^k f^{ij} + \hat{A}^{ij}_{TT} \f$, where \f$ X^i \f$ 
	 * represents the longitudinal part and \f$ \hat{A}^{ij}_{TT} \f$ is 
	 * the transverse-traceless part.
         */        
	mutable Evolution_std<Sym_tensor> hata_evol ; 

	/** Potential \e A associated with the symmetric tensor 
	 * \f$ \hat{A}^{ij}_{TT} \f$. (see documentation of \c Sym_tensor::p_aaa).
	 */
	mutable Evolution_std<Scalar> A_hata_evol ;
        
	/** Potential \f$ \tilde{B} \f$ associated with the symmetric tensor 
	 * \f$ \hat{A}^{ij}_{TT} \f$. (see documentation of \c Sym_tensor::p_tilde_b).
	 */
	mutable Evolution_std<Scalar> B_hata_evol ;
        
    // Derived data : 
    // ------------
    protected:
        /** Pointer on the conformal metric \f$ \tilde\gamma_{ij} \f$
         * at the current time step (\c jtime)
         */
	    mutable Metric* p_tgamma ; 
        
        /// Pointer on the factor \f$ \Psi^4 \f$ at the current time step (\c jtime)
	    mutable Scalar* p_psi4 ; 
        
        /// Pointer on the logarithm of \f$ \Psi \f$ at the current time step (\c jtime)
	    mutable Scalar* p_ln_psi ; 
        
        /** Pointer on the vector \f$ H^i = {\cal D}_j \tilde\gamma^{ij} \f$ 
         * (which vanishes in Dirac gauge), at the current time step (\c jtime).
         */
        mutable Vector* p_hdirac ; 

	/** Pointer on the vector \f$ X^i \f$ representing the longitudinal 
	 * part of \f$ \hat{A}^{ij} \f$.
	 * (see the documentation of \c hata_evol)
	 */
	mutable Vector* p_vec_X ;

    // Constructors - Destructor
    // -------------------------
    public:
    
    /** Constructor from conformal decomposition.
     *
     *  @param lapse_in lapse function \e N
     *  @param shift_in shift vector
     *  @param ff_in reference flat metric with respect to which the
     *           conformal decomposition is performed
     *  @param psi_in conformal factor \f$\Psi\f$ relating the
     *       physical metric \f$ \gamma_{ij} \f$ to the conformal one:
     *      \f$ \gamma_{ij} = \Psi^4 \tilde\gamma_{ij} \f$
     *  @param hh_in deviation \f$ h^{ij} \f$
     *      of the conformal metric \f$ \tilde\gamma^{ij} \f$ from 
     *      the flat metric \f$ f^{ij} \f$: 
     *      \f$\tilde\gamma^{ij} = f^{ij} + h^{ij} \f$.
     *      \f$ h^{ij} \f$ must be such that 
     *      \f$\det\tilde\gamma^{ij} = f^{-1} \f$.
     *  @param hata_in conformal representation \f$ A^{ij} \f$
     *      of the traceless part of the extrinsic curvature:
     *   \f$ \hat{A}^{ij} = \Psi^{10} \left( K^{ij} - \frac{1}{3} K \gamma^{ij} \right) \f$
     *  @param trk_in trace \e K of the extrinsic curvature 
     *  @param depth_in  number of stored time slices; this parameter is used
     *                   to set the \c scheme_order member with \c scheme_order
     *                   = \c depth_in - 1. \c scheme_order can be changed 
     *                   afterwards by the method \c set_scheme_order(int).
     */
    Time_slice_conf(const Scalar& lapse_in, const Vector& shift_in,
            const Metric_flat& ff_in, const Scalar& psi_in, 
            const Sym_tensor& hh_in, const Sym_tensor& hata_in, 
            const Scalar& trk_in, int depth_in = 3) ; 
    
    
    /** Constructor from physical metric.
     *  The conformal decomposition is performed by the constructor. 
     *
     *  @param lapse_in lapse function \e N
     *  @param shift_in shift vector
     *  @param gamma_in induced metric (covariant or contravariant components) 
     *  @param kk_in extrinsic curvature (covariant or contravariant components)
     *  @param ff_in reference flat metric with respect to which the
     *           conformal decomposition is performed
     *  @param depth_in  number of stored time slices; this parameter is used
     *                   to set the \c scheme_order member with \c scheme_order
     *                   = \c depth_in - 1. \c scheme_order can be changed 
     *                   afterwards by the method \c set_scheme_order(int).
     */
    Time_slice_conf(const Scalar& lapse_in, const Vector& shift_in,
               const Sym_tensor& gamma_in, const Sym_tensor& kk_in,
               const Metric_flat& ff_in, int depth_in = 3) ; 
               
    /** Constructor as standard time slice of flat spacetime (Minkowski). 
     *
     *  @param mp Mapping on which the various Lorene fields will be constructed
     *  @param triad vector basis with respect to which the various tensor
     *      components will be defined
     *  @param ff_in reference flat metric with respect to which the
     *           conformal decomposition is performed
     *  @param depth_in  number of stored time slices; this parameter is used
     *                   to set the \c scheme_order member with \c scheme_order
     *                   = \c depth_in - 1. \c scheme_order can be changed 
     *                   afterwards by the method \c set_scheme_order(int).
     */
    Time_slice_conf(const Map& mp, const Base_vect& triad, 
                    const Metric_flat& ff_in, int depth_in = 3) ; 
    
    /** Constructor from binary file.
     *
     *  The binary file must have been created by method \c save.
     *
     *  @param mp Mapping on which the various Lorene fields will be constructed
     *  @param triad vector basis with respect to which the various tensor
     *      components will be defined
     *  @param ff_in reference flat metric with respect to which the
     *           conformal decomposition is performed
     *  @param fich file containing the saved \c Time_slice_conf
     *  @param partial_read indicates whether the full object must
     *      be read in file or whether the final construction is
     *      devoted to a constructor of a derived class
     *  @param depth_in  number of stored time slices; the given must
     *      coincide with that stored in the file.
     */
    Time_slice_conf(const Map& mp, const Base_vect& triad, 
                    const Metric_flat& ff_in, FILE* fich, 
                    bool partial_read, int depth_in = 3) ; 

    Time_slice_conf(const Time_slice_conf& ) ;	///< Copy constructor

    virtual ~Time_slice_conf() ;			///< Destructor
 

    // Memory management
    // -----------------
    protected:
	    
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Time_slice_conf
	void operator=(const Time_slice_conf&) ;

	/// Assignment to a \c Time_slice
	void operator=(const Time_slice&) ;
        
        /** Sets the conformal factor \f$ \Psi \f$ relating the
         * physical metric \f$ \gamma_{ij} \f$ to the conformal one:
         * \f$ \gamma_{ij} = \Psi^4 \tilde\gamma_{ij} \f$. 
         * \f$ \Psi \f$ is defined by
         *  \f[ \Psi := \left( \frac{\det\gamma_{ij}}{\det f_{ij}} 
         *      \right) ^{1/12} \f] 
         * Sets the value at the current time step (\c jtime ) and
         * deletes the value of \f$N\Psi\f$.
         *
         */
        virtual void set_psi_del_npsi(const Scalar& psi_in) ; 
        
        /** Sets the conformal factor \f$ \Psi \f$ relating the
         * physical metric \f$ \gamma_{ij} \f$ to the conformal one:
         * \f$ \gamma_{ij} = \Psi^4 \tilde\gamma_{ij} \f$. 
         * \f$ \Psi \f$ is defined by
         *  \f[ \Psi := \left( \frac{\det\gamma_{ij}}{\det f_{ij}} 
         *      \right) ^{1/12} \f] 
         * Sets the value at the current time step (\c jtime ) and
         * deletes the value of \e N.
         *
         */
        virtual void set_psi_del_n(const Scalar& psi_in) ; 
        
        /** Sets the factor \f$ N\Psi \f$ at the 
         *  current time step (\c jtime ) and deletes the value
         *  of \f$\Psi\f$.
         */
        virtual void set_npsi_del_psi(const Scalar& npsi_in) ; 
        
        /** Sets the factor \f$ N\Psi \f$ at the 
         *  current time step (\c jtime ) and deletes the value
         *  of \e N.
         */
        virtual void set_npsi_del_n(const Scalar& npsi_in) ; 
        
        /** Sets the deviation \f$ h^{ij} \f$ 
         * of the conformal metric \f$ \tilde\gamma^{ij} \f$ from 
         * the flat metric \f$ f^{ij} \f$: 
         * \f$\tilde\gamma^{ij} = f^{ij} + h^{ij} \f$.
         *      \f$ h^{ij} \f$ must be such that 
         *      \f$\det\tilde\gamma^{ij} = f^{-1} \f$.
         * Sets the value at the current time step (\c jtime ).
         */        
        virtual void set_hh(const Sym_tensor& hh_in) ; 

        /** Sets the conformal representation \f$ \hat{A}{ij} \f$ of the traceless part
         * of the extrinsic curvature:
         * \f$ \hat{A}^{ij} = \Psi^{10} \left( K^{ij} - \frac{1}{3} K \gamma^{ij} \right) \f$.
         * Sets the value at the current time step (\c jtime ), and updates the
	 * potentials \c A_hata_evol, \c B_hata_evol and \c p_vec_X accordingly.
         */        
        virtual void set_hata(const Sym_tensor& hata_in) ; 

	/** Sets the TT part of \f$ \hat{A}^{ij} \f$ (see member \c hata_evol ).
	 * Sets the value at current time-step (\c jtime ) and updates the potentials
	 * \e A and \f$ \tilde{B} \f$. */
	virtual void set_hata_TT(const Sym_tensor_tt& hata_tt) ;

        /** Sets the conformal representation \f$ \hat{A}{ij} \f$ of the traceless part
         * of the extrinsic curvature from its potentials \e A, \f$ \tilde{B} \f$
	 * and \f$ X^i \f$.
         * These potentials must be up-to-date. It sets the value at the current 
	 * time step (\c jtime ). 
         */        
	virtual void set_hata_from_XAB(Param* par_bc=0x0, Param* par_mat=0x0) ;

    // Accessors
    // ---------
    public:

        // Virtual functions from base class Time_slice:
        // ---------------------------------------------

	/// Lapse function \e N at the current time step (\c jtime )
	virtual const Scalar& nn() const ;
	
	/** Induced metric (covariant components \f$ \gamma_{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& gam_dd() const ;
	
	/** Induced metric (contravariant components \f$ \gamma^{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& gam_uu() const ;
	
	/** Extrinsic curvature tensor (covariant components \f$ K_{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& k_dd() const ;
	
	/** Extrinsic curvature tensor (contravariant components \f$ K^{ij} \f$) 
	 * at the current time step (\c jtime )
	 */
	virtual const Sym_tensor& k_uu() const ;
	
        // Virtual functions from this class:
        // ----------------------------------

	/** Returns the potential \e A of \f$ \hat{A}^{ij} \f$.
	 * See the documentation of \c Sym_tensor for details.
	 * Returns the value at the current time step (\c jtime ).
	 */
	virtual const Scalar& A_hata() const ; 
	
	/** Returns the potential \f$\tilde{B}\f$ of \f$ \hat{A}^{ij} \f$.
	 * See the documentation of \c Sym_tensor_tt for details.
	 * Returns the value at the current time step (\c jtime ).
	 */
	virtual const Scalar& B_hata() const ;
    
        /** Conformal factor \f$ \Psi \f$ relating the
         * physical metric \f$ \gamma_{ij} \f$ to the conformal one:
         * \f$ \gamma_{ij} = \Psi^4 \tilde\gamma_{ij} \f$. 
         * \f$ \Psi \f$ is defined by
         *  \f[ \Psi := \left( \frac{\det\gamma_{ij}}{\det f_{ij}} \right) ^{1/12} \f] 
         * Returns the value at the current time step (\c jtime ).
         */
        virtual const Scalar& psi() const ; 
        
        /// Factor \f$ \Psi^4 \f$ at the current time step (\c jtime ).
	    const Scalar& psi4() const ; 
        
        /// Logarithm of \f$ \Psi \f$ at the current time step (\c jtime ).
	    const Scalar& ln_psi() const ; 
        
        /** Factor \f$ N\Psi \f$ at the current time step (\c jtime ).
         */
        virtual const Scalar& npsi() const ; 
        
        /** Conformal metric 
         * \f$ \tilde\gamma_{ij} = \Psi^{-4} \gamma_{ij} \f$
         * Returns the value at the current time step (\c jtime ).
         */        
        virtual const Metric& tgam() const ; 

        /** Deviation \f$ h^{ij} \f$ 
         * of the conformal metric \f$ \tilde\gamma^{ij} \f$ from 
         * the flat metric \f$ f^{ij} \f$: 
         * \f$\tilde\gamma^{ij} = f^{ij} + h^{ij} \f$.
         * Returns the value at the current time step (\c jtime ).
         */        
        virtual const Sym_tensor& hh(Param* = 0x0, Param* = 0x0) const ; 

        /** Conformal representation \f$ \hat{A}^{ij} \f$ of the traceless part
         * of the extrinsic curvature:
         * \f$ \hat{A}^{ij} = \Psi^{10} \left( K^{ij} - \frac{1}{3} K \gamma^{ij} \right) \f$.
         * Returns the value at the current time step (\c jtime ).
         */        
        virtual const Sym_tensor& hata() const ; 

        /** Conformal representation \f$ A^{ij} \f$ of the traceless part
         * of the extrinsic curvature:
         * \f$ A^{ij} = \Psi^4 \left( K^{ij} - \frac{1}{3} K \gamma^{ij} \right) \f$.
         * Returns the value at the current time step (\c jtime ).
         */        
        virtual Sym_tensor aa() const ; 

        /** Trace \e K of the extrinsic curvature 
         *  at the current time step (\c jtime )
         */        
        virtual const Scalar& trk() const ; 
        
        /** Vector \f$ H^i = {\cal D}_j \tilde\gamma^{ij} \f$ 
         * which vanishes in Dirac gauge.
         */
        virtual const Vector& hdirac() const ;

	/** Vector \f$ X^i \f$ representing the longitudinal part of 
	 * \f$ \hat{A}^{ij} \f$.(see the documentation of \c hata_evol)
	 */
        virtual const Vector& vec_X(int method_poisson=6) const ;
        
    // Computational methods
    // ---------------------
    public:
	/** Computes the vector \f$ X^i \f$ from the conformally-rescaled 
	 * momentum \f$ \hat{S}^i = \Psi^6 S^i \f$, using the momentum constraint.
	 */
	void compute_X_from_momentum_constraint
	    (const Vector& hat_S, const Sym_tensor_tt& hata_tt, 
	      int iter_max = 200, double precis = 1.e-12,
	     double relax = 0.8, int methode_poisson = 6) ;

	/** Sets the potentials \e A and \f$\tilde{B}\f$
         * of the TT part \f$ \hat{A}^{ij}
	 * (see the documentation of \c Sym_tensor for details).
         * Sets the value at the current time step (\c jtime ).
	 */
	virtual void set_AB_hata(const Scalar& A_in, const Scalar& B_in) ; 

        /** Computes valid initial data by solving the constraint 
         *  equations in the conformal thin-sandwich approach.
         *
         *  @param uu value of 
         *    \f$ {\tilde u}^{ij} = \partial h^{ij} /\partial t \f$ 
         *                  (freely specifiable data).
         *  This quantity must be trace-free with respect to the conformal
         *  metric \f$\tilde\gamma_{ij}\f$, reflecting the unimodular 
         *  character of \f$\tilde\gamma_{ij}\f$.
         *  @param trk_in value of \f$ K = K_i^{\ i} \f$ 
         *      (freely specifiable data)
         *  @param trk_point value of \f$ \partial K / \partial t \f$ 
         *      (freely specifiable data)
         *  @param pdt time step, to be used in order to fill \c depth
         *      slices
         *  @param precis convergence threshold required to stop the 
         *          iteration
         *  @param method_poisson_vect method to be used for solving 
         *      vector Poisson equation (for the shift), see 
         *      \c Vector::poisson(double, const Metric_flat&, int) const.  
         *  @param graph_device  name of type of graphical device: 0x0 
         *      (default value) will result in interactive choice; 
         *      "/xwin" in X-Window display and "/n" in no output.
         *  @param ener_dens matter energy density \e E as measured by the 
         *      Eulerian observer; this quantity is passed as a pointer,
         *      the null value of which (default) meaning \e E=0.
         *  @param mom_dens matter momentum density \e J as measured by the 
         *      Eulerian observer; this quantity is passed as a pointer,
         *      the null value of which (default) meaning \e J=0.
         *  @param trace_stress trace of the matter stress \e S as measured 
         *      by the Eulerian observer; this quantity is passed as a pointer,
         *      the null value of which (default) meaning \e S=0.
         */
         virtual void initial_data_cts(const Sym_tensor& uu, const Scalar& trk_in, 
                const Scalar& trk_point, double pdt, double precis = 1.e-12,
                int method_poisson_vect = 6, const char* graph_device = 0x0, 
                const Scalar* ener_dens = 0x0, const Vector* mom_dens = 0x0, 
                const Scalar* trace_stress = 0x0 ) ; 
        
        /** Returns the ADM mass (geometrical units) at the current step.
         * Moreover this method updates \c adm_mass_evol if
         * necessary. 
         */
        virtual double adm_mass() const ; 
        
        /** Checks the \f$\frac{\partial}{\partial t} \ln\Psi \f$ relation.
         *  @param tlnpsi_dot [output] maximun value in each domain of 
         *      \f$ \left| \frac{\partial}{\partial t} \ln\Psi \right| \f$
         *  @param tdiff [output] maximum value in each domain of  \f$ \left| 
         *          \frac{\partial}{\partial t} \ln\Psi -
         *          \beta^i {\cal D}_i \ln \Psi - \frac{1}{6} ( 
         *          {\cal D}_i \beta^i - N K) \right| \f$
         *  @param tdiff_rel [output] relative error on the above relation
         *      in each domain.
         *
         */
        void check_psi_dot(Tbl& tlnpsi_dot, Tbl& tdiff, Tbl& tdiff_rel) const ; 
        
    // Outputs
    // -------
    protected:
    /// Operator >> (virtual function called by the operator<<). 
    virtual ostream& operator>>(ostream& ) const ; 
	
    /** Total or partial saves in a binary file.
     *  This protected method is to be called either from public method 
     *  \c save or from method \c sauve of a derived class. 
     *  
     *  @param fich binary file 
     *  @param partial_save indicates whether the whole object must be
     *      saved.
     */
    virtual void sauve(FILE* fich, bool partial_save) const ; 

} ;	
                    //----------------------------//
                    //   class Tslice_dirac_max   //
                    //----------------------------//

/**
 * Spacelike time slice of a 3+1 spacetime with conformal decomposition
 * in the maximal slicing and Dirac gauge.
 * \ingroup (evol)
 * 
 */
class Tslice_dirac_max : public Time_slice_conf {

  // Data : 
  // -----
 protected: 
  /** The \e A potential of \f$ \bar{h}^{ij} \f$.
   *
   * (see the documentation of \c Sym_tensor::p_aaa for details).
   */
  mutable Evolution_std<Scalar> A_hh_evol ;
  
  /** The \f$\tilde{B} \f$ potential of \f$ \bar{h}^{ij} \f$.
   *
   * (see the documentation of \c Sym_tensor::p_tilde_b for details).
   */
  mutable Evolution_std<Scalar> B_hh_evol ;

  /** The \e A potential of the source of equation for \f$ \bar{h}^{ij} \f$.
   *
   * (see the documentation of \c Sym_tensor::p_aaa for details).
   */
  mutable Evolution_std<Scalar> source_A_hh_evol ;
  
  /** The \f$\tilde{B} \f$ potential of the source of equation for \f$ \bar{h}^{ij} \f$.
   *
   * (see the documentation of \c Sym_tensor::p_tilde_b for details).
   */
  mutable Evolution_std<Scalar> source_B_hh_evol ;

  /** The potential \e A of the source of equation for \f$ \hat{A}^{ij} \f$.
   *
   * (see the documentation of \c Sym_tensor::p_aaa for details).
   */
  mutable Evolution_std<Scalar> source_A_hata_evol ;
  
  /** The potential \f$\tilde{B}\f$ of the source of equation for \f$ \hat{A}^{ij} \f$.
   *
   * (see the documentation of \c Sym_tensor::p_tilde_b for details).
   */
  mutable Evolution_std<Scalar> source_B_hata_evol ;

  /// The trace, with respect to the flat metric \c ff , of \f$ h^{ij} \f$.
  mutable Evolution_std<Scalar> trh_evol ;


    // Constructors - Destructor
    // -------------------------
    public:
    /** Constructor from conformal decomposition.
     *
     *  @param lapse_in lapse function \e N
     *  @param shift_in shift vector
     *  @param ff_in reference flat metric with respect to which the
     *           conformal decomposition is performed
     *  @param psi_in conformal factor \f$\Psi\f$ relating the
     *       physical metric \f$ \gamma_{ij} \f$ to the conformal one:
     *      \f$ \gamma_{ij} = \Psi^4 \tilde\gamma_{ij} \f$
     *  @param hh_in deviation \f$ h^{ij} \f$
     *      of the conformal metric \f$ \tilde\gamma^{ij} \f$ from 
     *      the flat metric \f$ f^{ij} \f$: 
     *      \f$\tilde\gamma^{ij} = f^{ij} + h^{ij} \f$.
     *      \f$ h^{ij} \f$ must be such that 
     *      \f$\det\tilde\gamma^{ij} = f^{-1} \f$.
     *  @param hata_in conformal representation \f$ \hat{A}^{ij} \f$
     *      of the traceless part of the extrinsic curvature:
     *      \f$ \hat{A}^{ij} = \Psi^{10} \left( K^{ij} - \frac{1}{3} K \gamma^{ij} \right) \f$,
     *       with \e K = 0 in the present case
     *  @param depth_in  number of stored time slices; this parameter is used
     *                   to set the \c scheme_order member with \c scheme_order
     *                   = \c depth_in - 1. \c scheme_order can be changed 
     *                   afterwards by the method \c set_scheme_order(int).
     */
	Tslice_dirac_max(const Scalar& lapse_in, const Vector& shift_in,
            const Metric_flat& ff_in, const Scalar& psi_in, 
            const Sym_tensor_trans& hh_in, const Sym_tensor& hata_in, 
            int depth_in = 3) ;	
	
    /** Constructor as standard time slice of flat spacetime (Minkowski). 
     *
     *  @param mp Mapping on which the various Lorene fields will be constructed
     *  @param triad vector basis with respect to which the various tensor
     *      components will be defined
     *  @param ff_in reference flat metric with respect to which the
     *           conformal decomposition is performed
     *  @param depth_in  number of stored time slices; this parameter is used
     *                   to set the \c scheme_order member with \c scheme_order
     *                   = \c depth_in - 1. \c scheme_order can be changed 
     *                   afterwards by the method \c set_scheme_order(int).
     */
    Tslice_dirac_max(const Map& mp, const Base_vect& triad, 
                     const Metric_flat& ff_in, int depth_in = 3) ; 
    
    /** Constructor from binary file.
     *
     *  The binary file must have been created by method \c save.
     *
     *  @param mp Mapping on which the various Lorene fields will be constructed
     *  @param triad vector basis with respect to which the various tensor
     *      components will be defined
     *  @param ff_in reference flat metric with respect to which the
     *           conformal decomposition is performed
     *  @param fich file containing the saved \c Tslice_dirac_max
     *  @param partial_read indicates whether the full object must
     *      be read in file or whether the final construction is
     *      devoted to a constructor of a derived class
     *  @param depth_in  number of stored time slices; the given must
     *      coincide with that stored in the file.
     */
    Tslice_dirac_max(const Map& mp, const Base_vect& triad, 
                     const Metric_flat& ff_in, FILE* fich, 
                     bool partial_read, int depth_in = 3) ; 

    /// Construnction of a stationary slice from a rotating star
    Tslice_dirac_max(const Star_rot_Dirac& star, double pdt, int depth_in = 3) ;

    Tslice_dirac_max(const Tslice_dirac_max& ) ;   ///< Copy constructor

    virtual ~Tslice_dirac_max() ;			///< Destructor
 

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Tslice_dirac_max
	void operator=(const Tslice_dirac_max&) ;	
	
        // Virtual functions from base class Time_slice_conf:
        // -------------------------------------------------

        /** Sets the deviation \f$ h^{ij} \f$ 
         * of the conformal metric \f$ \tilde\gamma^{ij} \f$ from 
         * the flat metric \f$ f^{ij} \f$: 
         * \f$\tilde\gamma^{ij} = f^{ij} + h^{ij} \f$.
         *      \f$ h^{ij} \f$ must be such that 
         *      \f$\det\tilde\gamma^{ij} = f^{-1} \f$.
         * Sets the value at the current time step (\c jtime ).
         */        
        virtual void set_hh(const Sym_tensor& hh_in) ; 

        /** Computes valid initial data by solving the constraint 
         *  equations in the conformal thin-sandwich approach.
         *
         *  @param uu value of 
         *    \f$ {\tilde u}^{ij} = \partial h^{ij} /\partial t \f$ 
         *                  (freely specifiable data).
         *  This quantity must be trace-free with respect to the conformal
         *  metric \f$\tilde\gamma_{ij}\f$, reflecting the unimodular 
         *  character of \f$\tilde\gamma_{ij}\f$.
         *  @param trk_in value of \f$ K = K_i^{\ i} \f$ 
         *      (freely specifiable data)
         *  @param trk_point value of \f$ \partial K / \partial t \f$ 
         *      (freely specifiable data)
         *  @param pdt time step, to be used in order to fill \c depth
         *      slices
         *  @param precis convergence threshold required to stop the 
         *          iteration
         *  @param method_poisson_vect method to be used for solving 
         *      vector Poisson equation (for the shift), see 
         *      \c Vector::poisson(double, const Metric_flat&, int) const.  
         *  @param graph_device  name of type of graphical device: 0x0 
         *      (default value) will result in interactive choice; 
         *      "/xwin" in X-Window display and "/n" in no output.
         *  @param ener_dens matter energy density \e E as measured by the 
         *      Eulerian observer; this quantity is passed as a pointer,
         *      the null value of which (default) meaning \e E=0.
         *  @param mom_dens matter momentum density \e J as measured by the 
         *      Eulerian observer; this quantity is passed as a pointer,
         *      the null value of which (default) meaning \e J=0.
         *  @param trace_stress trace of the matter stress \e S as measured 
         *      by the Eulerian observer; this quantity is passed as a pointer,
         *      the null value of which (default) meaning \e S=0.
         */
         virtual void initial_data_cts(const Sym_tensor& uu, const Scalar& trk_in, 
                const Scalar& trk_point, double pdt, double precis = 1.e-12,
                int method_poisson_vect = 6, const char* graph_device = 0x0, 
                const Scalar* ener_dens = 0x0, const Vector* mom_dens = 0x0, 
                const Scalar* trace_stress = 0x0 ) ; 
        
        // Virtual functions from this class:
        // ----------------------------------

	/** Sets the potentials \f$\chi \f$ and \f$\mu\f$
         * of the TT part \f$ \bar{h}^{ij} \f$ of \f$ h^{ij} \f$
	 * (see the documentation of \c Sym_tensor_tt for details).
         * The value of \f$ h^{ij} \f$ is then deduced from the
         * unimodularity condition on the conformal metric.
         * Sets the value at the current time step (\c jtime ).
	 */
	virtual void set_khi_mu(const Scalar& khi_in, const Scalar& mu_in) ; 

	/** Sets the potentials \e A and \f$\tilde{B}\f$
         * of the TT part \f$ \bar{h}^{ij} \f$ of \f$ h^{ij} \f$
	 * (see the documentation of \c Sym_tensor for details).
	 * \f$ h^{ij} \f$ is not modified.
         * Sets the value at the current time step (\c jtime ).
	 */
	virtual void set_AB_hh(const Scalar& A_in, const Scalar& B_in) ; 

	/** Sets the trace, with respect to the flat metric 
	 * \c ff , of \f$ h^{ij} \f$.
         * Sets the value at the current time step (\c jtime ).
         * Note that this method does not ensure that the conformal
         * metric is unimodular. 
	 */
	virtual void set_trh(const Scalar& trh_in) ;
    
	/** Solves the elliptic equation for the conformal factor $\Psi$
	 *  (Hamiltonian constraint).
	 *  @param ener_dens conformal matter energy density \f$ \hat{E} = \Psi^6 E \f$,
	 *         where \e E is measured by the Eulerian observer; this quantity is 
	 *         passed as a pointer, the null value of which (default) meaning \e E=0.
	 *  @return solution \f$\Psi_{\rm new}\f$ of the elliptic equation 
	 *  (flat Laplacian) for the lapse with the source computed from the
	 *  quantities at the current time step. 
	 *  
	 */
	virtual Scalar solve_psi(const Scalar* ener_dens=0x0) const ; 
        
	/** Solves the elliptic equation for \f$ N\Psi \f$ (maximal
	 *  slicing condition + Hamiltonian constraint)
	 *  @param ener_dens conformal matter energy density \f$ \hat{E} = \Psi^6 E \f$,
	 *         where \e E is measured by the Eulerian observer; this quantity is 
	 *         passed as a pointer, the null value of which (default) meaning \e E=0.
	 *  @param trace_stress trace of the conformal matter stress 
	 *         \f$ S^* = \Psi^6 S \f$,
	 *                      where \e S is measured by the Eulerian observer; 
	 *                      this quantity is passed as a pointer,
	 *                      the null value of which (default) meaning \e S=0.
	 *  @return solution \f$(N\Psi)_{\rm new}\f$ of the elliptic equation 
	 *  (flat Laplacian) for \f$ N\Psi \f$ with the source computed from the
	 *  quantities at the current time step. 
	 *  
	 */
	virtual Scalar solve_npsi(const Scalar* ener_dens=0x0,
				  const Scalar* trace_stress=0x0) const ; 
        
	/** Solves the elliptic equation for the shift vector \f$\beta^i\f$ 
	 *  from \f$ A^{ij} \f$ (Eq. (73) of Bonazzola et al. 2004).
	 *  @param method method to be used for solving 
	 *      vector Poisson equation (for the shift), see 
	 *      \c Vector::poisson(double, const Metric_flat&, int) const.  
	 *  @return solution \f$\beta^i_{\rm new}\f$ of the elliptic equation 
	 *  (flat vector Laplacian) for the shift with the source computed from the
	 *  quantities at the current time step. 
	 *  
	 */
	virtual Vector solve_beta(int method = 6) const ; 
        
    /** Time evolution by resolution of Einstein equations.
     *  
     *  @param pdt  time step \e dt.
     *  @param nb_time_steps  number of time steps for the evolution
     *  @param niter_elliptic  number of iterations if the resolution 
     *      of elliptic equations
     *  @param relax_elliptic  relaxation factor for the elliptic
     *      equations
     *  @param check_mod determines the frequency of check of the 
     *      constraint equations: they are checked every \c check_mod time step
     *  @param save_mod determines the frequency of writing to file
     *      the monotoring quantities: they are written to file every
     *      \c save_mod time step
     *  @param method method_poisson_vect to be used for solving 
     *      vector Poisson equation (for the shift), see 
     *      \c Vector::poisson(double, const Metric_flat&, int) const.  
     *  @param nopause  = 1 if no pause between each time step, 0 otherwise
     *  @param graph_device  name of type of graphical device: 0x0 
     *      (default value) will result in interactive choice; 
     *      "/xwin" in X-Window display and "/n" in no output.
     */
    void evolve(double pdt, int nb_time_steps, int niter_elliptic,
                double relax_elliptic, int check_mod, int save_mod,
                int method_poisson_vect = 6, int nopause = 1, 
                const char* graph_device = 0x0, bool verbose=true, 
		const Scalar* ener_euler = 0x0,
		const Vector* mom_euler = 0x0, const Scalar* s_euler = 0x0,
		const Sym_tensor* strain_euler = 0x0) ; 
        
    /** Returns the ADM mass at (geometrical units) the current step.
     * Moreover this method updates \c adm_mass_evol if
     * necessary. 
     */
    virtual double adm_mass() const ; 
        
    protected:
    /** Computes the sources \c source_A_XXX_evol and \c source_B_XXX_evol ,
     * for the solution of the evolution equation for \f$ h^{ij} \f$ and
     * \f$ \hat{A}^{ij} \f$.
     *  @param strain_tensor [input] : a pointer on the strain_tensor 
     *      \f$ S^{ij} \f$ measured by the Eulerian observer of 4-velocity 
     *      \f$\mbox{\boldmath{$n $}}\f$ ; if this is the null pointer
     *      (default value), it is assumed that \f$ S_{ij} \f$ = 0 (vacuum).  
     */
    void compute_sources(const Sym_tensor* strain_tensor = 0x0) const ;

    /// Copy the sources \c source_A_XXX_evol and \c source_B_XXX_evol to all time-steps.
    void initialize_sources_copy() const ;

    /** Computes \f$ h^{ij} \f$ from the values of \e A and 
     * \f$\tilde{B}\f$ and using the condition 
     * \f$\det\tilde\gamma^{ij} = \det f^{ij} \f$, which fixes the
     * trace of \f$ h^{ij} \f$.
     * @param j time step at which the computation of \f$ h^{ij} \f$
     *      is required.
     */
    void hh_det_one(int j, Param* par_bc = 0x0, Param* par_mat = 0x0) const ; 
    
    /** Computes \f$ h^{ij} \f$ from the TT part using the condition 
     * \f$\det\tilde\gamma^{ij} = \det f^{ij} \f$, which fixes the
     * trace of \f$ h^{ij} \f$.
     * @param hijtt : the TT part.
     */
    void hh_det_one(const Sym_tensor_tt& hijtt, Param* par_mat = 0x0) const ;

    /** Computes the flat Laplace operator \f$\Delta\f$ acting on 
     * \f$h^{ij}\f$, using the potentials \e A and \f$\tilde{B}\f$. 
     * The trace of teh result is set to zero (to modify).
     */
    Sym_tensor_tt laplacian_h_tt(Param* par_mat = 0x0) const ;

    // Accessors
    // ---------
    public:
    // Virtual functions from base class Time_slice_conf:
    // -------------------------------------------------
    
    /** Deviation \f$ h^{ij} \f$ 
     * of the conformal metric \f$ \tilde\gamma^{ij} \f$ from 
     * the flat metric \f$ f^{ij} \f$: 
     * \f$\tilde\gamma^{ij} = f^{ij} + h^{ij} \f$.
     * Returns the value at the current time step (\c jtime ).
     */        
    virtual const Sym_tensor& hh(Param* par_bc = 0x0, Param* par_mat = 0x0) const ; 
    
    /** Trace \e K of the extrinsic curvature 
     *  at the current time step (\c jtime ).
     * It is null in the present case (maximal slicing)
     */        
    virtual const Scalar& trk() const ; 
    
    /** Vector \f$ H^i = {\cal D}_j \tilde\gamma^{ij} \f$ 
     * which vanishes in Dirac gauge.
     * It is null in the present case...
     */
    virtual const Vector& hdirac() const ; 
    
    // Virtual functions from this class:
    // ----------------------------------
    
    /** Returns the potential \e A of \f$ \bar{h}^{ij} \f$.
     * See the documentation of \c Sym_tensor for details.
     * Returns the value at the current time step (\c jtime ).
     */
    virtual const Scalar& A_hh() const ; 
    
    /** Returns the potential \f$\tilde{B}\f$ of \f$ \bar{h}^{ij} \f$.
     * See the documentation of \c Sym_tensor_tt for details.
     * Returns the value at the current time step (\c jtime ).
     */
    virtual const Scalar& B_hh() const ;
    
    /** Computes the trace \c h, with respect to the flat metric 
     * \c ff , of \f$ h^{ij} \f$. 
     * Returns the value at the current time step (\c jtime ).
     */
    virtual const Scalar& trh() const ;
    
    
    // Outputs
    // -------
    protected:
    /// Operator >> (virtual function called by the operator<<). 
    virtual ostream& operator>>(ostream& ) const ;	
    
    /** Total or partial saves in a binary file.
     *  This protected method is to be called either from public method 
     *  \c save or from method \c sauve of a derived class. 
     *  
     *  @param fich binary file 
     *  @param partial_save indicates whether the whole object must be
     *      saved.
     */
    virtual void sauve(FILE* fich, bool partial_save) const ; 
  
};
#endif
