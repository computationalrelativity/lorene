/*
 *  Definition of Lorene class Isol_hole
 *				 
 */

/*
 *   Copyright (c) 2009 Nicolas Vasset
 */


#ifndef __ISOLHOLE_H_ 
#define __ISOLHOLE_H_ 


// Headers Lorene
#include "tensor.h"
#include "metric.h"
#include "spheroid.h"
#include "excision_surf.h"
#include "excision_hor.h"

 class Eos ;

			//---------------------------//
			//    Class Isol_hole        //
			//---------------------------//

namespace Lorene {
/** Class to compute quasistationary single black hole spacetimes in vacuum.
 *  It takes as arguments:
 * - A Mapping \c mp
 * - A rotation rate \f$ \Omega\f$ of the horizon in the \f$ \varphi \f$
 * direction, assumed to be the direction of the rotational symmetry. 
 * - A Scalar NoK, storing the boundary value to be verified either for the
 * lapse or the surface gravity at the excised boundary (which is a MOTS in 
 * the philosophy of this class). Whether it corresponds to the lapse or the
 * surface gravity is controlled by the boolean NorKappa. (Warning: setting a
 * boundary on Kappa is not implemented yet).
 * - A boolean isCF, to determine whether or not one adopts the conformally 
 * flat approximation in the spacetime construction. Set to FALSE by default.
 *
 *   The main goal of this class is the computation of metric data from those
 *  parameters, as well as some global quantities related to the physical 
 * characteristics of the spacetime. Those metric data are given on a 
 * spacelike 3-slice, assuming a global timelike Killing field and using
 *  adapted coordinates.
 */


class Isol_hole {

    // Data : 
    // -----
    protected:
	const Map& mp ;	    ///< Mapping associated with the star

        double Omega ;     /** Rotation rate of the horizon in the azimuthal
                            *direction.
			    */
	
        bool NorKappa ;    /** Indicates if the boundary value for the lapse
			    * or the surface gravity is used.
			    */

        Scalar boundNoK ; /** Stores the boundary value of the lapse or surface
			   * gravity.
			   */
	
	bool isCF ;       ///< Indicates if the CF approximation is used.

	
	// Metric data
	// -----------------

	/// Lapse function 
	Scalar lapse ; 
	
	// Conformal factor
	Scalar conf_fact;

	/// Shift vector
	Vector shift ;
	
	/** Deviation tensor( non-flat part of the conformal 3-metric on
	 * the slice; see Bonazzola et al. (2004)).
	 */
        Sym_tensor hij ;

	/**  Rescaled tracefree extrinsic curvature tensor: see Cordero et al.(          *    2009).
	  */
	
	Sym_tensor hatA ; 

   // Derived data : 
    // ------------
    protected:

    
	/// Computation of the spheroid associated with the black hole horizon

	mutable Spheroid* p_hor ;
	///   Computation of the ADM mass of the BH spacetime.
	mutable double* p_adm_mass ; 

	/** Computation of the Komar angular momentum w.r.t. assumed
	*rotational symmetry
	  */
	mutable double* p_komar_angmom ;

	/**Computation of the Virial residual, as difference at infinity
	 *between the ADM mass and the Komar integral associated to the mass.
	 */
	mutable double* p_virial_residue ;
	

	


    // Constructors - Destructor
    // -------------------------
    public:
    
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the black hole slice will be defined
	 * @param Omega_i rotation rate of the horizon
	 * @param NorKappa_i FALSE: use of boundary value for the lapse
	 *                   TRUE use of boundary value for the surface gravity
	 *                   (Warning! not implemented yet!)
	 * @param NoK_i value either for the lapse or surface gravity.
	 * @param isCF_i FALSE: Full GR 3+1 equations to be verified
	 *               TRUE: IWM approximation used in determination of the 
	 *               spacetime geometry;
	 */
	Isol_hole(const Map& mp_i, double Omega_i, bool NorKappa_i, Scalar NoK_i, bool isCF_i = false) ;			
	
	
	Isol_hole(const Isol_hole& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE* )). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param Omega_i rotation rate of the horizon
	 * @param NorKappa_i FALSE: use of boundary value for the lapse
	 *                   TRUE use of boundary value for the surface gravity
	 *                   (Warning! not implemented yet!)
	 * @param NoK_i value either for the lapse or surface gravity.
	 * @param isCF_i FALSE: Full GR 3+1 equations to be verified
	 *               TRUE: IWM approximation used in determination of the 
	 *               spacetime geometry;
	 * @param fich	input file (must have been created by the function
	 *	\c sauve)
	 */
	Isol_hole(const Map& mp_i, double Omega_i, bool NorKappa_i, Scalar NoK_i, bool isCF_i, FILE* fich) ;    		

	virtual ~Isol_hole() ;			///< Destructor

	
    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
	virtual void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Isol_hole
	void operator=(const Isol_hole&) ;	

	/** Computes a quasi-stationary 3-slice from the chosen parameters. 
	 *  
	 *  @param precis [input] threshold in the relative difference between 
	 *	the radial shift for two consecutive steps
	 *	to stop the iterative procedure (default value: 1.e-11)
	 */
	void compute_stat_metric(double precis, double relax, int mer_max, int mer_max2, bool isvoid = true) ; 

	/** Computes the rhs of hyperbolic equation for conformal metric
	 *assuming statioarity; 
	 *WARNING; up to now, we are only able to handle void spacetimes.
	 */

	void secmembre_kerr(Sym_tensor& source_hh);
 
 

    // Accessors
    // ---------
    public:
	/// Returns the mapping
	const Map& get_mp() const {return mp; } ; 

	/// Returns the rotation rate
	double get_Omega() const {return Omega ;} ;

	/// Returns the boundary value used for the lapse (if it is the one used)
	const Scalar& get_boundN() const{
	if (NorKappa == false){
	  return boundNoK;
	}
	else cout << "The boundary condition is imposed on the surface gravity!" << endl;
	}
	/// Returns the surface gravity value at the boundary (if it is the one used)
	
	const Scalar& get_Kappa() const{
	  if(NorKappa == true){
	  return boundNoK;
	}
	else cout << "The boundary condition is imposed on the lapse!" <<endl;
	}


	/// Returns the lapse function \e N
	const Scalar& get_lapse() const {return lapse;} ;
	
	/// Returns the conformal factor 
	const Scalar& get_conf_fact() const {return conf_fact;};

	/// Returns the shift vector \f$\beta^i\f$.
	const Vector& get_shift() const {return shift;} ;

	/// Returns the deviation tensor \f$ h^{ij} \f$.
	const Sym_tensor& get_hij() const {return hij;} ;

	/// Returns the rescaled tracefree extrinsic curvature \f$\hat{A}^{ij}\f$.
	const Sym_tensor& get_hatA() const {return hatA;} ;




    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	    ///< Save in a file

        void Einstein_errors();                     ///< Prints out errors in Einstein equations for the data obtained.

 
    // Global quantities
    // -----------------
    public:

      
	/** Spheroid at the excised boundary associated with the black hole
	* MOTS on the slice. Set by default at the position \f$ r=1 \f$.
	*/

	Spheroid hor() ;
 
	///   Computation of the ADM mass of the BH spacetime.
	double adm_mass() ; 

	/** Computation of the Komar angular momentum w.r.t. assumed
	*rotational symmetry
	  */
	double komar_angmom() ;

	/**Computation of the Virial residual, as difference at infinity
	 *between the ADM mass and the Komar integral associated to the mass.
	 */
	double virial_residue() ;

	

};

}
#endif
