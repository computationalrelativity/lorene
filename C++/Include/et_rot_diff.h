/*
 *  Definition of Lorene class Et_rot_diff
 *
 */

/*
 *   Copyright (c) 2001 Eric Gourgoulhon
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


#ifndef __ET_ROT_DIFF_H_ 
#define __ET_ROT_DIFF_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.2  2001/10/25  09:20:35  eric
 * Ajout de la fonction virtuelle display_poly.
 *
 * Revision 1.1  2001/10/19  08:17:59  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

/**
 * Class for differentially rotating stars.
 * 
 *
 * @version #$Id$#
 */

// Headers Lorene
#include "etoile.h"

class Et_rot_diff : public Etoile_rot {

    // Data : 
    // -----
    protected:
	/** Function $F(\Omega)$ defining the rotation profile.
	 *  This function is linked to the components of the fluid 4-velocity
	 *  by 
	 *  \begin{equation}
	 *	F(\Omega) = u^t u_\varphi \ . 
	 *  \end{equation}
	 *  The first argument of {\tt frot} must be $\Omega$; the second
	 *  argument contains the parameters, the first of which being
	 *  necessarily the central value of $\Omega$.   
	 */
	double (*frot)(double, const Tbl&) ;
	
	/** Primitive of the function $F(\Omega)$, which vanishes at the
	 *  stellar center. 
	 *  The first argument of {\tt primfrot} must be $\Omega$; the second
	 *  argument contains the parameters, the first of which being
	 *  necessarily the central value of $\Omega$.   
	 */
	double (*primfrot)(double, const Tbl&) ;
	
	/** Parameters of the function $F(\Omega)$.
	 * 
	 * To be used as the second argument of functions {\tt frot}
	 * and {\tt primfrot}. 
	 * The parameter {\tt par\_frot(0)} must always be the central angular
	 * velocity. 
	 */
	Tbl par_frot ; 
	
	/// Field $\Omega(r,\theta)$
	Tenseur omega_field ; 
	
	double omega_min ;   /// Minimum value of $\Omega$
	double omega_max ;   /// Maximum value of $\Omega$
	
	/// Field $\int_{\Omega_{\rm c}}^\Omega F(\Omega') \, d\Omega' $
	Tenseur prim_field ; 
	
    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param nzet_i Number of domains occupied by the star
	 * @param relat should be {\tt true} for a relativistic
	 *			star,  {\tt false} for a Newtonian one
	 * @param eos_i Equation of state of the stellar matter
	 * @param frot_i Function $F(\Omega)$ defining the rotation profile.
	 * @param primfrot_i Primitive of $F(\Omega)$ which vanishes at the
	 *   stellar center
	 * @param par_frot_i Parameters of functions {\tt frot\_i}
	 *     and {\tt primfrot\_i}, 
	 *	    {\tt par\_frot\_i(0)} being the central value of
	 *	    $\Omega$. 
	 */
	Et_rot_diff(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i, 
		    double (*frot_i)(double, const Tbl&), 
		    double (*primfrot_i)(double, const Tbl&), 
		    const Tbl& par_frot_i) ;			

	Et_rot_diff(const Et_rot_diff& ) ;		/// Copy constructor

	/** Constructor from a file (see {\tt sauve(FILE* )}). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param fich	input file (must have been created by the function
	 *	{\tt sauve})
	 * @param frot_i Function $F(\Omega)$ defining the rotation profile.
	 * @param primfrot_i Primitive of $F(\Omega)$ which vanishes at the
	 *   stellar center
	 */
	Et_rot_diff(Map& mp_i, const Eos& eos_i, FILE* fich,
		    double (*frot_i)(double, const Tbl&), 
		    double (*primfrot_i)(double, const Tbl&) ) ;			
	
	virtual ~Et_rot_diff() ;			/// Destructor


    // Memory management
    // -----------------

    // Everything is inherited from Etoile_rot

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another {\tt Et\_rot\_diff}
	void operator=(const Et_rot_diff& ) ;	
	
    // Accessors
    // ---------
    public:
	/// Returns the angular velocity field $\Omega$
	const Tenseur& get_omega_field() const {return omega_field;} ; 

	/** Returns the central value of the rotation angular velocity 
	 *  ({\tt [f\_unit]})
	 */ 
	virtual double get_omega_c() const ;	    

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    /// Save in a file
    
	/// Display in polytropic units
	virtual void display_poly(ostream& ) const ; 

    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    


    // Computational routines
    // ----------------------
    public: 

	virtual double tsw() const ;		/// Ratio T/W

	/** Computes the hydrodynamical quantities relative to the Eulerian
	 *  observer from those in the fluid frame.
	 *
	 *  The calculation is performed starting from the quantities
	 *  {\tt omega\_field}, {\tt ent}, {\tt ener}, {\tt press}, 
	 *  and {\tt a\_car},  
	 *  which are supposed to be up to date.  
	 *  From these,  the following fields are updated:
	 *  {\tt gam\_euler}, {\tt u\_euler}, {\tt ener\_euler}, {\tt s\_euler}. 
	 * 
	 */
	virtual void hydro_euler() ; 

	/** Computes $\Omega(r,\theta)$ (member {\tt omega\_field}).
	 *  
	 *  The computation amounts to solving the equation
	 *  \begin{equation}
	 *	F(\Omega) - {B^2 r^2\sin^2\theta (\Omega - N^\varphi)
	 *	    \over N^2 - B^2 r^2 \sin^2(\Omega-N^\varphi)^2} = 0
	 *  \end{equation}
	 *  for $\Omega$. 
	 *
	 *  @param omeg_min [input] Lower bound of the interval for
	 *			     searching omega
	 *  @param omeg_max [input] Higher bound of the interval for
	 *			     searching omega
	 *  @param precis [input] Required precision in the determination of
	 *			the zero by the secant method
	 *  @param nitermax [input] Maximum number of iterations in the secant 
	 *			    method to compute the zero. 
	 */
	void fait_omega_field(double omeg_min, double omeg_max,
			      double precis, int nitermax) ;
	
	/// Computes the member {\tt prim\_field} from {\tt omga\_field}. 
	void fait_prim_field() ;
	
	/** Evaluates $F(\Omega)$, where $F$ is the function
	 *  defining the rotation profile.
	 *  This function is linked to the components of the fluid 4-velocity
	 *  by
	 *  \begin{equation}
	 *	F(\Omega) = u^t u_\varphi \ .
	 *  \end{equation}
	 *
	 *  {\tt funct\_omega} calls {\tt frot} with the parameters
	 *  {\tt par\_frot}.
	 *
	 *  @param omeg [input] value of $\Omega$
	 *  @return value of $F(\Omega)$
	 *
	 */
	double funct_omega(double omeg) const ;  	

	/** Evaluates the primitive of $F(\Omega)$, where $F$ is the function
	 *  defining the rotation profile.
	 *
	 *  @param omeg [input] value of $\Omega$
	 *  @return value of 
	 *	$\int_{\Omega_{\rm c}}^\Omega F(\Omega') \, d\Omega' $
	 *
	 */
	double prim_funct_omega(double omeg) const ;  	

	/** Computes an equilibrium configuration.
	 *  
	 *  @param ent_c  [input] Central enthalpy 
	 *  @param omega0  [input] Requested central angular velocity 
	 *			     (if {\tt fact\_omega=1.})
	 *  @param fact_omega [input] 1.01 = search for the Keplerian frequency,
	 *			      1. = otherwise.
	 *  @param nzadapt  [input] Number of (inner) domains where the mapping 
	 *			    adaptation to an iso-enthalpy surface
	 *			    should be performed
	 *  @param ent_limit [input] 1-D {\tt Tbl} of dimension {\tt nzet} which
	 *				defines the enthalpy at the outer boundary
	 *				of each domain
	 *  @param icontrol [input] Set of integer parameters (stored as a
	 *			    1-D {\tt Itbl} of size 8) to control the 
	 *			    iteration: \\
	 *	{\tt icontrol(0) = mer\_max} : maximum number of steps \\
	 *	{\tt icontrol(1) = mer\_rot} : step at which the rotation is 
	 *				      switched on \\
	 *	{\tt icontrol(2) = mer\_change\_omega} : step at which the rotation
	 *			  velocity is changed to reach the final one  \\
	 *	{\tt icontrol(3) = mer\_fix\_omega} :  step at which the final
	 *			    rotation velocity must have been reached  \\
	 *	{\tt icontrol(4) = mer\_mass} : the absolute value of 
	 *			    {\tt mer\_mass} is the step from which the 
	 *			    baryon mass is forced to converge, 
	 *			    by varying the central enthalpy 
	 *			    ({\tt mer\_mass > 0}) or the angular 
	 *			    velocity ({\tt mer\_mass < 0}) \\
	 *	{\tt icontrol(5) = mermax\_poisson} : maximum number of steps in 
	 *				{\tt Map\_et::poisson} \\
	 *	{\tt icontrol(6) = mer\_triax} : step at which the 3-D 
	 *				perturbation is switched on \\
	 *	{\tt icontrol(7) = delta\_mer\_kep} : number of steps
	 *			    after {\tt mer\_fix\_omega} when {\tt omega}
	 *			    starts to be increased by {\tt fact\_omega}
	 *			    to search for the Keplerian velocity
	 * 	 
	 *  @param control [input] Set of parameters (stored as a 
	 *			    1-D {\tt Tbl} of size 7) to control the 
	 *			    iteration: \\
	 *	{\tt control(0) = precis} : threshold on the enthalpy relative 
	 *				change for ending the computation \\
	 *	{\tt control(1) = omega\_ini} : initial central angular velocity, 
	 *			    switched on only if {\tt mer\_rot < 0}, 
	 *			    otherwise 0 is used \\ 
	 *	{\tt control(2) = relax} : relaxation factor in the main 
	 *				   iteration \\ 
	 *	{\tt control(3) = relax\_poisson} : relaxation factor in 
	 *				   {\tt Map\_et::poisson}\\ 
	 *	{\tt control(4) = thres\_adapt} :  threshold on dH/dr for 
	 *			    freezing the adaptation of the mapping \\
	 *	{\tt control(5) = ampli\_triax} :  relative amplitude of 
	 *			    the 3-D perturbation \\
	 *	{\tt control(6) = precis\_adapt} : precision for 
	 *			    {\tt Map\_et::adapt}
	 *
	 *  @param mbar_wanted [input] Requested baryon mass (effective only 
	 *				if {\tt mer\_mass > mer\_max})
	 *  @param aexp_mass [input] Exponent for the increase factor of the 
	 *			      central enthalpy to converge to the 
	 *			      requested baryon mass
	 *  @param diff [output]   1-D {\tt Tbl} of size 7 for the storage of 
	 *			    some error indicators : \\
	 *	    {\tt diff(0)} : Relative change in the enthalpy field
	 *			      between two successive steps \\
	 *	    {\tt diff(1)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt nuf} \\  
	 *	    {\tt diff(2)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt nuq} \\  
	 *	    {\tt diff(3)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt dzeta} \\  
	 *	    {\tt diff(4)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt tggg} \\  
	 *	    {\tt diff(5)} : Relative error in the resolution of the
	 *			    equation for {\tt shift} (x comp.) \\  
	 *	    {\tt diff(6)} : Relative error in the resolution of the
	 *			    equation for {\tt shift} (y comp.) \\  
	 */
	virtual void equilibrium(double ent_c, double omega0, double fact_omega, 
			 int nzadapt, const Tbl& ent_limit,
			 const Itbl& icontrol, const Tbl& control,
			 double mbar_wanted, double aexp_mass, 
			 Tbl& diff) ;
	
};

#endif
