/*
 *  Definition of Lorene class Gravastar
 *				 
 */

/*
 *   Copyright (c) 2010 Frederic Vincent
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


#ifndef __GRAVASTAR_H_ 
#define __GRAVASTAR_H_ 



// Headers Lorene
#include "star_rot.h"

namespace Lorene {
class Eos ;

			//--------------------------//
			//    class Gravastar       //
			//--------------------------//

/**
 * Class for perfect fluid rotating gravastar. \ingroup (star)
 * 
 
 * 
 */
class Gravastar : public Star_rot {

    // Data : 
    // -----
 protected:
  /**
   *  Energy density in gravastar's core
   */
  double rho_core;

    // Derived data : 
    // ------------
    protected:

	 

    // Constructors - Destructor
    // -------------------------
    public:
  /** Standard constructor. 
   * 
   * @param mp_i Mapping on which the gravastar is contructed
   * @param nzet_i Number of domains occupied by the gravastar
   * @param eos_i Equation of state of the crust matter
   * @param rho_core_i Energy density in gravastar's core (constant, =-p_core)
   *                   NB:pas de choix pour l'EOS en fait : a virer??
   */
  Gravastar(Map& mp_i, int nzet_i, const Eos& eos_i, const double rho_core_i) ;
	
	
  Gravastar(const Gravastar& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE*) ). 
	 * 
	 * @param mp_i Mapping on which the gravastar is constructed
	 * @param eos_i Equation of state of the crust matter
	 * @param fich	input file (must have been created by the function
	 *	\c Gravastar::sauve )
	 */
	Gravastar(Map& mp_i, const Eos& eos_i, FILE* fich) ;    		

	~Gravastar() ;			///< Destructor


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Gravastar 
	void operator=(const Gravastar& ) ;

	/** Allows to computes the proper baryon and energy density, as well as
	 *  pressure from the enthalpy, in the gravastar's crust only (in the core,
	 *  p=-rho=cst)
	 */
	void equation_of_state() ; 	

	// Computational routines
	// ----------------------
 public: 

	/** Computes an equilibrium configuration.
	 *  
	 *  @param omega0  [input] Requested angular velocity 
	 *			     (if \c fact_omega=1. )
	 *  @param fact_omega [input] 1.01 = search for the Keplerian frequency,
	 *			      1. = otherwise.
	 *  @param nzadapt  [input] Number of (inner) domains where the mapping 
	 *			    adaptation to an iso-enthalpy surface
	 *			    should be performed
	 *  @param ent_limit [input] 1-D \c Tbl  of dimension \c nzet  which
	 *				defines the enthalpy at the outer boundary
	 *				of each domain
	 *  @param icontrol [input] Set of integer parameters (stored as a
	 *			    1-D \c Itbl  of size 8) to control the 
	 *			    iteration: 
	 *	\li \c icontrol(0) = mer_max  : maximum number of steps 
	 *	\li \c icontrol(1) = mer_rot  : step at which the rotation is 
	 *				      switched on 
	 *	\li \c icontrol(2) = mer_change_omega  : step at which the rotation
	 *			  velocity is changed to reach the final one  
	 *	\li \c icontrol(3) = mer_fix_omega  :  step at which the final
	 *			    rotation velocity must have been reached  
	 *	\li \c icontrol(4) = mer_mass  : the absolute value of 
	 *			    \c mer_mass  is the step from which the 
	 *			    baryon mass is forced to converge, 
	 *			    by varying the central enthalpy 
	 *			    (\c mer_mass>0 ) or the angular 
	 *			    velocity (\c mer_mass<0 ) 
	 *	\li \c icontrol(5) = mermax_poisson  : maximum number of steps in 
	 *				\c Map_et::poisson  
	 *	\li \c icontrol(6) = mer_triax  : step at which the 3-D 
	 *				perturbation is switched on 
	 *	\li \c icontrol(7) = delta_mer_kep  : number of steps
	 *			    after \c mer_fix_omega  when \c omega 
	 *			    starts to be increased by \c fact_omega 
	 *			    to search for the Keplerian velocity
	 * 	 
	 *  @param control [input] Set of parameters (stored as a 
	 *			    1-D \c Tbl  of size 7) to control the 
	 *			    iteration: 
	 *	\li \c control(0) = precis  : threshold on the enthalpy relative 
	 *				change for ending the computation 
	 *	\li \c control(1) = omega_ini  : initial angular velocity, 
	 *			    switched on only if \c mer_rot<0 , 
	 *			    otherwise 0 is used  
	 *	\li \c control(2) = relax  : relaxation factor in the main 
	 *				   iteration  
	 *	\li \c control(3) = relax_poisson  : relaxation factor in 
	 *				   \c Map_et::poisson  
	 *	\li \c control(4) = thres_adapt  :  threshold on dH/dr for 
	 *			    freezing the adaptation of the mapping 
	 *	\li \c control(5) = ampli_triax  :  relative amplitude of 
	 *			    the 3-D perturbation 
	 *	\li \c control(6) = precis_adapt  : precision for 
	 *			    \c Map_et::adapt 
	 *
	 *  @param diff [output]   1-D \c Tbl  of size 7 for the storage of 
	 *			    some error indicators : 
	 *	    \li \c diff(0)  : Relative change in the enthalpy field
	 *			      between two successive steps 
	 *	    \li \c diff(1)  : Relative error in the resolution of the
	 *			    Poisson equation for \c nuf    
	 *	    \li \c diff(2)  : Relative error in the resolution of the
	 *			    Poisson equation for \c nuq    
	 *	    \li \c diff(3)  : Relative error in the resolution of the
	 *			    Poisson equation for \c dzeta    
	 *	    \li \c diff(4)  : Relative error in the resolution of the
	 *			    Poisson equation for \c tggg    
	 *	    \li \c diff(5)  : Relative error in the resolution of the
	 *			    equation for \c shift  (x comp.)   
	 *	    \li \c diff(6)  : Relative error in the resolution of the
	 *			    equation for \c shift  (y comp.)   
	 */
	void equilibrium(double omega0, double fact_omega, 
			 int nzadapt, const Tbl& ent_limit,
			 const Itbl& icontrol, const Tbl& control,
			 Tbl& diff, Param* = 0x0) ;

    // Outputs
    // -------
    public:


    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    
	
};

}
#endif
