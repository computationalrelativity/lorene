/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char map_af_elliptic_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/02/11 09:47:46  p_grandclement
 * Addition of a new elliptic solver, matching with the homogeneous solution
 * at the outer shell and not solving in the external domain (more details
 * coming soon ; check your local Lorene dealer...)
 *
 * Revision 1.2  2004/01/28 16:46:23  p_grandclement
 * Addition of the sol_elliptic_fixe_der_zero stuff
 *
 * Revision 1.1  2003/12/11 14:48:48  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * 
 * $Header$
 *
 */

// Header C : 
#include <stdlib.h>
#include <math.h>

// Headers Lorene :
#include "tbl.h"
#include "mtbl_cf.h"
#include "map.h"
#include "param_elliptic.h"
         
            //----------------------------------------------
	   //		General elliptic solver
	  //----------------------------------------------

void Map_af::sol_elliptic(const Param_elliptic& ope_var, const Scalar& source, 
			  Scalar& pot) const {
    
  assert(source.get_etat() != ETATNONDEF) ; 
  assert(source.get_mp().get_mg() == mg) ; 
  assert(pot.get_mp().get_mg() == mg) ; 
  
  assert(source.check_dzpuis(2) || source.check_dzpuis(3) || 
	 source.check_dzpuis(4)) ; 
  // Spherical harmonic expansion of the source
  // ------------------------------------------
  
  const Valeur& sourva = source.get_spectral_va() ; 
  
  if (sourva.get_etat() == ETATZERO) {
    pot.set_etat_zero() ;
    return ;  
    }
  
  // Spectral coefficients of the source
  assert(sourva.get_etat() == ETATQCQ) ; 
  
  Valeur rho(sourva.get_mg()) ; 
  sourva.coef() ; 
  rho = *(sourva.c_cf) ;	// copy of the coefficients of the source
  
  rho.ylm() ;			// spherical harmonic transforms 
  
  // Call to the Mtbl_cf version
  // ---------------------------
  Mtbl_cf resu = elliptic_solver (ope_var, *(rho.c_cf)) ;
  
  // Final result returned as a Scalar
  // ---------------------------------
  
  pot.set_etat_zero() ;  // to call Scalar::del_t().
  
  pot.set_etat_qcq() ; 
  
  pot.set_spectral_va() = resu ;
  pot.set_spectral_va().ylm_i() ; // On repasse en base standard.	    
  
  pot.set_dzpuis(0) ; 
}

            //----------------------------------------------
	   //	   General elliptic solver with no ZEC
	  //----------------------------------------------

void Map_af::sol_elliptic_no_zec(const Param_elliptic& ope_var, const Scalar& source, 
			  Scalar& pot) const {
    
  assert(source.get_etat() != ETATNONDEF) ; 
  assert(source.get_mp().get_mg() == mg) ; 
  assert(pot.get_mp().get_mg() == mg) ; 
  
  assert(source.check_dzpuis(2) || source.check_dzpuis(3) || 
	 source.check_dzpuis(4)) ; 
  // Spherical harmonic expansion of the source
  // ------------------------------------------
  
  const Valeur& sourva = source.get_spectral_va() ; 
  
  if (sourva.get_etat() == ETATZERO) {
    pot.set_etat_zero() ;
    return ;  
    }
  
  // Spectral coefficients of the source
  assert(sourva.get_etat() == ETATQCQ) ; 
  
  Valeur rho(sourva.get_mg()) ; 
  sourva.coef() ; 
  rho = *(sourva.c_cf) ;	// copy of the coefficients of the source
  
  rho.ylm() ;			// spherical harmonic transforms 
  
  // Call to the Mtbl_cf version
  // ---------------------------
  Mtbl_cf resu = elliptic_solver_no_zec (ope_var, *(rho.c_cf)) ;
  
  // Final result returned as a Scalar
  // ---------------------------------
  
  pot.set_etat_zero() ;  // to call Scalar::del_t().
  
  pot.set_etat_qcq() ; 
  
  pot.set_spectral_va() = resu ;
  pot.set_spectral_va().ylm_i() ; // On repasse en base standard.	    
  
  pot.set_dzpuis(0) ; 
}


            //----------------------------------------------
	   //	   General elliptic solver with no ZEC
           //  and a mtaching with sin(r)/r
	  //----------------------------------------------

void Map_af::sol_elliptic_sin_zec(const Param_elliptic& ope_var, 
				  const Scalar& source, Scalar& pot, 
				  double freq, 
				  int nbr_phase, double& ampli_min, 
				  double& phase_min) const {
    
  assert(source.get_etat() != ETATNONDEF) ; 
  assert(source.get_mp().get_mg() == mg) ; 
  assert(pot.get_mp().get_mg() == mg) ; 
  
  assert(source.check_dzpuis(2) || source.check_dzpuis(3) || 
	 source.check_dzpuis(4)) ; 
  // Spherical harmonic expansion of the source
  // ------------------------------------------
  
  const Valeur& sourva = source.get_spectral_va() ; 
  
  if (sourva.get_etat() == ETATZERO) {
    pot.set_etat_zero() ;
    return ;  
    }
  
  // Spectral coefficients of the source
  assert(sourva.get_etat() == ETATQCQ) ; 
  
  Valeur rho(sourva.get_mg()) ; 
  sourva.coef() ; 
  rho = *(sourva.c_cf) ;	// copy of the coefficients of the source
  
  rho.ylm() ;			// spherical harmonic transforms 
  
  // Call to the Mtbl_cf version
  // ---------------------------
  Mtbl_cf resu = elliptic_solver_sin_zec (ope_var, *(rho.c_cf), freq, 
					  nbr_phase, ampli_min, phase_min) ;
  
  // Final result returned as a Scalar
  // ---------------------------------
  
  pot.set_etat_zero() ;  // to call Scalar::del_t().
  
  pot.set_etat_qcq() ; 
  
  pot.set_spectral_va() = resu ;
  pot.set_spectral_va().ylm_i() ; // On repasse en base standard.	    
  
  pot.set_dzpuis(0) ; 
}


            //----------------------------------------------
	   //	   General elliptic solver with no ZEC
	  //----------------------------------------------

void Map_af::sol_elliptic_fixe_der_zero (double valeur, 
					 const Param_elliptic& ope_var, 
					 const Scalar& source, 
					 Scalar& pot) const {
    
  assert(source.get_etat() != ETATNONDEF) ; 
  assert(source.get_mp().get_mg() == mg) ; 
  assert(pot.get_mp().get_mg() == mg) ; 
  
  assert(source.check_dzpuis(2) || source.check_dzpuis(3) || 
	 source.check_dzpuis(4)) ; 
  // Spherical harmonic expansion of the source
  // ------------------------------------------
  
  const Valeur& sourva = source.get_spectral_va() ; 
  
  if (sourva.get_etat() == ETATZERO) {
    pot.set_etat_zero() ;
    return ;  
    }
  
  // Spectral coefficients of the source
  assert(sourva.get_etat() == ETATQCQ) ; 
  
  Valeur rho(sourva.get_mg()) ; 
  sourva.coef() ; 
  rho = *(sourva.c_cf) ;	// copy of the coefficients of the source
  
  rho.ylm() ;			// spherical harmonic transforms 
  
  // Call to the Mtbl_cf version
  // ---------------------------
  valeur *= alpha[0] ;
  Mtbl_cf resu = elliptic_solver_fixe_der_zero (valeur, ope_var, *(rho.c_cf)) ;
  
  // Final result returned as a Scalar
  // ---------------------------------
  
  pot.set_etat_zero() ;  // to call Scalar::del_t().
  
  pot.set_etat_qcq() ; 
  
  pot.set_spectral_va() = resu ;
  pot.set_spectral_va().ylm_i() ; // On repasse en base standard.	    
  
  pot.set_dzpuis(0) ; 
}

