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


 

/*
 * $Id$
 * $Log$
 * Revision 1.14  2016/12/05 16:17:56  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.13  2014/10/13 08:53:02  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.12  2014/10/06 15:13:12  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.11  2007/05/06 10:48:11  p_grandclement
 * Modification of a few operators for the vorton project
 *
 * Revision 1.10  2007/01/16 15:08:07  n_vasset
 * New constructor, usn Scalar on mono-domain angular grid for boundary,
 * for function sol_elliptic_boundary()
 *
 * Revision 1.9  2005/11/30 11:09:07  p_grandclement
 * Changes for the Bin_ns_bh project
 *
 * Revision 1.8  2005/08/26 14:02:40  p_grandclement
 * Modification of the elliptic solver that matches with an oscillatory exterior solution
 * small correction in Poisson tau also...
 *
 * Revision 1.7  2005/06/09 07:57:31  f_limousin
 * Implement a new function sol_elliptic_boundary() and
 * Vector::poisson_boundary(...) which solve the vectorial poisson
 * equation (method 6) with an inner boundary condition.
 *
 * Revision 1.6  2004/08/24 09:14:42  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.5  2004/06/22 08:49:58  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * Revision 1.4  2004/03/17 15:58:48  p_grandclement
 * Slight modification of sol_elliptic_no_zec
 *
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
#include <cstdlib>
#include <cmath>

// Headers Lorene :
#include "tbl.h"
#include "mtbl_cf.h"
#include "map.h"
#include "param_elliptic.h"
#include "proto.h"         

            //----------------------------------------------
	   //		General elliptic solver
	  //----------------------------------------------

namespace Lorene {
void Map_af::sol_elliptic(Param_elliptic& ope_var, const Scalar& source, 
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
  
  // On met les bonnes bases dans le changement de variable 
  // de ope_var :
  ope_var.var_F.set_spectral_va().coef() ;
  ope_var.var_F.set_spectral_va().ylm() ;
  ope_var.var_G.set_spectral_va().coef() ;
  ope_var.var_G.set_spectral_va().ylm() ;

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


 
            //-----------------------------------------------
	   //		General elliptic solver with boundary as Mtbl-cf
	  //-------------------------------------------------


void Map_af::sol_elliptic_boundary(Param_elliptic& ope_var, const Scalar& source, 
				   Scalar& pot,  const Mtbl_cf& bound, 
				   double fact_dir, double fact_neu ) const {
    
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
  
  // On met les bonnes bases dans le changement de variable 
  // de ope_var :
  ope_var.var_F.set_spectral_va().coef() ;
  ope_var.var_F.set_spectral_va().ylm() ;
  ope_var.var_G.set_spectral_va().coef() ;
  ope_var.var_G.set_spectral_va().ylm() ;

  // Call to the Mtbl_cf version
  // ---------------------------
  Mtbl_cf resu = elliptic_solver_boundary (ope_var, *(rho.c_cf),  bound, 
					   fact_dir, fact_neu) ;
  
  // Final result returned as a Scalar
  // ---------------------------------
  
  pot.set_etat_zero() ;  // to call Scalar::del_t().
  
  pot.set_etat_qcq() ; 
  
  pot.set_spectral_va() = resu ;
  pot.set_spectral_va().ylm_i() ; // On repasse en base standard.	    
  
  pot.set_dzpuis(0) ; 
}





            //-----------------------------------------------
	   //		General elliptic solver with boundary as Scalar
	  //-------------------------------------------------


void Map_af::sol_elliptic_boundary(Param_elliptic& ope_var, const Scalar& source, 
				   Scalar& pot,  const Scalar& bound, 
				   double fact_dir, double fact_neu ) const {
    
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
  
  // On met les bonnes bases dans le changement de variable 
  // de ope_var :
  ope_var.var_F.set_spectral_va().coef() ;
  ope_var.var_F.set_spectral_va().ylm() ;
  ope_var.var_G.set_spectral_va().coef() ;
  ope_var.var_G.set_spectral_va().ylm() ;

  // Call to the Mtbl_cf version
  // ---------------------------
 
// REMINDER : The scalar bound must be defined on a mono-domain angular grid corresponding with the full grid of the scalar source (following assert()) 
 
  Scalar bbound = bound;
  bbound.set_spectral_va().ylm() ;
  const Map& mapp = bbound.get_mp();
 
  const Mg3d& gri2d = *mapp.get_mg();

  assert(&gri2d == source.get_mp().get_mg()->get_angu_1dom()) ;  // Attention à cet assert !! 
  
  Mtbl_cf bound2 (gri2d , bbound.get_spectral_base()) ;
  bound2.annule_hard() ;  

  if (bbound.get_etat() != ETATZERO){ 
 
      int nr = gri2d.get_nr(0) ;
      int nt = gri2d.get_nt(0) ; 
      int np = gri2d.get_np(0) ; 
       
	  for(int k=0; k<np+2; k++)
	      for (int j=0; j<=nt-1; j++)
		  for(int xi=0; xi<= nr-1; xi++)
		  {

  bound2.set(0, k , j , xi) = (*bbound.get_spectral_va().c_cf)(0, k, j, xi) ;   
  }
  }
  Mtbl_cf resu = elliptic_solver_boundary (ope_var, *(rho.c_cf),  bound2, 
					   fact_dir, fact_neu) ;
  
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

void Map_af::sol_elliptic_no_zec(Param_elliptic& ope_var, const Scalar& source, 
			  Scalar& pot, double val) const {
    
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
    
  // On met les bonnes bases dans le changement de variable 
  // de ope_var :
  ope_var.var_F.set_spectral_va().coef() ;
  ope_var.var_F.set_spectral_va().ylm() ;
  ope_var.var_G.set_spectral_va().coef() ;
  ope_var.var_G.set_spectral_va().ylm() ;

  // Call to the Mtbl_cf version
  // ---------------------------
  Mtbl_cf resu = elliptic_solver_no_zec (ope_var, *(rho.c_cf), val) ;
  
  // Final result returned as a Scalar
  // ---------------------------------
  
  pot.set_etat_zero() ;  // to call Scalar::del_t().
  
  pot.set_etat_qcq() ; 
  
  pot.set_spectral_va() = resu ;
  pot.set_spectral_va().ylm_i() ; // On repasse en base standard.	    
  
  pot.set_dzpuis(0) ; 
}

          //----------------------------------------------
	   //	   General elliptic solver only in the ZEC
	  //----------------------------------------------

void Map_af::sol_elliptic_only_zec(Param_elliptic& ope_var, 
				   const Scalar& source, 
				   Scalar& pot, double val) const {
    
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
    
  // On met les bonnes bases dans le changement de variable 
  // de ope_var :
  ope_var.var_F.set_spectral_va().coef() ;
  ope_var.var_F.set_spectral_va().ylm() ;
  ope_var.var_G.set_spectral_va().coef() ;
  ope_var.var_G.set_spectral_va().ylm() ;

  // Call to the Mtbl_cf version
  // ---------------------------
  Mtbl_cf resu = elliptic_solver_only_zec (ope_var, *(rho.c_cf), val) ;
  
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

void Map_af::sol_elliptic_sin_zec(Param_elliptic& ope_var, 
				  const Scalar& source, Scalar& pot, double* amplis, double* phases) const {
    
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
    
  // On met les bonnes bases dans le changement de variable 
  // de ope_var :
  ope_var.var_F.set_spectral_va().coef() ;
  ope_var.var_F.set_spectral_va().ylm() ;
  ope_var.var_G.set_spectral_va().coef() ;
  ope_var.var_G.set_spectral_va().ylm() ;

  // Call to the Mtbl_cf version
  // ---------------------------
  Mtbl_cf resu = elliptic_solver_sin_zec (ope_var, *(rho.c_cf), amplis, phases) ;
  
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
					 Param_elliptic& ope_var, 
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
    
  // On met les bonnes bases dans le changement de variable 
  // de ope_var :
  ope_var.var_F.set_spectral_va().coef() ;
  ope_var.var_F.set_spectral_va().ylm() ;
  ope_var.var_G.set_spectral_va().coef() ;
  ope_var.var_G.set_spectral_va().ylm() ;

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

}
