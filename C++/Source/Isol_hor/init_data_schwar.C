/*
 *  Method of class Hor_isol to compute valid initial data for standard boundary 
 *   conditions
 *
 *    (see file isol_hor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Jose Luis Jaramillo
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

char init_data_schwar_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/09/17 13:34:50  f_limousin
 * Introduction of relaxation
 *
 * Revision 1.1  2004/09/09 16:41:50  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>
#include <assert.h>

// Lorene headers
#include "time_slice.h"
#include "isol_hor.h"
#include "metric.h"
#include "evolution.h"
#include "unites.h"
#include "graphique.h"
#include "cmp.h"
#include "tenseur.h"
#include "utilitaires.h"

void Isol_hor::init_data_schwar(const Sym_tensor& uu, 
                const Scalar& trk_in, const Scalar& trk_point, 
		double precis, double relax,
                const Scalar* p_ener_dens, const Vector* p_mom_dens, 
                const Scalar* p_trace_stress) {

    using namespace Unites ;
   
    // Verifications
    // -------------
    double tr_uu = max(maxabs(uu.trace(tgam()), "trace tgam_{ij} u^{ij}")) ; 
    if (tr_uu > 1.e-8) {
        cerr << 
        "Time_slice_conf::initial_data_cts : the trace of u^{ij} with respect\n"
        << "  to the conformal metric tgam_{ij} is not zero !\n" 
        << "  error = " << tr_uu << endl ; 
        abort() ; 
    }

    assert(trk_in.check_dzpuis(2)) ; 
    assert(trk_point.check_dzpuis(4)) ; 

    // Initialisations
    // ---------------
    double ttime = the_time[jtime] ; 
         
    trk_evol.update(trk_in, jtime, ttime) ; 


    const Map& map = uu.get_mp() ; 
    const Base_vect& triad = *(uu.get_triad()) ;
    
    Scalar tmp(map) ;
    Scalar tmp_scal(map) ; 
    Sym_tensor tmp_sym(map, CON, triad) ;

    // Reset of quantities depending on K:
    k_dd_evol.downdate(jtime) ; 
    k_uu_evol.downdate(jtime) ; 

    set_aa( ( beta().ope_killing_conf(tgam()) +  uu) / (2.* nn()) ) ; 


    // Iteration
    // ---------
    int imax = 100 ; 
    for (int i=0; i<imax; i++) {
	
	//============================================
	//          Boundary conditions
	//============================================
      
	// Boundary for N	
	//---------------
	
	Valeur nn_bound (boundary_nn_Dir_eff(0)) ;
	
	// Boundary for Psi (Neumann)
	//---------------------------
	
	Valeur psi_bound (boundary_psi_Dir_spat()) ;

      //=============================================
      // Resolution of elliptic equations
      //=============================================
      
      // Resolution of the Poisson equation for the lapse
      // ------------------------------------------------
      
	
	des_profile(source_nn_hor(trk_in), 1.00001, 10, 0., 0., "source") ;
	des_coef_xi(source_nn_hor(trk_in).get_spectral_va(), 1, 0, 0) ;
	des_coef_theta(source_nn_hor(trk_in).get_spectral_va(), 1, 0, 0) ;

	Scalar nn_jp1 = source_nn_hor(trk_in).poisson_dirichlet(nn_bound, 0) 
	    + 1. ; 

	
      // Relaxation (relax=1 -> new ; relax=0 -> old )  
      //-----------
      
	double relax = 0.5 ;      
	nn_jp1 = relax * nn_jp1 + (1 - relax) * nn() ;
      

      // Test:
      maxabs(nn_jp1.laplacian() - source_nn_hor(trk_in),
	     "Absolute error in the resolution of the equation for N") ;  
      
      //        des_meridian(nn_jp1, 0., 5., "N", 2) ; 
      
      // Resolution of the Poisson equation for Psi
      // ------------------------------------------
      
      Scalar psi_jp1 = source_psi_hor().poisson_dirichlet(psi_bound, 0) + 1. ; 

      
      // Relaxation (relax=1 -> new ; relax=0 -> old )  
      //-----------
      
      relax = 0.5 ;      
      psi_jp1 = relax * psi_jp1 + (1 - relax) * psi() ;
     
      
      // Test:
      maxabs(psi_jp1.laplacian() - source_psi_hor(),
	     "Absolute error in the resolution of the equation for Psi") ;  
      
      //        des_meridian(psi_jp1, 0., 5., "Psi", 1) ; 

/*        
      // Resolution of the vector Poisson equation for the shift
      //---------------------------------------------------------
      
      Vector beta_jp1(ff.get_mp(), CON, *(ff.get_triad()) ) ;

      Vector beta_cartesian = beta_bound_cart() ;
   
      Vector beta_cart_old (beta_cartesian) ;
      
      Tenseur beta_t (map, 1, CON, map.get_bvect_cart() ) ;

      beta_t.set_etat_qcq() ;

      Cmp  beta_0 ( beta_cartesian(1) ) ;
      Cmp  beta_1 ( beta_cartesian(2) ) ;
      Cmp  beta_2 ( beta_cartesian(3) ) ;

      beta_t.set(0) =   beta_0 ;
      beta_t.set(1) =   beta_1 ;
      beta_t.set(2) =   beta_2 ;



  
      // Source
      //-------
      Vector source_vector ( source_beta_hor() ) ;

      source_vector.change_triad(map.get_bvect_cart() ) ;
      
      Tenseur source (map, 1, CON, map.get_bvect_cart() ) ;

      source.set_etat_qcq() ;

      Cmp  source_0 (source_vector(1) ) ;
      Cmp  source_1 (source_vector(2) ) ;
      Cmp  source_2 (source_vector(3) ) ;
      
      source.set(0) =   source_0 ;
      source.set(1) =   source_1 ;
      source.set(2) =   source_2 ;



      // Boundary values
      //----------------

      Valeur boundary_x  ( boundary_beta_x() ) ;
      
      Valeur boundary_y  ( boundary_beta_y() ) ;
      
      Valeur boundary_z  ( boundary_beta_z() ) ;
   

      // Resolution
      //-----------
      
      double precision = precis ;
      poisson_vect_frontiere(1./3., source, beta_t, boundary_x, boundary_y, 
			     boundary_z, 0, precision, 20) ;
      
      
      beta_cartesian.set(1) = beta_t(0) ; 
      beta_cartesian.set(2) = beta_t(1) ;
      beta_cartesian.set(3) = beta_t(2) ;

      // Relaxation (relax=1 -> new ; relax=0 -> old )  
      //-----------
      
      double relax = 0. ;         //si se pone 1 no converge! posible problema con la 
                                  //division por cero?
      
      beta_cartesian = relax * beta_cartesian + (1 - relax) * beta_cart_old ;
      
      
      // Change to spherical basis
      //--------------------------
      
      beta_cartesian.change_triad(map.get_bvect_spher() ) ;
      
      
      beta_jp1 = beta_cartesian  ;

*/
      





      //===========================================
      //      Convergence control
      //===========================================
      
      double diff_psi = max( diffrel(psi(), psi_jp1) ) ; 
      double diff_nn = max( diffrel(nn(), nn_jp1) ) ; 
        
      /*
      Vector beta_diff = beta() - beta_jp1 ;
      Scalar mod_diff_beta = contract( beta_diff.down(0, ff), 0,  beta_diff, 0 ) ;
      tmp.set_etat_zero() ;
      double diff_beta = max( diffrel(mod_diff_beta, tmp) ) ; 
      */

      cout << "step = " << i << " :  diff_psi = " << diff_psi 
	   << "  diff_nn = " << diff_nn << endl << endl ;
//	   << "  diff_beta = " << diff_beta << endl ; 
      cout << "----------------------------------------------" << endl ;
	  if ( (diff_psi < precis) && (diff_nn < precis) )//&& (diff_beta < precis) )
	break ; 
      
      //=============================================
      //      Updates for next step 
      //=============================================
      
      set_psi_del_q(psi_jp1) ; 
      
      n_evol.update(nn_jp1, jtime, ttime) ; 
      
//      beta_evol.update(beta_jp1, jtime, ttime) ; 
      
      // New value of A^{ij}:
      Sym_tensor aa_jp1 = ( beta().ope_killing_conf(tgam()) + uu ) 
	/ (2.* nn()) ; 
      
      //## Alternative formula:
      // Sym_tensor aa_jp1 = ( beta().ope_killing_conf(ff) 
      //                      - hh().derive_lie(beta())
      //                      - 0.6666666666666666 * beta.divergence() * hh()
      //                      + uu ) / (2.* nn()) ; 
      
      set_aa(aa_jp1) ; 
      
    }

} 



