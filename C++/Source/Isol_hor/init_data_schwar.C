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
 * Revision 1.4  2004/10/01 16:47:51  f_limousin
 * Case \alpha=0 included
 *
 * Revision 1.3  2004/09/28 16:10:05  f_limousin
 * Many improvements. Now the resolution for the shift is working !
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
				double precis, double relax, int niter,
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

    assert(trk_point.check_dzpuis(4)) ;
    assert(trk_in.check_dzpuis(2)) ; 
    assert(trk_point.check_dzpuis(4)) ; 

    // Initialisations
    // ---------------
    double ttime = the_time[jtime] ;    
    //   trk_evol.update(trk_in, jtime, ttime) ; 


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
    for (int i=0; i<niter; i++) {
	
	//============================================
	//          Boundary conditions
	//============================================
      
	// Boundary for N	
	//---------------
	
	Valeur nn_bound (boundary_nn_Dir(0.2)) ;
	
	// Boundary for Psi (Neumann)
	//---------------------------
	
	Valeur psi_bound (boundary_psi_Neu_spat()) ;

      //=============================================
      // Resolution of elliptic equations
      //=============================================
      
      // Resolution of the Poisson equation for the lapse
      // ------------------------------------------------
     
	
	Scalar source_nn (source_nn_hor(trk_point)) ;
	source_nn.dec_dzpuis(4) ;
//	des_profile(source_nn, 1.00001, 10, 0., 0., "source nn") ;

	Scalar nn_jp1 = source_nn_hor(trk_point)
	    .poisson_dirichlet(nn_bound, 0) + 1. ; 

     // Test:
      maxabs(nn_jp1.laplacian() - source_nn_hor(trk_point),
	     "Absolute error in the resolution of the equation for N") ;  
      
      // Relaxation (relax=1 -> new ; relax=0 -> old )  
      //-----------
	if (i>0)
	    nn_jp1 = relax * nn_jp1 + (1 - relax) * nn() ;
	double diff_nn = max( diffrel(nn(), nn_jp1) ) ;   
	n_evol.update(nn_jp1, jtime, ttime) ; 
          
      // Resolution of the Poisson equation for Psi
      // ------------------------------------------

//     des_profile(psi(), 1.00001, 10, 0., 0., "psi") ;
	Scalar source_psi (source_psi_hor()) ;
	source_psi.dec_dzpuis(4) ;
//	des_profile(source_psi, 1.00001, 10, 0., 0., "source psi") ;

       Scalar psi_jp1 = source_psi_hor().poisson_neumann(psi_bound, 0) + 1. ;


      // Relaxation (relax=1 -> new ; relax=0 -> old )  
      //-----------
      
       psi_jp1 = relax * psi_jp1 + (1 - relax) * psi() ;
     
      
      // Test:
      maxabs(psi_jp1.laplacian() - source_psi_hor(),
	     "Absolute error in the resolution of the equation for Psi") ;  


              
      // Resolution of the vector Poisson equation for the shift
      //---------------------------------------------------------
	  
      Vector beta_jp1(beta()) ;
      
      // Source
      //-------

//      des_profile(beta()(1), 1., 10, 0., 0., "beta") ;


      Vector source_vector ( source_beta_hor() ) ;
      double lambda = 0. ;
      Vector source_reg = - (1./3. - lambda) * beta().divergence(ff)
	  .derive_con(ff) ;
      source_reg.inc_dzpuis() ;
      source_vector = source_vector + source_reg ;

      // Boundary values
      //----------------

      Valeur boundary_x ( boundary_beta_x() ) ;
      Valeur boundary_y ( boundary_beta_y() ) ;
      Valeur boundary_z ( boundary_beta_z() ) ;
   
      // Resolution
      //-----------

      double precision = 1e-8 ;
      poisson_vect_boundary(lambda, source_vector, beta_jp1, boundary_x, 
			    boundary_y, boundary_z, 0, precision, 20) ;

        
      // Check of the resolution
      // ------------------------
      
      source_vector.dec_dzpuis() ;
      maxabs(beta_jp1.derive_con(ff).divergence(ff) 
	     + lambda * beta_jp1.divergence(ff)
	     .derive_con(ff) - source_vector,
	     "Absolute error in the resolution of the equation for beta") ;  
      cout << endl ;

//      des_profile(beta_jp1(1), 1.00001, 10, 0., 0., "beta (1)") ;


      // Relaxation (relax=1 -> new ; relax=0 -> old )  
      //-----------
            
      beta_jp1 = relax * beta_jp1 + (1 - relax) * beta() ;
      //     beta_jp1.set(3) = 0. ; 
      //   des_profile(beta_jp1(1), 1.00001, 10, 0., 0., "beta") ;
     
      
      
      //===========================================
      //      Convergence control
      //===========================================
      
      //     double diff_nn = max( diffrel(nn(), nn_jp1) ) ;   
      double diff_psi = max( diffrel(psi(), psi_jp1) ) ; 
      double diff_beta = max( maxabs(beta_jp1 - beta()) ) ; 
       
      cout << "step = " << i << " :  diff_psi = " << diff_psi 
	   << "  diff_nn = " << diff_nn 
	   << "  diff_beta = " << diff_beta << endl ;
      cout << "----------------------------------------------" << endl ;
      if ( (diff_psi < precis) && (diff_nn < precis) && (diff_beta < precis) )
	break ; 
      
      //=============================================
      //      Updates for next step 
      //=============================================
      
      set_psi_del_q(psi_jp1) ; 
      
      n_evol.update(nn_jp1, jtime, ttime) ; 
      
      beta_evol.update(beta_jp1, jtime, ttime) ; 
      
      // New value of A^{ij}:

      Sym_tensor aa_jp1 (map, CON, map.get_bvect_spher()) ;
      int nnr = map.get_mg()->get_nr(1) ;
      int nnp = map.get_mg()->get_np(1) ;
      int nnt = map.get_mg()->get_nt(1) ;

      int check ;
      check = 0 ;
      for (int k=0; k<nnp; k++)
	  for (int j=0; j<nnt; j++){
	      if (nn().val_grid_point(1, k, j , 0) < 1e-1){
		  check = 1 ;
		  break ;
	      }
	  }

      if (check == 0)
	  aa_jp1 = ( 0*beta().ope_killing_conf(tgam()) + uu ) 
	      / (2.* nn()) ;            
      else {
	  Scalar nn_sxpun (division_xpun (Cmp(nn()), 0)) ;
	  
	  Sym_tensor aa_sxpun = beta().ope_killing_conf(tgam()) + uu ;
	  
	  for(int i=1; i<=3; i++)
	      for(int j=1; j<=i; j++){
		  aa_sxpun.set(i,j).set_inner_boundary(1, 0.) ;
		  aa_sxpun.set(i,j) = Scalar (division_xpun 
					    (Cmp(aa_sxpun(i,j)), 0)) ;
//		  aa_jp1.set(i,j) = aa_sxpun / (2*nn_sxpun) ;
	      }
 		  aa_jp1 = aa_sxpun / (2*nn_sxpun) ;
     }
      
      cout << "aa" << endl << norme(aa_jp1(1,1)/(nnr*nnt*nnp)) << endl ;
      cout << "beta" << endl << norme(beta()(1)/(nnr*nnt*nnp)) << endl ;
      //     des_profile(aa_jp1(1,1), 0.9, 10, 0., 0., "aa(1,1)") ;
      
      
      cout << "check = " << check << endl ;

      //## Alternative formula:
      // Sym_tensor aa_jp1 = ( beta().ope_killing_conf(ff) 
      //                      - hh().derive_lie(beta())
      //                      - 0.6666666666666666 * beta.divergence() * hh()
      //                      + uu ) / (2.* nn()) ; 
      
      set_aa(aa_jp1) ; 
      

    }

} 



