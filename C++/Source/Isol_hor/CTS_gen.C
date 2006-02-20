/*
 *  Method of class Isol_hor to compute intital data using  generalized
 *  Conformal Thni Sandwich equations with N = Ntilde \Psi ^a
 *  and K_ij = \Psi ^\zeta A_ij + 1/3 K \gamma_ij
 *
 *    (see file isol_hor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Jose Luis Jaramillo
 *                       Francois Limousin
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

char CTS_gen[] = "$Header$" ;

/*
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
#include "utilitaires.h"
#include "param.h"
#include "vector.h"

void Isol_hor::init_data_CTS_gen(int bound_nn, double lim_nn, int bound_psi, int bound_beta,
				int solve_lapse, int solve_psi, int solve_shift, 
				double precis, double relax_nn ,
				double relax_psi, double relax_beta , 
				int niter, double a, double zeta ) {

     using namespace Unites ;

     // Initialisations
     // ===============

     double lambda = 0. ;              // to be used in the beta source

     double ttime = the_time[jtime] ;    
     const Base_vect& triad = *(ff.get_triad()) ;

     // Output file
     ofstream conv("resconv.d") ; 
     ofstream kss("kss.d") ;
     conv << " # diff_nn   diff_psi   diff_beta " << endl ;

     // Dynamical relaxation
     double relax_nn_fin = relax_nn ;
     double relax_psi_fin = relax_psi ;
     double relax_beta_fin = relax_beta ;

     // Initial values of the unknowns
     Scalar nntilde_j = nn() * pow(psi(), a) ;
     nntilde_j.std_spectral_base() ;
     Scalar psi_j = psi() ;
     Vector beta_j = beta() ;


     // Iteration
     //==========

     for (int mer=0; mer<niter; mer++) {
	

       // Useful functions
       //-----------------
       Scalar temp_scal_0 (mp) ;
       Scalar temp_scal_2 (mp) ;
       Scalar temp_scal_3 (mp) ;
       Scalar temp_scal_4 (mp) ;
       
       Vector temp_vect_2 (mp, CON, triad) ; 
       Vector temp_vect_3 (mp, CON, triad) ; 
       Vector temp_vect_4 (mp, CON, triad) ; 

       Scalar psi_j_a = pow(psi_j, a) ;
       psi_j_a.std_spectral_base() ;
       Scalar psi_j_ma = pow(psi_j, -a) ;
       psi_j_ma.std_spectral_base() ;
       Vector beta_j_d = beta_j.down(0, met_gamt) ;
       Scalar nnpsia_j = nntilde_j *  psi_j_a /( psi_j*psi_j*psi_j*psi_j*psi_j*psi_j) ;    // to be used in the beta source


       //Dynamical relaxation
       //--------------------
       //       double relax_init = 0.05 ;
       //       double relax_speed = 0.005 ;

       //      relax_nn = relax_nn_fin - ( relax_nn_fin - relax_init ) * exp (- relax_speed *  mer) ;
       //      relax_psi = relax_psi_fin - ( relax_psi_fin - relax_init ) * exp (- relax_speed *  mer) ;
       //      relax_beta = relax_beta_fin - ( relax_beta_fin - relax_init ) * exp (- relax_speed *  mer) ;       
       cout << " relax_nn = " << relax_nn << endl ;
       cout << " relax_psi = " << relax_psi << endl ;
       cout << " relax_beta = " << relax_beta << endl ;
  
       //========
       // Sources
       //========
       
       // Source for psi
       //---------------

       Scalar sour_psi (mp) ;

       // Source physique
       
       temp_scal_3 = 1./8. * met_gamt.ricci_scal() * psi_j ;
       temp_scal_3.inc_dzpuis() ;
       temp_scal_4 =  1./12. * trK*trK * psi_j*psi_j*psi_j*psi_j*psi_j
	             - 1./32. * psi_j*psi_j*psi_j*psi_j*psi_j * psi_j_ma*psi_j_ma / (nntilde_j*nntilde_j) *
	               contract(beta_j_d.ope_killing_conf(met_gamt), 0, 1, beta_j.ope_killing_conf(met_gamt), 0, 1) ;
       
       sour_psi = temp_scal_3 + temp_scal_4 ;


       // Correction Source 

       temp_scal_3 = - contract (hh(), 0, 1,  psi_j.derive_cov(ff).derive_cov(ff), 0, 1) ;
       temp_scal_3.inc_dzpuis() ;
       temp_scal_4 = - contract (hdirac(), 0, psi_j.derive_cov(ff), 0) ;
       
       sour_psi +=  temp_scal_3 + temp_scal_4 ;

       sour_psi.annule_domain(0) ;



       /*
       cout << " dzpuis temp_scal_2 =" << temp_scal_2.get_dzpuis() << endl;
       cout << " dzpuis temp_scal_3 =" << temp_scal_3.get_dzpuis() << endl;
       cout << " dzpuis temp_scal_4 =" << temp_scal_4.get_dzpuis() << endl;
       */
       
      


       // Source for nntilde
       //--------------------

       Scalar sour_nntilde (mp) ;

       // Source physique       
  
       temp_scal_2 = -  psi_j*psi_j*psi_j*psi_j * psi_j_ma * trK_point ;
       temp_scal_2.inc_dzpuis(2) ;

       temp_scal_3 =  psi_j*psi_j*psi_j*psi_j * psi_j_ma * contract (beta_j, 0, trK.derive_cov(ff), 0) 
	              - a/8. * nntilde_j * met_gamt.ricci_scal() ;
       temp_scal_3.inc_dzpuis() ;
       
       temp_scal_4 =  nntilde_j * (4-a)/12. * psi_j*psi_j*psi_j*psi_j * trK*trK 
	              - 2 * (a+1) * contract(psi_j.derive_cov(ff), 0, nntilde_j.derive_con(met_gamt), 0) / psi_j 
	              - a*(a+1) * nntilde_j * contract(psi_j.derive_cov(met_gamt), 0, psi_j.derive_con(met_gamt), 0) / (psi_j * psi_j) 
	              + (a+8)/32. * psi_j*psi_j*psi_j*psi_j * psi_j_ma*psi_j_ma/nntilde_j 
	                          * contract(beta_j_d.ope_killing_conf(met_gamt), 0, 1, beta_j.ope_killing_conf(met_gamt), 0, 1) ;

       sour_nntilde = temp_scal_2 + temp_scal_3 + temp_scal_4 ;


       // Correction Source 
       temp_scal_3 = - contract (hh(), 0, 1,  nntilde_j.derive_cov(ff).derive_cov(ff), 0, 1) ;
       temp_scal_3.inc_dzpuis() ;

       temp_scal_4 = - contract (hdirac(), 0, nntilde_j.derive_cov(ff), 0) ;
       
       sour_nntilde +=  temp_scal_3 + temp_scal_4 ;
       sour_nntilde.annule_domain(0) ;


       // Source for beta
       //----------------

       Vector sour_beta(mp, CON, triad) ; 

       // Source physique

       temp_vect_3 =  4./3. * psi_j_a * nntilde_j * trK.derive_con(met_gamt) 
	              - contract(met_gamt.ricci(), 1, beta_j, 0).up(0, met_gamt) ;
       temp_vect_3.inc_dzpuis() ;

       temp_vect_4 = contract(beta_j.ope_killing_conf(met_gamt), 1, nnpsia_j.derive_cov(ff), 0) / nnpsia_j ;
                     
       sour_beta = temp_vect_3 + temp_vect_4 ;

       // Correction Source 


       temp_vect_3 = (lambda - 1./3.) * contract (beta_j.derive_cov(ff).derive_con(ff), 0, 1) 
	             - contract(contract(met_gamt.connect().get_delta(), 1, beta_j, 0).derive_con(ff), 1, 2)
                     - contract(hh(), 0, 1, beta_j.derive_cov(met_gamt).derive_cov(ff), 1, 2) 
                     - 1./3. * contract (hh(), 1, beta_j.divergence(ff).derive_cov(ff), 0) ;
       temp_vect_3.inc_dzpuis() ;

       temp_vect_4 = - contract(hdirac(), 0, beta_j.derive_cov(met_gamt), 1) 
	             - contract(met_gamt.connect().get_delta(), 1, 2, beta_j.derive_con(met_gamt), 0, 1) ;
 
       sour_beta += temp_vect_3 + temp_vect_4 ;


       cout << " dzpuis sour_beta =" << sour_beta(1).get_dzpuis() << endl;

       //====================
       // Boundary conditions
       //====================


       // BC Psi
       //-------
       Scalar tmp = psi_j * psi_j * psi_j * trK
	 - psi_j * met_gamt.radial_vect().divergence(met_gamt) 
         - psi_j_ma/(2. * nntilde_j * psi_j) * contract( beta_j_d.ope_killing_conf(met_gamt), 0, 1,  met_gamt.radial_vect()* met_gamt.radial_vect(), 0, 1)
	 - 1./3. * trK/psi_j * contract( met_gamt.cov(), 0, 1,  met_gamt.radial_vect()* met_gamt.radial_vect(), 0, 1) 
	 - 4 * ( met_gamt.radial_vect()(2) * psi_j.derive_cov(ff)(2)  + met_gamt.radial_vect()(3) * psi_j.derive_cov(ff)(3) ) ;
       
       tmp = tmp / (4 *  met_gamt.radial_vect()(1)) ;

       // in this case you don't have to substract any value
       Valeur psi_bound (mp.get_mg()->get_angu() )  ;
       
       int nnp = mp.get_mg()->get_np(1) ;
       int nnt = mp.get_mg()->get_nt(1) ;
       
       psi_bound = 1 ;
			
       for (int k=0 ; k<nnp ; k++)
	 for (int j=0 ; j<nnt ; j++)
	   psi_bound.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
       
       psi_bound.std_base_scal() ;
       


       // BC beta
       //-------
       const Coord& y = mp.y ;
       const Coord& x = mp.x ;

       Scalar xx(mp) ;
       xx = x ;
       xx.std_spectral_base() ;

       Scalar yy(mp) ;
       yy = y ;
       yy.std_spectral_base() ;

       Vector tmp_vect = nntilde_j * psi_j_a/(psi_j * psi_j) * met_gamt.radial_vect() ;
       tmp_vect.change_triad(mp.get_bvect_cart() ) ;

       // Beta-x
     
       Valeur lim_x (mp.get_mg()->get_angu()) ;
    
       lim_x = 1 ;  // Juste pour affecter dans espace des configs ;
       
       for (int k=0 ; k<nnp ; k++)
	 for (int j=0 ; j<nnt ; j++)
	  lim_x.set(0, k, j, 0) =  omega * yy.val_grid_point(1, k, j, 0) 
	    + tmp_vect(1).val_grid_point(1, k, j, 0) ;
       
       lim_x.set_base(*(mp.get_mg()->std_base_vect_cart()[0])) ;
       
       
      // Beta-y
     
       Valeur lim_y (mp.get_mg()->get_angu()) ;
    
       lim_y = 1 ;  // Juste pour affecter dans espace des configs ;
       
       for (int k=0 ; k<nnp ; k++)
	 for (int j=0 ; j<nnt ; j++)
	  lim_y.set(0, k, j, 0) =  - omega * xx.val_grid_point(1, k, j, 0) 
	    + tmp_vect(2).val_grid_point(1, k, j, 0) ;
       
       lim_y.set_base(*(mp.get_mg()->std_base_vect_cart()[1])) ;
       
     // Beta-z
     
       Valeur lim_z (mp.get_mg()->get_angu()) ;
    
       lim_z = 1 ;  // Juste pour affecter dans espace des configs ;
       
       for (int k=0 ; k<nnp ; k++)
	 for (int j=0 ; j<nnt ; j++)
	  lim_z.set(0, k, j, 0) = tmp_vect(3).val_grid_point(1, k, j, 0) ;
       
       lim_z.set_base(*(mp.get_mg()->std_base_vect_cart()[2])) ;
     
       // BC nn_tilde
       //------------
       
       // BC kappa=const
       
       Scalar  kappa (mp) ;
       kappa = 0.15 ;
       kappa.std_spectral_base() ;
       kappa.inc_dzpuis(2) ;

       tmp = - a / psi_j * nntilde_j *contract(met_gamt.radial_vect(), 0, psi_j.derive_cov(met_gamt), 0) 
	     + 1./2. * psi_j * psi_j * psi_j_ma * contract( beta_j_d.ope_killing_conf(met_gamt), 0, 1,  met_gamt.radial_vect()* met_gamt.radial_vect(), 0, 1)
	     - 1./3. * nntilde_j * psi_j * psi_j * trK * contract( met_gamt.cov(), 0, 1,  met_gamt.radial_vect()* met_gamt.radial_vect(), 0, 1) 
             + psi_j * psi_j * psi_j_ma * kappa 
  	     - met_gamt.radial_vect()(2) * nntilde_j.derive_cov(ff)(2)  - met_gamt.radial_vect()(3) * nntilde_j.derive_cov(ff)(3)  ;
       
       tmp = tmp / met_gamt.radial_vect()(1) ;

       // in this case you don't have to substract any value
       Valeur nn_bound_kappa (mp.get_mg()->get_angu() )  ;
       
       nn_bound_kappa = 1 ;
			
       for (int k=0 ; k<nnp ; k++)
	 for (int j=0 ; j<nnt ; j++)
	   nn_bound_kappa.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
       
       nn_bound_kappa.std_base_scal() ;
       

       
       // BC nn = const
       
       tmp = 0.2 * psi_j_ma - 1;

       // in this case you don't have to substract any value
       Valeur nn_bound_const (mp.get_mg()->get_angu() )  ;
       
       nn_bound_const = 1 ;
			
       for (int k=0 ; k<nnp ; k++)
	 for (int j=0 ; j<nnt ; j++)
	   nn_bound_const.set(0, k, j, 0) = tmp.val_grid_point(1, k, j, 0) ;
       
       nn_bound_const.std_base_scal() ;
       


       //============================
       // Resolution of the equations
       //============================

       Scalar psi_jp1(mp) ;
       Scalar nntilde_jp1(mp) ;
       Vector beta_jp1(beta_j) ;

       
       if (solve_psi == 1) {
	  psi_jp1 = sour_psi.poisson_neumann(psi_bound, 0) + 1. ;
       
       // Test:
	  maxabs(psi_jp1.laplacian() - sour_psi,
		 "Absolute error in the resolution of the equation for Psi") ;  
	  // Relaxation (relax=1 -> new ; relax=0 -> old )  
	  psi_jp1 = relax_psi * psi_jp1 + (1 - relax_psi) * psi_j ;
       }
       

       if (solve_lapse == 1) {
	 //nntilde_jp1 = sour_nntilde.poisson_neumann(nn_bound_kappa, 0) + 1. ;
	  nntilde_jp1 = sour_nntilde.poisson_dirichlet(nn_bound_const, 0) + 1. ;
       
       // Test:
	  maxabs(nntilde_jp1.laplacian() - sour_nntilde,
		 "Absolute error in the resolution of the equation for nntilde") ;  
	  // Relaxation (relax=1 -> new ; relax=0 -> old )  
	  nntilde_jp1 = relax_nn * nntilde_jp1 + (1 - relax_nn) * nntilde_j ;
       }

       if (solve_shift == 1) {
        double precision = 1e-8 ;
	poisson_vect_boundary(lambda, sour_beta, beta_jp1, lim_x, 
				  lim_y, lim_z, 0, precision, 20) ;
	    
	// Test
	sour_beta.dec_dzpuis() ;
	maxabs(beta_jp1.derive_con(ff).divergence(ff) 
	       + lambda * beta_jp1.divergence(ff)
	       .derive_con(ff) - sour_beta,
	       "Absolute error in the resolution of the equation for beta") ;  
	    
	cout << endl ;
	
	// Relaxation (relax=1 -> new ; relax=0 -> old )  
	beta_jp1 = relax_beta * beta_jp1 + (1 - relax_beta) * beta_j ;
       }



       //====================
       // Convergence control
       //====================

       double diff_nn, diff_psi, diff_beta ;
       diff_nn = 1.e-16 ;
       diff_psi = 1.e-16 ;
       diff_beta = 1.e-16 ;
       if (solve_lapse == 1)
	 diff_nn = max( diffrel(nntilde_j, nntilde_jp1) ) ;   
       if (solve_psi == 1)
	 diff_psi = max( diffrel(psi_j, psi_jp1) ) ; 
       if (solve_shift == 1)
	  diff_beta = max( maxabs(beta_jp1 - beta_j) ) ; 
       
       cout << "step = " << mer << " :  diff_psi = " << diff_psi 
	    << "  diff_nn = " << diff_nn 
	    << "  diff_beta = " << diff_beta << endl ;
       cout << "----------------------------------------------" << endl ;
       if ((diff_psi<precis) && (diff_nn<precis) && (diff_beta<precis))
	 break ; 
       
       if (mer>0) {conv << mer << "  " << log10(diff_nn) << " " << log10(diff_psi) 
			<< " " << log10(diff_beta) << endl ; } ;



       //======================
       // Updates for next step
       //======================

       psi_j = psi_jp1 ;
       nntilde_j = nntilde_jp1 ;
       beta_j = beta_jp1 ;


       //========================
       // Savings of output files
       //========================

       
	Scalar kkss (   psi_j_ma/(2. * nntilde_j ) * contract( beta_j_d.ope_killing_conf(met_gamt), 0, 1,  met_gamt.radial_vect()* met_gamt.radial_vect(), 0, 1)
			+ 1./3. * trK * contract( met_gamt.cov(), 0, 1,  met_gamt.radial_vect()* met_gamt.radial_vect(), 0, 1) ) ; 

	double max_kss = kkss.val_grid_point(1, 0, 0, 0) ;
	double min_kss = kkss.val_grid_point(1, 0, 0, 0) ;
	
	Scalar hh_tilde (contract(met_gamt.radial_vect().derive_cov(met_gamt), 0, 1)) ;
	double max_hh_tilde = hh_tilde.val_grid_point(1, 0, 0, 0) ;
	double min_hh_tilde = hh_tilde.val_grid_point(1, 0, 0, 0) ;
	
	
	for (int k=0 ; k<nnp ; k++)
	  for (int j=0 ; j<nnt ; j++){
	    if (kkss.val_grid_point(1, k, j, 0) > max_kss)
	      max_kss = kkss.val_grid_point(1, k, j, 0) ;
	    if (kkss.val_grid_point(1, k, j, 0) < min_kss)
	      min_kss = kkss.val_grid_point(1, k, j, 0) ;

	    if (hh_tilde.val_grid_point(1, k, j, 0) > max_hh_tilde)
	      max_hh_tilde = hh_tilde.val_grid_point(1, k, j, 0) ;
	    if (hh_tilde.val_grid_point(1, k, j, 0) < min_hh_tilde)
	      min_hh_tilde = hh_tilde.val_grid_point(1, k, j, 0) ;
	  }
     }
	       
     


     // Global updates
     //--------------
     Scalar psi_j_a = pow(psi_j, a) ;
     psi_j_a.std_spectral_base() ;
       
     
     if (solve_psi == 1)
       set_psi(psi_j) ; 
     if (solve_lapse == 1)
       n_evol.update(nntilde_j * psi_j_a, jtime, ttime) ; 
     if (solve_shift == 1)
	    beta_evol.update(beta_j, jtime, ttime) ;	

     if (solve_shift == 1)
       update_aa() ;
   
     Vector beta_j_d = beta().down(0, met_gamt) ;
     Scalar check ( psi() * psi() * psi() * trK
		    - psi() * met_gamt.radial_vect().divergence(met_gamt) 
		    - pow(psi(), -a)/(2. * nntilde_j * psi()) 
		    * contract( beta_j_d.ope_killing_conf(met_gamt), 0, 1,  met_gamt.radial_vect()* met_gamt.radial_vect(), 0, 1)
		    - 1./3. * trK/psi() * contract( met_gamt.cov(), 0, 1,  met_gamt.radial_vect()* met_gamt.radial_vect(), 0, 1) 
		    - 4 *contract (met_gamt.radial_vect(), 0, psi().derive_cov(met_gamt),0)) ;
     
     
       des_meridian(check, 1, 4., "CHECK", 1) ;
       arrete() ;

     conv.close() ;   
     kss.close() ;

}


