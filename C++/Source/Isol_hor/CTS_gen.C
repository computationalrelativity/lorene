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

     double lambda = 1./3. ;              // to be used in the beta source

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
       Scalar psi_j_ma = pow(psi_j, -a) ;
       psi_j_ma.std_spectral_base() ;
       Vector beta_j_d = beta_j.down(0, met_gamt) ;
       Scalar nnpsia_j = nntilde_j *  psi_j_ma /( psi_j*psi_j*psi_j*psi_j*psi_j*psi_j) ;    // to be used in the beta source


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
       
       temp_scal_0 = 1./12. * trK*trK * psi_j*psi_j*psi_j*psi_j*psi_j ;
       temp_scal_0.inc_dzpuis(4) ;
       temp_scal_3 = 1./8. * met_gamt.ricci_scal() * psi_j ;
       temp_scal_3.inc_dzpuis() ;
       temp_scal_4 = -1./32. * psi_j*psi_j*psi_j*psi_j*psi_j * psi_j_ma*psi_j_ma / (nntilde_j*nntilde_j) *
	 contract(beta_j_d.ope_killing_conf(met_gamt), 0, 1, beta_j.ope_killing_conf(met_gamt), 0, 1) ;
       
       sour_psi = temp_scal_0 + temp_scal_3 + temp_scal_4 ;


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
       temp_scal_0 = - nntilde_j * (a - 4)/12. * psi_j*psi_j*psi_j*psi_j * trK*trK 
	             -  psi_j*psi_j*psi_j*psi_j * psi_j_ma * trK_point ;
       temp_scal_0.inc_dzpuis(4) ;

       temp_scal_2 =  psi_j*psi_j*psi_j*psi_j * psi_j_ma * contract (beta_j, 0, trK.derive_cov(ff), 0) ;
       temp_scal_2.inc_dzpuis(2) ;

       temp_scal_3 = - a/8. * nntilde_j * met_gamt.ricci_scal() ;
       temp_scal_3.inc_dzpuis() ;
       
       temp_scal_4 = - 2 * (a+1) * contract(psi_j.derive_cov(ff), 0, nntilde_j.derive_con(ff), 0) / psi_j 
	 - a*(a+1) * nntilde_j * contract(psi_j.derive_cov(ff), 0, psi_j.derive_con(ff), 0) / (psi_j * psi_j) 
	 + (a+8)/32. * psi_j*psi_j*psi_j*psi_j * psi_j_ma*psi_j_ma/nntilde_j 
	 * contract(beta_j_d.ope_killing_conf(met_gamt), 0, 1, beta_j.ope_killing_conf(met_gamt), 0, 1) ;

       sour_nntilde = temp_scal_0 + temp_scal_2 + temp_scal_3 + temp_scal_4 ;


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
       temp_vect_2 = 4./3. * psi_j_a * nntilde_j * trK.derive_con(met_gamt) ;
       temp_vect_2.inc_dzpuis(2) ;

       temp_vect_3 = - contract(met_gamt.ricci(), 1, beta_j, 0).up(0, met_gamt) 
	             + (lambda - 1./3.) * contract (beta_j.derive_cov(ff).derive_con(met_gamt), 0, 1) ;
       temp_vect_3.inc_dzpuis() ;

       temp_vect_4 = contract(beta_j.ope_killing_conf(met_gamt), 1, nnpsia_j.derive_cov(ff), 0) / nnpsia_j ;
                     
       sour_beta = temp_vect_2 + temp_vect_3 + temp_vect_4 ;

       // Correction Source 



       cout << " dzpuis sour_beta =" << sour_beta(1).get_dzpuis() << endl;


       arrete() ;





       

       //====================
       // Boundary conditions
       //====================



       //============================
       // Resolution of the equations
       //============================


       //====================
       // Convergence control
       //====================



       //======================
       // Updates for next step
       //======================

       

       //========================
       // Savings of output files
       //========================





     }

     conv.close() ;   
     kss.close() ;

}


