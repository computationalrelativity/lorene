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

char init_data_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2004/12/22 18:15:30  f_limousin
 * Many different changes.
 *
 * Revision 1.5  2004/11/08 14:51:21  f_limousin
 * A regularisation for the computation of A^{ij } is done in the
 * case lapse equal to zero on the horizon.
 *
 * Revision 1.1  2004/10/29 12:54:53  jl_jaramillo
 * First version
 *
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

void Isol_hor::init_data( double precis, double relax, int niter, 
			  double ang_vel) {

    using namespace Unites ;
   
    // Initialisations
    // ---------------
    double ttime = the_time[jtime] ;    

    ofstream conv("resconv.d") ; 
    conv << " # diff_nn   diff_psi   diff_beta " << endl ;


    // Iteration
    // ---------
    for (int mer=0; mer<niter; mer++) {
	
	//============================================
	//          Boundary conditions
	//============================================
      
	// Boundary for N	
	//---------------
	
	Valeur nn_bound (boundary_nn_Neu_kk()) ;
	
	// Boundary for Psi (Neumann)
	//---------------------------
	
	Valeur psi_bound (boundary_psi_Neu_spat()) ;

      //=============================================
      // Resolution of elliptic equations
      //=============================================
      
      // Resolution of the Poisson equation for the lapse
      // ------------------------------------------------

	Scalar nn_jp1 = source_nn().poisson_neumann(nn_bound, 0) + 1. ;
 
	// Test:
	maxabs(nn_jp1.laplacian() - source_nn(),
	       "Absolute error in the resolution of the equation for N") ;  
	
	// Relaxation (relax=1 -> new ; relax=0 -> old )  
	if (mer==0)
	n_evol.update(nn_jp1, jtime, ttime) ; 
	else
	nn_jp1 = relax * nn_jp1 + (1 - relax) * nn() ;

		
	// Resolution of the Poisson equation for Psi
	// ------------------------------------------

	Scalar psi_jp1 = source_psi().poisson_neumann(psi_bound, 0) + 1. ;
           
	// Test:
	maxabs(psi_jp1.laplacian() - source_psi(),
	       "Absolute error in the resolution of the equation for Psi") ;  

	// Relaxation (relax=1 -> new ; relax=0 -> old )  
	psi_jp1 = relax * psi_jp1 + (1 - relax) * psi() ;

	
	// Resolution of the vector Poisson equation for the shift
	//---------------------------------------------------------
	
	// Source
	
	Vector beta_jp1(beta()) ;
	Vector source_vector ( source_beta() ) ;
	double lambda = 0. ;
	Vector source_reg = - (1./3. - lambda) * beta().divergence(ff)
	    .derive_con(ff) ;
	source_reg.inc_dzpuis() ;
	source_vector = source_vector + source_reg ;
	
	// Boundary values
	
	Valeur boundary_x ( boundary_beta_x(ang_vel) ) ;
	Valeur boundary_y ( boundary_beta_y(ang_vel) ) ;
	Valeur boundary_z ( boundary_beta_z(ang_vel) ) ;
	
	// Resolution
	//-----------
	
	double precision = 1e-8 ;
	poisson_vect_boundary(lambda, source_vector, beta_jp1, boundary_x, 
			      boundary_y, boundary_z, 0, precision, 20) ;
	
	// Test
	source_vector.dec_dzpuis() ;
	maxabs(beta_jp1.derive_con(ff).divergence(ff) 
	       + lambda * beta_jp1.divergence(ff)
	       .derive_con(ff) - source_vector,
	       "Absolute error in the resolution of the equation for beta") ;  
	cout << endl ;
	
	// Relaxation (relax=1 -> new ; relax=0 -> old )  
	beta_jp1 = relax * beta_jp1 + (1 - relax) * beta() ;

	
	//===========================================
	//      Convergence control
	//===========================================
	
	double diff_nn = max( diffrel(nn(), nn_jp1) ) ;   
	double diff_psi = max( diffrel(psi(), psi_jp1) ) ; 
	double diff_beta = max( maxabs(beta_jp1 - beta()) ) ; 
	
	cout << "step = " << mer << " :  diff_psi = " << diff_psi 
	     << "  diff_nn = " << diff_nn 
	     << "  diff_beta = " << diff_beta << endl ;
	cout << "----------------------------------------------" << endl ;
	if ( (diff_psi < precis) && (diff_nn < precis) && (diff_beta < precis))
	    break ; 
	
	conv << mer << "  " << log(diff_nn) << " " << log(diff_psi) 
	     << " " << log(diff_beta) << endl ;
	
	//=============================================
	//      Updates for next step 
	//=============================================
	
	set_psi_del_q(psi_jp1) ; 
	n_evol.update(nn_jp1, jtime, ttime) ; 
	beta_evol.update(beta_jp1, jtime, ttime) ; 
	
	// New value of A^{ij}:
	
	Sym_tensor aa_jp1 (mp, CON, mp.get_bvect_spher()) ;
	int nnr = mp.get_mg()->get_nr(1) ;
	int nnt = mp.get_mg()->get_nt(1) ;
	int nnp = mp.get_mg()->get_np(1) ;
	
	int check ;
	check = 0 ;
	for (int k=0; k<nnp; k++)
	    for (int j=0; j<nnt; j++){
		if (nn().val_grid_point(1, k, j , 0) < 1e-12){
		    check = 1 ;
		    break ;
		}
	    }
	
	if (check == 0)
	    aa_jp1 = ( beta().ope_killing_conf(tgam()) + gamt_point ) 
		/ (2.* nn()) ;            
	else {
	    Scalar nn_sxpun (division_xpun (Cmp(nn()), 0)) ;
	    
	    Sym_tensor aa_sxpun = beta().ope_killing_conf(tgam())
		+ gamt_point ;

	    const Coord& r = mp.r ;        // r field
	    Mtbl r_mtbl = r ;
	    Scalar rr (mp) ;
	    rr = r_mtbl ;

	    double r1, r2 ;
	    r1 = mp.val_r(1, -1., 0., 0.) ;
	    r2 = mp.val_r(1, 1., 0., 0.) ;

	    for(int k=0; k<nnp; k++)
		for(int j=0; j<nnt; j++)
		    for(int m=1; m<=3; m++)
			for(int n=1; n<=m; n++){
			    double aa_mn_jk = 
				aa_sxpun(m,n).val_grid_point(1, k, j, 0) ;     
			    for(int i=0; i<nnr; i++){
				aa_sxpun.set(m,n).set_grid_point(1, k, j, i) -=
				aa_mn_jk * 
				(- 2 * rr.val_grid_point(1, k, j, i) + 3 * r1
				 - r2) * (rr.val_grid_point(1, k, j, i) - r2) *
				(rr.val_grid_point(1, k, j, i) - r2) /
				(r1 - r2) / (r1 - r2) / (r1 - r2) ;
			    }
			}

	    for(int i=1; i<=3; i++)
		for(int j=1; j<=i; j++){
		    aa_sxpun.set(i,j).set_inner_boundary(1, 0.) ;
		    aa_sxpun.set(i,j) = Scalar (division_xpun 
						(Cmp(aa_sxpun(i,j)), 0)) ;
		}
	    aa_jp1 = aa_sxpun / (2*nn_sxpun) ;
	}
	cout << "check = " << check << endl ;
	
	set_aa(aa_jp1) ; 	
    }   
    conv.close() ;   
} 



