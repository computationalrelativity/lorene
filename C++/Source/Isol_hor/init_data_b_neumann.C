/*
 *  Method of class Isol_hor to compute valid initial data for 
 *  Berlin boundary conditions
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

char init_data_b_neumann_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2005/03/04 17:05:34  jl_jaramillo
 * Addition of boost to the shift after solving the shift equation
 *
 * Revision 1.3  2005/03/03 10:04:16  f_limousin
 * The boundary conditions for the lapse, psi and shift are now
 * parameters (in file par_hor.d).
 *
 * Revision 1.2  2004/12/22 18:15:43  f_limousin
 * Many different changes.
 *
 * Revision 1.1  2004/11/24 19:29:51  jl_jaramillo
 * construction of initial data with Berlin boundary conditions
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

void Isol_hor::init_data_b_neumann(int bound_psi, 
				   int bound_beta,  double precis, 
				   double relax, int niter) {

    using namespace Unites ;
   
    // Initialisations
    // ---------------
    double ttime = the_time[jtime] ;    

    const Base_vect& triad = *(nn().get_triad()) ;
    
    Scalar tmp(mp) ;
    Scalar tmp_scal(mp) ; 
    Sym_tensor tmp_sym(mp, CON, triad) ;

    // Reset of quantities depending on K:
    k_dd_evol.downdate(jtime) ; 
    k_uu_evol.downdate(jtime) ; 

    set_aa( (beta().ope_killing_conf(met_gamt) +  gamt_point) / (2.* nn()) ) ; 

    ofstream conv("resconv.d") ; 
    conv << " # diff_psi   diff_beta " << endl ;


    // Iteration
    // ---------
    for (int mer=0; mer<niter; mer++) {
	
	//=========================================================
	// Boundary conditions and resolution of elliptic equations
	//=========================================================


	Valeur psi_bound (mp.get_mg()-> get_angu()) ;
	Scalar psi_jp1 (mp) ;

	switch (bound_psi) {
	    
	    case 0 : {
		psi_bound = boundary_psi_app_hor() ;
		psi_jp1 = source_psi().poisson_neumann(psi_bound, 0) + 1. ;
		break ;
	    }
	    case 1 : {
		psi_bound = boundary_psi_Neu_spat() ;
		psi_jp1 = source_psi().poisson_neumann(psi_bound, 0) + 1. ;
		break ;
	    }
	    case 2 : {
		psi_bound = boundary_psi_Dir_spat() ;
		psi_jp1 = source_psi().poisson_dirichlet(psi_bound, 0) + 1. ;
		break ;
	    }
	    case 3 : {
		psi_bound = boundary_psi_Neu_evol() ;
		psi_jp1 = source_psi().poisson_neumann(psi_bound, 0) + 1. ;
		break ;
	    }
	    case 4 : {
		psi_bound = boundary_psi_Dir_evol() ;
		psi_jp1 = source_psi().poisson_dirichlet(psi_bound, 0) + 1. ;
		break ;
	    }
	    default : {
		cout <<"Unexpected type of boundary conditions for psi!" 
		     << endl 
		     << "  bound_psi = " << bound_psi << endl ; 
		abort() ;
		break ; 
	    }
		
	} // End of switch  
           
	// Test:
	maxabs(psi_jp1.laplacian() - source_psi(),
	       "Absolute error in the resolution of the equation for Psi") ;  
	
	// Relaxation (relax=1 -> new ; relax=0 -> old )  
	//-----------
	psi_jp1 = relax * psi_jp1 + (1 - relax) * psi() ;
           

	// Resolution of the vector Poisson equation for the shift
	//---------------------------------------------------------
	
	Vector beta_jp1(beta()) ;

	// Source
	//-------
	
	Vector source_vector ( source_beta() ) ;
	double lambda = 0. ;
	Vector source_reg = - (1./3. - lambda) * beta().divergence(ff)
	    .derive_con(ff) ;
	source_reg.inc_dzpuis() ;
	source_vector = source_vector + source_reg ;
	
	// Boundary values
	//----------------

	
	Valeur boundary_x (mp.get_mg()-> get_angu()) ;
	Valeur boundary_y (mp.get_mg()-> get_angu()) ;
	Valeur boundary_z (mp.get_mg()-> get_angu()) ;

	switch (bound_beta) {
	    
	    case 1 : {
		boundary_x = boundary_vv_x(omega) ;
		boundary_y = boundary_vv_y(omega) ;
		boundary_z = boundary_vv_z(omega) ;
		break ;
	    }
	    default : {
		cout <<"Unexpected type of boundary conditions for psi!" 
		     << endl 
		     << "  bound_psi = " << bound_psi << endl ; 
		abort() ;
		break ; 
	    }
		
	} // End of switch  
	
	if (boost_x != 0.) 
	    boundary_x -= beta_boost_x() ;
	if (boost_z != 0.) 
	    boundary_z -= beta_boost_z() ;

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
	
	Vector boost_vect(mp, CON, mp.get_bvect_cart()) ;
	if (boost_x != 0.) {
	  boost_vect.set(1) = boost_x ;
	  boost_vect.set(2) = 0. ;
	  boost_vect.set(3) = 0. ;
	  boost_vect.std_spectral_base() ;
	  boost_vect.change_triad(mp.get_bvect_spher()) ;
	  beta_jp1 = beta_jp1 + boost_vect ;
	}
	  
	if (boost_z != 0.) {
	  boost_vect.set(1) = boost_z ;
	  boost_vect.set(2) = 0. ;
	  boost_vect.set(3) = 0. ;
	  boost_vect.std_spectral_base() ;
	  boost_vect.change_triad(mp.get_bvect_spher()) ;
	  beta_jp1 = beta_jp1 + boost_vect ;
	}



	// Relaxation (relax=1 -> new ; relax=0 -> old )  
	//-----------
	
	beta_jp1 = relax * beta_jp1 + (1 - relax) * beta() ;
	
	//===========================================
	//      Convergence control
	//===========================================
	
	double diff_psi = max( diffrel(psi(), psi_jp1) ) ; 
	double diff_beta = max( maxabs(beta_jp1 - beta()) ) ; 
	
	cout << "step = " << mer << " :  diff_psi = " << diff_psi 
	     << "  diff_beta = " << diff_beta << endl ;
	cout << "----------------------------------------------" << endl ;
	if ( (diff_psi < precis) && (diff_beta < precis))
	    break ; 

	
	conv << mer << "  " << log10(diff_psi) 
	     << " " << log10(diff_beta) << endl ;
	
	//=============================================
	//      Updates for next step 
	//=============================================
	
	set_psi_del_q(psi_jp1) ; 
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
	    aa_jp1 = ( beta().ope_killing_conf(met_gamt) + gamt_point ) 
		/ (2.* nn()) ;            
	else {
	    Scalar nn_sxpun (division_xpun (Cmp(nn()), 0)) ;
	    
	    Sym_tensor aa_sxpun = beta().ope_killing_conf(met_gamt)
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



