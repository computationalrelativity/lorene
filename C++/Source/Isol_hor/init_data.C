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
 * Revision 1.18  2005/06/09 08:05:32  f_limousin
 * Implement a new function sol_elliptic_boundary() and
 * Vector::poisson_boundary(...) which solve the vectorial poisson
 * equation (method 6) with an inner boundary condition.
 *
 * Revision 1.17  2005/05/12 14:48:07  f_limousin
 * New boundary condition for the lapse : boundary_nn_lapl().
 *
 * Revision 1.16  2005/04/08 12:16:52  f_limousin
 * Function set_psi(). And dependance in phi.
 *
 * Revision 1.15  2005/04/03 19:48:22  f_limousin
 * Implementation of set_psi(psi_in). And minor changes to avoid warnings.
 *
 * Revision 1.14  2005/04/02 15:49:21  f_limousin
 * New choice (Lichnerowicz) for aaquad. New member data nz.
 *
 * Revision 1.13  2005/03/31 09:45:31  f_limousin
 * New functions compute_ww(...) and aa_kerr_ww().
 *
 * Revision 1.12  2005/03/24 16:50:28  f_limousin
 * Add parameters solve_shift and solve_psi in par_isol.d and in function
 * init_dat(...). Implement Isolhor::kerr_perturb().
 *
 * Revision 1.11  2005/03/22 13:25:36  f_limousin
 * Small changes. The angular velocity and A^{ij} are computed
 * with a differnet sign.
 *
 * Revision 1.10  2005/03/09 10:18:08  f_limousin
 * Save K_{ij}s^is^j in a file. Add solve_lapse in a file
 *
 * Revision 1.9  2005/03/06 16:56:13  f_limousin
 * The computation of A^{ij} is no more necessary here thanks to the new
 * function Isol_hor::aa().
 *
 * Revision 1.8  2005/03/04 17:04:57  jl_jaramillo
 * Addition of boost to the shift after solving the shift equation
 *
 * Revision 1.7  2005/03/03 10:03:55  f_limousin
 * The boundary conditions for the lapse, psi and shift are now
 * parameters (in file par_hor.d).
 *
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
#include "isol_hor.h"
#include "metric.h"
#include "unites.h"
#include "graphique.h"
#include "cmp.h"
#include "tenseur.h"
#include "utilitaires.h"
#include "param.h"

void Isol_hor::init_data(int bound_nn, double lim_nn, int bound_psi, 
			 int bound_beta, int solve_lapse, int solve_psi,
			 int solve_shift, double precis, 
			 double relax, int niter) {

    using namespace Unites ;
   
    // Initialisations
    // ---------------
    double ttime = the_time[jtime] ;    

    ofstream conv("resconv.d") ; 
    ofstream kss("kss.d") ;
    conv << " # diff_nn   diff_psi   diff_beta " << endl ;

    // Iteration
    // ---------
    for (int mer=0; mer<niter; mer++) {
	
	//=========================================================
	// Boundary conditions and resolution of elliptic equations
	//=========================================================

	// Resolution of the Poisson equation for the lapse
	// ------------------------------------------------



	Scalar sou_nn (source_nn()) ;
	Scalar nn_jp1 (mp) ;
	if (solve_lapse == 1) {
	    Valeur nn_bound (mp.get_mg()-> get_angu()) ;

	    switch (bound_nn) {
		
		case 0 : {
		    nn_bound = boundary_nn_Dir(lim_nn) ;
		    nn_jp1 = sou_nn.poisson_dirichlet(nn_bound, 0) + 1. ;
		    break ;
		}
		case 1 : {
		    nn_bound = boundary_nn_Neu_eff(lim_nn) ;
		    nn_jp1 = sou_nn.poisson_neumann(nn_bound, 0) + 1. ;
		    break ;
		}
		case 2 : {
		    nn_bound = boundary_nn_Dir_eff(lim_nn) ;
		    nn_jp1 = sou_nn.poisson_dirichlet(nn_bound, 0) + 1. ;
		    break ;
		}
		case 3 : {
		    nn_bound = boundary_nn_Neu_kk() ;
		    nn_jp1 = sou_nn.poisson_neumann(nn_bound, 0) + 1. ;
		    break ;
		}
		case 4 : {
		    nn_bound = boundary_nn_Dir_kk() ;
		    nn_jp1 = sou_nn.poisson_dirichlet(nn_bound, 0) + 1. ;
		    break ;
		}
		case 5 : {
		    nn_bound = boundary_nn_Dir_lapl() ;
		    nn_jp1 = sou_nn.poisson_dirichlet(nn_bound, 0) + 1. ;
		    break ;
		}
		default : {
		    cout <<"Unexpected type of boundary conditions for the lapse!" 
			 << endl 
			 << "  bound_nn = " << bound_nn << endl ; 
		    abort() ;
		    break ; 
	    }
		    
	    } // End of switch  
	    
	// Test:
	    maxabs(nn_jp1.laplacian() - sou_nn,
		   "Absolute error in the resolution of the equation for N") ;
  	    
	    // Relaxation (relax=1 -> new ; relax=0 -> old )  
	    if (mer==0)
		n_evol.update(nn_jp1, jtime, ttime) ; 
	    else
		nn_jp1 = relax * nn_jp1 + (1 - relax) * nn() ;
	}
	    
	    
	// Resolution of the Poisson equation for Psi
	// ------------------------------------------
	
	Scalar sou_psi (source_psi()) ;
	Scalar psi_jp1 (mp) ;
	if (solve_psi == 1) {
	    Valeur psi_bound (mp.get_mg()-> get_angu()) ;
	    
	    switch (bound_psi) {
		
		case 0 : {
		    psi_bound = boundary_psi_app_hor() ;
		    psi_jp1 = sou_psi.poisson_neumann(psi_bound, 0) + 1. ;
		    break ;
		}
		case 1 : {
		    psi_bound = boundary_psi_Neu_spat() ;
		    psi_jp1 = sou_psi.poisson_neumann(psi_bound, 0) + 1. ;
		    break ;
		}
		case 2 : {
		    psi_bound = boundary_psi_Dir_spat() ;
		    psi_jp1 = sou_psi.poisson_dirichlet(psi_bound, 0) + 1. ;
		    break ;
		}
		case 3 : {
		    psi_bound = boundary_psi_Neu_evol() ;
		    psi_jp1 = sou_psi.poisson_neumann(psi_bound, 0) + 1. ;
		    break ;
		}
		case 4 : {
		    psi_bound = boundary_psi_Dir_evol() ;
		    psi_jp1 = sou_psi.poisson_dirichlet(psi_bound, 0) + 1. ;
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
	    maxabs(psi_jp1.laplacian() - sou_psi,
		   "Absolute error in the resolution of the equation for Psi") ;  
	    // Relaxation (relax=1 -> new ; relax=0 -> old )  
	    psi_jp1 = relax * psi_jp1 + (1 - relax) * psi() ;
	}
	
	// Resolution of the vector Poisson equation for the shift
	//---------------------------------------------------------	

	// Source

	Vector beta_jp1(beta()) ;

	if (solve_shift == 1) {
	    Vector source_vector ( source_beta() ) ;
	    double lambda = 0. ;
	    Vector source_reg = - (1./3. - lambda) * beta().divergence(ff)
		.derive_con(ff) ;
	    source_reg.inc_dzpuis() ;
	    source_vector = source_vector + source_reg ;

/*	    
	   // CARTESIAN CASE 
	   // #################################

	    // Boundary values
	    	    
	    Valeur boundary_x (mp.get_mg()-> get_angu()) ;
	    Valeur boundary_y (mp.get_mg()-> get_angu()) ;
	    Valeur boundary_z (mp.get_mg()-> get_angu()) ; 
	    
	    switch (bound_beta) {
		
		case 0 : {
		  boundary_x = boundary_beta_x(omega) ;
		  boundary_y = boundary_beta_y(omega) ;
		  boundary_z = boundary_beta_z() ;
		    break ;
		}
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
	    
*/
	    
	    // SPHERICAL CASE 
	    // #################################

	    // Boundary values

	    Valeur boundary_r (mp.get_mg()-> get_angu()) ;
	    Valeur boundary_t (mp.get_mg()-> get_angu()) ;
	    Valeur boundary_p (mp.get_mg()-> get_angu()) ; 
	    
	    switch (bound_beta) {
		
		case 0 : {
		  boundary_r = boundary_beta_r() ;
		  boundary_t = boundary_beta_theta() ;
		  boundary_p = boundary_beta_phi(omega) ;
		    break ;
		}
		case 1 : {
		    boundary_r = boundary_vv_x(omega) ;
		    boundary_t = boundary_vv_y(omega) ;
		    boundary_p = boundary_vv_z(omega) ;
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

	    // Resolution
	    //-----------
	    
	    beta_jp1 = source_vector.poisson_dirichlet(lambda, boundary_r,
				  boundary_t, boundary_p, 0) ;
	    

	    des_meridian(beta_jp1(1), 1.0000001, 10., "beta_r", 0) ;
	    des_meridian(beta_jp1(2), 1.0000001, 10., "beta_t", 1) ;
	    des_meridian(beta_jp1(3), 1.0000001, 10., "beta_p", 2) ;
	    arrete() ;

	    // Test
	    source_vector.dec_dzpuis() ;
	    maxabs(beta_jp1.derive_con(ff).divergence(ff) 
		   + lambda * beta_jp1.divergence(ff)
		   .derive_con(ff) - source_vector,
		   "Absolute error in the resolution of the equation for beta") ;  
	    
	    cout << endl ;
		    
	    // Boost
	    // -----
	    
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
	    beta_jp1 = relax * beta_jp1 + (1 - relax) * beta() ;
	}
	    
	//===========================================
	//      Convergence control
	//===========================================
	
	double diff_nn, diff_psi, diff_beta ;
	diff_nn = 1.e-16 ;
	diff_psi = 1.e-16 ;
	diff_beta = 1.e-16 ;
	if (solve_lapse == 1)
	  diff_nn = max( diffrel(nn(), nn_jp1) ) ;   
	if (solve_psi == 1)
	  diff_psi = max( diffrel(psi(), psi_jp1) ) ; 
	if (solve_shift == 1)
	  diff_beta = max( maxabs(beta_jp1 - beta()) ) ; 
	
	cout << "step = " << mer << " :  diff_psi = " << diff_psi 
	     << "  diff_nn = " << diff_nn 
	     << "  diff_beta = " << diff_beta << endl ;
	cout << "----------------------------------------------" << endl ;
	if ((diff_psi<precis) && (diff_nn<precis) && (diff_beta<precis))
	    break ; 
	
	conv << mer << "  " << log10(diff_nn) << " " << log10(diff_psi) 
	     << " " << log10(diff_beta) << endl ;
    
	//=============================================
	//      Updates for next step 
	//=============================================
	
	
	if (solve_psi == 1)
	    set_psi(psi_jp1) ; 
	if (solve_lapse == 1)
	    n_evol.update(nn_jp1, jtime, ttime) ; 
	if (solve_shift == 1)
	    beta_evol.update(beta_jp1, jtime, ttime) ;	

	if (solve_shift == 1)
	  update_aa() ;

	// Saving ok K_{ij}s^is^j
	// -----------------------
	
	Scalar kkss (contract(k_dd(), 0, 1, gam().radial_vect()*
		     gam().radial_vect(), 0, 1)) ;
	double max_kss = kkss.val_grid_point(1, 0, 0, 0) ;
	int nnp = mp.get_mg()->get_np(1) ;
	int nnt = mp.get_mg()->get_nt(1) ;
	for (int k=0 ; k<nnp ; k++)
	    for (int j=0 ; j<nnt ; j++)
		if (kkss.val_grid_point(1, k, j, 0) > max_kss)
		    max_kss = kkss.val_grid_point(1, k, j, 0) ;

	kss << mer << " " << max_kss << endl ;
    }
    
    conv.close() ;   
    kss.close() ;

} 



