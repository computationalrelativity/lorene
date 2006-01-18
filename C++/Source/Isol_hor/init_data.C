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
 * Revision 1.25  2006/01/18 12:28:35  jl_jaramillo
 * new method (chi, Theta) for the spherical symmetric case
 *
 * Revision 1.24  2006/01/18 09:04:27  f_limousin
 * Minor modifs (warnings and errors at the compilation with gcc-3.4)
 *
 * Revision 1.23  2006/01/16 17:13:40  jl_jaramillo
 * function for solving the spherical case
 *
 * Revision 1.22  2005/11/02 16:09:44  jl_jaramillo
 * changes in boundary_nn_Dir_lapl
 *
 * Revision 1.21  2005/10/24 16:44:40  jl_jaramillo
 * Cook boundary condition ans minot bound of kss
 *
 * Revision 1.20  2005/10/21 16:20:55  jl_jaramillo
 * Version for the paper JaramL05
 *
 * Revision 1.19  2005/07/08 13:15:23  f_limousin
 * Improvements of boundary_vv_cart(), boundary_nn_lapl().
 * Add a fonction to compute the departure of axisymmetry.
 *
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
		    nn_bound = boundary_nn_Dir_lapl(mer) ;
		    nn_jp1 = sou_nn.poisson_dirichlet(nn_bound, 0) + 1. ;
		    break ;
		}
		case 6 : {
		    nn_bound = boundary_nn_Neu_Cook() ;
		    nn_jp1 = sou_nn.poisson_neumann(nn_bound, 0) + 1. ;
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
		case 5 : {
		    psi_bound = boundary_psi_Dir() ;
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
	    double lambda = 0.;
	    Vector source_reg = - (1./3. - lambda) * beta().divergence(ff)
		.derive_con(ff) ;
	    source_reg.inc_dzpuis() ;
	    source_vector = source_vector + source_reg ;

	    
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
	    

/*	    
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
	    // #########################
	    // End of spherical case
	    // #########################


*/
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
	
	if (mer>0) {conv << mer << "  " << log10(diff_nn) << " " << log10(diff_psi) 
			 << " " << log10(diff_beta) << endl ; } ;
    
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
	double min_kss = kkss.val_grid_point(1, 0, 0, 0) ;
	
	Scalar aaLss (pow(psi(), 6) * kkss) ;
	double max_aaLss = aaLss.val_grid_point(1, 0, 0, 0) ;
	double min_aaLss = aaLss.val_grid_point(1, 0, 0, 0) ;
	
	Scalar hh_tilde (contract(met_gamt.radial_vect().derive_cov(met_gamt), 0, 1)) ;
	double max_hh_tilde = hh_tilde.val_grid_point(1, 0, 0, 0) ;
	double min_hh_tilde = hh_tilde.val_grid_point(1, 0, 0, 0) ;
	
	
	int nnp = mp.get_mg()->get_np(1) ;
	int nnt = mp.get_mg()->get_nt(1) ;
	for (int k=0 ; k<nnp ; k++)
	  for (int j=0 ; j<nnt ; j++){
	    if (kkss.val_grid_point(1, k, j, 0) > max_kss)
	      max_kss = kkss.val_grid_point(1, k, j, 0) ;
	    if (kkss.val_grid_point(1, k, j, 0) < min_kss)
	      min_kss = kkss.val_grid_point(1, k, j, 0) ;

	    if (aaLss.val_grid_point(1, k, j, 0) > max_aaLss)
	      max_aaLss = aaLss.val_grid_point(1, k, j, 0) ;
	    if (aaLss.val_grid_point(1, k, j, 0) < min_aaLss)
	      min_aaLss = aaLss.val_grid_point(1, k, j, 0) ;
	    
	    if (hh_tilde.val_grid_point(1, k, j, 0) > max_hh_tilde)
	      max_hh_tilde = hh_tilde.val_grid_point(1, k, j, 0) ;
	    if (hh_tilde.val_grid_point(1, k, j, 0) < min_hh_tilde)
	      min_hh_tilde = hh_tilde.val_grid_point(1, k, j, 0) ;
	    
	  }
	
	
	kss << mer << " " << max_kss << " " << min_kss << " " << max_aaLss << " " << min_aaLss
	    << " " <<  -1 * max_hh_tilde << " "  << -1 * min_hh_tilde << endl ;
    }
    
    conv.close() ;   
    kss.close() ;

} 





void Isol_hor::init_data_spher(int bound_nn, double lim_nn, int bound_psi, 
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

      
      //       des_meridian(psi(), 1, 10., "psi", 0) ;
      //       des_meridian(b_tilde(), 1, 10., "b_tilde", 1) ;
      //       des_meridian(nn(), 1, 10., "nn", 2) ;
      //       arrete() ;
      

      //========
      // Sources
      //========

      // Useful functions
      // ----------------
       Vector tem_vect (beta() ) ;
       Scalar dif_b = b_tilde() - tem_vect.set(1) ;
       //       cout << "dif_b = " << dif_b << endl ;	
       //       arrete() ;
       
       Scalar dbdr ( b_tilde().dsdr() ) ;

       Scalar bsr (b_tilde()) ;
       bsr.div_r() ;
       bsr.inc_dzpuis(2) ;

       Scalar bsr2 ( bsr) ;
       bsr2.div_r() ;
       bsr2.inc_dzpuis(2)  ;
      
       Scalar psisr (psi()) ;
       psisr.div_r() ;
       psisr.inc_dzpuis(2) ;

      
     
       // Source Psi
       // ----------
       Scalar source_psi_spher(mp) ;
       source_psi_spher = -1./12. * psi4()*psi()/(nn() * nn()) * (dbdr - bsr) * (dbdr - bsr)   ;
       
       // Source N
       //---------
       Scalar source_nn_spher(mp) ;
       source_nn_spher = 2./3. * psi4() /nn() * (dbdr - bsr) * (dbdr - bsr)   
	 - 2 * ln_psi().dsdr() * nn().dsdr() ;
       
       // Source b_tilde
      //---------------
       Scalar source_btilde_spher(mp) ;
       
       Scalar tmp ( -1./3. * (dbdr + 2 * bsr).dsdr() ) ;
       tmp.std_spectral_base() ;
       tmp.inc_dzpuis() ;       

       source_btilde_spher = tmp + 2 * bsr2  
       	                     + 4./3. * (dbdr - bsr) * ( nn().dsdr()/nn()  - 6 * psi().dsdr()/psi() ) ;
    
       Scalar source_btilde_trun(mp) ;
       
       source_btilde_trun = tmp +
	 4./3. * (dbdr - bsr) * ( nn().dsdr()/nn()  - 6 * psi().dsdr()/psi() ) ;
       

       //       Scalar diff_dbeta ( (dbdr + 2 * bsr).dsdr() -  beta().divergence(ff).derive_con(ff)(1) ) ;
       

	
       // Parallel calculation
       //---------------------
       
       Scalar sourcepsi (source_psi()) ;
       Scalar sourcenn (source_nn()) ;
       
       Vector sourcebeta (source_beta()) ;
       Vector source_reg =  1./3. * beta().divergence(ff).derive_con(ff) ;
       source_reg.inc_dzpuis() ;
       sourcebeta -=  source_reg ;
       Scalar source_btilde (sourcebeta(1) ) ;
       
       //       Scalar diff_div =  source_reg(1) + tmp ;   ;
       
       Scalar mag_sou_psi ( source_psi_spher ) ;
       mag_sou_psi.dec_dzpuis(4) ;
       Scalar mag_sou_nn ( source_nn_spher ) ;
       mag_sou_nn.dec_dzpuis(4) ;
       Scalar mag_sou_btilde ( source_btilde_trun ) ;
       mag_sou_btilde.dec_dzpuis(4) ;
    
       Scalar diff_sou_psi ( source_psi_spher - sourcepsi) ;
       diff_sou_psi.dec_dzpuis(4) ;
       Scalar diff_sou_nn ( source_nn_spher - sourcenn) ;
       diff_sou_nn.dec_dzpuis(4) ;
       Scalar diff_sou_btilde ( source_btilde_trun - source_btilde) ;
       diff_sou_btilde.dec_dzpuis(4) ;
       
       /*
       cout << "dzpuis mag_btilde =" << mag_sou_btilde.get_dzpuis()<<endl  ;
       des_meridian(diff_sou_psi, 1, 10., "diff_psi", 0) ;
       des_meridian(diff_sou_nn, 1, 10., "diff_nn", 1) ;
       des_meridian(diff_sou_btilde, 1, 10., "diff_btilde", 2) ;
       des_meridian(mag_sou_psi, 1, 10., "mag_psi", 3) ;
       des_meridian(mag_sou_nn, 1, 10., "mag_nn", 4) ;
       des_meridian(mag_sou_btilde, 1, 10., "mag_btilde", 5) ;
       //       des_meridian(diff_dbeta, 1, 10., "diff_dbeta", 6) ;
       
       arrete() ;
       */

      
       //====================
       // Boundary conditions
       //====================
       
       // To avoid warnings;
       bound_nn = 1 ; lim_nn = 1. ; bound_psi = 1 ; bound_beta = 1 ;

       double kappa_0 = 0.2 - 1. ;
       
       Scalar kappa (mp) ;
       kappa = kappa_0 ;
       kappa.std_spectral_base() ;
       kappa.inc_dzpuis(2) ;
       
       
       int nnp = mp.get_mg()->get_np(1) ;
       int nnt = mp.get_mg()->get_nt(1) ;
       
       
       Valeur psi_bound (mp.get_mg()-> get_angu()) ;
       Valeur nn_bound (mp.get_mg()-> get_angu()) ;
       Valeur btilde_bound (mp.get_mg()-> get_angu()) ;
       psi_bound = 1. ; // Juste pour affecter dans espace des configs ;
       nn_bound = 1. ; // Juste pour affecter dans espace des configs ;
       btilde_bound = 1. ; // Juste pour affecter dans espace des configs ;
       
       Scalar tmp_psi = -1./4. * (2 * psisr +  
				  2./3. * psi4()/(psi() * nn()) * (dbdr - bsr) ) ;
       
       Scalar tmp_nn = kappa ; //+ 2./3. * psi() * psi() * (dbdr - bsr)   ;
       
       Scalar tmp_btilde = nn() / (psi() * psi()) ;
        

       for (int k=0 ; k<nnp ; k++)
	 for (int j=0 ; j<nnt ; j++){
	   psi_bound.set(0, k, j, 0) = tmp_psi.val_grid_point(1, k, j, 0) ;       // BC Psi
	   nn_bound.set(0, k, j, 0) = tmp_nn.val_grid_point(1, k, j, 0) ;         // BC N
	   btilde_bound.set(0, k, j, 0) = tmp_btilde.val_grid_point(1, k, j, 0) ; // BC b_tilde
	 }
       
       psi_bound.std_base_scal() ;
       nn_bound.std_base_scal() ;
       btilde_bound.std_base_scal() ;
       
       
       //=================================
       // Resolution of elliptic equations
       //=================================
       
       // Resolution of the Poisson equation for Psi
       // ------------------------------------------
       Scalar psi_jp1 (mp) ;
       if (solve_psi == 1) {
	 
	psi_jp1 = source_psi_spher.poisson_neumann(psi_bound, 0) + 1. ;
	
	// Test:
	maxabs(psi_jp1.laplacian() -  source_psi_spher,
	       "Absolute error in the resolution of the equation for Psi") ;  
	// Relaxation (relax=1 -> new ; relax=0 -> old )  
	psi_jp1 = relax * psi_jp1 + (1 - relax) * psi() ;
       }
       
      // Resolution of the Poisson equation for the lapse
      // ------------------------------------------------
       Scalar nn_jp1 (mp) ;
       if (solve_lapse == 1) {
	 
	 nn_jp1 = source_nn_spher.poisson_dirichlet(nn_bound, 0) + 1. ;
	
	 // Test:
	 maxabs(nn_jp1.laplacian() - source_nn_spher,
		"Absolute error in the resolution of the equation for N") ;
	 
	 // Relaxation (relax=1 -> new ; relax=0 -> old )  
	 if (mer==0)
	  n_evol.update(nn_jp1, jtime, ttime) ; 
	 else
	   nn_jp1 = relax * nn_jp1 + (1 - relax) * nn() ;
	 
       }
       
       // Resolution of the Poisson equation for b_tilde
      // ----------------------------------------------
       Scalar btilde_jp1 (mp) ;
       if (solve_shift == 1) {
	 
	btilde_jp1 = source_btilde_spher.poisson_dirichlet(btilde_bound, 0)  ;
	
	// Test:
	maxabs(btilde_jp1.laplacian() - source_btilde_spher,
	       "Absolute error in the resolution of the equation for btilde") ;  
	// Relaxation (relax=1 -> new ; relax=0 -> old )  
	btilde_jp1 = relax * btilde_jp1 + (1 - relax) * b_tilde() ;
       }
       
	    
       //===========================================
       //      Convergence control
       //===========================================
       
       double diff_nn, diff_psi, diff_btilde ;
       diff_nn = 1.e-16 ;
	diff_psi = 1.e-16 ;
	diff_btilde = 1.e-16 ;
	if (solve_lapse == 1)
	  diff_nn = max( diffrel(nn(), nn_jp1) ) ;   
	if (solve_psi == 1)
	  diff_psi = max( diffrel(psi(), psi_jp1) ) ; 
	if (solve_shift == 1)
	  diff_btilde = max( diffrel(btilde_jp1, b_tilde()) ) ; 
	
	cout << "step = " << mer << " :  diff_psi = " << diff_psi 
	     << "  diff_nn = " << diff_nn 
	     << "  diff_btilde = " << diff_btilde << endl ;
	cout << "----------------------------------------------" << endl ;
	if ((diff_psi<precis) && (diff_nn<precis) && (diff_btilde<precis))
	  break ; 
	
	if (mer>0) {conv << mer << "  " << log10(diff_nn) << " " << log10(diff_psi) 
			 << " " << log10(diff_btilde) << endl ; } ;
	
	//=============================================
	//      Updates for next step 
	//=============================================
	
	
	if (solve_psi == 1)
	  set_psi(psi_jp1) ; 
	if (solve_lapse == 1)
	  n_evol.update(nn_jp1, jtime, ttime) ; 
	if (solve_shift == 1)
	 { 
	   Vector beta_jp1 (btilde_jp1 * tgam().radial_vect()) ;
	   cout <<  tgam().radial_vect() << endl ;
	   beta_evol.update(beta_jp1, jtime, ttime) ;	
	 }
	if (solve_shift == 1 || solve_lapse == 1)
	  {
	    update_aa() ;
	  }

	// Saving ok K_{ij}s^is^j
	// -----------------------
	
	Scalar kkss (contract(k_dd(), 0, 1, gam().radial_vect()*
			      gam().radial_vect(), 0, 1)) ;
	double max_kss = kkss.val_grid_point(1, 0, 0, 0) ;
	double min_kss = kkss.val_grid_point(1, 0, 0, 0) ;
	
	Scalar aaLss (pow(psi(), 6) * kkss) ;
	double max_aaLss = aaLss.val_grid_point(1, 0, 0, 0) ;
	double min_aaLss = aaLss.val_grid_point(1, 0, 0, 0) ;
	
	Scalar hh_tilde (contract(met_gamt.radial_vect().derive_cov(met_gamt), 0, 1)) ;
	double max_hh_tilde = hh_tilde.val_grid_point(1, 0, 0, 0) ;
	double min_hh_tilde = hh_tilde.val_grid_point(1, 0, 0, 0) ;
	
	

	for (int k=0 ; k<nnp ; k++)
	  for (int j=0 ; j<nnt ; j++){
	    if (kkss.val_grid_point(1, k, j, 0) > max_kss)
	      max_kss = kkss.val_grid_point(1, k, j, 0) ;
	    if (kkss.val_grid_point(1, k, j, 0) < min_kss)
	      min_kss = kkss.val_grid_point(1, k, j, 0) ;

	    if (aaLss.val_grid_point(1, k, j, 0) > max_aaLss)
	      max_aaLss = aaLss.val_grid_point(1, k, j, 0) ;
	    if (aaLss.val_grid_point(1, k, j, 0) < min_aaLss)
	      min_aaLss = aaLss.val_grid_point(1, k, j, 0) ;
	    
	    if (hh_tilde.val_grid_point(1, k, j, 0) > max_hh_tilde)
	      max_hh_tilde = hh_tilde.val_grid_point(1, k, j, 0) ;
	    if (hh_tilde.val_grid_point(1, k, j, 0) < min_hh_tilde)
	      min_hh_tilde = hh_tilde.val_grid_point(1, k, j, 0) ;
	    
	  }
	
	
	kss << mer << " " << max_kss << " " << min_kss << " " << max_aaLss << " " << min_aaLss
	    << " " <<  -1 * max_hh_tilde << " "  << -1 * min_hh_tilde << endl ;
    }
    
    conv.close() ;   
    kss.close() ;

} 



void Isol_hor::init_data_alt(int bound_nn, double lim_nn, int bound_psi, 
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

    Scalar psi_j (psi()) ;
    Scalar nn_j (nn()) ;
    Scalar btilde_j (b_tilde()) ;    
    Scalar diffb (  btilde_j - b_tilde() ) ;
    Scalar theta_j (beta().divergence(ff)) ;
    theta_j.dec_dzpuis(2) ;
    Scalar chi_j (b_tilde()) ;
    chi_j.mult_r() ;

   
    // Iteration
    // ---------

    for (int mer=0; mer<niter; mer++) {

      
      des_meridian(psi_j, 1, 10., "psi", 0) ;
      des_meridian(nn_j, 1, 10., "nn", 1) ;
      des_meridian(theta_j, 1, 10., "Theta", 2) ;
      des_meridian(chi_j, 1, 10., "chi", 3) ;
      arrete() ;


      //========
      // Sources
      //========

      // Useful functions
      // ----------------
       Vector tem_vect (beta() ) ;
       //       Scalar dif_b = btilde_j - tem_vect.set(1) ;
       //       cout << "dif_b = " << dif_b << endl ;	
       //       arrete() ;
       
       Scalar dbdr ( btilde_j.dsdr() ) ;

       Scalar bsr (btilde_j) ;
       bsr.div_r() ;
       bsr.inc_dzpuis(2) ;

       Scalar bsr2 ( bsr) ;
       bsr2.div_r() ;
       bsr2.inc_dzpuis(2)  ;
      
       Scalar psisr (psi_j) ;
       psisr.div_r() ;
       psisr.inc_dzpuis(2) ;

       Scalar dchidr ( chi_j.dsdr() ) ;
       
       Scalar chisr (chi_j) ;
       chisr.div_r() ;
       chisr.inc_dzpuis(2) ;
      
       Scalar rdthetadr (theta_j.dsdr() ) ;
       rdthetadr.mult_r() ;
       rdthetadr.inc_dzpuis(2) ;

       Scalar theta_dz4 (theta_j) ;
       theta_dz4.inc_dzpuis(4) ; 

       Scalar dbmb (dchidr - 2 * chisr) ;
       dbmb.div_r() ;
     

       // Source Psi
       // ----------
       Scalar source_psi_spher(mp) ;

       source_psi_spher = -1./12. * psi_j*psi_j*psi_j*psi_j*psi_j/(nn_j * nn_j) 
       	 * dbmb *dbmb ; 

       
       // Source N
       //---------
       Scalar source_nn_spher(mp) ;
       source_nn_spher = 2./3. * psi_j*psi_j*psi_j*psi_j/nn_j * dbmb *dbmb
	 - 2 * psi_j.dsdr()/psi_j * nn_j.dsdr() ;



       // Source_theta
       //-------------
       Scalar source_theta_spher(mp) ;
       
       source_theta_spher =  (dbmb * (nn_j.dsdr()/nn_j - 6 * psi_j.dsdr()/psi_j)).dsdr() ; 
       source_theta_spher.dec_dzpuis() ;
       
       
       
       // Source chi
       //-----------
       Scalar source_chi_spher(mp) ;
       
       source_chi_spher = 4./3. * (dchidr - 2 * chisr) * ( nn_j.dsdr()/nn_j  - 6 * psi_j.dsdr()/psi_j ) 
                      	  - 1./3. * rdthetadr + 2 * theta_dz4 ;
			  

      
       //====================
       // Boundary conditions
       //====================
       double kappa_0 = 0.2 - 1. ;
       
       Scalar kappa (mp) ;
       kappa = kappa_0 ;
       kappa.std_spectral_base() ;
       kappa.inc_dzpuis(2) ;
       
       
       int nnp = mp.get_mg()->get_np(1) ;
       int nnt = mp.get_mg()->get_nt(1) ;
       
       
       Valeur psi_bound (mp.get_mg()-> get_angu()) ;
       Valeur nn_bound (mp.get_mg()-> get_angu()) ;
       Valeur theta_bound (mp.get_mg()-> get_angu()) ;
       Valeur chi_bound (mp.get_mg()-> get_angu()) ;

       psi_bound = 1. ; // Juste pour affecter dans espace des configs ;
       nn_bound = 1. ; // Juste pour affecter dans espace des configs ;
       theta_bound = 1. ; // Juste pour affecter dans espace des configs ;
       chi_bound = 1. ; // Juste pour affecter dans espace des configs ;

       //psi

       Scalar tmp_psi = -1./4. * (2 * psisr +  
				  2./3. * psi_j*psi_j*psi_j/ nn_j * dbmb ) ;

       //       tmp_psi = 2./3. * psi_j*psi_j*psi_j/ nn_j * (dchidr - 2 * chisr) ;
       //       tmp_psi.div_r() ; 
       //       tmp_psi =  -1./4. * (2 * psisr + tmp_psi) ;
       
       //nn
       Scalar tmp_nn = kappa ; //+ 2./3. * psi_j*psi_j * dbmb ;
      

       //theta
       Scalar tmp_theta = dchidr ;
       tmp_theta.dec_dzpuis(2) ;
       tmp_theta += nn_j/(psi_j*psi_j)  ;
       tmp_theta.div_r() ;
      

       //chi
       Scalar tmp_chi = nn_j/(psi_j*psi_j) ;
       tmp_chi.mult_r() ;
      

       for (int k=0 ; k<nnp ; k++)
	 for (int j=0 ; j<nnt ; j++){
	   psi_bound.set(0, k, j, 0) = tmp_psi.val_grid_point(1, k, j, 0) ;       // BC Psi
	   nn_bound.set(0, k, j, 0) = tmp_nn.val_grid_point(1, k, j, 0) ;         // BC N
	   theta_bound.set(0, k, j, 0) = tmp_theta.val_grid_point(1, k, j, 0) ;       // BC Theta
	   chi_bound.set(0, k, j, 0) = tmp_chi.val_grid_point(1, k, j, 0) ;       // BC chi	   
	 }
       
       psi_bound.std_base_scal() ;
       nn_bound.std_base_scal() ;
       theta_bound.std_base_scal() ;
       chi_bound.std_base_scal() ;       


       
       //=================================
       // Resolution of elliptic equations
       //=================================
       
       // Resolution of the Poisson equation for Psi
       // ------------------------------------------
       Scalar psi_jp1 (mp) ;
       if (solve_psi == 1) {
	 
	psi_jp1 = source_psi_spher.poisson_neumann(psi_bound, 0) + 1. ;
	
	// Test:
	maxabs(psi_jp1.laplacian() -  source_psi_spher,
	       "Absolute error in the resolution of the equation for Psi") ;  
	// Relaxation (relax=1 -> new ; relax=0 -> old )  
	psi_jp1 = relax * psi_jp1 + (1 - relax) * psi_j ;
       }
       
      // Resolution of the Poisson equation for the lapse
      // ------------------------------------------------
       Scalar nn_jp1 (mp) ;
       if (solve_lapse == 1) {
	 
	 nn_jp1 = source_nn_spher.poisson_dirichlet(nn_bound, 0) + 1. ;
	
	 // Test:
	 maxabs(nn_jp1.laplacian() - source_nn_spher,
		"Absolute error in the resolution of the equation for N") ;
	 
	 // Relaxation (relax=1 -> new ; relax=0 -> old )  
	 if (mer==0)
	  n_evol.update(nn_jp1, jtime, ttime) ; 
	 else
	   nn_jp1 = relax * nn_jp1 + (1 - relax) * nn_j ;
	 
       }


       // Resolution for chi and Theta
       // ----------------------------
       Scalar theta_jp1 (mp) ;
       Scalar chi_jp1 (mp) ;

       if (solve_shift == 1) {
			  
	 theta_jp1 = source_theta_spher.poisson_dirichlet(theta_bound, 0)  ;
	 chi_jp1 = source_chi_spher.poisson_dirichlet(chi_bound, 0)  ;
	 
	 // Test:
	 maxabs(theta_jp1.laplacian() - source_theta_spher,
		"Absolute error in the resolution of the equation for Theta") ;  
	 maxabs(chi_jp1.laplacian() - source_chi_spher,
		"Absolute error in the resolution of the equation for chi") ;
	 
	 // Relaxation (relax=1 -> new ; relax=0 -> old )  
       theta_jp1 = relax * theta_jp1 + (1 - relax) * theta_j ;
       chi_jp1 = relax * chi_jp1 + (1 - relax) * chi_j ;
       
       }
    

	    
       //===========================================
       //      Convergence control
       //===========================================
       
       double diff_nn, diff_psi, diff_btilde, diff_theta, diff_chi ;
       diff_nn = 1.e-16 ;
       diff_psi = 1.e-16 ;
       diff_btilde = 1.e-16 ;
       diff_theta = 1.e-16 ;
       diff_chi = 1.e-16 ;

	if (solve_lapse == 1)
	  diff_nn = max( diffrel(nn_j, nn_jp1) ) ;   
	if (solve_psi == 1)
	  diff_psi = max( diffrel(psi_j, psi_jp1) ) ; 
	if (solve_shift == 1)
	  diff_theta = max( diffrel(theta_jp1, theta_j) ) ; 
	if (solve_shift == 1)
	  diff_chi = max( diffrel(chi_jp1, chi_j) ) ; 

	
	cout << "step = " << mer << " :  diff_psi = " << diff_psi 
	     << "  diff_nn = " << diff_nn 
	     << "  diff_theta = " << diff_theta  
	     << "  diff_chi = " << diff_chi <<  endl ;
	cout << "----------------------------------------------" << endl ;
	if ((diff_psi<precis) && (diff_nn<precis) && (diff_theta<precis) &&  (diff_chi<precis))
	  break ; 
	
	if (mer>0) {conv << mer << "  " << log10(diff_nn) << " " << log10(diff_psi) 
			 << " " << log10(diff_theta) 
		         << " " << log10(diff_chi) << endl ;  } ;
	
	//=============================================
	//      Updates for next step 
	//=============================================
	
	
	if (solve_psi == 1)
	  set_psi(psi_jp1) ; 
	  psi_j = psi_jp1 ; 
	if (solve_lapse == 1)
	  n_evol.update(nn_jp1, jtime, ttime) ; 
	  nn_j = nn_jp1 ;
	if (solve_shift == 1)
	  {
	    theta_j = theta_jp1 ;
	    chi_j = chi_jp1 ;
	    chi_jp1.mult_r() ;
	    Vector beta_jp1 (chi_jp1 * tgam().radial_vect()) ;
	  }

	if (solve_shift == 1 || solve_lapse == 1)
	  {
	    update_aa() ;
	  }

	// Saving ok K_{ij}s^is^j
	// -----------------------
	
	Scalar kkss (contract(k_dd(), 0, 1, gam().radial_vect()*
			      gam().radial_vect(), 0, 1)) ;
	double max_kss = kkss.val_grid_point(1, 0, 0, 0) ;
	double min_kss = kkss.val_grid_point(1, 0, 0, 0) ;
	
	Scalar aaLss (pow(psi(), 6) * kkss) ;
	double max_aaLss = aaLss.val_grid_point(1, 0, 0, 0) ;
	double min_aaLss = aaLss.val_grid_point(1, 0, 0, 0) ;
	
	Scalar hh_tilde (contract(met_gamt.radial_vect().derive_cov(met_gamt), 0, 1)) ;
	double max_hh_tilde = hh_tilde.val_grid_point(1, 0, 0, 0) ;
	double min_hh_tilde = hh_tilde.val_grid_point(1, 0, 0, 0) ;
	
	

	for (int k=0 ; k<nnp ; k++)
	  for (int j=0 ; j<nnt ; j++){
	    if (kkss.val_grid_point(1, k, j, 0) > max_kss)
	      max_kss = kkss.val_grid_point(1, k, j, 0) ;
	    if (kkss.val_grid_point(1, k, j, 0) < min_kss)
	      min_kss = kkss.val_grid_point(1, k, j, 0) ;

	    if (aaLss.val_grid_point(1, k, j, 0) > max_aaLss)
	      max_aaLss = aaLss.val_grid_point(1, k, j, 0) ;
	    if (aaLss.val_grid_point(1, k, j, 0) < min_aaLss)
	      min_aaLss = aaLss.val_grid_point(1, k, j, 0) ;
	    
	    if (hh_tilde.val_grid_point(1, k, j, 0) > max_hh_tilde)
	      max_hh_tilde = hh_tilde.val_grid_point(1, k, j, 0) ;
	    if (hh_tilde.val_grid_point(1, k, j, 0) < min_hh_tilde)
	      min_hh_tilde = hh_tilde.val_grid_point(1, k, j, 0) ;
	    
	  }
	
	
	kss << mer << " " << max_kss << " " << min_kss << " " << max_aaLss << " " << min_aaLss
	    << " " <<  -1 * max_hh_tilde << " "  << -1 * min_hh_tilde << endl ;
    }
    
    conv.close() ;   
    kss.close() ;

} 



