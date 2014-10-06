/*
 *   Copyright (c) 2004-2005 Francois Limousin
 *                           Jose-Luis Jaramillo
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


char binhor_coal_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.14  2014/10/06 15:13:01  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.13  2007/04/13 15:28:55  f_limousin
 * Lots of improvements, generalisation to an arbitrary state of
 * rotation, implementation of the spatial metric given by Samaya.
 *
 * Revision 1.12  2006/08/01 14:37:19  f_limousin
 * New version
 *
 * Revision 1.11  2006/06/28 13:36:09  f_limousin
 * Convergence to a given irreductible mass
 *
 * Revision 1.10  2006/05/24 16:56:37  f_limousin
 * Many small modifs.
 *
 * Revision 1.9  2005/09/13 18:33:15  f_limousin
 * New function vv_bound_cart_bin(double) for computing binaries with
 * berlin condition for the shift vector.
 * Suppress all the symy and asymy in the importations.
 *
 * Revision 1.8  2005/07/11 08:21:57  f_limousin
 * Implementation of a new boundary condition for the lapse in the binary
 * case : boundary_nn_Dir_lapl().
 *
 * Revision 1.7  2005/03/10 17:21:52  f_limousin
 * Add the Berlin boundary condition for the shift.
 * Some changes to avoid warnings.
 *
 * Revision 1.6  2005/03/10 17:09:05  f_limousin
 * Display the logarithm of viriel and convergence.
 *
 * Revision 1.5  2005/03/10 16:57:00  f_limousin
 * Improve the convergence of the code coal_bh.
 *
 * Revision 1.4  2005/02/24 17:24:26  f_limousin
 * The boundary conditions for psi, N and beta are now parameters in
 * par_init.d and par_coal.d.
 *
 * Revision 1.3  2005/02/07 10:43:36  f_limousin
 * Add the printing of the regularisation of the shift in the case N=0
 * on the horizon.
 *
 * Revision 1.2  2004/12/31 15:40:21  f_limousin
 * Improve the initialisation of several quantities in set_statiques().
 *
 * Revision 1.1  2004/12/29 16:11:19  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

//standard
#include <cstdlib>

// Lorene
#include "tensor.h"
#include "isol_hor.h"
#include "graphique.h"


void Bin_hor::set_statiques (double precis, double relax, int bound_nn,
			     double lim_nn, int bound_psi) {
    
  int nz = hole1.mp.get_mg()->get_nzone() ;
    
  set_omega(0) ;
  hole1.init_met_trK() ;
  hole2.init_met_trK() ;
  init_bin_hor() ;
  extrinsic_curvature() ;
      
  int indic = 1 ;
  int conte = 0 ;
 
  cout << "Static black holes : " << endl ;
  while (indic == 1) {
    Scalar lapse_un_old (hole1.n_auto) ;

    solve_psi (precis, relax, bound_psi) ;
    solve_lapse (precis, relax, bound_nn, lim_nn) ;

    //	des_profile(hole1.nn(), 0, 20, M_PI/2, M_PI) ;

    double erreur = 0 ;
    Tbl diff (diffrelmax (lapse_un_old, hole1.n_auto)) ;
    for (int i=1 ; i<nz ; i++)
      if (diff(i) > erreur)
	erreur = diff(i) ;
	
    cout << "Step : " << conte << " Difference : " << erreur << endl ;
	
    if (erreur < precis)
      indic = -1 ;
    conte ++ ;
  }
}

double Bin_hor::coal (double angu_vel, double relax, int nb_ome,
		      int nb_it, int bound_nn, double lim_nn, 
		      int bound_psi, int bound_beta, double omega_eff,
		      double alpha,
		      ostream& fich_iteration, ostream& fich_correction,
		      ostream& fich_viriel, ostream& fich_kss, 
		      int step, int search_mass, double mass_irr, 
		      const int sortie) {
    
  int nz = hole1.mp.get_mg()->get_nzone() ;

  double precis = 1e-7 ;
    
  // LOOP INCREASING OMEGA  : 
  cout << "OMEGA INCREASED STEP BY STEP." << endl ;
  double homme = get_omega() ;
  double inc_homme = (angu_vel - homme)/nb_ome ;
  for (int pas = 0 ; pas <nb_ome ; pas ++) {
      
    bool verif = false ;
    if (omega_eff == alpha*homme ) verif = true ;
      
    homme += inc_homme ;
    set_omega (homme) ;
    if (verif)
      omega_eff = alpha*homme ;
    Scalar beta_un_old (hole1.beta_auto(1)) ;
      
    solve_shift (precis, relax, bound_beta, omega_eff) ;
    extrinsic_curvature() ;

    solve_psi (precis, relax, bound_psi) ;
    solve_lapse (precis, relax, bound_nn, lim_nn) ;
	
    // Convergence to the given irreductible mass 
    if (search_mass == 1 && step >= 30) {
      double mass_area = sqrt(hole1.area_hor()/16/M_PI) + 
	sqrt(hole2.area_hor()/16/M_PI) ;
      double error_m = (mass_irr - mass_area) / mass_irr ;
      double scaling_r = pow((2-error_m)/(2-2*error_m), 1.) ;
      hole1.mp.homothetie_interne(scaling_r) ;
      hole1.radius = hole1.radius *scaling_r ;
      hole2.mp.homothetie_interne(scaling_r) ;
      hole2.radius = hole2.radius *scaling_r ;
	
      // Update of the different metrics (another possibility would 
      // be to set all derived quantities to 0x0, especially
      // the connection p_connect
      hole1.ff = hole1.mp.flat_met_spher() ;
      hole1.tgam = hole1.mp.flat_met_spher() ;
      hole2.ff = hole2.mp.flat_met_spher() ;
      hole2.tgam = hole1.mp.flat_met_spher() ;
	
    }
      
    cout << "Angular momentum computed at the horizon : " << ang_mom_hor()
	 << endl ;
      
    double erreur = 0 ;
    Tbl diff (diffrelmax (beta_un_old, hole1.beta_auto(1))) ;
    for (int i=1 ; i<nz ; i++)
      if (diff(i) > erreur)
	erreur = diff(i) ;
      
    // Saving ok K_{ij}s^is^j
    // -----------------------
	
    Scalar kkss (contract(hole1.get_k_dd(), 0, 1, 
			  hole1.get_gam().radial_vect()*
			  hole1.get_gam().radial_vect(), 0, 1)) ;
    double max_kss = kkss.val_grid_point(1, 0, 0, 0) ;
    double min_kss = kkss.val_grid_point(1, 0, 0, 0) ;
    int nnp = hole1.mp.get_mg()->get_np(1) ;
    int nnt = hole2.mp.get_mg()->get_nt(1) ;
    for (int k=0 ; k<nnp ; k++)
      for (int j=0 ; j<nnt ; j++){
	if (kkss.val_grid_point(1, k, j, 0) > max_kss)
	  max_kss = kkss.val_grid_point(1, k, j, 0) ;
	if (kkss.val_grid_point(1, k, j, 0) < min_kss)
	  min_kss = kkss.val_grid_point(1, k, j, 0) ;
      }

    if (sortie != 0) {
      fich_iteration << step << " " << log10(erreur) << " " << homme << endl ;
      fich_correction << step << " " << log10(hole1.regul) << " " << homme << endl ;
      //	  fich_viriel << step << " " << log10(fabs(viriel())) << " " << homme << endl ;
      fich_viriel << step << " " << viriel() << " " << homme << " " << hole1.omega_hor() - alpha*homme << " " << omega_eff << endl ;
      fich_kss << step << " " << max_kss << " " << min_kss << endl ;
    }
	    
    cout << "STEP : " << step << " DIFFERENCE : " << erreur << endl ;
    step ++ ;
  }
    
  // LOOP WITH FIXED OMEGA :

  if (nb_it !=0)
    cout << "OMEGA FIXED" << endl ;
  double erreur ;

  for (int pas = 0 ; pas <nb_it ; pas ++) {
	
    Scalar beta_un_old (hole1.beta_auto(1)) ;

    solve_shift (precis, relax, bound_beta, omega_eff) ;
    extrinsic_curvature() ;

    solve_psi (precis, relax, bound_psi) ;
    solve_lapse (precis, relax, bound_nn, lim_nn) ;

    // Convergence to the given irreductible mass 
    if (search_mass == 1 && step >= 30) {
      double mass_area = sqrt(hole1.area_hor()/16/M_PI) + 
	sqrt(hole2.area_hor()/16/M_PI) ;
      double error_m = (mass_irr - mass_area) / mass_irr ;
      double scaling_r = pow((2-error_m)/(2-2*error_m), 1.) ;
	  
      hole1.mp.homothetie_interne(scaling_r) ;
      hole1.radius = hole1.radius *scaling_r ;
      hole2.mp.homothetie_interne(scaling_r) ;
      hole2.radius = hole2.radius *scaling_r ;
    }

    erreur = 0 ;
    Tbl diff (diffrelmax (beta_un_old, hole1.beta_auto(1))) ;
    for (int i=1 ; i<nz ; i++)
      if (diff(i) > erreur)
	erreur = diff(i) ;

    // Saving ok K_{ij}s^is^j
    // -----------------------
	
    Scalar kkss (contract(hole1.get_k_dd(), 0, 1, 
			  hole1.get_gam().radial_vect()*
			  hole1.get_gam().radial_vect(), 0, 1)) ;
    double max_kss = kkss.val_grid_point(1, 0, 0, 0) ;
    double min_kss = kkss.val_grid_point(1, 0, 0, 0) ;
    int nnp = hole1.mp.get_mg()->get_np(1) ;
    int nnt = hole2.mp.get_mg()->get_nt(1) ;
    for (int k=0 ; k<nnp ; k++)
      for (int j=0 ; j<nnt ; j++){
	if (kkss.val_grid_point(1, k, j, 0) > max_kss)
	  max_kss = kkss.val_grid_point(1, k, j, 0) ;
	if (kkss.val_grid_point(1, k, j, 0) < min_kss)
	  min_kss = kkss.val_grid_point(1, k, j, 0) ;
      }

	
    if (sortie != 0) {
      fich_iteration << step << " " << log10(erreur) << " " << homme << endl ;
      fich_correction << step << " " << log10(hole1.regul) << " " << homme << endl ;
      //	  fich_viriel << step << " " << log10(fabs(viriel())) << " " << homme << endl ;
      fich_viriel << step << " " << viriel() << " " << homme << " " << hole1.omega_hor() - alpha*homme << " " << omega_eff << endl ;
      fich_kss << step << " " << max_kss << " " << min_kss << endl ;
    }

    cout << "STEP : " << step << " DIFFERENCE : " << erreur << endl ;
    step ++ ;
  }

  if (nb_it != 0){
    fich_iteration << "#----------------------------"  << endl ;
    fich_correction << "#-----------------------------" << endl ;
    fich_viriel << "#------------------------------"  << endl ;
  }

  return viriel() ;
}
