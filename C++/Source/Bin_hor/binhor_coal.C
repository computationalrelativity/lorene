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
#include <stdlib.h>

// Lorene
#include "tensor.h"
#include "isol_hor.h"


void Bin_hor::set_statiques (double precis, double relax, int bound_nn,
			     double lim_nn, int bound_psi) {
    
    int nz = hole1.mp.get_mg()->get_nzone() ;
    
    set_omega(0) ;
    init_bin_hor() ;

    hole1.init_met_trK() ;
    hole2.init_met_trK() ;
    extrinsic_curvature() ;
      
    int indic = 1 ;
    int conte = 0 ;
 
    cout << "Static black holes : " << endl ;
    while (indic == 1) {
	Scalar lapse_un_old (hole1.n_auto()) ;
	
	solve_psi (precis, relax, bound_psi) ;
	solve_lapse (precis, relax, bound_nn, lim_nn) ;
	
	double erreur = 0 ;
	Tbl diff (diffrelmax (lapse_un_old, hole1.n_auto())) ;
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
		      int bound_psi, int bound_beta,
		      ostream& fich_iteration, ostream& fich_correction,
		      ostream& fich_viriel, int step, const int sortie) {
    
  int nz = hole1.mp.get_mg()->get_nzone() ;
    
  double precis = 1e-7 ;
    
    // LOOP INCREASING OMEGA  : 
    cout << "OMEGA INCREASED STEP BY STEP." << endl ;
    double homme = get_omega() ;
    double inc_homme = (angu_vel - homme)/nb_ome ;
    for (int pas = 0 ; pas <nb_ome ; pas ++) {
	
	homme += inc_homme ;
	set_omega (homme) ;
	Scalar beta_un_old (hole1.beta_auto()(1)) ;
	
	solve_shift (precis, relax, bound_beta) ;
	extrinsic_curvature() ;
	
	solve_psi (precis, relax, bound_psi) ;
	solve_lapse (precis, relax, bound_nn, lim_nn) ;
	
	double erreur = 0 ;
	Tbl diff (diffrelmax (beta_un_old, hole1.beta_auto()(1))) ;
	for (int i=1 ; i<nz ; i++)
	    if (diff(i) > erreur)
		erreur = diff(i) ;
	if (sortie != 0) {
	  fich_iteration << step << " " << erreur << " " << homme << endl ;
	  fich_correction << step << " " << hole1.regul << " " << homme << endl ;
	  fich_viriel << step << " " << viriel() << " " << homme << endl ;
	    }
	    
	cout << "STEP : " << step << " DIFFERENCE : " << erreur << endl ;
	step ++ ;
    }
    
    // LOOP WITH FIXED OMEGA :
    cout << "OMEGA FIXED" << endl ;
    double erreur ;

    for (int pas = 0 ; pas <nb_it ; pas ++) {
	
	Scalar beta_un_old (hole1.beta_auto()(1)) ;
	solve_shift (precis, relax, bound_beta) ;
	extrinsic_curvature() ;
	
	solve_psi (precis, relax, bound_psi) ;
	solve_lapse (precis, relax, bound_nn, lim_nn) ;

	erreur = 0 ;
	Tbl diff (diffrelmax (beta_un_old, hole1.beta_auto()(1))) ;
	for (int i=1 ; i<nz ; i++)
	    if (diff(i) > erreur)
		erreur = diff(i) ;
	
	if (sortie != 0) {
	  fich_iteration << step << " " << erreur << " " << homme << endl ;
	  fich_correction << step << " " << hole1.regul << " " << homme << endl ;
	  fich_viriel << step << " " << viriel() << " " << homme << endl ;
	    }

	cout << "STEP : " << step << " DIFFERENCE : " << erreur << endl ;
 	step ++ ;
   }

    fich_iteration << "#----------------------------"  << endl ;
    fich_correction << "#-----------------------------" << endl ;
    fich_viriel << "#------------------------------"  << endl ;

    return viriel() ;
}
