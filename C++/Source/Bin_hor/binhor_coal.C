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


void Bin_hor::set_statiques (double precis, double relax) {
    
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
	
	solve_psi (precis, relax) ;
	solve_lapse (precis, relax) ;
	
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

double Bin_hor::coal (double angu_vel, double precis, double relax, 
			double nb_ome, const int sortie) {
    
    assert (omega == 0) ;
    int nz = hole1.mp.get_mg()->get_nzone() ;
    
    int indic = 1 ;
    int conte = 0 ;
    
    char name_iteration[40] ;
    char name_correction[40] ;
    char name_viriel[40] ;
     
    sprintf(name_iteration, "ite_%e.dat", angu_vel) ;
    sprintf(name_correction, "cor_%e.dat", angu_vel) ;
    sprintf(name_viriel, "vir_%e.dat", angu_vel) ;
    
    ofstream fiche_iteration(name_iteration) ;
    fiche_iteration.precision(8) ; 

    ofstream fiche_correction(name_correction) ;
    fiche_correction.precision(8) ; 
    
    ofstream fiche_viriel(name_viriel) ;
    fiche_viriel.precision(8) ; 

    
    // LA BOUCLE EN AUGMENTANT OMEGA  : 
    cout << "OMEGA AUGMENTE A LA MAIN." << endl ;
    double homme = 0 ;
    for (int pas = 0 ; pas <nb_ome ; pas ++) {
	
	homme += angu_vel/nb_ome ;
	set_omega (homme) ;
	Scalar beta_un_old (hole1.beta_auto()(1)) ;
	
	solve_shift (precis, relax) ;
	extrinsic_curvature() ;
	
	solve_psi (precis, relax) ;
	solve_lapse (precis, relax) ;
	
	double erreur = 0 ;
	Tbl diff (diffrelmax (beta_un_old, hole1.beta_auto()(1))) ;
	for (int i=1 ; i<nz ; i++)
	    if (diff(i) > erreur)
		erreur = diff(i) ;
	if (sortie != 0) {
	    fiche_iteration << conte << " " << erreur << endl ;
	    fiche_correction << conte << " " << hole1.regul << endl ;
	    fiche_viriel << conte << " " << viriel() << endl ;
	    }
	    
	cout << "PAS TOTAL : " << conte << " DIFFERENCE : " << erreur << endl ;
	conte ++ ;
    }
    
    // BOUCLE AVEC OMEGA BLOQUE :
    cout << "OMEGA BLOQUE" << endl ;
    indic = 1 ; 
    double erreur ;
    while (indic == 1) {
	
	Scalar beta_un_old (hole1.beta_auto()(1)) ;
	solve_shift (precis, relax) ;
	extrinsic_curvature() ;
	
	solve_psi (precis, relax) ;
	solve_lapse (precis, relax) ;

	erreur = 0 ;
	Tbl diff (diffrelmax (beta_un_old, hole1.beta_auto()(1))) ;
	for (int i=1 ; i<nz ; i++)
	    if (diff(i) > erreur)
		erreur = diff(i) ;
	
	if (sortie != 0) {
	    fiche_iteration << conte << " " << erreur << endl ;
	    fiche_correction << conte << " " << hole1.regul << endl ;
	    fiche_viriel << conte << " " << viriel() << endl ;
	    }
	    
	cout << "PAS TOTAL : " << conte << " DIFFERENCE : " << erreur << endl ;
	if (erreur < precis)
	    indic = -1 ;
	conte ++ ;
    }
    
    fiche_iteration.close() ;
    fiche_correction.close() ;
    fiche_viriel.close() ;
      
    return viriel() ;
}
