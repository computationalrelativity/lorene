/*
 * Main code for computation of binary black hole quasiequilibrium
 *  configuration
 *
 *
 */

/*
 *   Copyright (c) 2004-2005 Francois Limousin
 *                           Jose Luis Jaramillo
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

char coal_bh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2005/02/25 12:32:35  f_limousin
 * The boundary conditions for psi, N and beta are now parameters in
 * par_init.d and par_coal.d.
 *
 * Revision 1.1  2004/12/31 15:42:50  f_limousin
 * First version
 *
 * 
 * $Header$
 *
 */


//standard
#include <stdlib.h>
#include <math.h>
//#include <fstream.h>

// LORENE
#include "type_parite.h"
#include "nbr_spx.h"
#include "proto.h"
#include "coord.h"
#include "tenseur.h"
#include "tensor.h"
#include "isol_hor.h"
#include "utilitaires.h"
#include "graphique.h"

int main() {
        
    char blabla [120] ;
    char nomini[120] ;
    double omega_inf, omega_sup, precis, relax, precis_viriel, lim_nn ;
    int nb_om, bound_nn, bound_psi, bound_beta ;
    
    ifstream param("par_coal.d") ;
	if ( !param.good() ) {
		cout << "Problem with opening the file par_coal.d ! " << endl ;
		abort() ;
	}
    param.ignore(1000, '\n') ;
    param.ignore(1000, '\n') ;
    param.getline(nomini, 80) ; 
    param >> omega_inf ; param >> omega_sup ; param.getline(blabla, 120) ;
    param >> precis ; param.getline(blabla, 120) ;
    param >> precis_viriel ;  param.getline(blabla, 120) ;
    param >> relax ; param.getline(blabla, 120) ;
    param >> nb_om ; param.getline(blabla, 120) ;
    param >> bound_beta ; param.getline(blabla, 120) ;    

    param.close() ;
    
    FILE* fich = fopen(nomini, "r") ;
    Mg3d grid (fich) ;
    Map_af map_un (grid, fich) ;
    Map_af map_deux (grid, fich) ;
    Isol_hor hole_un (map_un, fich, true) ;
    Isol_hor hole_deux (map_deux, fich, true) ;
    fread_be(&bound_nn, sizeof(int), 1, fich) ;	
    fread_be(&lim_nn, sizeof(double), 1, fich) ;
    fread_be(&bound_psi, sizeof(int), 1, fich) ;	
    fclose(fich) ;

    // Le fichier sortie pour la recherche de omega :
    char name_omega[20] ;
    sprintf(name_omega, "omega.dat") ;
    ofstream fiche_omega(name_omega) ;
    fiche_omega.precision(8) ;

    
    int depth = 3 ;
    Bin_hor bin (map_un, map_deux, depth) ;
    bin.set(1) = hole_un ;
    bin.set(2) = hole_deux ;
    bin.set_omega(0) ;
    bin.set(1).n_comp (bin(2)) ;
    bin.set(1).psi_comp (bin(2)) ;
    bin.set(2).n_comp (bin(1)) ;
    bin.set(2).psi_comp (bin(1)) ;
    bin.decouple() ;
    bin.extrinsic_curvature() ;
    
    double beta = bin(1).get_mp().get_ori_x() - bin(2).get_mp().get_ori_x() ;
    beta /= bin(1).get_radius() ;
    
    cout << "CALCUL AVEC BETA = " << beta << endl ;
    
    Bin_hor courant(map_un, map_deux, depth) ;
    
    courant = bin ;
    double omega_min = omega_inf ;
    double erreur_min = courant.coal(omega_min, precis, relax, nb_om, bound_nn,
				     lim_nn, bound_psi, bound_beta, 1) ;
 
    fiche_omega << omega_min << " " << erreur_min << endl ;
    if (erreur_min < 0) {
      cout << "Borne inf. too big" << endl ;
      abort() ;
    }

    courant = bin ;
    double omega_max = omega_sup ;
    double erreur_max = courant.coal(omega_max, precis, relax, nb_om, bound_nn,
				     lim_nn, bound_psi, bound_beta, 1) ;
    if (erreur_max > 0) {
      cout << "Borne max. too small" << endl ;
      abort() ;
    }
    fiche_omega << omega_max << " " << erreur_max << endl ;
	 
    bool boucle = true ;
    double erreur, omega ;

    while (boucle) {
      
      courant = bin ;
      omega = omega_min - erreur_min * (omega_max-omega_min)
	  /(erreur_max-erreur_min) ;
      erreur = courant.coal (omega, precis, relax, nb_om, bound_nn,
			     lim_nn, bound_psi, bound_beta,1) ;
      
      fiche_omega << omega << " " << erreur << endl ;
      
      if (fabs(erreur) < precis_viriel)
	boucle = false ;

      
      if (erreur > 0) {
	omega_min = omega ;
	erreur_min = erreur ;
      }
      else {
	omega_max = omega ;
	erreur_max = erreur ;
      }
    }

    char name[20] ;
    sprintf(name, "bin_%e.dat", omega) ;
    FILE* fich_sortie = fopen(name, "w") ;
    grid.sauve(fich_sortie) ;
    map_un.sauve(fich_sortie) ;
    map_deux.sauve(fich_sortie) ;
    courant(1).sauve(fich_sortie, true) ;
    courant(2).sauve(fich_sortie, true) ;
    fclose(fich_sortie) ;
    
    fiche_omega.close() ;
    return 1 ;
}
