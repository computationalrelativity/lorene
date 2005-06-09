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
 * Revision 1.6  2005/06/09 16:17:21  f_limousin
 * Many different changes.
 *
 * Revision 1.4  2005/03/10 16:57:02  f_limousin
 * Improve the convergence of the code coal_bh.
 *
 * Revision 1.3  2005/03/04 09:41:34  f_limousin
 * New construction of the object Bin_hor
 *
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
    double omega_init, relax, precis_viriel, lim_nn ;
    int nb_om, nb_it, bound_nn, bound_psi, bound_beta ;
    
    ifstream param("par_coal.d") ;
	if ( !param.good() ) {
		cout << "Problem with opening the file par_coal.d ! " << endl ;
		abort() ;
	}
    param.ignore(1000, '\n') ;
    param.ignore(1000, '\n') ;
    param.getline(nomini, 80) ; 
    param >> omega_init ; param.getline(blabla, 120) ;
    param >> precis_viriel ;  param.getline(blabla, 120) ;
    param >> relax ; param.getline(blabla, 120) ;
    param >> nb_om ; param.getline(blabla, 120) ;
    param >> nb_it ; param.getline(blabla, 120) ;
    param >> bound_beta ; param.getline(blabla, 120) ;    

    param.close() ;
    
    int depth = 3 ;

    FILE* fich = fopen(nomini, "r") ;
    Mg3d grid (fich) ;
    Map_af map_un (grid, fich) ;
    Map_af map_deux (grid, fich) ;
    Bin_hor bin(map_un, map_deux, fich, true, depth) ;
    fread_be(&bound_nn, sizeof(int), 1, fich) ;	
    fread_be(&lim_nn, sizeof(double), 1, fich) ;
    fread_be(&bound_psi, sizeof(int), 1, fich) ;	
    fclose(fich) ;
    
    // Le fichier sortie pour la recherche de omega :
    char name_omega[20] ;
    sprintf(name_omega, "omega.dat") ;
    ofstream fiche_omega(name_omega) ;
    fiche_omega.precision(8) ;
      
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
    
    ofstream fich_iteration("iteration.dat") ;
    fich_iteration.precision(8) ; 

    ofstream fich_correction("correction.dat") ;
    fich_correction.precision(8) ; 
    
    ofstream fich_viriel("viriel.dat") ;
    fich_viriel.precision(8) ; 

    fich_iteration << "# step  precision  omega"  << endl ;
    fich_correction << "# step  regularisation  omega"  << endl ;
    fich_viriel << "# step  viriel  omega"  << endl ;

    int step = 0 ;
    double omega_jp1, erreur_jp1 ;
    double omega_j = omega_init ;
    
    cout << "step = " << step << endl ;
    double erreur_j = bin.coal(omega_j, relax, nb_om, nb_it, bound_nn,
			       lim_nn, bound_psi, bound_beta, 
			       fich_iteration, fich_correction,
			       fich_viriel, step, 1) ;
    step += nb_om + nb_it ;
 
    fiche_omega << omega_j << " " << erreur_j << endl ;
    
    //    nb_om = (nb_om + nb_om%2) / 2 ;
    if (erreur_j < 0) {
      omega_jp1 = 0.8 * omega_j ;
      erreur_jp1 = bin.coal(omega_jp1, relax, nb_om, nb_it, bound_nn,
			    lim_nn, bound_psi, bound_beta, 
			    fich_iteration, fich_correction,
			    fich_viriel, step, 1) ;
      fiche_omega << omega_jp1 << " " << erreur_jp1 << endl ;
    }
    else {
      omega_jp1 = 1.25 * omega_j ;
      erreur_jp1 = bin.coal(omega_jp1, relax, nb_om, nb_it, bound_nn,
			    lim_nn, bound_psi, bound_beta, 
			    fich_iteration, fich_correction,
			    fich_viriel, step, 1) ;
      fiche_omega << omega_jp1 << " " << erreur_jp1 << endl ;
    }
    step += nb_om + nb_it ;

    bool boucle = true ;
    double erreur, omega ;

    while (boucle) {
      
      omega = omega_j - erreur_j * (omega_jp1-omega_j)
	  /(erreur_jp1-erreur_j) ;
      erreur = bin.coal (omega, relax, nb_om, nb_it, bound_nn,
			 lim_nn, bound_psi, bound_beta, 
			 fich_iteration, fich_correction,
			 fich_viriel, step, 1) ;
     step += nb_om + nb_it ;
     
      fiche_omega << omega << " " << erreur << endl ;
      
      if (fabs(erreur) < precis_viriel)
	boucle = false ;

      omega_j = omega_jp1 ;
      erreur_j = erreur_jp1 ;
      omega_jp1 = omega ;
      erreur_jp1 = erreur ;
    }


    fich_iteration.close() ;
    fich_correction.close() ;
    fich_viriel.close() ;

    FILE* fich_sortie = fopen("bin.dat", "w") ;
    grid.sauve(fich_sortie) ;
    map_un.sauve(fich_sortie) ;
    map_deux.sauve(fich_sortie) ;
    bin.sauve(fich_sortie, true) ;
    fclose(fich_sortie) ;
    
    fiche_omega.close() ;

    ofstream seqfich("resformat.dat") ; 
    if ( !seqfich.good() ) {
	cout << "coal_bh : problem with opening the file resformat.d !" 
	     << endl ;
	abort() ;
    }
    bin.write_global(seqfich) ; 
    seqfich.close() ; 

    return EXIT_SUCCESS ;
}
