/*
 * Main code for computation of binary black hole quasiequilibrium
 *  configurations
 *
 *
 */

/*
 *   Copyright (c) 2002 Philippe Grandclement
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
 * Revision 1.6  2014/10/13 08:53:52  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:09:41  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2004/03/25 12:35:35  j_novak
 * now using namespace Unites
 *
 * Revision 1.3  2003/11/14 08:14:52  p_grandclement
 * Correction of the codes for binary black holes in circular orbit
 *
 * Revision 1.2  2003/04/09 15:01:14  e_gourgoulhon
 *
 * Added copyright and log messages.
 *
 *
 * revision 1.1
 * date: 2003/03/31 16:08:23;  author: p_grandclement;  state: Exp;
 * Add codes for the binary black holes computation :
 *       - init_bh.C : computes the Misner-Lindquist solution (used as initial guess for the coal code)
 *       - param_statiques.d : example of parameter file for init_bh
 *       - coal_bh.C : computation of the binary
 *       - param_coal.d : example of the paramater file for coal_bh
 *       - lit_holes_bin.C : reads a configuration.
 * 
 * $Header$
 *
 */


//standard
#include <cstdlib>
#include <cmath>

// LORENE
#include "type_parite.h"
#include "nbr_spx.h"
#include "proto.h"
#include "coord.h"
#include "tenseur.h"
#include "bhole.h"
#include "utilitaires.h"
#include "graphique.h"

using namespace Lorene ;

int main(int argc, char** argv) {
    
    // Lecture du fichier de parametres :
     if (argc <3) {
	cout <<" Passer nom des fichiers en arguments SVP !" << endl ;
	abort() ;
    }

    char* name_fich = argv[2] ;
   
    FILE* fich = fopen(name_fich, "r") ;
    Mg3d grille (fich) ;
    Map_af map_un (grille, fich) ;
    Map_af map_deux (grille, fich) ;
    Bhole hole_un (map_un, fich) ;
    Bhole hole_deux (map_deux, fich) ;
    fclose(fich) ;
    
    char blabla [120] ;
    name_fich = argv[1] ;
    ifstream param(name_fich) ;
    
    double omega_inf, omega_sup, precis, 
	relax, precis_viriel ;
    int nbre_homme ;
    
    param >> omega_inf ; param >> omega_sup ; param.getline(blabla, 120) ;
    param >> precis ; param.getline(blabla, 120) ;
    param >> precis_viriel ;  param.getline(blabla, 120) ;
    param >> relax ; param.getline(blabla, 120) ;
    param >> nbre_homme ; param.getline(blabla, 120) ;
    
    param.close() ;
    
    // Le fichier sortie pour la recherche de omega :
    char name_omega[20] ;
    sprintf(name_omega, "omega.dat") ;
    ofstream fiche_omega(name_omega) ;
    fiche_omega.precision(8) ;

    Bhole_binaire bin (map_un, map_deux) ;
    bin.set(1) = hole_un ;
    bin.set(2) = hole_deux ;
    bin.set_omega(0) ;
    bin.set(1).fait_n_comp (bin(2)) ;
    bin.set(1).fait_psi_comp (bin(2)) ;
    bin.set(2).fait_n_comp (bin(1)) ;
    bin.set(2).fait_psi_comp (bin(1)) ;
    bin.fait_decouple() ;
    bin.fait_tkij() ;
    
    double beta = bin(1).get_mp().get_ori_x() - bin(2).get_mp().get_ori_x() ;
    beta /= bin(1).get_rayon() ;
    
    cout << "CALCUL AVEC BETA = " << beta << endl ;
    
    Bhole_binaire courant(map_un, map_deux) ;
    
    courant = bin ;
    double omega_min = omega_inf ;
    double erreur_min = courant.coal (omega_min, precis, relax, nbre_homme, 1) ; 
    fiche_omega << omega_min << " " << erreur_min << endl ;
    if (erreur_min < 0) {
      cout << "Borne inf. too big" << endl ;
      abort() ;
    }

    courant = bin ;
    double omega_max = omega_sup ;
    double erreur_max = courant.coal (omega_max, precis, relax, nbre_homme, 1) ;
    if (erreur_max > 0) {
      cout << "Borne max. too small" << endl ;
      abort() ;
    }
    fiche_omega << omega_max << " " << erreur_max << endl ;
	 
    bool boucle = true ;
    double erreur, omega ;

    while (boucle) {
      
      courant = bin ;
      omega = omega_min - erreur_min * (omega_max-omega_min)/(erreur_max-erreur_min) ;
      erreur = courant.coal (omega, precis, relax, nbre_homme, 1) ;
      
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
    grille.sauve(fich_sortie) ;
    map_un.sauve(fich_sortie) ;
    map_deux.sauve(fich_sortie) ;
    courant(1).sauve(fich_sortie) ;
    courant(2).sauve(fich_sortie) ;
    fclose(fich_sortie) ;
    
    fiche_omega.close() ;
    return 1 ;
}
