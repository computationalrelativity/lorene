/*
 * Main code for computation of binary black hole quasiequilibrium
 *  configurations
 *
 *
 */

/*
 *   Copyright (c) 2005 Philippe Grandclement
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

char coal_bh_mdiff_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:09:41  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2005/08/29 15:10:19  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * $Header$
 *
 */


 #include "headcpp.h"
 
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

int main(int argc, char** argv) {
    // Lecture du fichier de parametres :
     if (argc <3) {
	cout <<" Passer nom des fichiers en arguments SVP !" << endl ;
	abort() ;
    }

   char* name_fich = argv[2] ;
    
    FILE* fich = fopen(name_fich, "r") ;
    Mg3d grille_un (fich) ;
    Mg3d grille_deux (fich) ;
    Map_af map_un (grille_un, fich) ;
    Map_af map_deux (grille_deux, fich) ;
    Bhole hole_un (map_un, fich) ;
    Bhole hole_deux (map_deux, fich) ;
    fclose(fich) ;
    
    char blabla [120] ;
    name_fich = argv[1] ;
    ifstream param(name_fich) ;
    
    double precis, relax, search_ome, m1, m2 ;
    int nbr_ome ;
    
    param >> m1 ; param >> m2 ; param.getline(blabla, 120) ;
    param >> search_ome ; param.getline(blabla, 120) ;
    param >> precis ; param.getline(blabla, 120) ;
    param >> relax ; param.getline(blabla, 120) ;
    param >> nbr_ome ; param.getline(blabla, 120) ;
    
    param.close() ;
    
    // Le fichier sortie pour la recherche de omega :
   
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
    
    double distance = bin(1).get_mp().get_ori_x() - bin(2).get_mp().get_ori_x() ;  
    cout << "CALCUL AVEC L = " << distance << endl ;
   
    bin.coal (precis, relax, nbr_ome, search_ome, m1, m2, 1) ;

    char name[20] ;
    sprintf(name, "bin.dat") ;
    FILE* fich_sortie = fopen(name, "w") ;
    grille_un.sauve(fich_sortie) ;
    grille_deux.sauve(fich_sortie) ;
    map_un.sauve(fich_sortie) ;
    map_deux.sauve(fich_sortie) ;
    bin(1).sauve(fich_sortie) ;
    bin(2).sauve(fich_sortie) ;
    fclose(fich_sortie) ;
    
    return 0 ;
}
