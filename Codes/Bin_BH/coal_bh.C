//standard
#include <stdlib.h>
#include <math.h>
#include <fstream.h>

// LORENE
#include "type_parite.h"
#include "nbr_spx.h"
#include "proto.h"
#include "coord.h"
#include "tenseur.h"
#include "bhole.h"
#include "utilitaires.h"
#include "graphique.h"
#include "unites.h"

void main(int argc, char** argv) {
    
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
    
    int nbre = 10 ;
    double omega = omega_inf ;
    
    for (int i=0 ; i<nbre ; i++) {
	
      courant = bin ;
      double erreur = courant.coal (omega, precis, relax, nbre_homme, 1) ;
      fiche_omega << omega << " " << erreur << endl ;
	 
      char name[20] ;
      sprintf(name, "bin_%e.dat", omega) ;
      FILE* fich_sortie = fopen(name, "w") ;
      grille.sauve(fich_sortie) ;
      map_un.sauve(fich_sortie) ;
      map_deux.sauve(fich_sortie) ;
      courant(1).sauve(fich_sortie) ;
      courant(2).sauve(fich_sortie) ;
      fclose(fich_sortie) ;
	 
      omega += (omega_sup-omega_inf)/(nbre-1) ;
      
    }
    
    
   
    fiche_omega.close() ;
}
