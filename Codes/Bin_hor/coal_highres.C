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

char coal_highres_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2006/05/24 16:59:08  f_limousin
 * New version
 *
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
    ifstream param_init("par_init.d") ;

    double  precis, relax, radius, beta, lim_nn ;
    int nz, nt, np, nr1, nrp1, bound_nn, bound_psi ;

    param_init.getline(blabla, 120) ;
    param_init.getline(blabla, 120) ;
    param_init >> beta ; param_init.getline(blabla, 120) ;
    param_init >> nz ; param_init.getline(blabla, 120) ;
    param_init >> nt; param_init.ignore(1000, '\n');
    param_init >> np; param_init.ignore(1000, '\n');
    param_init >> nr1; param_init.ignore(1000, '\n');
    param_init >> nrp1; param_init.ignore(1000, '\n');

    double* bornes = new double[nz+1] ;
    int* nr_tab = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];

    for (int l=0 ; l<nz ; l++){
      if (l==1) nr_tab[1] = 17 ;
      else nr_tab[l] = 17 ;
      np_tab[l] = 12 ;
      nt_tab[l] = 13 ;
      param_init >> bornes[l] ;

    }
    radius = bornes[1] ;
    param_init.getline(blabla, 120) ;
    bornes[nz] = __infinity ;

    param_init >> precis ; param_init.getline(blabla, 120) ;
    param_init >> relax ; param_init.getline(blabla, 120) ;
    double distance = radius*beta ;

    param_init >> bound_nn ;
    param_init >> lim_nn ;  param_init.ignore(1000, '\n');
    param_init >> bound_psi ;  param_init.ignore(1000, '\n');

    param_init.close() ;

    // Type of sampling in theta and phi :
    int type_t = SYM ;
    int type_p = NONSYM ;


    int* type_r = new int[nz] ;
    type_r[0] = RARE ;
    for (int l=1 ; l<nz-1 ; l++)
      type_r[l] = FIN ;
    type_r[nz-1] = UNSURR ;

    Mg3d grid (nz, nr_tab, type_r, nt_tab, type_t, np_tab, type_p) ;
    Map_af map_un (grid, bornes) ;
    Map_af map_deux (grid, bornes) ;

    map_un.set_ori (distance/2.,0, 0) ;
    map_deux.set_ori (-distance/2., 0, 0) ;
    map_deux.set_rot_phi (M_PI) ;

    int depth = 3 ;
    Bin_hor bin (map_un, map_deux, depth) ;
    bin.set_statiques(precis, relax, bound_nn, lim_nn, bound_psi) ;

    FILE* fich = fopen("static.d", "w") ;
    grid.sauve(fich) ;
    map_un.sauve(fich) ;
    map_deux.sauve(fich) ;
    bin.sauve(fich, true) ;
    fwrite_be(&bound_nn, sizeof(int), 1, fich) ;
    fwrite_be (&lim_nn, sizeof(double), 1, fich) ;
    fwrite_be(&bound_psi, sizeof(int), 1, fich) ;
    fclose(fich) ;

    // Part of coal
    // ------------

    char nomini[120] ;
    double omega_init, precis_viriel;
    int nb_om, nb_it, bound_beta ;
    
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
        
    // Le fichier sortie pour la recherche de omega :
    char name_omega[20] ;
    sprintf(name_omega, "omega.dat") ;
    ofstream fiche_omega(name_omega) ;
    fiche_omega.precision(8) ;
      
    bin.set_omega(0) ;
    bin.set(1).n_comp (bin(2)) ;
    bin.set(1).psi_comp (bin(2)) ;
    bin.set(1).beta_comp (bin(2)) ;
    bin.set(2).n_comp (bin(1)) ;
    bin.set(2).psi_comp (bin(1)) ;
    bin.set(2).beta_comp (bin(1)) ;
    bin.decouple() ;
    bin.extrinsic_curvature() ;
        
    cout << "CALCUL AVEC BETA = " << beta << endl ;
    
    ofstream fich_iteration("iteration.dat") ;
    fich_iteration.precision(8) ; 

    ofstream fich_correction("correction.dat") ;
    fich_correction.precision(8) ; 
    
    ofstream fich_viriel("viriel.dat") ;
    fich_viriel.precision(8) ; 

    ofstream fich_kss("kss.dat") ;
    fich_kss.precision(8) ; 

    fich_iteration << "# step  precision  omega"  << endl ;
    fich_correction << "# step  regularisation  omega"  << endl ;
    fich_viriel << "# step  viriel  omega"  << endl ;
    fich_kss << "# step  kss  omega"  << endl ;

    int step = 0 ;
     
    cout << "step = " << step << endl ;
    double erreur = bin.coal(omega_init, relax, nb_om, nb_it, bound_nn,
			       lim_nn, bound_psi, bound_beta, 
			       fich_iteration, fich_correction,
			       fich_viriel, fich_kss, step, 1) ;
    step += nb_om + nb_it ;
 
    fiche_omega << omega_init << " " << erreur << endl ;
    
    // Convergence to the true Omega
    // ------------------------------

    bool boucle = true ;
    double omega = omega_init ;

    while (boucle) {
      
      omega = omega * pow((2-erreur)/(2-2*erreur), 1.) ;
      erreur = bin.coal (omega, relax, 1, 0, bound_nn,
			 lim_nn, bound_psi, bound_beta, 
			 fich_iteration, fich_correction,
			 fich_viriel, fich_kss, step, 1) ;
         
      fiche_omega << omega << " " << erreur << endl ;
      
      if (fabs(erreur) < precis_viriel)
	boucle = false ;

      step += 1 ;

    }

    char name[40] ;
    sprintf(name, "bin_%i.dat", grid.get_nr(1)) ;

    FILE* fich_sortie = fopen(name, "w") ;
    grid.sauve(fich_sortie) ;
    map_un.sauve(fich_sortie) ;
    map_deux.sauve(fich_sortie) ;
    bin.sauve(fich_sortie, true) ;
    fclose(fich_sortie) ;
    
    sprintf(name, "resformat_%i.dat", grid.get_nr(1)) ;
    ofstream seqfich(name) ; 
    if ( !seqfich.good() ) {
	cout << "coal_bh : problem with opening the file resformat.d !" 
	     << endl ;
	abort() ;
    }
    bin.write_global(seqfich) ; 
    seqfich.close() ; 


    // New mapping (with more collocation points)
    // -------------------------------------------

    // ??? TO BE IMPROVED WHEN Bin_hor::operator= will be possible... ????
    

    // Resolution 25x17x16
    // ------------------------

    for (int l=0 ; l<nz ; l++){
      if (l==1) nr_tab[1] = 25 ;
      else nr_tab[l] = 25 ;
      np_tab[l] = 16 ;
      nt_tab[l] = 17 ;
    }

    Mg3d grid_25 (nz, nr_tab, type_r, nt_tab, type_t, np_tab, type_p) ;
    
    Map_af map_un_25 (grid_25, bornes) ;
    Map_af map_deux_25 (grid_25, bornes) ;
      
    map_un_25.set_ori (distance/2.,0, 0) ;
    map_deux_25.set_ori (-distance/2., 0, 0) ;
    map_deux_25.set_rot_phi (M_PI) ;
    
    Bin_hor bin_25 (map_un_25, map_deux_25, depth) ;
    
    bin_25.import_bh (bin) ;
    
    bin_25.set(1).n_comp (bin_25(2)) ;
    bin_25.set(1).psi_comp (bin_25(2)) ;
    bin_25.set(1).beta_comp (bin_25(2)) ;
    bin_25.set(2).n_comp (bin_25(1)) ;
    bin_25.set(2).psi_comp (bin_25(1)) ;
    bin_25.set(2).beta_comp (bin_25(1)) ;
    bin_25.decouple() ;
    bin_25.extrinsic_curvature() ;
      
      
    fich_iteration << "#------------- New resolution"  << endl ;
    fich_correction << "#-------------- New resolution" << endl ;
    fich_viriel << "#------------------ New resolution"  << endl ;
    
    erreur = bin_25.coal (omega, relax, 1, 0, bound_nn,
			  lim_nn, bound_psi, bound_beta,
			  fich_iteration, fich_correction,
			  fich_viriel, fich_kss, step, 1) ;
    
    boucle = true ;
    
    while (boucle) {
      
      omega = omega * pow((2-erreur)/(2-2*erreur), 1.) ;
      erreur = bin_25.coal (omega, relax, 1, 0, bound_nn,
			    lim_nn, bound_psi, bound_beta,
			    fich_iteration, fich_correction,
			    fich_viriel, fich_kss, step, 1) ;
      
      fiche_omega << omega << " " << erreur << endl ;
      
      if (fabs(erreur) < precis_viriel)
	boucle = false ;
      
      step += 1 ;
    }
    
    sprintf(name, "bin_%i.dat", grid_25.get_nr(1)) ;
    
    FILE* fich_sortie_25 = fopen(name, "w") ;
    grid_25.sauve(fich_sortie_25) ;
    map_un_25.sauve(fich_sortie_25) ;
    map_deux_25.sauve(fich_sortie_25) ;
    bin_25.sauve(fich_sortie_25, true) ;
    fclose(fich_sortie_25) ;
    
    sprintf(name, "resformat_%i.dat", grid_25.get_nr(1)) ;
    ofstream seqfich_25 (name) ;
    if ( !seqfich_25.good() ) {
      cout << "coal_bh : problem with opening the file resformat.d !"
	   << endl ;
      abort() ;
    }
    bin_25.write_global(seqfich_25) ;
    seqfich_25.close() ;
    
    
    

    // Resolution 33x21x20
    // ------------------------

    for (int l=0 ; l<nz ; l++){
      if (l==1) nr_tab[1] = 33 ;
      else nr_tab[l] = 33 ;
      np_tab[l] = 20 ;
      nt_tab[l] = 21 ;
    }

    Mg3d grid_33 (nz, nr_tab, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_af map_un_33 (grid_33, bornes) ;
    Map_af map_deux_33 (grid_33, bornes) ;

    map_un_33.set_ori (distance/2.,0, 0) ;
    map_deux_33.set_ori (-distance/2., 0, 0) ;
    map_deux_33.set_rot_phi (M_PI) ;

    Bin_hor bin_33 (map_un_33, map_deux_33, depth) ;

    bin_33.import_bh (bin_25) ;

    bin_33.set(1).n_comp (bin_33(2)) ;
    bin_33.set(1).psi_comp (bin_33(2)) ;
    bin_33.set(1).beta_comp (bin_33(2)) ;
    bin_33.set(2).n_comp (bin_33(1)) ;
    bin_33.set(2).psi_comp (bin_33(1)) ;
    bin_33.set(2).beta_comp (bin_33(1)) ;
    bin_33.decouple() ;
    bin_33.extrinsic_curvature() ;


    fich_iteration << "#------------- New resolution"  << endl ;
    fich_correction << "#-------------- New resolution" << endl ;
    fich_viriel << "#------------------ New resolution"  << endl ;

    erreur = bin_33.coal (omega, relax, 1, 0, bound_nn,
                          lim_nn, bound_psi, bound_beta,
                          fich_iteration, fich_correction,
                          fich_viriel, fich_kss, step, 1) ;

    boucle = true ;

    while (boucle) {

      omega = omega * pow((2-erreur)/(2-2*erreur), 1.) ;
      erreur = bin_33.coal (omega, relax, 1, 0, bound_nn,
                            lim_nn, bound_psi, bound_beta,
                            fich_iteration, fich_correction,
                            fich_viriel, fich_kss, step, 1) ;

      fiche_omega << omega << " " << erreur << endl ;

      if (fabs(erreur) < precis_viriel)
        boucle = false ;

      step += 1 ;
    }

    sprintf(name, "bin_%i.dat", grid_33.get_nr(1)) ;


    FILE* fich_sortie_33 = fopen(name, "w") ;
    grid_33.sauve(fich_sortie_33) ;
    map_un_33.sauve(fich_sortie_33) ;
    map_deux_33.sauve(fich_sortie_33) ;
    bin_33.sauve(fich_sortie_33, true) ;
    fclose(fich_sortie_33) ;

    sprintf(name, "resformat_%i.dat", grid_33.get_nr(1)) ;
    ofstream seqfich_33 (name) ;
    if ( !seqfich_33.good() ) {
      cout << "coal_bh : problem with opening the file resformat.d !"
           << endl ;
      abort() ;
    }
    bin_33.write_global(seqfich_33) ;
    seqfich_33.close() ;



    // Resolution 41x29x28
    // ------------------------

    for (int l=0 ; l<nz ; l++){
      if (l==1) nr_tab[1] = 41 ;
      else nr_tab[l] = 41 ;
      np_tab[l] = 28 ;
      nt_tab[l] = 29 ;
    }

    Mg3d grid_41 (nz, nr_tab, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_af map_un_41 (grid_41, bornes) ;
    Map_af map_deux_41 (grid_41, bornes) ;

    map_un_41.set_ori (distance/2.,0, 0) ;
    map_deux_41.set_ori (-distance/2., 0, 0) ;
    map_deux_41.set_rot_phi (M_PI) ;

    Bin_hor bin_41 (map_un_41, map_deux_41, depth) ;

    bin_41.import_bh (bin_33) ;

    bin_41.set(1).n_comp (bin_41(2)) ;
    bin_41.set(1).psi_comp (bin_41(2)) ;
    bin_41.set(1).beta_comp (bin_41(2)) ;
    bin_41.set(2).n_comp (bin_41(1)) ;
    bin_41.set(2).psi_comp (bin_41(1)) ;
    bin_41.set(2).beta_comp (bin_41(1)) ;
    bin_41.decouple() ;
    bin_41.extrinsic_curvature() ;


    fich_iteration << "#------------- New resolution"  << endl ;
    fich_correction << "#-------------- New resolution" << endl ;
    fich_viriel << "#------------------ New resolution"  << endl ;

    erreur = bin_41.coal (omega, relax, 1, 0, bound_nn,
                          lim_nn, bound_psi, bound_beta,
                          fich_iteration, fich_correction,
                          fich_viriel, fich_kss, step, 1) ;

    boucle = true ;

    while (boucle) {

      omega = omega * pow((2-erreur)/(2-2*erreur), 1.) ;
      erreur = bin_41.coal (omega, relax, 1, 0, bound_nn,
                            lim_nn, bound_psi, bound_beta,
                            fich_iteration, fich_correction,
                            fich_viriel, fich_kss, step, 1) ;

      fiche_omega << omega << " " << erreur << endl ;

      if (fabs(erreur) < precis_viriel)
        boucle = false ;

      step += 1 ;
    }

    sprintf(name, "bin_%i.dat", grid_41.get_nr(1)) ;
    FILE* fich_sortie_41 = fopen(name, "w") ;
    grid_41.sauve(fich_sortie_41) ;
    map_un_41.sauve(fich_sortie_41) ;
    map_deux_41.sauve(fich_sortie_41) ;
    bin_41.sauve(fich_sortie_41, true) ;
    fclose(fich_sortie_41) ;

    sprintf(name, "resformat_%i.dat", grid_41.get_nr(1)) ;
    ofstream seqfich_41 (name) ;
    if ( !seqfich_41.good() ) {
      cout << "coal_bh : problem with opening the file resformat.d !"
           << endl ;
      abort() ;
    }
    bin_41.write_global(seqfich_41) ;
    seqfich_41.close() ;



    // Resolution 49x33x32
    // ------------------------

    for (int l=0 ; l<nz ; l++){
      if (l==1) nr_tab[1] = 49 ;
      else nr_tab[l] = 49 ;
      np_tab[l] = 32 ;
      nt_tab[l] = 33 ;
    }

    Mg3d grid_49 (nz, nr_tab, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_af map_un_49 (grid_49, bornes) ;
    Map_af map_deux_49 (grid_49, bornes) ;

    map_un_49.set_ori (distance/2.,0, 0) ;
    map_deux_49.set_ori (-distance/2., 0, 0) ;
    map_deux_49.set_rot_phi (M_PI) ;

    Bin_hor bin_49 (map_un_49, map_deux_41, depth) ;

    bin_49.import_bh (bin_41) ;

    bin_49.set(1).n_comp (bin_49(2)) ;
    bin_49.set(1).psi_comp (bin_49(2)) ;
    bin_49.set(1).beta_comp (bin_49(2)) ;
    bin_49.set(2).n_comp (bin_49(1)) ;
    bin_49.set(2).psi_comp (bin_49(1)) ;
    bin_49.set(2).beta_comp (bin_49(1)) ;
    bin_49.decouple() ;
    bin_49.extrinsic_curvature() ;

    fich_iteration << "#------------- New resolution"  << endl ;
    fich_correction << "#-------------- New resolution" << endl ;
    fich_viriel << "#------------------ New resolution"  << endl ;

    erreur = bin_49.coal (omega, relax, 1, 0, bound_nn,
                          lim_nn, bound_psi, bound_beta,
                          fich_iteration, fich_correction,
                          fich_viriel, fich_kss, step, 1) ;

    boucle = true ;

    while (boucle) {

      omega = omega * pow((2-erreur)/(2-2*erreur), 1.) ;
      erreur = bin_49.coal (omega, relax, 1, 0, bound_nn,
                            lim_nn, bound_psi, bound_beta,
                            fich_iteration, fich_correction,
                            fich_viriel, fich_kss, step, 1) ;

      fiche_omega << omega << " " << erreur << endl ;

      if (fabs(erreur) < precis_viriel)
        boucle = false ;

      step += 1 ;
    }

    sprintf(name, "bin_%i.dat", grid_49.get_nr(1)) ;
    FILE* fich_sortie_49 = fopen(name, "w") ;
    grid_49.sauve(fich_sortie_49) ;
    map_un_49.sauve(fich_sortie_49) ;
    map_deux_49.sauve(fich_sortie_49) ;
    bin_49.sauve(fich_sortie_49, true) ;
    fclose(fich_sortie_49) ;

    sprintf(name, "resformat_%i.dat", grid_49.get_nr(1)) ;
    ofstream seqfich_49 (name) ;
    if ( !seqfich_49.good() ) {
      cout << "coal_bh : problem with opening the file resformat.d !"
           << endl ;
      abort() ;
    }
    bin_49.write_global(seqfich_49) ;
    seqfich_49.close() ;


    // End of the computation

    fich_iteration.close() ;
    fich_correction.close() ;
    fich_viriel.close() ;
    fich_kss.close() ;
    fiche_omega.close() ;

    delete [] nr_tab ;
    delete [] nt_tab ;
    delete [] np_tab ;
    delete [] type_r ;
    delete [] bornes ;

    return EXIT_SUCCESS ;
}
