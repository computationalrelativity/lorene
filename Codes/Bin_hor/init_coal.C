/*
 * Main code for computing initial configuration
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

char init_coal_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2006/08/01 14:13:41  f_limousin
 * New version...
 *
 * Revision 1.4  2006/06/29 08:54:52  f_limousin
 * Boundary conditions and grid writen in resformat.dat
 *
 * Revision 1.3  2006/06/28 13:36:52  f_limousin
 * Convergence to a given irreductible mass
 *
 * Revision 1.2  2006/05/24 16:59:08  f_limousin
 * New version
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
    
    double  precis, relax, radius, separation, lim_nn, mass_irr ;
    int nz, nt, np, nr1, nrp1, bound_nn, bound_psi, search_mass ;
    
    param_init.getline(blabla, 120) ;
    param_init.getline(blabla, 120) ;
    param_init >> separation ; param_init.getline(blabla, 120) ;
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
      if (l==1) nr_tab[1] = nr1 ;
      else nr_tab[l] = nrp1 ;
      np_tab[l] = np ; 
      nt_tab[l] = nt ; 
      param_init >> bornes[l] ;

    }
    radius = bornes[1] ;
    param_init.getline(blabla, 120) ;
    bornes[nz] = __infinity ; 

    param_init >> precis ; param_init.getline(blabla, 120) ;
    param_init >> relax ; param_init.getline(blabla, 120) ;
    
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
    
    map_un.set_ori (separation/2.,0, 0) ;
    map_deux.set_ori (-separation/2., 0, 0) ;
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

    ifstream param_coal("par_coal.d") ;
    if ( !param_coal.good() ) {
      cout << "Problem with opening the file par_coal.d ! " << endl ;
      abort() ;
    }
 
    param_coal.ignore(1000, '\n') ;
    param_coal.ignore(1000, '\n') ;
    param_coal.getline(nomini, 80) ;
    param_coal >> omega_init ; param_coal.getline(blabla, 120) ;
    param_coal >> search_mass ;
    param_coal >> mass_irr ;  param_coal.ignore(1000, '\n');
    param_coal >> precis_viriel ;  param_coal.getline(blabla, 120) ;
    param_coal >> relax ; param_coal.getline(blabla, 120) ;
    param_coal >> nb_om ; param_coal.getline(blabla, 120) ;
    param_coal >> nb_it ; param_coal.getline(blabla, 120) ;
    param_coal >> bound_beta ; param_coal.getline(blabla, 120) ;

    param_coal.close() ;

    bin.set_omega(0) ;
    bin.set(1).n_comp (bin(2)) ;
    bin.set(1).psi_comp (bin(2)) ;
    bin.set(1).beta_comp (bin(2)) ;
    bin.set(2).n_comp (bin(1)) ;
    bin.set(2).psi_comp (bin(1)) ;
    bin.set(2).beta_comp (bin(1)) ;
    bin.decouple() ;
    bin.extrinsic_curvature() ;

    cout << "CALCUL AVEC SEPARATION = " << separation << endl ;

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
			     fich_viriel, fich_kss, step, search_mass,
			     mass_irr, 1) ;
    step += nb_om + nb_it ;

    // Convergence to the true Omega
    // ------------------------------

    bool boucle = true ;
    double omega = omega_init ;

    while (boucle) {

      omega = omega * pow((2-erreur)/(2-2*erreur), 1.) ;
  
      Scalar beta_old (bin(1).beta_auto()(1)) ;

      erreur = bin.coal (omega, relax, 1, 0, bound_nn,
                         lim_nn, bound_psi, bound_beta,
                         fich_iteration, fich_correction,
                         fich_viriel, fich_kss, step, search_mass,
			 mass_irr, 1) ;
 
     double erreur_it = 0 ;
      Tbl diff (diffrelmax (beta_old, bin(1).beta_auto()(1))) ;
      for (int i=1 ; i<bin(1).get_mp().get_mg()->get_nzone() ; i++)
        if (diff(i) > erreur_it)
          erreur_it = diff(i) ;


      if (fabs(erreur) < precis_viriel && erreur_it < precis_viriel)
        boucle = false ;

      step += 1 ;

    }

    fich_iteration.close() ;
    fich_correction.close() ;
    fich_viriel.close() ;
    fich_kss.close() ;

    FILE* fich_sortie = fopen("bin.dat", "w") ;
    grid.sauve(fich_sortie) ;
    map_un.sauve(fich_sortie) ;
    map_deux.sauve(fich_sortie) ;
    bin.sauve(fich_sortie, true) ;
    fclose(fich_sortie) ;

    ofstream seqfich("resformat.dat") ;
    if ( !seqfich.good() ) {
      cout << "coal_bh : problem with opening the file resformat.d !"
	   << endl ;
      abort() ;
    }
    bin.write_global(seqfich, lim_nn, bound_nn, bound_psi, bound_beta) ;
    seqfich.close() ;

    delete [] nr_tab ;
    delete [] nt_tab ;
    delete [] np_tab ;
    delete [] type_r ;
    delete [] bornes ;

    return 1 ;
}
