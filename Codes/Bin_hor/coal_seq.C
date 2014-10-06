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

char coal_seq_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2014/10/06 15:09:41  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2007/04/13 15:30:58  f_limousin
 * Lots of improvements, generalisation to an arbitrary state of
 * rotation, implementation of the spatial metric given by Samaya.
 *
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
#include <cstdlib>
#include <cmath>

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
    ifstream param_seq("par_seq.d") ;
    
    double  precis, relax, radius, separ, lim_nn ;
    int nz, nt, np, nr1, nrp1, bound_nn, bound_psi ;
    double omega_init, precis_viriel, separ_min, separ_max, mass_irr ;
    int nb_om, nb_it, bound_beta, nb_conf, search_mass ;
    
    param_seq.getline(blabla, 120) ;
    param_seq.getline(blabla, 120) ;
    param_seq >> nt; param_seq.ignore(1000, '\n');
    param_seq >> np; param_seq.ignore(1000, '\n');
    param_seq >> nr1; param_seq.ignore(1000, '\n');
    param_seq >> nrp1; param_seq.ignore(1000, '\n');
    param_seq >> precis ; param_seq.getline(blabla, 120) ;
    param_seq >> relax ; param_seq.getline(blabla, 120) ;
    param_seq >> bound_nn ;
    param_seq >> lim_nn ;  param_seq.ignore(1000, '\n');
    param_seq >> bound_psi ;  param_seq.ignore(1000, '\n');
    param_seq >> omega_init ; param_seq.getline(blabla, 120) ;
    param_seq >> search_mass ;
    param_seq >> mass_irr ;  param_seq.ignore(1000, '\n');
    param_seq >> precis_viriel ;  param_seq.getline(blabla, 120) ;
    param_seq >> nb_om ; param_seq.getline(blabla, 120) ;
    param_seq >> nb_it ; param_seq.getline(blabla, 120) ;
    param_seq >> bound_beta ; param_seq.getline(blabla, 120) ;
    param_seq >> nb_conf ; param_seq.getline(blabla, 120) ;
    param_seq >> separ_min ; param_seq.getline(blabla, 120) ;
    param_seq >> separ_max ; param_seq.getline(blabla, 120) ;
    
    if (nb_conf <=1) {
      cout << "too small value of nb_conf" << endl ;
      abort() ;
    }

    param_seq.close() ;

    for (int conf=0; conf<nb_conf; conf++){
      separ = separ_min + conf * (separ_max-separ_min)/(nb_conf-1) ;

      // To improve ???
      if (separ < 16) omega_init = 0.03 ;
      else omega_init = 0.015 ;
      

      if (separ <= 7){ 
	cout << "too small value of separation" << endl ;
	abort() ;
      }
      
      if (separ > 30){
        cout << "too large value of separation" << endl ;
        abort() ;
      }
      
      if (separ <= 16) nz = 6 ;
      else nz = 7 ;
      
      
      double* bornes = new double[nz+1] ;
      int* nr_tab = new int[nz];
      int* nt_tab = new int[nz];
      int* np_tab = new int[nz];
      
      for (int l=0 ; l<nz ; l++){
	if (l==1) nr_tab[1] = nr1 ;
	else nr_tab[l] = nrp1 ;
	np_tab[l] = np ; 
	nt_tab[l] = nt ; 
      }
      bornes[0] = 0 ;
      bornes[1] = 1 ;
      bornes[2] = 2 ;
      bornes[3] = 4 ;
      if (separ >=10) bornes[4] = 8 ;
      else{ 
	bornes[4] = 6 ;
	bornes[5] = 12 ;
      }
      if (separ >=10 && separ <= 14) bornes[5] = 16 ;
      if (separ > 14 && separ <= 16) bornes[5] = 18 ;
      if (separ > 16 && separ <= 18){
	bornes[5] = 14 ;
	bornes[6] = 28 ;
      }
      if (separ > 18){
	bornes[5] = 16 ;
	bornes[6] = 32 ;
      }
      bornes[nz] = __infinity ; 
      
      radius = 1. ;
 
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
      
      map_un.set_ori (separ/2.,0, 0) ;
      map_deux.set_ori (-separ/2., 0, 0) ;
      map_deux.set_rot_phi (M_PI) ;

      cout << map_un << endl ;
      
      Bin_hor bin (map_un, map_deux) ;
      bin.set_statiques(precis, relax, bound_nn, lim_nn, bound_psi) ;
      
      char name[40] ;
      sprintf(name, "static_%e.d", separ) ;
 
      FILE* fich = fopen(name, "w") ;
      grid.sauve(fich) ;
      map_un.sauve(fich) ;
      map_deux.sauve(fich) ;
      bin.sauve(fich) ;
      fwrite_be(&bound_nn, sizeof(int), 1, fich) ;
      fwrite_be (&lim_nn, sizeof(double), 1, fich) ;
      fwrite_be(&bound_psi, sizeof(int), 1, fich) ;
      fclose(fich) ;
      
      // Drawings
      const Coord& r = bin(1).get_mp().r ;        // r field 
      Mtbl usr = 1 / r ;
      Scalar unsr(bin(1).get_mp()) ;
      unsr = usr ;
      
      Scalar temp = 1. + unsr ;
      temp.std_spectral_base() ;
      
      // Part of coal
      // ------------
      
      bin.set_omega(0) ;
      bin.set(1).n_comp_import (bin(2)) ;
      bin.set(1).psi_comp_import (bin(2)) ;
      bin.set(1).beta_comp_import (bin(2)) ;
      bin.set(2).n_comp_import (bin(1)) ;
      bin.set(2).psi_comp_import (bin(1)) ;
      bin.set(2).beta_comp_import (bin(1)) ;
      bin.decouple() ;
      bin.extrinsic_curvature() ;
      
      cout << "CALCUL AVEC SEPARATION = " << separ << endl ;
      
      sprintf(name, "iteration_%e.dat", separ) ;
      ofstream fich_iteration(name) ;
      fich_iteration.precision(8) ;

      sprintf(name, "correction_%e.dat", separ) ;
      ofstream fich_correction(name) ;
      fich_correction.precision(8) ;
      
      sprintf(name, "viriel_%e.dat", separ) ;
      ofstream fich_viriel(name) ;
      fich_viriel.precision(8) ;
      

      sprintf(name, "kss_%e.dat", separ) ;
      ofstream fich_kss(name) ;
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

	Scalar beta_old (bin(1).get_beta_auto()(1)) ;

	erreur = bin.coal (omega, relax, 1, 0, bound_nn,
			   lim_nn, bound_psi, bound_beta,
			   fich_iteration, fich_correction,
			   fich_viriel, fich_kss, step, search_mass,
			   mass_irr, 1) ;
			   
	double erreur_it = 0 ;
	Tbl diff (diffrelmax (beta_old, bin(1).get_beta_auto()(1))) ;
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
      
      sprintf(name, "bin_%e.dat", separ) ;
      FILE* fich_sortie = fopen(name, "w") ;
      grid.sauve(fich_sortie) ;
      map_un.sauve(fich_sortie) ;
      map_deux.sauve(fich_sortie) ;
      bin.sauve(fich_sortie) ;
      fclose(fich_sortie) ;
      
   
      sprintf(name, "resformat_%e.dat", separ) ;
      ofstream seqfich(name) ;
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
    }
    
    return 1 ;
}
