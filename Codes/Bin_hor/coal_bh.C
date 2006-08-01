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
 * Revision 1.11  2006/08/01 14:13:41  f_limousin
 * New version...
 *
 * Revision 1.10  2006/06/29 08:54:51  f_limousin
 * Boundary conditions and grid writen in resformat.dat
 *
 * Revision 1.9  2006/06/28 13:36:52  f_limousin
 * Convergence to a given irreductible mass
 *
 * Revision 1.8  2006/05/24 16:59:08  f_limousin
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
    char nomini[120] ;
    double omega_init, relax, precis_viriel, lim_nn, mass_irr ;
    int nb_om, nb_it, bound_nn, bound_psi, bound_beta, search_mass ;
    
    ifstream param("par_coal.d") ;
	if ( !param.good() ) {
		cout << "Problem with opening the file par_coal.d ! " << endl ;
		abort() ;
	}
    param.ignore(1000, '\n') ;
    param.ignore(1000, '\n') ;
    param.getline(nomini, 80) ; 
    param >> omega_init ; param.getline(blabla, 120) ;
    param >> search_mass ;
    param >> mass_irr ;  param.ignore(1000, '\n');
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
    
    bin.set_omega(0) ;
    bin.set(1).n_comp (bin(2)) ;
    bin.set(1).psi_comp (bin(2)) ;
    bin.set(1).beta_comp (bin(2)) ;
    bin.set(2).n_comp (bin(1)) ;
    bin.set(2).psi_comp (bin(1)) ;
    bin.set(2).beta_comp (bin(1)) ;
    bin.decouple() ;
    bin.extrinsic_curvature() ;
    
    double separation = bin(1).get_mp().get_ori_x() - 
	bin(2).get_mp().get_ori_x() ;
        
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

    return EXIT_SUCCESS ;
}
