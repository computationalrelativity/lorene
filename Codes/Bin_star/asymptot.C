/*
 * Subroutine for analysis of asymptotic behavior of a Cmp
 *
 */

/*
 *   Copyright (c) 2002 Francois Limousin
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
char asymptot_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/06/18 14:07:26  f_limousin
 * Analysis of asymptotic behavior of binary NS and BH
 *
 *
 *
 * $Header$
 *
 */

// Headers standard du C++
#include <iostream.h>
#include <fstream.h>
 
// Headers standard du C
//  (par exemple definit la macro EXIT_SUCCESS)
#include <stdlib.h>
#include <stdio.h>

 
// Headers Lorene
#include "tenseur.h"
#include "bhole.h"
#include "binaire.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"

#include "unites.h"

void asymptot(const Cmp& nn, const char* coment, bool graphics, ostream& fich) {

  // Multi-grid
  const Mg3d& mg =  *( nn.get_mp()->get_mg() ) ;
  int nz = mg.get_nzone() ;
  int nzm1 = nz - 1 ; // index of the last domain

 // Number of coefficients
  int np = mg.get_np(nzm1) ;
  int nt = mg.get_nt(nzm1) ;
 
    // Asymptotic behavior of N :
        Valeur** nn_asymp = nn.asymptot(3,1) ; 
       

	fich << "Asymptotic behavior of " << coment << endl << endl ;
	
	for (int i=0; i<4; i++) {
	// Value (on the angular grid) containing the coef of 1/r^i
	Valeur& nn_i = *(nn_asymp[i]) ;
	
	// Computation of spectral expansions 
	nn_i.coef() ;

	cout << "Spectral coefficients of " << coment << " (1/r^" << i 
	     << ") : " << endl ;

	nn_i.affiche_seuil(cout,0,4,1e-5) ;

	fich << "Spectral coefficients of " << coment << " (1/r^" << i 
	     << ") : " << endl ;

		
	double nni_max = 0. ;
	for (int k=0; k<np+1; k++) {
	    if (k==1) continue ; 
	    for (int j=0; j<nt; j++) {
		double cf = (*nn_i.c_cf)(nzm1,k,j,0) ;
		if (fabs(cf) > fabs (nni_max)) {
		    nni_max = fabs(cf) ;
		}
	    }
	}

	fich << "nn" << i << "_max =" << nni_max << endl ;
	//fich << "nn" << i << "_min =" << nni_min << endl ;	    
	
	if ( nni_max > 1e-3 ) {
	for (int k=0; k<np+1; k++) {
	    if (k==1) continue ; 
	    for (int j=0; j<nt; j++) {
		double cf = (*nn_i.c_cf)(nzm1,k,j,0) ; 
		if ( fabs(cf) > 0.01 * nni_max && nni_max !=0) {
		    fich << "k= " << k << " j= " << j << " : " << cf << endl ;
		}
	    }
	}
	}
	
	    if (graphics && i==3) {
	  des_coef_theta(nn_i, 1, 0, 0, 1e-10) ; 
	  des_coef_theta(nn_i, 1, 2, 0, 1e-10) ; 
	  des_coef_theta(nn_i, 1, 4, 0, 1e-10) ; 

	  des_coef_theta(nn_i, 1, 6, 0, 1e-10) ; 
	  des_coef_theta(nn_i, 1, 8, 0, 1e-10) ; 
    
	  des_coef_phi(nn_i, 1, 0, 0, 1e-10) ; 
	  des_coef_phi(nn_i, 1, 1, 0, 1e-10) ; 
	  des_coef_phi(nn_i, 1, 2, 0, 1e-10) ; 
	  des_coef_phi(nn_i, 1, 3, 0, 1e-10) ; 
	  des_coef_phi(nn_i, 1, 4, 0, 1e-10) ; 
	    
	    }
	    fich << endl ;
	
      //double mass_nn = - nn_i(nz -1,0,0,0) / ggrav ; 
      //double mass_nn_sol = mass_nn / msol ;
      //cout << "Mass read from lapse : " << mass_nn_sol << " M_sol" << endl ; 

        //nn_3.ylm() ; // Spherical harmonics

	}
	    fich << "------------------------------------------------------------------------" << endl ;

}





  









