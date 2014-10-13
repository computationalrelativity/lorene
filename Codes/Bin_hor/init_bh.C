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

char init_bh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.9  2014/10/13 08:53:53  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2014/10/06 15:09:42  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2007/04/13 15:30:58  f_limousin
 * Lots of improvements, generalisation to an arbitrary state of
 * rotation, implementation of the spatial metric given by Samaya.
 *
 * Revision 1.6  2006/06/28 13:36:52  f_limousin
 * Convergence to a given irreductible mass
 *
 * Revision 1.5  2005/03/04 09:42:02  f_limousin
 * New construction of the object Bin_hor.
 *
 * Revision 1.4  2005/02/25 12:31:59  f_limousin
 * The boundary conditions for psi, N and beta are now parameters in
 * par_init.d and par_coal.d.
 *
 * Revision 1.3  2005/01/03 07:59:02  f_limousin
 * Drawings.
 *
 * Revision 1.2  2004/12/31 15:45:26  f_limousin
 * Change the parameters in par_init.d
 *
 * Revision 1.1  2004/12/29 18:00:20  f_limousin
 * First version
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


using namespace Lorene ;

int main() {
    
    char blabla [120] ;
    ifstream param("par_init.d") ;
    
    double  precis, relax, radius, separation, lim_nn ;
    int nz, nt, np, nr1, nrp1, bound_nn, bound_psi ;
    
    param.getline(blabla, 120) ;
    param.getline(blabla, 120) ;
    param >> separation ; param.getline(blabla, 120) ;
    param >> nz ; param.getline(blabla, 120) ;
    param >> nt; param.ignore(1000, '\n');
    param >> np; param.ignore(1000, '\n');
    param >> nr1; param.ignore(1000, '\n');
    param >> nrp1; param.ignore(1000, '\n');

    double* bornes = new double[nz+1] ;
    int* nr_tab = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];

    for (int l=0 ; l<nz ; l++){
      if (l==1) nr_tab[1] = nr1 ;
      else nr_tab[l] = nrp1 ;
      np_tab[l] = np ; 
      nt_tab[l] = nt ; 
      param >> bornes[l] ;

    }
    radius = bornes[1] ;
    param.getline(blabla, 120) ;
    bornes[nz] = __infinity ; 

    param >> precis ; param.getline(blabla, 120) ;
    param >> relax ; param.getline(blabla, 120) ;    
    param >> bound_nn ;
    param >> lim_nn ;  param.ignore(1000, '\n');
    param >> bound_psi ;  param.ignore(1000, '\n');
    
    param.close() ;
    
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

    Bin_hor bin (map_un, map_deux) ;
    bin.set_statiques(precis, relax, bound_nn, lim_nn, bound_psi) ;
    
    FILE* fich = fopen("static.d", "w") ;
    grid.sauve(fich) ;
    map_un.sauve(fich) ;
    map_deux.sauve(fich) ;
    bin.sauve(fich) ;
    fwrite_be(&bound_nn, sizeof(int), 1, fich) ;
    fwrite_be (&lim_nn, sizeof(double), 1, fich) ;
    fwrite_be(&bound_psi, sizeof(int), 1, fich) ;
    fclose(fich) ;

//    cout << bin << endl ;


    // Drawings
    const Coord& r = bin(1).get_mp().r ;        // r field 
    Mtbl usr = 1 / r ;
    Scalar unsr(bin(1).get_mp()) ;
    unsr = usr ;
    
    Scalar temp = 1. + unsr ;
    temp.std_spectral_base() ;
/*
    des_profile(bin(1).nn(), 1.00001, 10, M_PI/2., 0., "bin(1).nn()") ;
    des_profile(bin(1).psi(), 1.00001, 10, M_PI/2., 0., "bin(1).psi()") ;
    des_profile(temp, 1.00001, 10, M_PI/2., 0., "psi ana()") ;
    des_profile(temp-bin(1).psi(), 1.00001, 10, M_PI/2., 0., "diff psi") ;
*/
    delete [] nr_tab ;
    delete [] nt_tab ;
    delete [] np_tab ;
    delete [] type_r ;
    delete [] bornes ;

    return 1 ;
}
