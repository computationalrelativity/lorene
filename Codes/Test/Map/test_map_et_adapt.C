/*
 * Test program for the adapt and reevaluate methods of the Map_et class
 * 
 */
 
/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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

char test_map_et_adapt_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2015/08/10 15:32:27  j_novak
 * Better calls to Param::add_int(), to avoid weird problems (e.g. with g++ 4.8).
 *
 * Revision 1.3  2014/10/06 15:12:53  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/01/09 11:07:53  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2000/02/15  15:20:08  eric
 * *** empty log message ***
 *
 * Revision 1.2  2000/02/15  15:18:29  eric
 * Changement du parametre Param dans Map_et::adapt : fact_echelle est
 * desormais passe en double_mod.
 *
 * Revision 1.1  2000/01/04  16:48:56  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// headers C
#include <cstdlib>

// headers Lorene
#include "type_parite.h"
#include "cmp.h"
#include "param.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"

//******************************************************************************

void main(){
    
    // Identification of all the subroutines called by the code : 
    
    system("ident test_map_et_adapt") ; 

    //-----------------------------------------------------------------------
    //		Input data : number of points, types of sampling, etc...
    //-----------------------------------------------------------------------

    int nt = 9 ;    // Number of points in theta
    int np = 8 ;    // Number of points in phi
    int nz = 3 ;    // Number of domains

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
    int* type_r = new int[nz];

    for (int l=0; l<nz; l++) {
	nr[l] = 13 ; 
	np_tab[l] = np ; 
	nt_tab[l] = nt ; 
    }
    
    bornes[0] = 0. ; 
    bornes[1] = 2. ;
    bornes[2] = 3.5 ;
    bornes[3] = __infinity ; 
    
    type_r[0] = RARE ; 
    type_r[1] = FIN ; 
    type_r[2] = UNSURR ; 
    
    int type_t = SYM ; 

    int type_p = NONSYM ; 
    

    //-----------------------------------------------------------------------
    //		Construction of a multi-grid
    //-----------------------------------------------------------------------
    
    Mg3d mg(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    cout << endl << "Grid mg : " << mg << endl ; 
        

    //----------------------------
    // Construction of the Map_et
    //----------------------------
    
    Map_et mp(mg, bornes) ;
    
    cout << endl << "Mapping mp (before adaptation) : " << endl ; 
    cout << mp << endl ; 
    
//    const Coord& r = mp.r ; 
    const Coord& x = mp.x ; 
    const Coord& y = mp.y ; 
    const Coord& z = mp.z ; 
//    const Coord& cost = mp.cost ; 
//    const Coord& sint = mp.sint ; 
//    const Coord& cosp = mp.cosp ; 
//    const Coord& sinp = mp.sinp ; 

    //-----------------------------------------------------------------------
    //		Construction of a Cmp
    //-----------------------------------------------------------------------

    Cmp ent(mp) ; 

    ent = 1 - x*x/9 - y*y/4 - z*z ;
    ent.annule(nz-1) ; 
    ent.std_base_scal() ; 

    
    des_coupe_x(ent, 0., -3., 3., -3., 3., "Before mapping adaptation (x=0)", 
		&ent, 25) ; 
    des_coupe_y(ent, 0., -3., 3., -3., 3., "Before mapping adaptation (y=0)", 
		&ent, 25) ; 
    des_coupe_z(ent, 0., -3., 3., -3., 3., "Before mapping adaptation (z=0)", 
		&ent, 25) ; 

   
    //-----------------------------------------------------------------------
    //		Adaptation of the Map_et to the Cmp
    //-----------------------------------------------------------------------
    
    // Save the old value of the mapping : 
    Map_et mp_prev = mp ; 

    int nzet = 1 ;
    int nitermax = 100 ;  
    int niter ; 
    double precis = 1.e-15 ; 
    int j_bord = mg.get_nt(0) - 1 ; 
    int k_bord = 0 ; 
    double fact_echelle = 1 ; 
    double fact_lamu = 1 ;
    
    Tbl ent_limit(nzet) ; 
    ent_limit.set_etat_qcq() ; 
    ent_limit.set(0) = 0 ; 
    
    int i_one = 1 ;

    Param par ; 
    par.add_int(nitermax, 0) ;  
    par.add_int(nzet, 1) ;    
    par.add_int(i_one, 2) ; // 1 = full computation
    par.add_int(j_bord, 3) ; 
    par.add_int(k_bord, 4) ; 
    par.add_int_mod(niter) ;  
    par.add_double(precis, 0) ; 
    par.add_double(fact_lamu, 1) ; 
    par.add_double_mod(fact_echelle, 0) ; 
    par.add_tbl(ent_limit, 0) ; 
    
    mp.adapt(ent, par) ;
    
    cout << endl << "Mapping mp (after adaptation) : " << endl ; 
    cout << mp << endl ; 
    arrete() ; 
    

    cout << endl << "Coefficients of the mapping functions F and G" << endl ; 
    des_map_et(mp, 0) ;      
    des_map_et(mp, 1) ;      


    des_coupe_x(ent, 0., -3., 3., -3., 3., 
	"After mapping adaptation and before reevaluate (x=0)",  
		&ent, 25) ; 
    des_coupe_y(ent, 0., -3., 3., -3., 3., 
	"After mapping adaptation and before reevaluate (y=0)",  
		&ent, 25) ; 
    des_coupe_z(ent, 0., -3., 3., -3., 3., 
	"After mapping adaptation and before reevaluate (z=0)",  
		&ent, 25) ; 

    
    // Computation of the field ent at the new collocation points
    // ----------------------------------------------------------
    
    
    mp.reevaluate(&mp_prev, nz-1, ent) ;  


    des_coupe_x(ent, 0., -3., 3., -3., 3., "After reevaluate (x=0)",  
		&ent, 25) ; 
    des_coupe_y(ent, 0., -3., 3., -3., 3., "After reevaluate (y=0)",  
		&ent, 25) ; 
    des_coupe_z(ent, 0., -3., 3., -3., 3., "After reevaluate (z=0)",  
		&ent, 25) ; 
  

    Cmp ent0(mp) ; 
    ent0 = 1 - x*x/9 - y*y/4 - z*z ;
    ent0.annule(nz-1) ; 
    ent0.std_base_scal() ; 

    cout << "Relative difference between ent and ent0 : " << endl ; 
    cout << "  diffrel : " << diffrel(ent, ent0) << endl ; 
    cout << "  diffrelmax : " << diffrelmax(ent, ent0) << endl ; 
    
    
    
    
    // Cleaning
    // --------
    
    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ; 
    delete [] bornes ; 
    delete [] type_r ; 
    
    exit(EXIT_SUCCESS) ; 
    
}

