/*
 * Test program the derivative functions of the Cmp class with a Map_af mapping
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

char test_map_af_deriv_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/01/22 16:18:52  e_gourgoulhon
 * void main() --> int main().
 *
 * Revision 1.2  2003/01/09 11:07:53  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  2000/01/28  08:48:00  eric
 * Identification du code au debut
 *
 * Revision 1.1  2000/01/27  15:25:29  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// headers C
#include <stdlib.h>
#include <math.h>

// headers Lorene
#include "type_parite.h"
#include "cmp.h"
#include "utilitaires.h"
#include "param.h"
#include "graphique.h"
#include "nbr_spx.h"

//******************************************************************************

int main(){
    
    // Identification of all the subroutines called by the code : 
    
    system("ident test_map_af_deriv") ; 

    //-----------------------------------------------------------------------
    //		Input data : number of points, types of sampling, etc...
    //-----------------------------------------------------------------------

    int nt = 7 ;    // Number of points in theta
    int np = 8 ;    // Number of points in phi
    int nz = 3 ;    // Number of domains

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
    int* type_r = new int[nz];

    for (int l=0; l<nz; l++) {
	nr[l] = 9 ;		// Number of points in r
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
    
    type_t = SYM ; 
    type_p = NONSYM ; 
    
    const Mg3d mg(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    //----------------------------
    // Construction of the Map_af
    //----------------------------
    
    Map_af mp(mg, bornes) ;
        
    const Coord& r = mp.r ; 
    const Coord& x = mp.x ; 
    const Coord& y = mp.y ; 
    const Coord& z = mp.z ; 
    const Coord& cost = mp.cost ; 
    const Coord& sint = mp.sint ; 
    const Coord& cosp = mp.cosp ; 
    const Coord& sinp = mp.sinp ; 
    
    //-----------------------------------------------------------------------
    //		Construction of a Cmp 
    //-----------------------------------------------------------------------

    Cmp aa(mp) ; 
    Cmp aa_ext(mp) ;   
    Cmp da_ana(mp) ; 
    Cmp da_ana_ext(mp) ;   
    Cmp da(mp) ; 

    aa = x*y*z*z + y  ;		
    aa_ext = 1 / (r*r) ; 
    
    aa.annule(nz-1) ; 
    aa_ext.annule(0, nz-2) ; 
    
    aa = aa + aa_ext ; 
    
    aa.std_base_scal() ; // Sets the standard basis for spectral expansions
         
    cout << "Coefficients of aa = x y z^2 + y  (in ZEC : aa = 1/r^2) : " << endl ; 
    aa.affiche_seuil(cout) ;
     
    arrete() ; 

    // Computation of da/dx 
    // --------------------
    
    da = aa.dsdx() ; 
    
    // Analytical expression for da/dx :     
    da_ana = y*z*z ;  
    da_ana_ext = - 2 * sint * cosp / r ;

    da_ana.annule(nz-1) ; 
    da_ana_ext.annule(0, nz-2) ; 
    
    da_ana = da_ana + da_ana_ext ;
    da_ana.set_dzpuis(2) ;  
    da_ana.std_base_scal() ; 
    
    cout << "Coefficients of da/dx : " << endl ; 
    da.affiche_seuil(cout) ;
     
    cout << "Coefficients of da/dx analytical : " << endl ; 
    da_ana.affiche_seuil(cout) ;
     
    cout << "Relative difference da/dx <-> da/dx analytical : " << endl ; 
    cout << diffrel(da, da_ana) << endl ; 
    
    cout << "Difference coefficients da/dx <-> da/dx analytical : " << endl ; 
    (da.va).coef() ; 
    (da_ana.va).coef() ; 
    Mtbl_cf diff = *((da.va).c_cf) - *((da_ana.va).c_cf) ; 
    diff.affiche_seuil(cout) ; 
        
    arrete() ; 
    
    // Computation of da/dy 
    // --------------------
    
    da = aa.dsdy() ; 
    
    // Analytical expression for da/dy :     
    da_ana = x*z*z + 1 ;  
    da_ana_ext = - 2 * sint * sinp / r ;

    da_ana.annule(nz-1) ; 
    da_ana_ext.annule(0, nz-2) ; 
    
    da_ana = da_ana + da_ana_ext ;
    da_ana.set_dzpuis(2) ;  
    da_ana.std_base_scal() ; 
    
    cout << "Coefficients of da/dy : " << endl ; 
    da.affiche_seuil(cout) ;
     
    cout << "Coefficients of da/dy analytical : " << endl ; 
    da_ana.affiche_seuil(cout) ;
     
    cout << "Relative difference da/dy <-> da/dy analytical : " << endl ; 
    cout << diffrel(da, da_ana) << endl ; 
    
    cout << "Difference coefficients da/dy <-> da/dy analytical : " << endl ; 
    (da.va).coef() ; 
    (da_ana.va).coef() ; 
    diff = *((da.va).c_cf) - *((da_ana.va).c_cf) ; 
    diff.affiche_seuil(cout) ; 
    
    arrete() ; 
    
    // Computation of da/dz 
    // --------------------
    
    da = aa.dsdz() ; 
    
    // Analytical expression for da/dz :     
    da_ana = 2*z*x*y ;  
    da_ana_ext = - 2 * cost / r ;

    da_ana.annule(nz-1) ; 
    da_ana_ext.annule(0, nz-2) ; 
    
    da_ana = da_ana + da_ana_ext ;
    da_ana.set_dzpuis(2) ;  
    (da_ana.va).set_base( (da.va).base )  ; 
    
    cout << "Coefficients of da/dz : " << endl ; 
    da.affiche_seuil(cout) ;
     
    cout << "Coefficients of da/dz analytical : " << endl ; 
    da_ana.affiche_seuil(cout) ;
     
    cout << "Relative difference da/dz <-> da/dz analytical : " << endl ; 
    cout << diffrel(da, da_ana) << endl ; 
        
    cout << "Difference coefficients da/dz <-> da/dz analytical : " << endl ; 
    (da.va).coef() ; 
    (da_ana.va).coef() ; 
    diff = *((da.va).c_cf) - *((da_ana.va).c_cf) ; 
    diff.affiche_seuil(cout) ; 

    // Cleaning
    // --------
    
    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ; 
    delete [] bornes ; 
    delete [] type_r ; 
    
    return EXIT_SUCCESS ; 
    
}

