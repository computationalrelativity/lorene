/*
 * Test program the derivative functions of the Cmp class with a Map_et mapping
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

char test_map_et_deriv_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:31  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.2  2000/01/28  08:43:01  eric
 * Identification du code au debut.
 *
 * Revision 1.1  2000/01/27  15:16:16  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */


// headers C++
#include <iostream.h>

// headers C
#include <stdlib.h>
#include <math.h>

// headers Lorene
#include "type_parite.h"
#include "cmp.h"
#include "utilitaires.h"
#include "param.h"
#include "nbr_spx.h"

//******************************************************************************

void main(){
    
    // Identification of all the subroutines called by the code : 
    
    system("ident test_map_et_deriv") ; 

    //-----------------------------------------------------------------------
    //		Input data : number of points, types of sampling, etc...
    //-----------------------------------------------------------------------

    int nt = 25 ;    // Number of points in theta
    int np = 24 ;    // Number of points in phi
    int nz = 3 ;    // Number of domains

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
    int* type_r = new int[nz];

    for (int l=0; l<nz; l++) {
	nr[l] = 33 ;		// Number of points in r
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

    //-----------------------------------------------------------------------
    //	Construction of a pure angular mapping of class Map_af
    //-----------------------------------------------------------------------
    
    Map_af mp0( *(mg.get_angu()), bornes ) ;

    cout << endl << "Mapping mp0 : " << endl ; 
    cout << mp0 << endl ; 

    const Coord& cost0 = mp0.cost ; 
    const Coord& sint0 = mp0.sint ; 
    const Coord& cosp0 = mp0.cosp ; 
    const Coord& sinp0 = mp0.sinp ; 
    
    //-----------------------------------------------------------------------
    //   Functions F_l and G_l corresponding to a triaxial ellipsoid
    //-----------------------------------------------------------------------

    //---------------------------
    //  Triaxial ellipsoid		x^2/a^2 + y^2/b^2 + z^2/c^2 = 1
    //---------------------------

    Valeur ff0( *(mg.get_angu()) ) ; 
    Valeur gg0( *(mg.get_angu()) ) ; 

    double bsa = 0.8 ;	    // axis ratio b/a 
    double csa = 0.6 ;	    // axis ratio c/a 

    gg0 = 1. / sqrt( sint0*sint0 * ( cosp0*cosp0 + sinp0*sinp0 / (bsa*bsa) )
			+ cost0*cost0 / (csa*csa) ) - 1. ;
    gg0.annule(nz-2, nz-1) ;	// 0 in the external domains

    gg0.std_base_scal() ; 
    
    // F in domain 0 : 
    ff0.set_etat_c_qcq() ; 
    ff0.c->set_etat_qcq() ; 
    ff0.annule(0) ;	// the ellipsoid contains only even m
    
    // F in other domains is deduced from G by maching domain l and domain l+1 :

    const double* alpha = mp0.get_alpha() ; 
    
    for (int l=1; l<nz-1; l++) {
	*(ff0.c->t[l]) = alpha[l-1] / alpha[l] * (*(gg0.c->t[l-1])) ;
    } 

    ff0.annule(nz-1) ; 

    ff0.std_base_scal() ; 

    //----------------------------
    // Construction of the Map_et
    //----------------------------
    
    Map_et mp(mg, bornes) ;
        
    mp.set_ff(ff0) ; 
    mp.set_gg(gg0) ; 
    
    cout << endl << "Mapping mp : " << endl ; 
    cout << mp << endl ; 
    
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
         
//    cout << "Coefficients of aa = x y z^2 + y  (in ZEC : aa = 1/r^2) : " << endl ; 
//    aa.affiche_seuil(cout) ;     
//    arrete() ; 

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
    
//    cout << "Coefficients of da/dx : " << endl ; 
//    da.affiche_seuil(cout, 0, 4, 1e-4) ;
     
//    cout << "Coefficients of da/dx analytical : " << endl ; 
//    da_ana.affiche_seuil(cout, 0, 4, 1e-4) ;
     
    cout << "Relative difference da/dx <-> da/dx analytical : " << endl ; 
    cout << diffrel(da, da_ana) << endl ; 
    
    cout << "Difference coefficients da/dx <-> da/dx analytical : " << endl ; 
    (da.va).coef() ; 
    (da_ana.va).coef() ; 
    Mtbl_cf diff = *((da.va).c_cf) - *((da_ana.va).c_cf) ; 
    cout << norme(diff) / norme( *((da_ana.va).c_cf) ) << endl ; 
        
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
    
//    cout << "Coefficients of da/dy : " << endl ; 
//    da.affiche_seuil(cout, 0, 4, 1e-4) ;
     
//    cout << "Coefficients of da/dy analytical : " << endl ; 
//    da_ana.affiche_seuil(cout, 0, 4, 1e-4) ;
     
    cout << "Relative difference da/dy <-> da/dy analytical : " << endl ; 
    cout << diffrel(da, da_ana) << endl ; 
    
    cout << "Difference coefficients da/dy <-> da/dy analytical : " << endl ; 
    (da.va).coef() ; 
    (da_ana.va).coef() ; 
    diff = *((da.va).c_cf) - *((da_ana.va).c_cf) ; 
    cout << norme(diff) / norme( *((da_ana.va).c_cf) ) << endl ; 
    
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
    
//    cout << "Coefficients of da/dz : " << endl ; 
//    da.affiche_seuil(cout, 0, 4, 1e-4) ;
     
//    cout << "Coefficients of da/dz analytical : " << endl ; 
//    da_ana.affiche_seuil(cout, 0, 4, 1e-4) ;
     
    cout << "Relative difference da/dz <-> da/dz analytical : " << endl ; 
    cout << diffrel(da, da_ana) << endl ; 
        
    cout << "Difference coefficients da/dz <-> da/dz analytical : " << endl ; 
    (da.va).coef() ; 
    (da_ana.va).coef() ; 
    diff = *((da.va).c_cf) - *((da_ana.va).c_cf) ; 
    cout << norme(diff) / norme( *((da_ana.va).c_cf) ) << endl ; 

    // Cleaning
    // --------
    
    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ; 
    delete [] bornes ; 
    delete [] type_r ; 
    
    exit(EXIT_SUCCESS) ; 
    
}

