/*
 * Test program for Map_radial::poisson_compact 
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

char test_poisson_compact_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:31  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.4  2000/02/25  09:56:04  eric
 * Nouveau constructeur de Tenseur (triade passee en argument)
 * fact_echelle est desormais passe en double_mod dans le Param.
 *
 * Revision 1.3  2000/01/27  16:04:50  eric
 * Suppression lecture des parametres dans un fichier.
 *
 * Revision 1.2  2000/01/14  17:35:33  eric
 * Annulation de la source avant l'appel a Map_radial::poisson_compact.
 *
 * Revision 1.1  2000/01/13  16:45:10  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */


// version of 13.01.2000

// headers C++
#include <iostream.h>
#include <fstream.h>

// headers C
#include <stdlib.h>
#include <math.h>

// headers Lorene
#include "tenseur.h"
#include "param.h"
#include "utilitaires.h"
#include "nbr_spx.h"

//******************************************************************************

void main(){
    
    // Identification of all the subroutines called by the code : 
    
    system("ident test_poisson_compact") ; 

    //-----------------------------------------------------------------------
    //		Input data : number of points, types of sampling, etc...
    //-----------------------------------------------------------------------

    int nt = 17 ;    // Number of points in theta
    int np = 16 ;    // Number of points in phi
    int nz = 3 ;    // Number of domains

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
    int* type_r = new int[nz];

    for (int l=0; l<nz; l++) {
	nr[l] = 25 ;		// Number of points in r
	np_tab[l] = np ; 
	nt_tab[l] = nt ; 
    }
    
    bornes[0] = 0. ; 
    bornes[1] = 2. ;
    bornes[2] = 3. ;
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
    
    Mg3d mg(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    cout << endl << "Grid mg : " << mg << endl ; 
        
    //----------------------------
    // Construction of the Map_et
    //----------------------------
    
    Map_et mp(mg, bornes) ;
        
    const Coord& r = mp.r ; 
    const Coord& x = mp.x ; 
    const Coord& y = mp.y ; 
    const Coord& z = mp.z ; 
//    const Coord& cost = mp.cost ; 
    const Coord& sint = mp.sint ; 
//    const Coord& cosp = mp.cosp ; 
//    const Coord& sinp = mp.sinp ; 

    //-----------------------------------------------------------------------
    //		Construction of a Cmp
    //-----------------------------------------------------------------------

    Cmp ent(mp) ; 
    
    double bsa = 0.9 ; 
    double csa = 0.8 ; 
    double xtri = - 0. ; 
    double ray0 = mp.val_r(0, 1., M_PI/2, 0) ;

    double a = sqrt( ray0 * (ray0 + xtri) ) ; 
    double b = bsa * a ; 
    double c = csa * a ; 
    
    cout << "a, b, c, xtri : " << a << "  " << b << "  " << c << "  "
	  << xtri << endl ; 
	  
	     
    ent = 1 - (x+xtri)*x/(a*a) - (y*y)/(b*b) - (z*z)/(c*c) ; 

    ent.annule(nz-1) ; 
    ent.std_base_scal() ; 
    
       
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
    
    Param par ; 
    par.add_int(nitermax, 0) ;  
    par.add_int(nzet, 1) ;    
    par.add_int(1, 2) ; // 1 = full computation
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

    cout << endl << "Coefficients of the mapping functions F and G" << endl ; 
//    des_map_et(mp, 0) ;      
//    des_map_et(mp, 1) ;      


    // Computation of the field ent at the new collocation points
    // ----------------------------------------------------------
    
    mp.reevaluate(&mp_prev, nz-1, ent) ;  


//    des_coupe_x(ent, 0., -3., 3., -3., 3., "Enthalpy (x=0)",  
//		&ent, 25) ; 
//    des_coupe_y(ent, 0., -3., 3., -3., 3., "Enthalpy (y=0)",  
//		&ent, 25) ; 
//    des_coupe_z(ent, 0., -3., 3., -3., 3., "Enthalpy (z=0)",  
//		&ent, 25) ; 
  
    //-----------------------------------------------------------------------
    //		Resolution of the continuity equation
    //-----------------------------------------------------------------------
    
    Cmp psi(mp) ;	// velocity scalar potential
    
    double gamma = 2 ; 
    
    Cmp aa = (gamma-1) * ent ; 
    
    Tenseur ent1(ent) ; 
    Tenseur bb = ent1.gradient_spher() ; 
    bb.annule(1, nz-1) ; 
    
    double omega = 1. ; 
    
    Tenseur vit_rot(mp, 1, CON, mp.get_bvect_spher()) ; 
    vit_rot.set_etat_qcq() ; 
    vit_rot.set(0) = 0 ; 
    vit_rot.set(1) = 0 ; 
    vit_rot.set(2) = omega * r * sint ; 
    vit_rot.annule(1, nz-1) ; 
    
    Tenseur source1 = contract(vit_rot, 0, bb, 0) ; 
    
    double relax = 0.5 ; 
    nitermax = 100 ; 
    precis = 1.e-8 ; 
    
    Param par_com ; 
    par_com.add_int(nitermax) ;  
    par_com.add_int_mod(niter) ;  
    par_com.add_double(precis, 0) ; 
    par_com.add_double(relax, 1) ; 

    psi = 0 ; 

    mp.poisson_compact(source1(), aa, bb, par_com, psi) ;
    
    cout << "Coef of psi : " << endl ; 
    (psi.va).affiche_seuil(cout) ; 
    
    //-----------------------------------------------------------------------
    //		Test: does the solution satisfies to the equation ?
    //-----------------------------------------------------------------------
    
    Tenseur psi1(psi) ; 
    Tenseur dpsi = psi1.gradient_spher() ; 
    Cmp bdpsi = bb(0)*dpsi(0) + bb(1)*dpsi(1) + bb(2)*dpsi(2) ; 
    
    cout << "Coef of bdpsi : " << endl ; 
    (bdpsi.va).affiche_seuil(cout) ; 
    


    Cmp source = source1() ; 
    source.set_dzpuis(4) ; 
    source.annule(1, nz-1) ; 
    
    Cmp op_psi = aa * psi.laplacien() + bdpsi ; 
    
    cout << "Coef of op_psi : " << endl ; 
    (op_psi.va).affiche_seuil(cout) ; 
    
    Cmp diff = op_psi - source ; 
    
    cout << "Coef of diff : " << endl ; 
    (diff.va).affiche_seuil(cout) ; 
    
    cout << "diffrel(op_psi, source) : " << endl ; 
    cout << diffrel(op_psi, source) << endl ; 
    
    cout << "diffrelmax(op_psi, source) : " << endl ; 
    cout << diffrelmax(op_psi, source) << endl ; 
    
    
    
    // Cleaning
    // --------
    
    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ; 
    delete [] bornes ; 
    delete [] type_r ; 
    
    exit(EXIT_SUCCESS) ; 
    
}

