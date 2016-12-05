/*
 * Test program the contract and flat_scalar_prod operations on Tenseur's
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2016/12/05 16:18:29  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:54:03  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:12:55  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2008/08/27 08:54:07  jl_cornou
 * Tests for angular potential A added
 *
 * Revision 1.3  2003/01/09 11:07:59  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.2  2002/08/08 15:10:45  j_novak
 * The flag "plat" has been added to the class Metrique to show flat metrics.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2000/02/01  14:56:14  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// headers C
#include <cstdlib>
#include <cmath>

// headers Lorene
#include "type_parite.h"
#include "tenseur.h"
#include "utilitaires.h"
#include "nbr_spx.h"

//******************************************************************************

using namespace Lorene ;

int main(){
    
    // Identification of all the subroutines called by the code : 
    
    system("ident test_tenseur_contract") ; 

    //-----------------------------------------------------------------------
    //		Input data : number of points, types of sampling, etc...
    //-----------------------------------------------------------------------

    int nt = 7 ;    // Number of points in theta
    int np = 6 ;    // Number of points in phi
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

    //-----------------------------------------------------------------------
    //		Construction of a mapping
    //-----------------------------------------------------------------------
    
    const Map_af mp(mg, bornes) ;

    const Coord& r = mp.r ; 
//    const Coord& x = mp.x ; 
//    const Coord& y = mp.y ; 
//    const Coord& z = mp.z ; 
//    const Coord& cost = mp.cost ; 
//    const Coord& sint = mp.sint ; 
//    const Coord& cosp = mp.cosp ; 
//    const Coord& sinp = mp.sinp ; 
    
    //-----------------------------------------------------------------------
    //		Construction of a Cmp
    //-----------------------------------------------------------------------

    Cmp aa(mp) ; 
    
    aa = r*r  ;   
    aa.annule(nz-1) ;	
    aa.std_base_scal() ; 
    
    //-----------------------------------------------------------------------
    //		Scalar product of two vectors
    //-----------------------------------------------------------------------
    
    Tenseur uu(mp, 1, COV, mp.get_bvect_cart()) ; 
    Tenseur vv(mp, 1, CON, mp.get_bvect_cart()) ; 

    uu.set_etat_qcq() ; 
 
    uu.set(0) = aa ; 
    uu.set(1) = - aa ; 
    uu.set(2) = 4 * aa ; 

    vv.set_etat_qcq() ; 
 
    vv.set(0) = -3 * aa ; 
    vv.set(1) = aa ; 
    vv.set(2) = 2 * aa ; 


    cout << endl << "Test scalar product of two vectors" << endl 
		 << "==================================" << endl ; 

    Tenseur scal0( 4 * aa * aa ) ; 
    Tenseur scal1( uu(0) * vv(0) + uu(1) * vv(1) + uu(2) * vv(2) ) ; 
    Tenseur scal2 = contract(uu, 0, vv, 0) ; 
    Tenseur scal3 = flat_scalar_prod(uu, vv) ; 
    
    cout << "diff scal1 / scal0 : " << diffrel(scal1(), scal0()) << endl ; 
    cout << "diff scal2 / scal0 : " << diffrel(scal2(), scal0()) << endl ; 
    cout << "diff scal3 / scal0 : " << diffrel(scal3(), scal0()) << endl ; 
    cout << "diff scal2 / scal3 : " << diffrel(scal2(), scal3()) << endl ; 
 
    Tenseur diff = scal2 - scal3 ; 
    diff().affiche_seuil(cout, 1) ; 
     
    //-----------------------------------------------------------------------
    //		Scalar product of a rank 2 tensor and a vector
    //-----------------------------------------------------------------------
    
    Tenseur kk(mp, 2, COV, mp.get_bvect_cart()) ; 

    kk.set_etat_qcq() ; 
 
    kk.set(0, 0) = 1 ; 
    kk.set(0, 1) = 2 ; 
    kk.set(0, 2) = 3 ; 
    kk.set(1, 0) = 4 ; 
    kk.set(1, 1) = 5 ; 
    kk.set(1, 2) = 6 ; 
    kk.set(2, 0) = 7 ; 
    kk.set(2, 1) = 8 ; 
    kk.set(2, 2) = 9 ; 

    cout << endl << "Test scalar product of a rank 2 tensor and a vector" << endl
		 << "===================================================" << endl ; 

    Tenseur res0(mp, 1, COV, mp.get_bvect_cart()) ; 
    res0.set_etat_qcq() ; 
    
    res0.set(0) = 5 * aa ; 
    res0.set(1) = 5 * aa ; 
    res0.set(2) = 5 * aa ; 

    Tenseur res2 = contract(kk, 1, vv, 0) ; 
    Tenseur res3 = flat_scalar_prod(kk, vv) ; 
    
    cout << "diff res2 / res0 : " << endl
	 << diffrel(res2(0), res0(0)) << endl 
	 << diffrel(res2(1), res0(1)) << endl 
	 << diffrel(res2(2), res0(2)) << endl ;  
	 
    cout << "diff res3 / res0 : " << endl
	 << diffrel(res3(0), res0(0)) << endl 
	 << diffrel(res3(1), res0(1)) << endl 
	 << diffrel(res3(2), res0(2)) << endl ;  
	 
    cout << "diff res3 / res2 : " << endl 
	 << diffrel(res3(0), res2(0)) << endl 
	 << diffrel(res3(1), res2(1)) << endl 
	 << diffrel(res3(2), res2(2)) << endl ;  
	 
    Tenseur vdiff = res2 - res3 ; 
    vdiff(0).affiche_seuil(cout, 1) ; 
    vdiff(1).affiche_seuil(cout, 1) ; 
    vdiff(2).affiche_seuil(cout, 1) ; 
     

    // Cleaning
    // --------
    
    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ; 
    delete [] bornes ; 
    delete [] type_r ; 
    
    exit(EXIT_SUCCESS) ; 
    
}

