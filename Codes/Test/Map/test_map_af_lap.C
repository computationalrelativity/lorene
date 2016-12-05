/*
 * Test program the Laplacian operator of the Cmp class
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
 * Revision 1.4  2016/12/05 16:18:28  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
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
 * Revision 1.1  2000/02/25  09:21:36  eric
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
#include "nbr_spx.h"

//******************************************************************************

void main(){
    
    // Identification of all the subroutines called by the code : 
    
    system("ident test_map_af_lap") ; 

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
    
    Mg3d mg(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    //----------------------------
    // Construction of the Map_af
    //----------------------------
    
    Map_af mp(mg, bornes) ;
        
    const Coord& r = mp.r ; 
    const Coord& x = mp.x ; 
    const Coord& y = mp.y ; 
    const Coord& z = mp.z ; 
    
    //-----------------------------------------------------------------------
    //		Construction of a Cmp
    //-----------------------------------------------------------------------

    Cmp aa(mp) ; 
    Cmp aa_ext(mp) ;     
    
    aa = x*y*z*z + x*x  ;   // aa = x*y*z*z + x*x in all space
    aa.annule(nz-1) ;	// aa = x*y*z*z + x*x in internal domains, 
			// aa = 0 in extern. dom.

    aa_ext = 1 / (r*r) ;    // aa_ext = 1/r^2 in all space
    aa_ext.annule(0, nz-2) ;	// aa_ext = 0 in intern. dom., =1/r^2 in ext. dom.
    
    aa = aa + aa_ext ;	// aa = x*y*z*z + x*x in internal domains, 
			// aa = 1/r^2 in external domain 

    aa.std_base_scal() ; // Sets the standard basis for spectral expansions

    cout << "Coef de aa : " << endl ; 
    aa.affiche_seuil(cout) ; 
    
		// Computation of Laplacian of a with zec_mult_r = 4
		// -------------------------------------------------
    
    Cmp lap = aa.laplacien() ;  // lap = Delta(a)	 in internal domains
				// lap = r^4 Delta(a)   in the external domain
    
    cout << "Coef of lap (zec_mult_r=4) : " << endl ; 
    lap.affiche_seuil(cout) ; 

    // Analytical expression for the Laplacian : 
    Cmp lap_ana(mp) ; 
    Cmp lap_ana_ext(mp) ;
    
    lap_ana = 2 + 2*x*y ;  
    lap_ana.annule(nz-1) ;  

    lap_ana_ext = 2  ;           // 2/r^4 (dzpuis = 4)
    lap_ana_ext.annule(0, nz-2) ; 
    
    lap_ana = lap_ana + lap_ana_ext ; 
    (lap_ana.va).c->dzpuis = 4 ;  
    lap_ana.set_dzpuis(4) ; 
    lap_ana.std_base_scal() ; // Sets the standard base for spectral expansions
    
    cout << "Coef de lap_ana (dzpuis=4) : " << endl ; 
    lap_ana.affiche_seuil(cout) ; 

    cout << "diffrel lap <-> lap_ana : " << endl ; 
    cout << diffrel(lap, lap_ana) << endl ; 
    
    cout << "diffrelmax lap <-> lap_ana : " << endl ; 
    cout << diffrelmax(lap, lap_ana) << endl ; 
    
    
		// Computation of Laplacian of a with zec_mult_r = 2
		// -------------------------------------------------
    
    lap = aa.laplacien(2) ;  // lap = Delta(a)	 in internal domains
				// lap = r^2 Delta(a)   in the external domain
    
    cout << "Coef of lap (zec_mult_r=2) : " << endl ; 
    lap.affiche_seuil(cout) ; 

    // Analytical expression for the Laplacian : 
    lap_ana = 2 + 2*x*y ;  
    lap_ana.annule(nz-1) ;  

    lap_ana_ext = 2/(r*r)  ;           // 2/r^4 (dzpuis = 2)
    lap_ana_ext.annule(0, nz-2) ; 
    
    lap_ana = lap_ana + lap_ana_ext ; 
    (lap_ana.va).c->dzpuis = 2 ;  
    lap_ana.set_dzpuis(2) ; 
    lap_ana.std_base_scal() ; // Sets the standard base for spectral expansions
    
    cout << "Coef de lap_ana (dzpuis=2) : " << endl ; 
    lap_ana.affiche_seuil(cout) ; 

    cout << "diffrel lap <-> lap_ana : " << endl ; 
    cout << diffrel(lap, lap_ana) << endl ; 
    
    cout << "diffrelmax lap <-> lap_ana : " << endl ; 
    cout << diffrelmax(lap, lap_ana) << endl ; 
    
    
		// Computation of Laplacian of a with zec_mult_r = 0
		// -------------------------------------------------
    
    lap = aa.laplacien(0) ;  // lap = Delta(a)	 in internal domains
				//  lap = r^0 Delta(a)   in the external domain
    
    cout << "Coef of lap (zec_mult_r=0) : " << endl ; 
    lap.affiche_seuil(cout) ; 

    // Analytical expression for the Laplacian : 
    lap_ana = 2 + 2*x*y ;  
    lap_ana.annule(nz-1) ;  

    lap_ana_ext = 2/pow(r, 4)  ;           // 2/r^4 (dzpuis = 0)
    lap_ana_ext.annule(0, nz-2) ; 
    
    lap_ana = lap_ana + lap_ana_ext ; 
    (lap_ana.va).c->dzpuis = 0 ;  
    lap_ana.set_dzpuis(0) ; 
    lap_ana.std_base_scal() ; // Sets the standard base for spectral expansions
    
    cout << "Coef de lap_ana (dzpuis=0) : " << endl ; 
    lap_ana.affiche_seuil(cout) ; 

    cout << "diffrel lap <-> lap_ana : " << endl ; 
    cout << diffrel(lap, lap_ana) << endl ; 
    
    cout << "diffrelmax lap <-> lap_ana : " << endl ; 
    cout << diffrelmax(lap, lap_ana) << endl ; 
    
    
    exit(-1) ;   
        
}

