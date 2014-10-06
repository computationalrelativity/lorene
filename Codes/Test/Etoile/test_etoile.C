/*
 * Test program for the Etoile class
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

char test_etoile_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:12:52  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2002/10/16 14:37:18  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  2000/01/28  17:19:14  eric
 * *** empty log message ***
 *
 * Revision 1.1  2000/01/27  16:48:01  eric
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
#include "etoile.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"

//******************************************************************************

void main(){
    

    {	    // Main block

    // Identification of all the subroutines called by the code : 
    
    // system("ident test_etoile > identification.d") ; 

    //-----------------------------------------------------------------------
    //		Input data : number of points, types of sampling, etc...
    //-----------------------------------------------------------------------

    int nt = 7 ;    // Number of points in theta
    int np = 4 ;    // Number of points in phi
    int nz = 3 ;    // Number of domains

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
    int* type_r = new int[nz];

    for (int l=0; l<nz; l++) {
	nr[l] = 17 ; 
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
    
    Mg3d mg(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    cout << endl << "Grid mg : " << mg << endl ; 
        
    // Cleaning
    // --------
    
    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ; 
    delete [] type_r ; 

    //-------------------------------------------
    // Construction of a mapping of class Map_et
    //-------------------------------------------
    
    Map_et mp(mg, bornes) ;
    
    // Cleaning
    // --------

    delete [] bornes ; 
    
//    cout << endl << "Mapping mp (before adaptation) : " << endl ; 
//    cout << mp << endl ; 
    
    const Coord& r = mp.r ; 
//    const Coord& x = mp.x ; 
//    const Coord& y = mp.y ; 
//    const Coord& z = mp.z ; 
//    const Coord& cost = mp.cost ; 
//    const Coord& sint = mp.sint ; 
//    const Coord& cosp = mp.cosp ; 
//    const Coord& sinp = mp.sinp ; 


    //-----------------------------------------------------------------------
    //		Construction of an EOS
    //-----------------------------------------------------------------------

    double gamma = 2. ;
    double kappa = 0.0332 ;
     
    Eos_poly eos0(gamma, kappa) ; 



    //-----------------------------------------------------------------------
    //		Construction of a star
    //-----------------------------------------------------------------------

    int nzet = nz - 2 ;	    // Number of domains occupied by the star
    
    bool relat = true ;	    // true = relativisitic star
    
    Etoile star(mp, nzet, relat, eos0) ;	    
    
    //-----------------------------------------------------------------------
    //		Initialization of the stellar enthalpy field
    //-----------------------------------------------------------------------

    Cmp ent0(mp) ; 

    double ray = mp.val_r(nzet-1, 1, M_PI/2., 0.) ; 
    
    ent0 = 1 - r*r/(ray*ray) ;
    
    ent0.annule(nz-1) ; 
    ent0.std_base_scal() ; 

    star.set_enthalpy(ent0) ; 
    
    cout << star << endl ; 
    
    //-----------------------------------------------------------------------
    //		File output / input
    //-----------------------------------------------------------------------

    FILE* fich = fopen("res.d", "w") ; 
    
    mg.sauve(fich) ; 
    mp.sauve(fich) ; 
    eos0.sauve(fich) ; 
    star.sauve(fich) ; 
    
    fclose(fich) ; 
    
    fich = fopen("res.d", "r") ; 

    Mg3d mg1(fich) ; 
    Map_et mp1(mg1, fich) ; 
    Eos* peos1 = Eos::eos_from_file(fich) ; 
    
    Etoile star1(mp1, *peos1, fich) ; 

    fclose(fich) ; 

    
    if (*peos1 != eos0) {
	cout << "Problem : *peos1 != eos0 !" << endl ; 
	abort() ; 
    }    
    else {
	cout << "*peos1 == eos0 : OK " << endl ; 	
    }

    star1.equation_of_state() ;

//    cout << endl << "star1 : " << star1 << endl ; 


    delete peos1 ; 
    
    //-----------------------------------------------------------------------
    //		Computation of a static equilibrium configuration
    //-----------------------------------------------------------------------

    cout << endl << "Computation of a static equilibrium configuration"
	 << endl << "-------------------------------------------------"
	 << endl ; 

    double ent_c = 0.29920935 ;    // central enthalpy
    double precis = 1.e-12 ; 
    
    star.equilibrium_spher(ent_c, precis) ; 

    cout << "Final mapping : " << endl ;
    cout << mp << endl ; 
    cout << endl << "Final star : " << star << endl ; 
    
    des_profile(star.get_ent()(), 0., 1.2*ray, M_PI/2., 0., "Enthalpy") ; 


    //-----------------------------------------------------------------------
    //		The End
    //-----------------------------------------------------------------------

    } // End of main block : everything should be properly destroyed here
    

    exit(EXIT_SUCCESS) ; 
    
}

