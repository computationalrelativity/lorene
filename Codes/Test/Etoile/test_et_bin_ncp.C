/*
 *  Test code for class Et_bin_ncp
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

char test_et_bin_ncp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/12/09 10:29:38  f_limousin
 * Test code for class Et_bin_ncp
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>

// Lorene headers
#include "et_bin_ncp.h"
#include "nbr_spx.h"
#include "eos.h"

int main() {

    {	    // Main block

    // Identification of all the subroutines called by the code : 
    
    // system("ident test_etoile > identification.d") ; 

    //-----------------------------------------------------------------------
    //		Input data : number of points, types of sampling, etc...
    //-----------------------------------------------------------------------

    int nz = 3 ;    // Number of domains
    int nr = 17 ;    // Number of points in theta
    int nt = 7 ;    // Number of points in theta
    int np = 4 ;    // Number of points in phi

    double* bornes = new double[nz+1];

    bornes[0] = 0. ; 
    bornes[1] = 2. ;
    bornes[2] = 3. ;
    bornes[3] = __infinity ; 
    
    int typ_t = SYM ; 

    int typ_p = NONSYM ; 
    

    //-----------------------------------------------------------------------
    //		Construction of a multi-grid
    //-----------------------------------------------------------------------
    
    Mg3d mg(nz, nr, nt, np, typ_t, typ_p, true) ;   

    cout << endl << "Grid mg : " << mg << endl ; 
        
    //-------------------------------------------
    // Construction of a mapping of class Map_et
    //-------------------------------------------
    
    Map_et mp(mg, bornes) ;
    
    // Cleaning
    // --------

    delete [] bornes ; 
    
    cout << mp << endl ; 
    
    //-----------------------------------------------------------------------
    //		Construction of an EOS
    //-----------------------------------------------------------------------

    double gamma = 2. ;
    double kappa = 0.0332 ;
     
    Eos_poly eos0(gamma, kappa) ; 


    //-----------------------------------------------------------------------
    //         Construction of a flat metric
    //-----------------------------------------------------------------------

    Tenseur_sym plat(mp, 2, COV, mp.get_bvect_cart() ) ;
    plat.set_etat_qcq() ;
    for (int i=0; i<3; i++) {
      for (int j=0; j<i; j++) {
	plat.set(i,j) = 0 ;
      }
      plat.set(i,i) = 1 ;
    }
    plat.set_std_base() ;

    Metrique flat(plat, true) ; 

    //-----------------------------------------------------------------------
    //		Construction of a star
    //-----------------------------------------------------------------------

    int nzet = nz - 2 ;	    // Number of domains occupied by the star
    
    bool relat = true ;	    // true = relativisitic star
    bool irrot = true ; 

    Et_bin_ncp star(mp, nzet, relat, eos0, irrot, mp.get_bvect_cart(),
		    flat) ;

    //-----------------------------------------------------------------------
    //		Initialization of the stellar enthalpy field
    //-----------------------------------------------------------------------

    Cmp ent0(mp) ; 

    double ray = mp.val_r(nzet-1, 1, M_PI/2., 0.) ; 
    
    const Coord& r = mp.r ; 

    ent0 = 1 - r*r/(ray*ray) ;
    
    ent0.annule(nz-1) ; 
    ent0.std_base_scal() ; 

    star.set_enthalpy(ent0) ; 
    
    cout << star << endl ; 
    

    }

  return EXIT_SUCCESS ; 

}

