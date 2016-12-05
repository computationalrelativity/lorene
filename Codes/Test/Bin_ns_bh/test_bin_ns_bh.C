/*
 *  Test code for class Bin_ns_bh
 *
 */

/*
 *   Copyright (c) 2002  Philippe Grandclement, Keisuke Taniguchi,
 *              Eric Gourgoulhon
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:18:26  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:54:00  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:12:51  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2002/12/17 13:09:21  e_gourgoulhon
 * Test code for class Bin_ns_bh
 *
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>

// Lorene headers
#include "bin_ns_bh.h"
#include "nbr_spx.h"
#include "eos.h"

using namespace Lorene ;

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
    // Construction of a mapping of class Map_et  for the NS
    //-------------------------------------------

    Map_et mp_ns(mg, bornes) ;

    mp_ns.set_ori(-10.,0.,0.) ;

    //-------------------------------------------
    // Construction of a mapping of class Map_af for the BH
    //-------------------------------------------

    Map_af mp_bh(mg, bornes) ;

    mp_bh.set_ori(10.,0.,0.) ;
    mp_bh.set_rot_phi(M_PI) ;

    // Cleaning
    // --------

    delete [] bornes ;

    cout << "NS mapping : " << mp_ns << endl ;
    cout << "BH mapping : " << mp_bh << endl ;

    //-----------------------------------------------------------------------
    //		Construction of an EOS
    //-----------------------------------------------------------------------

    double gamma = 2. ;
    double kappa = 0.0332 ;

    Eos_poly eos(gamma, kappa) ;


    //-----------------------------------------------------------------------
    //		Construction of a NS-BH binary system
    //-----------------------------------------------------------------------

    int nzet = nz - 2 ;	    // Number of domains occupied by the star
    bool irrot_ns = true ;

    Bin_ns_bh bibi(mp_ns, nzet, eos, irrot_ns, mp_bh) ;


    //-----------------------------------------------------------------------
    //		Initialization of the stellar enthalpy field
    //-----------------------------------------------------------------------

    Cmp ent0(mp_ns) ;

    double ray = mp_ns.val_r(nzet-1, 1, M_PI/2., 0.) ;

    const Coord& r = mp_ns.r ;

    ent0 = 1 - r*r/(ray*ray) ;
    
    ent0.annule(nz-1) ; 
    ent0.std_base_scal() ; 

    (bibi.set_ns()).set_enthalpy(ent0) ;

    cout << bibi << endl ;


    }

  return EXIT_SUCCESS ;

}

