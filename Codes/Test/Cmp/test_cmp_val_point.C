/*
 *  Test code for the function Cmp::val_point
 *
 */

/*
 *   Copyright (c) 2002 Eric Gourgoulhon
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

char test_cmp_val_point_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/06 15:12:52  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2002/10/16 14:37:16  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/05/11 12:40:55  e_gourgoulhon
 * Added test for the basis T_COSSIN_SI.
 *
 * Revision 1.1  2002/05/05 16:26:47  e_gourgoulhon
 * Initial commit
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>

// Lorene headers
#include "cmp.h"
#include "nbr_spx.h"

int main() {

    //-----------------------------------------------------------------------
    //		Input data : number of points, types of sampling, etc...
    //-----------------------------------------------------------------------

    int nz = 3 ;    // Number of domains

    int nr = 5 ;    // Number of points in r (in each domain)
    int nt = 7 ;    // Number of points in theta (in each domain)
    int np = 8 ;    // Number of points in phi (in each domain)

    // Type of sampling in theta:
    //    SYM : theta in [0,pi/2]  (symmetry with respect equatorial plane)
    //	  NONSYM : theta in [0,pi] (no symmetry)
    int type_t = SYM ;

    // Type of sampling in phi:
    //    SYM : phi in [0,pi[  (symmetry phi --> phi + pi)
    //	  NONSYM : phi in [0,2 pi[ (no symmetry)
    int type_p = NONSYM ;

    // Shall the outermost domain be compactified ?
    bool compact = true ;

    // Boundaries of the domains:
    double* bornes = new double[nz+1];

    for (int l=0; l<nz; l++) {
	bornes[l] = 2*l ;
    }
    if (compact) {
    	bornes[nz] = __infinity ;
    }
    else {
    	bornes[nz] = 2*nz ;
    }

    //-----------------------------------------------------------------------
    //		Construction of a multi-grid
    //-----------------------------------------------------------------------

    const Mg3d mg(nz, nr, nt, np, type_t, type_p, compact) ;

    //-----------------------------------------------------------------------
    //		Construction of a mapping
    //-----------------------------------------------------------------------

    const Map_af mp(mg, bornes) ;

    const Coord& r = mp.r ;
    const Coord& tet = mp.tet ;
    const Coord& phi = mp.phi ;
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

    Cmp aa(mp) ;
    Cmp aa_aux(mp) ;   // auxilliary Cmp to construct aa

    /*
     *********** Test of R_CHEBPIM_I and T_COSSIN_SP  ************

    aa_aux = x + z*z + x*y  ;
    aa_aux.annule(nz-1) ;	

    // Spectral expansion bases :
    aa_aux.std_base_scal() ;   // standard basis for a symmetric scalar field

    aa = aa_aux.srdsdt() ;

    *
    ***************************************************************
    */

    /*
     *********** Test of R_CHEBPIM_P and T_COSSIN_SI **************
     */
    aa_aux =  ( x + z*z + x*y ) * z * x ;
    aa_aux.std_base_scal() ;
    aa_aux.va.set_base_r(0, R_CHEBPIM_I) ;
    aa_aux.va.set_base_t(T_COSSIN_CI) ;

    aa = aa_aux.srdsdt() ;

    /*
     ***************************************************************
     */

    aa.va.coef() ;

    aa.affiche_seuil(cout) ;

    //----------------------------------------------------------------------
    // Call to the val_point function
    //----------------------------------------------------------------------

    double r0, theta0, phi0 ;
    do {
  	cout << endl << "r ? (negative value to quit) " << endl ;
    	cin >> r0 ;
    	if (r0 >=0 ) {
    		cout << "theta ?" << endl ;
    		cin >> theta0 ;
    		cout << "phi ?" << endl ;
    		cin >> phi0 ;

    		double aa0 = aa.val_point(r0, theta0, phi0) ;
    		cout << "value at the point (" << r0 << "," << theta0 << ","
    			<< phi0 << ") : " << aa0 << endl ;
    	}
    }
    while (r0 >= 0.) ;


    // Test on every grid point
    // ------------------------
    Cmp bb(mp) ;
    bb.allocate_all() ;
    Mtbl mr = r ;
    Mtbl mtet = tet ;
    Mtbl mphi = phi ;
    for (int l=0; l<nz; l++) {
    	for (int k=0; k<np; k++) {
    		for (int j=0; j<nt; j++) {
    			for (int i=0; i<nr; i++) {
    				double r1 = mr(l,k,j,i) ;
    				double t1 = mtet(l,k,j,i) ;
    				double p1 = mphi(l,k,j,i) ;
    				bb.set(l,k,j,i) = aa.val_point(r1,t1,p1) ;
    			}
    		}
    	}
   }
   bb.annule(nz-1) ;	

   Cmp diff = aa - bb ;
   cout << "Test at every grid point:" << endl ;
   // cout << diff << endl ;
   diff.affiche_seuil(cout, 1) ;
   cout << endl << "Maximum value of the difference in each domain: " << endl ;
   cout << max(diff) << endl ;

    //----------------------------------------------------------------------
    // clean exit
    //----------------------------------------------------------------------

    delete [] bornes ;

    return EXIT_SUCCESS ;
}

