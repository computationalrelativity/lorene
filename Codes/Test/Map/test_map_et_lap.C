/*
 * Test program the Laplacian operator of the Cmp class with Map_et mapping
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
 * Revision 1.1  2000/02/25  09:34:18  eric
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
    
    system("ident test_map_et_lap") ; 

    //-----------------------------------------------------------------------
    //		Input data : number of points, types of sampling, etc...
    //-----------------------------------------------------------------------

    int nt = 21 ;    // Number of points in theta
    int np = 20 ;    // Number of points in phi
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
    double csa = 0.7 ;	    // axis ratio c/a 

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
    

    //-----------------------------------------------------------------------
    //		Construction of a Cmp
    //-----------------------------------------------------------------------

    Cmp aa(mp) ; 
    Cmp aa_ext(mp) ;     
    Cmp lap_ana(mp) ; 
    Cmp lap_ana_ext(mp) ;     
    
    const Coord& r = mp.r ; 
    const Coord& x = mp.x ; 
    const Coord& y = mp.y ; 
    const Coord& z = mp.z ; 
    const Coord& cost = mp.cost ; 
    const Coord& sint = mp.sint ; 
    const Coord& cosp = mp.cosp ; 
    const Coord& sinp = mp.sinp ; 

    aa = 0.3 * y - x + (x + z*z) * ( x*x + x*y + pow(z, 6) ) ;

    lap_ana = 2*(3*x + x*x + y + x*y + z*z + 15*x*pow(+z, 4) 
		    + 28*pow(+z,6) ) ;
    
    int nzm1 = nz-1 ; 
    aa.annule(nzm1) ; 
    lap_ana.annule(nzm1) ; 
    
    Mtbl x1 = sint * cosp ; 
    Mtbl y1 = sint * sinp ; 
    Mtbl z1 = cost ; 
    
    aa_ext =  1 / pow(+r, 8) + (x1*x1 - 3*x1*y1 + z1*z1) / pow(+r, 6)
    	  + (x1*z1*z1) / pow(+r, 5) ;
    
    lap_ana_ext = 2. / pow(+r, 3) * ( 28 / pow(+r, 3) 
		+ ( 14 *x1*x1 - 36 *x1*y1 + 2*y1*y1 + 14*z1*z1 ) / pow(+r, 1)
		+ ( x1*x1*x1 + x1*y1*y1 + 5*x1*z1*z1 ) ) ;
    
    aa_ext.annule(0, nzm1-1) ; 
    lap_ana_ext.annule(0, nzm1-1) ; 

    aa = aa + aa_ext ; 
    lap_ana = lap_ana + lap_ana_ext ; 
    
    aa.std_base_scal() ; 
    lap_ana.std_base_scal() ; 
    
    lap_ana.set_dzpuis(4) ; 
    
    Cmp lap = aa.laplacien(4) ; 
    
    cout << "Diffrel lap <-> lap_ana : " << endl ; 
    cout << diffrel(lap, lap_ana) << endl ; 
    
    cout << "Diffrelmax lap <-> lap_ana : " << endl ; 
    cout << diffrelmax(lap, lap_ana) << endl ; 
    
    exit(-1) ;   
        
}

