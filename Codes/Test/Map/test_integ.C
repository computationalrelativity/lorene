/*
 * Test program the integrale function of the Cmp class (Map_et mapping)
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
 * Revision 1.5  2016/12/05 16:18:27  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2015/08/10 15:32:27  j_novak
 * Better calls to Param::add_int(), to avoid weird problems (e.g. with g++ 4.8).
 *
 * Revision 1.3  2014/10/06 15:12:52  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/01/09 11:07:52  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  2000/02/25  13:04:40  eric
 * fact_echelle est desormais passe en double_mod dans le Param pour adapt.
 *
 * Revision 1.1  2000/01/17  13:59:10  eric
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
#include "cmp.h"
#include "utilitaires.h"
#include "param.h"
#include "graphique.h"
#include "nbr_spx.h"

//******************************************************************************

void main(){
    
    // Identification of all the subroutines called by the code : 
    
    system("ident test_integ") ; 

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
	nr[l] = 49 ;		// Number of points in r
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
    
    const Mg3d mg(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    //----------------------------
    // Construction of the Map_et
    //----------------------------
    
    Map_et mp(mg, bornes) ;
        
//    const Coord& r = mp.r ; 
    const Coord& x = mp.x ; 
    const Coord& y = mp.y ; 
    const Coord& z = mp.z ; 
//    const Coord& cost = mp.cost ; 
//    const Coord& sint = mp.sint ; 
//    const Coord& cosp = mp.cosp ; 
//    const Coord& sinp = mp.sinp ; 

    //-----------------------------------------------------------------------
    //		Construction of a Cmp to define a triaxial ellipsoid
    //-----------------------------------------------------------------------

    Cmp ent(mp) ; 
    
    double bsa = 0.7 ; 
    double csa = 0.5 ; 
    double xtri = - 0. ;    // xtri = 0  <=> ellipsoid
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
    // Map_et mp_prev = mp ; 

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
    
    int i_one = 1 ;

    Param par ; 
    par.add_int(nitermax, 0) ;  
    par.add_int(nzet, 1) ;    
    par.add_int(i_one, 2) ; // 1 = full computation
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
    

//    cout << endl << "Coefficients of the mapping functions F and G" << endl ; 
//    des_map_et(mp, 0) ;      
//    des_map_et(mp, 1) ;      

  

    //-----------------------------------------------------------------------
    //		Construction of a Cmp
    //-----------------------------------------------------------------------

    Cmp aa(mp) ; 
//    Cmp aa_ext(mp) ;     
    
    
    aa = 1 ;   
    aa.annule(nzet, nz-1) ;	
    
//    aa_ext = 1 / (r*r) ;    // aa_ext = 1/r^2 in all space
//    aa_ext.annule(0, nz-2) ;	// aa_ext = 0 in intern. dom., =1/r^2 in ext. dom.
    
//    aa = aa + aa_ext ;	// aa = x*y*z*z + x*x in internal domains, 
			// aa = 1/r^2 in external domain 

    aa.set_dzpuis(4) ; 
    
    aa.std_base_scal() ; // Sets the standard basis for spectral expansions

 
    des_coupe_y(aa, 0., -3., 3., -3., 3., 
		"Field the integral of which is to be computed (y=0)") ; 
 
 
    cout << "Coef of aa : " << endl ; 
    aa.affiche_seuil(cout) ; 
    
    cout << endl << "Integral of aa in each domain : " << endl ; 
    cout << aa.integrale_domains() << endl ;  
    
    cout << "Integral of aa in all space : "  << aa.integrale() << endl ;  
    
    cout << endl << "Once more : Integral of aa in each domain : " << endl ; 
    cout << aa.integrale_domains() << endl ;  
    
    cout << "Integral of aa in all space : "  << aa.integrale() << endl ;  
    
    double int_ana = double(4)/double(3) * M_PI * a * b * c ; 
        
    cout << endl <<
	"Analytical value of the integral (volume of the ellipsoid): " 
	 << int_ana << endl ; 
    cout << "Difference computed/analytical : " << aa.integrale()/int_ana -1
	 << endl ; 
    
    // Cleaning
    // --------
    
    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ; 
    delete [] bornes ; 
    delete [] type_r ; 
    
    exit(EXIT_SUCCESS) ; 
    
}

