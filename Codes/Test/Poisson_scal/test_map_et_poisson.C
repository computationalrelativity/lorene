/*
 * Test program for the resolution of Poisson equation
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

char test_map_et_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/01/09 11:07:55  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2000/02/09  15:44:50  eric
 * Modif affichage.
 *
 * Revision 1.2  2000/02/09  15:35:18  eric
 * Suppression de l'appel a rho.set_dzpuis(4) dans le cas d'une
 * source a support compact.
 *
 * Revision 1.1  2000/01/05  10:11:19  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// headers C
#include <stdlib.h>
#include <math.h>

// headers Lorene
#include "type_parite.h"
#include "cmp.h"
#include "utilitaires.h"
#include "graphique.h"
#include "param.h"
#include "nbr_spx.h"

//******************************************************************************

void main(){
    
    // Identification of all the subroutines called by the code : 
    
    system("ident test_map_et_poisson") ; 


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

    double bsa = 0.9 ;	    // axis ratio b/a 
    double csa = 0.8 ;	    // axis ratio c/a 

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
    
    for (int l=0; l<nz; l++) {
//    	des_map_et(mp, l) ;
    } 

    //-----------------------------------------------------------------------
    //		Construction of a Cmp
    //-----------------------------------------------------------------------

    Cmp uu(mp) ; 
    Cmp rho(mp) ; 
    Cmp rho_ext(mp) ; 
    
    const Coord& r = mp.r ; 
//    const Coord& x = mp.x ; 
//    const Coord& y = mp.y ; 
//    const Coord& z = mp.z ; 
//    const Coord& cost = mp.cost ; 
//    const Coord& sint = mp.sint ; 
//    const Coord& cosp = mp.cosp ; 
//    const Coord& sinp = mp.sinp ; 

    cout << "================================" << endl ; 
    cout << "  Source with compact support" << endl ; 
    cout << "================================" << endl ; 

    rho = 1 ; 
    rho.annule(nz-1) ; 

    rho.std_base_scal() ; // Sets the standard basis for spectral expansions
   
    
		// Resolution of Poisson equation 
		//   source with compact support
		// ------------------------------
    
    int nitermax = 100 ; 
    int niter ; 
    double lambda = 0.5 ; 
    double precis = 1.e-14 ; 
    Cmp ssj(mp) ;
    
    Param par ; 
    par.add_int(nitermax) ; 
    par.add_int_mod(niter) ; 
    par.add_double(lambda) ; 
    par.add_double(precis, 1) ; 
    par.add_cmp_mod(ssj) ; 
    
    // Initialization for Map_et::poisson  
    ssj = 0 ;     
    uu = 0 ; 
    ssj.std_base_scal() ; 
    uu.std_base_scal() ; 
    
    rho.poisson(par, uu) ;  
    
    
		// Comparison with the analytical solution 
		// ---------------------------------------
    
    Cmp uu_ana(mp) ; 
    Cmp uu_ana_ext(mp) ; 
    
    Mtbl rr = r ; 
    double ray = rr(nz-1, 0, 0, 0) ; 
    cout << "ray = " << ray << endl ; 
    
    uu_ana = r*r/6 - ray*ray/2 ; 
    uu_ana_ext = - pow(ray, 3) / (3*r) ; 
    
    uu_ana.annule(nz-1) ; 
    uu_ana_ext.annule(0, nz-2) ; 
    uu_ana = uu_ana + uu_ana_ext ; 
    
    cout << "diffrel uu <-> uu_ana : " << endl ; 
    cout << diffrel(uu, uu_ana) << endl ; 
    
    cout << "diffrelmax uu <-> uu_ana : " << endl ; 
    cout << diffrelmax(uu, uu_ana) << endl ; 

    arrete() ; 
    

		// Resolution of Poisson equation 
		//   source with non-compact support
		// ----------------------------------
    
    cout << "================================" << endl ; 
    cout << "  Source with extended support" << endl ; 
    cout << "================================" << endl ; 
    
    rho = 1  ;		
    rho.annule(nz-1) ;	
			
    rho_ext = pow(ray, 4) ;   
    rho_ext.annule(0, nz-2) ;	
    
    rho = rho + rho_ext ;	
    rho.set_dzpuis(4) ; 
   
    rho.std_base_scal() ; // Sets the standard basis for spectral expansions

    // Initialization for Map_et::poisson  
    ssj = 0 ;     
    uu = 0 ; 
    ssj.std_base_scal() ; 
    uu.std_base_scal() ; 
    
    rho.poisson(par, uu) ;  
    
        
		// Comparison with the analytical solution 
		// ---------------------------------------
    
    uu_ana = r*r/6 - ray*ray ; 
    uu_ana_ext = - 4*pow(ray, 3)/(3*r) + pow(ray, 4) / (2*r*r) ; 
    
    uu_ana.annule(nz-1) ; 
    uu_ana_ext.annule(0, nz-2) ; 
    uu_ana = uu_ana + uu_ana_ext ; 
    
    cout << "diffrel uu <-> uu_ana : " << endl ; 
    cout << diffrel(uu, uu_ana) << endl ; 
    
    cout << "diffrelmax uu <-> uu_ana : " << endl ; 
    cout << diffrelmax(uu, uu_ana) << endl ; 

    arrete() ; 

		// Comparison of the Laplacian of the solution 
		//  with the source
		// -------------------------------------------


    cout << endl << 
    "The Laplacian of the solution is computed and compared with the source:"
	<< endl ; 
	
    Cmp lap = uu.laplacien() ; 
    
    cout << "max( |lap(uu) - rho| ) " << max(abs(lap - rho)) << endl ; 
    cout << "norme( lap(uu) - rho)/norme(rho) " 
	  << norme( lap - rho)/norme(rho) << endl ; 

    
    exit(EXIT_SUCCESS) ; 
            
}

