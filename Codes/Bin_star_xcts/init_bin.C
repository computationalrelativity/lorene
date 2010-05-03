/*
 * Construction of initial conditions for a binary star computation
 * in XCTS formalism (see star_xcts.h and binary_xcts.h for details)
 * 
 */

/*
 *   Copyright (c) 2010 Michal Bejger 
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

char init_bin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2010/05/03 13:39:14  m_bejger
 * File handler fixed
 *
 * Revision 1.1  2010/04/29 15:05:17  m_bejger
 * Initial version
 *
 *
 * $Header$
 *
 */
 
// headers C
#include <stdlib.h>
#include <math.h>

// headers Lorene
#include "unites.h"
#include "binary_xcts.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"

int  main(){
    
    //------------------------------------------------------------------
    //		Input data for the multi-grid no. 1
    //------------------------------------------------------------------

    int nt, np, nz ;
    char blabla[80] ;

    ifstream fich("par_grid1.d") ;
    fich.getline(blabla, 80);
    fich.getline(blabla, 80);
    fich >> nz; fich.getline(blabla, 80) ;
    int nzet1 ; 
    fich >> nzet1; fich.getline(blabla, 80) ;
    fich >> nt; fich.getline(blabla, 80) ;
    fich >> np; fich.getline(blabla, 80) ;

    cout << "total number of domains :   nz = " << nz << endl ;
    cout << "number of points in phi :   np = " << np << endl ;
    cout << "number of points in theta : nt = " << nt << endl ;

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
     
    fich.getline(blabla, 80);
    for (int l=0; l<nz; l++) {
	fich >> nr[l]; 
	fich >> bornes[l]; fich.getline(blabla, 80) ;
	np_tab[l] = np ; 
	nt_tab[l] = nt ; 
    }
    bornes[nz] = __infinity ; 

    fich.close();


    // Type of r sampling : 
    int* type_r = new int[nz];
    type_r[0] = RARE ; 
    for (int l=1; l<nz-1; l++) {
	type_r[l] = FIN ; 
    }
    type_r[nz-1] = UNSURR ; 
    
    // Type of sampling in theta and phi :
    int type_t = SYM ; 
    int type_p = NONSYM ; 
    
    //------------------------------------------------------------------
    //		Construction of multi-grid 1 and mapping 1
    //------------------------------------------------------------------
   
    Mg3d mg1(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_et mp1(mg1, bornes) ;
   
    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ; 
    delete [] type_r ; 
    delete [] bornes ; 
       
    //------------------------------------------------------------------
    //		Input data for the multi-grid no. 2
    //------------------------------------------------------------------


    fich.open("par_grid2.d") ;
    fich.getline(blabla, 80);
    fich.getline(blabla, 80);
    fich >> nz; fich.getline(blabla, 80) ;
    int nzet2 ; 
    fich >> nzet2; fich.getline(blabla, 80) ;
    fich >> nt; fich.getline(blabla, 80) ;
    fich >> np; fich.getline(blabla, 80) ;

    cout << "total number of domains :   nz = " << nz << endl ;
    cout << "number of points in phi :   np = " << np << endl ;
    cout << "number of points in theta : nt = " << nt << endl ;

    nr = new int[nz];
    nt_tab = new int[nz];
    np_tab = new int[nz];
    bornes = new double[nz+1];
     
    fich.getline(blabla, 80);
    for (int l=0; l<nz; l++) {
	fich >> nr[l]; 
	fich >> bornes[l]; fich.getline(blabla, 80) ;
	np_tab[l] = np ; 
	nt_tab[l] = nt ; 
    }
    bornes[nz] = __infinity ; 

    fich.close();


    // Type of r sampling : 
    type_r = new int[nz];
    type_r[0] = RARE ; 
    for (int l=1; l<nz-1; l++) {
	type_r[l] = FIN ; 
    }
    type_r[nz-1] = UNSURR ; 
        
    //------------------------------------------------------------------
    //		Construction of multi-grid 2 and mapping 2
    //------------------------------------------------------------------
   
    Mg3d mg2(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_et mp2(mg2, bornes) ;
   
    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ; 
    delete [] type_r ; 
    delete [] bornes ; 
       
    cout << endl << "Multi-grid 1 : " 
	 << endl << "============   " << endl << mg1 << endl ; 
    cout << "Mapping 1 : " 
	 << endl << "=========   " << endl << mp1 << endl ; 

    cout << endl << "Multi-grid 2 : " 
	 << endl << "============   " << endl << mg2 << endl ; 
    cout << "Mapping 2 : " 
	 << endl << "=========   " << endl << mp2 << endl ; 

    //------------------------------------------------------------------
    //		Equation of state for star 1
    //------------------------------------------------------------------

    fich.open("par_eos1.d") ;

    Eos* peos1 = Eos::eos_from_file(fich) ; 
    Eos& eos1 = *peos1 ; 

    fich.close() ; 

    //------------------------------------------------------------------
    //		Equation of state for star 2
    //------------------------------------------------------------------

    fich.open("par_eos2.d") ;

    Eos* peos2 = Eos::eos_from_file(fich) ; 
    Eos& eos2 = *peos2 ; 

    fich.close() ; 

    cout << endl << "Equation of state of star 1 : " 
	 << endl << "===========================   " << endl << eos1 << endl ; 
    cout << endl << "Equation of state of star 2 : " 
	 << endl << "===========================   " << endl << eos2 << endl ; 


    //------------------------------------------------------------------
    //		Physical parameters imput
    //------------------------------------------------------------------

    using namespace Unites ;

    fich.open("par_init.d") ;
    fich.getline(blabla, 80) ;
    fich.getline(blabla, 80) ;

    double separ ; 
    fich >> separ; fich.getline(blabla, 80) ;
    separ *= km ;	// translation in machine units

    int irrot1_i, irrot2_i ; 
    double ent_c1, ent_c2 ; 
    fich >> ent_c1 ; fich.getline(blabla, 80) ;
    fich >> irrot1_i ; fich.getline(blabla, 80) ;
    bool irrot1 = (irrot1_i == 1) ; 
    fich >> ent_c2 ; fich.getline(blabla, 80) ;
    fich >> irrot2_i ; fich.getline(blabla, 80) ;
    bool irrot2 = (irrot2_i == 1) ; 
    
    fich.close() ; 
  
    cout << endl << "Requested orbital separation : " << separ / km 
	 << " km" << endl ; 

    //------------------------------------------------------------------
    //		Construction of a binary system
    //------------------------------------------------------------------

    Binary_xcts star(mp1, nzet1, eos1, irrot1, 
		             mp2, nzet2, eos2, irrot2) ;			
    
    //------------------------------------------------------------------
    //		Computation of two static configurations
    //------------------------------------------------------------------

    double precis = 1.e-12 ; 
    
    cout << endl << "Computation of a static configuration for star 1"
	 << endl << "================================================" << endl ;
    (star.set(1)).equilibrium_spher(ent_c1, precis) ; 

    (star.set(1)).set_Psi_auto() = exp(0.5*(star(1).get_lnq() 
    							 - star(1).get_logn()));
    (star.set(1)).set_Psi_auto().std_spectral_base() ;

    (star.set(1)).set_chi_auto() = exp(0.5*(star(1).get_lnq() 
    							 + star(1).get_logn())) ;
    (star.set(1)).set_chi_auto().std_spectral_base() ; 
    
    cout << endl << "Computation of a static configuration for star 2"
	 << endl << "================================================" << endl ; 

    (star.set(2)).equilibrium_spher(ent_c2, precis) ; 

    (star.set(2)).set_Psi_auto() = exp(0.5*(star(2).get_lnq() 
    							 - star(2).get_logn())) ;
    (star.set(2)).set_Psi_auto().std_spectral_base() ;

    (star.set(2)).set_chi_auto() = exp(0.5*(star(2).get_lnq() 
    							 + star(2).get_logn())) ;
    (star.set(2)).set_chi_auto().std_spectral_base() ;    

    //-----------------------------------------------------------------------
    //		Sets the stars at Newtonian (Keplerian) position 
    //-----------------------------------------------------------------------

    // Omega and rotation axis
    // -----------------------

    double total_mass =  star(1).mass_g() + star(2).mass_g() ; 

    star.set_omega() = sqrt( g_si/g_unit * total_mass / pow(separ, 3.) ) ;
 
    star.set_x_axe() = 0 ; 


    // Position of the two stars
    // -------------------------
    
    for (int i=1 ; i<=2 ; i++) {

	double xa_et = (star(3-i).mass_g()) / total_mass * separ ;
	if (i == 1) xa_et = - xa_et ; 

	((star.set(i)).set_mp()).set_ori(xa_et, 0., 0.) ; 	
    }

    // Orientation of the two stars
    // ----------------------------
    
    // Star 1 aligned with the absolute frame : 
    ((star.set(1)).set_mp()).set_rot_phi(0.) ; 
    
    // Star 2 anti-aligned with the absolute frame : 
    ((star.set(2)).set_mp()).set_rot_phi(M_PI) ; 
    
        
    cout << endl 
    << "=============================================================" << endl 
    << "=============================================================" << endl ;
    cout << endl << "Final characteristics of the computed system : " << endl ; 
    cout.precision(16) ; 
    cout << star << endl ; 
    
    //-----------------------------------------------------------------------
    //		The result is written in a file
    //-----------------------------------------------------------------------

    FILE* fresu = fopen("ini.d", "w") ; 
    
    int mer = 0 ; 
    fwrite(&mer, sizeof(int), 1, fresu) ;	// mer

    mg1.sauve(fresu) ; 
    mp1.sauve(fresu) ; 
    eos1.sauve(fresu) ; 

    mg2.sauve(fresu) ; 
    mp2.sauve(fresu) ; 
    eos2.sauve(fresu) ; 

    star.sauve(fresu) ;     

    fclose(fresu) ;     

    // Cleaning
    // --------

    delete peos1 ;    
    delete peos2 ;    

    return EXIT_SUCCESS ; 
    
}
