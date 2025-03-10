/*
 * Construction of initial conditions for a NS-BH binary computation.
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
 * Revision 1.12  2016/12/05 16:18:22  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.11  2014/10/13 08:53:53  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.10  2014/10/06 15:09:42  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.9  2008/09/26 08:44:04  p_grandclement
 * Mixted binaries with non vanishing spin
 *
 * Revision 1.6  2006/04/27 09:12:34  p_grandclement
 * First try at irrotational black holes
 *
 * Revision 1.5  2006/04/25 07:22:00  p_grandclement
 * Various changes for the NS_BH project
 *
 * Revision 1.4  2005/08/29 15:10:19  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.3  2004/03/25 12:35:35  j_novak
 * now using namespace Unites
 *
 * Revision 1.2  2002/12/19 14:57:15  e_gourgoulhon
 * New prototype for Bin_ns_bh::set_omega and set_x_axe.
 * The BH is now of the left (xa <0) and the NS on the right (xa>0).
 *
 * Revision 1.1  2002/12/18 10:33:10  e_gourgoulhon
 * Computations of NS - BH binaries
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
#include <cmath>

// Lorene headers
#include "bin_ns_bh.h"
#include "nbr_spx.h"
#include "eos.h"
#include "unites.h"
#include "graphique.h"

using namespace Lorene ;

int main(){

  using namespace Unites ;
    // Identification of all the subroutines called by the code :

    // system("ident init_bin") ;

    //-----------------------------------------------------------------------
    //		Input data for the multi-grid no. 1
    //-----------------------------------------------------------------------

    int nt, np, nz ;
    char blabla[80] ;

    ifstream fich("par_grid_ns.d") ;
    fich.getline(blabla, 80);
    fich.getline(blabla, 80);
    fich >> nz; fich.getline(blabla, 80) ;
    int nzet ;
    fich >> nzet; fich.getline(blabla, 80) ;
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

    //-----------------------------------------------------------------------
    //		Construction of multi-grid  and mapping for the NS
    //-----------------------------------------------------------------------

    Mg3d mg_ns(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_et mp_ns(mg_ns, bornes) ;

    delete [] nr ;
    delete [] nt_tab ;
    delete [] np_tab ;
    delete [] type_r ;
    delete [] bornes ;

    //-----------------------------------------------------------------------
    //		Input data for the BH multi-grid
    //-----------------------------------------------------------------------


    fich.open("par_grid_bh.d") ;
    fich.getline(blabla, 80);
    fich.getline(blabla, 80);
    fich >> nz; fich.getline(blabla, 80) ;
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

    //-----------------------------------------------------------------------
    //		Construction of multi-grid  and mapping for the BH
    //----------------------------------------------------------------------- 
    
    Mg3d mg_bh(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_af mp_bh(mg_bh, bornes) ;

    delete [] nr ;
    delete [] nt_tab ;
    delete [] np_tab ;
    delete [] type_r ;
    delete [] bornes ;

    cout << endl << "NS multi-grid : "
	 << endl << "==============   " << endl << mg_ns << endl ;
    cout << "NS mapping  : "
	 << endl << "===========   " << endl << mp_ns << endl ;

    cout << endl << "BH multi-grid : "
	 << endl << "=============   " << endl << mg_bh << endl ;
    cout << "BH mapping : "
	 << endl << "===========   " << endl << mp_bh << endl ;

    //-----------------------------------------------------------------------
    //		Equation of state for the star
    //-----------------------------------------------------------------------

    fich.open("par_eos.d") ;

    Eos* peos = Eos::eos_from_file(fich) ;
    Eos& eos = *peos ;

    fich.close() ;

    cout << endl << "Equation of state for the star : "
	 << endl << "==============================   " << endl << eos << endl ;

    //-----------------------------------------------------------------------
    //		Physical parameters imput
    //-----------------------------------------------------------------------

    fich.open("par_init.d") ;
    fich.getline(blabla, 80) ;
    fich.getline(blabla, 80) ;

    double separ ;
    fich >> separ; fich.getline(blabla, 80) ;
    separ *= km ;	// translation in Lorene units

    int irrot_i ;
    int state_rot_bh ;
    double ent_c ;
    fich >> ent_c ; fich.getline(blabla, 80) ;
    fich >> irrot_i ; fich.getline(blabla, 80) ;
    fich >> state_rot_bh ; fich.getline(blabla, 80) ;
    bool irrot_ns = (irrot_i == 1) ;
    double lapse_hori ;
    fich >> lapse_hori ; fich.getline(blabla, 80) ;
    
    fich.close() ;

    cout << endl << "Requested orbital separation : " << separ / km
	 << " km" << endl ;

    //-----------------------------------------------------------------------
    //		Construction of the binary system
    //-----------------------------------------------------------------------

    Bin_ns_bh bibi(mp_ns, nzet, eos, irrot_ns, mp_bh) ;
    bibi.set_bh().set_rot_state(state_rot_bh) ;
    bibi.set_bh().set_lapse_hori(lapse_hori) ;
    
    //-----------------------------------------------------------------------
    //		Computation of two static configurations
    //-----------------------------------------------------------------------

    double precis = 1.e-6 ;

    cout << endl << "Computation of a static configuration for the star"
	 << endl << "=================================================" << endl ;

    (bibi.set_ns()).equilibrium_spher(ent_c, precis) ;
    
  cout << endl << "Computation of a static configuration for the black hole"
	 << endl << "=======================================================" << endl ;

    (bibi.set_bh()).init_bhole_seul() ;      // Initialization to some kind of Schwarzschild


    //-----------------------------------------------------------------------
    //		Sets the stars at Newtonian (Keplerian) position
    //-----------------------------------------------------------------------

    // Omega and rotation axis
    // -----------------------

    double mass_ns =  (bibi.get_ns()).mass_g() ;
    double mass_bh =  (bibi.get_bh()).masse_adm_seul() / ggrav ;
    double total_mass =  mass_ns + mass_bh  ;

    cout << "NS mass : " << mass_ns / msol << " M_sol" << endl ;
    cout << "BH mass : " << mass_bh / msol << " M_sol" << endl ;
    bibi.set_x_axe(0.) ;

    // Position of the two objects
    // ---------------------------
    
    double xa_ns = mass_bh / total_mass * separ ;
    ((bibi.set_ns()).set_mp()).set_ori(xa_ns, 0., 0.) ;

    double xa_bh = - mass_ns /  total_mass * separ ;
    ((bibi.set_bh()).set_mp()).set_ori(xa_bh, 0., 0.) ;
    
    // Orientation of the two stars
    // ----------------------------

    // NS aligned with the absolute frame :
    ((bibi.set_ns()).set_mp()).set_rot_phi(M_PI) ;

    // BH anti-aligned with the absolute frame :
    ((bibi.set_bh()).set_mp()).set_rot_phi(0) ;

    // On decouple  
    
    bibi.init_auto() ;   
    int ite ;
    bibi.pseudo_misner (ite, 200, 0.7, precis, 0, lapse_hori) ;
 
  
    cout << endl
    << "=============================================================" << endl
    << "=============================================================" << endl ;
    cout << endl << "Final characteristics of the computed system : " << endl ;
    cout.precision(16) ;
    cout << bibi << endl ;
   
    //-----------------------------------------------------------------------
    //		The result is written in a file
    //-----------------------------------------------------------------------
    FILE* fresu = fopen("statiques.dat", "w") ;

    mg_ns.sauve(fresu) ;
    mp_ns.sauve(fresu) ;
    eos.sauve(fresu) ;

    mg_bh.sauve(fresu) ;
    mp_bh.sauve(fresu) ;

    bibi.sauve(fresu) ;

    fclose(fresu) ;

    // Cleaning
    // --------

    delete peos ;

    return EXIT_SUCCESS ;

}

