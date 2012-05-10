/*
 * Main code for computing stationary axisymmetric rotating stars in Dirac  
 * gauge and maximal slicing.
 *
 */

/*
 *   Copyright (c) 2005 Lap-Ming Lin & Jerome Novak
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

char rotstar_dirac_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2012/05/10 09:05:30  j_novak
 * New code examrot_dirac for reading the results of
 * rotstar_dirac. Simplification of the parrot.d parameter file for
 * rotstar_dirac.
 *
 * Revision 1.7  2009/10/26 10:42:30  j_novak
 * Changed dzpuis handling in final test.
 *
 * Revision 1.6  2008/08/22 06:34:13  j_novak
 * Now displaying \hat{A}^{ij} and \hat{A}^{ij}_{TT} (see paper by Cordero,
 * Cerda, Dimmelmeier et al.).
 *
 * Revision 1.5  2008/02/25 10:40:54  j_novak
 * Added the flag mer_hij to control the step from which the equation for h^{ij}
 * is being solved.
 *
 * Revision 1.4  2007/11/06 16:24:00  j_novak
 * Added the flag spectral_filter giving the order of possible spectral filtering
 * of the hydro sources of metric equations (some members *_euler). The filtering
 * is done in strot_dirac_hydro, if this flag is non-zero.
 *
 * Revision 1.3  2005/03/25 13:53:34  j_novak
 * Added drawings of various fields.
 *
 * Revision 1.2  2005/03/10 09:38:35  j_novak
 * *** empty log message ***
 *
 * Revision 1.1  2005/01/31 09:10:06  j_novak
 * Main code for rotating stars in Dirac gauge.
 *
 *
 * $Header$
 *
 */

// headers C
#include <stdlib.h>

// headers Lorene
#include "star_rot_dirac.h"
#include "cmp.h"
#include "eos.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"

// Local prototype (for drawings only)
Cmp raccord_c1(const Cmp& uu, int l1) ; 



//******************************************************************************

int main(){

  using namespace Unites ; 

    // Identification of all the subroutines called by the code : 
    
    system("ident rotstar_dirac > identif.d") ; 

    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    char blabla[120] ;

    int mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
	delta_mer_kep, mer_mass, graph, nz, nzet, nzadapt,
	nt, np, filter_order, mer_hij ; 
    double ent_c, freq_si, fact_omega, mbar_wanted, precis, freq_ini_si, 
	   aexp_mass, relax ;  
    
    ifstream fich("parrot.d") ;
    fich.getline(blabla, 120) ; fich.getline(blabla, 120) ;
    fich >> ent_c ; fich.getline(blabla, 120) ;
    fich >> freq_si ; fich.getline(blabla, 120) ;
    fich >> fact_omega ; fich.getline(blabla, 120) ;
    fich >> mbar_wanted ; fich.getline(blabla, 120) ;
    mbar_wanted *= msol ; 
    fich.getline(blabla, 120) ;
    fich >> mer_max ; fich.getline(blabla, 120) ;
    fich >> precis ; fich.getline(blabla, 120) ;
    fich >> mer_rot ; fich.getline(blabla, 120) ;
    fich >> freq_ini_si ; fich.getline(blabla, 120) ;
    fich >> mer_change_omega ; fich.getline(blabla, 120) ;
    fich >> mer_fix_omega ; fich.getline(blabla, 200) ;
    fich >> delta_mer_kep ; fich.getline(blabla, 120) ;
    fich >> mer_mass ; fich.getline(blabla, 120) ;
    fich >> aexp_mass ; fich.getline(blabla, 120) ;
    fich >> relax ; fich.getline(blabla, 120) ;
    fich >> graph ; fich.getline(blabla, 120) ;
    fich.getline(blabla, 120) ;
    fich >> nz ; fich.getline(blabla, 120) ;
    fich >> nzet; fich.getline(blabla, 120) ;
    fich >> nzadapt; fich.getline(blabla, 120) ;
    fich >> nt; fich.getline(blabla, 120) ;
    fich >> np; fich.getline(blabla, 120) ;

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
     
    fich.getline(blabla, 120);
    for (int l=0; l<nz; l++) {
	fich >> nr[l]; 
	fich >> bornes[l]; fich.getline(blabla, 120) ;
	np_tab[l] = np ; 
	nt_tab[l] = nt ; 
    }
    bornes[nz] = __infinity ;

    Tbl ent_limit(nzet) ;
    ent_limit.set_etat_qcq() ;
    ent_limit.set(nzet-1) = 0 ; 	// enthalpy at the stellar surface
    for (int l=0; l<nzet-1; l++) {
    	fich >> ent_limit.set(l) ; fich.getline(blabla, 120) ;
    }
    fich >> filter_order ; fich.getline(blabla, 120) ;
    fich >> mer_hij ;fich.getline(blabla, 120) ;

    fich.close();

    // Particular cases
    // ----------------

    // Initial frequency = final frequency
    if ( freq_ini_si < 0 ) {
	freq_ini_si = freq_si ; 
	mer_change_omega = mer_rot ; 
	mer_fix_omega = mer_rot + 1 ;  
    }

    
    //-----------------------------------------------------------------------
    //		Equation of state
    //-----------------------------------------------------------------------

    fich.open("par_eos.d") ;

    Eos* peos = Eos::eos_from_file(fich) ;
    Eos& eos = *peos ;

    fich.close() ;

    //-----------------------------------------------------------------------
    //		Construction of the multi-grid and the mapping
    //-----------------------------------------------------------------------

    // Rescale of bornes in the case where there more than 1 domain inside
    //   the star

    for (int l=0; l<nzet-1; l++) {
      
      bornes[l+1] = bornes[nzet] * sqrt(1 - ent_limit(l) / ent_c) ;

    }

    // Type of r sampling :
    int* type_r = new int[nz];
    type_r[0] = RARE ; 
    for (int l=1; l<nz-1; l++) {
	type_r[l] = FIN ; 
    }
    type_r[nz-1] = UNSURR ; 
    
    // Type of sampling in theta and phi :
    int type_t = SYM ; 
    int type_p = SYM ; 
    
    Mg3d mg(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_af mp(mg, bornes) ;
   
    // Cleaning
    // --------

    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ; 
    delete [] type_r ; 
    delete [] bornes ; 
       


    cout << endl 
	 << "==========================================================" << endl
	 << "                    Physical parameters                   " << endl
	 << "=========================================================="
	 << endl ; 
    cout << endl ;

    cout << endl << "Equation of state : " 
	 << endl << "=================   " << endl ;
    cout << eos << endl ; 

    cout << "Central enthalpy : " << ent_c << " c^2" << endl ; 
    cout << "Rotation frequency : " << freq_si << " Hz" << endl ; 
    if ( abs(mer_mass) < mer_max ) {
	cout << "Required Baryon mass [M_sol] : " 
	     << mbar_wanted / msol << endl ; 
    }
    
    cout << endl 
	 << "==========================================================" << endl
	 << "               Computational parameters                   " << endl
	 << "=========================================================="
	 << endl << endl ; 

    cout << "Maximum number of steps in the main iteration : " 
	 << mer_max << endl ; 
    cout << "Relaxation factor in the main iteration  : " 
	 << relax << endl ; 
    cout << "Threshold on the enthalpy relative change for ending the computation : " 
	 << precis << endl ; 
    cout << "Step from which the baryon mass is forced to converge : " 
	 << mer_mass << endl ; 
    cout << "Exponent for the increase factor of the central enthalpy : " 
	 << aexp_mass << endl ; 

    cout << endl << "Multi-grid : " 
	 << endl << "==========" << endl << mg << endl ; 
    cout << "Mapping : " 
	 << endl << "=======" << endl << mp << endl ; 


    //-----------------------------------------------------------------------
    //		Construction of the star
    //-----------------------------------------------------------------------
    
    Star_rot_Dirac star(mp, nzet, eos, filter_order) ; 
    

    //-----------------------------------------------------------------------
    //		Initialization of the enthalpy field
    //-----------------------------------------------------------------------


    const Coord& r = mp.r ;
    double ray0 = mp.val_r(nzet-1, 1., 0., 0.) ;  
    Scalar ent0(mp) ; 
    ent0 = ent_c * ( 1 - r*r / (ray0*ray0) ) ; 
    ent0.annule_domain(nz-1) ; 
    ent0.std_spectral_base() ; 
    star.set_enthalpy(ent0) ;  

    // Initialization of (n,e,p) from H
    star.equation_of_state() ; 

    // Initialization of (E,S,U,etc...) (quantities relative to the Eulerian obs)
    star.hydro_euler() ; 

    cout << endl << "Initial star : " 
	 << endl << "============   " << endl ;

    cout << star.get_gamma() ;

    cout << star << endl ; 
     
    //-----------------------------------------------------------------------
    //		Computation of the rotating equilibrium
    //-----------------------------------------------------------------------

    double omega = 2 * M_PI * freq_si / f_unit ; 
    double omega_ini = 2 * M_PI * freq_ini_si / f_unit ; 

    Itbl icontrol(7) ;
    icontrol.set_etat_qcq() ; 
    icontrol.set(0) = mer_max ; 
    icontrol.set(1) = mer_rot ; 
    icontrol.set(2) = mer_change_omega ; 
    icontrol.set(3) = mer_fix_omega ; 
    icontrol.set(4) = mer_mass ; 
    icontrol.set(5) = delta_mer_kep ; 
    icontrol.set(6) = mer_hij ;
    
    Tbl control(3) ; 
    control.set_etat_qcq() ; 
    control.set(0) = precis ; 
    control.set(1) = omega_ini ; 
    control.set(2) = relax ; 

    Tbl diff(8) ;     

    star.equilibrium(ent_c, omega, fact_omega, nzadapt, ent_limit, icontrol, control,
		     mbar_wanted, aexp_mass, diff) ;

     
    cout << endl << "Final star : " 
	 << endl << "==========   " << endl ;

    cout.precision(10) ; 
    cout << star << endl ;

    double rho_c = star.get_ener().val_grid_point(0,0,0,0) ;

    cout << "r_p/r_eq :" << star.aplat() << endl ;
    cout << "Omega rho0^{-1/2} : " << star.get_omega() /
      sqrt( ggrav * rho_c ) << endl ;

    cout << "M rho0^{1/2} : " << star.mass_g() * pow(ggrav,1.5) *
      sqrt( rho_c ) << endl ;
    
    cout << "M_B rho0^{1/2} : " << star.mass_b() * pow(ggrav,1.5) *
      sqrt( rho_c ) << endl ;
    
    cout << "R_circ rho0^{1/2} : " << star.r_circ() *
      sqrt( ggrav * rho_c ) << endl ;
    
    cout << "J rho0 : " << star.angu_mom() * ggrav * ggrav * rho_c  << endl ;
    
    cout << "GRV2: " << star.grv2() << endl ;
    cout << "GRV3: " << star.grv3() << endl ;

    double vit_triax = diff(6) ;

    //-----------------------------------------------
    //  General features of the final configuration
    //  saved in a file
    //-----------------------------------------------

    ofstream fichfinal("calcul.d") ;
    fichfinal.precision(10) ; 
    
    
    fichfinal << star.get_eos() << endl ;
    
    fichfinal << endl << "Total CPU time  : " << endl ;
    fichfinal << "Memory size : " << endl << endl ; 

    fichfinal << endl << endl ; 
    fichfinal << "Grid : " << endl ; 
    fichfinal << "------ " << endl ; 
    fichfinal << *(star.get_mp().get_mg()) << endl ; 
    fichfinal << endl << "Physical characteristics : " << endl ; 
    fichfinal	  << "-------------------------" << endl ; 
    fichfinal << star << endl ;
    fichfinal << "Growing rate of triaxial perturbation: " << vit_triax 
	      << endl ; 

    fichfinal << endl <<
    "===================================================================" 
    << endl ; 
    //    fichfinal << "Diff_ent : " << diff(0) << endl ; 
    //    fichfinal << "Relative error on the virial theorem GRV2 : "
    //	      << star.grv2() << endl ;   
    //    fichfinal << "Relative error on the virial theorem GRV3 : "
    //	      << star.grv3() << endl ;   
    
    fichfinal << endl <<
    "================================================================" << endl ;
    fichfinal <<
    "   PARAMETERS USED FOR THE COMPUTATION (file parrot.d) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    system("cat parrot.d >> calcul.d") ; 

    fichfinal.open("calcul.d", ios::app) ;
    fichfinal << endl <<
    "================================================================" << endl ;
    fichfinal <<
    "	           EOS PARAMETERS (file par_eos.d) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    system("cat par_eos.d >> calcul.d") ;

    // Identification du code et de ses sous-routines (no. de version RCS) :     	
    fichfinal.open("calcul.d", ios::app) ; 
    fichfinal << endl <<
    "================================================================" << endl ; 
    fichfinal << "	    IDENTIFICATION OF THE CODE : " << endl ; 
    fichfinal << 
    "================================================================" << endl ; 
    fichfinal.close() ; 
//    system("ident rotstar_dirac >> calcul.d") ; 


    // Saveguard of the whole configuration
    // ------------------------------------

	FILE* fresu = fopen("resu.d", "w") ;

	star.get_mp().get_mg()->sauve(fresu) ;		// writing of the grid
	star.get_mp().sauve(fresu) ;                // writing of the mapping
	star.get_eos().sauve(fresu) ;  				// writing of the EOS
	star.sauve(fresu) ;                         // writing of the star
	
	fclose(fresu) ;

 	Sym_tensor hat_aij = star.get_psi4()*star.get_psi2()*star.get_aa() ;
	const Metric_flat& mets = mp.flat_met_spher() ;
	cout << "Looking at \\hat{A}^{ij} and \\hat{A}^{ij}_{TT}: " << endl ;
 	Vector wwws = hat_aij.divergence(mets) ;
	wwws.inc_dzpuis() ;
	Vector www = wwws.poisson(1./3., 6) ;
 	Sym_tensor prueba = hat_aij - www.ope_killing_conf(mets) ;
	Tbl maxa = maxabs(hat_aij, "Max \\hat{A}^ij") ;
 	Tbl maxb = maxabs(prueba, "Max. \\hat{A}^{ij}_{TT}: ") ;
 	Tbl maxh = maxabs(star.get_hh(), "Max. h^ij: ") ;
	cout << "Max(A_TT)/Max(A): " << max(maxb)/max(maxa) <<endl ;
	cout << "Max(h): " << max(maxh) <<endl ;


    // Drawings
    // --------
    
	if (graph == 1) {

// 	for (int i=1; i<=3; i++) 
// 	    for (int j=i; j<=3; j++) {
// 		des_profile(star.get_hh()(i,j), 0., 3*rr, 1, 1) ;
// 	    }

	// Cmp defining the surface of the star (via the enthalpy field)
	Cmp surf(star.get_ent()) ; 
	Cmp surf_ext(mp) ; 
	surf_ext = - 0.2 * surf(0, 0, 0, 0) ; 
	surf_ext.annule(0, star.get_nzet()-1) ; 
	surf.annule(star.get_nzet(), mg.get_nzone()-1) ; 
	surf = surf + surf_ext ;
	surf = raccord_c1(surf, star.get_nzet()) ; 

	int nzdes = star.get_nzet() ; 

	des_coupe_y(star.get_ent(), 0., nzdes, "Enthalpy", &surf) ; 

 	des_coupe_y(star.get_logn(), 0., nzdes, "Gravitational potential \\gn", &surf) ; 
 	des_coupe_y(star.get_beta()(3), 0., nzdes, "Azimuthal shift \\gb\\u\\gf", 
		    &surf) ; 
 	des_coupe_y(star.get_lnq(), 0., nzdes, "Potential ln(Q)", &surf) ; 
	
 	des_coupe_y(star.get_hh()(1,1), 0., nzdes, "Metric potential h\\urr", &surf) ; 
	
	}

 
    // Cleaning
    // --------

    delete peos ;    

    exit(EXIT_SUCCESS) ; 
    
    return EXIT_SUCCESS ; 
   
}
