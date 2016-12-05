/*
 * Main code for computing stationary axisymmetric differentially rotating 
 * stars in Dirac gauge and maximal slicing.
 * (*** Under development ***)
 *
 */

/*
 *   Copyright (c) 2005 Motoyuki Saijo
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
 *
 * $Header$
 *
 */

// headers C
#include <cstdlib>

// headers Lorene
#include "star_rot_dirac_diff.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"
#include "cmp.h"

// Function defining the rotation profile
double frotlin(double omega, const Lorene::Tbl& par) ; 
double primfrotlin(double omega, const Lorene::Tbl& par) ; 

namespace Lorene {
// Local prototype (for drawings only)
//Scalar raccord_c1(const Scalar& uu, int l1) ; 
Cmp raccord_c1(const Cmp& uu, int l1) ; 
}


//******************************************************************************

using namespace Lorene ;

int main(){

  using namespace Unites ;

    // Identification of all the subroutines called by the code : 
    
    system("ident rotstar_dirac_diff > identif.d") ; 

    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    char blabla[120] ;

    int mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
	delta_mer_kep, mer_mass, mermax_poisson, graph, nz, nzet, nzadapt,
	nt, np, mer_triax, type_rot ; 
    double ent_c, freq_si, fact_omega, mbar_wanted, precis, freq_ini_si, 
	   thres_adapt, aexp_mass, relax, relax_poisson, ampli_triax, 
	   precis_adapt, rrot_si, arot ;  
    
    ifstream fich("parrotdiff.d") ;
    if ( !(fich.is_open()) ) {
	  cerr << "Problem in openning the file parrotdiff.d !" << '\n' ;
	  abort() ;
	}
    fich.getline(blabla, 120) ; 
    fich.getline(blabla, 120) ;
    fich >> ent_c ; fich.getline(blabla, 120) ;
    fich >> freq_si ; fich.getline(blabla, 120) ;
    fich >> type_rot ; fich.getline(blabla, 120) ;
    if (type_rot == 0) {
    	fich >> rrot_si ; fich.getline(blabla, 120) ;
	arot = 0 ; 
    }
    else {
    	assert (type_rot == 1) ; 
    	fich >> arot ; fich.getline(blabla, 120) ;
	rrot_si = 0 ; 	
    }
    fich >> fact_omega ; fich.getline(blabla, 120) ;
    fich >> mbar_wanted ; fich.getline(blabla, 120) ;
    mbar_wanted *= msol ; 
    fich.getline(blabla, 120) ;
    fich >> mer_max ; fich.getline(blabla, 120) ;
    fich >> precis ; fich.getline(blabla, 120) ;
    fich >> mer_rot ; fich.getline(blabla, 120) ;
    fich >> freq_ini_si ; fich.getline(blabla, 120) ;
    fich >> mer_change_omega ; fich.getline(blabla, 120) ;
    fich >> mer_fix_omega ; fich.getline(blabla, 120) ;
    fich >> delta_mer_kep ; fich.getline(blabla, 120) ;
    fich >> thres_adapt ; fich.getline(blabla, 120) ;
    fich >> mer_triax ; fich.getline(blabla, 120) ;
    fich >> ampli_triax ; fich.getline(blabla, 120) ;
    fich >> mer_mass ; fich.getline(blabla, 120) ;
    fich >> aexp_mass ; fich.getline(blabla, 120) ;
    fich >> relax ; fich.getline(blabla, 120) ;
    fich >> mermax_poisson ; fich.getline(blabla, 120) ;
    fich >> relax_poisson ; fich.getline(blabla, 120) ;
    fich >> precis_adapt ; fich.getline(blabla, 120) ;
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
       


    cout << '\n'
	 << "==========================================================" << '\n'
	 << "                    Physical parameters                   " << '\n'
	 << "=========================================================="
	 << '\n'; 
    cout << '\n';

    cout << '\n'<< "Equation of state : " 
	 << '\n'<< "=================   " << '\n';
    cout << eos << '\n'; 

    cout << "Central enthalpy : " << ent_c << " c^2" << '\n'; 
    cout << "Rotation frequency : " << freq_si << " Hz" << '\n'; 
    if ( abs(mer_mass) < mer_max ) {
	cout << "Required Baryon mass [M_sol] : " 
	     << mbar_wanted / msol << '\n'; 
    }
    
    cout << '\n'
	 << "==========================================================" << '\n'
	 << "               Computational parameters                   " << '\n'
	 << "=========================================================="
	 << '\n'<< '\n'; 

    cout << "Maximum number of steps in the main iteration : " 
	 << mer_max << '\n'; 
    cout << "Relaxation factor in the main iteration  : " 
	 << relax << '\n'; 
    cout << "Threshold on the enthalpy relative change for ending the computation : " 
	 << precis << '\n'; 
    cout << "Maximum number of steps in Map_et::poisson : " 
	 << mermax_poisson << '\n'; 
    cout << "Relaxation factor in Map_et::poisson : " 
	 << relax_poisson << '\n'; 
    cout << "Step from which the baryon mass is forced to converge : " 
	 << mer_mass << '\n'; 
    cout << "Exponent for the increase factor of the central enthalpy : " 
	 << aexp_mass << '\n'; 
    cout << 
    "Threshold on |dH/dr|_eq / |dH/dr|_pole for the adaptation of the mapping"
    << '\n'<< thres_adapt << '\n'; 


    cout << '\n'<< "Multi-grid : " 
	 << '\n'<< "==========" << '\n'<< mg << '\n'; 
    cout << "Mapping : " 
	 << '\n'<< "=======" << '\n'<< mp << '\n'; 

    //-----------------------------------------------------------------------
    //		Parameters for the function defining the differential rotation
    //-----------------------------------------------------------------------
    
    double omega_c = 2 * M_PI * freq_si / f_unit ; 
    double omega_c_ini = 2 * M_PI * freq_ini_si / f_unit ; 
    double rrot = rrot_si * km ; 

    Tbl parfrot(3) ;
    parfrot.set_etat_qcq() ; 
    parfrot.set(0) = omega_c_ini ;  
    parfrot.set(1) = rrot ; 
    parfrot.set(2) = arot ;
    
    //-----------------------------------------------------------------------
    //		Construction of the star
    //-----------------------------------------------------------------------
    
    Star_rot_Dirac_diff star(mp, nzet, eos, frotlin, primfrotlin, parfrot) ; 
    

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

    cout << '\n'<< "Initial star : " 
	 << '\n'<< "============   " << '\n';

    cout << star.get_gamma() ;

    cout << star << '\n'; 
     
    //-----------------------------------------------------------------------
    //		Computation of the rotating equilibrium
    //-----------------------------------------------------------------------

    //    double omega = 2 * M_PI * freq_si / f_unit ; 
    //    double omega_ini = 2 * M_PI * freq_ini_si / f_unit ; 

    Itbl icontrol(6) ;
    icontrol.set_etat_qcq() ; 
    icontrol.set(0) = mer_max ; 
    icontrol.set(1) = mer_rot ; 
    icontrol.set(2) = mer_change_omega ; 
    icontrol.set(3) = mer_fix_omega ; 
    icontrol.set(4) = mer_mass ; 
    icontrol.set(5) = delta_mer_kep ; 
    
    Tbl control(3) ; 
    control.set_etat_qcq() ; 
    control.set(0) = precis ; 
    control.set(1) = omega_c_ini ; 
    control.set(2) = relax ; 
    //    control.set(3) = relax_poisson ; 
    //    control.set(4) = thres_adapt ; 
    //    control.set(5) = ampli_triax ; 
    // control.set(6) = precis_adapt ; 

    Tbl diff(8) ;     

    star.equilibrium(ent_c, omega_c, fact_omega, nzadapt, ent_limit, icontrol, control,
		     mbar_wanted, aexp_mass, diff) ;

     
    cout << '\n'<< "Final star : " 
	 << '\n'<< "==========   " << '\n';

    cout.precision(10) ;

    cout << star << '\n';

    double rho_c = star.get_ener().val_grid_point(0,0,0,0) ;
    
    cout << "r_p/r_eq :" << star.aplat() << '\n';

    cout << "Omega rho0^{-1/2} : " << star.get_omega() /
             sqrt( ggrav * rho_c ) << '\n';

    cout << "M rho0^{1/2} : " << star.mass_g() * pow(ggrav,1.5) *
                                 sqrt( rho_c ) << '\n';
    
    cout << "M_B rho0^{1/2} : " << star.mass_b() * pow(ggrav,1.5) *
                                   sqrt( rho_c ) << '\n';
    
    cout << "R_circ rho0^{1/2} : " << star.r_circ() *
                                      sqrt( ggrav * rho_c ) << '\n';
    
    cout << "J rho0 : " << star.angu_mom() * ggrav * ggrav * rho_c 
                        << '\n';
    
    cout << "GRV2: " << star.grv2() << '\n';
    cout << "GRV3: " << star.grv3() << '\n';

    double vit_triax = diff(6) ;

    //-----------------------------------------------
    //  General features of the final configuration
    //  saved in a file
    //-----------------------------------------------

    ofstream fichfinal("calcul.d") ;
    fichfinal.precision(10) ; 
    
    
    fichfinal << star.get_eos() << '\n';
    
    fichfinal << '\n'<< "Total CPU time  : " << '\n';
    fichfinal << "Memory size : " << '\n'<< '\n'; 

    fichfinal << '\n'<< '\n'; 
    fichfinal << "Grid : " << '\n'; 
    fichfinal << "------ " << '\n'; 
    fichfinal << *(star.get_mp().get_mg()) << '\n'; 
    fichfinal << '\n'<< "Physical characteristics : " << '\n'; 
    fichfinal	  << "-------------------------" << '\n'; 
    fichfinal << star << '\n';
    fichfinal << "Growing rate of triaxial perturbation: " << vit_triax 
	      << '\n'; 

    fichfinal << '\n'<<
    "===================================================================" 
    << '\n'; 
    fichfinal << "Diff_ent : " << diff(0) << '\n'; 
    fichfinal << "Relative error on the virial theorem GRV2 : "
    	      << star.grv2() << '\n';   
    fichfinal << "Relative error on the virial theorem GRV3 : "
    	      << star.grv3() << '\n';   
    
    fichfinal << '\n'<<
    "================================================================" << '\n';
    fichfinal <<
    "   PARAMETERS USED FOR THE COMPUTATION (file parrotdiff.d) : " << '\n';
    fichfinal <<
    "================================================================" << '\n';
    fichfinal.close() ;
    system("cat parrotdiff.d >> calcul.d") ; 

    fichfinal.open("calcul.d", ios::app) ;
    fichfinal << '\n'<<
    "================================================================" << '\n';
    fichfinal <<
    "	           EOS PARAMETERS (file par_eos.d) : " << '\n';
    fichfinal <<
    "================================================================" << '\n';
    fichfinal.close() ;
    system("cat par_eos.d >> calcul.d") ;

    // Identification du code et de ses sous-routines (no. de version RCS) :     	
    fichfinal.open("calcul.d", ios::app) ; 
    fichfinal << '\n'<<
    "================================================================" << '\n'; 
    fichfinal << "	    IDENTIFICATION OF THE CODE : " << '\n'; 
    fichfinal << 
    "================================================================" << '\n'; 
    fichfinal.close() ; 
    system("ident rotstar_dirac_diff >> calcul.d") ; 


    // Saveguard of the whole configuration
    // ------------------------------------

	FILE* fresu = fopen("resu.d", "w") ;

	star.get_mp().get_mg()->sauve(fresu) ;		// writing of the grid
	star.get_mp().sauve(fresu) ;                // writing of the mapping
	star.get_eos().sauve(fresu) ;  				// writing of the EOS
	star.sauve(fresu) ;                         // writing of the star
	
	fclose(fresu) ;
	
	

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

	Cmp tmpdes = star.get_omega_field() / (2*M_PI) * f_unit ; 
	des_profile(tmpdes, 0., star.ray_eq(),
		    M_PI/2., 0., "\\gW/2\\gp  [Hz]", 
		    "Angular velocity in equatorial plane") ; 

	des_coupe_y(star.get_omega_field(), 0., nzdes, "\\gW", &surf) ; 

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

		//------------------------------//
		//	Function F(Omega)	//
		//------------------------------//

double frotlin(double omega, const Tbl& par){
    
	double omega_c = par(0) ; 
	double rrot = par(1) ; 
    
	return rrot*rrot* (omega_c - omega) ;  
    
}

		//------------------------------//
		//	Primitive of  F(Omega)	//
		//------------------------------//

double primfrotlin(double omega, const Tbl& par){
    
	double omega_c = par(0) ; 
	double rrot = par(1) ; 
    
	return - 0.5 * rrot*rrot* (omega_c - omega)*(omega_c - omega) ;  
    
}
