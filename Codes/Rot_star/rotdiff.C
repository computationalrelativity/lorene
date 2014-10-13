/*
 * Main code for computing stationary axisymmetric differentially 
 * rotating stars
 */

/*
 *   Copyright (c) 2001-2003 Eric Gourgoulhon
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

char rotdiff_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2014/10/13 08:53:58  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.7  2014/10/06 15:09:45  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2008/05/08 12:02:35  j_novak
 * saving parameters for rotation profiles definitions
 *
 * Revision 1.5  2004/03/25 12:35:44  j_novak
 * now using namespace Unites
 *
 * Revision 1.4  2003/05/25 19:56:49  e_gourgoulhon
 *
 * Added the possibility to choose the factor a = R_eq / R0, instead of R0
 * in the differential rotation law.
 *
 * Revision 1.3  2003/05/14 20:06:09  e_gourgoulhon
 * Suppressed the qualifier ios::nocreate in call to fstream::open
 * (not supported by gcc 3.2).
 *
 * Revision 1.2  2003/01/09 11:07:50  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2001/10/19  08:19:21  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */


// headers C
#include <cstdlib>
#include <cmath>
#include <cstring>

// headers Lorene
#include "et_rot_diff.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"	    

// Function defining the rotation profile
double frotlin(double omega, const Lorene::Tbl& par) ; 
double primfrotlin(double omega, const Lorene::Tbl& par) ; 

namespace Lorene {
// Local prototype (for drawings only)
Cmp raccord_c1(const Cmp& uu, int l1) ; 
}
//******************************************************************************

using namespace Lorene ;

int main(){

  using namespace Unites ;

    // Identification of all the subroutines called by the code : 
    
    system("ident rotdiff > identif.d") ; 

    // For the display : 
    char display_bold[]="x[1m" ; display_bold[0] = 27 ;

    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    char blabla[120] ;

    int relat_i, mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
	delta_mer_kep, mer_mass, mermax_poisson, graph, nz, nzet, nzadapt,
	nt, np, mer_triax, type_rot ; 
    double ent_c, freq_si, fact_omega, mbar_wanted, precis, freq_ini_si, 
	   thres_adapt, aexp_mass, relax, relax_poisson, ampli_triax, 
	   precis_adapt, rrot_si, arot ;  
    
    ifstream fich("parrotdiff.d") ;
    fich.getline(blabla, 120) ;
    fich >> relat_i ; fich.getline(blabla, 120) ;
    bool relat = (relat_i == 1) ; 
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


    // Special treatment of crust - liquid core boundary in the case
    //  of Eos_strange
    if (eos.identify() == 6) {
    	assert( nzet == 2 ) ;    	
    	const Eos_strange_cr* peos_cr = dynamic_cast<const Eos_strange_cr*>(peos) ;
    	if (peos_cr == 0x0) {
    	       cout << "rotdiff: problem : peos is not of type Eos_strange_cr !" << endl ;
    	       abort() ;
    	}
    	
    	ent_limit.set(0) = peos_cr->get_ent_nd() ;  // enthalpy at core/crust transition

    }

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

    Map_et mp(mg, bornes) ;
   
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
    cout << "Maximum number of steps in Map_et::poisson : " 
	 << mermax_poisson << endl ; 
    cout << "Relaxation factor in Map_et::poisson : " 
	 << relax_poisson << endl ; 
    cout << "Step from which the baryon mass is forced to converge : " 
	 << mer_mass << endl ; 
    cout << "Exponent for the increase factor of the central enthalpy : " 
	 << aexp_mass << endl ; 
    cout << 
    "Threshold on |dH/dr|_eq / |dH/dr|_pole for the adaptation of the mapping"
    << endl << thres_adapt << endl ; 


    cout << endl << "Multi-grid : " 
	 << endl << "==========" << endl << mg << endl ; 
    cout << "Mapping : " 
	 << endl << "=======" << endl << mp << endl ; 


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
    
    Et_rot_diff star(mp, nzet, relat, eos, frotlin, primfrotlin, parfrot) ; 
    
    if ( star.is_relativistic() ) {
	cout << "========================" << endl ;
	cout << "Relativistic computation" << endl ;
	cout << "========================" << endl ;
    }
    else {
	cout << "=====================" << endl ;
	cout << "Newtonian computation" << endl ;
	cout << "=====================" << endl ;
    }

    //-----------------------------------------------------------------------
    //		Initialization of the enthalpy field
    //-----------------------------------------------------------------------


    const Coord& r = mp.r ;
    double ray0 = mp.val_r(nzet-1, 1., 0., 0.) ;  
    Cmp ent0(mp) ; 
    ent0 = ent_c * ( 1 - r*r / (ray0*ray0) ) ; 
    ent0.annule(nz-1) ; 
    ent0.std_base_scal() ; 
    star.set_enthalpy(ent0) ;  
    
    // Initialization of (n,e,p) from H
    star.equation_of_state() ; 

    // Initialization of (E,S,U,etc...) (quantities relative to the Eulerian obs)
    star.hydro_euler() ; 

    cout << endl << "Initial star : " 
	 << endl << "============   " << endl ;

    cout << star << endl ; 
     
    //-----------------------------------------------------------------------
    //		Computation of the rotating equilibrium
    //-----------------------------------------------------------------------


    Itbl icontrol(8) ;
    icontrol.set_etat_qcq() ; 
    icontrol.set(0) = mer_max ; 
    icontrol.set(1) = mer_rot ; 
    icontrol.set(2) = mer_change_omega ; 
    icontrol.set(3) = mer_fix_omega ; 
    icontrol.set(4) = mer_mass ; 
    icontrol.set(5) = mermax_poisson ; 
    icontrol.set(6) = mer_triax ; 
    icontrol.set(7) = delta_mer_kep ; 
    
    Tbl control(7) ; 
    control.set_etat_qcq() ; 
    control.set(0) = precis ; 
    control.set(1) = omega_c_ini ; 
    control.set(2) = relax ; 
    control.set(3) = relax_poisson ; 
    control.set(4) = thres_adapt ; 
    control.set(5) = ampli_triax ; 
    control.set(6) = precis_adapt ; 

    Tbl diff(8) ;     

    star.equilibrium(ent_c, omega_c, fact_omega, nzadapt, ent_limit, icontrol, control,
		     mbar_wanted, aexp_mass, diff) ;

     
    cout << endl << "Final star : " 
	 << endl << "==========   " << endl ;

    cout.precision(10) ; 
    cout << star << endl ; 

    double vit_triax = diff(7) ; 

    //-----------------------------------------------
    //  General features of the final configuration
    //  saved in a file
    //-----------------------------------------------

    ofstream fichfinal("calcul.d") ;
    fichfinal.precision(10) ; 
    
    if ( star.is_relativistic() ) {
	fichfinal << "Relativistic computation" << endl ;
    }
    else {
	fichfinal << "Newtonian computation" << endl ;
    }
    
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
    fichfinal << "Diff_ent : " << diff(0) << endl ; 
    fichfinal << "Relative error on the virial theorem GRV2 : "
	      << star.grv2() << endl ;   
    fichfinal << "Relative error on the virial theorem GRV3 : "
	      << star.grv3() << endl ;   
    
    fichfinal << endl <<
    "================================================================" << endl ;
    fichfinal <<
    "   PARAMETERS USED FOR THE COMPUTATION (file parrotdiff.d) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    system("cat parrotdiff.d >> calcul.d") ; 

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
    system("ident rotdiff >> calcul.d") ; 


    // Saveguard of the whole configuration
	// ------------------------------------

	FILE* fresu = fopen("resu.d", "w") ;

	star.get_mp().get_mg()->sauve(fresu) ;		// writing of the grid
	star.get_mp().sauve(fresu) ;                // writing of the mapping
	star.get_eos().sauve(fresu) ;  				// writing of the EOS
	parfrot.sauve(fresu) ;             // writing parameters for rotation function
	star.sauve(fresu) ;                         // writing of the star
	
	fclose(fresu) ;
	
	

    // Drawings
    // --------
    
    if (graph == 1) {

	des_map_et(mp, 0) ; 

	// Cmp defining the surface of the star (via the enthalpy field)
	Cmp surf = star.get_ent()() ; 
	Cmp surf_ext(mp) ; 
	surf_ext = - 0.2 * surf(0, 0, 0, 0) ; 
	surf_ext.annule(0, star.get_nzet()-1) ; 
	surf.annule(star.get_nzet(), mg.get_nzone()-1) ; 
	surf = surf + surf_ext ;
	surf = raccord_c1(surf, star.get_nzet()) ; 

	int nzdes = star.get_nzet() ; 

	des_coupe_y(star.get_ent()(), 0., nzdes, "Enthalpy", &surf) ; 
	
	Cmp tmpdes = star.get_omega_field()() / (2*M_PI) * f_unit ; 
	des_profile(tmpdes, 0., star.ray_eq(),
		    M_PI/2., 0., "\\gW/2\\gp  [Hz]", 
		    "Angular velocity in equatorial plane") ; 

	des_coupe_y(star.get_omega_field()(), 0., nzdes, "\\gW", &surf) ; 

	if (mer_triax < mer_max) { 
	    des_coupe_z(star.get_ent()(), 0., nzdes, "Enthalpy (equatorial plane)", 
			&surf) ; 
	}
	    
	des_coupe_y(star.get_logn()(), 0., nzdes, 
		    "Gravitational potential \\gn", &surf) ; 
	
	if (star.is_relativistic()) {

	    des_coupe_y(star.get_nphi()(), 0., nzdes, 
		    "Azimuthal shift N\\u\\gf", &surf) ; 
	
	    des_coupe_y(star.get_dzeta()(), 0., nzdes, 
		    "Metric potential \\gz", &surf) ; 
	
	    des_coupe_y(star.get_tggg()(), 0., nzdes, 
		    "Metric potential (NB-1) r sin\\gh", &surf) ; 
	
	    des_coupe_y(star.get_ak_car()(), 0., nzdes, 
		    "A\\u2\\d K\\dij\\u K\\uij\\d", &surf) ; 
	}
	
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
