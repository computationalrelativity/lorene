/*
 * Main code for computing stationary axisymmetric rotating stars. 
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

char rotstar_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2003/01/09 11:07:51  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.2  2002/03/28 08:58:33  e_gourgoulhon
 * New Makefile, to produce any of the man codes
 * Updated README file
 * Template parameter files for rotseq
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.16  2001/10/10  13:58:04  eric
 * Modif Joachim: traitement des legendes graphiques
 *
 * Revision 1.14  2000/11/29  13:47:41  eric
 * Correction erreur initialisation du Tbl ent_limit.
 *
 * Revision 1.13  2000/11/24  13:28:57  eric
 * Initialisation de bornes dans le cas nzet >= 2.
 *
 * Revision 1.12  2000/11/23  15:44:41  eric
 * Version pour calcul avec 2 zones ou plus dans l'etoile
 *  --> Ajout de l'argument ent_limit a Etoile_rot::equilibrium
 *
 * Revision 1.11  2000/11/19  22:35:12  eric
 * Sauvegarde configuration complete a la fin
 * Ecriture parametres EOS dans calcul.d
 *
 * Revision 1.10  2000/11/10  15:22:57  eric
 * Ajout des parametres de calcul delta_mer_kep et precis_adapt
 *  (lus dans le fichier parrot.d)
 * Ajout de graphiques a la fin.
 *
 * Revision 1.9  2000/11/08  15:27:43  eric
 * Ajout de cout.precision(10) avant l'affichage final.
 *
 * Revision 1.8  2000/10/25  15:43:27  eric
 * Creation du fichier calcul.d a la fin.
 *
 * Revision 1.7  2000/10/23  13:50:17  dorota
 * Sortie de vit_triax.
 *
 * Revision 1.6  2000/10/20  13:58:07  eric
 * Ajout de l'argument nzadapt a Etoile_rot::equilibrium
 * Divers dessins.
 *
 * Revision 1.5  2000/10/17  16:01:26  eric
 * Ajout de la perturbation triaxiale.
 *
 * Revision 1.4  2000/10/16  09:25:01  eric
 * *** empty log message ***
 *
 * Revision 1.3  2000/08/25  12:29:12  eric
 * Dessins a la fin.
 *
 * Revision 1.2  2000/08/18  14:02:39  eric
 * *** empty log message ***
 *
 * Revision 1.1  2000/07/21  16:33:09  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// headers C
#include <stdlib.h>
#include <math.h>
#include <string.h>

// headers Lorene
#include "etoile.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"

// Local prototype (for drawings only)
Cmp raccord_c1(const Cmp& uu, int l1) ; 

//******************************************************************************

int main(){

    // Identification of all the subroutines called by the code : 
    
    system("ident rotstar > identif.d") ; 

    // For the display : 
    char display_bold[]="x[1m" ; display_bold[0] = 27 ;

    #include "unites.h"	    
    // To avoid some compilation warnings
    if (display_bold == 0x0) {
	cout << qpig << f_unit << km << mevpfm3 << endl ; 
    }    
    
    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    char blabla[120] ;

    int relat_i, mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
	delta_mer_kep, mer_mass, mermax_poisson, graph, nz, nzet, nzadapt,
	nt, np, mer_triax ; 
    double ent_c, freq_si, fact_omega, mbar_wanted, precis, freq_ini_si, 
	   thres_adapt, aexp_mass, relax, relax_poisson, ampli_triax, 
	   precis_adapt ;  
    
    ifstream fich("parrot.d") ;
    fich.getline(blabla, 120) ;
    fich >> relat_i ; fich.getline(blabla, 120) ;
    bool relat = (relat_i == 1) ; 
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
    	       cout << "rotstar: problem : peos is not of type Eos_strange_cr !" << endl ;
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
    //		Construction of the star
    //-----------------------------------------------------------------------
    
    Etoile_rot star(mp, nzet, relat, eos) ; 
    
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

    double omega = 2 * M_PI * freq_si / f_unit ; 
    double omega_ini = 2 * M_PI * freq_ini_si / f_unit ; 

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
    control.set(1) = omega_ini ; 
    control.set(2) = relax ; 
    control.set(3) = relax_poisson ; 
    control.set(4) = thres_adapt ; 
    control.set(5) = ampli_triax ; 
    control.set(6) = precis_adapt ; 

    Tbl diff(8) ;     

    star.equilibrium(ent_c, omega, fact_omega, nzadapt, ent_limit, icontrol, control,
		     mbar_wanted, aexp_mass, diff) ;

     
    cout << endl << "Final star : " 
	 << endl << "==========   " << endl ;

    cout.precision(10) ; 
    cout << star << endl ;

    double rho_c = star.get_ener()()(0,0,0,0) ;

    cout << "r_p/r_eq :" << star.aplat() << endl ;
    cout << "Omega rho0^{-1/2} : " << star.get_omega_c() /
    	sqrt( ggrav * rho_c ) << endl ;

    cout << "M rho0^{1/2} : " << star.mass_g() * pow(ggrav,1.5) *
    	sqrt( rho_c ) << endl ;

    cout << "M_B rho0^{1/2} : " << star.mass_b() * pow(ggrav,1.5) *
    	sqrt( rho_c ) << endl ;

    cout << "R_circ rho0^{1/2} : " << star.r_circ() *
    	sqrt( ggrav * rho_c ) << endl ;
    	
    cout << "J rho0 : " << star.angu_mom() * ggrav * ggrav * rho_c  << endl ;
    	
    cout << "Z_p :      " << star.z_pole() << endl ;
    cout << "Z_eq^f :   " << star.z_eqf() << endl ;
    cout << "Z_eq^b :   " << star.z_eqb() << endl ;

    cout << "GRV2: " << star.grv2() << endl ;
    cout << "GRV3: " << star.grv3() << endl ;

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
    "   PARAMETERS USED FOR THE COMPUTATION (file parrot.d) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    system("cat parrot.d >> calcul.d") ; 

    fichfinal.open("calcul.d", ios::app | ios::nocreate) ;
    fichfinal << endl <<
    "================================================================" << endl ;
    fichfinal <<
    "	           EOS PARAMETERS (file par_eos.d) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    system("cat par_eos.d >> calcul.d") ;

    // Identification du code et de ses sous-routines (no. de version RCS) :     	
    fichfinal.open("calcul.d", ios::app | ios::nocreate) ; 
    fichfinal << endl <<
    "================================================================" << endl ; 
    fichfinal << "	    IDENTIFICATION OF THE CODE : " << endl ; 
    fichfinal << 
    "================================================================" << endl ; 
    fichfinal.close() ; 
    system("ident rotstar >> calcul.d") ; 


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

	des_map_et(mp, 0) ; 


	char title[80] ;
	char bslash[2] = {92, '\0'} ;  // 92 is the ASCII code for backslash 

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

	if (mer_triax < mer_max) { 
	    des_coupe_z(star.get_ent()(), 0., nzdes, "Enthalpy (equatorial plane)", 
			&surf) ; 
	}
	    
	char partit[] = {92, 'g', 'n', '\0'} ; 
	strcpy(title, "Gravitational potential ") ; 
	strcat(title, partit) ; 

	des_coupe_y(star.get_logn()(), 0., nzdes, title, &surf) ; 
	
	
	strcpy(title, "Azimuthal shift N") ; 
	strcat(title, bslash) ; 
	strcat(title, "u") ; 
	strcat(title, bslash) ; 
	strcat(title, "gf") ; 
	des_coupe_y(star.get_nphi()(), 0., nzdes, title, &surf) ; 
	
	strcpy(title, "Metric potential ") ; 
	strcat(title, bslash) ; 
	strcat(title, "gz") ; 
	des_coupe_y(star.get_dzeta()(), 0., nzdes, title, &surf) ; 
	
	strcpy(title, "Metric potential (NB-1) r sin") ; 
	strcat(title, bslash) ; 
	strcat(title, "gh") ; 
	des_coupe_y(star.get_tggg()(), 0., nzdes, title, &surf) ; 
	

	char debtit[] = {'A', 92, 'u', '2', 92, 'd', ' ', 'K', 92, 'u', '\0'} ; 
	strcpy(title, debtit) ; 
	strcat(title, "ij") ; 
	strcat(title, bslash) ; 
	strcat(title, "d K") ; 
	strcat(title, bslash) ; 
	strcat(title, "dij") ; 
	strcat(title, bslash) ; 
	strcat(title, "u") ; 

	des_coupe_y(star.get_ak_car()(), 0., nzdes, title, &surf) ; 

    }

 
    // Cleaning
    // --------

    delete peos ;    

    exit(EXIT_SUCCESS) ; 
    
    return EXIT_SUCCESS ; 
   
}
