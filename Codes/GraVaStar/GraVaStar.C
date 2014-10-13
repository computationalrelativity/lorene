/*
 * Main code for computing stationary axisymmetric gravastar. 
 * 
 */

/*
 *   Copyright (c) 2010 Frederic Vincent
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

char gravastar_C[] = "$Header$" ;

/*
 *
 *
 * "$Header$"
 *
 */

// C headers
#include <cstdlib>
#include <cmath>
#include <cstring>

// Lorene headers
#include "gravastar.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"	    

namespace Lorene {
// Local prototype (for drawings only)
Scalar raccord_c1(const Scalar& uu, int l1) ; 
/* ***Comment: SERT A QUOI?? */
}
//******************************************************************************

using namespace Lorene ;

int main(){

    using namespace Unites ; 

    // Identification of all the subroutines called by the code : 
    
    //system("ident nrotstar > identif.d") ; 
    /* ***Comment: SERT A QUOI?? */

    // For the display : 
    char display_bold[]="x[1m" ; display_bold[0] = 27 ;

    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    int mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
      delta_mer_kep, mermax_poisson, graph, nz, nzet, nzadapt,
      nt, np; 
    double freq_si, fact_omega, precis, freq_ini_si, 
      thres_adapt, relax, relax_poisson,
      precis_adapt, rho_core ;  
    
    ifstream fpar("par_grav.d") ;
    if ( !fpar.good() ) {
        cerr << "Problem in opening the file par_grav.d ! " << endl ;
        abort() ;
    }
    

    fpar.ignore(1000,'\n') ;    // skip title
    fpar >> freq_si ; fpar.ignore(1000,'\n') ;
    fpar >> fact_omega ; fpar.ignore(1000,'\n') ;
    fpar >> rho_core ; fpar.ignore(1000,'\n') ;
    fpar.ignore(1000,'\n') ;	// skip title
    fpar >> mer_max ; fpar.ignore(1000,'\n') ;
    fpar >> precis ; fpar.ignore(1000,'\n') ;
    fpar >> mer_rot ; fpar.ignore(1000,'\n') ;
    fpar >> freq_ini_si ; fpar.ignore(1000,'\n') ;
    fpar >> mer_change_omega ; fpar.ignore(1000,'\n') ;
    fpar >> mer_fix_omega ; fpar.ignore(1000,'\n') ;
    fpar >> delta_mer_kep ; fpar.ignore(1000,'\n') ;
    fpar >> thres_adapt ; fpar.ignore(1000,'\n') ;
    fpar >> relax ; fpar.ignore(1000,'\n') ;
    fpar >> mermax_poisson ; fpar.ignore(1000,'\n') ;
    fpar >> relax_poisson ; fpar.ignore(1000,'\n') ;
    fpar >> precis_adapt ; fpar.ignore(1000,'\n') ;
    fpar >> graph ; fpar.ignore(1000,'\n') ;
    fpar.ignore(1000,'\n') ; // skip title
    fpar >> nz ; fpar.ignore(1000,'\n') ;
    fpar >> nzet; fpar.ignore(1000,'\n') ;
    fpar >> nzadapt; fpar.ignore(1000,'\n') ;
    fpar >> nt; fpar.ignore(1000,'\n') ;
    fpar >> np; fpar.ignore(1000,'\n') ;

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
     
    fpar.ignore(1000,'\n'); 	// skip title
    for (int l=0; l<nz; l++) {
	fpar >> nr[l]; 
	fpar >> bornes[l]; fpar.ignore(1000,'\n') ;
	np_tab[l] = np ; 
	nt_tab[l] = nt ; 
    }
    bornes[nz] = __infinity ;

    fpar.ignore(1000,'\n'); 	// skip title
    Tbl ent_limit(nzet) ;
    ent_limit.set_etat_qcq() ;
    ent_limit.set(nzet-1) = 0 ;  // enthalpy at the gravastar surface = 0
    for (int l=0; l<nzet-1; l++) {
    	fpar >> ent_limit.set(l) ; fpar.ignore(1000,'\n') ;
    }

    fpar.close();

    // Particular cases
    // ----------------

    // Initial frequency = final frequency
    if ( freq_ini_si < 0 ) {
	freq_ini_si = freq_si ; 
	mer_change_omega = mer_rot ; 
	mer_fix_omega = mer_rot + 1 ;  
    }

    
    //-----------------------------------------------------------------------
    //		Equation of state of the gravastar crust
    //-----------------------------------------------------------------------

    fpar.open("par_grav_eos.d") ;
    if ( !fpar.good() ) {
        cerr << "Problem in opening the file par_grav_eos.d ! " << endl ;
        abort() ;
    }

    Eos* peos = Eos::eos_from_file(fpar) ;
    Eos& eos = *peos ;

    fpar.close() ;

    //-----------------------------------------------------------------------
    //		Construction of the multi-grid and the mapping
    //-----------------------------------------------------------------------

    /* Remark:
       1st domain inside gravastar = vacuum core from bornes[0] to bornes[1]
       2nd domain = first domain inside the crust from bornes[1] to bornes[2]
       ...
       nth domain = last domain inside the crust from bornes[nzet-1] to bornes[nzet]
       (n+1)th domain = first domain outside the gravastar
       (n+2)th domain = last domain (compactified) outside the gravastar

       If there are more than one domain inside the crust, there is a need to rescale bornes[2,...,nzet-1] inside the crust, in order to allow enthalpy(r) to be linear there from H(inner crust boundary) to H(outer crust boundry)=0
    */



    for (int l=1; l<nzet-1; l++) {

      bornes[l+1] = bornes[nzet] + ent_limit(l)/ent_limit(0)*(bornes[1]-bornes[nzet]) ;

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

    cout << "Rotation frequency : " << freq_si << " Hz" << endl ; 
    
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
    cout << 
    "Threshold on |dH/dr|_eq / |dH/dr|_pole for the adaptation of the mapping"
    << endl << thres_adapt << endl ; 


    cout << endl << "Multi-grid : " 
	 << endl << "==========" << endl << mg << endl ; 
    cout << "Mapping : " 
	 << endl << "=======" << endl << mp << endl ; 


     //-----------------------------------------------------------------------
     //		Construction of the gravastar
     //-----------------------------------------------------------------------
    
    Gravastar gravastar(mp,nzet,eos,rho_core);
    
    //-----------------------------------------------------------------------
    //		Initialization of the enthalpy field
    //-----------------------------------------------------------------------


    const Coord& r = mp.r ;
    double ray1 = mp.val_r(0, 1., 0., 0.), ray2 = mp.val_r(nzet-1, 1., 0., 0.) ;
    Scalar ent0(mp) ; 
    ent0 = ent_limit(0) * (r-ray2)/(ray1-ray2) ; //linear profile of enthropy inside the gravastar
    for (int l=nzet;l<nz;l++) ent0.annule_domain(l) ; //sets enth to zero outside the gravastar 
    ent0.std_spectral_base() ; 
    gravastar.set_enthalpy(ent0) ;  

    // des_profile(ent0,0.,10.,0.,0.,"enthalpy0"); 
    
    // Initialization of (n,e,p) from H
    gravastar.equation_of_state() ; 

    // Initialization of (E,S,U,etc...) (quantities relative to the Eulerian obs)
    gravastar.hydro_euler() ; 

    cout << endl << "Initial gravastar : " 
	 << endl << "============   " << endl ;

    cout << gravastar << endl ; 
     
    //-----------------------------------------------------------------------
    //		Computation of the rotating equilibrium
    //-----------------------------------------------------------------------

    double omega = 2 * M_PI * freq_si / f_unit ; 
    double omega_ini = 2 * M_PI * freq_ini_si / f_unit ; 

    Itbl icontrol(6) ;
    icontrol.set_etat_qcq() ; 
    icontrol.set(0) = mer_max ; 
    icontrol.set(1) = mer_rot ; 
    icontrol.set(2) = mer_change_omega ; 
    icontrol.set(3) = mer_fix_omega ; 
    icontrol.set(4) = mermax_poisson ; 
    icontrol.set(5) = delta_mer_kep ; 
    
    Tbl control(6) ; 
    control.set_etat_qcq() ; 
    control.set(0) = precis ; 
    control.set(1) = omega_ini ; 
    control.set(2) = relax ; 
    control.set(3) = relax_poisson ; 
    control.set(4) = thres_adapt ; 
    control.set(5) = precis_adapt ; 

    Tbl diff(8) ;     

    gravastar.equilibrium(omega, fact_omega, nzadapt, ent_limit, icontrol, control,
			  diff) ; 

    cout << endl << "Final gravastar : " 
	 << endl << "==========   " << endl ;

    cout.precision(10) ; 
    cout << gravastar << endl ;

    double rho_c = gravastar.get_ener().val_grid_point(0,0,0,0) ;

    cout << "r_p/r_eq :" << gravastar.aplat() << endl ;
    cout << "Omega rho0^{-1/2} : " << gravastar.get_omega_c() /
    	sqrt( ggrav * rho_c ) << endl ;

    /*cout << "M rho0^{1/2} : " << gravastar.mass_g() * pow(ggrav,1.5) *
    	sqrt( rho_c ) << endl ;

    cout << "M_B rho0^{1/2} : " << gravastar.mass_b() * pow(ggrav,1.5) *
    	sqrt( rho_c ) << endl ;

    cout << "R_circ rho0^{1/2} : " << gravastar.r_circ() *
    sqrt( ggrav * rho_c ) << endl ;*/
    	
    cout << "J rho0 : " << gravastar.angu_mom() * ggrav * ggrav * rho_c  << endl ;
    	
    cout << "Z_p :      " << gravastar.z_pole() << endl ;
    cout << "Z_eq^f :   " << gravastar.z_eqf() << endl ;
    cout << "Z_eq^b :   " << gravastar.z_eqb() << endl ;

    cout << "GRV2: " << gravastar.grv2() << endl ;
    cout << "GRV3: " << gravastar.grv3() << endl ;

    //    double vit_triax = diff(7) ;

    //-----------------------------------------------
    //  General features of the final configuration
    //  saved in a file
    //-----------------------------------------------

    ofstream fichfinal("result.txt") ;
    fichfinal.precision(10) ; 
    
    /* //n'a plus de sens pour le gravastar
      if ( star.is_relativistic() ) {
	fichfinal << "Relativistic computation" << endl ;
    }
    else {
	fichfinal << "Newtonian computation" << endl ;
	}*/
    
    fichfinal << gravastar.get_eos() << endl ;
    
    fichfinal << endl << "Total CPU time  : " << endl ;
    fichfinal << "Memory size : " << endl << endl ; 

    fichfinal << endl << endl ; 
    fichfinal << "Grid : " << endl ; 
    fichfinal << "------ " << endl ; 
    fichfinal << *(gravastar.get_mp().get_mg()) << endl ; 
    fichfinal << endl << "Physical characteristics : " << endl ; 
    fichfinal	  << "-------------------------" << endl ; 
    fichfinal << gravastar << endl ;
    //    fichfinal << "Growing rate of triaxial perturbation: " << vit_triax 
    //      << endl ; 

    fichfinal << endl <<
    "===================================================================" 
    << endl ; 
    fichfinal << "Diff_ent : " << diff(0) << endl ; 
    fichfinal << "Relative error on the virial theorem GRV2 : "
	      << gravastar.grv2() << endl ;   
    fichfinal << "Relative error on the virial theorem GRV3 : "
	      << gravastar.grv3() << endl ;   
    
    fichfinal << endl <<
    "================================================================" << endl ;
    fichfinal <<
    "   PARAMETERS USED FOR THE COMPUTATION (file par_grav.d) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    system("cat par_grav.d >> result.txt") ; 

    fichfinal.open("result.txt", ios::app) ;
    fichfinal << endl <<
    "================================================================" << endl ;
    fichfinal <<
    "	           EOS PARAMETERS (file par_grav_eos.d) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    system("cat par_grav_eos.d >> result.txt") ;

    // Identification du code et de ses sous-routines (no. de version RCS) :     	
    fichfinal.open("result.txt", ios::app) ; 
    fichfinal << endl <<
    "================================================================" << endl ; 
    fichfinal << "	    IDENTIFICATION OF THE CODE : " << endl ; 
    fichfinal << 
    "================================================================" << endl ; 
    fichfinal.close() ; 
    //    system("ident nrotstar >> result.txt") ; 
    /* ***Comment: SERT A QUOI?? */

    // Saveguard of the whole configuration
    // ------------------------------------

	FILE* fresu = fopen("resu.d", "w") ;

	gravastar.get_mp().get_mg()->sauve(fresu) ;		// writing of the grid
	gravastar.get_mp().sauve(fresu) ;                // writing of the mapping
	gravastar.get_eos().sauve(fresu) ;  				// writing of the EOS
	gravastar.sauve(fresu) ;                         // writing of the gravastar
	
	fclose(fresu) ;
	
    // Drawings
    // --------

    ofstream fichdes("prof_n.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       N(theta=0)	N(theta=pi/2)" << endl ; 
    int npd = 300 ;
    double r_max = 3 ;
    double h = r_max/double(npd-1) ;
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << gravastar.get_nn().val_point(rr,0.,0.) << "   " 
	  << gravastar.get_nn().val_point(rr, M_PI/2.,0.) << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_omega.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       omega(theta=0)	  omega(theta=pi/2)" << endl ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << gravastar.get_nphi().val_point(rr,0.,0.) << "   " 
	  << gravastar.get_nphi().val_point(rr, M_PI/2.,0.) << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_omega.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       omega(theta=0)/Omega	  omega(theta=pi/2)/Omega" << endl ; 
    double omega_rot = gravastar.get_omega_c() ; 
    if (omega_rot == double(0)) omega_rot = 1 ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << gravastar.get_nphi().val_point(rr,0.,0.) / omega_rot << "   " 
	  << gravastar.get_nphi().val_point(rr, M_PI/2.,0.) / omega_rot << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_a.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       A(theta=0)	  A(theta=pi/2)" << endl ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << sqrt( gravastar.get_a_car().val_point(rr,0.,0.) ) << "   " 
	  << sqrt( gravastar.get_a_car().val_point(rr, M_PI/2.,0.) ) << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_bma.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       B-A(theta=0)	  B-A(theta=pi/2)" << endl ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << sqrt( gravastar.get_b_car().val_point(rr,0.,0.) ) 
	   - sqrt( gravastar.get_a_car().val_point(rr,0.,0.) ) << "   " 
	  << sqrt( gravastar.get_b_car().val_point(rr, M_PI/2.,0.) ) 
	   - sqrt( gravastar.get_a_car().val_point(rr, M_PI/2.,0.) ) << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_ener.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       ener(theta=0)	  ener(theta=pi/2)" << endl ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << gravastar.get_ener().val_point(rr,0.,0.) << "   " 
	  << gravastar.get_ener().val_point(rr, M_PI/2.,0.) << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_h.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       H(theta=0)	  H(theta=pi/2)" << endl ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << gravastar.get_ent().val_point(rr,0.,0.) << "   " 
	  << gravastar.get_ent().val_point(rr, M_PI/2.,0.) << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_u.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       U(theta=0)/c	  U(theta=pi/2)/c" << endl ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << gravastar.get_uuu().val_point(rr,0.,0.) << "   " 
	  << gravastar.get_uuu().val_point(rr, M_PI/2.,0.) << endl ; 
    }
    fichdes.close() ; 

    
   if (graph == 1) {

	char title[80] ;
	char bslash[2] = {92, '\0'} ;  // 92 is the ASCII code for backslash 

	// Scalar defining the surface of the gravastar (via the enthalpy field)
	Scalar surf( gravastar.get_ent() ); 
	Scalar surf_ext(mp) ; 
	surf_ext = - 0.2 * surf.val_grid_point(0, 0, 0, 0) ; 
	surf_ext.annule(0, gravastar.get_nzet()-1) ; 
	surf.annule(gravastar.get_nzet(), mg.get_nzone()-1) ; 
	surf = surf + surf_ext ;
	surf = raccord_c1(surf, gravastar.get_nzet()) ; 

	int nzdes = gravastar.get_nzet() ; 

	des_coupe_y(gravastar.get_ent(), 0., nzdes, "Log-enthalpy", &surf) ; 

	cout << endl << "Plot of the coefficients of the cos(j theta) expansion of the function G(theta) defining the stellar surface:" << endl ; 
	des_map_et(mp, gravastar.get_nzet()-1) ; 

	des_coupe_y(gravastar.get_ener(), 0., nzdes, "Fluid Proper energy density", &surf) ; 

	/*if (mer_triax < mer_max) { 
	    des_coupe_z(gravastar.get_ent(), 0., nzdes, "Log-enthalpy (equatorial plane)", 
			&surf) ; 
			}*/
	    
	char partit[] = {92, 'g', 'n', '\0'} ; 
	strcpy(title, "Gravitational potential ") ; 
	strcat(title, partit) ; 

	des_coupe_y(gravastar.get_logn(), 0., nzdes, title, &surf) ; 
	
	
	strcpy(title, "Azimuthal shift N") ; 
	strcat(title, bslash) ; 
	strcat(title, "u") ; 
	strcat(title, bslash) ; 
	strcat(title, "gf") ; 
	des_coupe_y(gravastar.get_nphi(), 0., nzdes, title, &surf) ; 
	
	strcpy(title, "Metric potential ") ; 
	strcat(title, bslash) ; 
	strcat(title, "gz") ; 
	des_coupe_y(gravastar.get_dzeta(), 0., nzdes, title, &surf) ; 
	
	strcpy(title, "Metric potential (NB-1) r sin") ; 
	strcat(title, bslash) ; 
	strcat(title, "gh") ; 
	des_coupe_y(gravastar.get_tggg(), 0., nzdes, title, &surf) ; 
	
	char debtit[] = {'A', 92, 'u', '2', 92, 'd', ' ', 'K', 92, 'u', '\0'} ; 
	strcpy(title, debtit) ; 
	strcat(title, "ij") ; 
	strcat(title, bslash) ; 
	strcat(title, "d K") ; 
	strcat(title, bslash) ; 
	strcat(title, "dij") ; 
	strcat(title, bslash) ; 
	strcat(title, "u") ; 

	des_coupe_y(gravastar.get_ak_car(), 0., nzdes, title, &surf) ; 

   }

   /*  ofstream seq("seq.d", ios::app) ; 
    seq << gravastar.get_ent().val_grid_point(0,0,0,0) << "   " << gravastar.mass_g() / msol 
    << "    " << qpig/(4.*M_PI) * gravastar.mass_g() / gravastar.r_circ() << endl ; 
   
    seq.close() ; */
	
    // Cleaning
    // --------

    delete peos ;    

    exit(EXIT_SUCCESS) ; 
    
    return EXIT_SUCCESS ; 
   
}
