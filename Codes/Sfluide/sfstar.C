
// headers C++
#include <iostream.h>
#include <fstream.h>

// headers C
#include <stdlib.h>
#include <math.h>
#include <string.h>

// headers Lorene
#include "et_rot_bifluid.h"
#include "eos_bifluid.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"

// Local prototype (for drawings only)
Cmp raccord_c1(const Cmp& uu, int l1) ; 

//******************************************************************************

int main(){

    // For the display : 
    char display_bold[]="x[1m" ; display_bold[0] = 27 ;

    #include "unites.h"	    
    // To avoid some compilation warnings
    if (display_bold == 0x0) {
      cout << qpig << f_unit << km << msol << mevpfm3 << endl ; 
    }    
    
    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    char blabla[120] ;

    int relat_i, mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
	mermax_poisson, graph, nz, nzet, nt, np; 
    double ent_c, ent2_c, freq_si, freq2_si, precis, freq_ini_si, 
      freq2_ini_si,relax, relax_poisson ;  
    
    ifstream fich("parrot.d") ;
    fich.getline(blabla, 120) ;
    fich >> relat_i ; fich.getline(blabla, 120) ;
    bool relat = (relat_i == 1) ; 
    fich >> ent_c ; fich.getline(blabla, 120) ;
    fich >> ent2_c ; fich.getline(blabla, 120) ;
    fich >> freq_si ; fich.getline(blabla, 120) ;
    fich >> freq2_si ; fich.getline(blabla, 120) ;
    fich.getline(blabla, 120) ;
    fich >> mer_max ; fich.getline(blabla, 120) ;
    fich >> precis ; fich.getline(blabla, 120) ;
    fich >> mer_rot ; fich.getline(blabla, 120) ;
    fich >> freq_ini_si ; fich.getline(blabla, 120) ;
    fich >> freq2_ini_si ; fich.getline(blabla, 120) ;
    fich >> mer_change_omega ; fich.getline(blabla, 120) ;
    fich >> mer_fix_omega ; fich.getline(blabla, 120) ;
    fich >> relax ; fich.getline(blabla, 120) ;
    fich >> mermax_poisson ; fich.getline(blabla, 120) ;
    fich >> relax_poisson ; fich.getline(blabla, 120) ;
    fich >> graph ; fich.getline(blabla, 120) ;
    fich.getline(blabla, 120) ;
    fich >> nz ; fich.getline(blabla, 120) ;
    fich >> nzet; fich.getline(blabla, 120) ;
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
    Tbl ent2_limit(ent_limit) ;
    for (int l=0; l<nzet-1; l++) {
    	fich >> ent_limit.set(l) ; fich.getline(blabla, 120) ;
	ent2_limit.set(l) = ent_limit(l) ; 
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

    Eos_bifluid* peos = Eos_bifluid::eos_from_file(fich) ;
    Eos_bifluid& eos = *peos ;

    fich.close() ;

    //-----------------------------------------------------------------------
    //		Construction of the multi-grid and the mapping
    //-----------------------------------------------------------------------

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

    cout << "Central enthalpy 1 : " << ent_c << " c^2" << endl ; 
    cout << "Central enthalpy 2 : " << ent2_c << " c^2" << endl ; 
    cout << "Rotation frequency : " << freq_si << " Hz" << endl ; 
    cout << "Rotation frequency 2 : " << freq2_si << " Hz" << endl ; 
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

    cout << endl << "Multi-grid : " 
	 << endl << "==========" << endl << mg << endl ; 
    cout << "Mapping : " 
	 << endl << "=======" << endl << mp << endl ; 


    //-----------------------------------------------------------------------
    //		Construction of the star
    //-----------------------------------------------------------------------
    
    Et_rot_bifluid star(mp, nzet, relat, eos) ;

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
    Cmp enta(mp) ; 
    enta = ent_c * ( 1 - r*r / (ray0*ray0) ) ; 
    enta.annule(nz-1) ; 
    enta.std_base_scal() ; 
    Cmp entb = enta * ent2_c/ent_c ; 
    star.set_enthalpies(enta, entb) ;  
    
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

    double omega2 = 2 * M_PI * freq2_si / f_unit ; 
    double omega2_ini = 2 * M_PI * freq2_ini_si / f_unit ; 

    Itbl icontrol(5) ;
    icontrol.set_etat_qcq() ; 
    icontrol.set(0) = mer_max ; 
    icontrol.set(1) = mer_rot ; 
    icontrol.set(2) = mer_change_omega ; 
    icontrol.set(3) = mer_fix_omega ; 
    icontrol.set(4) = mermax_poisson ; 
    
    Tbl control(5) ; 
    control.set_etat_qcq() ; 
    control.set(0) = precis ; 
    control.set(1) = omega_ini ;
    control.set(2) = omega2_ini ;
    control.set(3) = relax ; 
    control.set(4) = relax_poisson ; 

    Tbl diff(8) ;     

    star.equilibrium_bi(ent_c, ent2_c, omega, omega2, ent_limit, 
		     ent2_limit, icontrol, control, diff) ;

     
    cout << endl << "Final star : " 
	 << endl << "==========   " << endl ;

    cout.precision(10) ; 
    cout << star << endl ; 

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
    star.get_eos().sauve(fresu) ;  		// writing of the EOS
    star.sauve(fresu) ;                         // writing of the star
    
    fclose(fresu) ;
    
    
    
    // Drawings
    // --------
    
    if (graph == 1) {
      
      char title[80] ;
      
      int nzdes = star.get_nzet() ; 
      
      // Cmp defining the surface of the star (via the density fields)
      // 
      Cmp t1(enta) ;
      Cmp t2(entb) ;
      peos->nbar_ent(star.get_ent()(), star.get_ent2()(), star.get_xxx2()(),
		     t1, t2, star.get_nzet(), 0, false) ;
      Cmp surf(t1) ;
      Cmp surf2(t2) ;

      des_bi_coupe_y(star.get_nbar()(), 0., nzdes, "Fluid 1 baryonic density", 
		     &surf, &surf) ; 

      des_bi_coupe_y(star.get_nbar2()(), 0., nzdes, "Fluid 2 baryonic density", 
      	     &surf2, &surf2) ; 
      
      des_bi_coupe_y(star.get_ent()(), 0., nzdes, "Fluid 1 Enthalpy", 
		     &surf, &surf2) ; 

      des_bi_coupe_y(star.get_ent2()(), 0., nzdes, "Fluid 2 Enthalpy", 
      	     &surf, &surf2) ; 
      strcpy(title, "Gravitational potential \\gn") ; 
      des_bi_coupe_y(star.get_logn()(), 0., nzdes, title, &surf, &surf2) ; 
      
      strcpy(title, "Azimuthal shift N\\u\\gf") ; 
      des_bi_coupe_y(star.get_nphi()(), 0., nzdes, title, &surf, &surf2) ; 
      
      strcpy(title, "Metric potential \\gz") ; 
      des_bi_coupe_y(star.get_dzeta()(), 0., nzdes, title, &surf, &surf2) ; 
      
      strcpy(title, "Metric potential (NB-1) r sin\\gh") ; 
      des_bi_coupe_y(star.get_tggg()(), 0., nzdes, title, &surf, &surf2) ; 
      
      
      strcpy(title, "A\\u2\\d K\\uij\\d K\\dij\\u") ; 
      des_bi_coupe_y(star.get_ak_car()(), 0., nzdes, title, &surf, &surf2) ; 
	
    }

 
    // Cleaning
    // --------
    
    delete peos ;    

    exit(EXIT_SUCCESS) ; 
    
    return EXIT_SUCCESS ; 
    
}
