/*
 * Main code for computing stationary axisymmetric rotating stars in CFC.
 *
 */

/*
 *   Copyright (c) 2024 Jerome Novak
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

 
// headers Lorene
#include "star_rot_cfc.h"
#include "cmp.h"
#include "eos.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"

namespace Lorene {
  // Local prototype (for drawings only)
  Cmp raccord_c1(const Cmp& uu, int l1) ; 
}

//******************************************************************************

using namespace Lorene ;

int main(){

  using namespace Unites ; 

  //------------------------------------------------------------------
  //	    Parameters of the computation 
  //------------------------------------------------------------------

  int relat_i, mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
    delta_mer_kep, mer_mass, graph, nz, nzet, nzadapt,
    nt, np, filter_order ; 
  double ent_c, freq_si, fact_omega, mbar_wanted, precis, freq_ini_si, 
    aexp_mass, relax ;  
    
  ifstream fich("parrot.d") ;
  if ( !fich.good() ) {
    cerr << "Problem in opening the file parrot.d ! " << endl ;
    abort() ;
  }
  // Number of characters to gnore at the end of parameter file lines
  int n_ignore = 120 ;
  
  fich.ignore(n_ignore, '\n') ;
  fich >> relat_i; fich.ignore(n_ignore, '\n') ;
  fich >> ent_c ; fich.ignore(n_ignore, '\n') ;
  fich >> freq_si ; fich.ignore(n_ignore, '\n') ;
  fich >> fact_omega ; fich.ignore(n_ignore, '\n') ;
  fich >> mbar_wanted ; fich.ignore(n_ignore, '\n') ;
  mbar_wanted *= msol ; 
  fich.ignore(n_ignore, '\n') ;
  fich >> mer_max ; fich.ignore(n_ignore, '\n') ;
  fich >> precis ; fich.ignore(n_ignore, '\n') ;
  fich >> mer_rot ; fich.ignore(n_ignore, '\n') ;
  fich >> freq_ini_si ; fich.ignore(n_ignore, '\n') ;
  fich >> mer_change_omega ; fich.ignore(n_ignore, '\n') ;
  fich >> mer_fix_omega ; fich.ignore(n_ignore, '\n') ;
  fich >> delta_mer_kep ; fich.ignore(n_ignore, '\n') ;
  fich >> mer_mass ; fich.ignore(n_ignore, '\n') ;
  fich >> aexp_mass ; fich.ignore(n_ignore, '\n') ;
  fich >> relax ; fich.ignore(n_ignore, '\n') ;
  fich >> graph ; fich.ignore(n_ignore, '\n') ;
  fich.ignore(n_ignore, '\n') ;
  fich >> nz ; fich.ignore(n_ignore, '\n') ;
  fich >> nzet; fich.ignore(n_ignore, '\n') ;
  fich >> nzadapt; fich.ignore(n_ignore, '\n') ;
  fich >> nt; fich.ignore(n_ignore, '\n') ;
  fich >> np; fich.ignore(n_ignore, '\n') ;
  
  int* nr = new int[nz];
  int* nt_tab = new int[nz];
  int* np_tab = new int[nz];
  double* bornes = new double[nz+1];
  
  fich.ignore(n_ignore, '\n');
  for (int l=0; l<nz; l++) {
    fich >> nr[l]; 
    fich >> bornes[l]; fich.ignore(n_ignore, '\n') ;
    np_tab[l] = np ; 
    nt_tab[l] = nt ; 
  }
  bornes[nz] = __infinity ;

  Tbl ent_limit(nzet) ;
  ent_limit.set_etat_qcq() ;
  ent_limit.set(nzet-1) = 0 ; 	// enthalpy at the stellar surface
  for (int l=0; l<nzet-1; l++) {
    fich >> ent_limit.set(l) ; fich.ignore(n_ignore, '\n') ;
  }
  fich >> filter_order ; fich.ignore(n_ignore, '\n') ;
  
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
  if ( !fich.good() ) {
    cerr << "Problem in opening the file par_eos.d ! " << endl ;
    abort() ;
  }
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
    
    Star_rot_CFC star(mp, nzet, eos, relat_i, filter_order) ; 
    

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
    control.set(1) = omega_ini ; 
    control.set(2) = relax ; 

    Tbl diff(8) ;     

    star.equilibrium(ent_c, omega, fact_omega, nzadapt, ent_limit, icontrol,
		     control, mbar_wanted, aexp_mass, diff) ;

     
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
    
    fichfinal << endl <<
    "================================================================" << endl ;
    fichfinal <<
    "   PARAMETERS USED FOR THE COMPUTATION (file parrot.d) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    int ret = system("cat parrot.d >> calcul.d") ; 

    fichfinal.open("calcul.d", ios::app) ;
    fichfinal << endl <<
    "================================================================" << endl ;
    fichfinal <<
    "	           EOS PARAMETERS (file par_eos.d) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    ret = system("cat par_eos.d >> calcul.d") ;

    if (ret == -1) cout << "Warning: problems in identification of the parameters!"
			<< endl ;

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
      
      //      double rr = star.ray_eq() ; 
      
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
      des_coupe_y(star.get_psi(), 0., nzdes, "Conformal factor \\gQ", &surf) ; 
      
    }
    
 
    // Cleaning
    // --------

    delete peos ;    

    exit(EXIT_SUCCESS) ; 
    
    return EXIT_SUCCESS ; 
   
}
