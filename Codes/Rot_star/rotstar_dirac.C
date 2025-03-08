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

 

/*
 * $Id$
 * $Log$
 * Revision 1.14  2023/12/20 10:16:19  j_novak
 * New options in Star_rot_Dirac/rotstar_dirac: a flag to choose between GR, CFC and CFC+ A^{ij}_{TT}=0 (as described n Cordero_Carrion et al. 2009).
 *
 * Revision 1.13  2023/03/13 11:18:41  j_novak
 * Added new plots of themro quantities (H, e, nB) + minor changes.
 *
 * Revision 1.12  2016/12/05 16:18:26  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.11  2015/01/27 14:22:39  j_novak
 * New methods in Eos_tabul to correct for EoS themro consistency (optional).
 *
 * Revision 1.10  2014/10/13 08:53:59  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.9  2014/05/09 14:44:25  j_novak
 * Tests are now done for the opening of parameter files.
 *
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
#include <cstdlib>

// headers Lorene
#include "star_rot_dirac.h"
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
    nt, np, filter_order, mer_hij ; 
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
    fich >> mer_hij ;fich.ignore(n_ignore, '\n') ;

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
    
    Star_rot_Dirac star(mp, nzet, eos, relat_i, filter_order) ; 
    

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
    if (max(maxa) > 0.0)
      cout << "Max(A_TT)/Max(A): " << max(maxb)/max(maxa) <<endl ;
    cout << "Max(h): " << max(maxh) <<endl ;
    
    
    // Drawings
    // --------
    
    if (graph == 1) {
      
      double rr = star.ray_eq() ; 
      // des_meridian_log(star.get_ent(), 0., rr, "H", 1) ; 
      // des_meridian_log(star.get_ener(), 0., rr, "e", 2) ; 
      // des_meridian_log(star.get_nbar(), 0., rr, "n\\dB\\u", 3) ;
      // cout << "Type <RETURN> for next page:" << endl ;
      // cin.get() ;
      
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
