/*
 * Main code for computing a sequence of stationary axisymmetric rigidly
 * rotating stars
 */

/*
 *   Copyright (c) 2001-2003 Eric Gourgoulhon
 *             (c) 2019 Jerome Novak
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
 * Revision 1.1  2019/10/02 13:56:30  j_novak
 * New code for computing sequences of rotating stars
 *
 *
 *
 * $Header$
 *
 */


// headers C++
#include <sstream>

// headers Lorene
#include "star_rot.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"	    

namespace Lorene {
// Local prototype (for drawings only)
Scalar raccord_c1(const Scalar& uu, int l1) ;
}
//******************************************************************************


using namespace Lorene ;

int main(){

  using namespace Unites ;
  
  //------------------------------------------------------------------
  //	    Parameters of the computation 
  //------------------------------------------------------------------
  
  int relat_i, mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
    delta_mer_kep, mer_mass, mermax_poisson, graph, nz, nzet, nzadapt,
    nt, np, mer_triax, n_conf ;
  double entc_min, entc_max, fact_omega, mbar_wanted, precis, freq_min_si,
    freq_max_si, thres_adapt, aexp_mass, relax, relax_poisson,
    ampli_triax, precis_adapt, rrot_km, arot ;  
  
  ifstream fpar("par_seq.d") ;
  if ( !fpar.good() ) {
    cerr << "Problem with opening the file par_seq.d ! " << endl ;
    abort() ;
  }
  
  fpar.ignore(1000,'\n') ;
  fpar >> relat_i ; fpar.ignore(1000,'\n') ;
  bool relat = (relat_i == 1) ; 
  fpar >> entc_min ; fpar.ignore(1000,'\n') ;
  fpar >> entc_max ; fpar.ignore(1000,'\n') ;
  fpar >> freq_min_si ; fpar.ignore(1000,'\n') ;
  fpar >> freq_max_si ; fpar.ignore(1000,'\n') ;
  fpar >> n_conf ; fpar.ignore(1000,'\n') ;    
  if (n_conf > 10000) {
    cout << "nrotseq: n_conf must be smaller than 10000" << endl ;
    abort() ; 
  }
  
  double omega_c_min = 2 * M_PI * freq_min_si / f_unit ; 
  double omega_c_max = 2 * M_PI * freq_max_si / f_unit ; 
  fact_omega = 1 ; 
  fpar >> mbar_wanted ; fpar.ignore(1000,'\n');
  mbar_wanted *= msol ; 
  fpar.ignore(1000,'\n');
  fpar >> mer_max ; fpar.ignore(1000,'\n');
  fpar >> precis ; fpar.ignore(1000,'\n');
  mer_rot = 0 ; 
  mer_change_omega = 0 ; 
  mer_fix_omega = 1 ; 
  delta_mer_kep = 0 ; 
  fpar >> thres_adapt ; fpar.ignore(1000,'\n');
  mer_triax = 100000 ; 
  ampli_triax = 0 ; 
  fpar >> mer_mass ; fpar.ignore(1000,'\n');
  fpar >> aexp_mass ; fpar.ignore(1000,'\n');
  fpar >> relax ; fpar.ignore(1000,'\n');
  fpar >> mermax_poisson ; fpar.ignore(1000,'\n');
  fpar >> relax_poisson ; fpar.ignore(1000,'\n');
  fpar >> precis_adapt ; fpar.ignore(1000,'\n');
  fpar >> graph ; fpar.ignore(1000,'\n');
  fpar.ignore(1000,'\n');
  fpar >> nz ; fpar.ignore(1000,'\n') ;
  fpar >> nzet; fpar.ignore(1000,'\n') ;
  fpar >> nzadapt; fpar.ignore(1000,'\n') ;
  fpar >> nt; fpar.ignore(1000,'\n') ;
  fpar >> np; fpar.ignore(1000,'\n') ;
  
  int* nr = new int[nz];
  int* nt_tab = new int[nz];
  int* np_tab = new int[nz];
  double* bornes = new double[nz+1];
  
  fpar.ignore(1000,'\n') ;
  for (int l=0; l<nz; l++) {
    fpar >> nr[l]; 
    fpar >> bornes[l]; fpar.ignore(1000,'\n');
    np_tab[l] = np ; 
    nt_tab[l] = nt ; 
  }
  bornes[nz] = __infinity ;
  
  Tbl ent_limit(nzet) ;
  ent_limit.set_etat_qcq() ;
  ent_limit.set(nzet-1) = 0 ; 	// enthalpy at the stellar surface
  for (int l=0; l<nzet-1; l++) {
    fpar >> ent_limit.set(l) ; fpar.ignore(1000,'\n');
  }
  
  
  fpar.close();

    
  //-----------------------------------------------------------------------
  //		Equation of state
  //-----------------------------------------------------------------------
  
  fpar.open("par_eos.d") ;
  if ( !fpar.good() ) {
    cerr << "Problem in opening the file par_eos.d ! " << endl ;
    abort() ;
  }
  
  Eos* peos = Eos::eos_from_file(fpar) ;
  Eos& eos = *peos ;

  fpar.close() ;
  
  
  // Special treatment of crust - liquid core boundary in the case
  //  of Eos_strange
  if (eos.identify() == 6) {
    assert( nzet == 2 ) ;    	
    const Eos_strange_cr* peos_cr = dynamic_cast<const Eos_strange_cr*>(peos) ;
    if (peos_cr == 0x0) {
      cout << "nrotseq: problem : peos is not of type Eos_strange_cr !" << endl ;
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
    
    bornes[l+1] = bornes[nzet] * sqrt(1 - ent_limit(l) / entc_min) ;
    
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
  
  cout << "Central enthalpy range: " << entc_min << " - " << entc_max << " c^2" << endl ;
  cout << "Rotation frequency range : " << freq_min_si 
       << " - " << freq_max_si << " Hz" << endl ; 
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
  cout << "Exponent for the increase factor of the central enthalpy or frequency : "
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

  Star_rot star(mp, nzet, relat, eos) ;
  
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
  Scalar ent0(mp) ;
  ent0 = entc_min * ( 1 - r*r / (ray0*ray0) ) ;
  ent0.annule_domain(nz-1) ;
  ent0.std_spectral_base() ;
  star.set_enthalpy(ent0) ;
  
  // Initialization of (n,e,p) from H
  star.equation_of_state() ;
  
  // Initialization of (E,S,U,etc...) (quantities relative to the Eulerian obs)
  star.hydro_euler() ;
  
  cout << endl << "Initial star : "
       << endl << "============   " << endl ;
  
  cout << star << endl ;
  
  //---------------------------------------------------
  //  Sequence parameters saved in file "calcul_seq.d"
  //---------------------------------------------------
  ofstream fichseq("calcul_seq.d") ;
  fichseq << endl <<
    "================================================================" << endl ;
  fichseq <<
    "   PARAMETERS USED FOR THE COMPUTATION (file par_seq.d) : " << endl ;
  fichseq <<
    "================================================================" << endl ;
  fichseq.close() ;
  system("cat par_seq.d >> calcul_seq.d") ;
  
  fichseq.open("calcul_seq.d", ios::app) ;
  fichseq << endl <<
    "================================================================" << endl ;
  fichseq <<
    "	           EOS PARAMETERS (file par_eos.d) : " << endl ;
  fichseq <<
    "================================================================" << endl ;
  fichseq.close() ;
  system("cat par_eos.d >> calcul_seq.d") ;
  
  //-----------------------------------------------------------------------
  //		Computation of the equilibrium sequence
  //-----------------------------------------------------------------------

  
  Itbl icontrol(8) ;
  icontrol.set_etat_qcq() ;
  icontrol.set(0) = mer_max ;
  icontrol.set(1) = mer_rot ;
  icontrol.set(2) = mer_change_omega ;
  icontrol.set(3) = mer_fix_omega ;
  // icontrol.set(4) set later
  icontrol.set(5) = mermax_poisson ;
  icontrol.set(6) = mer_triax ;
  icontrol.set(7) = delta_mer_kep ;
  
  Tbl control(7) ;
  control.set_etat_qcq() ;
  control.set(0) = precis ;
  control.set(2) = relax ;
  control.set(3) = relax_poisson ;
  control.set(4) = thres_adapt ;
  control.set(5) = ampli_triax ;
  control.set(6) = precis_adapt ;

  Tbl diff(8) ;
  
  ofstream fichresu("outseq.d") ;
  fichresu << "# R_circ [km]  R_mean [km]  M_G [M_sol]  M_B [M_sol]  J [G M_sol^2/c] f_c [Hz]     H_c        r_p/r_e"
	   << "         T/W       GRV2         GRV3   " << endl ;
  
  // Loop on the configurations
  // --------------------------
  
  bool seq_freq = ( entc_min == entc_max ) ;
  if ( !seq_freq ) {
    if ( omega_c_min != omega_c_max ) {
      cout << "nrotseq : one must have freq_min == freq_max "
	   << "when entc_min != entc_max !" << endl ;
      abort() ;
    }
  }

  double domega = (n_conf > 1) ? 
    (omega_c_max - omega_c_min) / double(n_conf-1) : 0 ;
  double dent = (n_conf > 1) ? (entc_max - entc_min) / double(n_conf-1) : 0 ;
  
  for (int jj = 0; jj < n_conf; jj++) {

    double omega_c, ent_c ;
    
    if ( seq_freq ) {  
      icontrol.set(4) = mer_mass ;
      ent_c = star.get_ent().val_grid_point(0, 0, 0, 0) ;
      omega_c = omega_c_min + jj * domega ;
    }
    else {
      icontrol.set(4) = - mer_mass ;
      ent_c = entc_min + jj * dent ;
      omega_c = star.get_omega_c() ;
      if (omega_c == double(0)) {
	omega_c = omega_c_min ;
	icontrol.set(1) = 10 ; 	// mer_rot = 10
	icontrol.set(2) = 10 ;
	icontrol.set(3) = 11 ;
      }
      else {
	icontrol.set(1) = 0 ; 	// mer_rot = 0
	icontrol.set(2) = 0 ;
	icontrol.set(3) = 1 ;
      }
    }
    control.set(1) = omega_c ;
    
    star.equilibrium(ent_c, omega_c, fact_omega, nzadapt, ent_limit,
		     icontrol, control, mbar_wanted, aexp_mass, diff) ;
    
    int precisaff = 8 ;
    int tailleaff = precisaff + 4 ;
    fichresu.precision(precisaff) ;
    fichresu << setw(tailleaff)
	     << star.r_circ() / km << " "
	     << setw(tailleaff)
	     << star.mean_radius() / km << " "
	     << setw(tailleaff)
	     << star.mass_g() / msol << " "
	     << setw(tailleaff)
	     << setw(tailleaff)
	     << star.mass_b() / msol << " "
	     << setw(tailleaff)
	     << star.angu_mom()/( qpig / (4* M_PI) * msol*msol) << " "
	     << setw(tailleaff)
	     << star.get_omega_c() / (2.*M_PI) * f_unit << " "
	     << setw(tailleaff)
	     << star.get_ent().val_grid_point(0, 0, 0, 0) << " "
	     << setw(tailleaff)
	     << star.aplat() << " "
	     << setw(tailleaff)
	     << star.tsw() << " "
	     << setw(tailleaff)
	     << star.grv2() << " "
	     << setw(tailleaff)
	     << star.grv3() << endl ;

    cout << endl << "Configuration " << jj << " : "
	 << endl << "===========================" << endl ;

    cout.precision(10) ;
    cout << star << endl ;

    //-----------------------------------------------
    //  General features of the final configuration
    //  saved in a file
    //-----------------------------------------------

    ostringstream oo ;
    oo << "calcul" << setw(4) << setfill('0') << jj << oo.str() << ".d" ;
    cout << endl << "File name : " << oo.str() << endl ;

    ofstream fichfinal(oo.str().data()) ;
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
    fichfinal << "Growing rate of triaxial perturbation: 0 " << endl ; 
    
    fichfinal << endl <<
      "===================================================================" 
	      << endl ; 
    fichfinal << "Diff_ent : " << diff(0) << endl ; 
    fichfinal << "Relative error on the virial theorem GRV2 : "
	      << star.grv2() << endl ;   
    fichfinal << "Relative error on the virial theorem GRV3 : "
	      << star.grv3() << endl ;   
    
    fichfinal.close() ;

    // Drawings
    // --------
    
    if (graph == 1) {
      
      des_map_et(mp, 0) ; 
      
      // Scalar defining the surface of the star (via the enthalpy field)
      Scalar surf = star.get_ent() ; 
      Scalar surf_ext(mp) ; 
      surf_ext = - 0.2 * surf(0, 0, 0, 0) ; 
      surf_ext.annule(0, star.get_nzet()-1) ; 
      surf.annule(star.get_nzet(), mg.get_nzone()-1) ; 
      surf = surf + surf_ext ;
      surf = raccord_c1(surf, star.get_nzet()) ; 
      
      int nzdes = star.get_nzet() ; 
      
      des_coupe_y(star.get_ent(), 0., nzdes, "Enthalpy", &surf) ; 
      
      
      des_coupe_y(star.get_logn(), 0., nzdes,
		  "Gravitational potential \\gn", &surf) ;
      
      if (star.is_relativistic()) {
	
	des_coupe_y(star.get_nphi(), 0., nzdes,
		    "Azimuthal shift N\\u\\gf", &surf) ;
	
	des_coupe_y(star.get_dzeta(), 0., nzdes,
		    "Metric potential \\gz", &surf) ;
	
	des_coupe_y(star.get_tggg(), 0., nzdes,
		    "Metric potential (NB-1) r sin\\gh", &surf) ;
	
	des_coupe_y(star.get_ak_car(), 0., nzdes, 
		    "A\\u2\\d K\\dij\\u K\\uij\\d", &surf) ; 
      }
	
    }


  } // End of the loop on the configurations

  fichresu.close() ; 
 
  // Cleaning
  // --------

  delete peos ;
  
  exit(EXIT_SUCCESS) ; 
  
  return EXIT_SUCCESS ; 
  
}

