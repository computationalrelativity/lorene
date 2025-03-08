/*
 * Main code for computing magnetized stationary axisymmetric rotating stars. 
 * 
 */

/*
 *   Copyright (c) 2002 Jerome Novak & Emanuel Marcq
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
 * Revision 1.15  2022/02/23 13:21:35  j_novak
 * Removed output file identification.d
 *
 * Revision 1.14  2014/10/13 08:53:57  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.13  2014/09/03 15:34:02  j_novak
 * Filtering of Maxwell sources is now optional.
 *
 * Revision 1.12  2014/07/04 12:42:41  j_novak
 * Parameter file for mass-shedding configurations.
 *
 * Revision 1.11  2012/08/10 12:52:42  j_novak
 * Minor changes.
 *
 * Revision 1.10  2006/03/09 08:11:59  j_novak
 * Output of the divergence of B at the end of the iteration.
 *
 * Revision 1.9  2004/03/25 12:35:42  j_novak
 * now using namespace Unites
 *
 * Revision 1.8  2003/05/14 20:06:09  e_gourgoulhon
 * Suppressed the qualifier ios::nocreate in call to fstream::open
 * (not supported by gcc 3.2).
 *
 *
 * $Header$
 *
 */

// headers Lorene
#include "et_rot_mag.h"
#include "metric.h"
#include "eos.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"


namespace Lorene{
  // Local prototype (for drawings only)
  Cmp raccord_c1(const Cmp& uu, int l1) ; 
}

using namespace Lorene ;

// Definition de f_j et de M_j
Cmp f_j(const Cmp& x, double a_j){ 
  Cmp resu(x.get_mp() ) ;
  resu = a_j ;
  return resu ;
}
Cmp M_j(const Cmp& x, double a_j){ 
  return - a_j*x ;
}
//******************************************************************************

int main(){

  using namespace Unites_mag ;
  
  //------------------------------------------------------------------
  //	    Parameters of the computation 
  //------------------------------------------------------------------
  
  int relat_i, mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
    delta_mer_kep, mer_mass, mermax_poisson, graph, nz, nzet, nzadapt,
    nt, np ; 
  double ent_c, freq_si, fact_omega, mbar_wanted, precis, freq_ini_si, 
    thres_adapt, aexp_mass, relax, relax_poisson, precis_adapt ;  
  
  double Q0, a_j0, Q_ini, a_j_ini ;
  int mer_mag, mer_change_mag, mer_fix_mag, mag_filter, conduct ;
  
  ifstream fich("parrot.d") ;
  fich.ignore(1000, '\n') ;
  fich >> relat_i ; fich.ignore(1000, '\n') ;
  bool relat = (relat_i == 1) ; 
  fich >> ent_c ; fich.ignore(1000, '\n') ;
  fich >> freq_si ; fich.ignore(1000, '\n') ;
  fich >> fact_omega ; fich.ignore(1000, '\n') ;
  fich >> mbar_wanted ; fich.ignore(1000, '\n') ;
  mbar_wanted *= msol ; 
  fich >> conduct ; fich.ignore(1000, '\n') ;
  fich.ignore(1000, '\n') ;
  fich >> Q0 ; fich.ignore(1000, '\n') ;
  fich >> a_j0 ; fich.ignore(1000, '\n') ;
  fich >> Q_ini; fich.ignore(1000, '\n') ;
  fich >> a_j_ini ; fich.ignore(1000, '\n') ;
  fich >> mer_mag ; fich.ignore(1000, '\n') ;
  fich >> mer_change_mag ; fich.ignore(1000, '\n') ;
  fich >> mer_fix_mag ; fich.ignore(1000, '\n');
  fich >> mag_filter ; fich.ignore(1000, '\n') ;
  fich.ignore(1000, '\n');
  fich >> mer_max ; fich.ignore(1000, '\n') ;
  fich >> precis ; fich.ignore(1000, '\n') ;
  fich >> mer_rot ; fich.ignore(1000, '\n') ;
  fich >> freq_ini_si ; fich.ignore(1000, '\n') ;
  fich >> mer_change_omega ; fich.ignore(1000, '\n') ;
  fich >> mer_fix_omega ; fich.ignore(1000, '\n') ;
  fich >> delta_mer_kep ; fich.ignore(1000, '\n') ;
  fich >> thres_adapt ; fich.ignore(1000, '\n') ;
  fich >> mer_mass ; fich.ignore(1000, '\n') ;
  fich >> aexp_mass ; fich.ignore(1000, '\n') ;
  fich >> relax ; fich.ignore(1000, '\n') ;
  fich >> mermax_poisson ; fich.ignore(1000, '\n') ;
  fich >> relax_poisson ; fich.ignore(1000, '\n') ;
  fich >> precis_adapt ; fich.ignore(1000, '\n') ;
  fich >> graph ; fich.ignore(1000, '\n') ;
  fich.ignore(1000, '\n') ;
  fich >> nz ; fich.ignore(1000, '\n') ;
  fich >> nzet; fich.ignore(1000, '\n') ;
  fich >> nzadapt; fich.ignore(1000, '\n') ;
  fich >> nt; fich.ignore(1000, '\n') ;
  fich >> np; fich.ignore(1000, '\n') ;
  
  int* nr = new int[nz];
  int* nt_tab = new int[nz];
  int* np_tab = new int[nz];
  double* bornes = new double[nz+1];
  
  fich.ignore(1000, '\n');
  for (int l=0; l<nz; l++) {
    fich >> nr[l]; 
    fich >> bornes[l]; fich.ignore(1000, '\n') ;
    np_tab[l] = np ; 
    nt_tab[l] = nt ; 
  }
  bornes[nz] = __infinity ;
  
  Tbl ent_limit(nzet) ;
  ent_limit.set_etat_qcq() ;
  ent_limit.set(nzet-1) = 0 ; 	// enthalpy at the stellar surface
  for (int l=0; l<nzet-1; l++) {
    fich >> ent_limit.set(l) ; fich.ignore(1000, '\n') ;
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
  
  Et_rot_mag star(mp, nzet, relat, eos, conduct) ; 
  
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
  
  // Initialization of eulerian electromagnetic quantities.
  star.MHD_comput() ;
  

  cout << endl << "Initial star : " 
       << endl << "============   " << endl ;
  
  cout << star << endl ; 
  
  //-----------------------------------------------------------------------
  //		Computation of the rotating equilibrium
  //-----------------------------------------------------------------------

  double omega = 2 * M_PI * freq_si / f_unit ; 
  double omega_ini = 2 * M_PI * freq_ini_si / f_unit ; 
  
  Itbl icontrol(11) ;
  icontrol.set_etat_qcq() ; 
  icontrol.set(0) = mer_max ; 
  icontrol.set(1) = mer_rot ; 
  icontrol.set(2) = mer_change_omega ; 
  icontrol.set(3) = mer_fix_omega ; 
  icontrol.set(4) = mer_mass ; 
  icontrol.set(5) = mermax_poisson ;  
  icontrol.set(6) = delta_mer_kep ; 
  icontrol.set(7) = mer_mag ;
  icontrol.set(8) = mer_change_mag ;
  icontrol.set(9) = mer_fix_mag ;
  icontrol.set(10) = mag_filter ;
  
  Tbl control(8) ; 
  control.set_etat_qcq() ; 
  control.set(0) = precis ; 
  control.set(1) = omega_ini ; 
  control.set(2) = relax ; 
  control.set(3) = relax_poisson ; 
  control.set(4) = thres_adapt ; 
  control.set(5) = precis_adapt ; 
  control.set(6) = Q_ini ;
  control.set(7) = a_j_ini ;
  
  Tbl diff(1) ;     
  
  star.equilibrium_mag(ent_c, omega, fact_omega, nzadapt, ent_limit, 
		       icontrol, control, mbar_wanted, aexp_mass, diff, 
		       Q0, a_j0, &f_j, &M_j) ;

     
  cout << endl << "Final star : " 
       << endl << "==========   " << endl ;
  
  cout.precision(10) ; 
  cout << star << endl ; 
    
  const Base_vect_spher& bspher = mp.get_bvect_spher() ;
  Sym_tensor gij(mp, COV, bspher) ;
  gij.set_etat_zero() ;
  gij.set(1,1) = star.get_a_car()() ;
  gij.set(2,2) = star.get_a_car()() ;
  gij.set(3,3) = star.get_b_car()() ;
  Metric gam(gij) ;
  
  Scalar fac = sqrt(star.get_a_car()()) ;
  fac.std_spectral_base() ;
  Vector Bmag(mp, CON, bspher) ;
  Bmag.set(1) = Scalar(star.Magn()(0)) / fac ;
  Bmag.set(1).dec_dzpuis(2) ;
  Bmag.set(2) = Scalar(star.Magn()(1)) / fac ;
  Bmag.set(2).dec_dzpuis(2) ;
  Bmag.set(3) = 0 ;
  Tbl maxB = 0.5*(max(abs(Bmag(1))) + max(abs(Bmag(2)))) ;
  cout << maxB ;
  cout << "div(B) / max(B) in each domain :  "
       << max(abs(Bmag.divergence(gam))) / maxB ; 
    
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
  int res = system("cat parrot.d >> calcul.d") ; 
  if (res==-1) {
    cout << "Problem in shell command:" << endl ;
    cout << "cat parrot.d >> calcul.d" << endl ;
    cout << "calcul.d may be altered." << endl ;
  }
  
  fichfinal.open("calcul.d", ios::app) ;
  fichfinal << endl <<
    "================================================================" << endl ;
  fichfinal <<
    "	           EOS PARAMETERS (file par_eos.d) : " << endl ;
  fichfinal <<
    "================================================================" << endl ;
  fichfinal.close() ;
  res = system("cat par_eos.d >> calcul.d") ;
  if (res==-1) {
    cout << "Problem in shell command:" << endl ;
    cout << "cat par_eos.d >> calcul.d" << endl ;
    cout << "calcul.d may be altered." << endl ;
  }
  
  // Saveguard of the whole configuration
  // ------------------------------------
  
  FILE* fresu = fopen("resu.d", "w") ;
  
  star.get_mp().get_mg()->sauve(fresu) ;      // writing of the grid
  star.get_mp().sauve(fresu) ;                // writing of the mapping
  star.get_eos().sauve(fresu) ;  	      // writing of the EOS
  star.sauve(fresu) ;                         // writing of the star
  
  fclose(fresu) ;
	
	
  // Drawings
  // --------
  
  if (graph == 1) {
    
    // Cmp defining the surface of the star (via the enthalpy field)
    Cmp surf = star.get_ent()() ; 
    Cmp surf_ext(mp) ; 
    surf_ext = - 0.2 * surf(0, 0, 0, 0) ; 
    surf_ext.annule(0, star.get_nzet()-1) ; 
    surf.annule(star.get_nzet(), mg.get_nzone()-1) ; 
    surf = surf + surf_ext ;
    surf = raccord_c1(surf, star.get_nzet()) ; 
    
    int nzdes = star.get_nzet() ; 
    
    des_coupe_y(star.get_At(), 0., nzdes, "A\\dt\\u Potential", &surf) ; 
    des_coupe_y(star.get_Aphi(), 0., nzdes, "Magnetic field", &surf) ; 
    des_coupe_y(star.get_ent()(), 0., nzdes, "Enthalpy", &surf) ; 
    
    des_coupe_y(star.get_logn()(), 0., nzdes, "Gravitational potential \\gn",
		&surf) ; 
       
    des_coupe_y(star.get_nphi()(), 0., nzdes, "Azimuthal shift N\\u\\gf", &surf) ; 
    
    des_coupe_y(star.get_dzeta()(), 0., nzdes, "Metric potential \\gz", &surf) ; 
    
    des_coupe_y(star.get_tggg()(), 0., nzdes, "Metric potential (NB-1) r sin\\gh",
		&surf) ; 
    
    des_coupe_y(star.get_ak_car()(), 0., nzdes, "A\\u2\\dK\\uij\\dK\\dij\\u",
		&surf) ;     
  }
  
  // Cleaning
  // --------
  
  delete peos ;    
  
  exit(EXIT_SUCCESS) ; 
  
  return EXIT_SUCCESS ; 
   
}
