/*
 *   Copyright (c) 2002, 2003 Jerome Novak
 *   Copyright (c) 2003, 2004 Reinhard Prix
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

/***********************************************************************
 *   compute expected value of the glitch rise time using Eq. (21) of 
 *   Sourie, Chamel, Novak and Oertel, MNRAS 464, 4641 (2017)
 ***********************************************************************/

/* The input parameters are very similar to those of sfstar. Some new parameters
 * relative to the computation of the partial moments of inertia and the mutual friction 
 * coefficient should be set in glitch_rise.par.
 *
 * The partial moments of inertia I_{XY} are computed at a fixed total baryon mass M^B
 * with beta-equilibrium or at fixed partial baryon masses M^B_n and M^B_p. To get the values
 * of these masses corresponding to the gravitational mass under interest, it is recommended 
 * to start by using the sfstar code (with convergence towards a given grav. mass).
 */ 

// headers C
#include <cmath>
#include <gsl/gsl_integration.h>
#include <ctime>
#include <gsl/gsl_matrix.h>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

// headers Lorene
#include "et_rot_bifluid.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"


namespace Lorene{
  // Local prototype (for drawings only)
  Cmp raccord_c1(const Cmp& uu, int l1) ; 

  void compare_analytic (Et_rot_bifluid& star, int adapt, const char *resdir); 
  string get_file_base (bool relat, double xp, double sig);
  string get_file_base (bool relat, double xp, double sig, double eps, double om1, double om2, int typeos, int nzet, int adapt);

  // some global EOS parameters 
  double kap1, kap2, kap3, beta, m1, m2, detA, k1, k2, kk, R;
  double sigma, eps, eps_n, xp;
  double nc, rhoc;  // central density

  //----------------------------------------------------------------------
  // stuff for GSL-integration of A(r)
  typedef struct {
    Cmp *AofR;
    double theta;
    double phi; 
  } AofR_params;

  double AofR (double r, void *params);

}
//----------------------------------------------------------------------
using namespace Lorene ;
//******************************************************************************

int main(){
  
  time_t tbegin1,tend1;           // to measure computation time
  double texec1 =0.;
  using namespace Unites ; 
  
  //------------------------------------------------------------------
  //	    Parameters of the computation 
  //------------------------------------------------------------------
  bool relat, graph;
  int mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
    mermax_poisson, nz, nzet, nt, np, id_eos; 
  double ent1_c, ent2_c, freq_si, freq2_si, precis, freq_ini_si;
  double freq2_ini_si,relax, relax_poisson ;  

  // adaptive-grid paramters + defaults
  int nzadapt = 0;
  double thres_adapt = 0.0;
  double precis_adapt = 1.e-14;

  const char *parrot = "settings.par"; // config-file (input parameters for the configuration)
  int res = 0;

  res += read_variable (parrot, "relat", relat);
  res += read_variable (NULL, "ent1_c", ent1_c);
  res += read_variable (NULL, "ent2_c",ent2_c);
  res += read_variable (NULL, "freq_si", freq_si);
  res += read_variable (NULL, "freq2_si", freq2_si);
  res += read_variable (NULL, "mer_max", mer_max);
  res += read_variable (NULL, "precis", precis);
  res += read_variable (NULL, "mer_rot", mer_rot);
  res += read_variable (NULL, "freq_ini_si", freq_ini_si);
  res += read_variable (NULL, "freq2_ini_si", freq2_ini_si);
  res += read_variable (NULL, "mer_change_omega", mer_change_omega);
  res += read_variable (NULL, "mer_fix_omega", mer_fix_omega);
  res += read_variable (NULL, "relax", relax);
  res += read_variable (NULL, "mermax_poisson", mermax_poisson);
  res += read_variable (NULL, "relax_poisson", relax_poisson);
  res += read_variable (NULL, "graph", graph);
  res += read_variable (NULL, "nz", nz);
  res += read_variable (NULL, "nzet", nzet);
  res += read_variable (NULL, "nt", nt);
  res += read_variable (NULL, "np", np);
  res += read_variable (NULL, "ident", id_eos);
    
  if ( res != 0 )
    {
      cerr << "An error ocurred in reading the parameter file 'settings.par'. Terminating...\n";
      exit (-1);
    }
      
  int* nr = new int[nz];
  int* nt_tab = new int[nz];
  int* np_tab = new int[nz];
  double* bornes = new double[nz+1];
    
  char cbuf[50]; // man, <ios> does not seem to exist here... 
  for (int l=0; l<nz; l++) {
    sprintf (cbuf, "nr%d", l);
    res += read_variable (NULL, cbuf, nr[l]);
    sprintf (cbuf, "rmin%d", l);
    res += read_variable (NULL, cbuf, bornes[l]);
    np_tab[l] = np ; 
    nt_tab[l] = nt ; 
  }

  bornes[nz] = __infinity ;

  Tbl ent_limit(nzet) ;
  ent_limit.set_etat_qcq() ;
  ent_limit.set(nzet-1) = 0 ; 	// enthalpy at the stellar surface
  Tbl ent2_limit(ent_limit) ;
  ent2_limit.set_etat_qcq() ;
  ent2_limit.set(nzet-1) = 0 ; // enthalpy 2 at the stellar surface

  double tmp;
  for (int l=0; l<nzet-1; l++) {
    sprintf (cbuf, "ent_limit%d", l);
    res += read_variable (NULL, cbuf, tmp);
    ent_limit.set(l) = tmp;
    ent2_limit.set(l) = tmp;
  }

  if ( res != 0 )
    {
      cerr << "An error ocurred in reading the parameter file " << parrot <<". Terminating...\n";
      exit (-1);
    }

  // Read paramters specific to adaptive grid (optional)
  // --------------------------------------------------
  if ( (read_variable (NULL, "nzadapt", nzadapt) == 0) && (nzadapt > 0) )  // do we want adaptive grid?
    {
      res += read_variable (NULL, "thres_adapt", thres_adapt);
      res += read_variable (NULL, "precis_adapt", precis_adapt);
	
      if (res != 0)
	cout << "WARNING: some adaptive-grid variables were not found! Using default ... \n";
    }

  // Read parameters specific to Kepler-limit (optional)
  // --------------------------------------------------
  // Default values: 
  int kepler_fluid	= 0; 	// Kepler limit for which fluid? 1,2; 3 = both
  int kepler_wait_steps = 1; 	// how many steps after mer_fix_omega shall we start?
  double kepler_factor = 1.01; // factor to increase omega in each step to approach Kepler (>1!)

  res = 0;
  if( (read_variable (NULL, "kepler_fluid", kepler_fluid) == 0) && (kepler_fluid > 0) )
    {
      res += read_variable (NULL, "kepler_wait_steps", kepler_wait_steps);
      res += read_variable (NULL, "kepler_factor", kepler_factor);
	
      if (res != 0)
	cout << "WARNING: some Kepler-limit parameters were not found in settings.par! Using default ... \n";
    }

  // Read parameters specific to triaxial perturbation (optional)
  // ------------------------------------------------------------
  // Default values: 
  int mer_triax = 2000 ;       //step at which the 3-D perturbation is switched on
  double ampli_triax = 1.e-3 ; // relative amplitude of the 3-D perturbation

  res = 0;
  if( (read_variable (NULL, "mer_triax", mer_triax) == 0) && (mer_triax < mer_max) )
    {
      res += read_variable (NULL, "ampli_triax", ampli_triax) ;
      if (res != 0)
	cout << "WARNING: some 3D perturbation parameters were not found in settings.par! Using default ... \n";
    }

  // read parameters related to mass-convergence (optional)
  // -------------------------------------------------------
  int mer_mass = -1;	/* default: no mass-adaption */
  /* if mer_mass > 0, different convergence scenarios are possible : 
   * - if mbar2_wanted = -1, the code converges towards a configuration in which M_G (grav. mass) = mbar1_wanted assuming beta equilibrium
   * - if mbar2_wanted = 0 , the code converges towards a configuration in which M_B (total bar. mass) = mbar1_wanted assuming beta equilibrium
   * - otherwise, the code converges towards a configuration in which M_B^n = mbar1_wanted and M_B^p = mbar2_wanted (no beta equilibrium in this case)
   */
  double mbar1_wanted, mbar2_wanted, aexp_mass = 0.5;
  res = 0;
  if ( (read_variable (parrot, "mer_mass", mer_mass) == 0) && (mer_mass > 0) )
    {
      res += read_variable (NULL, "mbar1_wanted", mbar1_wanted);
      res += read_variable (NULL, "mbar2_wanted", mbar2_wanted);
	
      if (res != 0)
	{
	  cout << "ERROR: mer_mass > 0 but 'mbar1_wanted' or 'mbar2_wanted' not set!\n";
	  exit (-1);
	}
      // correct units of baryon-masses: (input is in Solar masses)
      mbar1_wanted *= msol;
      mbar2_wanted *= msol;
      read_variable (NULL, "aexp_mass", aexp_mass);

    } /* if fixed baryon masses */

  /* Directory to store results in */
  string resdir = "Results";	/* default value */
  char *resdir_tmp = NULL;
  if ( read_variable (NULL, "resdir", &resdir_tmp) == 0)
    resdir = resdir_tmp ;

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
    
  /* eos-params are read from different parameter files depending on the EoS type ! */
    
  Eos_bifluid* peos = 0x0 ;
  
  // Analytical EoSs
  if (id_eos == 1 || id_eos == 2 ) {

    const char *fpar_anal = "eos_anal.par"; 
    peos = Eos_bifluid::eos_from_file(fpar_anal);	

  }
  // Tabulated EoSs
  else if (id_eos == 3) {
    ifstream fpar_tab;
    fpar_tab.open("eos_tab.par") ;
    if ( !fpar_tab.good() ) {
      cerr << "Problem in opening the file eos_tab.par ! " << endl ;
      abort() ;
    }    
    peos = Eos_bifluid::eos_from_file(fpar_tab) ;
    fpar_tab.close() ;
  }
    
  Eos_bifluid& eos = *peos  ;
    
    
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
   
  cout << endl 
       << "==========================================================" << endl
       << "                    Physical parameters                   " << endl
       << "=========================================================="
       << endl ; 
  cout << endl ;

  cout << endl << "Equation of state : " 
       << endl << "=================   " << endl ;
  cout << eos << endl ; 

  cout << "Central enthalpy 1 : " << ent1_c << " c^2" << endl ; 
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
	 
  tbegin1=time(NULL); 
    
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

  //------------------------------------------------------------------
  //	Read extra parameters necessary to compute the glitch rise time
  //------------------------------------------------------------------
	 
  double coeff_B ; 		// mutual friction parameter
  double eps_freq ;		// epsilon parameter used to compute the partial moment of inertia 
  // using a fourth-order finite difference method

  const char *par_glitch = "glitch_rise.par"; // input parameters for the computation of the glitch rise time
  int res_gl = 0;

  res_gl += read_variable (par_glitch, "B", coeff_B);
  res_gl += read_variable (NULL, "eps_freq", eps_freq);
    
  if ( res_gl != 0 )
    {
      cerr << "An error ocurred in reading the parameter file 'glitch_rise.par'. Terminating...\n";
      exit (-1);
    }

  //-------------------------------
  // To store intermediate results 
  //-------------------------------
  	
  ofstream fich_der1 ; // Store quantites computed for fixed Omega_p
  fich_der1.open("coupling_omega_p_fixed.d", ios::app); 

  ofstream fich_der2 ; // Store quantites computed for fixed Omega_n
  fich_der2.open("coupling_omega_n_fixed.d", ios::app); 

  fich_der1.precision(12);
  fich_der1 << "# Omega_n[Hz]   Omega_p[Hz]      Mb(Msol)	 Mbn(Msol)	Mbp(Msol)	Mg(Msol)	"
	    << "  GRV2		GRV3	  J[10^38 kg m^2 s^-1]  Jn[10^38 kg m^2 s^-1]	Jp[10^38 kg m^2 s^-1]	"
	    << "  tilde{I}_n[10^38 kg m^2]   tilde{I}_p[10^38 kg m^2]   tilde{omega}_n[rad s^-1]   tilde{omega}_p[rad s^-1]    tilde{eps}_n   tilde{eps}_p  "
	    << endl;
		  
  fich_der2.precision(12);
  fich_der2 << "# Omega_n[Hz]   Omega_p[Hz]      Mb(Msol)	 Mbn(Msol)	Mbp(Msol)	Mg(Msol)	"
	    << "  GRV2		GRV3	  J[10^38 kg m^2 s^-1]  Jn[10^38 kg m^2 s^-1]	Jp[10^38 kg m^2 s^-1]	"
	    << "  tilde{I}_n[10^38 kg m^2]   tilde{I}_p[10^38 kg m^2]   tilde{omega}_n[rad s^-1]   tilde{omega}_p[rad s^-1]    tilde{eps}_n   tilde{eps}_p  "
	    << endl;

  double Jn_minus_2eps_N = 0, Jn_minus_eps_N = 0, Jn_plus_eps_N = 0, Jn_plus_2eps_N  = 0;
  double Jn_minus_2eps_P = 0, Jn_minus_eps_P = 0, Jn_plus_eps_P  = 0, Jn_plus_2eps_P  = 0;
  double Jp_minus_2eps_N = 0, Jp_minus_eps_N = 0, Jp_plus_eps_N = 0, Jp_plus_2eps_N  = 0;
  double Jp_minus_2eps_P = 0, Jp_minus_eps_P = 0, Jp_plus_eps_P = 0, Jp_plus_2eps_P = 0 ;	
  double J_minus_2eps_N = 0,  J_minus_eps_N = 0,  J_plus_eps_N = 0,  J_plus_2eps_N  = 0;
  double J_minus_2eps_P = 0,  J_minus_eps_P = 0,  J_plus_eps_P = 0,  J_plus_2eps_P = 0 ;
 
  double tilde_omegan_minus_2eps_N = 0, tilde_omegan_minus_eps_N = 0, tilde_omegan_plus_eps_N = 0, tilde_omegan_plus_2eps_N = 0 ;
  double tilde_omegan_minus_2eps_P = 0, tilde_omegan_minus_eps_P = 0, tilde_omegan_plus_eps_P = 0, tilde_omegan_plus_2eps_P = 0 ;
  double tilde_omegap_minus_2eps_N = 0, tilde_omegap_minus_eps_N = 0, tilde_omegap_plus_eps_N = 0, tilde_omegap_plus_2eps_N = 0 ;
  double tilde_omegap_minus_2eps_P = 0, tilde_omegap_minus_eps_P = 0, tilde_omegap_plus_eps_P = 0, tilde_omegap_plus_2eps_P = 0 ;	

  // loop on the frequencies
  //-------------------------
	
  double freq_si_0 = freq_si ;
  double freq2_si_0 = freq2_si ;
      
  for (int i_freq = 1; i_freq <=9; i_freq++){

    freq_si = freq_si_0 ;
    freq2_si = freq2_si_0 ;	
	
    if (i_freq==1){ freq_si = freq_si_0* (1.+ 2.* eps_freq );}
    if (i_freq==2){ freq_si = freq_si_0* (1.+ eps_freq );}
    if (i_freq==3){ freq_si = freq_si_0* (1.- eps_freq );}
    if (i_freq==4){ freq_si = freq_si_0* (1.- 2.* eps_freq );}
    if (i_freq==5){ freq2_si = freq2_si_0* (1.+ 2.* eps_freq );}
    if (i_freq==6){ freq2_si = freq2_si_0* (1.+ eps_freq );}
    if (i_freq==7){ freq2_si = freq2_si_0* (1.- eps_freq );}
    if (i_freq==8){ freq2_si = freq2_si_0* (1.- 2.* eps_freq );}
	 	 

    //-----------------------------------------------------------------------
    //		Initialization of the enthalpy field
    //-----------------------------------------------------------------------

   
    const Coord& r = mp.r ;
    double ray0 = mp.val_r(nzet-1, 1., 0., 0.) ;  
    Cmp enta(mp) ; 
    enta = ent1_c * ( 1 - r*r / (ray0*ray0) ) ; 
    enta.annule(nz-1) ; 
    enta.std_base_scal() ; 
    Cmp entb = enta * ent2_c/ent1_c ; 
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

    Itbl icontrol(9) ;
    icontrol.set_etat_qcq() ; 
    icontrol.set(0) = mer_max ; 
    icontrol.set(1) = mer_rot ; 
    icontrol.set(2) = mer_change_omega ; 
    icontrol.set(3) = mer_fix_omega ; 
    icontrol.set(4) = mermax_poisson ; 
    icontrol.set(5) = nzadapt;  				// nb of domains for adaptive grid
    icontrol.set(6) = kepler_fluid;  		// index of fluid for Kepler-search (0=none)
    icontrol.set(7) = kepler_wait_steps;
    icontrol.set(8) = mer_triax ; 			// starting of triaxial perturbations

    Tbl control(9) ; 
    control.set_etat_qcq() ; 
    control.set(0) = precis ; 
    control.set(1) = omega_ini ;
    control.set(2) = omega2_ini ;
    control.set(3) = relax ; 
    control.set(4) = relax_poisson ; 
    control.set(5) = thres_adapt;
    control.set(6) = precis_adapt;
    control.set(7) = kepler_factor;
    control.set(8) = ampli_triax ;

    Tbl diff(9) ;     

    star.equilibrium_bi(ent1_c, ent2_c, omega, omega2, ent_limit, 
			ent2_limit, icontrol, control, diff, 
			mer_mass, mbar1_wanted, mbar2_wanted, aexp_mass);

     
    cout << endl << "Final star : " 
	 << endl << "==========   " << endl ;

    cout.precision(10) ; // modify output precision
    cout << star << endl ; 
    cout << "Total mass of neutron-fluid: " << star.mass_b1() / msol << " Msol\n";
    cout << "Total mass of proton-fluid: " << star.mass_b2() / msol << " Msol\n";
    cout << "Total baryon mass: " << star.mass_b()/msol << "Msol\n";
    double vit_triax = diff(8) ;
    cout << "Growing rate of triaxial perturbation: " << vit_triax 
	 << endl ; 

    //-----------------------------------------------
    //  General features of the final configuration
    //  saved in a file
    //-----------------------------------------------

    /*
      ostringstream nom_calcul ; 
      nom_calcul << "result_" ;
      nom_calcul << ent2_c ;
      nom_calcul << "_" ;
      nom_calcul << freq_si ;
      nom_calcul << ".d" ;
      string nom_fich = nom_calcul.str();
      ofstream fichfinal(nom_fich.c_str(), ios::out) ; */
    ofstream fichfinal("result.d") ;
    
    fichfinal.precision(16) ; 
    
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

    // total central densities
    /*   	m1 = eos.get_m1();
		m2 = eos.get_m2();
		double rhon_c = m1 * star.get_nbar()()(0,0,0,0);  
		double rhop_c = m2 * star.get_nbar2()()(0,0,0,0);
    */
    double rhon_c = star.get_nbar()()(0,0,0,0);  
    double rhop_c = star.get_nbar2()()(0,0,0,0);

    fichfinal << "\n\nCentral neutron density: " << rhon_c << " x 0.1 fm^-3 \n";
    fichfinal << "Central proton density: " << rhop_c << " x 0.1 fm^-3 \n";
    fichfinal << "Total central baryon density : " << rhon_c + rhop_c << " x 0.1 fm^-3 \n\n";

    fichfinal << endl <<
      "================================================================" << endl ;
    fichfinal <<
      "   PARAMETERS USED FOR THE COMPUTATION (file settings.par) : " << endl ;
    fichfinal <<
      "================================================================" << endl ;
    fichfinal.close() ;
	
    system("cat settings.par >> result.d") ; 
    fichfinal.open("result.d", ios::app) ;
    /* system(("cat settings.par >> "+nom_fich).c_str()) ;
       fichfinal.open(nom_fich.c_str(), ios::app) ; */
   
    fichfinal << endl <<
      "================================================================" << endl ; 
    fichfinal << "	    IDENTIFICATION OF THE CODE : " << endl ; 
    fichfinal << 
      "================================================================" << endl ; 
    fichfinal.close() ; 

    system("ident sfstar >> result.d") ; 
    // system(("ident result.d >> "+nom_fich).c_str()) ;

    tend1=time(NULL);
    	
    // Saveguard of the whole configuration
    // ------------------------------------
    
    FILE* fresu = fopen("resu.d", "w") ;
    
    star.get_mp().get_mg()->sauve(fresu) ;			// writing of the grid
    star.get_mp().sauve(fresu) ;                	// writing of the mapping
    star.get_eos().sauve(fresu) ;  					// writing of the EOS
    star.sauve(fresu) ;                         	// writing of the star
    
    fclose(fresu) ;
    

    // Store intermediate results :

    // freq_p = freq2_si_0 (fixed)
    if (i_freq==1){ // freq_n = freq_si_0 * (1 + 2 * eps_freq) 
      J_plus_2eps_N 	= star.angu_mom()	; 
      Jn_plus_2eps_N = star.angu_mom_1() ;
      Jp_plus_2eps_N = star.angu_mom_2() ;
      tilde_omegan_plus_2eps_N = star.coupling_LT_1()/ star.coupling_mominert_1() ;
      tilde_omegap_plus_2eps_N = star.coupling_LT_2()/ star.coupling_mominert_2() ;
    }
    if (i_freq==2){ // freq_n = freq_si_0 * (1 + eps_freq) 
      J_plus_eps_N  = star.angu_mom()	; 
      Jn_plus_eps_N = star.angu_mom_1() ;
      Jp_plus_eps_N = star.angu_mom_2() ;
      tilde_omegan_plus_eps_N = star.coupling_LT_1()/ star.coupling_mominert_1() ;
      tilde_omegap_plus_eps_N = star.coupling_LT_2()/ star.coupling_mominert_2() ;
    }
    if (i_freq==3){ // freq_n = freq_si_0 * (1 - eps_freq)
      J_minus_eps_N 	= star.angu_mom()	; 
      Jn_minus_eps_N = star.angu_mom_1() ;
      Jp_minus_eps_N = star.angu_mom_2() ;
      tilde_omegan_minus_eps_N = star.coupling_LT_1()/ star.coupling_mominert_1() ;
      tilde_omegap_minus_eps_N = star.coupling_LT_2()/ star.coupling_mominert_2() ;
    }
    if (i_freq==4){ // freq_n = freq_si_0 * (1 - 2 * eps_freq)
      J_minus_2eps_N  = star.angu_mom()	; 
      Jn_minus_2eps_N = star.angu_mom_1() ;
      Jp_minus_2eps_N = star.angu_mom_2() ;
      tilde_omegan_minus_2eps_N = star.coupling_LT_1()/ star.coupling_mominert_1() ;
      tilde_omegap_minus_2eps_N = star.coupling_LT_2()/ star.coupling_mominert_2() ;
    }		
    // freq_n = freq_si_0 (fixed)
    if (i_freq==5){ // freq_p = freq2_si_0 * (1 + 2 * eps_freq) 
      J_plus_2eps_P 	= star.angu_mom()	; 
      Jn_plus_2eps_P = star.angu_mom_1() ;
      Jp_plus_2eps_P = star.angu_mom_2() ;
      tilde_omegan_plus_2eps_P = star.coupling_LT_1()/ star.coupling_mominert_1() ;
      tilde_omegap_plus_2eps_P = star.coupling_LT_2()/ star.coupling_mominert_2() ;
    }
    if (i_freq==6){ // freq_p = freq2_si_0 * (1 + eps_freq) 
      J_plus_eps_P  = star.angu_mom()	; 
      Jn_plus_eps_P = star.angu_mom_1() ;
      Jp_plus_eps_P = star.angu_mom_2() ;
      tilde_omegan_plus_eps_P = star.coupling_LT_1()/ star.coupling_mominert_1() ;
      tilde_omegap_plus_eps_P = star.coupling_LT_2()/ star.coupling_mominert_2() ;
    }
    if (i_freq==7){ // freq_p = freq2_si_0 * (1 - eps_freq) 
      J_minus_eps_P 	= star.angu_mom()	; 
      Jn_minus_eps_P = star.angu_mom_1() ;
      Jp_minus_eps_P = star.angu_mom_2() ;
      tilde_omegan_minus_eps_P = star.coupling_LT_1()/ star.coupling_mominert_1() ;
      tilde_omegap_minus_eps_P = star.coupling_LT_2()/ star.coupling_mominert_2() ;
    }
    if (i_freq==8){ // freq_p = freq2_si_0 * (1 - 2 * eps_freq) 
      J_minus_2eps_P 	= star.angu_mom()	; 
      Jn_minus_2eps_P = star.angu_mom_1() ;
      Jp_minus_2eps_P = star.angu_mom_2() ;				
      tilde_omegan_minus_2eps_P = star.coupling_LT_1()/ star.coupling_mominert_1() ;
      tilde_omegap_minus_2eps_P = star.coupling_LT_2()/ star.coupling_mominert_2() ;
    }

    if ( (i_freq == 9) || (i_freq == 1) || (i_freq == 2) || (i_freq == 3) || (i_freq == 4) ) { 
      fich_der1  << freq_si  << "     " 
		 << freq2_si << "     "	
		 << star.mass_b() / msol		<< "     "
		 << star.mass_b1()/ msol 		<< "     "
		 << star.mass_b2() /msol 		<< "     "
		 << star.mass_g() /msol 		<< "     " 
		 << star.grv2() 			<< "     "
		 << star.grv3() 			<< "     "
		 << star.angu_mom() * rho_unit * pow(r_unit, double(4.)) * v_unit / double(1.e38) << "     " 
		 << star.angu_mom_1()  * rho_unit * pow(r_unit, double(4.)) * v_unit / double(1.e38) << "     " // Jn   
		 << star.angu_mom_2()	 * rho_unit * pow(r_unit, double(4.))	* v_unit / double(1.e38) << "     " // Jp  	  
		 << star.coupling_mominert_1()* rho_unit * pow(r_unit, double(5.)) / double(1.e38) << "     " //\tilde{I}_n
		 << star.coupling_mominert_2()	* rho_unit * pow(r_unit, double(5.)) / double(1.e38) << "     " //\tilde{I}_p 
		 << star.coupling_LT_1()/ star.coupling_mominert_1() * f_unit		<< "     " //\tilde{omega}_n  
		 << star.coupling_LT_2()/ star.coupling_mominert_2() * f_unit		<< "     " //\tilde{omega}_p  	
		 << star.coupling_entr()/ star.coupling_mominert_1()					<< "     " //\tilde{eps}_n
		 << star.coupling_entr()/ star.coupling_mominert_2()					<< "     " //\tilde{eps}_p
		 << endl;
    }
    
    if ( (i_freq== 9) || (i_freq == 5) || (i_freq == 6) || (i_freq == 7) || (i_freq == 8) ) { 
      fich_der2  << freq_si				<< "     " 
		 << freq2_si 				<< "     "	
		 << star.mass_b() / msol		<< "     "
		 << star.mass_b1()/ msol 		<< "     "
		 << star.mass_b2() /msol 		<< "     "
		 << star.mass_g() /msol			<< "     " 
		 << star.grv2() 			<< "     "
		 << star.grv3()				<< "     "
		 << star.angu_mom() * rho_unit * pow(r_unit, double(4.))	* v_unit / double(1.e38)	<< "     " 
		 << star.angu_mom_1()  * rho_unit * pow(r_unit, double(4.))	* v_unit / double(1.e38)	<< "     " // Jn   
		 << star.angu_mom_2()	 * rho_unit * pow(r_unit, double(4.))	* v_unit / double(1.e38)	<< "     " // Jp
		 << star.coupling_mominert_1()	* rho_unit * pow(r_unit, double(5.)) / double(1.e38) 	<< "     " //\tilde{I}_n
		 << star.coupling_mominert_2() * rho_unit * pow(r_unit, double(5.)) / double(1.e38) 	<< "     " //\tilde{I}_p 
		 << star.coupling_LT_1()/ star.coupling_mominert_1() * f_unit		<< "     " //\tilde{omega}_n  
		 << star.coupling_LT_2()/ star.coupling_mominert_2() * f_unit		<< "     " //\tilde{omega}_p  	
		 << star.coupling_entr()/ star.coupling_mominert_1() 			<< "     " //\tilde{eps}_n
		 << star.coupling_entr()/ star.coupling_mominert_2() 			<< "     " //\tilde{eps}_p
		 << endl;
    }

		
  }  // end of the loop on the frequencies

		
  fich_der1.close() ;
  fich_der2.close() ;

		
  // Drawings
  // --------
    	
  // To study the Newtonian limit:
  double unsurc2 ;
  if ( star.is_relativistic() ) { 
    unsurc2 = 1;     
  }
  else { 
    unsurc2 = 0.; 
  }


	  
  const Map_radial *map_colloc = dynamic_cast<const Map_radial*>(&star.get_mp());
  const Mg3d* mg_colloc = map_colloc->get_mg() ;	// Multi-grid
  const Coord& r_colloc = map_colloc->r;
   
  // Computation of the vorticity profile (and the mutual friction moment)
  // ---------------------------------------------------------------------
  if ( star.get_omega_c() != 0 ) {  	
    // if omega_n = 0 -> no vorticity
 
    // Scalar fields :
    Cmp nbar1 = star.get_nbar()(); 							// n_n
    Cmp nbar2 = star.get_nbar2()();							// n_p
    Cmp Knp = star.get_K_np()() ;  							// K^np
    Cmp Knn = star.get_K_nn()() ;   							// K^nn
    Cmp ent1 = star.get_ent()();								// H^n
    Cmp ent2 = star.get_ent2()();								// H^p
    Cmp mu1  =  eos.get_m1() * exp ( unsurc2 * ent1) ; //mu^n
    mu1.std_base_scal() ;

    Cmp Gamman = star.get_gam_euler()(); 					// Gamma_n
    Cmp Gammap = star.get_gam_euler2()();					// Gamma_p
    Cmp N_metric = star.get_nnn()() ;						// N 
    Cmp B2_metric = star.get_b_car()();						// B^2
    Cmp A2_metric = star.get_a_car()() ; 					// A^2
    Cmp omega_metric = star.get_nphi()() ;					// omega
    Cmp omegaT_metric = star.get_tnphi()();				// omega * r * sin(theta)
    Cmp delta_car 	= star.get_delta_car()(); 				// relative speed^2
    Cmp GammaDelta 	= 1./sqrt(1 - unsurc2 * delta_car) ;	// Relative Lorentz factor	
 
    Cmp g_tt = - N_metric * N_metric +  B2_metric * omegaT_metric * omegaT_metric ; 		// g_tt
    Cmp g_phph = B2_metric;																					// g_phiphi
    g_phph.std_base_scal();
    g_phph.mult_rsint() ;
    g_phph.mult_rsint() ;   // NB : wrong values in the last domain but this domain is unused in the vorticity calculation
    Cmp g_tph = - B2_metric * omega_metric ;
    g_tph.std_base_scal();
    g_tph.mult_rsint() ;
    g_tph.mult_rsint() ;  
		
    // Components of the neutron momentum : p^n_t and p^n_phi
    double omega_n = star.get_omega_c() ;
    double omega_p = star.get_omega2() ;
    Cmp p_t =  g_tt * ( mu1 * Gamman / N_metric  +   Knp * nbar2 / N_metric * (Gammap - GammaDelta * Gamman) ) 
      + g_tph * ( mu1 * Gamman / N_metric * omega_n  +   Knp * nbar2 / N_metric * (Gammap * omega_p - GammaDelta * Gamman *  omega_n ) ) ; 
    Cmp p_ph =  g_tph * ( mu1 * Gamman / N_metric  +   Knp * nbar2 / N_metric * (Gammap - GammaDelta * Gamman) ) 
      + g_phph * ( mu1 * Gamman / N_metric * omega_n  +   Knp * nbar2 / N_metric * (Gammap * omega_p - GammaDelta * Gamman *  omega_n ) ) ;   
	  	
    //derivatives :
    Cmp dp_tsdr = p_t.dsdr() ;			// d p^n_t / dr
    Cmp dp_phsdr = p_ph.dsdr() ;		// d p^n_phi / dr
    Cmp dp_tsdth = p_t.srdsdt() ;		// 1/r * d p^n_t / dtheta
    Cmp dp_phsdth = p_ph.srdsdt() ;	// 1/r * d p^n_phi / dtheta
 
    dp_tsdr.dec2_dzpuis ();
    dp_phsdr.dec2_dzpuis ();
    dp_tsdth.dec2_dzpuis ();
    dp_phsdth.dec2_dzpuis () ;
 
    // To treat properly the origin (r=0)
    Cmp num = ( N_metric * N_metric * ( dp_phsdr  * dp_phsdr  + dp_phsdth * dp_phsdth )		
		-  g_phph * (  pow(omega_metric * dp_phsdr +   dp_tsdr, 2.)  + pow( omega_metric * dp_phsdth  + dp_tsdth,2.)  ) ) ;
    num.div_rsint() ;
    num.div_rsint() ;
    Cmp Vorticity = num / ( B2_metric *  N_metric * N_metric * A2_metric); 
    Vorticity = sqrt(Vorticity) ; 
    Vorticity.std_base_scal();
		
    Cmp num2 = ( N_metric * N_metric * pow( dp_phsdr  * dp_phsdr  + dp_phsdth * dp_phsdth , 2. ) 
		 - g_phph * pow( omega_metric * (  dp_phsdr *  dp_phsdr + dp_phsdth * dp_phsdth ) + dp_phsdr * dp_tsdr + dp_tsdth *dp_phsdth, 2.) );
    num2.std_base_scal();
    num2.div_rsint() ;
    num2.div_rsint() ;	
    Cmp h_perp_2 = num2 	/ ( B2_metric  *  N_metric * N_metric * A2_metric * A2_metric * pow(Vorticity, 4.) ); 
    h_perp_2.std_base_scal();	 

    // The neutron vorticity should vanish in areas where the neutron fluid is absent: 
    for (int i=0 ; i<nz; i++){
      for (int l=0; l<np; l++){
	for (int k = 0; k<nt; k++){
	  for (int j=0; j< nr[i]; j++){
	    if (nbar1(i,l,k,j) <= 0) { Vorticity.set(i,l,k,j) =0. ; h_perp_2.set(i,l,k,j) =0.;}
	  }
	}
      }
    }
  
    // to plot the vorticity profiles and make some checks (collocation points)

    ofstream fichdes("prof_vorticity.d") ;
    fichdes.precision(16) ; 
    fichdes 		<< "#  r [km](theta=0)  	   r [km](theta=pi/2)  	 Newtonian_Vorticity	   "
			<< " Vorticity(theta=0)	   	h_perp_2(theta=0)  		Vorticity(theta=pi/2)	   	h_perp_2(theta=pi/2) " 
			<< endl ;	
    for (int i=0 ; i<nz; i++){
      for (int j=0; j< nr[i]; j++){
	int pi_2 = mg_colloc->get_nt(i)-1;	
	double position_0 = (+r_colloc)(i,0,0,j) ;
	double position_pi_2 =  (+r_colloc)(i,0,pi_2,j) ;
 	 
	// theta = 0
	double vorticity_0	= 	Vorticity(i,0,0,j);
	double hperp_2_0		= 	h_perp_2(i,0, 0, j) ;
	// theta = pi/2	
	double g_phph_pi2		= 	g_phph(i,0,pi_2,j); 
	double dp_tsdr_pi2 	= 	dp_tsdr(i,0,pi_2,j); 
	double dp_phsdr_pi2  =	dp_phsdr(i,0,pi_2,j); 
	double dp_tsdth_pi2  =	dp_tsdth(i,0,pi_2,j); 
	double dp_phsdth_pi2 =	dp_phsdth(i,0,pi_2,j); 	
	double vorticity_pi2	= 	Vorticity(i,0,pi_2,j);
	double hperp_2_pi2 	= 	h_perp_2(i,0, pi_2, j) ;
	double N_metric_pi2 		=	N_metric(i,0, pi_2, j) ;  
	//double B2_metric_pi2 	= 	B2_metric(i,0, pi_2, j) ;
	double A2_metric_pi2 	=	A2_metric(i,0, pi_2, j) ; 
	double omega_metric_pi2 =  omega_metric(i,0, pi_2, j) ; 
		 
	double vorticity_verif =  ( N_metric_pi2 * N_metric_pi2 * ( dp_phsdr_pi2  * dp_phsdr_pi2 ) 
				    -  g_phph_pi2 * (   (omega_metric_pi2 * dp_phsdr_pi2  + dp_tsdr_pi2) * (omega_metric_pi2 * dp_phsdr_pi2  + dp_tsdr_pi2) )) / 
	  ( A2_metric_pi2 * N_metric_pi2 * N_metric_pi2 * g_phph_pi2 )    ;
	vorticity_verif = sqrt(vorticity_verif) ; // normalisation constant = *f_unit  * rhonuc_si * 1e-44

	double hperp_2_verif =  N_metric_pi2 * N_metric_pi2 * pow( dp_phsdr_pi2  * dp_phsdr_pi2 +  dp_phsdth_pi2  * dp_phsdth_pi2, 2. ) 
	  -  g_phph_pi2 * pow(   omega_metric_pi2 * (dp_phsdr_pi2  * dp_phsdr_pi2   + dp_phsdth_pi2 * dp_phsdth_pi2 )  
				 + dp_phsdth_pi2 * dp_tsdth_pi2 + dp_phsdr_pi2  *  dp_tsdr_pi2, 2. )  ;
	hperp_2_verif = hperp_2_verif / (  A2_metric_pi2 *  A2_metric_pi2  * N_metric_pi2 * N_metric_pi2 * g_phph_pi2  * pow ( vorticity_verif, 4.) ) ;
	
	fichdes 	<< position_0 * 10 		<< "   " 
			<< position_pi_2 * 10 	<< "   " 
			<< 2.* omega_n * eos.get_m1()  << "   " // Newtonian expression for the vorticity (in the limit of corotation)
			<< vorticity_0						<< "   " 
			<< hperp_2_0						<< "   " 
			<< vorticity_pi2					<< "   " 
			<< hperp_2_pi2 					<< "   " 	
			<<	endl ;
 		
      }
    }
    fichdes.close() ; 
     

    /*
     * Computation of the partial moments of inertia and 
     * every quantity needed to compute the rise time
     * Note that the notations are similar to those
     * used in Sourie et al., MNRAS 464 (2017)
     */
     

    // the partial moments of inertia are computed from a fourth-order finite difference method

    double Inn = 1./(12. * eps_freq *  star.get_omega_c()) * (Jn_minus_2eps_N - 8.* Jn_minus_eps_N + 8. * Jn_plus_eps_N - Jn_plus_2eps_N );
    double Inp = 1./(12. * eps_freq *  star.get_omega2()) * (Jn_minus_2eps_P - 8.* Jn_minus_eps_P + 8. * Jn_plus_eps_P - Jn_plus_2eps_P );
    double Ipn = 1./(12. * eps_freq *  star.get_omega_c()) * (Jp_minus_2eps_N - 8.* Jp_minus_eps_N + 8. * Jp_plus_eps_N - Jp_plus_2eps_N );
    double Ipp = 1./(12. * eps_freq *  star.get_omega2()) * (Jp_minus_2eps_P - 8.* Jp_minus_eps_P + 8. * Jp_plus_eps_P - Jp_plus_2eps_P );
    double hat_In = 1./(12. * eps_freq *  star.get_omega_c()) * (J_minus_2eps_N - 8.* J_minus_eps_N + 8. * J_plus_eps_N - J_plus_2eps_N ); // = Inn + Inp
    double hat_Ip = 1./(12. * eps_freq *  star.get_omega2()) * (J_minus_2eps_P - 8.* J_minus_eps_P + 8. * J_plus_eps_P - J_plus_2eps_P ); // = Ipp + Ipn
    double hat_I = hat_In + hat_Ip ;

    // computation of the mutual friction torque (see Eq. 6 of Sourie, Chamel, Novak & Oertel, MNRAS, 2017) : 

    Cmp dens(mp) ; 
    dens = Gamman * nbar1 * Vorticity * h_perp_2 * A2_metric * sqrt(B2_metric)  ;
    dens.std_base_scal();
    //double Gamma_int = dens.integrale() * coeff_B * (omega_p - omega_n)  ; 
		
    // computation of the kappa and zeta quantities (see Eqs. 18 and 23 of Sourie, Chamel, Novak & Oertel, MNRAS, 2017 for definition) : 

    double kappa = dens.integrale() ;
    double zeta = kappa/ ( 2. * hat_In *  star.get_omega_c() ) ; 

    // to compute the theoretical expression of the rise time 
    // (see Eq. 21 of Sourie, Chamel, Novak & Oertel, MNRAS, 2017 for definition) :

    double tau_r = (Inn * Ipp - Inp * Ipn ) / ( kappa * hat_I * coeff_B ) ;

      
    // to compute the eps_XY^{LT} quantities 
    // (see Eq. C9 of Sourie, Chamel, Novak & Oertel, MNRAS, 2017 for definition) :
    // be careful, these quantities (as the \tilde{I}_n, \tilde{\omega}_n and \tilde{\epsilon}_n
    // only make sense in the slow-rotation approximation and at first order in the lag.
	
    double eps_nn_LT  = 1./(12. * eps_freq * star.get_omega_c()) * (tilde_omegan_minus_2eps_N - 8.* tilde_omegan_minus_eps_N + 8. * tilde_omegan_plus_eps_N - tilde_omegan_plus_2eps_N );
    double eps_pn_LT 	= 1./(12. * eps_freq * star.get_omega2())  * (tilde_omegan_minus_2eps_P - 8.* tilde_omegan_minus_eps_P + 8. * tilde_omegan_plus_eps_P - tilde_omegan_plus_2eps_P );
    double eps_np_LT  = 1./(12. * eps_freq * star.get_omega_c()) * (tilde_omegap_minus_2eps_N - 8.* tilde_omegap_minus_eps_N + 8. * tilde_omegap_plus_eps_N - tilde_omegap_plus_2eps_N );
    double eps_pp_LT  = 1./(12. * eps_freq * star.get_omega2())  * (tilde_omegap_minus_2eps_P - 8.* tilde_omegap_minus_eps_P + 8. * tilde_omegap_plus_eps_P - tilde_omegap_plus_2eps_P );

	
    // coupling coefficients 
    // (see Eq. 29 of Sourie, Chamel, Novak & Oertel, MNRAS, 2017 for definition) :

    double hat_eps_n = ( star.coupling_entr()/ star.coupling_mominert_1() - eps_pn_LT ) / (1. -  eps_pn_LT  - eps_nn_LT ) ;
    double hat_eps_p = ( star.coupling_entr()/ star.coupling_mominert_2() - eps_np_LT ) / (1. -  eps_np_LT  - eps_pp_LT ) ;

    // other expression for the rise time 
    // (see Eq. 30 of Sourie, Chamel, Novak & Oertel, MNRAS, 2017)

    double tau_r_bis = hat_Ip / hat_I * (1. - hat_eps_n - hat_eps_p) / ( kappa * coeff_B) * hat_In;
	
    ofstream fich_glitch ;
    fich_glitch.open("tau_rise.d", ios:: app); 
    fich_glitch.precision(12) ;
    fich_glitch 	<< " Omega_n(Hz)		Omega_p(Hz)		H^n_c(c^2)    H^p_c(c^2) 		"
			<< " M^B(Msol)			M^B_n(Msol)		M^B_p(Msol)   M_G(Msol)    	"
			<< " Rcirc(km) 		Rcirc2(km) 		Rmean(km) 	  Rmean2(km)      GRV2		    GRV3		"
			<< " Inn[10^38 kg m^2]  Inp[10^38 kg m^2]  Ipn[10^38 kg m^2]  Ipp[10^38 kg m^2] 	 "
			<< " hat_In[10^38 kg m^2]  hat_Ip[10^38 kg m^2]  hat_I[10^38 kg m^2] 		" 
			<< " kappa[10^38 kg m^2 s^-1]	   zeta       tau_r(s) " 
			<< " tilde{I}_n[10^38 kg m^2]		tilde{I}_p[10^38 kg m^2]    tilde{eps}_n	tilde{eps_p} "
			<< " eps_nn^LT		eps_pn^LT		eps_np^LT		eps_pp^LT	" 
			<< " hat_eps_n   hat_eps_p  tau_r_bis(s) " 
			<< endl;

    fich_glitch	<< freq_si 							<< "     "
		<< freq2_si 						<< "     "
		<< ent1_c 							<< "     "
		<< ent2_c 							<< "     "
		<< star.mass_b() /msol 			<< "     " 
		<< star.mass_b1()/ msol 		<< "     "
		<< star.mass_b2() /msol 		<< "     "
		<< star.mass_g() / msol			<< "     "
		<< star.r_circ() / km   		<< "     "
		<< star.r_circ2() / km  		<< "     "
		<< star.mean_radius() / km 	<< "     "
		<< star.mean_radius2() / km 	<< "     "
		<< star.grv2() 					<< "     "
		<< star.grv3() 					<< "     "
		<< Inn * rho_unit * (pow(r_unit, double(5.)) / double(1.e38) )			<< "     "
		<< Inp * rho_unit * (pow(r_unit, double(5.)) / double(1.e38) )			<< "     "
		<< Ipn * rho_unit * (pow(r_unit, double(5.)) / double(1.e38) )		   << "     "
		<< Ipp * rho_unit * (pow(r_unit, double(5.)) / double(1.e38) )			<< "     "
		<< hat_In * rho_unit * (pow(r_unit, double(5.)) / double(1.e38) )		<< "     "
		<< hat_Ip * rho_unit * (pow(r_unit, double(5.)) / double(1.e38) )		<< "     "
		<< hat_I  * rho_unit * (pow(r_unit, double(5.)) / double(1.e38) )		<< "     "			
		<< kappa	* rho_unit  * pow(r_unit, double(5.)) / t_unit	/ double(1.e38) 				<< "     " 
		<< zeta 										<< "     "
		<< tau_r* t_unit							<< "     "
		<< star.coupling_mominert_1()* rho_unit * (pow(r_unit, double(5.)) / double(1.e38) )						<< "     " //\tilde{I}_n
		<< star.coupling_mominert_2()* rho_unit * (pow(r_unit, double(5.)) / double(1.e38) )						<< "     " //\tilde{I}_p 
		<< star.coupling_entr()/ star.coupling_mominert_1()		<< "     " //\tilde{eps}_n
		<< star.coupling_entr()/ star.coupling_mominert_2()		<< "     " //\tilde{eps}_p
		<< eps_nn_LT		<< "    " 
		<< eps_pn_LT		<< "    " 
		<< eps_np_LT		<< "    " 
		<< eps_pp_LT		<< "    " 
		<< hat_eps_n		<< "    " 
		<< hat_eps_p		<< "    " 
		<< tau_r_bis * t_unit 		<< "    " 
		<< endl ;
    fich_glitch.close() ;

  } // end of the if condition on Omega_n



  texec1=difftime(tend1,tbegin1);
  cout << " Time of computation of the star : " << texec1 << endl;
  //-----------------------------------------------------------------//
    
  // Cleaning
  // --------
  delete [] nr ; 
  delete [] nt_tab ; 
  delete [] np_tab ; 
  delete [] type_r ; 
  delete [] bornes ; 
    

  if (graph == 1) {
    char title[80] ;
    char bslash[2] = {92, '\0'} ;  // 92 is the ASCII code for backslash 
  
    int nzdes = star.get_nzet() ; 
        
    // Cmp defining the surface of the star (via the density fields)
    Cmp surf(mp) ;
    surf = -0.2*star.get_nbar()()(0,0,0,0) ;
    surf.annule(0, star.get_nzet()-1) ;
    surf += star.get_nbar()() ; ; 
    surf.std_base_scal();
    //	     des_profile(surf, 0, 1.5, M_PI/2, 0, "surf before prolonge");
    surf = prolonge_c1(surf, star.get_nzet()) ;
    //	     des_profile(surf, 0, 1.5, M_PI/2, 0, "surf after prolonge");
  
    Cmp surf2(mp) ;
    surf2 = -0.2*star.get_nbar2()()(0,0,0,0) ;
    surf2.annule(0, star.get_nzet()-1) ;
    surf2 += star.get_nbar2()() ; ;
    surf2 = prolonge_c1(surf2, star.get_nzet()) ;
  
    des_bi_coupe_y(star.get_nbar()(), 0., nzdes, "Fluid 1 baryonic density", &surf, &surf) ; 
    //      des_bi_coupe_y(star.get_nbar()(), 0., nzdes, "Fluid 1 baryonic density", &surf, &surf, 1.2, true, 15, 5000, 5000) ; 
    des_bi_coupe_y(star.get_nbar2()(), 0., nzdes, "Fluid 2 baryonic density", &surf2, &surf2) ; 
    des_bi_coupe_y(star.get_press()(), 0., nzdes, "Pressure", &surf, &surf2) ;
    des_bi_coupe_y(star.get_logn()(), 0., nzdes, "Grav. potential", &surf, &surf2) ; 
  
    strcpy(title, "Azimuthal shift N") ; 
    strcat(title, bslash) ; 
    strcat(title, "u") ; 
    strcat(title, bslash) ; 
    strcat(title, "gf") ; 
    des_bi_coupe_y(star.get_nphi()(), 0., nzdes, title, &surf, &surf2) ; 
  	
    strcpy(title, "Metric potential ") ; 
    strcat(title, bslash) ; 
    strcat(title, "gz") ; 
    des_bi_coupe_y(star.get_dzeta()(), 0., nzdes, title, &surf, &surf2) ; 
  
  }
 
  // now print out key-values of the configuration in such a "translated" way
  // that we can compare the results to the analytic solution of PCA02:
  /*
    if (eos.identify() == 2)  // only do that if type = eos_bf_poly_newt
    compare_analytic (star, nzadapt, resdir);
  */

  // Cleaning
  // --------
  delete peos ;    
                    
  exit(EXIT_SUCCESS) ; 
  return EXIT_SUCCESS ; 
     
} // main()




// ----------------------------------------------------------------------
// compare_analytic()
// print out appropriate "translations" of parameters and results such that
// we can compare them to the analytic solution of PCA02
//----------------------------------------------------------------------
void 
compare_analytic (Et_rot_bifluid& star, int adapt, const char *resdir)
{
  using namespace Unites ;
  using Lorene::beta ;

  Eos_bf_poly eos = dynamic_cast<const Eos_bf_poly&>(star.get_eos());

  //   // let's see if it's really what we called  "analytic" EOS in PCA02:
  //   if ( (eos.identify() != 2) || (eos.get_gam1() != 2) || (eos.get_gam2() != 2) ||
  //        (eos.get_gam3() != 1) || (eos.get_gam4() != 1) )
  //     {
  //       cout << "This EOS is not of type Newtonian analytic EOS, compare_analytic() useless here!\n";
  //       return;
  //     }
  //   else
  //     {
  //       cout << "\n\n----------------------------------------------------------------------" << endl;
  //       cout << " compare_analytic() called on AnalyticEOS: now comparing... " << endl << endl;
  //     }

  double muc1 = star.get_ent()()(0,0,0,0);
  double muc2 = star.get_ent2()()(0,0,0,0);
  cout.precision (15);
  if (muc1 != muc2)
    cout << "\n!! WARNING !!: central chemical potentials differ..!!\n: mu1 = " << muc1 << "; mu2 = " << muc2 << endl;;

  // get "raw" EOS parameters
  kap1 = eos.get_kap1();
  kap2 = eos.get_kap2();
  kap3 = eos.get_kap3();
  beta = eos.get_beta();
  m1 = eos.get_m1();
  m2 = eos.get_m2();
  
  cout << "Raw EOS parameters: kappa1 = " <<  kap1 << " kappa2 = " << kap2; 
  cout << " kappa3 = " << kap3 << " beta = " << beta << endl;


  // Central densities
  double nn_c = star.get_nbar()()(0,0,0,0);
  double np_c = star.get_nbar2()()(0,0,0,0);
  nc = nn_c + np_c;
  rhoc = m1 * nn_c + m2 * np_c;

  // baryon densities in natural units  
  //  Cmp nn = star.get_nbar()() / nc;
  //  Cmp np = star.get_nbar2()() / nc;

  // translate EOS parameters into x_p, sigma and epsilon  
  detA = kap1*kap2 - kap3*kap3;
  k1 = m1 * (kap2 - kap3) / detA;
  k2 = m2 * (kap1 - kap3) / detA;
  kk = m1*k1 + m2* k2;
  R = M_PI /sqrt(qpig * kk) ;  // analytic prediction
  
  sigma = - kap3 / kap1;
  if ( fabs(sigma) < 1e-9 ) sigma = 0;
  eps = 2.0 * beta * nn_c / m2;
  eps_n = 2.0 * beta * np_c / m1;

  xp =  k2 / (k1 + k2);

  double xp_num = np_c / (nn_c + np_c) ;

  cout << setprecision(9);  
  cout << "Translated EOS parameters: sigma = " << sigma << ", epsilon_c = " << eps << ", xp = " << xp << endl;
  
  cout << "Central neutron density: " << nn_c << " rho_nuc" << endl;
  cout << "Central proton density: " << np_c << " rho_nuc" << endl;
  cout << "Central baryon density: " << nc << " rho_nuc" << endl;
  double rel = 2* g_si * star.mass_b()*m_unit / (c_si*c_si * R * r_unit);
  cout << "Relativity parameter: " << rel << endl;

  double om0 = sqrt ( 4.0 * M_PI * g_si * rhoc * rho_unit );
  cout << "Rotation-rate Unit: " << om0 << endl;
  double om_n = star.get_omega_c() * f_unit / om0;
  double om_p = star.get_omega2() * f_unit / om0;
  cout << "Natural rotation rates: Om_n = " << om_n << "; Om_p = " << om_p << endl;
  cout << "Analytic static radius: " << R << endl;
  cout << "eps_p(0) = " << eps << "  eps_n(0) = " << eps_n << endl;
  if ( eps >= 1.0 || eps_n >= 1.0 )
    cout << "*** WARNING **** negative effective masses if eps, eps_n >= 1 !! \n";

  // ******************************
  // now start with tests and output

  cout << "Total mass of neutron-fluid: " << star.mass_b1() / msol << " Msol\n";
  cout << "Total mass of proton-fluid: " << star.mass_b2() / msol << " Msol\n";
  cout << "Total baryon mass: " << star.mass_b()/msol << "Msol\n";

  // save some data about this configuration: 
  string fname;

  // make sure Results-directory exists, otherwise create it:
  DIR *thisdir;
  if ( (thisdir = opendir(resdir)) == NULL)
    if (mkdir (resdir, S_IREAD|S_IWRITE|S_IEXEC) == -1)
      {
	cout << "Failed to create results-dir: " << resdir << endl;
	exit (-1);
      }
  closedir (thisdir);
      
  // construct path+name for results-file
  string path = resdir;
  path += "/";
  fname = path + get_file_base (star.is_relativistic(), xp, sigma, eps, om_n, om_p, eos.get_typeos(), star.get_nzet(), adapt);


  const Map_radial *map = dynamic_cast<const Map_radial*>(&star.get_mp());
  const Mg3d* mg = map->get_mg() ;	// Multi-grid

  //----------------------------------------------------------------------
  // get radii at intermediate angle, ~pi/4

  const Coord& theta = map->tet ;
  int g = mg->get_nt(0)/2;   // theta close to pi/4

  double thetaI = (+theta)(0,0,g,0);
  double RnI = map->val_r_jk(star.l_surf()(0,g), star.xi_surf()(0,g), g, 0);
  double RpI = map->val_r_jk(star.l_surf2()(0,g), star.xi_surf2()(0,g), g, 0);

  cout << "theta = " << thetaI << "; RnI = " << RnI << "; RpI = " << RpI << endl;

  double mnat = rhoc * R * R * R;

  //----------------------------------------------------------------------
  // calculate magnitude of shift-vector = N^phi B r sin(th)
  Cmp shiftMag = star.get_tnphi()(); 	// N^phi r sin(th)
  shiftMag *= star.get_bbb()();


  //----------------------------------------------------------------------
  // calculate proper radii!
  Cmp aaa = sqrt( star.get_a_car()());		// a = sqrt(a^2)
  aaa.std_base_scal();
  // ok, not too sure about the spectral inner workings of this, so we
  // integrate this "by hand"... (i.e using GSL...)
  double rmax;
  double abserr;
  size_t neval;
  double RR_pol_n, RR_pol_p, RR_eq_n, RR_eq_p;	// the proper radii


  AofR_params params;
  params.AofR = &aaa;
  params.phi = 0;

  gsl_function func_AofR;
  func_AofR.function = AofR;
  func_AofR.params = &params;

  // RR_pol_n
  rmax = star.ray_pole();	/* upper limit of integration */  
  params.theta = 0;		/* pole */
  if (gsl_integration_qng (&func_AofR, 0, rmax, 0, 1e-7, &RR_pol_n, &abserr, &neval) ) 
    {
      cout << "Integration of proper radius RR_pol_n failed!" << endl;
      exit (-1);
    }
  cout << "Proper-radius RR_pol_n = "<<RR_pol_n<< ", abserr= "<<abserr<<", neval= "<<neval <<endl;
  // RR_pol_p
  rmax = star.ray_pole2();	/* upper limit of integration */  
  params.theta = 0;		/* pole */
  if (gsl_integration_qng (&func_AofR, 0, rmax, 0, 1e-7, &RR_pol_p, &abserr, &neval) ) 
    {
      cout << "Integration of proper radius RR_pol_p failed!" << endl;
      exit (-1);
    }
  cout << "Proper-radius RR_pol_p = "<<RR_pol_p<< ", abserr= "<<abserr<<", neval= "<<neval <<endl;
  // RR_eq_n
  rmax = star.ray_eq();		/* upper limit of integration */  
  params.theta = M_PI/2;	/* equator */
  if (gsl_integration_qng (&func_AofR, 0, rmax, 0, 1e-7, &RR_eq_n, &abserr, &neval) ) 
    {
      cout << "Integration of proper radius RR_eq_n failed!" << endl;
      exit (-1);
    }
  cout << "Proper-radius RR_eq_n = "<<RR_eq_n<< ", abserr= "<<abserr<<", neval= "<<neval <<endl;
  // RR_eq_p
  rmax = star.ray_eq2();	/* upper limit of integration */  
  params.theta = M_PI/2;	/* equator */
  if (gsl_integration_qng (&func_AofR, 0, rmax, 0, 1e-7, &RR_eq_p, &abserr, &neval) ) 
    {
      cout << "Integration of proper radius RR_eq_p failed!" << endl;
      exit (-1);
    }
  cout << "Proper-radius RR_eq_p = "<<RR_eq_p<< ", abserr= "<<abserr<<", neval= "<<neval <<endl;


  cout << "Flattening of neutrons, r_pole/r_eq = " << RR_pol_n / RR_eq_n << endl;
  cout << "Flattening of protons,  r_pole/r_eq = " << RR_pol_p / RR_eq_p << endl;

  //----------------------------------------------------------------------

  std::ofstream data((fname+".d").c_str());
  data << setprecision(17);

  data << "# EOS and stellar parameters: \n";

  data << "kappa1 = " <<  kap1 << endl;
  data << "kappa2 = " <<  kap2 << endl; 
  data << "kappa3 = " <<  kap3 << endl;
  data << "beta = "   <<  beta << endl;

  data << "muc_n = " <<  star.get_ent()()(0,0,0,0) << endl;
  data << "muc_p = " << star.get_ent2()()(0,0,0,0) << endl;

  data << "rhoc = " << nc << " rho_nuc" << endl;
  data << "Om0 = " << om0 << endl;
  data << "Relativity-parameter = " << rel << endl;

  data << "\n# in natural units:\n";
  data << "sigma = "   << sigma << endl;
  data << "epsilon = " << eps   << endl;
  data << "xp = "      << xp    << endl;
  data << "xp_num = "  << xp_num << endl;

  data << "Om_n = " << om_n << endl;
  data << "Om_p = " << om_p << endl;
    
  data << "\n# Global stellar quantities:\n";
  data << "Mn = " << star.mass_b1() / mnat << endl;
  data << "Mp = " << star.mass_b2() / mnat << endl;

  data << "Rn = " << star.ray_pole()/R  << "\t" << star.ray_eq()/ R  << endl;
  data << "Rp = " << star.ray_pole2()/R << "\t" << star.ray_eq2()/ R << endl;

  data << "\n# Intermediate radii: \n";
  data << "thetaI = " << thetaI << endl;
  data << "RXI = " << RnI/R << "\t" << RpI/R  << endl;

  data << "\n# Viriel identity violations:" << endl;
  data << "GRV3 = " << star.grv3() << endl;
  data << "GRV2 = " << star.grv2() << endl;

  data << "\n# relativistic stuff in physical units: \n";
  data << "Mbar_n = " << star.mass_b1() / msol << " Msol\n";
  data << "Mbar_p = " << star.mass_b2() / msol << " Msol\n";
  data << "Mbar = " << star.mass_b() / msol << " Msol\n";
  data << "Mgrav = " << star.mass_g() / msol << " Msol\n";
  data << "Rcirc_n = " << star.r_circ() * 10.0 << " km\n";
  data << "Rcirc_p = " << star.r_circ2() * 10.0 << " km\n";
  data << "lapse N(0)      = " << star.get_nnn()()(0,0,0,0) << endl;
  data << "lapse N(eq)     = " << star.get_nnn()().va.val_point(0,1,M_PI/2,0) << endl;
  data << "lapse N(pol)    = " << star.get_nnn()().va.val_point(0,1,0, 0) << endl;
  data << "|N^phi|(eq)     = " << shiftMag.va.val_point(0,1,M_PI/2,0) << endl;

  data << "RR_eq_n         = " << RR_eq_n << endl;
  data << "RR_pol_n        = " << RR_pol_n << endl;
  data << "RR_eq_p         = " << RR_eq_p << endl;
  data << "RR_pol_p        = " << RR_pol_p << endl;

  // in order to uniquely idenfity the run, we append the output of "result.d" to this file:
  data << "\n======================================================================\n";
  data << "               Identification of the run: result.d\n";
  data << "======================================================================\n";
  data.close();

  system(("cat result.d >> "+fname+".d").c_str()) ; 


  return;

} //  compare_analytic



namespace Lorene {

  //----------------------------------------------------------------------
  // function A(r) for numerical integration
  double AofR (double r, void *params)
  {
    double res;
    AofR_params *myparams = reinterpret_cast<AofR_params*>(params);

    res = myparams->AofR->val_point( r, myparams->theta, myparams->phi);

    return (res);
  
  } /* AofR() */
  //----------------------------------------------------------------------

  /*----------------------------------------------------------------------
   * get_file_base(): construct filename-base from EOS parameters
   *
   *----------------------------------------------------------------------*/
  string
  get_file_base (bool relat, double xp0, double sig0)
  {
    ostringstream s;
    char cbuf[50]; // man, <ios> does not seem to exist here... 
    const char *head;

    if (relat)
      head = "Rel";
    else
      head = "Newt";

    sprintf (cbuf, "%s_xp%4.2f_sig%4.2f", head, xp0, sig0);

    s << cbuf;
  
    return (s.str());
  }

  string
  get_file_base (bool relat, double xp0, double sig0, double eps0, double om1, double om2, int typeos, int nzet, int adapt)
  {
    ostringstream s;
    char cbuf[50]; // man, <ios> does not seem to exist here... 
    double relOm;

    if (om2 != 0)
      relOm = (om1 - om2)/om2;
    else
      relOm = 0;

    s << get_file_base (relat, xp0, sig0);

    sprintf (cbuf, "_eps%4.2f_Om%8.6f_R%4.2f%s%s%s", eps0, om2, relOm, 
	     (typeos==5)? "_sr" : "",
	     (nzet==2)? "_2dom" : "",
	     (adapt>0)? "_adapt" : "");

    s << cbuf;
  
    return (s.str());
  }

}
