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
 *   calculate stationary configuration of superfluid neutron star
 ***********************************************************************/

// // headers C
#include <math.h>
#include <gsl/gsl_integration.h>

#include <sstream>
// headers Lorene
#include "et_rot_bifluid.h"
#include "eos_bifluid.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"

#include "unites.h"

// Local prototype (for drawings only)
Cmp raccord_c1(const Cmp& uu, int l1) ; 

void compare_analytic (Et_rot_bifluid& star, int adapt); 
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
//----------------------------------------------------------------------

//******************************************************************************

int main(){

  using namespace Unites ; 

    // For the display : 
    char display_bold[]="x[1m" ; display_bold[0] = 27 ;

    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------
    bool relat, graph;
    int mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
	mermax_poisson, nz, nzet, nt, np; 
    double ent1_c, ent2_c, freq_si, freq2_si, precis, freq_ini_si;
    double freq2_ini_si,relax, relax_poisson ;  

    // adaptive-grid paramters + defaults
    int nzadapt = 0;
    double thres_adapt = 0.0;
    double precis_adapt = 1.e-14;

    char *parrot = "settings.par"; // config-file
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
      //      cout << "DEBUG: l= "<<l<< ", read in nr[l] = " << nr[l] << "; rmin[l] = " << bornes[l] << endl;
    }

    bornes[nz] = __infinity ;

    Tbl ent_limit(nzet) ;
    ent_limit.set_etat_qcq() ;
    ent_limit.set(nzet-1) = 0 ; 	// enthalpy at the stellar surface
    Tbl ent2_limit(ent_limit) ;
    ent2_limit.set_etat_qcq() ;
    ent2_limit.set(nzet-1) = 0 ;  // enthalpy 2 at the stellar surface

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

    // Read paramters specific to adaptive grid:
    // ----------------------------------------
    if ( (read_variable (NULL, "nzadapt", nzadapt) == 0) && (nzadapt > 0) )  // do we want adaptive grid?
      {
	res += read_variable (NULL, "thres_adapt", thres_adapt);
	res += read_variable (NULL, "precis_adapt", precis_adapt);
	
	if (res != 0)
	  cout << "WARNING: some adaptive-grid variables were not found! Using default ... \n";
      }

    // Read parameters specific to Kepler-limit
    // ----------------------------------------
    int kepler_fluid	= 0; 	// Kepler limit for which fluid? 1,2; 3 = both
    int kepler_wait_steps = 1; 	// how many steps after mer_fix_omega shall we start?
    double kepler_factor = 1.01; // factor to increase omega in each step to approach Kepler (>1!)

    res = 0;
    if( (read_variable (NULL, "kepler_fluid", kepler_fluid) == 0) && (kepler_fluid > 0) )
      {
	res += read_variable (NULL, "kepler_wait_steps", kepler_wait_steps);
	res += read_variable (NULL, "kepler_factor", kepler_factor);
	
	if (res != 0)
	  cout << "WARNING: some Kepler-limit paramters were not found in settings.par! Using default ... \n";
      }


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
    Eos_bifluid* peos = Eos_bifluid::eos_from_file("eos.par") ;
    Eos_bifluid& eos = *peos ;

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

    Itbl icontrol(8) ;
    icontrol.set_etat_qcq() ; 
    icontrol.set(0) = mer_max ; 
    icontrol.set(1) = mer_rot ; 
    icontrol.set(2) = mer_change_omega ; 
    icontrol.set(3) = mer_fix_omega ; 
    icontrol.set(4) = mermax_poisson ; 
    icontrol.set(5) = nzadapt;  	// nb of domains for adaptive grid
    icontrol.set(6) = kepler_fluid;  	// index of fluid for Kepler-search (0=none)
    icontrol.set(7) = kepler_wait_steps;

    Tbl control(8) ; 
    control.set_etat_qcq() ; 
    control.set(0) = precis ; 
    control.set(1) = omega_ini ;
    control.set(2) = omega2_ini ;
    control.set(3) = relax ; 
    control.set(4) = relax_poisson ; 
    control.set(5) = thres_adapt;
    control.set(6) = precis_adapt;
    control.set(7) = kepler_factor;


    Tbl diff(8) ;     

    star.equilibrium_bi(ent1_c, ent2_c, omega, omega2, ent_limit, 
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
    "   PARAMETERS USED FOR THE COMPUTATION (file settings.par) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    system("cat settings.par >> calcul.d") ; 

    fichfinal.open("calcul.d", ios::app) ;
    fichfinal << endl <<
    "================================================================" << endl ;
    fichfinal <<
    "	           EOS PARAMETERS (file eos.par) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    system("cat eos.par >> calcul.d") ;

    // Identification du code et de ses sous-routines (no. de version RCS) :     	
    fichfinal.open("calcul.d", ios::app) ; 
    fichfinal << endl <<
    "================================================================" << endl ; 
    fichfinal << "	    IDENTIFICATION OF THE CODE : " << endl ; 
    fichfinal << 
    "================================================================" << endl ; 
    fichfinal.close() ; 
    system("ident sfstar >> calcul.d") ; 


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
      char bslash[2] = {92, '\0'} ;  // 92 is the ASCII code for backslash 

      int nzdes = star.get_nzet() ; 
      
      // Cmp defining the surface of the star (via the density fields)
      // 
      Cmp surf(mp) ;
      surf = -0.2*star.get_nbar()()(0,0,0,0) ;
      surf.annule(0, star.get_nzet()-1) ;
      surf += star.get_nbar()() ; ;

      surf.std_base_scal();

//       des_profile(surf, 0, 1.5, M_PI/2, 0, "surf before prolonge");
      surf = prolonge_c1(surf, star.get_nzet()) ;
//       des_profile(surf, 0, 1.5, M_PI/2, 0, "surf after prolonge");

      Cmp surf2(mp) ;
      surf2 = -0.2*star.get_nbar2()()(0,0,0,0) ;
      surf2.annule(0, star.get_nzet()-1) ;
      surf2 += star.get_nbar2()() ; ;
      surf2 = prolonge_c1(surf2, star.get_nzet()) ;

      des_bi_coupe_y(star.get_nbar()(), 0., nzdes, "Fluid 1 baryonic density", &surf, &surf) ; 

      des_bi_coupe_y(star.get_nbar2()(), 0., nzdes, "Fluid 2 baryonic density", &surf2, &surf2) ; 
      
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
    //    if (eos.identify() == 2)  // only do that if type = eos_bf_poly_newt
    compare_analytic (star, nzadapt);

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
compare_analytic (Et_rot_bifluid& star, int adapt)
{
  using namespace Unites ; 

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
  Cmp nn = star.get_nbar()() / nc;
  Cmp np = star.get_nbar2()() / nc;

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
  string resdir = "Results/";
  string fname;

  fname = resdir + get_file_base (star.is_relativistic(), xp, sigma, eps, om_n, om_p, eos.get_typeos(), star.get_nzet(), adapt);


  Map_et &map = (Map_et&)(star.get_mp());
  const Mg3d* mg = map.get_mg() ;	// Multi-grid

  //----------------------------------------------------------------------
  // get radii at intermediate angle, ~pi/4

  const Coord& theta = map.tet ;
  int g = mg->get_nt(0)/2;   // theta close to pi/4

  double thetaI = (+theta)(0,0,g,0);
  double RnI = map.val_r_jk(star.l_surf()(0,g), star.xi_surf()(0,g), g, 0);
  double RpI = map.val_r_jk(star.l_surf2()(0,g), star.xi_surf2()(0,g), g, 0);

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

  // in order to uniquely idenfity the run, we append the output of "calcul.d" to this file:
  data << "\n======================================================================\n";
  data << "               Identification of the run: calcul.d\n";
  data << "======================================================================\n";
  data.close();

  system(("cat calcul.d >> "+fname+".d").c_str()) ; 


  return;

} //  compare_analytic

//----------------------------------------------------------------------
// function A(r) for numerical integration
double AofR (double r, void *params)
{
  double res;
  AofR_params *myparams = (AofR_params*)params;

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
  char *head;

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

