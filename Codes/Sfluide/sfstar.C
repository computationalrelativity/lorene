/***********************************************************************
 *   calculate stationary configuration of superfluid neutron star
 ***********************************************************************/

// // headers C
#include <math.h>
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

void compare_analytic (Et_rot_bifluid& star); 
string get_file_base (double xp, double sig);
string get_file_base (double xp, double sig, double eps, double om1, double om2);


// some global EOS parameters 
double kap1, kap2, kap3, beta, m1, m2, detA, k1, k2, kk, R;
double sigma, eps, eps_n, xp;
double nc, rhoc;  // central density

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
      
      int nzdes = star.get_nzet() ; 
      
      // Cmp defining the surface of the star (via the density fields)
      // 
      Cmp surf(mp) ;
      surf = -0.2*star.get_nbar()()(0,0,0,0) ;
      surf.annule(0, star.get_nzet()-1) ;
      surf += star.get_nbar()() ; ;
      surf = prolonge_c1(surf, star.get_nzet()) ;

      Cmp surf2(mp) ;
      surf2 = -0.2*star.get_nbar2()()(0,0,0,0) ;
      surf2.annule(0, star.get_nzet()-1) ;
      surf2 += star.get_nbar2()() ; ;
      surf2 = prolonge_c1(surf2, star.get_nzet()) ;

      des_bi_coupe_y(star.get_nbar()(), 0., nzdes, "Fluid 1 baryonic density", 
		     &surf, &surf) ; 

      des_bi_coupe_y(star.get_nbar2()(), 0., nzdes, "Fluid 2 baryonic density", 
      	     &surf2, &surf2) ; 
      
      des_bi_coupe_y(star.get_logn()(), 0., nzdes, "Grav. potential", 
		     &surf, &surf2) ; 
    }

    // now print out key-values of the configuration in such a "translated" way
    // that we can compare the results to the analytic solution of PCA02:
    if (eos.identify() == 2)  // only do that if type = eos_bf_poly_newt
      compare_analytic (star);
 
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
compare_analytic (Et_rot_bifluid& star)
{
  bool is_static = false;

  Eos_bf_poly_newt eos = dynamic_cast<const Eos_bf_poly_newt&>(star.get_eos());

  // let's see if it's really what we called  "analytic" EOS in PCA02:
  if ( (eos.identify() != 2) || (eos.get_gam1() != 2) || (eos.get_gam2() != 2) ||
       (eos.get_gam3() != 1) || (eos.get_gam4() != 1) )
    {
      cout << "This EOS is not of type Newtonian analytic EOS, compare_analytic() useless here!\n";
      return;
    }
  else
    {
      cout << "\n\n----------------------------------------------------------------------" << endl;
      cout << " compare_analytic() called on AnalyticEOS: now comparing... " << endl << endl;
    }


  if (star.get_ent()()(0,0,0,0) != star.get_ent2()()(0,0,0,0) ) 
    cout << "\n!! WARNING !!: central chemical potentials differ..!!\n";


  if ( star.get_omega_c() == 0  && star.get_omega2() == 0 )
    is_static = true;

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

  fname = resdir + get_file_base (xp, sigma, eps, om_n, om_p);

// #define NUM_POINTS 100 		// number of points in density profiles
//   Tbl rr (NUM_POINTS);
//   rr.set_etat_qcq();

//   std::ofstream prof((fname+".prof").c_str());
//   prof << setprecision(9);
//   for (int i=0; i < NUM_POINTS; i++)
//     {
//       rr.set(i) = 1.0 * i / (NUM_POINTS-1);
//       prof << rr(i) << "\t" ;
      
//       prof << nn.val_point(R*rr(i), 0, 0) << "\t";
//       prof << np.val_point(R*rr(i), 0, 0) << "\t";

//       prof << nn.val_point(R*rr(i), M_PI/4, 0) << "\t";
//       prof << np.val_point(R*rr(i), M_PI/4, 0) << "\t";
      
//       prof << nn.val_point(R*rr(i), M_PI/2, 0) << "\t";
//       prof << np.val_point(R*rr(i), M_PI/2, 0) << "\t";
      
//       prof << endl;
//     }

  double mnat = rhoc * R * R * R;
  std::ofstream data((fname+".d").c_str());
  data << setprecision(16);

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

  data << "Om_n = " << om_n << endl;
  data << "Om_p = " << om_p << endl;
    
  data << "\n# Global stellar quantities:\n";
  data << "Mn = " << star.mass_b1() / mnat << endl;
  data << "Mp = " << star.mass_b2() / mnat << endl;

  data << "Rn = " << star.ray_pole()/R  << "\t" << star.ray_eq()/ R  << endl;
  data << "Rp = " << star.ray_pole2()/R << "\t" << star.ray_eq2()/ R << endl;

  data << "# Viriel identity violations:" << endl;
  data << "GRV3 = " << star.grv3() << endl;
  data << "GRV2 = " << star.grv2() << endl;

  return;

} //  compare_analytic

/*----------------------------------------------------------------------
 * get_file_base(): construct filename-base from EOS parameters
 *
 *----------------------------------------------------------------------*/
string
get_file_base (double xp0, double sig0)
{
  ostringstream s;
  char cbuf[50]; // man, <ios> does not seem to exist here... 

  sprintf (cbuf, "Newt_xp%4.2f_sig%4.2f", xp0, sig0);

  s << cbuf;
  
  return (s.str());
}

string
get_file_base (double xp0, double sig0, double eps0, double om1, double om2)
{
  ostringstream s;
  char cbuf[50]; // man, <ios> does not seem to exist here... 
  double relOm;

  if (om2 != 0)
    relOm = (om1 - om2)/om2;
  else
    relOm = 0;

  s << get_file_base (xp0, sig0);

  sprintf (cbuf, "_eps%4.2f_Om%8.6f_R%4.2f", eps0, om2, relOm);

  s << cbuf;
  
  return (s.str());
}
