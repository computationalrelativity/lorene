/*
 * Main code for computation of equilibrium configuration of a NS-BH binary system.
 *
 */

/*
 *   Copyright (c) 2002  Philippe Grandclement, Keisuke Taniguchi,
 *              Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char coal_ns_bh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2004/03/25 12:35:35  j_novak
 * now using namespace Unites
 *
 * Revision 1.3  2003/11/29 07:01:25  k_taniguchi
 * Large modification. However, the code has not worked yet.
 *
 * Revision 1.2  2002/12/19 14:58:20  e_gourgoulhon
 * First non-empty version.
 *
 * Revision 1.1  2002/12/18 10:33:10  e_gourgoulhon
 * Computations of NS - BH binaries
 *
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>
#include <math.h>

// Lorene headers
#include "bin_ns_bh.h"
#include "nbr_spx.h"
#include "eos.h"
#include "utilitaires.h"
#include "unites.h"

int main(){
  
  using namespace Unites ;

    // Identification of all the subroutines called by the code :

    // system("ident coal_ns_bh > identif.d") ;

    // For the display :
    char display_bold[]="x[1m" ; display_bold[0] = 27 ;
    char display_normal[] = "x[0m" ; display_normal[0] = 27 ;

    //------------------------------------------------------------------
    //	    Parameters of the computation
    //------------------------------------------------------------------

    char blabla[120], nomini[80] ;
    int mermax, mermax_eqb, prompt, graph, fmer_stop, fmer_save, mermax_poisson ;
    int mermax_potvit, mer_masse, fmer_upd_met, ind_rel_met ;
    double seuil, relax_poisson, relax_potvit, relax, aexp_masse ;
    double mbar_voulue, mirr_voulue, fact_separ, relax_met, relax_omeg ;
    double fact_omeg_min, fact_omeg_max, thres_adapt ;

    ifstream fpar("par_coal.d") ;
    fpar.getline(blabla, 120) ;
    fpar.getline(blabla, 120) ;
    fpar.getline(nomini, 80) ;
    fpar >> fact_separ ; fpar.getline(blabla, 120);
    fpar >> mbar_voulue ; fpar.getline(blabla, 120) ;
    fpar >> mirr_voulue ; fpar.getline(blabla, 120) ;
    mbar_voulue *= msol ;
    mirr_voulue *= msol ;
    fpar.getline(blabla, 120) ;
    fpar >> mermax ; fpar.getline(blabla, 120) ;
    fpar >> relax ; fpar.getline(blabla, 120) ;
    fpar >> mermax_eqb ; fpar.getline(blabla, 120) ;
    fpar >> prompt ; fpar.getline(blabla, 120) ;
    fpar >> graph ; fpar.getline(blabla, 120) ;
    fpar >> seuil ; fpar.getline(blabla, 120) ;
    fpar >> fmer_stop ; fpar.getline(blabla, 120) ;
    fpar >> fmer_save ; fpar.getline(blabla, 120) ;
    fpar >> mermax_poisson ; fpar.getline(blabla, 120) ;
    fpar >> relax_poisson ; fpar.getline(blabla, 120) ;
    fpar >> mermax_potvit ; fpar.getline(blabla, 120) ;
    fpar >> relax_potvit ; fpar.getline(blabla, 120) ;
    fpar >> mer_masse ; fpar.getline(blabla, 120) ;
    fpar >> aexp_masse ; fpar.getline(blabla, 120) ;
    fpar >> fmer_upd_met ; fpar.getline(blabla, 120);
    fpar >> ind_rel_met ; fpar.getline(blabla, 120);
    fpar >> relax_met ; fpar.getline(blabla, 120);
    if (ind_rel_met == 0) relax_met = 1. ;
    fpar >> relax_omeg ; fpar.getline(blabla, 120);
    fpar >> fact_omeg_min ; fpar.getline(blabla, 120);
    fpar >> fact_omeg_max ; fpar.getline(blabla, 120);
    fpar >> thres_adapt ; fpar.getline(blabla, 120);
    fpar.close() ;

    cout << endl
	 << "==========================================================" << endl
	 << "                    Physical parameters                   " << endl
	 << "=========================================================="
	 << endl ;
    cout << endl << endl ;
    cout << "File containing the initial conditions : " << nomini << endl ;
    cout << "Factor by which the initial separation will be multiplied : "
	 << fact_separ << endl ;
    if ( abs(mer_masse) < mermax ) {
	cout << "Baryon mass required for the NS [M_sol] : "
	     << mbar_voulue / msol << endl ;
	cout << "Irreducible mass required for the BH  [M_sol] : "
	     << mirr_voulue / msol << endl ;
    }
    cout << endl
	 << "==========================================================" << endl
	 << "              Parameters of the computation               " << endl
	 << "=========================================================="
	 << endl ;
    cout << "Maximum number of steps in the main iteration : "
	 << mermax << endl ;
    cout << "Relaxation factor in the main iteration  : "
	 << relax << endl ;
    cout << "Maximum number of steps in Et_bin_nsbh::equilibrium_nsbh : "
	 << mermax_eqb << endl ;
    cout << "Threshold on the enthalpy relative change for ending the computation : "
	 << seuil << endl ;
    cout << "Step interval between safeguards of the whole configuration  : "
	 << fmer_save << endl ;
    cout << "Maximum number of steps in Map_et::poisson : "
	 << mermax_poisson << endl ;
    cout << "Relaxation factor in Map_et::poisson : "
	 << relax_poisson << endl ;
    cout << "Maximum number of steps in Map_radial::poisson_compact : "
	 << mermax_potvit << endl ;
    cout << "Relaxation factor in Map_radial::poisson_compact : "
	 << relax_potvit << endl ;
    cout << "Step from which the baryon mass is forced to converge : "
	 << mer_masse << endl ;
    cout << "Exponent for the increase factor of the central enthalpy : "
	 << aexp_masse << endl ;
    cout << "Step interval between metric updates : "
	 << fmer_upd_met << endl ;
    if (ind_rel_met == 1) {
	cout << "Relaxation factor of the metric : "
	     << relax_met << endl ;
    }
    else {
	cout << "No relaxation on the metric" << endl ;
    }
    cout << "Relaxation factor on Omega (orbital angular velocity) : "
	 << relax_omeg << endl ;
    cout << "Relative low bound in the omega search :  "
	 << fact_omeg_min << endl ;
    cout << "Relative high bound in the omega search : "
	 << fact_omeg_max << endl ;
    cout <<
    "Threshold on |dH/dr|_eq / |dH/dr|_pole for the adaptation of the mapping for the NS"
    << endl << thres_adapt << endl ;

    arrete(prompt) ;

    //------------------------------------------------------------------
    //	    Read of the initial conditions
    //------------------------------------------------------------------

    FILE* fich = fopen(nomini, "r") ;

    int mer_ini ;
    fread(&mer_ini, sizeof(int), 1, fich) ;

    Mg3d mg_ns(fich) ;
    Map_et mp_ns(mg_ns, fich) ;
    Eos* peos = Eos::eos_from_file(fich) ;

    Mg3d mg_bh(fich) ;
    Map_af mp_bh(mg_bh, fich) ;

    cout << "Make an object of Bin_ns_bh" << endl ;
    arrete(prompt) ;
    Bin_ns_bh bibi(mp_ns, *peos, mp_bh, fich) ;

    fclose(fich) ;

    //------------------------------------------------------------------
    //	    Modification of the separation between the two stars
    //------------------------------------------------------------------

    cout << "Modification of the orbital separation" << endl ;
    arrete(prompt) ;
    double ori_x = ((bibi.get_ns()).get_mp()).get_ori_x() ;
    ori_x *= fact_separ ;
    ((bibi.set_ns()).set_mp()).set_ori(ori_x, 0., 0.) ;

    ori_x = ((bibi.get_bh()).get_mp()).get_ori_x() ;
    ori_x *= fact_separ ;
    ((bibi.set_bh()).set_mp()).set_ori(ori_x, 0., 0.) ;


    //------------------------------------------------------------------
    //	    Update of the initial conditions
    //------------------------------------------------------------------

    // Initialisation of A^2, N, etc...
    // ---------------------------------
    cout << "Update the initial conditions for NS" << endl ;
    arrete(prompt) ;
    (bibi.set_ns()).update_metric( bibi.get_bh() ) ;
    // ????? arguments ?????
    cout << "Update the initial conditions for BH" << endl ;
    arrete(prompt) ;
    (bibi.set_bh()).update_metric( bibi.get_ns(), 1., 1.) ; // double, double

    // Initialisation of gradients of companion potentials
    // ---------------------------------------------------
    (bibi.set_ns()).update_metric_der_comp( bibi.get_bh() ) ;

    (bibi.set_bh()).fait_n_comp( bibi.get_ns() ) ;
    (bibi.set_bh()).fait_psi_comp( bibi.get_ns() ) ;

    bibi.fait_decouple() ;
    bibi.fait_tkij() ;

    // Initialisation of hydro quantities for NS
    // -----------------------------------------
    (bibi.set_ns()).equation_of_state() ;
    (bibi.set_ns()).kinematics(bibi.get_omega(), bibi.get_x_axe()) ;
    (bibi.set_ns()).fait_d_psi() ;

    // Save in a file
    // --------------
    FILE* fresu = fopen("resu.d", "w") ;

    int mer1 = 0 ;
    fwrite(&mer1, sizeof(int), 1, fresu) ;	// mer

    (bibi.get_ns()).get_mp().get_mg()->sauve(fresu) ;
    (bibi.get_ns()).get_mp().sauve(fresu) ;
    (bibi.get_ns()).get_eos().sauve(fresu) ;

    (bibi.get_bh()).get_mp().get_mg()->sauve(fresu) ;
    (bibi.get_bh()).get_mp().sauve(fresu) ;

    bibi.sauve(fresu) ;

    fclose(fresu) ;

    cout << endl
	 << "======================================================" << endl
	 << "                  Initial conditions                  " << endl
	 << "======================================================" << endl ;
    cout << bibi << endl ;

    if ( (bibi.get_ns()).is_relativistic() ) {
      cout << "========================" << endl ;
      cout << "Relativistic computation" << endl ;
      cout << "========================" << endl ;
    }
    else {
      cout << "===========================================" << endl ;
      cout << "!!! WARNING : NS should be relativistic !!!" << endl ;
      cout << "===========================================" << endl ;
      abort() ;
    }

    arrete(prompt) ;


    //----------------------------------------------------------------
    //			Auxiliary quantities
    //----------------------------------------------------------------

    double ent_c ;     // Central enthaply in NS
    double dentdx ;    // Central d/dx(enthalpy) in NS

    // Resizing factor of the first shell
    // Computation is shown just befor "equilibrium" of NS
    Tbl fact_resize[] = {Tbl(1)} ;
    fact_resize[0].set_etat_qcq() ;

    // Error indicators in NS
    Tbl differ[] = {Tbl(7)} ;
    differ[0].set_etat_qcq() ;

    ent_c = (bibi.get_ns()).get_ent()()(0, 0, 0, 0) ;
    differ[0].set(0) = 1 ;	    // diff_ent = 1
    differ[0].set(1) = 1 ;	    // err_psi = 1
    differ[0].set(2) = 1 ;	    // err_n = 1
    differ[0].set(3) = 1 ;	    // err_confpsi = 1
    differ[0].set(4) = 1 ;	    // err_shift_x = 1
    differ[0].set(5) = 1 ;	    // err_shift_y = 1
    differ[0].set(6) = 1 ;	    // err_shift_z = 1

    double relax_jm1 = 1. - relax ;
    double relax_omeg_jm1 = 1. - relax_omeg ;

    //----------------------------------------------------------------
    //	 Binary system at the previous step (for the relaxation)
    //----------------------------------------------------------------

    Bin_ns_bh bibi_jm1 = bibi ;

    double omega_jm1 = bibi_jm1.get_omega() ;

    // n_comp and pot_centri of NS are initialized to 0.5 and 0 on bibi_jm1
    // --------------------------------------------------------------------
    bibi_jm1.set_ns().set_n_comp() = 0.5 ;
    bibi_jm1.set_ns().set_pot_centri() = 0 ;

    //---------------------------
    //	 Openning of log files
    //---------------------------

    ofstream fichresu("resglob.d") ;
    fichresu.precision(16) ;

    ofstream fichrota("resrota.d") ;
    fichrota.precision(16) ;

    ofstream fichconv("resconv.d") ;
    fichconv.precision(16) ;

    ofstream fichet("resstar.d") ;
    fichet.precision(16) ;

    fichrota <<
      "#      Omega [rad/s]        x_axe [km]        x_g (NS) [km]        x_g (BH) [km]        M_grav [M_sol]        J [G M_sol^2/c]"
	     << endl ;

    fichconv <<
      "#      diff_ent        err_psi        err_n        err_confpsi        err_shift_x        err_shift_y        err_shift_z"
	     << endl ;

    fichet <<
      "#      ori_x [km]        ent_c        M_bar [M_sol]        R(theta=0) [km]        R(pi/2, 0) [km]        R(pi/2, pi/2) [km]        R(pi/2, pi) [km]"
	   << endl ;

    //    double omega_kep ;
    int mer ;

//=========================================================================
//		Start of iteration
//=========================================================================

    for (mer=0; (differ[0](0) > seuil) && (mer < mermax); mer++) {

      cout <<
	"================================================================"
	   << endl ;
      cout << "step = " << mer << "        diff_ent : "
	   << differ[0](0) << endl ;
      cout <<
	"================================================================"
	   << endl ;

      fichresu << mer;
      fichresu << "    step" << endl ;

      fichrota << mer ;
      fichconv << mer ;
      fichet << mer ;

      //------------------------------------------------
      //     Computation of the metric coefficients
      //------------------------------------------------

      if ( (mer % fmer_upd_met) == 0 ) {

	(bibi.set_ns()).update_metric(bibi.get_bh()) ;

	(bibi.set_ns()).update_metric_der_comp(bibi.get_bh()) ;

      }

      //------------------------------------------------------------
      //     Computation of the orbital angular velocity Omega
      //------------------------------------------------------------

      bibi.orbit_omega(fact_omeg_min, fact_omeg_max) ;

      // Translation of the stars in order to set the origin
      //  of the absolute frame on the rotation axis
      //-----------------------------------------------------
      // !!! The subroutine for computation of the position of ratation axis
      // !!!  has not been developped

      cout << display_bold << "New orbital velocity Omega : "
	   << bibi.get_omega() * f_unit << " rad/s" << display_normal << endl ;

      // Keplerian velocity (for comparison only)
      // ----------------------------------------
      /*
      omega_kep = 0. ;   // I do not know whether it is necessary or not.

      cout << "``Keplerian'' velocity (for comparison only) : "
	   << omega_kep * f_unit << " rad/s" << endl ;
      */
      /*
      cout << "New X coordinate of the rotation axis : "
	   << bibi.get_x_axe() / km << " km" << endl ;
      */

      arrete(prompt) ;

      fichresu << bibi.get_x_axe() / km ;
      fichresu << "    abscidia of the rotation axis [km] " << endl ;

      fichresu << bibi.get_omega() * f_unit ;
      fichresu << "    Orbital frequency Omega [rad/s] " << endl ;

      fichrota << "  " << bibi.get_omega() * f_unit ;
      fichrota << "  " << bibi.get_x_axe() / km ;

      //---------------------------------------------------------
      //     Computation of B^i/N (bsn) and pot_centri in NS
      //---------------------------------------------------------

      (bibi.set_ns()).kinematics( bibi.get_omega(), bibi.get_x_axe() ) ;

      //------------------------------------------------------------------
      //     Computation of gam_euler, u_euler, ener_euler, s_euler,
      //     wit_w and loggam in NS
      //------------------------------------------------------------------

      (bibi.set_ns()).fait_d_psi() ;
      (bibi.set_ns()).hydro_euler() ;

      // Check of the Bin_ns_bh::orbit_omega computation
      //------------------------------------------------
      // !!! These out-put has not been prepared.

      //---------------------------------------------------------------
      //     Computation of the stellar equilibrium configurations
      //---------------------------------------------------------------

      /*
      // !!! This rotine has not worked yet.
      // Computation of the resizing factor
      double ray_eq_auto = (bibi.get_ns()).ray_eq() ;
      double ray_eq_comp = (bibi.get_bh()).ray_eq() ;

      int num_resize ;

      if (mg_ns.get_nzone() > 3) {
	num_resize = mg_ns.get_nzone() - 3 ;
      }
      else {
	num_resize = (bibi.get_ns()).get_nzet() ;
      }

      double lambda_resize = 0.95 *
	(bibi.separation() - ray_eq_comp)/ray_eq_auto ;
      fact_resize[0].set(0) =
	(lambda_resize < 2.*num_resize) ? lambda_resize : 2.*num_resize ;
      */

      // Relaxation on n_comp (only if it has not been done by update_metric)
      // --------------------------------------------------------------------
      if ( (ind_rel_met == 0) || ( (mer % fmer_upd_met) != 0 ) ) {
	    bibi.set_ns().set_n_comp() = relax * (bibi.get_ns()).get_n_comp()
	      + relax_jm1 * (bibi_jm1.get_ns()).get_n_comp() ;
      }

      // Relaxation on pot_centri
      // ------------------------
      bibi.set_ns().set_pot_centri() = relax * (bibi.get_ns()).get_pot_centri()
	+ relax_jm1 * (bibi_jm1.get_ns()).get_pot_centri() ;

      // Call to Et_bin_nsbh::equilibrium_nsbh
      // -------------------------------------

      (bibi.set_ns()).equilibrium_nsbh(ent_c, mermax_eqb,
				       mermax_poisson, relax_poisson,
				       mermax_potvit, relax_potvit,
				       thres_adapt, fact_resize[0],
				       differ[0]) ;

    //------------------------------------------------------------------
    //	  Relaxations
    //------------------------------------------------------------------

    bibi.set_ns().relaxation( bibi_jm1.get_ns(), relax, relax_met, mer,
			      fmer_upd_met ) ;

    bibi.set_ns().hydro_euler() ;

    //------------------------------------------------------------------
    //	  Change in the central enthalpy to get a fixed baryon mass
    //------------------------------------------------------------------

    if (mer >= mer_masse) {

      double xx = (bibi.get_ns()).mass_b() / mbar_voulue - 1. ;

      cout << "Discrepancy M_b / wanted M_b : " << xx << endl ;

      double xprog = ( mer > 2*mer_masse) ? 1. :
	double(mer-mer_masse)/double(mer_masse) ;

      xx *= xprog ;
      double ax = .5 * ( 2. + xx ) / (1. + xx ) ;

      double fact_ent = pow(ax, aexp_masse) ;

      cout << "  xprog, xx, ax, fact : " << xprog << "  " << xx
	   << "  " << ax << "  " << fact_ent << endl ;

      ent_c *= fact_ent ;

    }

    // Updates for the next step
    // -------------------------

    bibi_jm1 = bibi ;

    cout << bibi << endl ;

    // Graphical output
    // ----------------
    // !!! This routine has not worked yet.

    //-----------------------------------------------------------------------
    //     The whole configuration is saved in a file
    //-----------------------------------------------------------------------

    if ( (mer % fmer_save) == 0 ) {

      FILE* fresud = fopen("resu.d", "w") ;

      fwrite(&mer, sizeof(int), 1, fresud) ;	// mer

      (bibi.get_ns()).get_mp().get_mg()->sauve(fresud) ;
      (bibi.get_ns()).get_mp().sauve(fresud) ;
      (bibi.get_ns()).get_eos().sauve(fresud) ;

      (bibi.get_bh()).get_mp().get_mg()->sauve(fresud) ;
      (bibi.get_bh()).get_mp().sauve(fresud) ;

      bibi.sauve(fresud) ;

      fclose(fresud) ;

    }

    //--------------------------------------------
    //  Writing of global quantities in log files
    //--------------------------------------------
    // !!! This routine has not worked yet.

    }	// End of the main loop (mer)

//=========================================================================
//		End of iteration
//=========================================================================

    fichresu.close() ;
    fichrota.close() ;
    fichconv.close() ;
    fichet.close() ;

    //-------------------------------------------------------------
    // General features of the final configuration saved in a file
    //-------------------------------------------------------------

    ofstream fichfinal("calcul.d") ;
    fichfinal.precision(6) ;

    time_t rawtime = time(0x0) ;
    fichfinal << "Date: " << asctime( localtime( &rawtime ) ) << endl ;
    char* hostname = getenv("HOST") ;
    if (hostname != 0x0) {
      fichfinal << "Computer: " << hostname << endl ;
    }
    fichfinal <<
      "==================================================================="
	      << endl << endl ;

    if ( (bibi.get_ns()).is_irrotational() ) {
      fichfinal << "                           Irrotational" << endl ;
    }
    else {
      fichfinal << "                           Co-rotating" << endl ;
    }

    fichfinal << (bibi.get_ns()).get_eos() << endl ;

    fichfinal << "Omega = " << bibi.get_omega() * f_unit << " rad/s"
	      << "                Orbital frequency f = "
	      << bibi.get_omega() / (2*M_PI) * f_unit << " Hz" << endl ;
    //    fichfinal << "Omega_kepler = " << omega_kep * f_unit << " rad/s" << endl ;
    fichfinal << "Coordinate separation : " << bibi.separation()  / km
	      << " km" << endl ;
    // ADM mass has not been defined.
    /*
    fichfinal << "1/2 ADM mass :        " << 0.5 * bibi.mass_adm() / msol
	      << " Mo" << endl ;
    */
    // Angular momentum has not been defined.
    /*
    fichfinal << "Total angular momentum : "
	      << bibi.angu_mom()(2)/ ( qpig / (4* M_PI) * msol*msol)
	      << " G M_sol^2 / c" << endl ;
    */

    fichfinal << endl << "Number of steps : " << mer << endl ;

    fichfinal << endl <<
      "==================================================================="
	      << endl ;
    fichfinal << "       Neutron Star" << endl ;
    fichfinal <<
      "==================================================================="
	      << endl ;
    fichfinal << "Grid : " << endl ;
    fichfinal << "------ " << endl ;
    fichfinal << *((bibi.get_ns()).get_mp().get_mg()) << endl ;
    fichfinal << endl << "Physical characteristics : " << endl ;
    fichfinal	  << "-------------------------" << endl ;
    fichfinal << bibi.get_ns() << endl ;

    // !!! It is necessary to save the information of EOS here.

    fichfinal << endl <<
      "==================================================================="
	      << endl ;
    fichfinal << "       Black Hole" << endl ;
    fichfinal <<
      "==================================================================="
	      << endl ;
    fichfinal << "Grid : " << endl ;
    fichfinal << "------ " << endl ;
    fichfinal << *((bibi.get_bh()).get_mp().get_mg()) << endl ;
    // No output for black hole
    /*
    fichfinal << endl << "Physical characteristics : " << endl ;
    fichfinal	  << "-------------------------" << endl ;
    fichfinal << bibi.get_bh() << endl ;
    */
    fichfinal << endl ;

    fichfinal << endl <<
      "==================================================================="
	      << endl ;
    fichfinal << "Diff_ent :  neutron star : " << differ[0](0) << endl ;
    fichfinal << "dH/dx at r = 0 :  star 1 : " << dentdx << endl ;

    fichfinal << endl <<
      "================================================================"
	      << endl ;
    fichfinal << "	    PARAMETERS USED FOR THE COMPUTATION : " << endl ;
    fichfinal <<
      "================================================================"
	      << endl ;
    fichfinal.close() ;
    system("cat parcoal.d >> calcul.d") ;

    // Identification du code et de ses sous-routines (no. de version RCS) :

    fichfinal.open("calcul.d", ios::app ) ;
    fichfinal << endl <<
      "================================================================"
	      << endl ;
    fichfinal << "	    IDENTIFICATION OF THE CODE : " << endl ;
    fichfinal <<
      "================================================================"
	      << endl ;
    fichfinal.close() ;
    system("ident coal >> calcul.d") ;


    // Cleaning
    // --------

    delete peos ;

    return EXIT_SUCCESS ;

}
