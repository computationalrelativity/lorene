/*
 * Computes the equilibrium configuration of a binary system. 
 * 
 */

/*
 *   Copyright (c) 2000-2003 Eric Gourgoulhon
 *   Copyright (c) 2001-2002 Keisuke Taniguchi
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

char coal_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.9  2003/09/18 15:06:59  e_gourgoulhon
 * Added printing of date and host name in the header of the file calcul.d.
 *
 * Revision 1.8  2003/09/17 15:40:33  e_gourgoulhon
 * Test of good opening of the initial data file.
 *
 * Revision 1.7  2003/09/16 13:36:50  e_gourgoulhon
 * -- initial analytical shift decreased by a factor reduce_shift which
 *    can be specified in the parcoal.d file (default value = 0.6)
 * -- replaced the fpar.getline(blabla, 120) by fpar.ignore(1000,'\n')
 * -- new output file "resformat.d" at the end of the computation for
 *    formatted output of global quantities (to be read by external code)
 *
 * Revision 1.6  2003/09/08 12:06:20  e_gourgoulhon
 * Added the printing of the virial errors in the output file "calcul.d".
 *
 * Revision 1.5  2003/02/20 15:08:42  e_gourgoulhon
 * Suppressed the qualifier ios::nocreate in call to fstream::open
 * (not supported by gcc 3.2).
 *
 * Revision 1.4  2003/01/09 11:07:48  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.3  2002/12/11 13:53:01  k_taniguchi
 * Add "thres_adapt" to star 2,
 *   and modify the part of "resize".
 *
 * Revision 1.2  2001/12/06 16:20:14  e_gourgoulhon
 * Return type of main changed from 'void' to 'int'
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 2.18  2001/08/07  09:51:33  keisuke
 * Change of the argument in Etoile_bin::equilibrium.
 * Addition of the procedure to calculate a factor for "resize".
 *
 * Revision 2.17  2001/03/20  16:03:22  keisuke
 * Translation of the origin of the absolute frame.
 *
 * Revision 2.16  2000/07/07  14:18:01  eric
 * Sortie des quantites en unites polytropiques dans calcul.d
 *  (appel de la nouvelle fonction  Binaire::display_poly).
 *
 * Revision 2.15  2000/07/06  10:06:35  eric
 * Ajout des parametres et de l'identification dans calcul.d
 * Creation du fichier identif.d pour l'identification du code.
 *
 * Revision 2.14  2000/05/25  15:19:55  eric
 * Adaptation du mapping gele (cusp).
 *
 * Revision 2.13  2000/03/29  09:13:35  eric
 * Suppression des appels aux routines de dessin (car probleme sur mesioc).
 *
 * Revision 2.12  2000/03/22  16:42:22  eric
 * Sortie des erreurs sur les equations de Poisson.
 *
 * Revision 2.11  2000/03/22  14:00:44  eric
 * Les sorties pour Explorer ne sont effectuees que si graph = 1.
 *
 * Revision 2.10  2000/03/22  13:56:41  eric
 * Changement de prototypage d'Etoile_bin:equilibrium: on passe desormais
 * le Tbl differ pour avoir en sortie des estimateurs d'erreur.
 * Reorganisation des fichiers log : introduction de resconv1 et resconv2
 *
 * Revision 2.9  2000/03/22  11:08:00  eric
 * Reintroduction des messages Log RCS des versions 2.4-2.6 qui avaient
 *  disparus.
 * /
 *
 * Revision 2.8  2000/03/20  14:25:33  eric
 * resviriel.d rebaptise resconver.d
 * Theoreme du viriel calcule a la fin et sorti dans calcul.d
 * Suppression de l'ecriture du mapping dans calcul.d
 *
 * Revision 2.7  2000/03/17  16:52:50  eric
 * Rajout des headers Id et Log qui avaient disparus.
 *
 *
 * Revision 2.6 
 * Mutliplication du shift analytique par un facteur 0.6
 * Omega_initial donne par Binaire::analytical_omega().
 *
 * Revision 2.5 
 * Introduction d'un shift analytique dans le cas ou le shift initial est nul
 *
 * Revision 2.4 
 * Sortie des log(diff_ent) et log(diff_masse).
 *
 * Revision 2.3  2000/02/22  16:58:15  eric
 * Ajout de la convergence vers une masse baryonique donnee.
 *
 * Revision 2.2  2000/02/21  16:47:27  eric
 * Appel de fait_d_psi avant hydro_euler.
 * Suppression dessin surface.
 *
 * Revision 2.1  2000/02/17  19:54:44  eric
 * Premiere version operationnelle !
 *
 * Revision 2.0  2000/02/12  18:41:28  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers C
#include <stdlib.h>
#include <math.h>
#include <time.h>

// headers Lorene
#include "binaire.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"

// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ; 

//******************************************************************************

int main(){

    // Identification of all the subroutines called by the code : 
    
    system("ident coal > identif.d") ; 

    // For the display : 
    char display_bold[]="x[1m" ; display_bold[0] = 27 ;
    char display_normal[] = "x[0m" ; display_normal[0] = 27 ;

    #include "unites.h"	    
    // To avoid some compilation warnings
    if (display_bold == 0x0) {
	cout << qpig << f_unit << mevpfm3 << endl ; 
    }    
    
    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    char nomini[80] ;
    int mermax, mermax_eqb, prompt, graph, fmer_stop, fmer_save, mermax_poisson ;
    int mermax_potvit, mer_masse, fmer_upd_met, ind_rel_met ;  
    double seuil, relax_poisson, relax_potvit, relax, aexp_masse ; 
    double mbar_voulue[2], fact_separ, relax_met, relax_omeg ;
    double fact_omeg_min, fact_omeg_max, thres_adapt[2], reduce_shift ; 
    
    ifstream fpar("parcoal.d") ;
	if ( !fpar.good() ) {
		cout << "Problem with opening the file parcoal.d ! " << endl ;
		abort() ;
	}
    fpar.ignore(1000, '\n') ;
    fpar.ignore(1000, '\n') ;
    fpar.getline(nomini, 80) ; 
    fpar >> fact_separ ; fpar.ignore(1000, '\n');
    fpar >> mbar_voulue[0] ; fpar.ignore(1000, '\n') ;  
    fpar >> mbar_voulue[1] ; fpar.ignore(1000, '\n') ;  
    mbar_voulue[0] *= msol ;
    mbar_voulue[1] *= msol ;
    fpar.ignore(1000, '\n') ;
    fpar >> mermax ; fpar.ignore(1000, '\n') ;   
    fpar >> relax ; fpar.ignore(1000, '\n') ;  
    fpar >> mermax_eqb ; fpar.ignore(1000, '\n') ;  
    fpar >> prompt ; fpar.ignore(1000, '\n') ;  
    fpar >> graph ; fpar.ignore(1000, '\n') ;  
    fpar >> seuil ; fpar.ignore(1000, '\n') ;  
    fpar >> fmer_stop ; fpar.ignore(1000, '\n') ;  
    fpar >> fmer_save ; fpar.ignore(1000, '\n') ;  
    fpar >> mermax_poisson ; fpar.ignore(1000, '\n') ;  
    fpar >> relax_poisson ; fpar.ignore(1000, '\n') ;  
    fpar >> mermax_potvit ; fpar.ignore(1000, '\n') ;  
    fpar >> relax_potvit ; fpar.ignore(1000, '\n') ;  
    fpar >> mer_masse ; fpar.ignore(1000, '\n') ;  
    fpar >> aexp_masse ; fpar.ignore(1000, '\n') ;  
    fpar >> fmer_upd_met ; fpar.ignore(1000, '\n');
    fpar >> ind_rel_met ; fpar.ignore(1000, '\n');
    fpar >> relax_met ; fpar.ignore(1000, '\n');
    if (ind_rel_met == 0) relax_met = 1. ; 
    fpar >> relax_omeg ; fpar.ignore(1000, '\n');
    fpar >> fact_omeg_min ; fpar.ignore(1000, '\n');
    fpar >> fact_omeg_max ; fpar.ignore(1000, '\n');
    fpar >> thres_adapt[0] ; fpar.ignore(1000, '\n');
    fpar >> thres_adapt[1] ; fpar.ignore(1000, '\n');
	fpar >> reduce_shift ; 
	if ( ! fpar.good() ) {	  // to ensure compatibility with old 
		reduce_shift = 0.6 ;  // parcoal.d files which did not had
	}						  // the reduce_shift line
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
	cout << "Baryon mass required for star 1 [M_sol] : " 
	     << mbar_voulue[0] / msol << endl ; 
	cout << "Baryon mass required for star 2 [M_sol] : " 
	     << mbar_voulue[1] / msol << endl ; 
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
    cout << "Maximum number of steps in Etoile_bin::equilibrium : " 
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
    "Threshold on |dH/dr|_eq / |dH/dr|_pole for the adaptation of the mapping for star 1"
    << endl << thres_adapt[0] << endl ;
    cout << 
    "Threshold on |dH/dr|_eq / |dH/dr|_pole for the adaptation of the mapping for star 2"
    << endl << thres_adapt[1] << endl ;
	cout << "Factor by which the initial analytical shift is reduced : "
		<< reduce_shift << endl ; 

    arrete(prompt) ; 
    
    //------------------------------------------------------------------
    //	    Read of the initial conditions 
    //------------------------------------------------------------------

    FILE* fich = fopen(nomini, "r") ; 
    if (fich == 0x0) {
    	cout << "Problem in opening the file " << nomini << " ! " << endl ; 
	perror(" reason") ; 
	abort() ; 
    }

    int mer_ini ; 
    fread(&mer_ini, sizeof(int), 1, fich) ;	
    
    Mg3d mg1(fich) ;
    Map_et mp1(mg1, fich) ; 
    Eos* peos1 = Eos::eos_from_file(fich) ; 
    
    Mg3d mg2(fich) ;
    Map_et mp2(mg2, fich) ; 
    Eos* peos2 = Eos::eos_from_file(fich) ; 
    
    Binaire star(mp1, *peos1, mp2, *peos2, fich) ; 

    fclose(fich) ; 

    //------------------------------------------------------------------
    //	    Modification of the separation between the two stars
    //------------------------------------------------------------------

    for (int i=1 ; i<=2 ; i++) {

	double ori_x = (star(i).get_mp()).get_ori_x() ; 
	ori_x *= fact_separ ; 
	((star.set(i)).set_mp()).set_ori(ori_x, 0., 0.) ; 	
    }

    //------------------------------------------------------------------
    //	    Update of the initial conditions 
    //------------------------------------------------------------------

    // Initialisation of A^2, N, etc...
    // ---------------------------------
    for (int i=1; i<=2; i++) {
	(star.set(i)).update_metric(star(3-i)) ; 
    }

    // Initialisation of gradients of companion potentials
    // ---------------------------------------------------

    for (int i=1; i<=2; i++) {
	(star.set(i)).update_metric_der_comp(star(3-i)) ; 
    }

    // Initialisation of hydro quantities
    // ----------------------------------

    for (int i=1; i<=2; i++) {
	(star.set(i)).equation_of_state() ; 
	(star.set(i)).kinematics(star.get_omega(), star.get_x_axe()) ; 
	(star.set(i)).fait_d_psi() ; 
	(star.set(i)).hydro_euler() ; 
    }

    // If the shift vector has not been set previously, it is set to
    //  some analytical value
    // -------------------------------------------------------------
    
    if (star(1).get_w_shift()(1).get_etat() == ETATZERO) {

	assert( star(1).get_w_shift()(0).get_etat() == ETATZERO ) ;
	assert( star(1).get_w_shift()(2).get_etat() == ETATZERO ) ;

	assert( star(2).get_w_shift()(0).get_etat() == ETATZERO ) ;
	assert( star(2).get_w_shift()(1).get_etat() == ETATZERO ) ;
	assert( star(2).get_w_shift()(2).get_etat() == ETATZERO ) ;


	// New initial of value Omega (taking into account the fact
	//  that the separation has changed)

	star.analytical_omega() ; 

	// Sets some value to w_shift and khi_shift corresponding to the
	//  new values of Omega and separation

	star.analytical_shift() ;
	
	//## Shift decreased by a factor 2 
	for (int i=1; i<=2; i++) {
	     star.set(i).set_w_shift() = 0.6 * star(i).get_w_shift() ; 
	     star.set(i).set_khi_shift() = 0.6 * star(i).get_khi_shift() ; 
	}
	

	// Computes shift_auto from w_shift and khi_shift
	for (int i=1; i<=2; i++) {
	    (star.set(i)).fait_shift_auto() ; 
	}
	
	// A second call to update_metric must be performed to update
	//  shift_comp, tkij_auto and akcar_auto. 
	for (int i=1; i<=2; i++) {
	    (star.set(i)).update_metric(star(3-i)) ; 
	}
	
	// Second update of gradients of companion potentials

	for (int i=1; i<=2; i++) {
	    (star.set(i)).update_metric_der_comp(star(3-i)) ; 
	}

	// Second update of hydro quantities

	for (int i=1; i<=2; i++) {
	    (star.set(i)).equation_of_state() ; 
	    (star.set(i)).kinematics(star.get_omega(), star.get_x_axe()) ; 
	    (star.set(i)).fait_d_psi() ; 
	    (star.set(i)).hydro_euler() ; 
	}
	
	
    }  // End of the setting of an analytical shift






//##
	FILE* fresu = fopen("resu.d", "w") ; 
    
	int mer1 = 0 ;
	fwrite(&mer1, sizeof(int), 1, fresu) ;	// mer

	star(1).get_mp().get_mg()->sauve(fresu) ; 
	star(1).get_mp().sauve(fresu) ; 
	star(1).get_eos().sauve(fresu) ; 

	star(2).get_mp().get_mg()->sauve(fresu) ; 
	star(2).get_mp().sauve(fresu) ; 
	star(2).get_eos().sauve(fresu) ; 

	star.sauve(fresu) ;     
	
	fclose(fresu) ;     
//##

    cout << endl 
	 << "==========================================================" << endl
	 << "                    Initial conditions                    " << endl
	 << "=========================================================="
	 << endl ; 
    cout << star << endl ; 

    if ( star(1).is_relativistic() ) {
	cout << "========================" << endl ;
	cout << "Relativistic computation" << endl ;
	cout << "========================" << endl ;
    }
    else {
	cout << "=====================" << endl ;
	cout << "Newtonian computation" << endl ;
	cout << "=====================" << endl ;
    }

    arrete(prompt) ; 
    

    //----------------------------------------------------------------
    //			Auxiliary quantities
    //----------------------------------------------------------------

    double ent_c[2] ;	    // Central enthalpy in each star
    double dentdx[2] ;	    // Central d/dx(enthalpy) in each star

    // Resizing factor of the first shell
    // Computation is shown just befor "equilibrium"
    Tbl fact_resize[] = {Tbl(1), Tbl(1)} ;
    fact_resize[0].set_etat_qcq() ; 
    fact_resize[1].set_etat_qcq() ; 

    // Error indicators in each star
    Tbl differ[] = {Tbl(7), Tbl(7)} ;
    differ[0].set_etat_qcq() ; 
    differ[1].set_etat_qcq() ; 

    for (int i=1 ; i<=2 ; i++) {
	ent_c[i-1] = star(i).get_ent()()(0, 0, 0, 0) ;  
	differ[i-1].set(0) = 1 ;	    // diff_ent = 1
	differ[i-1].set(1) = 1 ;	    // err_psi = 1
	differ[i-1].set(2) = 1 ;	    // err_logn = 1
	differ[i-1].set(3) = 1 ;	    // err_beta = 1
	differ[i-1].set(4) = 1 ;	    // err_shift_x = 1
	differ[i-1].set(5) = 1 ;	    // err_shift_y = 1
	differ[i-1].set(6) = 1 ;	    // err_shift_z = 1
    }
    
    double relax_jm1 = 1. - relax ; 
    double relax_omeg_jm1 = 1. - relax_omeg ; 
    
    //----------------------------------------------------------------
    //	 Binary system at the previous step (for the relaxation)
    //----------------------------------------------------------------

    Binaire star_jm1 = star ;
     
    double omega_jm1 = star_jm1.get_omega() ; 

    
    // logn_comp and pot_centri are initialized to 0 on star_jm1 : 
    // ---------------------------------------------------------
    for (int i=1 ; i<=2 ; i++) {
	star_jm1.set(i).set_logn_comp() = 0 ; 
	star_jm1.set(i).set_pot_centri() = 0 ; 
    }
    
    //----------------------------------------------------------------
    //	 Openning of log files
    //----------------------------------------------------------------

    ofstream fichresu("resglob.d") ;
    fichresu.precision(16) ; 

    ofstream fichrota("resrota.d") ; 
    fichrota.precision(16) ; 

    ofstream fichvir("resdiffm.d") ;
    fichvir.precision(16) ; 

    ofstream fichconv[2] ;
    fichconv[0].open("resconv1.d") ; 
    fichconv[0].precision(16) ; 
    fichconv[1].open("resconv2.d") ; 
    fichconv[1].precision(16) ; 

    ofstream fichet[2] ;
    fichet[0].open("resstar1.d") ; 
    fichet[0].precision(16) ; 
    fichet[1].open("resstar2.d") ;
    fichet[1].precision(16) ; 
 
    fichrota << 
      "#	Omega [rad/s]	    x_axe [km]		x_g (et 0) [km]		x_g (et 1) [km]	   M_grav [M_sol]    J [ G M_sol^2 / c]"
      << endl ;

    fichvir << 
      "#    diff_mass"
      << endl ; 

    for (int i=1; i<=2; i++) {
	fichconv[i-1] << 
      "#     diff_ent           err_psi         err_logn          err_beta        err_shift_x         err_shift_y         err_shift_z"
      << endl ; 

	fichet[i-1] << 
      "#     ori_x [km]      ent_c      M_bar [M_sol]	  R(theta=0) [km]   R(pi/2, 0) [km]   R(pi/2, pi/2) [km]   R(pi/2, pi) [km]"
      << endl ; 
    }      

    double omega_kep, diff_mass ; 
    int mer ; 
    
//============================================================================
//		Start of iteration 
//============================================================================

    for (mer=0; (differ[0](0) > seuil) && (mer < mermax); mer++) {

    cout << 
    "=========================================================================="
    << endl ;
    cout << "step = " << mer << "        diff_ent (1) (2) (1)<->(2) : "
	 << differ[0](0) << "  " << differ[1](0) << "  "
	 << differ[0](0) - differ[1](0) << endl ; 
    cout << 
    "=========================================================================="
    << endl ;

    fichresu << mer; 
    fichresu << "    step" << endl ; 

    fichrota << mer ; 
    fichvir << mer ; 
    fichconv[0] << mer ; 
    fichconv[1] << mer ; 
    fichet[0] << mer ; 
    fichet[1] << mer ; 


    //------------------------------------------------------------------
    //	    Computation of the metric coefficients
    //------------------------------------------------------------------

    if ( (mer % fmer_upd_met) == 0 ) {
    
	for (int i=1; i<=2; i++) {
	    (star.set(i)).update_metric(star(3-i), star_jm1(i), relax_met) ; 
	}

	for (int i=1; i<=2; i++) {
	    (star.set(i)).update_metric_der_comp(star(3-i)) ; 
	}

    }


    //------------------------------------------------------------------
    //	    Computation of the orbital angular velocity Omega
    //------------------------------------------------------------------

    double xgg[2] ; 

    star.orbit(fact_omeg_min, fact_omeg_max, xgg[0], xgg[1]) ; 

    // Translation of the stars in order to set the origin
    //  of the absolute frame on the rotation axis
    //-----------------------------------------------------

    double x_rot = star.get_x_axe() ;

    for (int i=1 ; i<=2 ; i++) {
      double ori_x_old = (star(i).get_mp()).get_ori_x() ;
      double ori_x_new = ori_x_old - x_rot ;
      ((star.set(i)).set_mp()).set_ori(ori_x_new, 0., 0.) ;
    }

    star.set_x_axe() = 0. ;

    // Relaxation on the orbital velocity
    // ----------------------------------

    double omega_j = star.get_omega() ; 
    omega_j = relax_omeg * omega_j + relax_omeg_jm1 * omega_jm1 ; 
    omega_jm1 = omega_j ; 

    star.set_omega() = omega_j ; 

    cout << display_bold << "New orbital velocity Omega : " 
	 << star.get_omega() * f_unit << " rad/s" << display_normal << endl ; 

    // Keplerian velocity (for comparison only)
    // ----------------------------------------
    omega_kep = sqrt( g_si/g_unit * (star(1).mass_g() + star(2).mass_g()) 
			    / pow( star.separation(), 3.) ) ; 

    cout << "``Keplerian'' velocity (for comparison only) : " 
	 << omega_kep * f_unit << " rad/s" << endl ; 
    cout << "New X coordinate of the rotation axis : " 
	 << star.get_x_axe() / km << " km" << endl ; 

    arrete(prompt) ; 

    fichresu << star.get_x_axe() / km ; 
    fichresu << "    abscidia of the rotation axis [km] " << endl ; 

    fichresu << star.get_omega() * f_unit ; 
    fichresu << "    Orbital frequency Omega [rad/s] " << endl ; 

    fichresu << xgg[0] / km ; 
    fichresu << "    Abscidia ``center of mass'' star 1 [km] " << endl ; 

    fichresu << xgg[1] / km ; 
    fichresu << "    Abscidia ``center of mass'' star 2 [km] " << endl ; 

    fichrota << "  " << star.get_omega() * f_unit ;
    fichrota << "  " << star.get_x_axe() / km ;
    fichrota << "  " << xgg[0] / km ;
    fichrota << "  " << xgg[1] / km ;

    //------------------------------------------------------------------
    //	    Computation of B^i/N (bsn) and pot_centri in each star
    //------------------------------------------------------------------

    for (int i=1; i<=2; i++) {
	
	(star.set(i)).kinematics( star.get_omega(), star.get_x_axe() ) ; 
	
    }

    //------------------------------------------------------------------
    //	    Computation of gam_euler, u_euler, ener_euler, s_euler, 
    //	    wit_w and loggam  in each star
    //------------------------------------------------------------------

    for (int i=1; i<=2; i++) {
	
	(star.set(i)).fait_d_psi() ; 
	(star.set(i)).hydro_euler() ; 
	
	// Check of the Binaire::orbit computation 
	//----------------------------------------
    
	Cmp tmp = star(i).get_logn_auto()() + star(i).get_logn_comp()()  
			+ star(i).get_loggam()() ;
	double grad1 = tmp.dsdx()(0, 0, 0, 0) ;

	double grad2 = star(i).get_pot_centri()().dsdx()(0, 0, 0, 0) ; 

	dentdx[i-1] = star(i).get_ent()().dsdx()(0, 0, 0, 0) ; 

	double grad3 = star(i).get_loggam()().dsdx()(0, 0, 0, 0) ; 

	cout << "Star " << i << " : " << endl ; 
	cout << "  central dH/dx  : " <<  dentdx[i-1] << endl ; 
	cout << "  central d(log(Gam))/dx  : " <<  grad3 << endl ; 
	cout << "  central d/dx(nu + log(Gam)) : " << grad1 << endl ; 
	cout << "  central d/dx(pot_centri) : " << grad2 << endl ; 
	cout << "  central d/dx(nu + log(Gam) + pot_centri) : " 
	     << grad1 + grad2 << endl ; 

	
    }

    
    //------------------------------------------------------------------
    //	  Computation of the stellar equilibrium configurations
    //------------------------------------------------------------------

    for (int i=1; i<=2; i++) {

	// Computation of the resizing factor
	double ray_eq_auto = star(i).ray_eq() ;
	double ray_eq_comp = star(3-i).ray_eq() ;

	int num_resize ;

      	if (mg1.get_nzone() > 3) {
	  assert( mg2.get_nzone() == mg1.get_nzone() ) ;
	  num_resize = mg1.get_nzone() - 3 ;
	}
	else {
	  num_resize = star(i).get_nzet() ;
	}

	double lambda_resize = 0.95 *
	  (star.separation() - ray_eq_comp)/ray_eq_auto ;
	fact_resize[i-1].set(0) =
	  (lambda_resize < 2.*num_resize) ? lambda_resize : 2.*num_resize ;

    }

    for (int i=1; i<=2; i++) {

	// Relaxation on logn_comp (only if it has not been done by
	//			    update_metric)
	// --------------------------------------------------------
	if ( (ind_rel_met == 0) || ( (mer % fmer_upd_met) != 0 ) ) {
	    star.set(i).set_logn_comp() = relax * star(i).get_logn_comp()
				+ relax_jm1 * star_jm1(i).get_logn_comp() ; 
	}

	// Relaxation on pot_centri
	// ------------------------
	star.set(i).set_pot_centri() = relax * star(i).get_pot_centri()
				+ relax_jm1 * star_jm1(i).get_pot_centri() ; 


	// Call to Etoile_bin::equilibrium
	// --------------------------------
	
	(star.set(i)).equilibrium(ent_c[i-1], mermax_eqb, mermax_poisson, 
				  relax_poisson, mermax_potvit, relax_potvit, 
				  thres_adapt[i-1], fact_resize[i-1], differ[i-1]) ;

    }


    //------------------------------------------------------------------
    //	  Relaxations
    //------------------------------------------------------------------

    for (int i=1; i<=2; i++) {

	star.set(i).relaxation( star_jm1(i), relax, relax_met, mer, 
				fmer_upd_met ) ; 

	star.set(i).hydro_euler() ; 
    }    
    
    //------------------------------------------------------------------
    //	  Change in the central enthalpy to get a fixed baryon mass
    //------------------------------------------------------------------

    if (mer >= mer_masse) {

	for (int i=1; i<=2; i++) {

	double xx = star(i).mass_b() / mbar_voulue[i-1] - 1. ;

	cout << "Discrepancy M_b / wanted M_b : " << xx << endl ; 
	
	double xprog = ( mer > 2*mer_masse) ? 1. : 
				 double(mer-mer_masse)/double(mer_masse) ; 
	xx *= xprog ; 
	double ax = .5 * ( 2. + xx ) / (1. + xx ) ; 

	double fact_ent = pow(ax, aexp_masse) ; 

	cout << "  xprog, xx, ax, fact : " << xprog << "  " <<
	xx << "  " << ax << "  " << fact_ent << endl ; 
	
	ent_c[i-1] *= fact_ent ; 

	}
    }

    // Updates for the next step
    // -------------------------
    
    star_jm1 = star ; 
    

    cout << star << endl ; 


    // Graphical output
    // ----------------
    
    if ( (graph==1) && (mer % 5 == 0) ) {
	double xdes_min = - star(1).ray_eq_pi() + star(1).get_mp().get_ori_x() ; 
	xdes_min *= 1.5 ; 
	double xdes_max = star(2).ray_eq_pi() + star(2).get_mp().get_ori_x() ; 
	xdes_max *= 1.5 ; 
	double ydes_min = - 2.5 * star(1).ray_eq_pis2() ; 
	double ydes_max =  2.5 * star(2).ray_eq_pis2() ; 

	Cmp surf1 = star(1).get_ent()() ; 
	Cmp surf1_ext(mp1) ; 
	surf1_ext = - 0.2 * surf1(0, 0, 0, 0) ; 
	surf1_ext.annule(0, star(1).get_nzet()-1) ; 
	surf1.annule(star(1).get_nzet(), mg1.get_nzone()-1) ; 
	surf1 = surf1 + surf1_ext ;
    
	surf1 = raccord_c1(surf1, star(1).get_nzet()) ; 

	Cmp surf2 = star(2).get_ent()() ; 
	Cmp surf2_ext(mp2) ; 
	surf2_ext = - 0.2 * surf2(0, 0, 0, 0) ; 
	surf2_ext.annule(0, star(2).get_nzet()-1) ; 
	surf2.annule(star(2).get_nzet(), mg2.get_nzone()-1) ; 
	surf2 = surf2 + surf2_ext ;

	surf2 = raccord_c1(surf2, star(2).get_nzet()) ; 

	des_coupe_bin_z(star(1).get_nbar()(), star(2).get_nbar()(), 0., 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    "Baryon density (z=0)", &surf1, &surf2 ) ; 

	des_coupe_bin_y(star(1).get_nbar()(), star(2).get_nbar()(), 0., 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    "Baryon density (y=0)", &surf1, &surf2 ) ; 
    }

    

    //-----------------------------------------------------------------------
    //		The whole configuration is saved in a file
    //-----------------------------------------------------------------------
 
    if ( (mer % fmer_save) == 0 ) {

	FILE* fresu = fopen("resu.d", "w") ; 
    
	fwrite(&mer, sizeof(int), 1, fresu) ;	// mer

	star(1).get_mp().get_mg()->sauve(fresu) ; 
	star(1).get_mp().sauve(fresu) ; 
	star(1).get_eos().sauve(fresu) ; 

	star(2).get_mp().get_mg()->sauve(fresu) ; 
	star(2).get_mp().sauve(fresu) ; 
	star(2).get_eos().sauve(fresu) ; 

	star.sauve(fresu) ;     
	
	fclose(fresu) ;     
    }

    //--------------------------------------------
    //  Writing of global quantities in log files
    //--------------------------------------------
    for (int i=1 ; i<=2 ; i++) {
	fichresu << star(i).mass_b() / msol ;
	fichresu << "   Baryon mass of star " << i << "  [M_sol] " << endl ;
	fichresu << differ[i-1](0) ;
	fichresu << "   relative variation enth. star " << i << endl ;

	fichresu << star(i).ray_pole() / km ;
	fichresu << "   R(theta=0) [km] " << endl ;
	fichresu << star(i).ray_eq() / km ;
	fichresu << "   R(theta=pi/2, phi=0) [km] " << endl ;
	fichresu << star(i).ray_eq_pis2() / km ;
	fichresu << "   R(theta=pi/2, phi=pi/2) [km] " << endl ;
	fichresu << star(i).ray_eq_pi() / km ; 
	fichresu << "   R(theta=pi/2, phi=pi) [km] " << endl ;

	fichconv[i-1] << "  " << log10( fabs(differ[i-1](0)) + 1.e-16 ) ;
	fichconv[i-1] << "  " << log10( fabs(differ[i-1](1)) + 1.e-16 ) ;
	fichconv[i-1] << "  " << log10( fabs(differ[i-1](2)) + 1.e-16 ) ;
	fichconv[i-1] << "  " << log10( fabs(differ[i-1](3)) + 1.e-16 ) ;
	fichconv[i-1] << "  " << log10( fabs(differ[i-1](4)) + 1.e-16 ) ;
	fichconv[i-1] << "  " << log10( fabs(differ[i-1](5)) + 1.e-16 ) ;
	fichconv[i-1] << "  " << log10( fabs(differ[i-1](6)) + 1.e-16 ) ;
    
	fichet[i-1] << "  " << star(i).get_mp().get_ori_x() / km ;
	fichet[i-1] << "  " << ent_c[i-1] ;
	fichet[i-1] << "  " << star(i).mass_b() / msol ;
	fichet[i-1] << "  " << star(i).ray_pole() / km ;
	fichet[i-1] << "  " << star(i).ray_eq() / km ;
	fichet[i-1] << "  " << star(i).ray_eq_pis2() / km ;
	fichet[i-1] << "  " << star(i).ray_eq_pi() / km ;

    } // End of loop on the stars

    diff_mass = ( star(2).mass_b() - star(1).mass_b() ) 
			/ star(1).mass_b() ;
    cout << "Relative difference between the baryon masses: " 
	 << diff_mass << endl ;
    fichresu << diff_mass ; 
    fichresu << "   Relative difference between the baryon masses" << endl ;

    fichvir << "  " << log10( fabs(diff_mass) + 1.e-16 )  ; 


    fichresu << "  " << endl ; 
    fichresu.flush() ; 
    fichrota << "  " << endl ; 
    fichrota.flush() ; 
    fichvir << "  " << endl ; 
    fichvir.flush() ; 
    fichconv[0] << "  " << endl ; 
    fichconv[0].flush() ; 
    fichconv[1] << "  " << endl ; 
    fichconv[1].flush() ; 
    fichet[0] << "  " << endl ; 
    fichet[0].flush() ; 
    fichet[1] << "  " << endl ; 
    fichet[1].flush() ; 

    }	// End of the main loop (mer)

//============================================================================
//		End of iteration 
//============================================================================

    fichresu.close() ; 
    fichrota.close() ; 
    fichvir.close() ; 
    fichconv[0].close() ; 
    fichconv[1].close() ; 
    fichet[0].close() ; 
    fichet[1].close() ; 
        
    
    //-----------------------------------------------
    //  General features of the final configuration
    //  saved in a file
    //-----------------------------------------------

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
	
    if ( star(1).is_relativistic() ) {
	fichfinal << "Relativistic computation"  ;
    }
    else {
	fichfinal << "Newtonian computation"  ;
    }
    
    if ( star(1).is_irrotational() ) {
	fichfinal << "                           Irrotational" << endl ; 
    }
    else {
	fichfinal << "                           Co-rotating" << endl ; 
    }   

    fichfinal << star(1).get_eos() << endl ;
    
    fichfinal << "Omega = " << star.get_omega() * f_unit << " rad/s" 
	<< "                Orbital frequency f = " 
	<< star.get_omega() / (2*M_PI) * f_unit << " Hz" << endl ; 
    fichfinal << "Omega_kepler = " << omega_kep * f_unit << " rad/s" << endl ; 
    fichfinal << "Coordinate separation : " << star.separation()  / km 
	      << " km" << endl ; 
    fichfinal << "1/2 ADM mass :        " << 0.5 * star.mass_adm() / msol 
	      << " Mo" << endl ;
    fichfinal << "Total angular momentum : "  
	      << star.angu_mom()(2)/ ( qpig / (4* M_PI) * msol*msol)
	      << " G M_sol^2 / c" << endl ;

    fichfinal << endl << "Number of steps : " << mer << endl ;

    for (int i=1 ; i<=2; i++) {
	fichfinal << endl <<
	"===================================================================" 
	<< endl ; 
	fichfinal << "       Star no. " << i << endl ; 
	fichfinal <<
	"===================================================================" 
	<< endl ; 
	fichfinal << "Grid : " << endl ; 
	fichfinal << "------ " << endl ; 
	fichfinal << *(star(i).get_mp().get_mg()) << endl ; 
	fichfinal << endl << "Physical characteristics : " << endl ; 
	fichfinal	  << "-------------------------" << endl ; 
	fichfinal << star(i) << endl ; 
    }
    fichfinal << endl ;
    
    star.display_poly(fichfinal) ; // Reduced quantities for polytropic EOS

    fichfinal << endl <<
    "===================================================================" 
    << endl ; 
    fichfinal << "Diff_ent :  star 1 : " << differ[0](0) << "   star 2 : "
	<< differ[1](0) << endl ; 
    fichfinal << 
    "Relative difference between the baryon masses of the two stars : " 
	      << diff_mass << endl ; 
    fichfinal << "dH/dx at r = 0 :  star 1 : " << dentdx[0] 
		<< "   star 2 : " << dentdx[1] << endl ; 


    if ( star(1).is_relativistic() ) {
	fichfinal << "Relative error on the virial theorem : " << endl ; 
	fichfinal << "   VE(M)= " << star.virial() 
	          << "   VE(GB)= "<< star.virial_gb()  
		  << "   VE(FUS)= " << star.virial_fus() << endl ;   
    }
    else {
	fichfinal << "Relative error on the virial theorem : " 
		  <<  star.virial() << endl ; 
    }
    
    fichfinal << endl <<
    "================================================================" << endl ; 
    fichfinal << "	    PARAMETERS USED FOR THE COMPUTATION : " << endl ; 
    fichfinal << 
    "================================================================" << endl ; 
    fichfinal.close() ; 
    system("cat parcoal.d >> calcul.d") ; 

    // Identification du code et de ses sous-routines (no. de version RCS) :     	
    fichfinal.open("calcul.d", ios::app ) ; 
    fichfinal << endl <<
    "================================================================" << endl ; 
    fichfinal << "	    IDENTIFICATION OF THE CODE : " << endl ; 
    fichfinal << 
    "================================================================" << endl ; 
    fichfinal.close() ; 
    system("ident coal >> calcul.d") ; 
        
    // Preparation for CPU infos printing :     	
    fichfinal.open("calcul.d", ios::app ) ; 
    fichfinal << endl <<
    "================================================================" << endl ; 
    fichfinal << "	    CPU TIME and MEMORY infos : " << endl ; 
    fichfinal << 
    "================================================================" << endl ; 
	fichfinal << endl ; 
    fichfinal.close() ; 
        

    //-----------------------------------------------
    //  General features of the final configuration
    //  saved in a file with scientific notation and
	//  14 digits for further reading by a code
    //-----------------------------------------------

	ofstream seqfich("resformat.d") ; 
	if ( !seqfich.good() ) {
		cout << "coal : problem with opening the file resformat.d !" << endl ;
		abort() ;
	}
	star.write_global(seqfich) ; 
	seqfich.close() ; 
	
    
    if (graph == 1) {

	des_explorer(star, "latbin") ; 

//##	double xdes_min = - star(1).ray_eq_pi() + star(1).get_mp().get_ori_x() ; 
//	xdes_min *= 1.5 ; 
//	double xdes_max = star(2).ray_eq_pi() + star(2).get_mp().get_ori_x() ; 
//	xdes_max *= 1.5 ; 
//	double ydes_min = - 2.5 * star(1).ray_eq_pis2() ; 
//	double ydes_max =  2.5 * star(2).ray_eq_pis2() ; 
//
//	des_coupe_bin_z(star(1).get_nbar()(), star(2).get_nbar()(), 0, 
//		    xdes_min, xdes_max, ydes_min, ydes_max, "density (z=0)",  
//		    &(star(1).get_ent()()), &(star(2).get_ent()()) ) ; 
//
//	des_coupe_bin_z(star(1).get_logn_auto()(), star(2).get_logn_auto()(), 0, 
//		    xdes_min, xdes_max, ydes_min, ydes_max, "ln(N) (z=0)",  
//		    &(star(1).get_ent()()), &(star(2).get_ent()()) ) ; 
//
//	des_coupe_z(star(1).get_logn_comp()(), 0., 1, "logn_comp (z=0)", 
//##		&(star(1).get_ent()()) ) ; 

    }

    // Cleaning
    // --------

    delete peos1 ;    
    delete peos2 ;    

    return EXIT_SUCCESS ; 

}
