/*
 * Main code for computing a sequence of stationary axisymmetric differentially
 * rotating stars
 */

/*
 *   Copyright (c) 2001-2003 Eric Gourgoulhon
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

char rotseq_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.9  2014/10/06 15:09:45  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2004/03/25 12:35:44  j_novak
 * now using namespace Unites
 *
 * Revision 1.7  2003/09/16 13:06:24  e_gourgoulhon
 * Replaced the fich.getline(blabla, 120) by fich.ignore(1000, '\n')
 *
 * Revision 1.6  2003/08/26 08:58:50  e_gourgoulhon
 *
 * Added M/R and quantities in polytropic units in the
 * output file res.d
 *
 * Revision 1.5  2003/05/25 19:56:49  e_gourgoulhon
 *
 * Added the possibility to choose the factor a = R_eq / R0, instead of R0
 * in the differential rotation law.
 *
 * Revision 1.4  2003/05/14 20:06:09  e_gourgoulhon
 * Suppressed the qualifier ios::nocreate in call to fstream::open
 * (not supported by gcc 3.2).
 *
 * Revision 1.3  2003/01/09 11:07:50  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.2  2002/03/27 22:00:23  e_gourgoulhon
 * Now can computes sequences of rigidly rotating stars (class Etoile_rot)
 * as well as differentially rotating stars (class Et_rot_diff)
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2001/10/26  17:02:18  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */


// headers C
#include <cstdlib>
#include <cmath>
#include <cstring>

// headers Lorene
#include "et_rot_diff.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"

// Function defining the rotation profile
double frotlin(double omega, const Tbl& par) ;
double primfrotlin(double omega, const Tbl& par) ;

// Local prototype (for drawings only)
Cmp raccord_c1(const Cmp& uu, int l1) ;

//******************************************************************************


int main(){

  using namespace Unites ;

    // For the display :
    char display_bold[]="x[1m" ; display_bold[0] = 27 ;

    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    int relat_i, mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
	delta_mer_kep, mer_mass, mermax_poisson, graph, nz, nzet, nzadapt,
	nt, np, mer_triax, n_conf, diffrot_i ;
    double entc_min, entc_max, fact_omega, mbar_wanted, precis, freq_min_si,
	    freq_max_si, thres_adapt, aexp_mass, relax, relax_poisson,
	    ampli_triax, precis_adapt, rrot_km, arot ;  
    
    ifstream fich("parrotseq.d") ;
	if ( !fich.good() ) {
		cout << "Problem with opening the file parrotseq.d ! " << endl ;
		abort() ;
	}
	
	fich.ignore(1000,'\n') ;
    fich >> relat_i ; fich.ignore(1000,'\n') ;
    bool relat = (relat_i == 1) ; 
    fich >> entc_min ; fich.ignore(1000,'\n') ;
    fich >> entc_max ; fich.ignore(1000,'\n') ;
    fich >> freq_min_si ; fich.ignore(1000,'\n') ;
    fich >> freq_max_si ; fich.ignore(1000,'\n') ;
    fich >> n_conf ; fich.ignore(1000,'\n') ;    
    if (n_conf > 1000) {
	cout << "rotseq: n_conf must be smaller than 1000" << endl ;
        abort() ; 
    }

    fich >> diffrot_i ;  fich.ignore(1000,'\n');
	if (diffrot_i <= 1) {
    	fich >> rrot_km ; fich.ignore(1000,'\n');
	arot = 0 ; 
    }
    else {
    	assert (diffrot_i == 2) ; 
    	fich >> arot ; fich.ignore(1000,'\n');
	rrot_km = 0 ; 	
    }
	
    fact_omega = 1 ; 
    fich >> mbar_wanted ; fich.ignore(1000,'\n');
    mbar_wanted *= msol ; 
    fich.ignore(1000,'\n');
    fich >> mer_max ; fich.ignore(1000,'\n');
    fich >> precis ; fich.ignore(1000,'\n');
    mer_rot = 0 ; 
    mer_change_omega = 0 ; 
    mer_fix_omega = 1 ; 
    delta_mer_kep = 0 ; 
    fich >> thres_adapt ; fich.ignore(1000,'\n');
    mer_triax = 100000 ; 
    ampli_triax = 0 ; 
    fich >> mer_mass ; fich.ignore(1000,'\n');
    fich >> aexp_mass ; fich.ignore(1000,'\n');
    fich >> relax ; fich.ignore(1000,'\n');
    fich >> mermax_poisson ; fich.ignore(1000,'\n');
    fich >> relax_poisson ; fich.ignore(1000,'\n');
    fich >> precis_adapt ; fich.ignore(1000,'\n');
    fich >> graph ; fich.ignore(1000,'\n');
    fich.ignore(1000,'\n');
    fich >> nz ; fich.ignore(1000,'\n') ;
    fich >> nzet; fich.ignore(1000,'\n') ;
    fich >> nzadapt; fich.ignore(1000,'\n') ;
    fich >> nt; fich.ignore(1000,'\n') ;
    fich >> np; fich.ignore(1000,'\n') ;

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
     
    fich.ignore(1000,'\n') ;
    for (int l=0; l<nz; l++) {
	fich >> nr[l]; 
	fich >> bornes[l]; fich.ignore(1000,'\n');
	np_tab[l] = np ; 
	nt_tab[l] = nt ; 
    }
    bornes[nz] = __infinity ;

    Tbl ent_limit(nzet) ;
    ent_limit.set_etat_qcq() ;
    ent_limit.set(nzet-1) = 0 ; 	// enthalpy at the stellar surface
    for (int l=0; l<nzet-1; l++) {
    	fich >> ent_limit.set(l) ; fich.ignore(1000,'\n');
    }


    fich.close();

    
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
    	       cout << "rotseq: problem : peos is not of type Eos_strange_cr !" << endl ;
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
    //		Parameters for the function defining the differential rotation
    //-----------------------------------------------------------------------
    
    double omega_c_min = 2 * M_PI * freq_min_si / f_unit ; 
    double omega_c_max = 2 * M_PI * freq_max_si / f_unit ; 

    double rrot = rrot_km * km ; 

    Tbl parfrot(3) ;
    parfrot.set_etat_qcq() ; 
    parfrot.set(0) = 0 ;  
    parfrot.set(1) = rrot ;  
    parfrot.set(2) = arot ;      
    
    //-----------------------------------------------------------------------
    //		Construction of the star
    //-----------------------------------------------------------------------

    bool diffrot = (diffrot_i >= 1) ;
    Etoile_rot* p_star ;
    Et_rot_diff* p_star_diff = 0x0 ;

    if ( diffrot ) {
        p_star_diff = new Et_rot_diff(mp, nzet, relat, eos, frotlin, primfrotlin, parfrot) ;
    	p_star = p_star_diff ;
    }
    else {
        p_star = new Etoile_rot(mp, nzet, relat, eos) ;
    }
    Etoile_rot& star = *p_star ;

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
    ent0 = entc_min * ( 1 - r*r / (ray0*ray0) ) ;
    ent0.annule(nz-1) ;
    ent0.std_base_scal() ;
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
    "   PARAMETERS USED FOR THE COMPUTATION (file parrotseq.d) : " << endl ;
    fichseq <<
    "================================================================" << endl ;
    fichseq.close() ;
    system("cat parrotseq.d >> calcul_seq.d") ;

    fichseq.open("calcul_seq.d", ios::app) ;
    fichseq << endl <<
    "================================================================" << endl ;
    fichseq <<
    "	           EOS PARAMETERS (file par_eos.d) : " << endl ;
    fichseq <<
    "================================================================" << endl ;
    fichseq.close() ;
    system("cat par_eos.d >> calcul_seq.d") ;

    // Identification du code et de ses sous-routines (no. de version RCS) :
    fichseq.open("calcul_seq.d", ios::app) ;
    fichseq << endl <<
    "================================================================" << endl ;
    fichseq << "	    IDENTIFICATION OF THE CODE : " << endl ;
    fichseq <<
    "================================================================" << endl ;
    fichseq.close() ;
    system("ident rotseq >> calcul_seq.d") ;



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

    ofstream fichresu("seq.d") ;
    fichresu << "# M/R  M [poly]  M_B [poly]  J [G M_sol^2/c] M [M_sol]   f_c [Hz]     H_c     r_p/r_e"
             << "       T/W       M_B [M_sol]    GRV2       GRV3   " << endl ;

    // Loop on the configurations
    // --------------------------

    bool seq_freq = ( entc_min == entc_max ) ;
    if ( !seq_freq ) {
          if ( omega_c_min != omega_c_max ) {
          	cout << "rotseq : one must have freq_min == freq_max "
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

		ent_c = star.get_ent()()(0, 0, 0, 0) ;
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


    	const Eos_poly* p_eos_poly = dynamic_cast<const Eos_poly*>(
		&(star.get_eos()) ) ; 	  

    	double mass_g_poly, mass_b_poly ; 
	if (p_eos_poly != 0x0) {

		double kappa = p_eos_poly->get_kap() ; 
		double gamma = p_eos_poly->get_gam() ;  ; 

		// kappa^{n/2}
		double kap_ns2 = pow( kappa,  0.5 /(gamma-1) ) ; 
    
		// Polytropic unit of length in terms of r_unit : 
		double r_poly = kap_ns2 / sqrt(ggrav) ; 
    
		// Polytropic unit of mass in terms of m_unit :
		double m_poly = r_poly / ggrav ; 
    
		mass_g_poly = star.mass_g() / m_poly ;
		mass_b_poly = star.mass_b() / m_poly ;

	}
	else{
		mass_g_poly = 0 ;  
		mass_b_poly = 0 ;  
	}
	
	
	int precisaff = 8 ;
	int tailleaff = precisaff + 3 ;
	fichresu.precision(precisaff) ;
	fichresu << setw(tailleaff)
		 << ggrav * star.mass_g() / star.r_circ() << " "
		 << setw(tailleaff)
		 << mass_g_poly << " "
		 << setw(tailleaff)
		 << mass_b_poly << " "
		 << setw(tailleaff)
		 << star.angu_mom()/( qpig / (4* M_PI) * msol*msol) << " "
	         << setw(tailleaff)
		 << star.mass_g() / msol << " "
		 << setw(tailleaff)
		 << star.get_omega_c() / (2.*M_PI) * f_unit << " "
		 << setw(tailleaff)
		 << star.get_ent()()(0, 0, 0, 0) << " "
		 << setw(tailleaff)
		 << star.aplat() << " "
		 << setw(tailleaff)
		 << star.tsw() << " "
		 << setw(tailleaff)
		 << star.mass_b() / msol << " "
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

	char nomfich[20] ;
	strcpy(nomfich, "calcul") ;

	char numero[3] ;
	sprintf(numero, "%3.3d", jj) ;
	strcat(nomfich, numero) ;
	strcat(nomfich, ".d") ;
	cout << endl << "File name : " << nomfich << endl ;

	ofstream fichfinal(nomfich) ;
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

	// Cmp defining the surface of the star (via the enthalpy field)
	Cmp surf = star.get_ent()() ; 
	Cmp surf_ext(mp) ; 
	surf_ext = - 0.2 * surf(0, 0, 0, 0) ; 
	surf_ext.annule(0, star.get_nzet()-1) ; 
	surf.annule(star.get_nzet(), mg.get_nzone()-1) ; 
	surf = surf + surf_ext ;
	surf = raccord_c1(surf, star.get_nzet()) ; 

	int nzdes = star.get_nzet() ; 

	des_coupe_y(star.get_ent()(), 0., nzdes, "Enthalpy", &surf) ; 
	
	if (mer_triax < mer_max) { 
	    des_coupe_z(star.get_ent()(), 0., nzdes, "Enthalpy (equatorial plane)",
			&surf) ;
	}

	if (diffrot) {
		Cmp tmpdes = p_star_diff->get_omega_field()() / (2*M_PI) * f_unit ;
		des_profile(tmpdes, 0., star.ray_eq(),
			    M_PI/2., 0., "\\gW/2\\gp  [Hz]",
			"Angular velocity in equatorial plane") ;
         	
		des_coupe_y(p_star_diff->get_omega_field()(), 0., nzdes,
				"\\gW", &surf) ;
	}

	des_coupe_y(star.get_logn()(), 0., nzdes,
		    "Gravitational potential \\gn", &surf) ;

	if (star.is_relativistic()) {

	    des_coupe_y(star.get_nphi()(), 0., nzdes,
		    "Azimuthal shift N\\u\\gf", &surf) ;

	    des_coupe_y(star.get_dzeta()(), 0., nzdes,
		    "Metric potential \\gz", &surf) ;

	    des_coupe_y(star.get_tggg()(), 0., nzdes,
		    "Metric potential (NB-1) r sin\\gh", &surf) ;
	
	    des_coupe_y(star.get_ak_car()(), 0., nzdes, 
		    "A\\u2\\d K\\dij\\u K\\uij\\d", &surf) ; 
	}
	
	}


    } // End of the loop on the configurations

    fichresu.close() ; 
 
    // Cleaning
    // --------

    delete p_star ;

    delete peos ;

    exit(EXIT_SUCCESS) ; 
    
    return EXIT_SUCCESS ; 
   
}

		//------------------------------//
		//	Function F(Omega)	//
		//------------------------------//

double frotlin(double omega, const Tbl& par){
    
	double omega_c = par(0) ; 
	double rrot = par(1) ; 
    
	return rrot*rrot* (omega_c - omega) ;  
    
}

		//------------------------------//
		//	Primitive of  F(Omega)	//
		//------------------------------//

double primfrotlin(double omega, const Tbl& par){
    
	double omega_c = par(0) ; 
	double rrot = par(1) ; 
    
	return - 0.5 * rrot*rrot* (omega_c - omega)*(omega_c - omega) ;  
    
}
