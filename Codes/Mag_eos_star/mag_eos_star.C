/*
 * Main code for computing magnetized stars with B-field dependent EoS. 
 * 
 */

/*
 *   Copyright (c) 2002 Jerome Novak & Emanuel Marcq
 *   Copyright (c) 2014 Debarati Chatterjee
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
 * Revision 1.6  2014/10/13 08:53:57  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/09/22 15:00:09  j_novak
 * Maximal magnetic field display.
 *
 * Revision 1.4  2014/07/04 13:01:47  j_novak
 * Parameter file for mass-shedding limit.
 *
 * Revision 1.3  2014/05/27 12:41:25  j_novak
 * Added the possibility to converge to a given mgnetic moment.
 *
 * Revision 1.2  2014/04/29 13:46:54  j_novak
 * Addition of swicthes 'use_B_in_eos' and 'include_magnetisation' to control the model.
 *
 * Revision 1.1  2014/04/28 14:41:15  j_novak
 * Files for magnetisation star main code.
 *
 *
 * $Header$
 *
 */

// headers C
#include <cmath>

// headers Lorene
#include "et_rot_mag.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"
#include "metric.h"	    

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

    char blabla[120] ;

    int relat_i, mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
      delta_mer_kep, mer_mass, mermax_poisson, graph, nz, nzet, nzadapt,
      nt, np ; 
    double ent_c, freq_si, fact_omega, mbar_wanted, precis, freq_ini_si, 
      thres_adapt, aexp_mass, relax, relax_poisson, precis_adapt ;  
    
    double Q0, a_j0, Q_ini, a_j_ini, magmom_wanted ;
    int mer_mag, mer_change_mag, mer_fix_mag, mag_in_eos, mer_magmom, 
      use_magnetisation ;

    ifstream fich("parrot.d") ;
    fich.getline(blabla, 120) ;
    fich >> relat_i ; fich.getline(blabla, 120) ;
    bool relat = (relat_i == 1) ; 
    fich >> ent_c ; fich.getline(blabla, 120) ;
    fich >> freq_si ; fich.getline(blabla, 120) ;
    fich >> fact_omega ; fich.getline(blabla, 120) ;
    fich >> mbar_wanted ; fich.getline(blabla, 120) ;
    mbar_wanted *= msol ; 
    fich >> magmom_wanted ; fich.getline(blabla, 120) ;
    fich.getline(blabla, 120) ;
    fich >> Q0 ; fich.getline(blabla,120) ;
    fich >> a_j0 ; fich.getline(blabla,120) ;
    fich >> Q_ini; fich.getline(blabla,120) ;
    fich >> a_j_ini ; fich.getline(blabla,120) ;
    fich >> mer_mag ; fich.getline(blabla,120) ;
    fich >> mer_change_mag ; fich.getline(blabla,120) ;
    fich >> mer_fix_mag ; fich.getline(blabla,120);
    fich >> mag_in_eos ; fich.getline(blabla,120);
    fich >> use_magnetisation ; fich.getline(blabla,120);
    fich.getline(blabla,120);
    fich >> mer_max ; fich.getline(blabla, 120) ;
    fich >> precis ; fich.getline(blabla, 120) ;
    fich >> mer_rot ; fich.getline(blabla, 120) ;
    fich >> freq_ini_si ; fich.getline(blabla, 120) ;
    fich >> mer_change_omega ; fich.getline(blabla, 120) ;
    fich >> mer_fix_omega ; fich.getline(blabla, 120) ;
    fich >> delta_mer_kep ; fich.getline(blabla, 120) ;
    fich >> thres_adapt ; fich.getline(blabla, 120) ;
    fich >> mer_mass ; fich.getline(blabla, 120) ;
    fich >> mer_magmom ; fich.getline(blabla, 120) ;
    fich >> aexp_mass ; fich.getline(blabla, 120) ;
    fich >> relax ; fich.getline(blabla, 120) ;
    fich >> mermax_poisson ; fich.getline(blabla, 120) ;
    fich >> relax_poisson ; fich.getline(blabla, 120) ;
    fich >> precis_adapt ; fich.getline(blabla, 120) ;
    fich >> graph ; fich.getline(blabla, 120) ;
    fich.getline(blabla, 120) ;
    fich >> nz ; fich.getline(blabla, 120) ;
    fich >> nzet; fich.getline(blabla, 120) ;
    fich >> nzadapt; fich.getline(blabla, 120) ;
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
    for (int l=0; l<nzet-1; l++) {
    	fich >> ent_limit.set(l) ; fich.getline(blabla, 120) ;
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
    if ( mer_magmom < mer_max ) {
	cout << "Required magnetic moment [SI] : " 
	     << magmom_wanted << endl ; 
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
    cout << "Step from which the magnetic moment is forced to converge : " 
	 << mer_magmom << endl ; 
    cout << "Exponent for the increase factor of the central enthalpy/CFA : " 
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
    
    Et_magnetisation star(mp, nzet, relat, eos, (use_magnetisation!=0), 
			  (mag_in_eos!=0) ) ; 
    
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
    icontrol.set(9)= mer_fix_mag ;
    icontrol.set(10) = mer_magmom ;
    
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
			 icontrol, control, mbar_wanted, magmom_wanted,
			 aexp_mass, diff, Q0, a_j0, &f_j, &M_j) ;

     
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

    Scalar norm_B_Euler 
      = sqrt( Scalar(star.get_a_car()())*( Bmag(1)*Bmag(1) + Bmag(2)*Bmag(2) ) ) ;
    norm_B_Euler.std_spectral_base() ;
    Tbl maxB = 0.5*(max(abs(Bmag(1))) + max(abs(Bmag(2)))) ;
    Scalar divB = Bmag.divergence(gam) ;
    cout << "div(B) / max(B) in each domain :  " 
	 << max(abs(divB)) / maxB ; 
    cout << "Maximal magnetic field in the Eulerian frame: " 
	 << max(max(norm_B_Euler))*mag_unit << " [T]" << endl ;

    // Saveguard of the whole configuration
    // ------------------------------------

    FILE* fresu = fopen("resu.d", "w") ;
    
    star.get_mp().get_mg()->sauve(fresu) ;	// writing of the grid
    star.get_mp().sauve(fresu) ;                // writing of the mapping
    star.get_eos().sauve(fresu) ;  		// writing of the EOS
    star.sauve(fresu) ;                         // writing of the star
    
    fclose(fresu) ;
    
	

    // Drawings
    // --------
    
    if (graph == 1) {

	char title[80] ;

	// Cmp defining the surface of the star (via the enthalpy field)
	Cmp surf = star.get_ent()() ; 
	Cmp surf_ext(mp) ; 
	surf_ext = - 0.2 * surf(0, 0, 0, 0) ; 
	surf_ext.annule(0, star.get_nzet()-1) ; 
	surf.annule(star.get_nzet(), mg.get_nzone()-1) ; 
	surf = surf + surf_ext ;
	surf = raccord_c1(surf, star.get_nzet()) ; 

	int nzdes = star.get_nzet() ; 

	des_coupe_y(star.get_At(), 0., nzdes, "A\\dt\\u Potential", &surf,
		    2.) ; 
	des_coupe_y(star.get_Aphi(), 0., nzdes, "Magnetic field", &surf) ; 
	des_coupe_y(norm_B_Euler, 0., nzdes, "Magnetic field norm", &surf) ; 
	des_coupe_y(star.get_ent()(), 0., nzdes, "Enthalpy", &surf) ; 

	strcpy(title, "Gravitational potential \\gn") ; 

	des_coupe_y(star.get_logn()(), 0., nzdes, title, &surf) ; 
	
	
	strcpy(title, "Azimuthal shift N\\u\\gf") ; 
	des_coupe_y(star.get_nphi()(), 0., nzdes, title, &surf) ; 
	
	strcpy(title, "Metric potential \\gz") ; 
	des_coupe_y(star.get_dzeta()(), 0., nzdes, title, &surf) ; 
	
	strcpy(title, "Metric potential (NB-1) r sin\\gh") ; 
	des_coupe_y(star.get_tggg()(), 0., nzdes, title, &surf) ; 
	
	strcpy(title, "A\\u2\\dK\\uij\\dK\\dij\\u") ; 
	des_coupe_y(star.get_ak_car()(), 0., nzdes, title, &surf) ; 

    }

 
    // Cleaning
    // --------

    delete peos ;    

    exit(EXIT_SUCCESS) ; 
    
    return EXIT_SUCCESS ; 
   
}
