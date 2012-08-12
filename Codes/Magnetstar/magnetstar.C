/*
 * Main code for computing magnetized axisymmetric stars without deformations. 
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
 * Revision 1.1  2012/08/12 17:57:41  p_cerda
 * Magnetstar code
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

// headers C
#include <stdlib.h>
#include <math.h>
#include <string.h>

// headers Lorene
#include "et_rot_mag.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"
#include "metric.h"	    

// Local prototype (for drawings only)
Cmp raccord_c1(const Cmp& uu, int l1) ; 

// Definition de f_j et de M_j

Cmp f_j(const Cmp& x, const double a_j){ 
  Cmp resu(x.get_mp() ) ;
  resu = a_j ;
  return resu ;
}
Cmp M_j(const Cmp& x, const double a_j){ 
  return - a_j*x ;
}

Cmp fn_j(const Cmp& x, const Tbl an_j){ 
  Cmp resu(x.get_mp() ) ;
  Cmp xsign(x.get_mp() ) ;
  resu = 0.;
  assert ( an_j.get_dim(0) == 7);  
  
  xsign = x / sqrt(x*x);

  // M = a0 * x + a1 * x^2 + a2 * x^3 ...
  // f = dM/dx = a0 + 2 * a1 * x + 3 * a2 * x^2 ...
  resu = resu + an_j(0);  
  resu = resu + 2 * an_j(1) * x;
  resu = resu + 3 * an_j(2) * x * x;

  // M = c0_j * log (1 + x / rc0_j) 
  // f = dM/dx = c0_j / (1 + x/rc0_j) / rc0_j 
  resu = resu + an_j (3) / (1.0 + x / an_j (4)) / an_j (4);

  // M = c1_j * rc1_j * (x/rc1_j - atan(x / rc1_j))
  // f = dM/dx = c1_j * (1 - 1 / (1 + (x/rc1_j)^2))
  resu = resu + an_j (5) * (1.0 - 1.0 / (x*x/an_j(6)/an_j(6) + 1.0));

  
  return resu ;
}
Cmp Mn_j(const Cmp& x, const Tbl an_j){ 
  Cmp resu(x.get_mp() ) ;
  resu = 0.;
  assert ( an_j.get_dim(0) == 7);

  // M =  a0 * x + a1 * x^2 + a2 * x^3 ...
  resu = resu - an_j(0) * x;  
  resu = resu - an_j(1) * x * x ;  
  resu = resu - an_j(2) * x * x * x;
  
  // M = c0_j * log (1 + x / rc0_j)  
  resu = resu - an_j (3)*log(1.0+x/an_j(4));

  // M = c1_j * rc1_j * (x/rc1_j - atan(x / rc1_j))
  resu = resu - an_j (5) * an_j (6) * ( (x/an_j (6)) - atan(x/an_j(6)));
  
  return resu ;
}

Cmp gn_j(const Cmp& x, const Tbl bn_j){ 
  Cmp resu(x.get_mp() ) ;
  Cmp N(x.get_mp() ) ;
  Cmp dNdx(x.get_mp() ) ;
  
  N = 0.;
  N = N + bn_j(0) * x ;
  N = N + bn_j(1) * x * x ;
  N = N + bn_j(2) * x * x * x;

  dNdx= 0.;
  dNdx = dNdx + bn_j (0);
  dNdx = dNdx + 2. * bn_j (1) * x;
  dNdx = dNdx + 3. * bn_j (2) * x * x;
    
  resu = N*dNdx;
  return resu ;
}
Cmp Nn_j(const Cmp& x, const Tbl bn_j){ 
  Cmp resu(x.get_mp() ) ;

  resu = 0.;
  resu = resu + bn_j(0) * x ;
  resu = resu + bn_j(1) * x * x ;
  resu = resu + bn_j(2) * x * x * x;
  
  return resu ;
}



//******************************************************************************

int main(){

  using namespace Unites_mag ;

    // Identification of all the subroutines called by the code : 
    
    system("ident magnetstar > identif.d") ; 

    // For the display : 
    //    char display_bold[]="x[1m" ; display_bold[0] = 27 ;

    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    char blabla[120] ;

    int relat_i, mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
	delta_mer_kep, mer_mass, mermax_poisson, graph, nz, nzet, nzadapt,
	nt, np ; 
    double ent_c, freq_si, fact_omega, mbar_wanted, precis, freq_ini_si, 
	   thres_adapt, aexp_mass, relax, relax_poisson, precis_adapt ;  
   
    double a0_j, a1_j, a2_j, b0_j, b1_j, b2_j, c0_j, c1_j;
    double rc0_j, rc1_j;

    double relax_mag;
    int initial_j;

    ifstream fich("parrot.d") ;
    fich.getline(blabla, 120) ;
    fich >> relat_i ; fich.getline(blabla, 120) ;
    bool relat = (relat_i == 1) ; 
    fich >> ent_c ; fich.getline(blabla, 120) ;
    fich >> freq_si ; fich.getline(blabla, 120) ;
    fich >> fact_omega ; fich.getline(blabla, 120) ;
    fich >> mbar_wanted ; fich.getline(blabla, 120) ;
    mbar_wanted *= msol ; 
    fich.getline(blabla, 120) ;
    fich >> a0_j ; fich.getline(blabla,120) ;
    fich >> a1_j ; fich.getline(blabla,120) ;
    fich >> a2_j ; fich.getline(blabla,120) ;
    fich >> c0_j ; fich.getline(blabla,120) ;
    fich >> rc0_j ; fich.getline(blabla,120) ;
    //    fich >> c1_j ; fich.getline(blabla,120) ;
    //    fich >> rc1_j ; fich.getline(blabla,120) ;
    c1_j = 0.;
    rc1_j = 1.;
    fich >> b0_j ; fich.getline(blabla,120);
    fich >> b1_j ; fich.getline(blabla,120);
    fich >> b2_j ; fich.getline(blabla,120);
    fich >> relax_mag ; fich.getline(blabla,120);
    fich >> initial_j ; fich.getline(blabla,120);
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

    Tbl an_j (7);
    an_j.set_etat_qcq();
    an_j.set(0) = a0_j;
    an_j.set(1) = a1_j;
    an_j.set(2) = a2_j;
    an_j.set(3) = c0_j;
    an_j.set(4) = rc0_j;
    an_j.set(5) = c1_j;
    an_j.set(6) = rc1_j;
    
    Tbl bn_j (3);
    bn_j.set_etat_qcq();
    bn_j.set(0) = b0_j;
    bn_j.set(1) = b1_j;
    bn_j.set(2) = b2_j;



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
    int type_t;
    if ( initial_j == 2 ) {
      type_t = NONSYM ;  
    }else{
      type_t = SYM ;  
    }
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
    
    Et_rot_mag star(mp, nzet, relat, eos, 1) ;  
    
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

    Itbl icontrol(10) ;
    icontrol.set_etat_qcq() ; 
    icontrol.set(0) = mer_max ; 
    icontrol.set(1) = mer_rot ; 
    icontrol.set(2) = mer_change_omega ; 
    icontrol.set(3) = mer_fix_omega ; 
    icontrol.set(4) = mer_mass ; 
    icontrol.set(5) = mermax_poisson ;  
    icontrol.set(6) = delta_mer_kep ; 
    icontrol.set(7) = 0 ;
    icontrol.set(8) = 0 ;
    icontrol.set(9)= 1 ;
    
    Tbl control(8) ; 
    control.set_etat_qcq() ; 
    control.set(0) = precis ; 
    control.set(1) = omega_ini ; 
    control.set(2) = relax ; 
    control.set(3) = relax_poisson ; 
    control.set(4) = thres_adapt ; 
    control.set(5) = precis_adapt ; 
    control.set(6) = 0. ;
    control.set(7) = 0. ;

    Tbl diff(1) ;     
    

    star.equilibrium_mag(ent_c, omega, fact_omega, nzadapt, ent_limit, 
    			 icontrol, control, mbar_wanted, aexp_mass, diff, 
    			 0., 0., &f_j, &M_j) ;

    cout << endl << "Non-magnetized star : " 
	 << endl << "==========   " << endl ;
    
    
    star.equilibrium_mag_plus(icontrol, control, diff, initial_j, an_j, &fn_j, &Mn_j, bn_j, &gn_j, &Nn_j,
			      relax_mag) ;

     
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
    cout << "max(B) in each domain          :  " << maxB ;
    cout << "div(B) / max(B) in each domain :  " << max(abs(Bmag.divergence(gam))) / maxB ; 
    
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
    system("ident magnetstar >> calcul.d") ; 


    // Saveguard of the whole configuration
	// ------------------------------------

	FILE* fresu = fopen("resu.d", "w") ;

	star.get_mp().get_mg()->sauve(fresu) ;		// writing of the grid
	star.get_mp().sauve(fresu) ;                // writing of the mapping
	star.get_eos().sauve(fresu) ;  				// writing of the EOS
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

	//	int nzdes = star.get_nzet() ; 

	//	des_coupe_y(star.get_At(), 0., nzdes, "A\\dt\\u Potential", &surf,
	//		    2.) ;
	/*
	for (int i = 1; i<mg.get_nzone()-1; i++) { 
	  double rmin = mp.val_r_jk(i, 1, 0, 0);
	  cout <<"rmin " << rmin << endl ;
	  	  des_coupe_y(star.get_Aphi(), 0., i, "Magnetic field", 0x0, rmin, true, 20,100,100) ; 
	}
	*/

	des_coupe_y(star.get_Aphi(), 0., 1, "Magnetic field lines (A_phi)", 0x0, 1.2, true, 20,100,100) ; 
	des_coupe_y(star.get_jphi(), 0., 1, "Current", 0x0, 1.2, true, 20,100,100) ; 
	des_coupe_y(star.get_Bphi(), 0., 1, "Toroidal magnetic field", 0x0, 1.2, true, 20,100,100) ; 

	//	des_coupe_y(star.get_ent()(), 0., nzdes, "Enthalpy", &surf) ; 

	//strcpy(title, "Gravitational potential \\gn") ; 

	//des_coupe_y(star.get_logn()(), 0., nzdes, title, &surf) ; 
	
	//char title[80] ;
	
	//strcpy(title, "Azimuthal shift N\\u\\gf") ; 
	//des_coupe_y(star.get_nphi()(), 0., nzdes, title, &surf) ; 
	
	//strcpy(title, "Metric potential \\gz") ; 
	//des_coupe_y(star.get_dzeta()(), 0., nzdes, title, &surf) ; 
	
	//strcpy(title, "Metric potential (NB-1) r sin\\gh") ; 
	//des_coupe_y(star.get_tggg()(), 0., nzdes, title, &surf) ; 
	
	//strcpy(title, "A\\u2\\dK\\uij\\dK\\dij\\u") ; 
	//des_coupe_y(star.get_ak_car()(), 0., nzdes, title, &surf) ; 

//output of function M
	if(1)
	  {
	   
	    Tbl maxA_phi = max(abs(star.get_Aphi()));
	    Cmp A_phi_scaled = star.get_Aphi() / maxA_phi.set(0);

	    des_coupe_y(-Mn_j(A_phi_scaled, an_j), 0.,
	    		1, "deformation M", 0x0, 1.2, true, 20,100,100) ;

	    /*
	    FILE* opf_M;
	    double rr, theta, phi;
	    int ii, jj;
	    opf_M=fopen("M.dat","w");
	    cout << maxA_phi.set(0) << endl;
	    for(ii=1;ii<=100;ii++)
	      {
		for(jj=0;jj<=50;jj++)
		  {
		    rr = 1.2e0 * ii/ 100.0;
		    theta = 3.14 * jj / 50.0;
		    phi = 0.0;
		    fprintf(opf_M,"%e %e %e\n",rr*sin(theta),rr*cos(theta),
			    -float(Mn_j(A_phi_scaled,an_j).val_point(rr, theta, phi)));
		  }
		fprintf(opf_M,"\n");
	      }
	    fclose(opf_M);
	    */
	  }

    }

 
    // Cleaning
    // --------

    delete peos ;    

    exit(EXIT_SUCCESS) ; 
    
    return EXIT_SUCCESS ; 
   
}
