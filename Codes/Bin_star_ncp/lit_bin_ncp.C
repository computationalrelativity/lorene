/*
 * Reads a file containing a binary configuration (class Bin_ns_ncp) and
 * performs various plots. 
 * 
 */

/*
 *   Copyright (c) 2003 Francois Limousin
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

char lit_bin_ncp_C[] = "$Header$" ;

/*
 * $Header$
 *
 */
// headers C
#include <stdlib.h>
#include <math.h>
#include <string.h>

// headers Lorene
#include "bin_ns_ncp.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"

// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ; 

//******************************************************************************

int main(int argc, char** argv){

    // Identification of all the subroutines called by the code : 
    
    // system("ident lit_bin") ; 
    
    if (argc < 2) {
	cout << 
	"lit_bin_ncp : the name of a file containing a binary configuration"
	<< endl << " must be given in argument !" << endl ; 
	abort() ; 
    }
    
    char* nomresu = argv[1] ; 
    cout << "Name of the file to be read : " << nomresu << endl ; 
    
    cout << endl << 
    "Do you want to draw the boundaries of the various domains (y/n) ? [y]"
	 << endl ; 
    char rep ; 
    cin.get(rep) ;
    bool draw_bound = !(rep == 'n') ; 
        
    #include "unites.h"	    
    // To avoid some compilation warnings
    if (nomresu == 0x0) {
	cout << qpig << f_unit << mevpfm3 << endl ; 
    }    
    
  
    FILE* fich = fopen(nomresu, "r") ; 

    int mer ; 
    fread(&mer, sizeof(int), 1, fich) ;	// mer
    
    Mg3d mg1(fich) ;
    Map_et mp1(mg1, fich) ; 
    Eos* peos1 = Eos::eos_from_file(fich) ; 
    
    Mg3d mg2(fich) ;
    Map_et mp2(mg2, fich) ; 
    Eos* peos2 = Eos::eos_from_file(fich) ; 

    Tenseur_sym plat1(mp1, 2, COV, mp1.get_bvect_cart() ) ;
    plat1.set_etat_qcq() ;
    for (int i=0; i<3; i++) {
      for (int j=0; j<i; j++) {
	plat1.set(i,j) = 0 ;
      }
      plat1.set(i,i) = 1 ;
    }
    plat1.set_std_base() ;

    Metrique flat1(plat1, true) ; 
    
    Tenseur_sym plat2(mp2, 2, COV, mp2.get_bvect_cart() ) ;
    plat2.set_etat_qcq() ;
    for (int i=0; i<3; i++) {
      for (int j=0; j<i; j++) {
	plat2.set(i,j) = 0 ;
      }
      plat2.set(i,i) = 1 ;
    }
    plat2.set_std_base() ;

    Metrique flat2(plat2, true) ; 
   
    Bin_ns_ncp star(mp1, *peos1, mp2, *peos2, flat1, flat2, fich) ; 

    fclose(fich) ; 
    
    bool relativistic = star(1).is_relativistic() ; 
    bool irrotational = star(1).is_irrotational() ; 

    cout << endl << "Grid on which star 1 is defined : " << endl ; 
    cout << "=============================== " << endl ; 
    cout << *((star(1).get_mp()).get_mg()) << endl ; 

    cout << endl << "Grid on which star 2 is defined : " << endl ; 
    cout << "=============================== " << endl ; 
    cout << *((star(2).get_mp()).get_mg()) << endl ; 

    cout << endl << "Mapping on which star 1 is defined : " << endl ; 
    cout << "================================== " << endl ; 
    cout << star(1).get_mp() << endl ; 

    des_map_et(mp1, 0) ; 
    des_map_et(mp1, 1) ; 

    cout << endl << "Mapping on which star 2 is defined : " << endl ; 
    cout << "================================== " << endl ; 
    cout << star(2).get_mp() << endl ; 

    for (int i=1; i<=2; i++) {
	(star.set(i)).update_metric(star(3-i)) ; 
    }

    for (int i=1; i<=2; i++) {
	(star.set(i)).update_metric_der_comp(star(3-i)) ; 
    }

    for (int i=1; i<=2; i++) {
	(star.set(i)).equation_of_state() ; 
	(star.set(i)).kinematics(star.get_omega(), star.get_x_axe()) ; 
	(star.set(i)).fait_d_psi() ; 
	(star.set(i)).hydro_euler() ; 
    }

    cout << "Binary system read in file : " << endl ;
    cout << star << endl ; 
    
    /*
    cout << "ADM mass [M_sol] : " << star.mass_adm() / msol  << endl ; 
    cout << "Total energy [M_sol c^2] : " 
	 << star.total_ener() / msol << endl ; 
    cout << "Total angular momentum [M_sol c km] : " 
	 << (star.angu_mom())(2) / msol / km << endl ; 
    if (!relativistic) {
	cout << "Relative error on the virial theorem : " 
	     << star.virial() << endl ; 
    }
    
    cout << "Relative error in the Hamiltonian constraint : " << endl ; 
    cout << star.ham_constr() << endl ; 
	 
    cout << "Relative error in the momentum constraint : " << endl ; 
    cout << " X component : " << star.mom_constr()(0) << endl ; 
    cout << " Y component : " << star.mom_constr()(1) << endl ; 
    cout << " Z component : " << star.mom_constr()(2) << endl ; 
    */

    //==============================================================
    //  Reduced quantities for polytropic EOS
    //==============================================================

    star.display_poly(cout) ; 

    //==============================================================
    //  Drawings
    //==============================================================

    des_explorer_symz(star, "latbin") ; 
    
    int nzdes1 = star(1).get_nzet() ; 
    
    double ori_x1 = star(1).get_mp().get_ori_x() ; 
    double ori_x2 = star(2).get_mp().get_ori_x() ; 

    double xdes_min = - 1.5 * star(1).ray_eq_pi() + ori_x1 ;
    xdes_min += 0.2 * xdes_min ;  
    double xdes_max = 1.5 * star(2).ray_eq_pi() + ori_x2 ; 
    xdes_max += 0.2 * fabs(xdes_min) ;  

    double ydes_min1 = - 4. * star(1).ray_eq_pis2() ; 
    double ydes_min2 = - 4. * star(2).ray_eq_pis2() ; 
    double ydes_min = (ydes_min1 < ydes_min2) ? ydes_min1 : ydes_min2 ; 

    double ydes_max1 =  4. * star(1).ray_eq_pis2() ; 
    double ydes_max2 =  4. * star(2).ray_eq_pis2() ; 
    double ydes_max = (ydes_max1 > ydes_max2) ? ydes_max1 : ydes_max2 ; 

    double zdes_min1 = - 4. * star(1).ray_pole() ; 
    double zdes_min2 = - 4. * star(2).ray_pole() ; 
    double zdes_min = (zdes_min1 < zdes_min2) ? zdes_min1 : zdes_min2 ; 

    double zdes_max1 =  4. * star(1).ray_pole() ; 
    double zdes_max2 =  4. * star(2).ray_pole() ; 
    double zdes_max = (zdes_max1 > zdes_max2) ? zdes_max1 : zdes_max2 ; 

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

    char title[80] ;
    char bslash[2] = {92,  '\0'} ;  // 92 is the ASCII code for backslash 

    ofstream fent("enthalpy.d") ; 
    
    fent << "Enthalpy field at the boundary of last inner domain of star 1 : " 
	 << endl ; 
    int lzet =  star(1).get_nzet() - 1 ; 
    for (int k=0; k<mg1.get_np(lzet); k++) {
	fent << "k = " << k << " : " ; 
	for (int j=0; j<mg1.get_nt(lzet); j++) {
	    fent 
		<< "  " << star(1).get_ent()()(lzet, k, j, mg1.get_nr(lzet)-1) ;
	}
	fent << endl ; 
    }
   
    fent << endl << "enthalpy field of star 1 : " << endl ;     
    fent << star(1).get_ent() << endl ; 
    fent.close() ; 
    
    des_coupe_z(star(1).get_ent()(), 0., 1,
		"Enthalpy (z=0)", &surf1, 1.2, draw_bound ) ; 

    des_coupe_y(star(1).get_ent()(), 0., 1,
		"Enthalpy (y=0)", &surf1, 1.2, draw_bound ) ; 


    des_profile (star(1).get_ent()(), 0., 2* star(1).ray_eq(), 0., 0.,  
	"H", "H (theta=0)" ) ; 

    des_profile (star(1).get_ent()(), 0., 2* star(1).ray_eq(), M_PI/2., M_PI/2.,  
	"H", "H (theta=pi/2,  phi=pi/2)" ) ; 



    //==========================================
    // Metric quantities
    //==========================================

    //----------------------------
    // ln(N)
    //----------------------------

    Tenseur logn1 = - star(1).get_logn_auto() ; 
    Tenseur logn2 = - star(2).get_logn_auto() ; 

    des_coupe_bin_y(logn1(), logn2(), 0, 
			xdes_min, xdes_max, zdes_min, zdes_max, 
		    "ln(N) (y=0)",  &surf1, &surf2, draw_bound ) ; 

    des_coupe_bin_z(logn1(), logn2(), 0, 
			xdes_min, xdes_max, ydes_min, ydes_max, 
		    "ln(N) (z=0)",  &surf1, &surf2, draw_bound ) ; 

    des_coupe_bin_x(logn1(), logn2(),
			ori_x1, ydes_min, ydes_max, zdes_min, zdes_max, 
		    "ln(N) (x=x1)",  &surf1, 0x0, draw_bound ) ; 

   
    if (relativistic) {

	//--------------
	// Shift vector
	//--------------
      
      Tenseur shift2_auto = star(2).get_shift_auto() ;
      shift2_auto.change_triad(mp1.get_bvect_cart()) ;

      des_vect_bin_z(star(1).get_shift_auto(), shift2_auto, 0., 
		    -2., 0.5, xdes_min, xdes_max, ydes_min, ydes_max,
		    "Shift vector  (z=0)", 
		    &surf1, &surf2, draw_bound ) ; 

	//----------------------------
	// Extrinsic curvature tensor
	//----------------------------
    
	Tenseur tkij1 = star(1).get_tkij_auto() ;
	Tenseur tkij2 = star(2).get_tkij_auto() ;
    
	// Division by r^2 in the external compactified domain in order to get
	//  K^{ij} : 
	tkij1.dec2_dzpuis() ;     
	tkij2.dec2_dzpuis() ;     
    
	char debtit[] = {'K', 92, 'u', '\0'} ; 

	strcpy(title, debtit) ; 
	strcat(title, "xx") ; 
	strcat(title, bslash) ; 
	strcat(title, "d (z=0)") ; 

	des_coupe_bin_z(tkij1(0, 0), tkij2(0, 0), 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 

	strcpy(title, debtit) ; 
	strcat(title, "xy") ; 
	strcat(title, bslash) ; 
	strcat(title, "d (z=0)") ; 
	
	des_coupe_bin_z(tkij1(0, 1), tkij2(0, 1), 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 

	strcpy(title, debtit) ; 
	strcat(title, "xz") ; 
	strcat(title, bslash) ; 
	strcat(title, "d (y=0)") ; 
	
	des_coupe_bin_y(tkij1(0, 2), tkij2(0, 2), 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 

	strcpy(title, debtit) ; 
	strcat(title, "yy") ; 
	strcat(title, bslash) ; 
	strcat(title, "d (z=0)") ; 
	
	des_coupe_bin_z(tkij1(1, 1), tkij2(1, 1), 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 

	strcpy(title, debtit) ; 
	strcat(title, "yz") ; 
	strcat(title, bslash) ; 
	strcat(title, "d (y=0)") ; 
	
	des_coupe_bin_y(tkij1(1, 2), tkij2(1, 2), 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 

	strcpy(title, debtit) ; 
	strcat(title, "zz") ; 
	strcat(title, bslash) ; 
	strcat(title, "d (z=0)") ; 
	
	des_coupe_bin_z(tkij1(2, 2), tkij2(2, 2), 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 


	//----------------------------
	// Term  K_ij K^ij
	//----------------------------
    
	Tenseur kcar1 = star(1).get_kcar_auto() ; 
	Tenseur kcar2 = star(2).get_kcar_auto() ; 

	// Division by r^4 in the external compactified domain in order to get
	//  A^2 K_ij K^{ij} : 
	kcar1.dec2_dzpuis() ; 
	kcar1.dec2_dzpuis() ; 
	kcar2.dec2_dzpuis() ; 
	kcar2.dec2_dzpuis() ; 
    
	char title1[80] ;
	strcpy(title1, debtit) ; 
	strcat(title1, "ij") ; 
	strcat(title1, bslash) ; 
	strcat(title1, "d K") ; 
	strcat(title1, bslash) ; 
	strcat(title1, "dij") ; 
	strcat(title1, bslash) ; 

	strcpy(title, title1) ; 
	strcat(title, "u (x=x1)") ; 
	
	des_coupe_bin_x(kcar1(), kcar2(), ori_x1, 
			ydes_min, ydes_max, zdes_min, zdes_max, 
			title,  &surf1, 0x0, draw_bound ) ; 

	strcpy(title, title1) ; 
	strcat(title, "u (y=0)") ; 
	
	des_coupe_bin_y(kcar1(), kcar2(), 0, 
			xdes_min, xdes_max, zdes_min, zdes_max, 
			title,  &surf1, &surf2, draw_bound ) ; 

	strcpy(title, title1) ; 
	strcat(title, "u (z=0)") ; 
	
	des_coupe_bin_z(kcar1(), kcar2(), 0, 
			xdes_min, xdes_max, ydes_min, ydes_max, 
			title,  &surf1, &surf2, draw_bound ) ; 
    
	//----------------------------
	// determinant gamma
	//----------------------------
    
	Tenseur gamma1 = - star(1).get_loggamma_auto() ; 
	Tenseur gamma2 = - star(2).get_loggamma_auto() ; 
	
	des_coupe_bin_x(gamma1(), gamma2(),
			ori_x1, ydes_min, ydes_max, zdes_min, zdes_max, 
		    "loggamma (x=x1)",  &surf1, 0x0, draw_bound ) ; 

	des_coupe_bin_y(gamma1(), gamma2(), 0, 
			xdes_min, xdes_max, zdes_min, zdes_max, 
		    "loggamma (y=0)",  &surf1, &surf2, draw_bound ) ; 

	des_coupe_bin_z(gamma1(), gamma2(), 0, 
			xdes_min, xdes_max, ydes_min, ydes_max, 
		    "loggamma (z=0)",  &surf1, &surf2, draw_bound ) ; 

    }

    //==========================================
    // Hydro quantities
    //==========================================


    des_coupe_bin_z(star(1).get_nbar()(), star(2).get_nbar()(), 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    "Baryon density (z=0)",  
		    &surf1, &surf2, draw_bound ) ; 

    des_coupe_bin_y(star(1).get_nbar()(), star(2).get_nbar()(), 0, 
		    xdes_min, xdes_max, zdes_min, zdes_max, 
		    "Baryon density (y=0)",  
		    &surf1, &surf2, draw_bound ) ; 

    des_coupe_z(star(1).get_nbar()(), 0., 1,
		"Baryon density (z=0)", &surf1, 1.2, draw_bound ) ; 

    des_coupe_y(star(1).get_nbar()(), 0., 1,
		"Baryon density (y=0)", &surf1, 1.2, draw_bound ) ; 

    if (irrotational) {
	Tenseur tmp_draw1 = star(1).get_wit_w() ; 
	tmp_draw1.annule(star(1).get_nzet(), mg1.get_nzone()-1) ; 
	Tenseur tmp_draw2 = star(2).get_wit_w() ; 
	tmp_draw2.annule(star(2).get_nzet(), mg2.get_nzone()-1) ; 

	des_vect_bin_z(tmp_draw1, tmp_draw2, 0., 
		    -3., 0.5, xdes_min, xdes_max, ydes_min, ydes_max,
		    "Velocity w.r.t corotating frame  (z=0)", 
		    &surf1, &surf2, draw_bound, 40, 40) ; 
    }
    		    
    Tenseur tmp_draw1 = star(1).get_u_euler() ; 
    tmp_draw1.annule(star(1).get_nzet(), mg1.get_nzone()-1) ; 

    Tenseur tmp_draw2 = star(2).get_u_euler() ; 
    tmp_draw2.annule(star(2).get_nzet(), mg2.get_nzone()-1) ; 
    tmp_draw2.change_triad(mp1.get_bvect_cart()) ;

    des_coupe_vect_x(tmp_draw1, mp1.get_ori_x(), -1., 0.5, nzdes1, 
		     "U (x=x1)", &surf1, 1.2, draw_bound ) ; 

    des_vect_bin_z(tmp_draw1, tmp_draw2, 0., 
		    -2., 0.5, xdes_min, xdes_max, ydes_min, ydes_max,
		    "U  (z=0)", &surf1, &surf2, draw_bound, 40, 40) ; 

    if (irrotational) {
	tmp_draw1 = star(1).get_d_psi() ; 
	tmp_draw1.annule(star(1).get_nzet(), mg1.get_nzone()-1) ; 
	des_coupe_vect_z(tmp_draw1, 0, -1., 0.5, nzdes1,
			 "Grad(psi) (z=0)", &surf1, 1.2, draw_bound ) ;

	des_coupe_z(star(1).get_psi0()(), 0., 1,
		    "psi0 (z=0)", &surf1, 1.2, draw_bound ) ; 
    
	Tenseur d_psi0 = star(1).get_psi0().gradient() ; 
	//##    d_psi0.change_triad(star.get_ref_triad()) ; 

	des_coupe_vect_z(d_psi0, 0, -3., 0.5, nzdes1, "Grad(psi0) (z=0)", 
			 &surf1, 1.2, draw_bound ) ; 

	tmp_draw1 = star(1).get_wit_w() ; 
	tmp_draw1.annule(star(1).get_nzet(), mg1.get_nzone()-1) ; 
	des_coupe_vect_z(tmp_draw1, 0, -3., 0.5, nzdes1, "W (z=0)", 
			 &surf1, 1.2, draw_bound ) ; 

    }

    // Cleaning
    // --------

    delete peos1 ;    
    delete peos2 ; 

    return EXIT_SUCCESS ;

}
