/*
 * Reads a file containing a binary configuration (class Binary) and
 * performs various plots. 
 * 
 */

/*
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
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

char lit_bin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2005/09/13 19:47:28  f_limousin
 * Reintroduction of the resolution of the equations in cartesian coordinates.
 *
 * Revision 1.1  2004/09/16 12:14:47  f_limousin
 * *** empty log message ***
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
#include "unites.h"
#include "binary.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "cmp.h"
#include "tenseur.h" 

// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ; 

//******************************************************************************

int main(int argc, char** argv){

  //    Identification of all the subroutines called by the code : 
    
  //     system("ident lit_bin") ; 
    
    if (argc < 2) {
		cout << 
		"lit_bin : the name of a file containing a binary configuration"
		<< endl << " must be given in argument !" << endl ; 
		abort() ; 
    }
    
    char* nomresu = argv[1] ; 
  

    //    char* nomresu = "resu.d" ;
    cout << "Name of the file to be read : " << nomresu << endl ; 

    cout << endl << 
    "Do you want to draw the boundaries of the various domains (y/n) ? [y]"
	 << endl ; 
    char rep ; 
    cin.get(rep) ;

    bool draw_bound = !(rep == 'n') ; 
        
    using namespace Unites ;
    
    FILE* fich = fopen(nomresu, "r") ; 
    if (fich == 0x0) {
    	cout << "Problem in opening the file " << nomresu << " ! " << endl ; 
		perror(" reason") ; 
		abort() ; 
    }

    int mer ; 
    fread(&mer, sizeof(int), 1, fich) ;	// mer
   
    Mg3d mg1(fich) ;
    Map_et mp1(mg1, fich) ; 
    Eos* peos1 = Eos::eos_from_file(fich) ; 

    Mg3d mg2(fich) ;
    Map_et mp2(mg2, fich) ; 
    Eos* peos2 = Eos::eos_from_file(fich) ; 

    Binary star(mp1, *peos1, mp2, *peos2, fich) ; 

    fclose(fich) ; 
    
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

    cout << endl << "Mapping on which star 2 is defined : " << endl ; 
    cout << "================================== " << endl ; 
    cout << star(2).get_mp() << endl ; 

    star.fait_decouple() ;

    for (int i=1; i<=2; i++) {
	(star.set(i)).update_metric(star(3-i), star.get_omega()) ; 
    }

    for (int i=1; i<=2; i++) {
	(star.set(i)).update_metric_der_comp(star(3-i), star.get_omega()) ; 
    }

    for (int i=1; i<=2; i++) {
	(star.set(i)).equation_of_state() ; 
	(star.set(i)).kinematics(star.get_omega(), star.get_x_axe()) ; 
	(star.set(i)).fait_d_psi() ; 
	(star.set(i)).hydro_euler() ; 
    }

    
    ofstream seqfich("resformat.d") ; 
    if ( !seqfich.good() ) {
	cout << "coal : problem with opening the file resformat.d !" << endl ;
	abort() ;
    }
    star.write_global(seqfich) ; 
    seqfich.close() ; 
    
    cout.precision(6) ;
    cout << "mass_bar = " << star(1).mass_b()/msol  << endl ;
    cout << "mass_vol = " << star.mass_adm_vol()/msol  << endl ;
    cout << "mass_adm = " << star.mass_adm()/msol << endl ;
    cout << "d = " << star.separation() << endl ;
    cout << "ray_eq = " << star(1).ray_eq() << endl ;
    cout << "ray_eq_pi = " << star(1).ray_eq_pi() << endl ;
    cout << "R = " << (star(1).ray_eq() + star(1).ray_eq_pi())/2. << endl ;
    cout << "d_milieu = " << star.separation()/2. + (star(1).ray_eq_pi()
						     - star(1).ray_eq())/2. 
	 << endl ;
    cout << "d/R = " << (star.separation() + (star(1).ray_eq_pi()
	 - star(1).ray_eq()))/(star(1).ray_eq() + star(1).ray_eq_pi())
	 << endl ;
    abort() ;
    

    cout << "Binary system read in file : " << endl ;
    cout << star << endl ; 
    star.display_poly(cout) ; //  Reduced quantities for polytropic EOS

    cout << "ADM mass [M_sol] : " << star.mass_adm() / msol  << endl ; 
/*    cout << "Total energy [M_sol c^2] : " 
	 << star.total_ener() / msol << endl ; 
    cout << "Total angular momentum [M_sol c km] : " 
	 << (star.angu_mom())(2) / msol / km << endl ; 

    
    cout << "Relative error in the Hamiltonian constraint : " << endl ; 
    cout << star.ham_constr() << endl ; 
	 
    cout << "Relative error in the momentum constraint : " << endl ; 
    cout << " X component : " << star.mom_constr()(0) << endl ; 
    cout << " Y component : " << star.mom_constr()(1) << endl ; 
    cout << " Z component : " << star.mom_constr()(2) << endl ; 

*/

    Vector beta (star(1).get_beta()) ;
    beta.change_triad(star(1).get_mp().get_bvect_cart()) ;

    Cmp beta_y (beta(2)) ;

    Sym_tensor gtilde (star(1).get_gtilde().cov()) ;
    //   Sym_tensor gtilde
    Sym_tensor flat (star(1).get_flat().cov()) ;
    Sym_tensor hij (gtilde-flat) ;
    hij.change_triad(star(1).get_mp().get_bvect_cart()) ;
    /*
    Sym_tensor hij1 (star(1).get_hij_auto()) ;
    hij1.change_triad(star(1).get_mp().get_bvect_cart()) ;

    Sym_tensor hij2 (star(2).get_hij_auto()) ;
    hij2.change_triad(star(2).get_mp().get_bvect_cart()) ;
    hij2.change_triad(star(1).get_mp().get_bvect_cart()) ;
    
    */
    Sym_tensor hij1(hij) ;
    Sym_tensor gtilde2 (star(2).get_gtilde().cov()) ;
    Sym_tensor flat2 (star(2).get_flat().cov()) ;
    Sym_tensor hij2 (gtilde2-flat2) ;
    hij2.change_triad(star(2).get_mp().get_bvect_cart()) ;
    //    hij2 = hij2 * star(2).get_decouple() ;
    

    Cmp hxx (hij(1,1)) ;
    Cmp hxy (hij(2,1)) ;
    Cmp hxz (hij(3,1)) ;
    Cmp hyy (hij(2,2)) ;
    Cmp hyz (hij(3,2)) ;
    Cmp hzz (hij(3,3)) ;

    Cmp hxx1 (hij1(1,1)) ;
    Cmp hxy1 (hij1(2,1)) ;
    Cmp hxz1 (hij1(3,1)) ;
    Cmp hyy1 (hij1(2,2)) ;
    Cmp hyz1 (hij1(3,2)) ;
    Cmp hzz1 (hij1(3,3)) ;
    Cmp hxx2 (hij2(1,1)) ;
    Cmp hxy2 (hij2(2,1)) ;
    Cmp hxz2 (hij2(3,1)) ;
    Cmp hyy2 (hij2(2,2)) ;
    Cmp hyz2 (hij2(3,2)) ;
    Cmp hzz2 (hij2(3,3)) ;
   
    Cmp logn (star(1).get_logn()) ;
    Cmp psi4 (star(1).get_psi4()) ;


    ofstream fichmetric("metric.d") ;    
    fichmetric.precision(6) ; 

    int n1, n2 ;
    n1 = 100 ;
    n2 = 10000 ;

    double ori ;
    double fact ;
    fact = star.get_omega() * f_unit / M_PI / c_si ;
    ori = star(1).get_mp().get_ori_x() ;
      
    double ii ;
    fichmetric << "# x/lambda   beta^y   hxx    hyy    hzz " << endl ;
    for (int i=1; i<n1; i++){
	ii = i ;
	fichmetric << (- ii/n1 * ori)*10000. * fact  << " " 
		   << beta_y.val_point(- ori + ii/n1 * ori, M_PI/2, 0.)<<" "
 		   << hxx.val_point(- ori + ii /n1 * ori, M_PI/2, 0.) << " " 
		   << hyy.val_point(- ori + ii /n1 * ori, M_PI/2, 0.) << " " 
		   << hzz.val_point(- ori + ii /n1 * ori, M_PI/2, 0.) << endl ; 
	  //		   << logn.val_point(- ori + ii /n1 * ori, M_PI/2, 0.) << " " 
	  //	   << psi4.val_point(- ori + ii /n1 * ori, M_PI/2, 0.) - 1. << endl ;
    }

    for (int i=0; i<=n2; i++){
	ii = i ;
	fichmetric << (-100.* ii/n2*ori - ori)*10000. * fact  << " " 
		   << beta_y.val_point(-100.*ii/n2*ori, M_PI/2, M_PI)<<" "
 		   << hxx.val_point(-100.*ii/n2*ori, M_PI/2, M_PI) << " " 
		   << hyy.val_point(-100.*ii/n2*ori, M_PI/2, M_PI) << " " 
		   << hzz.val_point(-100.*ii/n2*ori, M_PI/2, M_PI) << endl ;
	  //		   << logn.val_point(-100.*ii/n2*ori, M_PI/2, M_PI) << " " 
	  //	   << psi4.val_point(-100.*ii/n2*ori, M_PI/2, M_PI) - 1. << endl ;
    }
    
    fichmetric.close() ; 
    
     arrete() ;

    //==============================================================
    //  Drawings
    //==============================================================

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

    Cmp surf1 (star(1).get_ent()) ; 
    Cmp surf1_ext(mp1) ; 
    surf1_ext = - 0.2 * surf1(0, 0, 0, 0) ; 
    surf1_ext.annule(0, star(1).get_nzet()-1) ; 
    surf1.annule(star(1).get_nzet(), mg1.get_nzone()-1) ; 
    surf1 = surf1 + surf1_ext ;
    surf1 = raccord_c1(surf1, star(1).get_nzet()) ; 

    Cmp surf2 (star(2).get_ent()) ; 
    Cmp surf2_ext(mp2) ; 
    surf2_ext = - 0.2 * surf2(0, 0, 0, 0) ; 
    surf2_ext.annule(0, star(2).get_nzet()-1) ; 
    surf2.annule(star(2).get_nzet(), mg2.get_nzone()-1) ; 
    surf2 = surf2 + surf2_ext ;

    surf2 = raccord_c1(surf2, star(2).get_nzet()) ; 

    char title[80] ;
    char bslash[2] = {92,  '\0'} ;  // 92 is the ASCII code for backslash 

    ofstream fent("enthalpy.d") ; 
	if ( !fent.good() ) {
		cout << "lit_bin : problem with opening the file enthalpy.d !" << endl ;
		abort() ;
	}
   
    fent << "Enthalpy field at the boundary of last inner domain of star 1 : " 
	 << endl ; 
    int lzet =  star(1).get_nzet() - 1 ; 
    for (int k=0; k<mg1.get_np(lzet); k++) {
	fent << "k = " << k << " : " ; 
	for (int j=0; j<mg1.get_nt(lzet); j++) {
	    fent 
		<< "  " << star(1).get_ent().val_grid_point(lzet, k, j, mg1.get_nr(lzet)-1) ;
	}
	fent << endl ; 
    }

   
    fent << endl << "enthalpy field of star 1 : " << endl ;     
    fent << star(1).get_ent() << endl ; 
    fent.close() ; 
    
    Cmp ent1 (star(1).get_ent()) ;
/*
    des_coupe_z(ent1, 0., 1,
		"Enthalpy (z=0)", &surf1, 1.2, draw_bound ) ; 

    des_coupe_y(ent1, 0., 1,
		"Enthalpy (y=0)", &surf1, 1.2, draw_bound ) ; 


    des_profile (ent1, 0., 2* star(1).ray_eq(), 0., 0.,  
	"H", "H (theta=0)" ) ; 

    des_profile (ent1, 0., 2* star(1).ray_eq(), M_PI/2., M_PI/2.,  
	"H", "H (theta=pi/2,  phi=pi/2)" ) ; 

*/

    //==========================================
    // Metric quantities
    //==========================================

    //----------------------------
    // ln(N)
    //----------------------------

    Cmp logn1 (- star(1).get_logn_auto()) ; 
    Cmp logn2 (- star(2).get_logn_auto()) ; 


    cout << "logn xz plane" << endl ;

    des_coupe_bin_y(logn1, logn2, 0, 
			xdes_min, xdes_max, zdes_min, zdes_max, 
		    "ln(N) (y=0)",  &surf1, &surf2, draw_bound ) ; 

    cout << "logn xy plane" << endl ;

    des_coupe_bin_z(logn1, logn2, 0, 
			xdes_min, xdes_max, ydes_min, ydes_max, 
		    "ln(N) (z=0)",  &surf1, &surf2, draw_bound ) ; 
/*
	double xdes_min_large = 2 * xdes_min ; 
	double xdes_max_large = 2 * xdes_max ; 
	double ydes_min_large = 2 * ydes_min ; 
	double ydes_max_large = 2 * ydes_max ; 

    des_coupe_bin_z(logn1, logn2, 0, 
			xdes_min_large, xdes_max_large, ydes_min_large, ydes_max_large, 
		    "ln(N) (z=0)",  &surf1, &surf2, draw_bound ) ; 

    des_coupe_bin_x(logn1, logn2,
			ori_x1, ydes_min, ydes_max, zdes_min, zdes_max, 
		    "ln(N) (x=x1)",  &surf1, 0x0, draw_bound ) ; 
*/
    //--------------
    // Shift vector
    //--------------

    double dmax ;
    dmax = 8. * star(1).ray_eq() ;
    
    Vector beta1 (star(1).get_beta_auto()) ;
    Vector beta2 (star(2).get_beta_auto()) ;
    
    beta1.change_triad(star(1).get_mp().get_bvect_cart()) ;
    beta2.change_triad(star(2).get_mp().get_bvect_cart()) ;
    beta2.change_triad(star(1).get_mp().get_bvect_cart()) ;
    
    Tenseur tmp_beta1 (star(1).get_mp(), 1, CON, star(1).get_mp().get_bvect_cart()) ;
    tmp_beta1.set_etat_qcq() ;
    Cmp tmp_beta11 (beta1(1)) ;
    Cmp tmp_beta12 (beta1(2)) ;
    Cmp tmp_beta13 (beta1(3)) ;
    tmp_beta1.set(0) = tmp_beta11 ;
    tmp_beta1.set(1) = tmp_beta12 ;
    tmp_beta1.set(2) = tmp_beta13 ;

    Tenseur tmp_beta2 (star(2).get_mp(), 1, CON, star(1).get_mp().get_bvect_cart()) ;
    tmp_beta2.set_etat_qcq() ;
    Cmp tmp_beta21 (beta2(1)) ;
    Cmp tmp_beta22 (beta2(2)) ;
    Cmp tmp_beta23 (beta2(3)) ;
    tmp_beta2.set(0) = tmp_beta21 ;
    tmp_beta2.set(1) = tmp_beta22 ;
    tmp_beta2.set(2) = tmp_beta23 ;
    
    cout << "shift_x xy plane" << endl ;
    des_coupe_bin_z(tmp_beta11, tmp_beta21, 0., -dmax, dmax, -dmax, dmax,
		    "shift_x (z=0)" , &surf1, &surf2, draw_bound) ;
    cout << "shift_y xy plane" << endl ;
    des_coupe_bin_z(tmp_beta12, tmp_beta22, 0., -dmax, dmax, -dmax, dmax,
		    "shift_y (z=0)" , &surf1, &surf2, draw_bound) ;
    cout << "shift_y xz plane" << endl ;
    des_coupe_bin_y(tmp_beta12, tmp_beta22, 0., -dmax, dmax, -dmax, dmax,
		    "shift_y (y=0)" , &surf1, &surf2, draw_bound) ;


    cout << "shift xy plane" << endl ;

    des_vect_bin_z(tmp_beta1, tmp_beta2, 0., 
		   -2., 0.5, xdes_min, xdes_max, ydes_min, ydes_max,
		   "Shift vector  (z=0)", 
		   &surf1, &surf2, draw_bound ) ; 
 

    //---------------------------
    // Conformal factor psi4
    //---------------------------
    /*
     Cmp lnq1 (star(1).get_lnq_auto()) ;
     Cmp lnq2 (star(2).get_lnq_auto()) ;
    
     Cmp log_psi1 = 0.5*(lnq1 - logn1) ;
     Cmp log_psi2 = 0.5*(lnq2 - logn2) ;

     cout << "log_psi xz plane" << endl ;
     des_coupe_bin_y(log_psi1, log_psi2, 0, 
		     xdes_min, xdes_max, zdes_min, zdes_max, 
		     "log_psi (y=0)",  &surf1, &surf2, draw_bound ) ; 
     
     cout << "log_psi xy plane" << endl ;
     des_coupe_bin_z(log_psi1, log_psi2, 0, 
		     xdes_min, xdes_max, zdes_min, zdes_max, 
		     "log_psi (z=0)",  &surf1, &surf2, draw_bound ) ; 
    */


    //---------------------------
    // Metric coefficients hij
    //---------------------------

    cout << "hxx-hyy xy plane" << endl ;
    des_coupe_bin_z(hxx1- hyy1, hxx2-hyy2, 0., -dmax, dmax, -dmax, dmax,
		    "hxx-hyy (z=0)" , &surf1, &surf2, draw_bound, 20) ;
    cout << "hxx-hyy xz plane" << endl ;
    des_coupe_bin_y(hxx1- hyy1, hxx2-hyy2, 0., -dmax, dmax, -dmax, dmax,
		    "hxx-hyy (y=0)" , &surf1, &surf2, draw_bound, 20) ;
    cout << "hxx-hzz xy plane" << endl ;
    des_coupe_bin_z(hxx1- hzz1, hxx2-hzz2, 0., -dmax, dmax, -dmax, dmax,
		    "hxx-hzz (z=0)" , &surf1, &surf2, draw_bound, 20) ;
    cout << "hxx-hzz xz plane" << endl ;
    des_coupe_bin_y(hxx1- hzz1, hxx2-hzz2, 0., -dmax, dmax, -dmax, dmax,
		    "hxx-hzz (y=0)" , &surf1, &surf2, draw_bound, 20) ;




    cout << "hxx xy plane" << endl ;
    des_coupe_bin_z(hxx1, hxx2, 0., -dmax, dmax, -dmax, dmax,
		    "hxx (z=0)" , &surf1, &surf2, draw_bound, 20) ;
    cout << "hxy xy plane" << endl ;
    des_coupe_bin_z(hxy1, hxy2, 0., -dmax, dmax, -dmax, dmax,
		"hxy (z=0)", &surf1, &surf2, draw_bound, 20) ;
    cout << "hyy xy plane" << endl ;
    des_coupe_bin_z(hyy1, hyy2, 0., -dmax, dmax, -dmax, dmax,
		"hyy (z=0)", &surf1, &surf2, draw_bound, 20) ;
    cout << "hzz xy plane" << endl ;
    des_coupe_bin_z(hzz1, hzz2, 0., -dmax, dmax, -dmax, dmax,
		"hzz (z=0)", &surf1, &surf2, draw_bound, 20) ;
 
    cout << "hxx xz plane" << endl ;
    des_coupe_bin_y(hxx1, hxx2, 0., -dmax, dmax, -dmax, dmax,
		"hxx (y=0)", &surf1, &surf2, draw_bound, 20) ;
    cout << "hxz xz plane" << endl ;
    des_coupe_bin_y(hxz1, hxz2, 0., -dmax, dmax, -dmax, dmax,
		"hxz (y=0)", &surf1, &surf2, draw_bound, 20) ;
    cout << "hyy xz plane" << endl ;
    des_coupe_bin_y(hyy1, hyy2, 0., -dmax, dmax, -dmax, dmax,
		"hyy (y=0)", &surf1, &surf2, draw_bound, 20) ;
    cout << "hzz xz plane" << endl ;
    des_coupe_bin_y(hzz1, hzz2, 0., -dmax, dmax, -dmax, dmax,
		"hzz (y=0)", &surf1, &surf2, draw_bound, 20) ;
 
    
    /*
    des_profile (hxx, 0, 10, M_PI/2., 0, "Ampli", "hxx on x_axis") ;
    des_profile (hxy, 0, 10, M_PI/2., 0, "Ampli", "hxy on x_axis") ;
    des_profile (hxz, 0, 10, M_PI/2., 0, "Ampli", "hxz on x_axis") ;
    des_profile (hyy, 0, 10, M_PI/2., 0, "Ampli", "hyy on x_axis") ;
    des_profile (hyz, 0, 10, M_PI/2., 0, "Ampli", "hyz on x_axis") ;
    des_profile (hzz, 0, 10, M_PI/2., 0, "Ampli", "hzz on x_axis") ;
    */


    //----------------------------
    // Extrinsic curvature tensor
    //----------------------------
    
    Tensor tkij1 (star(1).get_tkij_auto()) ;
    Tensor tkij2 (star(2).get_tkij_auto()) ;
    tkij1.change_triad(mp1.get_bvect_cart()) ;
    tkij2.change_triad(mp2.get_bvect_cart()) ;
    tkij2.change_triad(mp1.get_bvect_cart()) ;
    
    // Division by r^2 in the external compactified domain in order to get
    //  A^2 K^{ij} : 
    tkij1.dec_dzpuis(2) ;     
    tkij2.dec_dzpuis(2) ;     
    
    char debtit[] = {'K', 92, 'u', '\0'} ; 

    strcpy(title, debtit) ; 
    strcat(title, "xx") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (z=0)") ; 
    
    Cmp tmp11 (tkij1(1, 1)) ;
    Cmp tmp21 (tkij2(1, 1)) ;

    cout << "K^xx (z=0)" << endl ;

    des_coupe_bin_z(tmp11, tmp21, 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 
    
    strcpy(title, debtit) ; 
    strcat(title, "xy") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (z=0)") ; 
    
    Cmp tmp12 (tkij1(1, 2)) ;
    Cmp tmp22 (tkij2(1, 2)) ;

    cout << "K^xy (z=0)" << endl ;

    des_coupe_bin_z(tmp12, tmp22, 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 
    
    strcpy(title, debtit) ; 
    strcat(title, "xz") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (y=0)") ; 
    
    Cmp tmp13 (tkij1(1, 3)) ;
    Cmp tmp23 (tkij2(1, 3)) ;

    des_coupe_bin_x(tmp13, tmp23, 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 
    
    strcpy(title, debtit) ; 
    strcat(title, "yy") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (z=0)") ; 
    

    cout << "K^yy (z=0)" << endl ;

    Cmp tmp14 (tkij1(2, 2)) ;
    Cmp tmp24 (tkij2(2, 2)) ;

    des_coupe_bin_z(tmp14, tmp24, 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 
    
    strcpy(title, debtit) ; 
    strcat(title, "yz") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (y=0)") ; 
    
    Cmp tmp15 (tkij1(2, 3)) ;
    Cmp tmp25 (tkij2(2, 3)) ;

    des_coupe_bin_y(tmp15, tmp25, 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 
    
    strcpy(title, debtit) ; 
    strcat(title, "zz") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (z=0)") ; 
    
    Cmp tmp16 (tkij1(3, 3)) ;
    Cmp tmp26 (tkij2(3, 3)) ;

    des_coupe_bin_z(tmp16, tmp26, 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 

   
    //==========================================
    // Hydro quantities
    //==========================================

    Cmp nbar1 (star(1).get_nbar()) ;
    Cmp nbar2 (star(2).get_nbar()) ;

    cout << "nbar xy plane" << endl ;
    des_coupe_bin_z(nbar1, nbar2, 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    "Baryon density (z=0)",  
		    &surf1, &surf2, draw_bound ) ; 
    
    cout << "nbar xz plane" << endl ;
    des_coupe_bin_y(nbar1, nbar2, 0, 
		    xdes_min, xdes_max, zdes_min, zdes_max, 
		    "Baryon density (y=0)",  
		    &surf1, &surf2, draw_bound ) ; 

    des_coupe_z(nbar1, 0., 1,
		"Baryon density (z=0)", &surf1, 1.2, draw_bound ) ; 

    des_coupe_y(nbar1, 0., 1,
		"Baryon density (y=0)", &surf1, 1.2, draw_bound ) ; 

    if (irrotational) {
	Vector tmp_draw_1 = star(1).get_wit_w() ; 
	tmp_draw_1.annule(star(1).get_nzet(), mg1.get_nzone()-1) ; 
	Vector tmp_draw_2 = star(2).get_wit_w() ; 
	tmp_draw_2.annule(star(2).get_nzet(), mg2.get_nzone()-1) ; 

	tmp_draw_1.change_triad(star(1).get_mp().get_bvect_cart()) ;
	tmp_draw_2.change_triad(star(2).get_mp().get_bvect_cart()) ;
	tmp_draw_2.change_triad(star(1).get_mp().get_bvect_cart()) ;

	Tenseur tmp_draw1 (star(1).get_mp(), 1, CON, star(1).get_mp().get_bvect_cart()) ;
	tmp_draw1.set_etat_qcq() ;
	Cmp tmp_draw11 (tmp_draw_1(1)) ;
	Cmp tmp_draw12 (tmp_draw_1(2)) ;
	Cmp tmp_draw13 (tmp_draw_1(3)) ;
	tmp_draw1.set(0) = tmp_draw11 ;
	tmp_draw1.set(1) = tmp_draw12 ;
	tmp_draw1.set(2) = tmp_draw13 ;

	Tenseur tmp_draw2 (star(2).get_mp(), 1, CON, star(1).get_mp().get_bvect_cart()) ;
	tmp_draw2.set_etat_qcq() ;
	Cmp tmp_draw21 (tmp_draw_2(1)) ;
	Cmp tmp_draw22 (tmp_draw_2(2)) ;
	Cmp tmp_draw23 (tmp_draw_2(3)) ;
	tmp_draw2.set(0) = tmp_draw21 ;
	tmp_draw2.set(1) = tmp_draw22 ;
	tmp_draw2.set(2) = tmp_draw23 ;


	des_vect_bin_z(tmp_draw1, tmp_draw2, 0., 
		    -3., 0.5, xdes_min, xdes_max, ydes_min, ydes_max,
		    "Velocity w.r.t corotating frame  (z=0)", 
		    &surf1, &surf2, draw_bound, 40, 40) ; 
    }
    		    
    Vector tmp_draw_1 = star(1).get_u_euler() ; 
    tmp_draw_1.annule(star(1).get_nzet(), mg1.get_nzone()-1) ; 

    Vector tmp_draw_2 = star(2).get_u_euler() ; 
    tmp_draw_2.annule(star(2).get_nzet(), mg2.get_nzone()-1) ; 

    tmp_draw_1.change_triad(star(1).get_mp().get_bvect_cart()) ;
    tmp_draw_2.change_triad(star(2).get_mp().get_bvect_cart()) ;
    tmp_draw_2.change_triad(star(1).get_mp().get_bvect_cart()) ;
   
    Tenseur tmp_draw1 (star(1).get_mp(), 1, CON, star(1).get_mp().get_bvect_cart()) ;
    tmp_draw1.set_etat_qcq() ;
    Cmp tmp_draw11 (tmp_draw_1(1)) ;
    Cmp tmp_draw12 (tmp_draw_1(2)) ;
    Cmp tmp_draw13 (tmp_draw_1(3)) ;
    tmp_draw1.set(0) = tmp_draw11 ;
    tmp_draw1.set(1) = tmp_draw12 ;
    tmp_draw1.set(2) = tmp_draw13 ;
    
    Tenseur tmp_draw2 (star(2).get_mp(), 1, CON, star(1).get_mp().get_bvect_cart()) ;
    tmp_draw2.set_etat_qcq() ;
    Cmp tmp_draw21 (tmp_draw_2(1)) ;
    Cmp tmp_draw22 (tmp_draw_2(2)) ;
    Cmp tmp_draw23 (tmp_draw_2(3)) ;
    tmp_draw2.set(0) = tmp_draw21 ;
    tmp_draw2.set(1) = tmp_draw22 ;
    tmp_draw2.set(2) = tmp_draw23 ;


    des_coupe_vect_x(tmp_draw1, mp1.get_ori_x(), -1., 0.5, nzdes1, 
		     "U (x=x1)", &surf1, 1.2, draw_bound ) ; 

    cout << "fluid velocity xy plane" << endl ;
    des_vect_bin_z(tmp_draw1, tmp_draw2, 0., 
		    -2., 0.5, xdes_min, xdes_max, ydes_min, ydes_max,
		    "U  (z=0)", &surf1, &surf2, draw_bound, 40, 40) ; 

    if (irrotational) {
	Vector tmp_dpsii = star(1).get_d_psi() ; 
	tmp_dpsii.annule(star(1).get_nzet(), mg1.get_nzone()-1) ; 

	tmp_dpsii.change_triad(star(1).get_mp().get_bvect_cart()) ;
	
	Tenseur tmp_dpsi (star(1).get_mp(), 1, CON, star(1).get_mp().get_bvect_cart()) ;
	tmp_dpsi.set_etat_qcq() ;
	Cmp tmp_dpsi1 (tmp_dpsii(1)) ;
	Cmp tmp_dpsi2 (tmp_dpsii(2)) ;
	Cmp tmp_dpsi3 (tmp_dpsii(3)) ;
	tmp_dpsi.set(0) = tmp_dpsi1 ;
	tmp_dpsi.set(1) = tmp_dpsi2 ;
	tmp_dpsi.set(2) = tmp_dpsi3 ;
	

	des_coupe_vect_z(tmp_dpsi, 0, -1., 0.5, nzdes1,
			 "Grad(psi) (z=0)", &surf1, 1.2, draw_bound ) ;

	Cmp psi00 (star(1).get_psi0()) ;
	des_coupe_z(psi00, 0., 1,
		    "psi0 (z=0)", &surf1, 1.2, draw_bound ) ; 
    
	Tenseur psi000 (psi00) ;
	Tenseur d_psi0 = psi000.gradient() ; 
	//##    d_psi0.change_triad(star.get_ref_triad()) ; 

	des_coupe_vect_z(d_psi0, 0, -3., 0.5, nzdes1, "Grad(psi0) (z=0)", 
			 &surf1, 1.2, draw_bound ) ; 

	Vector tmp_witt = star(1).get_wit_w() ; 
	tmp_witt.annule(star(1).get_nzet(), mg1.get_nzone()-1) ; 

	tmp_witt.change_triad(star(1).get_mp().get_bvect_cart()) ;

	Tenseur tmp_wit (star(1).get_mp(), 1, CON, star(1).get_mp().get_bvect_cart()) ;
	tmp_wit.set_etat_qcq() ;
	Cmp tmp_wit1 (tmp_witt(1)) ;
	Cmp tmp_wit2 (tmp_witt(2)) ;
	Cmp tmp_wit3 (tmp_witt(3)) ;
	tmp_wit.set(0) = tmp_wit1 ;
	tmp_wit.set(1) = tmp_wit2 ;
	tmp_wit.set(2) = tmp_wit3 ;

	des_coupe_vect_z(tmp_wit, 0, -3., 0.5, nzdes1, "W (z=0)", 
			 &surf1, 1.2, draw_bound ) ; 

    }

    // Cleaning
    // --------

    delete peos1 ;    
    delete peos2 ; 

    return EXIT_SUCCESS ;

}
