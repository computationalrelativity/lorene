/*
 * Reads a file containing a binary configuration (class Binary_xcts) and
 * performs various tests and plots. 
 * 
 */
 
/*
 *   Copyright (c) 2010 Michal Bejger 
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

char lit_bin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2011/03/25 10:01:17  m_bejger
 * Curvature tensor plots corrected
 *
 * Revision 1.5  2010/12/09 10:49:37  m_bejger
 * *** empty log message ***
 *
 * Revision 1.4  2010/07/20 19:58:41  m_bejger
 * Correcting diagnostic plots of logn_auto
 *
 * Revision 1.3  2010/06/04 19:51:26  m_bejger
 * Minor corrections
 *
 * Revision 1.2  2010/05/03 13:39:14  m_bejger
 * File handler fixed
 *
 * Revision 1.1  2010/04/29 15:10:47  m_bejger
 * First version: not working properly
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
#include "binary_xcts.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "cmp.h"
#include "tenseur.h" 
#include "nbr_spx.h"

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

    cout << "Name of the file to be read : " << nomresu << endl ; 

    bool draw_bound = true ;      
    using namespace Unites ;
    
    FILE* fich = fopen(nomresu, "r") ; 
    if (fich == 0x0) {
    	cout << "Problem in opening the file " << nomresu << " ! " << endl ; 
		perror(" reason") ; 
		abort() ; 
    }
    
    int mer ; 
    int f_status = fread(&mer, sizeof(int), 1, fich) ;	// mer
                      
    Mg3d mg1(fich) ;
    Map_et mp1(mg1, fich) ; 
    Eos* peos1 = Eos::eos_from_file(fich) ; 

    Mg3d mg2(fich) ;
    Map_et mp2(mg2, fich) ; 
    Eos* peos2 = Eos::eos_from_file(fich) ; 
    
    Binary_xcts star(mp1, *peos1, mp2, *peos2, fich) ; 
    
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

    for (int i=1; i<=2; i++) 
		(star.set(i)).update_metric(star(3-i)) ; 
	         
    for (int i=1; i<=2; i++) 
		(star.set(i)).update_metric_der_comp(star(3-i)) ; 
	 	 
    for (int i=1; i<=2; i++) {
		
		(star.set(i)).equation_of_state() ; 
		(star.set(i)).kinematics(star.get_omega(), star.get_x_axe()) ; 
		(star.set(i)).fait_d_psi() ; 
		(star.set(i)).hydro_euler() ; 
    }
    				
    // Writing of resformat.d
    ofstream seqfich("resformat.d") ; 
    if ( !seqfich.good() ) {
	cout << "coal : problem with opening the file resformat.d !" << endl ;
	abort() ;
    }

    star.write_global(seqfich) ; 
    seqfich.close() ; 

    // Some printings
    cout.precision(6) ;
	cout << "Central gamma determinant :   " 
		 << ((star(1).get_gamma()).determinant()).val_grid_point(0,0,0,0) << endl ;     
	cout << "Central value of N        :   "
		 << (star(1).get_nn()).val_grid_point(0,0,0,0) << endl ;    
	cout << "Central value of Psi^4    :   "
		 << (star(1).get_psi4()).val_grid_point(0,0,0,0) << endl ;    
//	cout << "Central value of Psi_auto    :   "
//		 << (star(1).get_Psi_auto()).val_grid_point(0,0,0,0) << endl ; 
//	cout << "Central value of Psi_comp    :   "
//		 << (star(1).get_Psi_comp()).val_grid_point(0,0,0,0) << endl ; 		          

	cout << "Star(1) mass_bar = " << (star(1).mass_b())/msol  << endl ; 
    cout << "mass_adm_vol     = " << star.mass_adm_vol()/msol  << endl ;
    cout << "mass_kom vol     = " << star.mass_kom_vol()/msol << endl ;
	cout << "mass_adm         = " << star.mass_adm()/msol << endl ;
    cout << "mass_kom         = " << star.mass_kom()/msol << endl ;
    cout << "d = " << star.separation() << endl ;
    cout << "ray_eq = " << star(1).ray_eq() << endl ;
    cout << "ray_eq_pi = " << star(1).ray_eq_pi() << endl ;
    cout << "R = " << (star(1).ray_eq() + star(1).ray_eq_pi())/2. << endl ;
    cout << "d_milieu = " << star.separation()/2. + (star(1).ray_eq_pi()
						     - star(1).ray_eq())/2. << endl ;
    cout << "d/R = " << (star.separation() + (star(1).ray_eq_pi()
	 - star(1).ray_eq()))/(star(1).ray_eq() + star(1).ray_eq_pi())
	 << endl ;

    arrete() ; 
    
    cout << "Binary system read in file : " << endl ;
    cout << star << endl ; 
     
   star.display_poly(cout) ; //  Reduced quantities for polytropic EOS    
   
    //////////////////////////////////////////////////////////
    //    Plot of different fields along X, Y and Z axis    //
    //////////////////////////////////////////////////////////

    Vector beta (star(1).get_beta()) ;
    beta.change_triad(star(1).get_mp().get_bvect_cart()) ;
    Scalar beta_y_aux (beta(2)) ;

    Scalar logn_aux (star(1).get_logn()) ;
    Scalar psi_aux (star(1).get_Psi()) ;
    psi_aux.std_spectral_base() ;
   
    // Construction of an auxiliar grid and mapping
    int nz = star(1).get_mp().get_mg()->get_nzone() ;
    assert (nz >= 4);
    double* bornes = new double [nz+1] ;
    double r_in = 0.95 * (-star(1).get_mp().get_ori_x()-star(1).ray_eq()) ;
    double r_ext = 1.05 * (-star(1).get_mp().get_ori_x()+star(1).ray_eq_pi()) ;

    bornes[0] = 0 ;
    bornes[1] = r_in ;
    bornes[2] = r_ext ;    
    for (int l=3; l<nz; l++)
      bornes[l] = r_ext * pow(2., l-2) ;    
    bornes[nz] = __infinity ;
    
    Map_af mapping (*(star(1).get_mp().get_mg()), bornes) ;
    delete [] bornes ; 
    
    // Importation of fields 
    // ----------------------
    
    assert (star(1).get_mp().get_rot_phi() == 0) ;
    
    Scalar logn (mapping) ;
    logn.import(logn_aux) ;
    logn.std_spectral_base() ;

    Scalar psi (mapping) ;
    psi.import(psi_aux) ;
    psi.std_spectral_base() ;

    Scalar beta_y (mapping) ;
    beta_y.import(beta_y_aux) ;
    beta_y.std_spectral_base() ;

     if (star(1).get_nzet() > 1) {
      cout.precision(10) ;

      for (int ltrans = 0; ltrans < star(1).get_nzet()-1; ltrans++) {
	cout << endl << "Star 1 : values at boundary between domains no. " 
	<< ltrans << " and " << 	ltrans+1 << " for theta = pi/2 and phi = 0 :" << endl ;

	double rt1 = star(1).get_mp().val_r(ltrans, 1., M_PI/2, 0.) ; 
	double rt2 = star(1).get_mp().val_r(ltrans+1, -1., M_PI/2, 0.) ; 
	double diff_rt = (rt2 - rt1)/rt1 ; 
	cout << "   Coord. r [km] (left, right, rel. diff) : " 
	  << rt1 / km << "  " << rt2 / km << "  " << diff_rt << endl ; 
  
	int ntm1 = star(1).get_mp().get_mg()->get_nt(ltrans) - 1; 
	int nrm1 = star(1).get_mp().get_mg()->get_nr(ltrans) - 1 ; 
	double ent1 = star(1).get_ent().val_grid_point(ltrans, 0, ntm1, nrm1) ; 
	double ent2 = star(1).get_ent().val_grid_point(ltrans+1, 0, ntm1, 0) ; 
	double diff_ent = (ent2-ent1)/ent1 ; 
	  cout << "   Enthalpy (left, right, rel. diff) : " 
	  << ent1 << "  " << ent2 << "  " << diff_ent << endl ; 

	double press1 = star(1).get_press().val_grid_point(ltrans, 0, ntm1, nrm1) ; 
	double press2 = star(1).get_press().val_grid_point(ltrans+1, 0, ntm1, 0) ; 
	double diff_press = (press2-press1)/press1 ; 
	  cout << "   Pressure (left, right, rel. diff) : " 
	  << press1 << "  " << press2 << "  " << diff_press << endl ; 

	double nb1 = star(1).get_nbar().val_grid_point(ltrans, 0, ntm1, nrm1) ; 
	double nb2 = star(1).get_nbar().val_grid_point(ltrans+1, 0, ntm1, 0) ; 
	double diff_nb = (nb2-nb1)/nb1 ; 
	  cout << "   Baryon density (left, right, rel. diff) : " 
	  << nb1 << "  " << nb2 << "  " << diff_nb << endl ; 
      }

    double r_max = 1.2 * star(1).ray_eq() ; 
    des_profile(star(1).get_nbar(), 0., r_max, M_PI/2, 0., "n", "Baryon density") ; 
    des_profile(star(1).get_ener(), 0., r_max, M_PI/2, 0., "e", "Energy density") ; 
    des_profile(star(1).get_press(), 0., r_max, M_PI/2, 0., "p", "Pressure") ; 
    des_profile(star(1).get_ent(), 0., r_max, M_PI/2, 0., "H", "Enthalpy") ; 
 
    }


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
    char bslash[2] = {92,'\0'} ;  // 92 is the ASCII code for backslash 

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
	  fent << "  " << star(1).get_ent().val_grid_point(lzet, k, j, mg1.get_nr(lzet)-1) ;
	}
	fent << endl ; 
    }

   
    fent << endl << "enthalpy field of star 1 : " << endl ;     
    fent << star(1).get_ent() << endl ; 
    fent.close() ; 
    
    Scalar ent1  (star(1).get_ent()) ;
    
    Scalar Psi1  (star(1).get_Psi()) ;
    Scalar chi1  (star(1).get_chi()) ;
          
    Scalar Psi2  (star(2).get_Psi()) ;     
    Scalar chi2  (star(2).get_chi()) ;
    
    Scalar Psi1_auto  (star(1).get_Psi_auto()) ;
    Scalar Psi2_auto  (star(2).get_Psi_auto()) ;     

    Scalar chi1_auto  (star(1).get_chi_auto()) ;
    Scalar chi2_auto  (star(2).get_chi_auto()) ; 
    
    des_coupe_z(Psi1_auto, 0., 1,
		"Psi1_auto (z=0)", &surf1, 5., draw_bound ) ; 
 
//  des_profile (Psi1_auto, 0., 15* star(1).ray_eq(),  M_PI/2., 0.,  
//	"Psi1_auto", "Psi1_auto (theta=pi/2,  phi=0)" ) ; 
 
//    des_profile (Psi1_comp, 0., 15* star(1).ray_eq(), M_PI/2., 0.,  
//	"Psi1_comp", "Psi1_comp (theta=pi/2,  phi=0)" ) ;  

    des_profile (Psi1, 0., 15* star(1).ray_eq(), M_PI/2., 0.,  
	"Psi1", "Psi1 (theta=pi/2,  phi=0)" ) ;  

/*
//    des_profile (chi1_auto, 0., 15* star(1).ray_eq(),  M_PI/2., 0.,  
//	"chi1_auto", "chi1_auto (theta=pi/2,  phi=0)" ) ; 
 
//    des_profile (chi1_comp, 0., 15* star(1).ray_eq(), M_PI/2., 0.,  
//	"chi1_comp", "chi1_comp (theta=pi/2,  phi=0)" ) ;  

    des_profile (chi1, 0., 15* star(1).ray_eq(), M_PI/2., 0.,  
	"chi1", "chi1 (theta=pi/2,  phi=0" ) ;  
	
    des_coupe_z(ent1, 0., 1,
		"Enthalpy (z=0)", &surf1, 1.2, draw_bound ) ; 

    des_coupe_y(ent1, 0., 1,
		"Enthalpy (y=0)", &surf1, 1.2, draw_bound ) ; 
	
    des_profile (ent1, 0., 2* star(1).ray_eq(), 0., 0.,  
	"H", "H (theta=0)" ) ; 

    des_profile (ent1, 0., 2* star(1).ray_eq(), M_PI/2., M_PI/2.,  
	"H", "H (theta=pi/2,  phi=pi/2)" ) ; 



    //==========================================
    // Metric quantities
    //==========================================

    //----------------------------
    // ln(N)
    //----------------------------

    Scalar logn1 (log((chi1_auto + 1.)/(Psi1_auto + 1.))) ;
    Scalar logn2 (log((chi2_auto + 1.)/(Psi2_auto + 1.))) ;     

	logn1.std_spectral_base() ; 
	logn2.std_spectral_base() ; 

    cout << "logn xz plane" << endl ;

    des_coupe_bin_y(logn1, logn2, 0, 
			xdes_min, xdes_max, zdes_min, zdes_max, 
		    "ln(N) (y=0)",  &surf1, &surf2, draw_bound ) ; 

    cout << "logn xy plane" << endl ;

    des_coupe_bin_z(logn1, logn2, 0, 
			xdes_min, xdes_max, ydes_min, ydes_max, 
		    "ln(N) (z=0)",  &surf1, &surf2, draw_bound ) ; 

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
*/     
    //----------------------------
    // Extrinsic curvature tensor
    //----------------------------
    
    Tensor tkij1 (star(1).get_haij_auto()) ;
    Tensor tkij2 (star(2).get_haij_auto()) ;
    
    // Division by r^2 in the external compactified domain 
    // in order to get \hat{A}^{ij} : 
    tkij1.dec_dzpuis(2) ;     
    tkij2.dec_dzpuis(2) ; 
              
    char debtit[] = {'A', 92, 'u', '\0'} ; 

    //----------------------------------- 
	cout << "\\hat{A}^xx (z=0)" << endl ;
    strcpy(title, debtit) ; 
    strcat(title, "xx") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (z=0)") ; 
    
    des_coupe_bin_z(tkij1(1, 1), tkij2(1, 1), 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title, &surf1, &surf2, draw_bound ) ; 

	//-----------------------------------
	cout << "\\hat{A}^xy (z=0)" << endl ;
    strcpy(title, debtit) ; 
    strcat(title, "xy") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (z=0)") ; 
    
    des_coupe_bin_z(tkij1(1, 2), tkij2(1, 2), 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 
    
	//-----------------------------------
	cout << "\\hat{A}^xz (y=0)" << endl ;
    strcpy(title, debtit) ; 
    strcat(title, "xz") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (y=0)") ; 
    
    des_coupe_bin_y(tkij1(1, 3), tkij2(1, 3), 0, 
		    xdes_min, xdes_max, zdes_min, zdes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 

	//-----------------------------------    
    cout << "\\hat{A}^yy (z=0)" << endl ;
    strcpy(title, debtit) ; 
    strcat(title, "yy") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (z=0)") ; 

    des_coupe_bin_z(tkij1(2, 2), tkij2(2, 2), 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 

	//-----------------------------------    
    cout << "\\hat{A}^yz (y=0)" << endl ;    
    strcpy(title, debtit) ; 
    strcat(title, "yz") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (y=0)") ; 

    des_coupe_bin_y(tkij1(2, 3), tkij2(2, 3), 0, 
		    xdes_min, xdes_max, zdes_min, zdes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 

	//-----------------------------------    
    cout << "\\hat{A}^zz (z=0)" << endl ;    
    strcpy(title, debtit) ; 
    strcat(title, "zz") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (z=0)") ; 

    des_coupe_bin_z(tkij1(3, 3), tkij2(3, 3), 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 
		    
	//----------------------------
	// \hat{A}_ij \hat{A}^ij term
 	//----------------------------
    
	Tenseur hacar1 ( star(1).get_hacar_auto() ) ; 
	Tenseur hacar2 ( star(2).get_hacar_auto() ) ; 

	// Division by r^4 in the external compactified domain 
	// in order to get \hat{A}_ij \hat{A}^{ij} : 
	hacar1.dec2_dzpuis() ; 
	hacar1.dec2_dzpuis() ; 
	hacar2.dec2_dzpuis() ; 
	hacar2.dec2_dzpuis() ; 
    
	char title1[80] ;
	strcpy(title1, debtit) ; 
	strcat(title1, "ij") ; 
	strcat(title1, bslash) ; 
	strcat(title1, "d A") ; 
	strcat(title1, bslash) ; 
	strcat(title1, "dij") ; 
	strcat(title1, bslash) ; 

	strcpy(title, title1) ; 
	strcat(title, "u (x=x1)") ; 
	
	des_coupe_bin_x(hacar1(), hacar2(), ori_x1, 
			ydes_min, ydes_max, zdes_min, zdes_max, 
			title,  &surf1, 0x0, draw_bound ) ; 

	strcpy(title, title1) ; 
	strcat(title, "u (y=0)") ; 
	
	des_coupe_bin_y(hacar1(), hacar2(), 0, 
			xdes_min, xdes_max, zdes_min, zdes_max, 
			title,  &surf1, &surf2, draw_bound ) ; 

	strcpy(title, title1) ; 
	strcat(title, "u (z=0)") ; 
	
	des_coupe_bin_z(hacar1(), hacar2(), 0, 
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
    
    //cout << "nbar xz plane" << endl ;
    //des_coupe_bin_y(nbar1, nbar2, 0, 
	//	    xdes_min, xdes_max, zdes_min, zdes_max, 
	//	    "Baryon density (y=0)",  
	//	    &surf1, &surf2, draw_bound ) ; 

    //des_coupe_z(nbar1, 0., 1,
	//	"Baryon density (z=0)", &surf1, 1.2, draw_bound ) ; 

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
