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
  
    //char* nomresu = "resu.d" ;
    cout << "Name of the file to be read : " << nomresu << endl ; 

   /* cout << endl << 
    "Do you want to draw the boundaries of the various domains (y/n) ? [y]"
	 << endl ; 
    char rep ; 
    cin.get(rep) ;
   
    bool draw_bound = !(rep == 'n') ; 
        
    */
    
    bool draw_bound = true ;      
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

    //star.fait_decouple() ;
            
    for (int i=1; i<=2; i++) {
		(star.set(i)).update_metric(star(3-i)) ; 
		             
        Scalar Psiauto (star(i).get_Psi_auto() ) ; 
        Psiauto.std_spectral_base() ; 
        Scalar Psicomp ( star(i).get_Psi_comp() ) ; 
        Psicomp.std_spectral_base() ;         

        Scalar chiauto (star(i).get_chi_auto() ) ; 
        chiauto.std_spectral_base() ; 
        Scalar chicomp ( star(i).get_chi_comp() ) ; 
        chicomp.std_spectral_base() ; 
                       
        Scalar psi_total = Psiauto + Psicomp + 1.; 
        psi_total.std_spectral_base() ; 
              
        Scalar psi4_test = pow(psi_total, 4.) ; 
        psi4_test.std_spectral_base() ; 
        
        Scalar logn_test ( star(i).get_chi()/star(i).get_Psi() ) ; 
        logn_test.std_spectral_base() ; 
        
        Scalar nvalue_total = log(chiauto + chicomp + 1.)/(Psiauto + Psicomp + 1.) ; 
        nvalue_total.std_spectral_base() ; 
              
	    cout << "For star(" << i << ") : " << endl ;  
        cout << "Central gamma determinant :   " << ((star(i).get_gamma()).determinant()).val_grid_point(0,0,0,0) << endl ;     
        cout << "Central value of log(N)   :   " << logn_test.val_grid_point(0,0,0,0) << endl ;
        cout << "Central value of log(N)   :   " << nvalue_total.val_grid_point(0,0,0,0) << endl ;
        cout << "get_psi4                      " << psi4_test.val_grid_point(0,0,0,0) << endl ;
        cout << "psi_total                     " << psi_total.val_grid_point(0,0,0,0) << endl ;        
	    cout << "Central value of Psiauto  :   " << Psiauto.val_grid_point(0,0,0,0) <<  endl ; 
	    cout << "Central value of Psicomp  :   " << Psicomp.val_grid_point(0,0,0,0) <<  endl ; 
	    cout << "Central value of chiauto  :   " << chiauto.val_grid_point(0,0,0,0) <<  endl ; 
	    cout << "Central value of chicomp  :   " << chicomp.val_grid_point(0,0,0,0) <<  endl ; 
    }
         
    arrete() ; 
         
    for (int i=1; i<=2; i++) {
		(star.set(i)).update_metric_der_comp(star(3-i)) ; 
    }
	 	 
    for (int i=1; i<=2; i++) {
		(star.set(i)).equation_of_state() ; 
		
		cout << "get omega : " << star.get_omega() << endl ;
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
    cout << "Star(1) mass_bar = " << (star(1).mass_b())/msol  << endl ; 
    cout << "mass_adm_vol     = " << star.mass_adm_vol()/msol  << endl ;
    cout << "mass_adm         = " << star.mass_adm()/msol << endl ;
    cout << "mass_kom vol     = " << star.mass_kom_vol()/msol << endl ;
    cout << "mass_kom         = " << star.mass_kom()/msol << endl ;
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

    Scalar Psi1_auto  (star(1).get_Psi_auto()) ;
    Scalar Psi2_auto  (star(2).get_Psi_auto()) ;     

    Scalar chi1_auto  (star(1).get_chi_auto()) ;
    Scalar chi2_auto  (star(2).get_chi_auto()) ; 
    
    des_coupe_z(Psi1_auto, 0., 1,
		"Psi1_auto (z=0)", &surf1, 5., draw_bound ) ; 
 
/*  des_profile (Psi1_auto, 0., 15* star(1).ray_eq(),  M_PI/2., 0.,  
	"Psi1_auto", "Psi1_auto (theta=pi/2,  phi=0)" ) ; 
 
    des_profile (Psi1_comp, 0., 15* star(1).ray_eq(), M_PI/2., 0.,  
	"Psi1_comp", "Psi1_comp (theta=pi/2,  phi=0)" ) ;  
*/
    des_profile (Psi1, 0., 15* star(1).ray_eq(), M_PI/2., 0.,  
	"Psi1", "Psi1 (theta=pi/2,  phi=0)" ) ;  

/*    des_profile (chi1_auto, 0., 15* star(1).ray_eq(),  M_PI/2., 0.,  
	"chi1_auto", "chi1_auto (theta=pi/2,  phi=0)" ) ; 
 
    des_profile (chi1_comp, 0., 15* star(1).ray_eq(), M_PI/2., 0.,  
	"chi1_comp", "chi1_comp (theta=pi/2,  phi=0)" ) ;  
*/
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
     
    //----------------------------
    // Extrinsic curvature tensor
    //----------------------------
    
    Tensor tkij1 (star(1).get_haij_auto()) ;
    Tensor tkij2 (star(2).get_haij_auto()) ;
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

    cout << "\\hat{A}^xx (z=0)" << endl ;

    des_coupe_bin_z(tmp11, tmp21, 0, 
		    xdes_min, xdes_max, ydes_min, ydes_max, 
		    title,  &surf1, &surf2, draw_bound ) ; 
    
    strcpy(title, debtit) ; 
    strcat(title, "xy") ; 
    strcat(title, bslash) ; 
    strcat(title, "d (z=0)") ; 
    
    Cmp tmp12 (tkij1(1, 2)) ;
    Cmp tmp22 (tkij2(1, 2)) ;

    cout << "\\hat{A}^xy (z=0)" << endl ;

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
    

    cout << "\\hat{A}^yy (z=0)" << endl ;

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
