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
 * Revision 1.12  2014/10/13 08:53:54  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.11  2014/10/06 15:09:43  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.10  2008/11/14 13:55:44  e_gourgoulhon
 * Added more outputs in the case of more than one domain inside the
 * stars.
 *
 * Revision 1.9  2004/03/25 12:35:37  j_novak
 * now using namespace Unites
 *
 * Revision 1.8  2003/09/18 07:31:49  e_gourgoulhon
 * -- Test of good opening of files
 * -- Suppressed the call to des_explorer_symz
 * -- Enlarged view of log(N) in the z=0 plane
 *
 * Revision 1.7  2003/09/16 09:17:15  e_gourgoulhon
 * Changed the name of the output file "seq.d" to "resformat0.d".
 *
 * Revision 1.6  2003/09/15 15:11:25  e_gourgoulhon
 * Added the creation of file "seq.d" by a call to Binaire::write_global.
 *
 * Revision 1.5  2003/09/08 08:21:43  e_gourgoulhon
 * Added printing out of the virial errors in the relativistic case.
 *
 * Revision 1.4  2003/09/07 22:11:24  e_gourgoulhon
 * Corrected a bug at the end (introduced new vectors tmp_dpsi and
 * tmp_wit for drawings in the irrotational case).
 *
 * Revision 1.3  2003/01/09 11:07:49  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.2  2002/06/10 10:06:36  f_limousin
 * "void main" changed to "int main".
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 * Revision 1.17  2001/04/17  11:34:26  eric
 * Modif determination x,y,zdes_min,max (adaptation au cas de
 *  deux etoiles differentes).
 *
 * Revision 1.16  2001/03/07  10:58:40  eric
 * Appel de des_explorer_symz (pour dessiner des etoiles completes)
 * Annulation des champs de vitesse dans les zones externes.
 *
 * Revision 1.15  2000/12/20  10:27:21  eric
 * Ajout de la sortie dans le fichier "enthalpy.d" de l'enthalpie pour
 *  un examen des valeurs point par point.
 *
 * Revision 1.14  2000/07/07  14:17:30  eric
 * Appel de la nouvelle fonction Binaire::display_poly pour l'affichage
 *   des quantites en unites polytropiques.
 *
 * Revision 1.13  2000/07/06  10:23:56  eric
 * Ajout affichage barycentre dans les quantites polytropiques.
 *
 * Revision 1.12  2000/06/22  16:13:14  eric
 * Ameliorations diverses.
 *
 * Revision 1.11  2000/06/21  16:09:39  keisuke
 * Suppress the lines for the boundaries of second domains.
 *
 * Revision 1.10  2000/06/11  15:52:15  keisuke
 * Add a few figures for the baryon density.
 *
 * Revision 1.9  2000/03/22  10:32:18  eric
 * Appel de raccord_c1 sur surf1 et surf2 avant le dessin.
 *
 * Revision 1.8  2000/03/15  13:16:07  eric
 * Changement de l'ordre des dessins.
 *
 * Revision 1.7  2000/03/13  17:05:51  eric
 * Ajout de la verification des equations de constrainte.
 *
 * Revision 1.6  2000/03/08  14:46:18  eric
 * Dessin des quantites relativistes.
 *
 * Revision 1.5  2000/03/02  10:32:04  eric
 * Affichage des mappings.
 * Introduction des champs surf1 et surf2 pour definir les surfaces des
 *  etoiles (plutot que de prendre directement l'enthalpie)
 * Dessins de champs de vecteurs.
 *
 * Revision 1.4  2000/02/18  15:50:06  eric
 * Affichage des quantites globales a la fin.
 *
 * Revision 1.3  2000/02/17  19:58:14  eric
 * Modif affichage
 *
 * Revision 1.2  2000/02/12  12:01:30  eric
 * *** empty log message ***
 *
 * Revision 1.1  2000/02/12  11:53:15  eric
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
#include "binaire.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "unites.h"	    

namespace Lorene {
// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ; 
}
//******************************************************************************

using namespace Lorene ;

int main(int argc, char** argv){

  using namespace Unites ;
    // Identification of all the subroutines called by the code : 
    
    // system("ident lit_bin") ; 
    
    if (argc < 2) {
		cout << 
		"lit_bin : the name of a file containing a binary configuration"
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
    
    Binaire star(mp1, *peos1, mp2, *peos2, fich) ; 

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
    star.display_poly(cout) ; //  Reduced quantities for polytropic EOS

    cout << "ADM mass [M_sol] : " << star.mass_adm() / msol  << endl ; 
    cout << "Total energy [M_sol c^2] : " 
	 << star.total_ener() / msol << endl ; 
    cout << "Total angular momentum [M_sol c km] : " 
	 << (star.angu_mom())(2) / msol / km << endl ; 

    if ( relativistic ) {
	cout << "Relative error on the virial theorem : " << endl ; 
	cout << "   VE(M)= " << star.virial() 
	          << "   VE(GB)= "<< star.virial_gb()  
		  << "   VE(FUS)= " << star.virial_fus() << endl ;   
    }
    else {
	cout << "Relative error on the virial theorem : " 
		  <<  star.virial() << endl ; 
    }
    
    cout << "Relative error in the Hamiltonian constraint : " << endl ; 
    cout << star.ham_constr() << endl ; 
	 
    cout << "Relative error in the momentum constraint : " << endl ; 
    cout << " X component : " << star.mom_constr()(0) << endl ; 
    cout << " Y component : " << star.mom_constr()(1) << endl ; 
    cout << " Z component : " << star.mom_constr()(2) << endl ; 


	ofstream seqfich("resformat0.d") ; 
	if ( !seqfich.good() ) {
		cout << "lit_bin : problem with opening the file resformat0.d !" << endl ;
		abort() ;
	}
	star.write_global(seqfich) ; 
	seqfich.close() ; 
	

     if (star(1).get_nzet() > 1) {
      cout.precision(10) ;

      for (int ltrans = 0; ltrans < star(1).get_nzet()-1; ltrans++) {
	cout << endl << "Star 1 : values at boundary between domains no. " << ltrans << " and " << 	ltrans+1 << " for theta = pi/2 and phi = 0 :" << endl ;

	double rt1 = star(1).get_mp().val_r(ltrans, 1., M_PI/2, 0.) ; 
	double rt2 = star(1).get_mp().val_r(ltrans+1, -1., M_PI/2, 0.) ; 
	double diff_rt = (rt2 - rt1)/rt1 ; 
	cout << "   Coord. r [km] (left, right, rel. diff) : " 
	  << rt1 / km << "  " << rt2 / km << "  " << diff_rt << endl ; 
  
	int ntm1 = star(1).get_mp().get_mg()->get_nt(ltrans) - 1; 
	int nrm1 = star(1).get_mp().get_mg()->get_nr(ltrans) - 1 ; 
	double ent1 = star(1).get_ent()()(ltrans, 0, ntm1, nrm1) ; 
	double ent2 = star(1).get_ent()()(ltrans+1, 0, ntm1, 0) ; 
	double diff_ent = (ent2-ent1)/ent1 ; 
	  cout << "   Enthalpy (left, right, rel. diff) : " 
	  << ent1 << "  " << ent2 << "  " << diff_ent << endl ; 

	double press1 = star(1).get_press()()(ltrans, 0, ntm1, nrm1) ; 
	double press2 = star(1).get_press()()(ltrans+1, 0, ntm1, 0) ; 
	double diff_press = (press2-press1)/press1 ; 
	  cout << "   Pressure (left, right, rel. diff) : " 
	  << press1 << "  " << press2 << "  " << diff_press << endl ; 

	double nb1 = star(1).get_nbar()()(ltrans, 0, ntm1, nrm1) ; 
	double nb2 = star(1).get_nbar()()(ltrans+1, 0, ntm1, 0) ; 
	double diff_nb = (nb2-nb1)/nb1 ; 
	  cout << "   Baryon density (left, right, rel. diff) : " 
	  << nb1 << "  " << nb2 << "  " << diff_nb << endl ; 
      }

    double r_max = 1.2 * star(1).ray_eq() ; 
    des_profile(star(1).get_nbar()(), 0., r_max, M_PI/2, 0., "n", "Baryon density") ; 
    des_profile(star(1).get_ener()(), 0., r_max, M_PI/2, 0., "e", "Energy density") ; 
    des_profile(star(1).get_press()(), 0., r_max, M_PI/2, 0., "p", "Pressure") ; 
    // des_profile(star(1).get_ent()(), 0., r_max, M_PI/2, 0., "H", "Enthalpy") ; 
 
    }


    //==============================================================
    //  Drawings
    //==============================================================



    // des_explorer_symz(star, "latbin") ; 
    
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

	double xdes_min_large = 2 * xdes_min ; 
	double xdes_max_large = 2 * xdes_max ; 
	double ydes_min_large = 2 * ydes_min ; 
	double ydes_max_large = 2 * ydes_max ; 

    des_coupe_bin_z(logn1(), logn2(), 0, 
			xdes_min_large, xdes_max_large, ydes_min_large, ydes_max_large, 
		    "ln(N) (z=0)",  &surf1, &surf2, draw_bound ) ; 

    des_coupe_bin_x(logn1(), logn2(),
			ori_x1, ydes_min, ydes_max, zdes_min, zdes_max, 
		    "ln(N) (x=x1)",  &surf1, 0x0, draw_bound ) ; 

   
    if (relativistic) {

	//--------------
	// Shift vector
	//--------------
	des_vect_bin_z(star(1).get_shift_auto(), star(2).get_shift_auto(), 0., 
		    -2., 0.5, xdes_min, xdes_max, ydes_min, ydes_max,
		    "Shift vector  (z=0)", 
		    &surf1, &surf2, draw_bound ) ; 

	//----------------------------
	// Extrinsic curvature tensor
	//----------------------------
    
	Tenseur tkij1 = star(1).get_tkij_auto() ;
	Tenseur tkij2 = star(2).get_tkij_auto() ;
    
	// Division by r^2 in the external compactified domain in order to get
	//  A^2 K^{ij} : 
	tkij1.dec2_dzpuis() ;     
	tkij2.dec2_dzpuis() ;     
    
	char debtit[] = {'A', 92, 'u', '2', 92, 'd', ' ', 'K', 92, 'u', '\0'} ; 

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
	// Term A^2 K_ij K^ij
	//----------------------------
    
	Tenseur akcar1 = star(1).get_akcar_auto() ; 
	Tenseur akcar2 = star(2).get_akcar_auto() ; 

	// Division by r^4 in the external compactified domain in order to get
	//  A^2 K_ij K^{ij} : 
	akcar1.dec2_dzpuis() ; 
	akcar1.dec2_dzpuis() ; 
	akcar2.dec2_dzpuis() ; 
	akcar2.dec2_dzpuis() ; 
    
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
	
	des_coupe_bin_x(akcar1(), akcar2(), ori_x1, 
			ydes_min, ydes_max, zdes_min, zdes_max, 
			title,  &surf1, 0x0, draw_bound ) ; 

	strcpy(title, title1) ; 
	strcat(title, "u (y=0)") ; 
	
	des_coupe_bin_y(akcar1(), akcar2(), 0, 
			xdes_min, xdes_max, zdes_min, zdes_max, 
			title,  &surf1, &surf2, draw_bound ) ; 

	strcpy(title, title1) ; 
	strcat(title, "u (z=0)") ; 
	
	des_coupe_bin_z(akcar1(), akcar2(), 0, 
			xdes_min, xdes_max, ydes_min, ydes_max, 
			title,  &surf1, &surf2, draw_bound ) ; 
    
	//----------------------------
	// beta = ln(AN)
	//----------------------------
    
	Tenseur beta1 = - star(1).get_beta_auto() ; 
	Tenseur beta2 = - star(2).get_beta_auto() ; 
	
	des_coupe_bin_x(beta1(), beta2(),
			ori_x1, ydes_min, ydes_max, zdes_min, zdes_max, 
		    "ln(AN) (x=x1)",  &surf1, 0x0, draw_bound ) ; 

	des_coupe_bin_y(beta1(), beta2(), 0, 
			xdes_min, xdes_max, zdes_min, zdes_max, 
		    "ln(AN) (y=0)",  &surf1, &surf2, draw_bound ) ; 

	des_coupe_bin_z(beta1(), beta2(), 0, 
			xdes_min, xdes_max, ydes_min, ydes_max, 
		    "ln(AN) (z=0)",  &surf1, &surf2, draw_bound ) ; 

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

    des_coupe_vect_x(tmp_draw1, mp1.get_ori_x(), -1., 0.5, nzdes1, 
		     "U (x=x1)", &surf1, 1.2, draw_bound ) ; 

    des_vect_bin_z(tmp_draw1, tmp_draw2, 0., 
		    -2., 0.5, xdes_min, xdes_max, ydes_min, ydes_max,
		    "U  (z=0)", &surf1, &surf2, draw_bound, 40, 40) ; 

    if (irrotational) {
	Tenseur tmp_dpsi = star(1).get_d_psi() ; 
	tmp_dpsi.annule(star(1).get_nzet(), mg1.get_nzone()-1) ; 
	des_coupe_vect_z(tmp_dpsi, 0, -1., 0.5, nzdes1,
			 "Grad(psi) (z=0)", &surf1, 1.2, draw_bound ) ;

	des_coupe_z(star(1).get_psi0()(), 0., 1,
		    "psi0 (z=0)", &surf1, 1.2, draw_bound ) ; 
    
	Tenseur d_psi0 = star(1).get_psi0().gradient() ; 
	//##    d_psi0.change_triad(star.get_ref_triad()) ; 

	des_coupe_vect_z(d_psi0, 0, -3., 0.5, nzdes1, "Grad(psi0) (z=0)", 
			 &surf1, 1.2, draw_bound ) ; 

	Tenseur tmp_wit = star(1).get_wit_w() ; 
	tmp_wit.annule(star(1).get_nzet(), mg1.get_nzone()-1) ; 
	des_coupe_vect_z(tmp_wit, 0, -3., 0.5, nzdes1, "W (z=0)", 
			 &surf1, 1.2, draw_bound ) ; 

    }

    // Cleaning
    // --------

    delete peos1 ;    
    delete peos2 ; 

    return EXIT_SUCCESS ;

}
