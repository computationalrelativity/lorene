/*
 * Code for reading a binary black hole or binary neutron star configuration
 *  and performing some analysis of the asymptotic
 *  behavior of the metric coefficients
 *
 */

/*
 *   Copyright (c) 2002 Francois Limousin
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
char analyse_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/06 15:09:42  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2004/03/25 12:35:36  j_novak
 * now using namespace Unites
 *
 * Revision 1.2  2003/01/09 11:07:48  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1  2002/06/18 14:07:26  f_limousin
 * Analysis of asymptotic behavior of binary NS and BH
 *
 *
 *
 * $Header$
 *
 */

// Headers standard du C
//  (par exemple definit la macro EXIT_SUCCESS)
#include <cstdlib>
#include <cstdio>

 
// Headers Lorene
#include "tenseur.h"
#include "bhole.h"
#include "binaire.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"

// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ; 

void asymptot(const Cmp&,  const char*, bool, ostream& );
 
int main() {

  using namespace Unites ;

  // Reads the parameters of the computation
  // ---------------------------------------

  char nom_config[100], blabla[100] ; 
  int type_objet, nr_a, nt_a, np_a, graphics_i ;
  double rayon_a ; 

  ifstream fichparam("par_ana.d") ;
  fichparam >> graphics_i ; fichparam.getline(blabla,100) ;
  fichparam >> type_objet ; fichparam.getline(blabla,100) ; 
  fichparam.getline(nom_config,100) ; 
  fichparam >> nr_a ; fichparam.getline(blabla,100) ; 
  fichparam >> nt_a ; fichparam.getline(blabla,100) ; 
  fichparam >> np_a ; fichparam.getline(blabla,100) ; 
  fichparam >> rayon_a ; fichparam.getline(blabla,100) ;  
  fichparam.close() ; 

  bool black_hole = (type_objet == 1) ;
  bool graphics = (graphics_i == 1) ;

  cout << "Name of the file containing the binary system :"
       << endl << nom_config << endl ;

  // Pointers on global objects (independent of BH or NS)
  Mg3d* grille_un ;
  Mg3d* grille_deux ;
  Map* map_un ;
  Map* map_deux ; 
  Cmp* p_nn1 ;
  Cmp* p_nn2 ; 
  Tenseur* p_shift1 ;
  Tenseur* p_shift2 ;

  FILE* fich = fopen(nom_config, "r") ;
  if (fich == 0x0) {
      cout << "Problem in opening the file " << nom_config << " !" << endl ;
      abort() ; 
  }

  //*******************************************************************
  // Reading of binary black hole data
  //*******************************************************************
  if (black_hole) {

    grille_un = new Mg3d(fich) ;
    grille_deux = grille_un ; 
    Map_af* map_un_af = new Map_af(*grille_un, fich) ;
    Map_af* map_deux_af = new Map_af(*grille_deux, fich) ;
    map_un = map_un_af ; 
    map_deux = map_deux_af ; 
    Bhole hole_un (*map_un_af, fich) ;
    Bhole hole_deux (*map_deux_af, fich) ;
    fclose(fich) ;
    
    assert (hole_un.get_omega() == hole_deux.get_omega()) ;

    cout << "Multi-grid read in file : " << endl << *grille_un << endl ;
    // arrete() ; 
    cout << "Mapping 1 read in file : " << endl << *map_un << endl ;
    cout << "Mapping 2 read in file : " << endl << *map_deux << endl ;
    // arrete() ; 

    // Construction of the binary system
    // ---------------------------------
    Bhole_binaire systeme (*map_un_af, *map_deux_af) ;
    systeme.set(1) = hole_un ;
    systeme.set(2) = hole_deux ;
    systeme.set_omega(hole_un.get_omega()) ;

    // On initialise les grandeurs derivees :
    systeme.set(1).fait_n_comp (systeme(2)) ;
    systeme.set(1).fait_psi_comp (systeme(2)) ;
    systeme.set(2).fait_n_comp (systeme(1)) ;
    systeme.set(2).fait_psi_comp (systeme(1)) ;
    systeme.fait_decouple() ;
    systeme.fait_tkij() ;

    // Initialisation of member data
    // -----------------------------

    // Unit of length:
    double aa = systeme(1).get_rayon() ;
    double aa2 = systeme(2).get_rayon() ;

    double omega = systeme.get_omega() * aa ;
    double dist = ( map_un->get_ori_x() - map_deux->get_ori_x() ) / aa ;

    cout << endl << "Binary system read in file : " << endl ;
    cout <<	    "---------------------------- " << endl ;
    cout << "  Separation d/a :       " << dist << endl ;
    cout << "  Omega :                " << omega << " / a" << endl ;
    cout << "  Size of black hole 2 : " << aa2 / aa << " a" << endl ;
    cout << "  ADM mass :             " << systeme.adm_systeme() / aa
         << " a" << endl ;
    cout << "  Komar-lile mass :      " << systeme.komar_systeme() / aa
         << " a" << endl ;
    cout << "  Angular momentum :     " << systeme.moment_systeme_inf() 
	    / (aa*aa) << " a^2" << endl ; 
    cout << "  Proper distance between the two throats : "
	  << systeme.distance_propre() / aa << " a" << endl ;
    cout << "  Area of black hole 1 apparent horizon : " << 
	      systeme(1).area() / (aa*aa) << " a^2" << endl ; 
    cout << "  Area of black hole 2 apparent horizon : " <<
	      systeme(2).area() / (aa*aa) << " a^2" << endl ; 

    //------------------------------------------------------
    //  Lapse function
    //------------------------------------------------------

    // nn1 defined as a Tenseur copy of systeme(1).n_auto
    //  (constructed by the copy constructor of class Tenseur) 
    // Tenseur nn1 = systeme(1).get_n_auto() ; 

    

    // *nn1 is a copy of the Cmp systeme(1).n_auto.c[0]
    p_nn1 = new Cmp( systeme(1).get_n_auto()() ) ;

    // cout << "N_1 : " << endl << *nn1 << endl ; 
    
    // nn2 defined as a reference to the Cmp systeme(2).n_auto.c[0]
    p_nn2 = new Cmp( systeme(2).get_n_auto()() ) ;

    // Shift vector
    p_shift1 = new Tenseur( systeme(1).get_shift_auto() ) ;
    p_shift2 = new Tenseur( systeme(2).get_shift_auto() ) ;


  } // End of the black hole case
  
  //*******************************************************************
  // Reading of binary neutron star data
  //*******************************************************************
  else{

    int mer ; 
    fread_be(&mer, sizeof(int), 1, fich) ;	// step 
    
    grille_un = new Mg3d(fich) ;
    map_un = new Map_et(*grille_un, fich) ; 
    Eos* peos_un = Eos::eos_from_file(fich) ; 
    
    grille_deux = new Mg3d(fich) ;
    map_deux = new Map_et(*grille_deux, fich) ; 
    Eos* peos_deux = Eos::eos_from_file(fich) ; 
    
    Binaire star(*map_un, *peos_un, *map_deux, *peos_deux, fich) ; 

    fclose(fich) ; 
    
    bool relativistic = star(1).is_relativistic() ; 

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

    star.display_poly(cout) ; //  Reduced quantities for polytropic EOS

    //==============================================================
    //  Drawings
    //==============================================================

    // int nzdes1 = star(1).get_nzet() ; 
    
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
    Cmp surf1_ext(*map_un) ; 
    surf1_ext = - 0.2 * surf1(0, 0, 0, 0) ; 
    surf1_ext.annule(0, star(1).get_nzet()-1) ; 
    surf1.annule(star(1).get_nzet(), grille_un->get_nzone()-1) ; 
    surf1 = surf1 + surf1_ext ;
    surf1 = raccord_c1(surf1, star(1).get_nzet()) ; 

    Cmp surf2 = star(2).get_ent()() ; 
    Cmp surf2_ext(*map_deux) ; 
    surf2_ext = - 0.2 * surf2(0, 0, 0, 0) ; 
    surf2_ext.annule(0, star(2).get_nzet()-1) ; 
    surf2.annule(star(2).get_nzet(), grille_deux->get_nzone()-1) ; 
    surf2 = surf2 + surf2_ext ;

    surf2 = raccord_c1(surf2, star(2).get_nzet()) ; 

    Tenseur logn1 = - star(1).get_logn_auto() ; 
    Tenseur logn2 = - star(2).get_logn_auto() ; 

    if (graphics) {
      des_coupe_bin_y(logn1(), logn2(), 0, 
			xdes_min, xdes_max, zdes_min, zdes_max, 
		    "ln(N) (y=0)",  &surf1, &surf2) ; 

      des_coupe_bin_z(logn1(), logn2(), 0, 
			xdes_min, xdes_max, ydes_min, ydes_max, 
		    "ln(N) (z=0)",  &surf1, &surf2) ; 


      des_coupe_z(star(1).get_ent()(), 0., 1,
		"Enthalpy (z=0)", &surf1, 1.2) ; 

      des_coupe_y(star(1).get_ent()(), 0., 1,
		"Enthalpy (y=0)", &surf1, 1.2) ; 
    }

    // *nn1 is a copy of the Cmp systeme(1).n_auto.c[0]
    p_nn1 = new Cmp( exp( star(1). get_logn_auto()() ) ) ;
    p_nn1->std_base_scal() ;

    

    // cout << "N_1 : " << endl << *nn1 << endl ; 
    
    // nn2 defined as a reference to the Cmp systeme(2).n_auto.c[0]
    p_nn2 = new Cmp( exp( star(2). get_logn_auto()() ) ) ;
    p_nn2->std_base_scal() ;

    // Shift vector
    p_shift1 = new Tenseur( star(1).get_shift_auto() ) ;
    p_shift2 = new Tenseur( star(2).get_shift_auto() ) ;

    delete peos_un ; 
    delete peos_deux ;

  }  // End of neutron star case

//***************************************************************************
// Begin of the asymptotic study
//***************************************************************************


  /*

    double z0 = 0 ; // cut by the plane z=z0
    double x_min = - 100 ; 
    double x_max = + 100 ; 
    double y_min = - 50 ; 
    double y_max = + 50 ; 

    des_coupe_bin_z(nn1, nn2, z0, x_min, x_max, y_min, y_max,
		    "Lapse function N") ; 

    // Asymptotic behavior of N_1 :
        Valeur** nn1_asymp = nn2.asymptot(3,1) ; 

	// Value (on the angular grid) containing the coef of 1/r^0
	Valeur& nn1_0 = *(nn1_asymp[0]) ;

	// Value (on the angular grid) containing the coef of 1/r
	Valeur& nn1_1 = *(nn1_asymp[1]) ;

	// Value (on the angular grid) containing the coef of 1/r^2
	Valeur& nn1_2 = *(nn1_asymp[2]) ;

	// Value (on the angular grid) containing the coef of 1/r^3
	Valeur& nn1_3 = *(nn1_asymp[3]) ;


	// Computation of spectral expansions
	nn1_0.coef() ; 
	cout << "Spectral coefficients of nn1_0 : " << endl ;
	nn1_0.affiche_seuil(cout,0,4,1e-3) ;
	
	nn1_1.coef() ; 
	cout << "Spectral coefficients of nn1_1 : " << endl ;
	nn1_1.affiche_seuil(cout,0,4,1e-3) ;
	
	nn1_2.coef() ; 
	cout << "Spectral coefficients of nn1_2 : " << endl ;
	nn1_2.affiche_seuil(cout,0,4,1e-3) ;

	nn1_3.coef() ; 
	cout << "Spectral coefficients of nn1_3 : " << endl ;
	nn1_3.affiche_seuil(cout,0,4,1e-5) ;
	
    
	*/
  
	//----------------------------------------------
	// Grid centered on the system "center of mass"
	//----------------------------------------------

	int nz_a = 2 ; // number of domains
	Mg3d mg(nz_a, nr_a, nt_a, np_a, grille_un->get_type_t(), 
		  grille_un->get_type_p(), true) ;

	double bornes_a[3] ;
	bornes_a[0] = 0 ; 
	bornes_a[1] = rayon_a ;
	//        bornes_a[2] = 2*rayon_a ;
	//        bornes_a[3] = 4*rayon_a ;
	//        bornes_a[4] = 8*rayon_a ;
	bornes_a[2] = __infinity ;
	Map_af mp(mg, bornes_a) ;
	cout << "Mapping on the \"centered\" grid : " << endl
	     << mp << endl ;

	// Evaluation of the lapse on the mapping mp
	Cmp nn1(mp) ;
	Cmp nn2(mp) ;

	nn1.import( *(p_nn1) ) ;
	nn2.import( *(p_nn2) ) ;

	nn1.std_base_scal() ;
	nn2.std_base_scal() ;

	Cmp nn(mp) ; 
	if (black_hole) {  // N = N_1 + N_2  for black holes
	    nn = nn1 + nn2  ;
	}
	else{             // ln(N) = ln(N_1) + ln(N_2) for neutron stars
	    nn = nn1*nn2 ;
	}
	
	// The lapse is set to 1 in the inner domain
	nn.annule(0) ;
	Cmp tmp(mp) ;
	tmp = 1 ; 
	tmp.annule(1) ;	
	nn = nn + tmp ; 

	// Logarithm of the lapse
	Cmp logn = log( nn ) ;
	logn.std_base_scal() ;
		

	// Evaluation of the shift vector on the mapping mp
	
	// Ensures that the components of shift1 and shift2 are
	//  defined with respect to the same triad (this was
	//  true for binary NS but not for BH)
	p_shift1->change_triad( mp.get_bvect_cart() ) ;
	p_shift2->change_triad( mp.get_bvect_cart() ) ;

	Tenseur shift1(mp, 1, CON, mp.get_bvect_cart()) ;
	Tenseur shift2(mp, 1, CON, mp.get_bvect_cart()) ;

	shift1.set_etat_qcq() ; 
	shift2.set_etat_qcq() ; 
	for (int i=0; i<3; i++) {
	    (shift1.set(i)).import( (*p_shift1)(i) ) ; 
	    (shift2.set(i)).import( (*p_shift2)(i) ) ; 
	}

	// Spectral bases for the x-component and y-component
	for (int i=0; i<2; i++) {
	    shift1.set(i).va.set_base_r(0,R_CHEBPIM_P) ;
	    shift1.set(i).va.set_base_r(1,R_CHEBU) ;
	    shift1.set(i).va.set_base_t(T_COSSIN_CP) ;
	    shift1.set(i).va.set_base_p(P_COSSIN) ;

	    shift2.set(i).va.set_base_r(0,R_CHEBPIM_P) ;
	    shift2.set(i).va.set_base_r(1,R_CHEBU) ;
	    shift2.set(i).va.set_base_t(T_COSSIN_CP) ;
	    shift2.set(i).va.set_base_p(P_COSSIN) ;
	}
	// Spectral bases for the z-component
	shift1.set(2).va.set_base_r(0,R_CHEBPIM_I) ;
	shift1.set(2).va.set_base_r(1,R_CHEBU) ;
	shift1.set(2).va.set_base_t(T_COSSIN_CI) ;
	shift1.set(2).va.set_base_p(P_COSSIN) ;
	
	shift2.set(2).va.set_base_r(0,R_CHEBPIM_I) ;
	shift2.set(2).va.set_base_r(1,R_CHEBU) ;
	shift2.set(2).va.set_base_t(T_COSSIN_CI) ;
	shift2.set(2).va.set_base_p(P_COSSIN) ;
	

	Tenseur shift = shift1 + shift2 ; 

	//--------------------------------------------
	// Graphical outputs
	//-------------------------------------------

	if (graphics) {
	  des_coupe_z(nn1, 0., 1, "nn1") ;
	 
	  des_coupe_vect_z(shift, 0., -3., 1., 1,"shift vector") ;

	}

	Cmp shift_x  = shift.set(0) ;
	Cmp shift_y  = shift.set(1) ;
	Cmp shift_z  = shift.set(2) ;
	
	if (graphics) {
	  des_coupe_z(shift_y, 0., 1, "shift_y") ;
	  des_coupe_z(shift_x, 0., 1, "shift_x") ;
	  des_coupe_z(shift_z, 0., 1, "shift_z") ;
	}

	cout << "msol" << 4.62/msol << endl ;

	ofstream fichresu("resu_asymptot.d") ; 
	fichresu << "nom de fichier :" << nom_config << endl
		<< "nr_a : number of points in r for the centered grid = "
		<< nr_a << endl
                <<"nt_a : number of points in theta for the center grid = " 
                << nt_a  << endl
		<< "np_a : number of points in phi for the centered grid = "
		<< np_a << endl
		<<"rayon_a (10 km) : inner radius of the last domain = " 
		 << rayon_a << endl <<endl << endl ; 


	asymptot(nn, "N", graphics, fichresu) ;
	asymptot(logn, "logn", graphics, fichresu) ;
	asymptot(shift_x, "shift_x", graphics, fichresu) ;
        asymptot(shift_y, "shift_y", graphics, fichresu) ; 	   
	asymptot(shift_z, "shift_z", graphics, fichresu) ;
 	
	fichresu.close() ; 

	

    // Freeing memory

	  //   delete [] nn1_asymp ;
	//  delete [] nn_asymp ;

      delete p_shift1 ;
      delete p_shift2 ;
      delete p_nn1 ;
      delete p_nn2 ;
      delete map_un ;
      delete map_deux ;
      delete grille_un ;
      if (!black_hole) delete grille_deux ;


//home/francois/EUNetwork/Meudon/Data/BinNS/GR/irrotation/G18vs18_g2_ir_M.d
  return EXIT_SUCCESS ; 

}



