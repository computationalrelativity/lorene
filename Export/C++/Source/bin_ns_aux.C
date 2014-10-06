/*
 * Constructor of class Bin_NS (binary neutron star exportation)
 * which depends explicitely on Lorene objects.
 *
 * (see file bin_ns.h for documentation).
 */

/*
 *   Copyright (c) 2002  Eric Gourgoulhon
 *   Copyright (c) 2002  Keisuke Taniguchi
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

char bin_ns_aux_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2014/10/06 15:13:25  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2014/02/24 16:08:11  e_gourgoulhon
 * Eulerian velocities set to zero outside the stars
 *
 * Revision 1.5  2006/09/12 08:04:07  j_novak
 * Removal of the include path Export/C++/Include, updating of the relevant
 * source files in Export/C++/Source.
 *
 * Revision 1.4  2004/04/29 20:29:18  e_gourgoulhon
 * Corrected a bug in the computation of ener_spec.
 *
 * Revision 1.3  2004/03/25 10:29:27  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.2  2002/01/15 09:11:04  e_gourgoulhon
 * Display of Lorene object (class binaire) read in file
 *
 * Revision 1.1  2002/01/11 17:03:02  e_gourgoulhon
 * Exportation of binary neutron stars configuration to a Cartesian grid
 *
 *
 * $Header$
 *
 */

#include "../Include/bin_ns.h"

// C headers
#include <cstring>
#include <cmath>

// Lorene headers
#include "tenseur.h"
#include "binaire.h"
#include "eos.h"
#include "unites.h"

		    //----------------------------------------//
		    //	    Constructor from LORENE data      //
		    //----------------------------------------//

Bin_NS::Bin_NS(int nbpoints, const double* xi, const double* yi,
	       const double* zi, const char* filename)
	       : np(nbpoints) {

  using namespace Unites ;

    // Reading of data
    // ---------------
    FILE* fich = fopen(filename, "r") ;

    int mer ;
    fread(&mer, sizeof(int), 1, fich) ; // mer

    Mg3d mg1(fich) ;
    Map_et mp1(mg1, fich) ; 
    Eos* peos1 = Eos::eos_from_file(fich) ; 

    Mg3d mg2(fich) ;
    Map_et mp2(mg2, fich) ; 
    Eos* peos2 = Eos::eos_from_file(fich) ; 
    
    Binaire star(mp1, *peos1, mp2, *peos2, fich) ; 

    fclose(fich) ;
    
    if ( (star(1).is_relativistic() == false) ||
         (star(2).is_relativistic() == false)    ) {
        cout << "Bin_NS::Bin_NS can import only relativistic stars !" << endl ;
        abort() ;
    }

    // Construction of the binary system
    // ---------------------------------
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

    // Initialisation of member data
    // -----------------------------

    const Eos_poly* p_eos_poly1 = dynamic_cast<const Eos_poly*>( peos1 ) ;
    const Eos_poly* p_eos_poly2 = dynamic_cast<const Eos_poly*>( peos2 ) ;

    if ( p_eos_poly1 != 0x0 ) {
        strcpy(eos_name1, "Polytropic EOS") ;
	gamma_poly1 = p_eos_poly1->get_gam() ;
	kappa_poly1 = p_eos_poly1->get_kap() ;
    }
    else {
        strncpy(eos_name1, (star(1).get_eos()).get_name(), 100) ;
        gamma_poly1 = 0 ;
        kappa_poly1 = 0 ;
    }

    if ( p_eos_poly2 != 0x0 ) {
        strcpy(eos_name2, "Polytropic EOS") ;
	gamma_poly2 = p_eos_poly2->get_gam() ;
	kappa_poly2 = p_eos_poly2->get_kap() ;
    }
    else {
        strncpy(eos_name2, (star(2).get_eos()).get_name(), 100) ;
        gamma_poly2 = 0 ;
        kappa_poly2 = 0 ;
    }

    omega = star.get_omega() * f_unit ;
    dist = star.separation() / km ;
    dist_mass = ( star(2).xa_barycenter()
		  - star(1).xa_barycenter() ) / km ;
    mass1_b = star(1).mass_b() / msol ;
    mass2_b = star(2).mass_b() / msol ;
    mass_adm = star.mass_adm() / msol ;
    angu_mom = star.angu_mom()(2) / ( ggrav * msol*msol) ;
    rad1_x_comp = star(1).ray_eq() / km ;
    rad1_y = star(1).ray_eq_pis2() / km ;
    rad1_z = star(1).ray_pole() / km ;
    rad1_x_opp = star(1).ray_eq_pi() / km ;
    rad2_x_comp = star(2).ray_eq() / km ;
    rad2_y = star(2).ray_eq_pis2() / km ;
    rad2_z = star(2).ray_pole() / km ;
    rad2_x_opp = star(2).ray_eq_pi() / km ;

    cout.precision(13) ;
    cout << endl << "Binary system read in file : " << endl ;
    cout <<	    "---------------------------- " << endl ;
    cout << star << endl ;
    cout << endl << "Summary : " << endl ;
    cout << "-------" << endl ;
    cout << "  Separation d   : " << dist << " km" << endl ;
    cout << "  Separation d_G : " << dist_mass << " km" << endl ;
    cout << "  Omega          : " << omega << " rad/s" << endl ;
    cout << "  Baryon mass of star 1  : " << mass1_b
	 << " M_sol" << endl ;
    cout << "  Baryon mass of star 2  : " << mass2_b
	 << " M_sol" << endl ;
    cout << "  ADM mass of the system : " << mass_adm
         << " M_sol" << endl ;
    cout << "  Total angular momentum : " << angu_mom
	 << " G M_sol^2 / c" << endl ;
    cout << "  Radius of star 1 (x_comp) : " << rad1_x_comp
	 << " km" << endl ;
    cout << "  Radius of star 1 (y)      : " << rad1_y
	 << " km" << endl ;
    cout << "  Radius of star 1 (z)      : " << rad1_z
	 << " km" << endl ;
    cout << "  Radius of star 1 (x_opp)  : " << rad1_x_opp
	 << " km" << endl ;
    cout << "  Radius of star 2 (x_comp) : " << rad2_x_comp
	 << " km" << endl ;
    cout << "  Radius of star 2 (y)      : " << rad2_y
	 << " km" << endl ;
    cout << "  Radius of star 2 (z)      : " << rad2_z
	 << " km" << endl ;
    cout << "  Radius of star 2 (x_opp)  : " << rad2_x_opp
	 << " km" << endl ;

    // Creation of the various arrays on the Cartesian grid
    // ----------------------------------------------------

    alloc_memory() ;

    // Initialisation of the Cartesian grid
    // ------------------------------------

    for (int i=0; i<np; i++) {
	xx[i] = xi[i] ;
    }
    for (int i=0; i<np; i++) {
	yy[i] = yi[i] ;
    }
    for (int i=0; i<np; i++) {
	zz[i] = zi[i] ;
    }

    // Computation of the values at the points of the Cartesian grid
    // -------------------------------------------------------------

    // Parameter 1/c^2
    double unsurc2 = star(1).is_relativistic() ? double(1) : double(0) ;

    Tenseur lapse1 = exp( unsurc2 * star(1).get_logn_auto() ) ;
    lapse1.set_std_base() ;
    Tenseur lapse2 = exp( unsurc2 * star(2).get_logn_auto() ) ;
    lapse2.set_std_base() ;

    Valeur vnn1 = lapse1().va ;
    Valeur vnn2 = lapse2().va ;
    vnn1.coef() ;
    vnn2.coef() ;

    const Valeur& vshiftx1 = (star(1).get_shift_auto()(0)).va ;
    const Valeur& vshiftx2 = (star(2).get_shift_auto()(0)).va ;
    const Valeur& vshifty1 = (star(1).get_shift_auto()(1)).va ;
    const Valeur& vshifty2 = (star(2).get_shift_auto()(1)).va ;
    const Valeur& vshiftz1 = (star(1).get_shift_auto()(2)).va ;
    const Valeur& vshiftz2 = (star(2).get_shift_auto()(2)).va ;
    vshiftx1.coef() ;
    vshiftx2.coef() ;
    vshifty1.coef() ;
    vshifty2.coef() ;
    vshiftz1.coef() ;
    vshiftz2.coef() ;

    Tenseur a_car1 = exp( 2.*unsurc2*( star(1).get_beta_auto()
				       - star(1).get_logn_auto() ) ) ;
    a_car1.set_std_base() ;
    Tenseur a_car2 = exp( 2.*unsurc2*( star(2).get_beta_auto()
				       - star(2).get_logn_auto() ) ) ;
    a_car2.set_std_base() ;

    Valeur vacar1 = a_car1().va ;
    Valeur vacar2 = a_car2().va ;
    vacar1.coef() ;
    vacar2.coef() ;

    Tenseur_sym k_one(star(1).get_tkij_auto()) ;
    k_one.set_std_base() ;
    k_one.dec2_dzpuis() ;
    Tenseur_sym k_two(star(2).get_tkij_auto()) ;
    k_two.set_std_base() ;
    k_two.dec2_dzpuis() ;

    Valeur vkxx1 = (k_one(0, 0)).va ;
    Valeur vkxx2 = (k_two(0, 0)).va ;
    Valeur vkxy1 = (k_one(0, 1)).va ;
    Valeur vkxy2 = (k_two(0, 1)).va ;
    Valeur vkxz1 = (k_one(0, 2)).va ;
    Valeur vkxz2 = (k_two(0, 2)).va ;
    Valeur vkyy1 = (k_one(1, 1)).va ;
    Valeur vkyy2 = (k_two(1, 1)).va ;
    Valeur vkyz1 = (k_one(1, 2)).va ;
    Valeur vkyz2 = (k_two(1, 2)).va ;
    Valeur vkzz1 = (k_one(2, 2)).va ;
    Valeur vkzz2 = (k_two(2, 2)).va ;

    vkxx1.set_base( (star(1).get_tkij_auto()(0, 0)).va.base ) ;
    vkxx2.set_base( (star(2).get_tkij_auto()(0, 0)).va.base ) ;
    vkxy1.set_base( (star(1).get_tkij_auto()(0, 1)).va.base ) ;
    vkxy2.set_base( (star(2).get_tkij_auto()(0, 1)).va.base ) ;
    vkxz1.set_base( (star(1).get_tkij_auto()(0, 2)).va.base ) ;
    vkxz2.set_base( (star(2).get_tkij_auto()(0, 2)).va.base ) ;
    vkyy1.set_base( (star(1).get_tkij_auto()(1, 1)).va.base ) ;
    vkyy2.set_base( (star(2).get_tkij_auto()(1, 1)).va.base ) ;
    vkyz1.set_base( (star(1).get_tkij_auto()(1, 2)).va.base ) ;
    vkyz2.set_base( (star(2).get_tkij_auto()(1, 2)).va.base ) ;
    vkzz1.set_base( (star(1).get_tkij_auto()(2, 2)).va.base ) ;
    vkzz2.set_base( (star(2).get_tkij_auto()(2, 2)).va.base ) ;
    vkxx1.coef() ;
    vkxx2.coef() ;
    vkxy1.coef() ;
    vkxy2.coef() ;
    vkxz1.coef() ;
    vkxz2.coef() ;
    vkyy1.coef() ;
    vkyy2.coef() ;
    vkyz1.coef() ;
    vkyz2.coef() ;
    vkzz1.coef() ;
    vkzz2.coef() ;

    const Valeur& vnbar1 = (star(1).get_nbar()()).va ;
    const Valeur& vnbar2 = (star(2).get_nbar()()).va ;
    vnbar1.coef() ;
    vnbar2.coef() ;

    const Valeur& vener1 = (star(1).get_ener()()).va ;
    const Valeur& vener2 = (star(2).get_ener()()).va ;
    vener1.coef() ;
    vener2.coef() ;

    Valeur vueulerx1 = (star(1).get_u_euler()(0)).va ;
    Valeur vueulerx2 = (star(2).get_u_euler()(0)).va ;
    Valeur vueulery1 = (star(1).get_u_euler()(1)).va ;
    Valeur vueulery2 = (star(2).get_u_euler()(1)).va ;
    Valeur vueulerz1 = (star(1).get_u_euler()(2)).va ;
    Valeur vueulerz2 = (star(2).get_u_euler()(2)).va ;

    // Velocities set to zero outside the stars:
    int nzet1 = star(1).get_nzet() ; 
    int nz1m1 = mg1.get_nzone() - 1 ; 
    int nzet2 = star(2).get_nzet() ; 
    int nz2m1 = mg2.get_nzone() - 1 ; 
    vueulerx1.annule(nzet1, nz1m1) ; 
    vueulery1.annule(nzet1, nz1m1) ; 
    vueulerz1.annule(nzet1, nz1m1) ; 
    vueulerx2.annule(nzet2, nz2m1) ; 
    vueulery2.annule(nzet2, nz2m1) ; 
    vueulerz2.annule(nzet2, nz2m1) ; 
    
    vueulerx1.coef() ;
    vueulerx2.coef() ;
    vueulery1.coef() ;
    vueulery2.coef() ;
    vueulerz1.coef() ;
    vueulerz2.coef() ;

    for (int i=0; i<np; i++) {

	double x0 = xx[i] * km ;    // x in Lorene's unit
	double y0 = yy[i] * km ;
	double z0 = zz[i] * km ;

	// Values of (l1, xi1, theta1, phi1) (grid 1)
	// corresponding to (x,y,z):
	// ------------------------------------------
	double r1, theta1, phi1 ;   // polar coordinates centered on NS 1
	mp1.convert_absolute(x0, y0, z0, r1, theta1, phi1) ;

	int l1 ;	    // domain index
	double xi1 ;	    // radial coordinate xi in [0,1] or [-1,1]
	mp1.val_lx(r1, theta1, phi1, l1, xi1) ;

	// Values of (l2, xi2, theta2, phi2) (grid 2) 
	// corresponding to (x,y,z):
	// ------------------------------------------
	double r2, theta2, phi2 ;   // polar coordinates centered on b.h. 2
	mp2.convert_absolute(x0, y0, z0, r2, theta2, phi2) ; 
	
	int l2 ;	    // domain index
	double xi2 ;	    // radial coordinate xi in [0,1] or [-1,1]
	mp2.val_lx(r2, theta2, phi2, l2, xi2) ;
	
	// Lapse function
	// --------------
	
	nnn[i] = vnn1.c_cf->val_point_symy(l1, xi1, theta1, phi1) 
		 * vnn2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ;
	
	// Shift vector
	// ------------
	
	beta_x[i] = vshiftx1.c_cf->val_point_asymy(l1, xi1, theta1, phi1)
		 +  vshiftx2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ;

	beta_y[i] = vshifty1.c_cf->val_point_symy(l1, xi1, theta1, phi1)
		 +  vshifty2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ;
	//		 + omega * xx[i] ;

	beta_z[i] = vshiftz1.c_cf->val_point_asymy(l1, xi1, theta1, phi1)
		 +  vshiftz2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ;

	// Conformal factor
	// ----------------

	double psi4 = vacar1.c_cf->val_point_symy(l1, xi1, theta1, phi1)
	           * vacar2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ;

	g_xx[i] = psi4 ;
	g_yy[i] = psi4 ;
	g_zz[i] = psi4 ;
	g_xy[i] = 0 ;
	g_xz[i] = 0 ;
	g_yz[i] = 0 ;

	// Extrinsic curvature
	// -------------------

	double pre = km * psi4 ;

	k_xx[i] = pre*(  vkxx1.c_cf->val_point_asymy(l1, xi1, theta1, phi1)
		       +vkxx2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ) ;

	k_xy[i] = pre*(  vkxy1.c_cf->val_point_symy(l1, xi1, theta1, phi1)
		       +vkxy2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ) ;

	k_xz[i] = pre*(  vkxz1.c_cf->val_point_asymy(l1, xi1, theta1, phi1)
		       +vkxz2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ) ;

	k_yy[i] = pre*(  vkyy1.c_cf->val_point_asymy(l1, xi1, theta1, phi1)
		       +vkyy2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ) ;

	k_yz[i] = pre*(  vkyz1.c_cf->val_point_symy(l1, xi1, theta1, phi1)
		       +vkyz2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ) ;

	k_zz[i] = pre*(  vkzz1.c_cf->val_point_asymy(l1, xi1, theta1, phi1)
		       +vkzz2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ) ;

	// Baryon density [kg/m^3]
	// --------------
	
	nbar[i] = rho_unit*(
			  vnbar1.c_cf->val_point_symy(l1, xi1, theta1, phi1)
			+ vnbar2.c_cf->val_point_symy(l2, xi2, theta2, phi2)
			) ;

	// Energy density
	// --------------

	double ener = vener1.c_cf->val_point_symy(l1, xi1, theta1, phi1)
		 + vener2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ;

        if ( nbar[i] == double(0) ) {
                ener_spec[i] = 0 ;
        }
        else {
                ener_spec[i] = ener / (nbar[i]/rho_unit) - double(1) ;
        }
        

        // 3-velocity with respect to the Eulerian observer
        // ------------------------------------------------
	u_euler_x[i] = vueulerx1.c_cf->val_point(l1, xi1, theta1, phi1)
	            +  vueulerx2.c_cf->val_point(l2, xi2, theta2, phi2) ;

	u_euler_y[i] = vueulery1.c_cf->val_point(l1, xi1, theta1, phi1)
	            +  vueulery2.c_cf->val_point(l2, xi2, theta2, phi2) ;
	//		 + omega * xx[i] ;

	u_euler_z[i] = vueulerz1.c_cf->val_point(l1, xi1, theta1, phi1)
	            +  vueulerz2.c_cf->val_point(l2, xi2, theta2, phi2) ;



    }	// End of loop on the points


     delete peos1 ;
     delete peos2 ;
}

