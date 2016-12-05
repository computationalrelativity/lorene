/*
 * Constructor of class Mag_NS (magnetized neutron star exportation)
 * which depends explicitely on Lorene objects.
 *
 * (see file mag_ns.h for documentation).
 */

/*
 *   Copyright (c) 2002  Eric Gourgoulhon
 *   Copyright (c) 2009  Jerome Novak
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:18:31  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2016/02/17 09:46:35  j_novak
 * u_euler output in units of c.
 *
 * Revision 1.3  2014/10/13 08:54:06  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/05/13 10:06:45  j_novak
 * Update to take unto account the change in Lorene magnetic units.
 *
 * Revision 1.1  2009/11/19 16:15:21  j_novak
 * Export class for magnetized neutron stars.
 *
 *
 * $Header$
 *
 */

#include "../Include/mag_ns.h"

// C headers
#include <cstring>
#include <cmath>

// Lorene headers
#include "tenseur.h"
#include "et_rot_mag.h"
#include "eos.h"
#include "unites.h"
#include "metric.h"
#include "utilitaires.h"

		    //----------------------------------------//
		    //	    Constructor from LORENE data      //
		    //----------------------------------------//

namespace Lorene {
Mag_NS::Mag_NS(int nbpoints, const double* xi, const double* yi,
	       const double* zi, const char* filename)
	       : np(nbpoints) {

  using namespace Unites ;
  using namespace Unites_mag ;

    // Reading of data
    // ---------------
    FILE* fich = fopen(filename, "r") ;

    Mg3d spectral_grid(fich) ;
    int nphi = spectral_grid.get_np(0) ;
    if (nphi < 4) {
	cout << "Contructor of Mag_NS: " << endl ;
	cout << "Fatal problem with initial data, "<< endl ;
	cout << "the number of points in the phi direction is lower than 4."
	     << endl ;
	cout << "Impossible to rotate tensors to Cartesian triad. Please give as"
	     << endl ;
	cout << "an input a result file from magstar with nphi >=4." << endl ;
	abort() ;
    }
    Map_et mapping(spectral_grid, fich) ;
    Eos* p_eos = Eos::eos_from_file(fich) ;

    Et_rot_mag star(mapping, *p_eos, fich) ;
    star.equation_of_state() ;
    star.update_metric() ;
    star.extrinsic_curvature() ;
    star.hydro_euler() ;

    const Map& mp = star.get_mp();
    
    const Base_vect_spher& bspher = mp.get_bvect_spher() ;
    Sym_tensor gij(mp, COV, bspher) ;
    gij.set_etat_zero() ;
    gij.set(1,1) = star.get_a_car()() ;
    gij.set(2,2) = star.get_a_car()() ;
    gij.set(3,3) = star.get_b_car()() ;
    Metric gam(gij) ; //the 3-metric in quasi-isotropic coordinates

    // Initialisation of member data
    // -----------------------------

    const Eos_poly* p_eos_poly = dynamic_cast<const Eos_poly*>( p_eos ) ;
    double m_0 = 1. ;
    if ( p_eos_poly != 0x0 ) {
        strcpy(eos_name, "Polytropic EOS") ;
	gamma_poly = p_eos_poly->get_gam() ;
	kappa_poly = p_eos_poly->get_kap() ;
	m_0 = p_eos_poly->get_m_0() ;
    }
    else {
        strncpy(eos_name, (star.get_eos()).get_name(), 100) ;
        gamma_poly = 0 ;
        kappa_poly = 0 ;
	cout << "Mag_NS::Mag_NS : not a ploytropic EOS" << endl ;
	cout << "WARNING!!" << endl ;
	cout << "exporting the baryon density instead of the mass density!" << endl ;
	arrete() ;
    }

    Scalar sp_density = m_0 * star.get_nbar()()*rho_unit ;
    Scalar sp_energy = star.get_ener()() * rho_unit;
    omega = star.get_omega_c() * f_unit ;
    rho_c = m_0*sp_density.val_grid_point(0,0,0,0) ;
    eps_c = sp_energy.val_grid_point(0,0,0,0) / rho_c - 1. ;
    mass_b = star.mass_b() / msol ;
    mass_g = star.mass_g() / msol ;
    r_eq = star.ray_eq() / km ;
    r_p = star.ray_pole() / km ;
    angu_mom = star.angu_mom() / ( ggrav * msol*msol) ;
    T_over_W =  star.tsw() ;
    magn_mom = star.MagMom() ;
    b_z_pole = 
      star.Magn()(0).va.val_point(star.l_surf()(0,0), star.xi_surf()(0,0),0.,0.) 
      * mag_unit / 1.e9 ;
    int theta_eq = mapping.get_mg()->get_nt(star.get_nzet() - 1) - 1 ;
    b_z_eq = star.Magn()(1).va.val_point(star.l_surf()(0,theta_eq),
					 star.xi_surf()(0,theta_eq),M_PI_2,0.)
      * mag_unit / 1.e9 ;
    
    cout.precision(13) ;
    cout << endl << "Magnetized star read in file : " << endl ;
    cout <<	    "------------------------------ " << endl ;
    cout << star << endl ;
    cout << endl << "Summary : " << endl ;
    cout << "-------" << endl ;
    cout << "  Omega          : " << omega << " rad/s" << endl ;
    cout << "  Baryon mass : " << mass_b
	 << " M_sol" << endl ;
    cout << "  Gravitational mass : " << mass_g
	 << " M_sol" << endl ;
    cout << "  Equatorial radius : " << r_eq
         << " km" << endl ;
    cout << "  Polar radius : " << r_p
         << " km" << endl ;
    cout << "  Total angular momentum : " << angu_mom
	 << " G M_sol^2 / c" << endl ;
    cout << "  T/W : " << T_over_W << endl ;
    cout << "  Magnetic momentum : " << magn_mom
	 << " A/m^2" << endl ;
    cout << " Radial magnetic field polar value  : " << b_z_pole
	 << " GT" << endl ;
    cout << "  Tangent magnetic field equatorial value  : " << b_z_eq
	 << " GT" << endl ;

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
    assert(star.is_relativistic()) ;

    Scalar sp_lapse = star.get_nnn()() ;
    Vector sp_shift(mapping, CON, mapping.get_bvect_cart()) ;
    sp_shift.set(1) = star.get_shift()(0) ; //contravariant representation!!
    sp_shift.set(2) = star.get_shift()(1) ;
    sp_shift.set(3) = star.get_shift()(2) ;

    Sym_tensor sp_gamma(mapping, COV, mapping.get_bvect_spher()) ;
    sp_gamma.set_etat_qcq() ;
    for (int i=1; i<=3; i++) {
	for (int j=1; j<=i; j++)
	    sp_gamma.set(i,j) = 0 ;
	if (i != 3) sp_gamma.set(i,i) = star.get_a_car()() ;
	else sp_gamma.set(3,3) = star.get_bbb()()*star.get_bbb()() ;
    }
    sp_gamma.change_triad(mapping.get_bvect_cart()) ;

    Sym_tensor sp_kij(mapping, COV, mapping.get_bvect_cart()) ;
    sp_kij.set_etat_qcq() ;
    for (int i=0; i<3; i++)
	for (int j=i; j<3; j++) {
	    sp_kij.set(i+1, j+1) =
		star.get_bbb()()*star.get_bbb()()*star.get_tkij()(i,j) / r_unit ;
	    sp_kij.set(i+1, j+1).dec_dzpuis(2) ;
	}
    
    Vector sp_u_euler(mapping, CON, mapping.get_bvect_cart()) ;
    sp_u_euler.set(1) = star.get_u_euler()(0) ;
    sp_u_euler.set(2) = star.get_u_euler()(1) ;
    sp_u_euler.set(3) = star.get_u_euler()(2) ;

    Vector sp_current(mp, CON, mapping.get_bvect_spher()) ;
    sp_current.set(1) = 0 ;
    sp_current.set(2) = 0 ;
    sp_current.set(3) = Scalar(star.get_jphi()*j_unit) ;
    sp_current.set(3).mult_rsint() ; 
    sp_current.change_triad(mapping.get_bvect_cart()) ;
    Scalar sp_jt = star.get_jt() * j_unit ;

    Scalar fac = 1.e9*sqrt(star.get_a_car()())/mag_unit ;//to transform B^(i) into B^i
    fac.std_spectral_base() ;
    Vector sp_Bmag(mp, CON, bspher) ;
    sp_Bmag.set(1) = Scalar(star.Magn()(0)) / fac ;
    sp_Bmag.set(1).dec_dzpuis(2) ;
    sp_Bmag.set(2) = Scalar(star.Magn()(1)) / fac ;
    sp_Bmag.set(2).dec_dzpuis(2) ;
    sp_Bmag.set(3) = 0 ;
    sp_Bmag.change_triad(mapping.get_bvect_cart()) ;

    for (int i=0; i<np; i++) {

	double x0 = xx[i] * km ;    // x, y, z in Lorene's unit
	double y0 = yy[i] * km ;
	double z0 = zz[i] * km ;

	// Values of (l, xi, theta, phi) corresponding to (x,y,z):
	// -------------------------------------------------------
	double r, theta, phi ;   // polar coordinates centered on the star
	mapping.convert_absolute(x0, y0, z0, r, theta, phi) ;

	// Lapse function
	// --------------
	
	nnn[i] = sp_lapse.val_point(r, theta, phi) ;
	
	// Shift vector
	// ------------
	
	beta_x[i] = sp_shift(1).val_point(r, theta, phi) ;
	beta_y[i] = sp_shift(2).val_point(r, theta, phi) ;
	beta_z[i] = sp_shift(3).val_point(r, theta, phi) ;

	// 3-metric
	// ---------

	g_xx[i] = sp_gamma(1,1).val_point(r, theta, phi) ;
	g_yy[i] = sp_gamma(2,2).val_point(r, theta, phi) ;
	g_zz[i] = sp_gamma(3,3).val_point(r, theta, phi) ;
	g_xy[i] = sp_gamma(1,2).val_point(r, theta, phi) ;
	g_xz[i] = sp_gamma(1,3).val_point(r, theta, phi) ;
	g_yz[i] = sp_gamma(2,3).val_point(r, theta, phi) ;

	// Extrinsic curvature
	// -------------------

	k_xx[i] = sp_kij(1,1).val_point(r, theta, phi) ;
	k_xy[i] = sp_kij(1,1).val_point(r, theta, phi) ;
	k_xz[i] = sp_kij(1,1).val_point(r, theta, phi) ;
	k_yy[i] = sp_kij(1,1).val_point(r, theta, phi) ;
	k_yz[i] = sp_kij(1,1).val_point(r, theta, phi) ;
	k_zz[i] = sp_kij(1,1).val_point(r, theta, phi) ;

	// Baryon density [kg/m^3]
	// --------------
	
	nbar[i] = sp_density.val_point(r, theta, phi) ;

	// Energy density
	// --------------

	double ener = sp_energy.val_point(r, theta, phi) ;

        if ( nbar[i] == double(0) ) {
                ener_spec[i] = 0 ;
        }
        else {
                ener_spec[i] = ener / nbar[i] - double(1) ;
        }
        
        // 3-velocity with respect to the Eulerian observer
        // ------------------------------------------------
	u_euler_x[i] = sp_u_euler(1).val_point(r, theta, phi) ;
	u_euler_y[i] = sp_u_euler(2).val_point(r, theta, phi) ;
	u_euler_z[i] = sp_u_euler(3).val_point(r, theta, phi) ;

	// Magnetic field
	//---------------
	bb_x[i] = sp_Bmag(1).val_point(r, theta, phi) ;
	bb_y[i] = sp_Bmag(2).val_point(r, theta, phi) ;
	bb_z[i] = sp_Bmag(3).val_point(r, theta, phi) ;

	// 4-current
	//----------
	jj_t[i] = sp_jt.val_point(r, theta, phi) ;
	jj_x[i] = sp_current(1).val_point(r, theta, phi) ;
	jj_y[i] = sp_current(2).val_point(r, theta, phi) ;
	jj_z[i] = sp_current(3).val_point(r, theta, phi) ;


    }	// End of loop on the points


     delete p_eos ;
}

}
