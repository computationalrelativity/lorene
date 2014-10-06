/*
 *  Method of class Bin_bhns_extr to compute the orbital angular velocity
 *  {\tt omega}
 *
 *    (see file bin_bhns_extr.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004-2005 Keisuke Taniguchi
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

char bin_bhns_extr_orbit_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:13:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2005/02/28 23:07:47  k_taniguchi
 * Addition of the code for the conformally flat case
 *
 * Revision 1.1  2004/11/30 20:46:57  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include <cmath>

// Lorene headers
#include "bin_bhns_extr.h"
#include "eos.h"
#include "param.h"
#include "utilitaires.h"
#include "unites.h"

double fonc_bhns_orbit_ks(double, const Param& ) ;
double fonc_bhns_orbit_cf(double, const Param& ) ;

//************************************************************************

void Bin_bhns_extr::orbit_omega_ks(double fact_omeg_min,
				   double fact_omeg_max) {

  using namespace Unites ;

    //------------------------------------------------------------------
    // Evaluation of various quantities at the center of a neutron star
    //------------------------------------------------------------------

    double dnulg, asn2, dasn2, nx, dnx, ny, dny, npn, dnpn ;
    double msr ;

    const Map& mp = star.get_mp() ;

    const Cmp& logn_auto_regu = star.get_logn_auto_regu()() ;
    const Cmp& loggam = star.get_loggam()() ;
    const Cmp& nnn = star.get_nnn()() ;
    const Cmp& a_car = star.get_a_car()() ;
    const Tenseur& shift = star.get_shift() ;
    const Tenseur& d_logn_auto_div = star.get_d_logn_auto_div() ;

    Tenseur dln_auto_div = d_logn_auto_div ;

    if ( *(dln_auto_div.get_triad()) != ref_triad ) {

        // Change the basis from spherical coordinate to Cartesian one
        dln_auto_div.change_triad( mp.get_bvect_cart() ) ;

	// Change the basis from mapping coordinate to absolute one
	dln_auto_div.change_triad( ref_triad ) ;

    }

    const Tenseur& d_logn_comp = star.get_d_logn_comp() ;

    Tenseur dln_comp = d_logn_comp ;

    if ( *(dln_comp.get_triad()) != ref_triad ) {

        // Change the basis from spherical coordinate to Cartesian one
        dln_comp.change_triad( mp.get_bvect_cart() ) ;

	// Change the basis from mapping coordinate to absolute one
	dln_comp.change_triad( ref_triad ) ;

    }

    //-------------------------------
    // Parameters with respect to BH
    //-------------------------------

    msr = ggrav * mass_bh / separ ;

    //-----------------------------------------------------------------
    // Calculation of d/dX( nu + ln(Gamma) ) at the center of the star
    //  ---> dnulg
    //-----------------------------------------------------------------

    // Factor for the coordinate transformation x --> X :
    double factx ;
    if (fabs(mp.get_rot_phi()) < 1.e-14) {
	factx = 1. ;
    }
    else {
	if (fabs(mp.get_rot_phi() - M_PI) < 1.e-14) {
	    factx = - 1. ;
	}
	else {
	    cout << "Bin_bhns_extr::orbit_omega : unknown value of rot_phi !"
		 << endl ;
	    abort() ;
	}
    }

    Cmp tmp = logn_auto_regu + loggam ;

    // ... gradient
    dnulg = dln_auto_div(0)(0, 0, 0, 0) + dln_comp(0)(0, 0, 0, 0)
      + factx * tmp.dsdx()(0, 0, 0, 0) ;

    //------------------------------------------------------------
    // Calculation of A^2/N^2 at the center of the star ---> asn2
    //------------------------------------------------------------
    double nc = nnn(0, 0, 0, 0) ;
    double a2c = a_car(0, 0, 0, 0) ;
    asn2 = a2c / (nc * nc) ;

    if ( star.is_relativistic() ) {

        //------------------------------------------------------------------
	// Calculation of d/dX(A^2/N^2) at the center of the star ---> dasn
	//------------------------------------------------------------------

        double da2c = factx * a_car.dsdx()(0, 0, 0, 0) ; 
	double dnc =  factx * nnn.dsdx()(0, 0, 0, 0) ;

	dasn2 = ( da2c - 2. * a2c / nc * dnc ) / (nc*nc) ;

	//------------------------------------------------------
	// Calculation of N^X at the center of the star ---> nx
	//------------------------------------------------------
	nx = shift(0)(0, 0, 0, 0) ;

	//-----------------------------------------------------------
	// Calculation of dN^X/dX at the center of the star ---> dnx
	//-----------------------------------------------------------
	dnx = factx * shift(0).dsdx()(0, 0, 0, 0) ;

	//------------------------------------------------------
	// Calculation of N^Y at the center of the star ---> ny
	//------------------------------------------------------
	ny = shift(1)(0, 0, 0, 0) ;

	//-----------------------------------------------------------
	// Calculation of dN^Y/dX at the center of the star ---> dny
	//-----------------------------------------------------------
	dny = factx * shift(1).dsdx()(0, 0, 0, 0) ;

	//--------------------------------------------
	// Calculation of (N^X)^2 + (N^Y)^2 + (N^Z)^2
	//  at the center of the star ---> npn
	//--------------------------------------------
	tmp = flat_scalar_prod(shift, shift)() ; // For the flat part
	npn = tmp(0, 0, 0, 0) ;

	//----------------------------------------------------
	// Calculation of d/dX( (N^X)^2 + (N^Y)^2 + (N^Z)^2 )
	//  at the center of the star ---> dnpn
	//----------------------------------------------------
	dnpn = factx * tmp.dsdx()(0, 0, 0, 0) ;

    }  // Finish of the relativistic case
    else {
	cout << "Bin_bhns_extr::orbit_omega : "
	     << "It should be the relativistic calculation !" << endl ;
	abort() ;
    }

    cout << "Bin_bhns_extr::orbit_omega: coord. ori. d(nu+log(Gam))/dX : "
	 << dnulg << endl ;
    cout << "Bin_bhns_extr::orbit_omega: coord. ori. A^2/N^2 : "
	 << asn2 << endl ;
    cout << "Bin_bhns_extr::orbit_omega: coord. ori. d(A^2/N^2)/dX : "
	 << dasn2 << endl ; 
    cout << "Bin_bhns_extr::orbit_omega: coord. ori. N^X : " << nx << endl ;
    cout << "Bin_bhns_extr::orbit_omega: coord. ori. dN^X/dX : "
	 << dnx << endl ;
    cout << "Bin_bhns_extr::orbit_omega: coord. ori. N^Y : " << ny << endl ;
    cout << "Bin_bhns_extr::orbit_omega: coord. ori. dN^Y/dX : "
	 << dny << endl ;
    cout << "Bin_bhns_extr::orbit_omega: coord. ori. N.N : " << npn << endl ;
    cout << "Bin_bhns_extr::orbit_omega: coord. ori. d(N.N)/dX : "
	 << dnpn << endl ; 

    //------------------------------------------------------
    // Start of calculation of the orbital angular velocity
    //------------------------------------------------------
    int relat = ( star.is_relativistic() ) ? 1 : 0 ;

    Param parf ;
    parf.add_int(relat) ;
    parf.add_double( (star.get_mp()).get_ori_x(), 0) ;
    parf.add_double( dnulg, 1) ;
    parf.add_double( asn2, 2) ;
    parf.add_double( dasn2, 3) ;
    parf.add_double( nx, 4) ;
    parf.add_double( dnx, 5) ;
    parf.add_double( ny, 6) ;
    parf.add_double( dny, 7) ;
    parf.add_double( npn, 8) ;
    parf.add_double( dnpn, 9) ;
    parf.add_double( msr, 10) ;

    double omega1 = fact_omeg_min * omega ;
    double omega2 = fact_omeg_max * omega ;

    cout << "Bin_bhns_extr::orbit_omega: omega1, omega2 [rad/s] : "
	 << omega1 * f_unit << "  " << omega2 * f_unit << endl ;

    // Search for the various zeros in the interval [omega1,omega2]
    // ------------------------------------------------------------
    int nsub = 50 ;  // total number of subdivisions of the interval
    Tbl* azer = 0x0 ;
    Tbl* bzer = 0x0 ;
    zero_list(fonc_bhns_orbit_ks, parf, omega1, omega2, nsub,
	      azer, bzer) ;

    // Search for the zero closest to the previous value of omega
    // ----------------------------------------------------------
    double omeg_min, omeg_max ;
    int nzer = azer->get_taille() ; // number of zeros found by zero_list
    cout << "Bin_bhns_extr::orbit_omega : " << nzer <<
      "zero(s) found in the interval [omega1,  omega2]." << endl ;
    cout << "omega, omega1, omega2 : " << omega << "  " << omega1
	 << "  " << omega2 << endl ; 
    cout << "azer : " << *azer << endl ;
    cout << "bzer : " << *bzer << endl ;

    if (nzer == 0) {
	cout << "Bin_bhns_extr::orbit_omega: WARNING : "
	     << "no zero detected in the interval" << endl
	     << "   [" << omega1 * f_unit << ", " 
	     << omega2 * f_unit << "]  rad/s  !" << endl ;
	omeg_min = omega1 ;
	omeg_max = omega2 ;
    }
    else {
	double dist_min = fabs(omega2 - omega1) ;
	int i_dist_min = -1 ;
	for (int i=0; i<nzer; i++) {
	    // Distance of previous value of omega from the center of the
	    //  interval [azer(i), bzer(i)]
	    double dist = fabs( omega - 0.5 * ( (*azer)(i) + (*bzer)(i) ) ) ;

	    if (dist < dist_min) {
		dist_min = dist ;
		i_dist_min = i ;
	    }
	}
	omeg_min = (*azer)(i_dist_min) ;
	omeg_max = (*bzer)(i_dist_min) ;
    }

    delete azer ; // Tbl allocated by zero_list
    delete bzer ; //

    cout << "Bin_bhns_extr:orbit_omega : "
	 << "interval selected for the search of the zero : "
	 << endl << "  [" << omeg_min << ", " << omeg_max << "] = ["
	 << omeg_min * f_unit << ", " << omeg_max * f_unit << "] rad/s "
	 << endl ;

    // Computation of the zero in the selected interval by the secant method
    // ---------------------------------------------------------------------

    int nitermax = 200 ;
    int niter ;
    double precis = 1.e-13 ;
    omega = zerosec_b(fonc_bhns_orbit_ks, parf, omeg_min, omeg_max,
		      precis, nitermax, niter) ;

    cout << "Bin_bhns_extr::orbit_omega : "
	 << "Number of iterations in zerosec for omega : "
	 << niter << endl ;

    cout << "Bin_bhns_extr::orbit_omega : omega [rad/s] : "
	 << omega * f_unit << endl ;

}


void Bin_bhns_extr::orbit_omega_cf(double fact_omeg_min,
				   double fact_omeg_max) {

  using namespace Unites ;

    //------------------------------------------------------------------
    // Evaluation of various quantities at the center of a neutron star
    //------------------------------------------------------------------

    double dnulg, asn2, dasn2, ny, dny, npn, dnpn ;

    const Map& mp = star.get_mp() ;

    const Cmp& logn_auto_regu = star.get_logn_auto_regu()() ;
    const Cmp& logn_comp = star.get_logn_comp()() ;
    const Cmp& loggam = star.get_loggam()() ;
    const Cmp& nnn = star.get_nnn()() ;
    const Cmp& a_car = star.get_a_car()() ;
    const Tenseur& shift = star.get_shift() ;
    const Tenseur& d_logn_auto_div = star.get_d_logn_auto_div() ;

    Tenseur dln_auto_div = d_logn_auto_div ;

    if ( *(dln_auto_div.get_triad()) != ref_triad ) {

        // Change the basis from spherical coordinate to Cartesian one
        dln_auto_div.change_triad( mp.get_bvect_cart() ) ;

	// Change the basis from mapping coordinate to absolute one
	dln_auto_div.change_triad( ref_triad ) ;

    }

    //-----------------------------------------------------------------
    // Calculation of d/dX( nu + ln(Gamma) ) at the center of the star
    //  ---> dnulg
    //-----------------------------------------------------------------

    // Factor for the coordinate transformation x --> X :
    double factx ;
    if (fabs(mp.get_rot_phi()) < 1.e-14) {
	factx = 1. ;
    }
    else {
	if (fabs(mp.get_rot_phi() - M_PI) < 1.e-14) {
	    factx = - 1. ;
	}
	else {
	    cout <<
	      "Bin_bhns_extr::orbit_omega_cfc : unknown value of rot_phi !"
		 << endl ;
	    abort() ;
	}
    }

    Cmp tmp = logn_auto_regu + logn_comp + loggam ;

    // ... gradient
    dnulg = dln_auto_div(0)(0, 0, 0, 0)
      + factx * tmp.dsdx()(0, 0, 0, 0) ;

    //------------------------------------------------------------
    // Calculation of A^2/N^2 at the center of the star ---> asn2
    //------------------------------------------------------------
    double nc = nnn(0, 0, 0, 0) ;
    double a2c = a_car(0, 0, 0, 0) ;
    asn2 = a2c / (nc * nc) ;

    if ( star.is_relativistic() ) {

        //------------------------------------------------------------------
	// Calculation of d/dX(A^2/N^2) at the center of the star ---> dasn
	//------------------------------------------------------------------

        double da2c = factx * a_car.dsdx()(0, 0, 0, 0) ; 
	double dnc =  factx * nnn.dsdx()(0, 0, 0, 0) ;

	dasn2 = ( da2c - 2. * a2c / nc * dnc ) / (nc*nc) ;

	//------------------------------------------------------
	// Calculation of N^Y at the center of the star ---> ny
	//------------------------------------------------------
	ny = shift(1)(0, 0, 0, 0) ;

	//-----------------------------------------------------------
	// Calculation of dN^Y/dX at the center of the star ---> dny
	//-----------------------------------------------------------
	dny = factx * shift(1).dsdx()(0, 0, 0, 0) ;

	//--------------------------------------------
	// Calculation of (N^X)^2 + (N^Y)^2 + (N^Z)^2
	//  at the center of the star ---> npn
	//--------------------------------------------
	tmp = flat_scalar_prod(shift, shift)() ; // For the flat part
	npn = tmp(0, 0, 0, 0) ;

	//----------------------------------------------------
	// Calculation of d/dX( (N^X)^2 + (N^Y)^2 + (N^Z)^2 )
	//  at the center of the star ---> dnpn
	//----------------------------------------------------
	dnpn = factx * tmp.dsdx()(0, 0, 0, 0) ;

    }  // Finish of the relativistic case
    else {
	cout << "Bin_bhns_extr::orbit_omega_cfc : "
	     << "It should be the relativistic calculation !" << endl ;
	abort() ;
    }

    cout << "Bin_bhns_extr::orbit_omega_cfc: coord. ori. d(nu+log(Gam))/dX : "
	 << dnulg << endl ;
    cout << "Bin_bhns_extr::orbit_omega_cfc: coord. ori. A^2/N^2 : "
	 << asn2 << endl ;
    cout << "Bin_bhns_extr::orbit_omega_cfc: coord. ori. d(A^2/N^2)/dX : "
	 << dasn2 << endl ; 
    cout << "Bin_bhns_extr::orbit_omega_cfc: coord. ori. N^Y : "
	 << ny << endl ;
    cout << "Bin_bhns_extr::orbit_omega_cfc: coord. ori. dN^Y/dX : "
	 << dny << endl ;
    cout << "Bin_bhns_extr::orbit_omega_cfc: coord. ori. N.N : "
	 << npn << endl ;
    cout << "Bin_bhns_extr::orbit_omega_cfc: coord. ori. d(N.N)/dX : "
	 << dnpn << endl ; 

    //------------------------------------------------------
    // Start of calculation of the orbital angular velocity
    //------------------------------------------------------
    int relat = ( star.is_relativistic() ) ? 1 : 0 ;

    Param parf ;
    parf.add_int(relat) ;
    parf.add_double( (star.get_mp()).get_ori_x(), 0) ;
    parf.add_double( dnulg, 1) ;
    parf.add_double( asn2, 2) ;
    parf.add_double( dasn2, 3) ;
    parf.add_double( ny, 4) ;
    parf.add_double( dny, 5) ;
    parf.add_double( npn, 6) ;
    parf.add_double( dnpn, 7) ;

    double omega1 = fact_omeg_min * omega ;
    double omega2 = fact_omeg_max * omega ;

    cout << "Bin_bhns_extr::orbit_omega_cfc: omega1, omega2 [rad/s] : "
	 << omega1 * f_unit << "  " << omega2 * f_unit << endl ;

    // Search for the various zeros in the interval [omega1,omega2]
    // ------------------------------------------------------------
    int nsub = 50 ;  // total number of subdivisions of the interval
    Tbl* azer = 0x0 ;
    Tbl* bzer = 0x0 ;
    zero_list(fonc_bhns_orbit_cf, parf, omega1, omega2, nsub,
	      azer, bzer) ;

    // Search for the zero closest to the previous value of omega
    // ----------------------------------------------------------
    double omeg_min, omeg_max ;
    int nzer = azer->get_taille() ; // number of zeros found by zero_list
    cout << "Bin_bhns_extr::orbit_omega_cfc : " << nzer <<
      "zero(s) found in the interval [omega1,  omega2]." << endl ;
    cout << "omega, omega1, omega2 : " << omega << "  " << omega1
	 << "  " << omega2 << endl ; 
    cout << "azer : " << *azer << endl ;
    cout << "bzer : " << *bzer << endl ;

    if (nzer == 0) {
	cout << "Bin_bhns_extr::orbit_omega_cfc: WARNING : "
	     << "no zero detected in the interval" << endl
	     << "   [" << omega1 * f_unit << ", " 
	     << omega2 * f_unit << "]  rad/s  !" << endl ;
	omeg_min = omega1 ;
	omeg_max = omega2 ;
    }
    else {
	double dist_min = fabs(omega2 - omega1) ;
	int i_dist_min = -1 ;
	for (int i=0; i<nzer; i++) {
	    // Distance of previous value of omega from the center of the
	    //  interval [azer(i), bzer(i)]
	    double dist = fabs( omega - 0.5 * ( (*azer)(i) + (*bzer)(i) ) ) ;

	    if (dist < dist_min) {
		dist_min = dist ;
		i_dist_min = i ;
	    }
	}
	omeg_min = (*azer)(i_dist_min) ;
	omeg_max = (*bzer)(i_dist_min) ;
    }

    delete azer ; // Tbl allocated by zero_list
    delete bzer ; //

    cout << "Bin_bhns_extr:orbit_omega_cfc : "
	 << "interval selected for the search of the zero : "
	 << endl << "  [" << omeg_min << ", " << omeg_max << "] = ["
	 << omeg_min * f_unit << ", " << omeg_max * f_unit << "] rad/s "
	 << endl ;

    // Computation of the zero in the selected interval by the secant method
    // ---------------------------------------------------------------------

    int nitermax = 200 ;
    int niter ;
    double precis = 1.e-13 ;
    omega = zerosec_b(fonc_bhns_orbit_cf, parf, omeg_min, omeg_max,
		      precis, nitermax, niter) ;

    cout << "Bin_bhns_extr::orbit_omega_cfc : "
	 << "Number of iterations in zerosec for omega : "
	 << niter << endl ;

    cout << "Bin_bhns_extr::orbit_omega_cfc : omega [rad/s] : "
	 << omega * f_unit << endl ;

}


//***********************************************************
//  Function used for search of the orbital angular velocity
//***********************************************************

double fonc_bhns_orbit_ks(double om, const Param& parf) {

    int relat = parf.get_int() ;

    double xc = parf.get_double(0) ;
    double dnulg = parf.get_double(1) ;
    double asn2 = parf.get_double(2) ;
    double dasn2 = parf.get_double(3) ;
    double nx = parf.get_double(4) ;
    double dnx = parf.get_double(5) ;
    double ny = parf.get_double(6) ;
    double dny = parf.get_double(7) ;
    double npn = parf.get_double(8) ;
    double dnpn = parf.get_double(9) ;
    double msr = parf.get_double(10) ;

    double om2 = om*om ;

    double dphi_cent ;

    if (relat == 1) {
	double bpb = om2 * xc*xc - 2.*om * ny * xc + npn + 2.*msr * nx*nx;

	dphi_cent = ( asn2* ( om* (ny + xc*dny) - om2*xc - 0.5*dnpn 
			      - 2.*msr * nx * dnx + msr * nx*nx / xc )
		      - 0.5*bpb* dasn2 ) / ( 1 - asn2 * bpb ) ;
    }
    else {
	cout << "Bin_bhns_extr::orbit_omega_ks : "
	     << "It should be the relativistic calculation !" << endl ;
	abort() ;
    }

    return dnulg + dphi_cent ;

}

double fonc_bhns_orbit_cf(double om, const Param& parf) {

    int relat = parf.get_int() ;

    double xc = parf.get_double(0) ;
    double dnulg = parf.get_double(1) ;
    double asn2 = parf.get_double(2) ;
    double dasn2 = parf.get_double(3) ;
    double ny = parf.get_double(4) ;
    double dny = parf.get_double(5) ;
    double npn = parf.get_double(6) ;
    double dnpn = parf.get_double(7) ;

    double om2 = om*om ;

    double dphi_cent ;

    if (relat == 1) {
	double bpb = om2 * xc*xc - 2.*om * ny * xc + npn ;

	dphi_cent = ( asn2* ( om* (ny + xc*dny) - om2*xc - 0.5*dnpn )
		      - 0.5*bpb* dasn2 ) / ( 1 - asn2 * bpb ) ;
    }
    else {
	cout << "Bin_bhns_extr::orbit_omega_cf : "
	     << "It should be the relativistic calculation !" << endl ;
	abort() ;
    }

    return dnulg + dphi_cent ;

}
