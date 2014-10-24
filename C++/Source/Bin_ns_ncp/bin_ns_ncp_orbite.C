/*
 * Method of class Bin_ns_ncp to compute the orbital angular velocity {\tt omega}
 * and the position of the rotation axis {\tt x_axe}.
 *
 * (See file bin_ns_ncp.h for documentation)
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


char bin_ns_ncp_orbite_C[] = "$Header$" ;

/*
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene 
#include "bin_ns_ncp.h"
#include "eos.h"
#include "param.h"
#include "utilitaires.h"
#include "unites.h"

namespace Lorene {
double  fonc_bin_ncp_axe(double , const Param& ) ;
double  fonc_bin_ncp_orbit(double , const Param& ) ;

//******************************************************************************

void Bin_ns_ncp::orbit(double fact_omeg_min, double fact_omeg_max, double& xgg1, 
		     double& xgg2) {
  using namespace Unites ;
    
    //-------------------------------------------------------------
    // Evaluation of various quantities at the center of each star
    //-------------------------------------------------------------


    double g00[2], g10[2], g20[2], g11[2], g21[2], g22[2], bx[2], by[2] ;
      
    double bz[2], d1sn2[2], unsn2[2] ;

    double dnulg[2], xgg[2], ori_x[2], dg00[2], dg10[2], dg20[2], dg11[2] ;
      
    double dg21[2], dg22[2], dbx[2], dby[2], dbz[2], dbymo[2] ;

    for (int i=0; i<2; i++) {
	
	const Cmp& logn_auto = et[i]->get_logn_auto()() ; 
	const Cmp& logn_comp = et[i]->get_logn_comp()() ; 
	const Cmp& loggam = et[i]->get_loggam()() ; 
	const Cmp& nnn = et[i]->get_nnn()() ; 
	const Tenseur& shift = et[i]->get_shift() ; 
	const Metrique& gtilde = et[i]->get_gtilde() ;
	const Tenseur& gamma = et[i]->get_gamma() ;
	const Tenseur a_car = et[i]->get_a_car() ;

	const Metrique met_gamma(pow(gamma, 1./3.)*gtilde.cov()) ;

	const Cmp& gg00 = met_gamma.cov()(0,0) ;
	const Cmp& gg10 = met_gamma.cov()(1,0) ;
	const Cmp& gg20 = met_gamma.cov()(2,0) ;
	const Cmp& gg11 = met_gamma.cov()(1,1) ;
	const Cmp& gg21 = met_gamma.cov()(2,1) ;
	const Cmp& gg22 = met_gamma.cov()(2,2) ;

	const Cmp& bbx = shift(0) ;
	const Cmp& bby = shift(1) ;
	const Cmp& bbz = shift(2) ;

	
	//----------------------------------
	// Calcul de d/dX( nu + ln(Gamma) ) au centre de l'etoile ---> dnulg[i]
	//----------------------------------

	Cmp tmp = logn_auto + logn_comp + loggam ;
	
	// ... gradient suivant X : 		
	dnulg[i] = tmp.dsdx()(0, 0, 0, 0) ; 

	//----------------------------------
	// Calcul de gij, lapse et shift au centre de l'etoile
	//----------------------------------

	g00[i] = gg00(0,0,0,0) ; 
	g10[i] = gg10(0,0,0,0) ; 
	g20[i] = gg20(0,0,0,0) ; 
	g11[i] = gg11(0,0,0,0) ; 
	g21[i] = gg21(0,0,0,0) ; 
	g22[i] = gg22(0,0,0,0) ; 

	bx[i] = bbx(0,0,0,0) ;
	by[i] = bby(0,0,0,0) ;
	bz[i] = bbz(0,0,0,0) ;

	unsn2[i] = 1/(nnn(0,0,0,0)*nnn(0,0,0,0)) ; 

	//----------------------------------
	// Calcul de d/dX(gij), d/dX(shift) au centre de l'etoile
	//----------------------------------
	
	dg00[i] = gg00.dsdx()(0,0,0,0) ;
	dg10[i] = gg10.dsdx()(0,0,0,0) ;
	dg20[i] = gg20.dsdx()(0,0,0,0) ;
	dg11[i] = gg11.dsdx()(0,0,0,0) ;
	dg21[i] = gg21.dsdx()(0,0,0,0) ;
	dg22[i] = gg22.dsdx()(0,0,0,0) ;
	
	dbx[i] = bbx.dsdx()(0,0,0,0) ;
	dby[i] = bby.dsdx()(0,0,0,0) ;
 	dbz[i] = bbz.dsdx()(0,0,0,0) ;

	dbymo[i] = bby.dsdx()(0,0,0,0) - omega ;


	d1sn2[i] = (1/(nnn*nnn)).dsdx()(0,0,0,0) ;


	cout << "Bin_ns_ncp::orbit: central d(nu+log(Gam))/dX : " 
	     << dnulg[i] << endl ; 
	cout << "Bin_ns_ncp::orbit: central g00 :" << g00[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central g10 :" << g10[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central g20 :" << g20[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central g11 :" << g11[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central g21 :" << g21[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central g22 :" << g22[i] << endl ;

	cout << "Bin_ns_ncp::orbit: central shift_x :" << bx[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central shift_y :" << by[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central shift_z :" << bz[i] << endl ;

	cout << "Bin_ns_ncp::orbit: central d/dX(g00) :" << dg00[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central d/dX(g10) :" << dg10[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central d/dX(g20) :" << dg20[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central d/dX(g11) :" << dg11[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central d/dX(g21) :" << dg21[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central d/dX(g22) :" << dg22[i] << endl ;

	cout << a_car()(0,0,0,0) << " " << a_car().dsdx()(0,0,0,0) << endl ;

	cout << "Bin_ns_ncp::orbit: central d/dX(shift_x) :" << dbx[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central d/dX(shift_y) :" <<dby[i] << endl ;
	cout << "Bin_ns_ncp::orbit: central d/dX(shift_z) :" << dbz[i] << endl ;

	//----------------------
	// Pour information seulement : 1/ calcul des positions des "centres de
	//				    de masse"
	//				2/ calcul de dH/dX en r=0
	//-----------------------

        ori_x[i] = (et[i]->get_mp()).get_ori_x() ;

	xgg[i] = (et[i]->xa_barycenter() - ori_x[i]) ;
		 
    } // fin de la boucle sur les etoiles 

    xgg1 = xgg[0] ;
    xgg2 = xgg[1] ;
    
//---------------------------------
//  Position de l'axe de rotation   
//---------------------------------

    int relat = ( et[0]->is_relativistic() ) ? 1 : 0 ;
    double ori_x1 = ori_x[0] ;
    double ori_x2 = ori_x[1] ;

    if ( et[0]->get_eos() == et[1]->get_eos() &&
	 et[0]->get_ent()()(0,0,0,0) == et[1]->get_ent()()(0,0,0,0) ) {

        x_axe = 0. ;

    }
    else {

	Param paraxe ;
	paraxe.add_int(relat) ;
	paraxe.add_double( ori_x1, 0) ;
	paraxe.add_double( ori_x2, 1) ;
	paraxe.add_double( dnulg[0], 2) ;
	paraxe.add_double( dnulg[1], 3) ;
	paraxe.add_double( g00[0], 4) ;
	paraxe.add_double( g00[1], 5) ;
	paraxe.add_double( g10[0], 6) ;
	paraxe.add_double( g10[1], 7) ;
	paraxe.add_double( g20[0], 8) ;
	paraxe.add_double( g20[1], 9) ;
	paraxe.add_double( g11[0], 10) ;
	paraxe.add_double( g11[1], 11) ;
	paraxe.add_double( g21[0], 12) ;
	paraxe.add_double( g21[1], 13) ;
	paraxe.add_double( g22[0], 14) ;
	paraxe.add_double( g22[1], 15) ;
	paraxe.add_double( bx[0], 16) ;
	paraxe.add_double( bx[1], 17) ;
	paraxe.add_double( by[0], 18) ;
	paraxe.add_double( by[1], 19) ;
	paraxe.add_double( bz[0], 20) ;
	paraxe.add_double( bz[1], 21) ;
	paraxe.add_double( dg00[0], 22) ;
	paraxe.add_double( dg00[1], 23) ;
	paraxe.add_double( dg10[0], 24) ;
	paraxe.add_double( dg10[1], 25) ;
	paraxe.add_double( dg20[0], 26) ;
	paraxe.add_double( dg20[1], 27) ;
	paraxe.add_double( dg11[0], 28) ;
	paraxe.add_double( dg11[1], 29) ;
	paraxe.add_double( dg21[0], 30) ;
	paraxe.add_double( dg21[1], 31) ;
	paraxe.add_double( dg22[0], 32) ;
	paraxe.add_double( dg22[1], 33) ;
	paraxe.add_double( dbx[0], 34) ;
	paraxe.add_double( dbx[1], 35) ;
	paraxe.add_double( dbz[0], 36) ;
	paraxe.add_double( dbz[1], 37) ;
	paraxe.add_double( dbymo[0], 38) ;
	paraxe.add_double( dbymo[1], 39) ;
	paraxe.add_double( d1sn2[0], 40) ;
	paraxe.add_double( d1sn2[1], 41) ;
	paraxe.add_double( unsn2[0], 42) ;
	paraxe.add_double( unsn2[1], 43) ;
	paraxe.add_double( omega, 44) ;

	int nitmax_axe = 200 ; 
	int nit_axe ; 
	double precis_axe = 1.e-13 ;

	x_axe = zerosec(fonc_bin_ncp_axe, paraxe, 0.9*ori_x1, 0.9*ori_x2,
			precis_axe, nitmax_axe, nit_axe) ;

	cout << "Bin_ns_ncp::orbit : Number of iterations in zerosec for x_axe : "
	     << nit_axe << endl ;
    }

    cout << "Bin_ns_ncp::orbit : x_axe [km] : " << x_axe / km << endl ; 

//-------------------------------------
//  Calcul de la vitesse orbitale    
//-------------------------------------

    Param parf ; 
    parf.add_int(relat) ; 
    parf.add_double( ori_x1, 0) ; 
    parf.add_double( dnulg[0], 1) ;
    parf.add_double( g00[0], 2) ;
    parf.add_double( g10[0], 3) ;
    parf.add_double( g20[0], 4) ;
    parf.add_double( g11[0], 5) ;
    parf.add_double( g21[0], 6) ;
    parf.add_double( g22[0], 7) ;
    parf.add_double( bx[0], 8) ;
    parf.add_double( by[0], 9) ;
    parf.add_double( bz[0], 10) ;
    parf.add_double( dg00[0], 11) ;
    parf.add_double( dg10[0], 12) ;
    parf.add_double( dg20[0], 13) ;
    parf.add_double( dg11[0], 14) ;
    parf.add_double( dg21[0], 15) ;
    parf.add_double( dg22[0], 16) ;
    parf.add_double( dbx[0], 17) ;
    parf.add_double( dbz[0], 18) ;
    parf.add_double( dby[0], 19) ;
    parf.add_double( d1sn2[0], 20) ;
    parf.add_double( unsn2[0], 21) ;
    parf.add_double( x_axe, 22) ;
 

    double omega1 = fact_omeg_min * omega  ; 
    double omega2 = fact_omeg_max * omega ; 
    cout << "Bin_ns_ncp::orbit: omega1,  omega2 [rad/s] : " 
	 << omega1 * f_unit << "  " << omega2 * f_unit << endl ; 


	// Search for the various zeros in the interval [omega1,omega2]
	// ------------------------------------------------------------
	int nsub = 50 ;  // total number of subdivisions of the interval
	Tbl* azer = 0x0 ;
	Tbl* bzer = 0x0 ; 
	zero_list(fonc_bin_ncp_orbit, parf, omega1, omega2, nsub,
		  azer, bzer) ; 
	
	// Search for the zero closest to the previous value of omega
	// ----------------------------------------------------------
	double omeg_min, omeg_max ; 
	int nzer = azer->get_taille() ; // number of zeros found by zero_list
	cout << "Bin_ns_ncp:orbit : " << nzer << 
	     " zero(s) found in the interval [omega1,  omega2]." << endl ; 
	cout << "omega, omega1, omega2 : " << omega << "  " << omega1
		<< "  " << omega2 << endl ; 
	cout << "azer : " << *azer << endl ;
	cout << "bzer : " << *bzer << endl ;
	
	if (nzer == 0) {
		cout << 
		"Bin_ns_ncp::orbit: WARNING : no zero detected in the interval"
		<< endl << "   [" << omega1 * f_unit << ", " 
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
    cout << "Bin_ns_ncp:orbit : interval selected for the search of the zero : "
	 << endl << "  [" << omeg_min << ", " << omeg_max << "]  =  [" 
	 << omeg_min * f_unit << ", " << omeg_max * f_unit << "] rad/s " << endl ; 
    
    // Computation of the zero in the selected interval by the secant method
    // ---------------------------------------------------------------------

    int nitermax = 200 ; 
    int niter ; 
    double precis = 1.e-13 ;
    omega = zerosec_b(fonc_bin_ncp_orbit, parf, omeg_min, omeg_max,
		    precis, nitermax, niter) ;
    
    cout << "Bin_ns_ncp::orbit : Number of iterations in zerosec for omega : "
	 << niter << endl ; 
	
    cout << "Bin_ns_ncp::orbit : omega [rad/s] : "
	 << omega * f_unit << endl ; 
          

}


//*************************************************
//  Function used for search of the rotation axis
//*************************************************

double  fonc_bin_ncp_axe(double x_rot, const Param& paraxe) {

    double ori_x1 = paraxe.get_double(0) ;
    double ori_x2 = paraxe.get_double(1) ;
    double dnulg_1 = paraxe.get_double(2) ;
    double dnulg_2 = paraxe.get_double(3) ;
    double g00_1 = paraxe.get_double(4) ;
    double g00_2 = paraxe.get_double(5) ;
    double g10_1 = paraxe.get_double(6) ;
    double g10_2 = paraxe.get_double(7) ;
    double g20_1 = paraxe.get_double(8) ;
    double g20_2 = paraxe.get_double(9) ;
    double g11_1 = paraxe.get_double(10) ;
    double g11_2 = paraxe.get_double(11) ;
    double g21_1 = paraxe.get_double(12) ;
    double g21_2 = paraxe.get_double(13) ;
    double g22_1 = paraxe.get_double(14) ;
    double g22_2 = paraxe.get_double(15) ;
    double bx_1 = paraxe.get_double(16) ;
    double bx_2 = paraxe.get_double(17) ;
    double by_1 = paraxe.get_double(18) ;
    double by_2 = paraxe.get_double(19) ;
    double bz_1 = paraxe.get_double(20) ;
    double bz_2 = paraxe.get_double(21) ;
    double dg00_1 = paraxe.get_double(22) ;
    double dg00_2 = paraxe.get_double(23) ;
    double dg10_1 = paraxe.get_double(24) ;
    double dg10_2 = paraxe.get_double(25) ;
    double dg20_1 = paraxe.get_double(26) ;
    double dg20_2 = paraxe.get_double(27) ;
    double dg11_1 = paraxe.get_double(28) ;
    double dg11_2 = paraxe.get_double(29) ;
    double dg21_1 = paraxe.get_double(30) ;
    double dg21_2 = paraxe.get_double(31) ;
    double dg22_1 = paraxe.get_double(32) ;
    double dg22_2 = paraxe.get_double(33) ;
    double dbx_1 = paraxe.get_double(34) ;
    double dbx_2 = paraxe.get_double(35) ;
    double dbz_1 = paraxe.get_double(36) ;
    double dbz_2 = paraxe.get_double(37) ;
    double dbymo_1 = paraxe.get_double(38) ;
    double dbymo_2 = paraxe.get_double(39) ;
    double d1sn2_1 = paraxe.get_double(40) ;
    double d1sn2_2 = paraxe.get_double(41) ;
    double unsn2_1 = paraxe.get_double(42) ;
    double unsn2_2 = paraxe.get_double(43) ;
    double omega = paraxe.get_double(44) ;

    double om2_star1 ;
    double om2_star2 ;

    double x1 = ori_x1 - x_rot ;
    double x2 = ori_x2 - x_rot ;

    double bymxo_1 = by_1-x1*omega ; 
    double bymxo_2 = by_2-x2*omega ;

    
    double beta1 = g00_1*bx_1*bx_1 + 2*g10_1*bx_1*bymxo_1 + 2*g20_1*bx_1*bz_1 ;
    double beta2 = g11_1*bymxo_1*bymxo_1 + 2*g21_1*bz_1*bymxo_1 
                   + g22_1*bz_1*bz_1 ;

    double beta_1 = beta1 + beta2 ;


    double delta1 = dg00_1*bx_1*bx_1 + 2*g00_1*dbx_1*bx_1 + 2*dg10_1*bx_1*bymxo_1 ;
    double delta2 = 2*g10_1*bymxo_1*dbx_1 + 2*g10_1*bx_1*dbymo_1 + 2*dg20_1*bx_1*bz_1 ;
    double delta3 = 2*g20_1*bx_1*dbz_1 +2*g20_1*bz_1*dbx_1 + dg11_1*bymxo_1*bymxo_1 ;
    double delta4 = 2*g11_1*bymxo_1*dbymo_1 + 2*dg21_1*bz_1*bymxo_1;
    double delta5 = 2*g21_1*bymxo_1*dbz_1 +2*g21_1*bz_1*dbymo_1 + dg22_1*bz_1*bz_1 + 2*g22_1*bz_1*dbz_1 ;

    double delta_1 = delta1 + delta2 + delta3 + delta4 + delta5 ;

    // Computation of omega for star 1
    //---------------------------------

    om2_star1 = dnulg_1 / (beta_1/(omega*omega)*(dnulg_1*unsn2_1 + d1sn2_1/2.) 
			   + unsn2_1*delta_1/(omega*omega)/2.) ;



    double beta3 = g00_2*bx_2*bx_2 + 2*g10_2*bx_2*bymxo_2 + 2*g20_2*bx_2*bz_2 ;
    double beta4 = g11_2*bymxo_2*bymxo_2 + 2*g21_2*bz_2*bymxo_2
                   + g22_2*bz_2*bz_2 ;

    double beta_2 = beta3 + beta4 ;

       
    double delta6 = dg00_2*bx_2*bx_2 + 2*g00_2*dbx_2*bx_2 + 2*dg10_2*bx_2*bymxo_2 ;
    double delta7 = 2*g10_2*bymxo_2*dbx_2 + 2*g10_2*bx_2*dbymo_2 + 2*dg20_2*bx_2*bz_2 ;
    double delta8 = 2*g20_2*bx_2*dbz_2 +2*g20_2*bz_2*dbx_2 + dg11_2*bymxo_2*bymxo_2 ;
    double delta9 = 2*g11_2*bymxo_2*dbymo_2 + 2*dg21_2*bz_2*bymxo_2;
    double delta10 = 2*g21_2*bymxo_2*dbz_2 +2*g21_2*bz_2*dbymo_2 + dg22_2*bz_2*bz_2 + 2*g22_2*bz_2*dbz_2 ;

    double delta_2 = delta6 + delta7 + delta8 + delta9 + delta10 ;

    // Computation of omega for star 2
    //---------------------------------

    om2_star2 = dnulg_2 / (beta_2/(omega*omega)*(dnulg_2*unsn2_2 + d1sn2_2/2.) 
			   + unsn2_2*delta_2/(omega*omega)/2.) ;
                                                                            ; 
  
    return om2_star1 - om2_star2 ;

}

//*****************************************************************************
//  Fonction utilisee pour la recherche de omega par la methode de la secante
//*****************************************************************************

double fonc_bin_ncp_orbit(double om, const Param& parf) {

    double xc = parf.get_double(0) ; 
    double dnulg = parf.get_double(1) ;
    double g00 = parf.get_double(2) ;
    double g10 = parf.get_double(3) ;
    double g20 = parf.get_double(4) ;
    double g11 = parf.get_double(5) ;
    double g21 = parf.get_double(6) ;
    double g22 = parf.get_double(7) ;
    double bx = parf.get_double(8) ;
    double by = parf.get_double(9) ;
    double bz = parf.get_double(10) ;
    double dg00 = parf.get_double(11) ;
    double dg10 = parf.get_double(12) ;
    double dg20 = parf.get_double(13) ;
    double dg11 = parf.get_double(14) ;
    double dg21 = parf.get_double(15) ;
    double dg22 = parf.get_double(16) ;
    double dbx = parf.get_double(17) ;
    double dbz = parf.get_double(18) ;
    double dby = parf.get_double(19) ;
    double d1sn2 = parf.get_double(20) ;
    double unsn2 = parf.get_double(21) ;
    double x_axe = parf.get_double(22) ;
 

    double dbymo = dby - om ;
    double xx = xc - x_axe ; 
    
    double bymxo = by-xx*om ;  


    double beta1 = g00*bx*bx + 2*g10*bx*bymxo + 2*g20*bx*bz ;
    double beta2 = g11*bymxo*bymxo + 2*g21*bz*bymxo+ g22*bz*bz ;
    double beta = beta1 + beta2 ;
       
    double alpha = 1 - unsn2*beta ;

    double delta1 = dg00*bx*bx + 2*g00*dbx*bx + 2*dg10*bx*bymxo ;
    double delta2 = 2*g10*bymxo*dbx + 2*g10*bx*dbymo + 2*dg20*bx*bz ;
    double delta3 = 2*g20*bx*dbz +2*g20*bz*dbx + dg11*bymxo*bymxo ;
    double delta4 = 2*g11*bymxo*dbymo + 2*dg21*bz*bymxo;
    double delta5 = 2*g21*bymxo*dbz +2*g21*bz*dbymo + dg22*bz*bz + 2*g22*bz*dbz ;

    double delta = delta1 + delta2 + delta3 + delta4 + delta5 ;

    // Difference entre les 2 termes de l'eq.(95) de Gourgoulhon et.al (2001)  
    //centre de l'etoile 
    //-----------------------------------------------------------------------

     double diff = dnulg + (1/(2.*alpha))*(-d1sn2*beta - unsn2*delta) ; 
 
     return diff ; 
       
}



}
