 /*
 * Method of class Binaire to compute the orbital angular velocity {\tt omega}
 * and the position of the rotation axis {\tt x_axe}.
 *
 * (See file binaire.h for documentation)
 *
 */

/*
 *   Copyright (c) 2000-2003 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Keisuke Taniguchi
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


char binaire_orbite_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2014/10/13 08:52:44  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:12:59  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2011/03/27 16:42:21  e_gourgoulhon
 * Added output via new function save_profile for graphics.
 *
 * Revision 1.4  2009/06/18 18:40:57  k_taniguchi
 * Added a slightly modified code to determine
 * the orbital angular velocity.
 *
 * Revision 1.3  2004/03/25 10:28:59  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.2  2003/09/16 13:41:27  e_gourgoulhon
 * Search of sub-intervals containing the zero(s) via zero_list
 * Selection of the sub-interval as the closest one to previous value of omega.
 * Replaced zerosec by zerosec_b.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.9  2001/08/14  12:36:45  keisuke
 * Change of the minimum and maximum values for searching a new position
 *  of the rotation axis.
 *
 * Revision 2.8  2001/03/20  16:06:55  keisuke
 * Addition of the method directly to calculate
 *  the orbital angular velocity (we do not use this method!).
 *
 * Revision 2.7  2001/03/19  17:10:28  keisuke
 * Change of the definition of the canter of mass.
 *
 * Revision 2.6  2001/03/19  13:37:04  keisuke
 * Set x_axe to be zero for the case of identical stars.
 *
 * Revision 2.5  2001/03/17  16:17:20  keisuke
 * Input a subroutine to determine the position of the rotation axis
 *  in the case of binary systems composed of different stars.
 *
 * Revision 2.4  2000/10/02  08:29:34  keisuke
 * Change the method to construct d_logn_auto in the calculation of dnulg[i].
 *
 * Revision 2.3  2000/09/22  15:56:32  keisuke
 * *** empty log message ***
 *
 * Revision 2.2  2000/09/22  15:52:12  keisuke
 * Changement calcul de dlogn (prise en compte d'une partie divergente
 * dans logn_auto par l'appel a d_logn_auto).
 *
 * Revision 2.1  2000/02/12  18:36:09  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/02/12  17:09:21  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene 
#include "binaire.h"
#include "eos.h"
#include "param.h"
#include "utilitaires.h"
#include "unites.h"

#include "scalar.h"
#include "graphique.h"

namespace Lorene {
double  fonc_binaire_axe(double , const Param& ) ;
double  fonc_binaire_orbit(double , const Param& ) ;

//******************************************************************************

void Binaire::orbit(double fact_omeg_min, double fact_omeg_max, double& xgg1, 
		     double& xgg2) {

  using namespace Unites ;
    
    //-------------------------------------------------------------
    // Evaluation of various quantities at the center of each star
    //-------------------------------------------------------------

    double dnulg[2], asn2[2], dasn2[2], ny[2], dny[2], npn[2], dnpn[2], xgg[2] ;     
    double nyso[2], dnyso[2], npnso2[2], dnpnso2[2], ori_x[2] ;

    for (int i=0; i<2; i++) {
	
	const Map& mp = et[i]->get_mp() ; 

	const Cmp& logn_auto_regu = et[i]->get_logn_auto_regu()() ; 
	const Cmp& logn_comp = et[i]->get_logn_comp()() ; 
	const Cmp& loggam = et[i]->get_loggam()() ; 
	const Cmp& nnn = et[i]->get_nnn()() ; 
	const Cmp& a_car = et[i]->get_a_car()() ; 
	const Tenseur& shift = et[i]->get_shift() ; 
	const Tenseur& d_logn_auto_div = et[i]->get_d_logn_auto_div() ; 

	Tenseur dln_auto_div = d_logn_auto_div ;

	if ( *(dln_auto_div.get_triad()) != ref_triad ) {

	  // Change the basis from spherical coordinate to Cartesian one
	  dln_auto_div.change_triad( mp.get_bvect_cart() ) ;

	  // Change the basis from mapping coordinate to absolute one
	  dln_auto_div.change_triad( ref_triad ) ;

	}

	//----------------------------------
	// Calcul de d/dX( nu + ln(Gamma) ) au centre de l'etoile ---> dnulg[i]
	//----------------------------------
	
	// Facteur de passage x --> X :
	double factx ;
	if (fabs(mp.get_rot_phi()) < 1.e-14) {
	    factx = 1. ; 
	}
	else {
	    if (fabs(mp.get_rot_phi() - M_PI) < 1.e-14) {
		factx = - 1. ; 
	    }
	    else {
		cout << "Binaire::orbit : unknown value of rot_phi !" << endl ;
		abort() ; 
	    }
	}
	    
	Cmp tmp = logn_auto_regu + logn_comp + loggam ;
	
	// ... gradient suivant X : 		
	dnulg[i] = dln_auto_div(0)(0, 0, 0, 0)
	  + factx * tmp.dsdx()(0, 0, 0, 0) ; 
	
	tmp = logn_auto_regu + logn_comp ; 
	cout << "dlnndx_div_c : " <<  dln_auto_div(0)(0, 0, 0, 0)  << endl ; 
	cout << "dlnndx_c : " <<  dln_auto_div(0)(0, 0, 0, 0) + factx*tmp.dsdx()(0, 0, 0, 0) << endl ; 
				  
	cout << "dloggamdx_c : " <<  factx*loggam.dsdx()(0, 0, 0, 0) << endl ; 
	Scalar stmp(logn_comp) ; 
        save_profile(stmp, 0., 10., 0.5*M_PI, 0., "prof_logn.d") ; 
	stmp = loggam ; 
        save_profile(stmp, 0., 1.8, 0.5*M_PI, 0., "prof_loggam.d") ; 
  
	//----------------------------------
	// Calcul de A^2/N^2 au centre de l'etoile ---> asn2[i]
	//----------------------------------

	double nc = nnn(0, 0, 0, 0) ;
	double a2c = a_car(0, 0, 0, 0) ;
	asn2[i] = a2c / (nc * nc) ;
	 
	if ( et[i]->is_relativistic() ) {

	    //----------------------------------
	    // Calcul de d/dX(A^2/N^2) au centre de l'etoile ---> dasn2[i]
	    //----------------------------------

	    double da2c = factx * a_car.dsdx()(0, 0, 0, 0) ; 
	    double dnc =  factx * nnn.dsdx()(0, 0, 0, 0) ;

	    dasn2[i] = ( da2c - 2 * a2c / nc * dnc ) / (nc*nc) ; 

	    //----------------------------------
	    // Calcul de N^Y au centre de l'etoile ---> ny[i]
	    //----------------------------------

	    ny[i] = shift(1)(0, 0, 0, 0) ; 
	    nyso[i] = ny[i] / omega ;
	    
	    //----------------------------------
	    // Calcul de dN^Y/dX au centre de l'etoile ---> dny[i]
	    //----------------------------------
	    
	    dny[i] = factx * shift(1).dsdx()(0, 0, 0, 0) ; 
	    dnyso[i] = dny[i] / omega ;

	    //----------------------------------
	    // Calcul de (N^X)^2 + (N^Y)^2 + (N^Z)^2 
	    //				     au centre de l'etoile ---> npn[i]
	    //----------------------------------

	    tmp = flat_scalar_prod(shift, shift)() ; 

	    npn[i] = tmp(0, 0, 0, 0) ; 
	    npnso2[i] = npn[i] / omega / omega ;

	    //----------------------------------
	    // Calcul de d/dX( (N^X)^2 + (N^Y)^2 + (N^Z)^2 )
	    //				     au centre de l'etoile ---> dnpn[i]
	    //----------------------------------
	    
	    dnpn[i] = factx * tmp.dsdx()(0, 0, 0, 0) ; 
	    dnpnso2[i] = dnpn[i] / omega / omega ;

	}	    // fin du cas relativiste 
	else {
	    dasn2[i] = 0 ; 
	    ny[i] = 0 ; 
	    nyso[i] = 0 ;
	    dny[i] = 0 ; 
	    dnyso[i] = 0 ;
	    npn[i] = 0 ; 
	    npnso2[i] = 0 ;
	    dnpn[i] = 0 ; 
	    dnpnso2[i] = 0 ;
	}

	cout << "Binaire::orbit: central d(nu+log(Gam))/dX : " 
	     << dnulg[i] << endl ; 
	cout << "Binaire::orbit: central A^2/N^2 : " << asn2[i] << endl ; 
	cout << "Binaire::orbit: central d(A^2/N^2)/dX : " << dasn2[i] << endl ; 
	cout << "Binaire::orbit: central N^Y : " << ny[i] << endl ; 
	cout << "Binaire::orbit: central dN^Y/dX : " << dny[i] << endl ; 
	cout << "Binaire::orbit: central N.N : " << npn[i] << endl ; 
	cout << "Binaire::orbit: central d(N.N)/dX : " << dnpn[i] << endl ; 

	//----------------------
	// Pour information seulement : 1/ calcul des positions des "centres de
	//				    de masse"
	//				2/ calcul de dH/dX en r=0
	//-----------------------

        ori_x[i] = (et[i]->get_mp()).get_ori_x() ;

	xgg[i] = factx * (et[i]->xa_barycenter() - ori_x[i]) ;
		 
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
	paraxe.add_double( asn2[0], 4) ;
	paraxe.add_double( asn2[1], 5) ;
	paraxe.add_double( dasn2[0], 6) ;
	paraxe.add_double( dasn2[1], 7) ;
	paraxe.add_double( nyso[0], 8) ;
	paraxe.add_double( nyso[1], 9) ;
	paraxe.add_double( dnyso[0], 10) ;
	paraxe.add_double( dnyso[1], 11) ;
	paraxe.add_double( npnso2[0], 12) ;
	paraxe.add_double( npnso2[1], 13) ;
	paraxe.add_double( dnpnso2[0], 14) ;
	paraxe.add_double( dnpnso2[1], 15) ;

	int nitmax_axe = 200 ; 
	int nit_axe ; 
	double precis_axe = 1.e-13 ;

	x_axe = zerosec(fonc_binaire_axe, paraxe, 0.9*ori_x1, 0.9*ori_x2,
			precis_axe, nitmax_axe, nit_axe) ;

	cout << "Binaire::orbit : Number of iterations in zerosec for x_axe : "
	     << nit_axe << endl ;
    }

    cout << "Binaire::orbit : x_axe [km] : " << x_axe / km << endl ; 

//-------------------------------------
//  Calcul de la vitesse orbitale    
//-------------------------------------

    Param parf ; 
    parf.add_int(relat) ; 
    parf.add_double( (et[0]->get_mp()).get_ori_x(), 0) ; 
    parf.add_double( dnulg[0], 1) ;  
    parf.add_double( asn2[0], 2) ;    
    parf.add_double( dasn2[0], 3) ;    
    parf.add_double( ny[0], 4) ;    
    parf.add_double( dny[0], 5) ;    
    parf.add_double( npn[0], 6) ;    
    parf.add_double( dnpn[0], 7) ;    
    parf.add_double( x_axe, 8) ; 

    double omega1 = fact_omeg_min * omega  ; 
    double omega2 = fact_omeg_max * omega ; 
    cout << "Binaire::orbit: omega1,  omega2 [rad/s] : " 
	 << omega1 * f_unit << "  " << omega2 * f_unit << endl ; 

	// Search for the various zeros in the interval [omega1,omega2]
	// ------------------------------------------------------------
	int nsub = 50 ;  // total number of subdivisions of the interval
	Tbl* azer = 0x0 ;
	Tbl* bzer = 0x0 ; 
	zero_list(fonc_binaire_orbit, parf, omega1, omega2, nsub,
		  azer, bzer) ; 
	
	// Search for the zero closest to the previous value of omega
	// ----------------------------------------------------------
	double omeg_min, omeg_max ; 
	int nzer = azer->get_taille() ; // number of zeros found by zero_list
	cout << "Binaire:orbit : " << nzer << 
	     " zero(s) found in the interval [omega1,  omega2]." << endl ; 
	cout << "omega, omega1, omega2 : " << omega << "  " << omega1
		<< "  " << omega2 << endl ; 
	cout << "azer : " << *azer << endl ;
	cout << "bzer : " << *bzer << endl ;
	
	if (nzer == 0) {
		cout << 
		"Binaire::orbit: WARNING : no zero detected in the interval"
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
    
	cout << "Binaire:orbit : interval selected for the search of the zero : "
		<< endl << "  [" << omeg_min << ", " << omeg_max << "]  =  [" 
		 << omeg_min * f_unit << ", " << omeg_max * f_unit << "] rad/s " << endl ; 
	
	// Computation of the zero in the selected interval by the secant method
	// ---------------------------------------------------------------------
    int nitermax = 200 ; 
    int niter ; 
    double precis = 1.e-13 ;
    omega = zerosec_b(fonc_binaire_orbit, parf, omeg_min, omeg_max,
		    precis, nitermax, niter) ;
    
    cout << "Binaire::orbit : Number of iterations in zerosec for omega : "
	 << niter << endl ; 
	
    cout << "Binaire::orbit : omega [rad/s] : "
	 << omega * f_unit << endl ; 
          

    /* We do not use the method below.
    // Direct calculation
    //--------------------

    double om2_et1 ;
    double om2_et2 ;

    double x_et1 = ori_x1 - x_axe ;
    double x_et2 = ori_x2 - x_axe ;

    if (relat == 1) {

      double andan_et1 = 0.5 * dasn2[0] + asn2[0] * dnulg[0] ;
      double andan_et2 = 0.5 * dasn2[1] + asn2[1] * dnulg[1] ;

      double bpb_et1 = x_et1 * x_et1 - 2. * nyso[0] * x_et1 + npnso2[0] ;
      double bpb_et2 = x_et2 * x_et2 - 2. * nyso[1] * x_et2 + npnso2[1] ;

      double cpc_et1 = 0.5 * dnpnso2[0] + x_et1 * (1. - dnyso[0]) - nyso[0] ;
      double cpc_et2 = 0.5 * dnpnso2[1] + x_et2 * (1. - dnyso[1]) - nyso[1] ;

      om2_et1 = dnulg[0] / (andan_et1 * bpb_et1 + asn2[0] * cpc_et1) ;
      om2_et2 = dnulg[1] / (andan_et2 * bpb_et2 + asn2[1] * cpc_et2) ;

    }
    else {

      om2_et1 = dnulg[0] / x_et1 ;
      om2_et2 = dnulg[1] / x_et2 ;

    }

    double ome_et1 = sqrt( om2_et1 ) ;
    double ome_et2 = sqrt( om2_et2 ) ;

    double diff_om1 = 1. - ome_et1 / ome_et2 ;
    double diff_om2 = 1. - ome_et1 / omega ;

    cout << "Binaire::orbit : omega (direct) [rad/s] : "
	 << ome_et1 * f_unit << endl ;

    cout << "Binaire::orbit : relative difference between " << endl
	 << " omega_star1 and omega_star2 (direct) : " << diff_om1 << endl
	 << " omega (direct) and omega (iretation) : " << diff_om2 << endl ;
	 */

}


void Binaire::orbit_eqmass(double fact_omeg_min, double fact_omeg_max,
			   double mass1, double mass2,
			   double& xgg1, double& xgg2) {

  using namespace Unites ;
    
    //-------------------------------------------------------------
    // Evaluation of various quantities at the center of each star
    //-------------------------------------------------------------

    double dnulg[2], asn2[2], dasn2[2], ny[2], dny[2], npn[2], dnpn[2], xgg[2] ;     
    double nyso[2], dnyso[2], npnso2[2], dnpnso2[2], ori_x[2] ;

    for (int i=0; i<2; i++) {
	
	const Map& mp = et[i]->get_mp() ; 

	const Cmp& logn_auto_regu = et[i]->get_logn_auto_regu()() ; 
	const Cmp& logn_comp = et[i]->get_logn_comp()() ; 
	const Cmp& loggam = et[i]->get_loggam()() ; 
	const Cmp& nnn = et[i]->get_nnn()() ; 
	const Cmp& a_car = et[i]->get_a_car()() ; 
	const Tenseur& shift = et[i]->get_shift() ; 
	const Tenseur& d_logn_auto_div = et[i]->get_d_logn_auto_div() ; 

	Tenseur dln_auto_div = d_logn_auto_div ;

	if ( *(dln_auto_div.get_triad()) != ref_triad ) {

	  // Change the basis from spherical coordinate to Cartesian one
	  dln_auto_div.change_triad( mp.get_bvect_cart() ) ;

	  // Change the basis from mapping coordinate to absolute one
	  dln_auto_div.change_triad( ref_triad ) ;

	}

	//----------------------------------
	// Calcul de d/dX( nu + ln(Gamma) ) au centre de l'etoile ---> dnulg[i]
	//----------------------------------
	
	// Facteur de passage x --> X :
	double factx ;
	if (fabs(mp.get_rot_phi()) < 1.e-14) {
	    factx = 1. ; 
	}
	else {
	    if (fabs(mp.get_rot_phi() - M_PI) < 1.e-14) {
		factx = - 1. ; 
	    }
	    else {
		cout << "Binaire::orbit : unknown value of rot_phi !" << endl ;
		abort() ; 
	    }
	}
	    
	Cmp tmp = logn_auto_regu + logn_comp + loggam ;
	
	// ... gradient suivant X : 		
	dnulg[i] = dln_auto_div(0)(0, 0, 0, 0)
	  + factx * tmp.dsdx()(0, 0, 0, 0) ; 
	
	//----------------------------------
	// Calcul de A^2/N^2 au centre de l'etoile ---> asn2[i]
	//----------------------------------

	double nc = nnn(0, 0, 0, 0) ;
	double a2c = a_car(0, 0, 0, 0) ;
	asn2[i] = a2c / (nc * nc) ;
	 
	if ( et[i]->is_relativistic() ) {

	    //----------------------------------
	    // Calcul de d/dX(A^2/N^2) au centre de l'etoile ---> dasn2[i]
	    //----------------------------------

	    double da2c = factx * a_car.dsdx()(0, 0, 0, 0) ; 
	    double dnc =  factx * nnn.dsdx()(0, 0, 0, 0) ;

	    dasn2[i] = ( da2c - 2 * a2c / nc * dnc ) / (nc*nc) ; 

	    //----------------------------------
	    // Calcul de N^Y au centre de l'etoile ---> ny[i]
	    //----------------------------------

	    ny[i] = shift(1)(0, 0, 0, 0) ; 
	    nyso[i] = ny[i] / omega ;
	    
	    //----------------------------------
	    // Calcul de dN^Y/dX au centre de l'etoile ---> dny[i]
	    //----------------------------------
	    
	    dny[i] = factx * shift(1).dsdx()(0, 0, 0, 0) ; 
	    dnyso[i] = dny[i] / omega ;

	    //----------------------------------
	    // Calcul de (N^X)^2 + (N^Y)^2 + (N^Z)^2 
	    //				     au centre de l'etoile ---> npn[i]
	    //----------------------------------

	    tmp = flat_scalar_prod(shift, shift)() ; 

	    npn[i] = tmp(0, 0, 0, 0) ; 
	    npnso2[i] = npn[i] / omega / omega ;

	    //----------------------------------
	    // Calcul de d/dX( (N^X)^2 + (N^Y)^2 + (N^Z)^2 )
	    //				     au centre de l'etoile ---> dnpn[i]
	    //----------------------------------
	    
	    dnpn[i] = factx * tmp.dsdx()(0, 0, 0, 0) ; 
	    dnpnso2[i] = dnpn[i] / omega / omega ;

	}	    // fin du cas relativiste 
	else {
	    dasn2[i] = 0 ; 
	    ny[i] = 0 ; 
	    nyso[i] = 0 ;
	    dny[i] = 0 ; 
	    dnyso[i] = 0 ;
	    npn[i] = 0 ; 
	    npnso2[i] = 0 ;
	    dnpn[i] = 0 ; 
	    dnpnso2[i] = 0 ;
	}

	cout << "Binaire::orbit: central d(nu+log(Gam))/dX : " 
	     << dnulg[i] << endl ; 
	cout << "Binaire::orbit: central A^2/N^2 : " << asn2[i] << endl ; 
	cout << "Binaire::orbit: central d(A^2/N^2)/dX : " << dasn2[i] << endl ; 
	cout << "Binaire::orbit: central N^Y : " << ny[i] << endl ; 
	cout << "Binaire::orbit: central dN^Y/dX : " << dny[i] << endl ; 
	cout << "Binaire::orbit: central N.N : " << npn[i] << endl ; 
	cout << "Binaire::orbit: central d(N.N)/dX : " << dnpn[i] << endl ; 

	//----------------------
	// Pour information seulement : 1/ calcul des positions des "centres de
	//				    de masse"
	//				2/ calcul de dH/dX en r=0
	//-----------------------

        ori_x[i] = (et[i]->get_mp()).get_ori_x() ;

	xgg[i] = factx * (et[i]->xa_barycenter() - ori_x[i]) ;
		 
    } // fin de la boucle sur les etoiles 

    xgg1 = xgg[0] ;
    xgg2 = xgg[1] ;
    
//---------------------------------
//  Position de l'axe de rotation   
//---------------------------------

    int relat = ( et[0]->is_relativistic() ) ? 1 : 0 ;
    double ori_x1 = ori_x[0] ;
    double ori_x2 = ori_x[1] ;

    if ( et[0]->get_eos() == et[1]->get_eos() && mass1 == mass2 ) {

        x_axe = 0. ;

    }
    else {

	Param paraxe ;
	paraxe.add_int(relat) ;
	paraxe.add_double( ori_x1, 0) ;
	paraxe.add_double( ori_x2, 1) ;
	paraxe.add_double( dnulg[0], 2) ;
	paraxe.add_double( dnulg[1], 3) ;
	paraxe.add_double( asn2[0], 4) ;
	paraxe.add_double( asn2[1], 5) ;
	paraxe.add_double( dasn2[0], 6) ;
	paraxe.add_double( dasn2[1], 7) ;
	paraxe.add_double( nyso[0], 8) ;
	paraxe.add_double( nyso[1], 9) ;
	paraxe.add_double( dnyso[0], 10) ;
	paraxe.add_double( dnyso[1], 11) ;
	paraxe.add_double( npnso2[0], 12) ;
	paraxe.add_double( npnso2[1], 13) ;
	paraxe.add_double( dnpnso2[0], 14) ;
	paraxe.add_double( dnpnso2[1], 15) ;

	int nitmax_axe = 200 ; 
	int nit_axe ; 
	double precis_axe = 1.e-13 ;

	x_axe = zerosec(fonc_binaire_axe, paraxe, 0.9*ori_x1, 0.9*ori_x2,
			precis_axe, nitmax_axe, nit_axe) ;

	cout << "Binaire::orbit : Number of iterations in zerosec for x_axe : "
	     << nit_axe << endl ;
    }

    cout << "Binaire::orbit : x_axe [km] : " << x_axe / km << endl ; 

//-------------------------------------
//  Calcul de la vitesse orbitale    
//-------------------------------------

    Param parf ; 
    parf.add_int(relat) ; 
    parf.add_double( (et[0]->get_mp()).get_ori_x(), 0) ; 
    parf.add_double( dnulg[0], 1) ;  
    parf.add_double( asn2[0], 2) ;    
    parf.add_double( dasn2[0], 3) ;    
    parf.add_double( ny[0], 4) ;    
    parf.add_double( dny[0], 5) ;    
    parf.add_double( npn[0], 6) ;    
    parf.add_double( dnpn[0], 7) ;    
    parf.add_double( x_axe, 8) ; 

    double omega1 = fact_omeg_min * omega  ; 
    double omega2 = fact_omeg_max * omega ; 
    cout << "Binaire::orbit: omega1,  omega2 [rad/s] : " 
	 << omega1 * f_unit << "  " << omega2 * f_unit << endl ; 

	// Search for the various zeros in the interval [omega1,omega2]
	// ------------------------------------------------------------
	int nsub = 50 ;  // total number of subdivisions of the interval
	Tbl* azer = 0x0 ;
	Tbl* bzer = 0x0 ; 
	zero_list(fonc_binaire_orbit, parf, omega1, omega2, nsub,
		  azer, bzer) ; 
	
	// Search for the zero closest to the previous value of omega
	// ----------------------------------------------------------
	double omeg_min, omeg_max ; 
	int nzer = azer->get_taille() ; // number of zeros found by zero_list
	cout << "Binaire:orbit : " << nzer << 
	     " zero(s) found in the interval [omega1,  omega2]." << endl ; 
	cout << "omega, omega1, omega2 : " << omega << "  " << omega1
		<< "  " << omega2 << endl ; 
	cout << "azer : " << *azer << endl ;
	cout << "bzer : " << *bzer << endl ;
	
	if (nzer == 0) {
		cout << 
		"Binaire::orbit: WARNING : no zero detected in the interval"
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
    
	cout << "Binaire:orbit : interval selected for the search of the zero : "
		<< endl << "  [" << omeg_min << ", " << omeg_max << "]  =  [" 
		 << omeg_min * f_unit << ", " << omeg_max * f_unit << "] rad/s " << endl ; 
	
	// Computation of the zero in the selected interval by the secant method
	// ---------------------------------------------------------------------
    int nitermax = 200 ; 
    int niter ; 
    double precis = 1.e-13 ;
    omega = zerosec_b(fonc_binaire_orbit, parf, omeg_min, omeg_max,
		    precis, nitermax, niter) ;
    
    cout << "Binaire::orbit : Number of iterations in zerosec for omega : "
	 << niter << endl ; 
	
    cout << "Binaire::orbit : omega [rad/s] : "
	 << omega * f_unit << endl ; 
          

    /* We do not use the method below.
    // Direct calculation
    //--------------------

    double om2_et1 ;
    double om2_et2 ;

    double x_et1 = ori_x1 - x_axe ;
    double x_et2 = ori_x2 - x_axe ;

    if (relat == 1) {

      double andan_et1 = 0.5 * dasn2[0] + asn2[0] * dnulg[0] ;
      double andan_et2 = 0.5 * dasn2[1] + asn2[1] * dnulg[1] ;

      double bpb_et1 = x_et1 * x_et1 - 2. * nyso[0] * x_et1 + npnso2[0] ;
      double bpb_et2 = x_et2 * x_et2 - 2. * nyso[1] * x_et2 + npnso2[1] ;

      double cpc_et1 = 0.5 * dnpnso2[0] + x_et1 * (1. - dnyso[0]) - nyso[0] ;
      double cpc_et2 = 0.5 * dnpnso2[1] + x_et2 * (1. - dnyso[1]) - nyso[1] ;

      om2_et1 = dnulg[0] / (andan_et1 * bpb_et1 + asn2[0] * cpc_et1) ;
      om2_et2 = dnulg[1] / (andan_et2 * bpb_et2 + asn2[1] * cpc_et2) ;

    }
    else {

      om2_et1 = dnulg[0] / x_et1 ;
      om2_et2 = dnulg[1] / x_et2 ;

    }

    double ome_et1 = sqrt( om2_et1 ) ;
    double ome_et2 = sqrt( om2_et2 ) ;

    double diff_om1 = 1. - ome_et1 / ome_et2 ;
    double diff_om2 = 1. - ome_et1 / omega ;

    cout << "Binaire::orbit : omega (direct) [rad/s] : "
	 << ome_et1 * f_unit << endl ;

    cout << "Binaire::orbit : relative difference between " << endl
	 << " omega_star1 and omega_star2 (direct) : " << diff_om1 << endl
	 << " omega (direct) and omega (iretation) : " << diff_om2 << endl ;
	 */

}


//*************************************************
//  Function used for search of the rotation axis
//*************************************************

double  fonc_binaire_axe(double x_rot, const Param& paraxe) {

    int relat = paraxe.get_int() ;

    double ori_x1 = paraxe.get_double(0) ;
    double ori_x2 = paraxe.get_double(1) ;
    double dnulg_1 = paraxe.get_double(2) ;
    double dnulg_2 = paraxe.get_double(3) ;
    double asn2_1 = paraxe.get_double(4) ;
    double asn2_2 = paraxe.get_double(5) ;
    double dasn2_1 = paraxe.get_double(6) ;
    double dasn2_2 = paraxe.get_double(7) ;
    double nyso_1 = paraxe.get_double(8) ;
    double nyso_2 = paraxe.get_double(9) ;
    double dnyso_1 = paraxe.get_double(10) ;
    double dnyso_2 = paraxe.get_double(11) ;
    double npnso2_1 = paraxe.get_double(12) ;
    double npnso2_2 = paraxe.get_double(13) ;
    double dnpnso2_1 = paraxe.get_double(14) ;
    double dnpnso2_2 = paraxe.get_double(15) ;

    double om2_star1 ;
    double om2_star2 ;

    double x1 = ori_x1 - x_rot ;
    double x2 = ori_x2 - x_rot ;

    if (relat == 1) {

      double andan_1 = 0.5 * dasn2_1 + asn2_1 * dnulg_1 ;
      double andan_2 = 0.5 * dasn2_2 + asn2_2 * dnulg_2 ;

      double bpb_1 = x1 * x1 - 2. * nyso_1 * x1 + npnso2_1 ;
      double bpb_2 = x2 * x2 - 2. * nyso_2 * x2 + npnso2_2 ;

      double cpc_1 = 0.5 * dnpnso2_1 + x1 * (1. - dnyso_1) - nyso_1 ;
      double cpc_2 = 0.5 * dnpnso2_2 + x2 * (1. - dnyso_2) - nyso_2 ;

      om2_star1 = dnulg_1 / (andan_1 * bpb_1 + asn2_1 * cpc_1) ;
      om2_star2 = dnulg_2 / (andan_2 * bpb_2 + asn2_2 * cpc_2) ;

    }
    else {

      om2_star1 = dnulg_1 / x1 ;
      om2_star2 = dnulg_2 / x2 ;

    }

    return om2_star1 - om2_star2 ;

}

//*****************************************************************************
//  Fonction utilisee pour la recherche de omega par la methode de la secante
//*****************************************************************************

double fonc_binaire_orbit(double om, const Param& parf) {

    int relat = parf.get_int() ;

    double xc = parf.get_double(0) ; 
    double dnulg = parf.get_double(1) ;
    double asn2 = parf.get_double(2) ;
    double dasn2 = parf.get_double(3) ;
    double ny = parf.get_double(4) ;
    double dny = parf.get_double(5) ;
    double npn = parf.get_double(6) ;
    double dnpn = parf.get_double(7) ;
    double x_axe = parf.get_double(8) ;

    double xx = xc - x_axe ; 
    double om2 = om*om ; 
    
    double dphi_cent ;
    
    if (relat == 1) {
	double bpb = om2 * xx*xx - 2*om * ny * xx + npn ; 
    
	dphi_cent = ( asn2* ( om* (ny + xx*dny) - om2*xx - 0.5*dnpn )
			 - 0.5*bpb* dasn2 ) 
		       / ( 1 - asn2 * bpb ) ; 
    }
    else {
	dphi_cent =  - om2*xx  ; 
    }
		    
    return dnulg + dphi_cent ; 
       
}


}
