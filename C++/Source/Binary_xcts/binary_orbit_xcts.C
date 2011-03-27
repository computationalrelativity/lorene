/*
 * Method of class Binary_xcts to compute the orbital angular velocity 
 * {\tt omega} and the position of the rotation axis {\tt x_axe}.
 * (See file binary_xcts.h for documentation)
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

char binary_orbit_xcts_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.9  2011/03/27 16:41:19  e_gourgoulhon
 * -- Corrected sign of ny and dny for star no. 2
 * -- Added output via new function save_profile for graphics
 *
 * Revision 1.8  2011/03/27 14:58:48  m_bejger
 * dnulg by means of dsdx(); rearrangements to use primary variables
 *
 * Revision 1.7  2011/03/25 16:28:36  e_gourgoulhon
 * Still in progress
 *
 * Revision 1.6  2010/12/09 10:41:20  m_bejger
 * For testing; not sure if working properly
 *
 * Revision 1.5  2010/10/26 19:45:45  m_bejger
 * Cleanup
 *
 * Revision 1.4  2010/07/16 16:27:19  m_bejger
 * This version is basically a copy of the one used by Binaire (binaire_orbite.C)
 *
 * Revision 1.3  2010/06/17 14:15:41  m_bejger
 * Using method get_Psi()
 *
 * Revision 1.2  2010/06/15 07:57:30  m_bejger
 * Minor corrections
 *
 * Revision 1.1  2010/05/04 07:35:54  m_bejger
 * Initial version
 *
 * $Header$
 *
 */

// Headers C
#include <math.h>

// Headers Lorene 
#include "binary_xcts.h"
#include "eos.h"
#include "param.h"
#include "utilitaires.h"
#include "unites.h"

#include "graphique.h"

double  fonc_binary_xcts_axe(double , const Param& ) ;
double  fonc_binary_xcts_orbit(double , const Param& ) ;

//******************************************************************************

void Binary_xcts::orbit(double fact_omeg_min, 
						double fact_omeg_max, 
						double& xgg1, 
						double& xgg2) {

  using namespace Unites ;
    
    //-------------------------------------------------------------
    // Evaluation of various quantities at the center of each star
    //-------------------------------------------------------------

    double dnulg[2], asn2[2], dasn2[2], ny[2], dny[2], npn[2], dnpn[2], xgg[2] ;     
    double nyso[2], dnyso[2], npnso2[2], dnpnso2[2], ori_x[2] ;

    for (int i=0; i<2; i++) {
	
	const Map& mp = et[i]->get_mp() ; 
	const Metric& flat = et[i]->get_flat() ;
		
	const Scalar& chi = et[i]->get_chi() ;
	const Scalar& Psi = et[i]->get_Psi() ;
	const Scalar& psi4 = et[i]->get_psi4() ;

	const Vector& dchi = et[i]->get_dcov_chi() ;
	const Vector& dPsi = et[i]->get_dcov_Psi() ;
	
	const Scalar& nn = et[i]->get_nn() ; 			
	const Scalar& loggam = et[i]->get_loggam() ; 
		
	Scalar dlnnlng = chi.dsdx()/chi - Psi.dsdx()/Psi + loggam.dsdx() ;  
	
	//Vector dlnnlng = dchi / Psi - chi/(Psi*Psi)*dPsi + loggam.derive_cov(flat); 
	// Change the basis from spherical coordinate to Cartesian one
	//dlnnlng.change_triad( mp.get_bvect_cart() ) ;

	// Sign convention for shift (beta^i = - N^i)
	Vector shift =  - ( et[i]->get_beta() ) ;
	shift.change_triad(et[i]->mp.get_bvect_cart()) ;

	//------------------------------------------------------------------
	// d/dX( log(N) + log(Gamma) ) 
	// in the center of the star ---> dnulg[i]
	//------------------------------------------------------------------
	
	// Facteur de passage x --> X :
	double factx ;
	if (fabs(mp.get_rot_phi()) < 1.e-14) {
	    factx = 1. ; 

	} else {
		
	  	  if (fabs(mp.get_rot_phi() - M_PI) < 1.e-14) {
			factx = - 1. ; 
			
	      } else {
			  
			cout << "Binary_xcts::orbit : unknown value of rot_phi !" << endl ;
			abort() ; 
	      }
	}
	    
	//dnulg[i] = factx*dlnnlng(1).val_grid_point(0, 0, 0, 0) ;
	dnulg[i] = factx*dlnnlng.val_grid_point(0, 0, 0, 0) ;

	//## Alternative : 
	// Scalar tmp = et[i]->get_logn() + et[i]->get_loggam() ; 
	// dnulg[i] = factx*tmp.dsdx().val_grid_point(0, 0, 0, 0) ; 
		
	// For graphical outputs: 
	Scalar tgraph = et[i]->get_logn() - log( (1. + et[i]->get_chi_auto()) / (1. + et[i]->get_Psi_auto()) ) ; 
	// tmp = log( (1. + et[i]->get_chi_comp()) / (1. + et[i]->get_Psi_comp()) ) ; 
	tgraph.std_spectral_base() ; 
        save_profile(tgraph, 0., 10., 0.5*M_PI, 0., "prof_logn.d") ; 
        save_profile(et[i]->get_loggam(), 0., 1.8, 0.5*M_PI, 0., "prof_loggam.d") ;
 
	//------------------------------------------------------------------
	// Psi^4/N^2 = in the center of the star ---> asn2[i]
	//------------------------------------------------------------------

	Scalar Psi6schi2 = pow(Psi, 6)/(chi*chi) ; 
	Psi6schi2.std_spectral_base() ; 
	asn2[i] = Psi6schi2.val_grid_point(0, 0, 0, 0) ;  
	 
	//------------------------------------------------------------------
	// d/dX(A^2/N^2) in the center of the star ---> dasn2[i]
	//------------------------------------------------------------------

	//Scalar dPsi6schi2 = Psi6schi2.dsdx() ; 
	dasn2[i] = Psi6schi2.dsdx().val_grid_point(0, 0, 0, 0) ; 
		
	//------------------------------------------------------------------
	// N^Y in the center of the star ---> ny[i]
	//------------------------------------------------------------------

	//## ny[i] = shift(2).val_grid_point(0, 0, 0, 0) ; 
	ny[i] = factx*shift(2).val_grid_point(0, 0, 0, 0) ; 

	nyso[i] = ny[i] / omega ;
	    
	//------------------------------------------------------------------
	// dN^Y/dX in the center of the star ---> dny[i]
	//------------------------------------------------------------------
	    
	//## dny[i] = factx * shift(2).dsdx().val_grid_point(0, 0, 0, 0) ; 
	dny[i] = shift(2).dsdx().val_grid_point(0, 0, 0, 0) ; 

	dnyso[i] = dny[i] / omega ;

	    //--------------------------------------------------------------
	    // (N^X)^2 + (N^Y)^2 + (N^Z)^2 
	    // in the center of the star ---> npn[i]
	    //--------------------------------------------------------------
 
        Scalar tmp = contract(shift, 0, shift.up_down(flat), 0) ; 

	    npn[i] = tmp.val_grid_point(0, 0, 0, 0) ; 
	    npnso2[i] = npn[i] / omega / omega ;

	    //--------------------------------------------------------------
	    // d/dX( (N^X)^2 + (N^Y)^2 + (N^Z)^2 )
	    // in the center of the star ---> dnpn[i]
	    //--------------------------------------------------------------
	    
	    dnpn[i] = factx * tmp.dsdx().val_grid_point(0, 0, 0, 0) ; 
	    dnpnso2[i] = dnpn[i] / omega / omega ;

	cout << "Binary_xcts::orbit: central d(nu+log(Gam))/dX : " 
	     << dnulg[i] << endl ; 
	cout << "Binary_xcts::orbit: central A^2/N^2 : " 
	     << asn2[i] << endl ; 
	cout << "Binary_xcts::orbit: central d(A^2/N^2)/dX : " 
	     << dasn2[i] << endl ; 
	cout << "Binary_xcts::orbit: central N^Y : " 
	     << ny[i] << endl ; 
	cout << "Binary_xcts::orbit: central dN^Y/dX : " 
	     << dny[i] << endl ; 
	cout << "Binary_xcts::orbit: central N.N : " 
	     << npn[i] << endl ; 
	cout << "Binary_xcts::orbit: central d(N.N)/dX : " 
	     << dnpn[i] << endl ; 


    ori_x[i] = (et[i]->get_mp()).get_ori_x() ;
	xgg[i] = factx * (et[i]->xa_barycenter() - ori_x[i]) ;
		 
    } 

    xgg1 = xgg[0] ;
    xgg2 = xgg[1] ;
    
//---------------------------------
//  axis of rotation    
//---------------------------------

    int relat = 1 ; 
     
    double ori_x1 = ori_x[0] ;
    double ori_x2 = ori_x[1] ;

    if ( et[0]->get_eos() == et[1]->get_eos() &&
	 fabs( et[0]->get_ent().val_grid_point(0,0,0,0) - 
		   et[1]->get_ent().val_grid_point(0,0,0,0) ) < 1.e-14 ) {

        x_axe = 0. ;

    } else {

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

	x_axe = zerosec(fonc_binary_xcts_axe, paraxe, 0.9*ori_x1, 0.9*ori_x2,
			precis_axe, nitmax_axe, nit_axe) ;

	cout << "Binary_xcts::orbit : Number of iterations in zerosec for x_axe : "
	     << nit_axe << endl ;
    }

    cout << "Binary_xcts::orbit: x_axe [km] : " << x_axe / km << endl ; 

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
    cout << "Binary_xcts::orbit: omega1,  omega2 [rad/s] : " 
	 << omega1 * f_unit << "  " << omega2 * f_unit << endl ; 

	// Search for the various zeros in the interval [omega1, omega2]
	// ------------------------------------------------------------
	int nsub = 50 ;  // total number of subdivisions of the interval
	Tbl* azer = 0x0 ;
	Tbl* bzer = 0x0 ; 
	zero_list(fonc_binary_xcts_orbit, parf, omega1, omega2, nsub,
		  azer, bzer) ; 
	
	// Search for the zero closest to the previous value of omega
	// ----------------------------------------------------------
	double omeg_min, omeg_max ; 
	int nzer = azer->get_taille() ; // number of zeros found by zero_list
	cout << "Binary_xcts:orbit : " << nzer << 
	     " zero(s) found in the interval [omega1,  omega2]." << endl ; 
	cout << "omega, omega1, omega2 : " << omega << "  " << omega1
		<< "  " << omega2 << endl ; 
	cout << "azer : " << *azer << endl ;
	cout << "bzer : " << *bzer << endl ;
	
	if (nzer == 0) {
		cout << 
		"Binary_xcts::orbit: WARNING : no zero detected in the interval"
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
    
	cout << "Binary_xcts::orbit : interval selected for the search of the zero : "
		<< endl << "  [" << omeg_min << ", " << omeg_max << "]  =  [" 
		 << omeg_min * f_unit << ", " << omeg_max * f_unit << "] rad/s " << endl ; 
	
	// Computation of the zero in the selected interval by the secant method
	// ---------------------------------------------------------------------
    int nitermax = 200 ; 
    int niter ; 
    double precis = 1.e-13 ;
    omega = zerosec_b(fonc_binary_xcts_orbit, parf, omeg_min, omeg_max,
		    precis, nitermax, niter) ;
    
    cout << "Binary_xcts::orbit : Number of iterations in zerosec for omega : "
	 << niter << endl ; 
	
    cout << "Binary_xcts::orbit : omega [rad/s] : "
	 << omega * f_unit << endl ; 
          
}

//*************************************************
//  Function used for search of the rotation axis
//*************************************************

double  fonc_binary_xcts_axe(double x_rot, const Param& paraxe) {

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

double fonc_binary_xcts_orbit(double om, const Param& parf) {

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


