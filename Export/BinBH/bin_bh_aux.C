/*
 * Constructor of class Bin_BH (binary black hole exportation)
 * which depends explicitely on Lorene objects.
 */

/*
 *   Copyright (c) 2001  Eric Gourgoulhon
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

char bin_bh_aux_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/12/14 08:59:18  e_gourgoulhon
 * Exportation of Lorene Bhole_binaire object to a Cartesian grid
 *
 * Revision 1.2  2001/12/11 06:44:41  e_gourgoulhon
 * template files
 *
 *
 *
 * $Header$
 *
 */


#include "bin_bh.h"

// C headers
#include <math.h>

// Lorene headers
#include "tenseur.h"
#include "bhole.h"


		    //----------------------------------------//
		    //	    Constructor from LORENE data      //
		    //----------------------------------------//

Bin_BH::Bin_BH(int nbpoints, const double* xi, const double* yi, 
	       const double* zi, const char* filename) 
	       : np(nbpoints) {

    // Reading of data
    // ---------------
    FILE* fich = fopen(filename, "r") ;
    Mg3d grille (fich) ;
    Map_af map_un (grille, fich) ;
    Map_af map_deux (grille, fich) ;
    Bhole hole_un (map_un, fich) ;
    Bhole hole_deux (map_deux, fich) ;
    fclose(fich) ;
    
    assert (hole_un.get_omega() == hole_deux.get_omega()) ;
    
    // Construction of the binary system
    // ---------------------------------
    Bhole_binaire systeme (map_un, map_deux) ;
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
    
    omega = systeme.get_omega() * aa ;
    dist = ( map_un.get_ori_x() - map_deux.get_ori_x() ) / aa ;
    radius2 = systeme(2).get_rayon() / aa ; 
   
    cout << endl << "Binary system read in file : " << endl ; 
    cout <<	    "---------------------------- " << endl ;
    cout << "  Separation d/a :       " << dist << endl ;
    cout << "  Omega :                " << omega << " / a" << endl ; 
    cout << "  Size of black hole 2 : " << radius2 << " a" << endl ; 
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
    
    const Map_af& mp1 = systeme(1).get_mp() ; 
    const Map_af& mp2 = systeme(2).get_mp() ; 
    
    const Valeur& vnn1 = (systeme(1).get_n_auto()()).va ;
    const Valeur& vnn2 = (systeme(2).get_n_auto()()).va ;
    vnn1.coef() ;		// The sprectral coefficients are required
    vnn2.coef() ;		
    
    const Valeur& vbetax1 = (systeme(1).get_shift_auto()(0)).va ;
    const Valeur& vbetax2 = (systeme(2).get_shift_auto()(0)).va ;
    const Valeur& vbetay1 = (systeme(1).get_shift_auto()(1)).va ;
    const Valeur& vbetay2 = (systeme(2).get_shift_auto()(1)).va ;
    const Valeur& vbetaz1 = (systeme(1).get_shift_auto()(2)).va ;
    const Valeur& vbetaz2 = (systeme(2).get_shift_auto()(2)).va ;
    vbetax1.coef() ; 
    vbetax2.coef() ; 
    vbetay1.coef() ; 
    vbetay2.coef() ; 
    vbetaz1.coef() ; 
    vbetaz2.coef() ; 

    const Valeur& vpsi1 = (systeme(1).get_psi_auto()()).va ;
    const Valeur& vpsi2 = (systeme(2).get_psi_auto()()).va ;
    vpsi1.coef() ;		
    vpsi2.coef() ;		
    
    Tenseur_sym k_un (systeme(1).get_tkij_auto()) ;
    k_un.set_std_base() ;
    k_un.dec2_dzpuis() ;
    Tenseur_sym k_deux (systeme(2).get_tkij_auto()) ;
    k_deux.set_std_base() ;
    k_deux.dec2_dzpuis() ;
    
    Valeur vkxx1 = (k_un(0, 0)).va ;
    Valeur vkxx2 = (k_deux(0, 0)).va ;
    Valeur vkxy1 = (k_un(0, 1)).va ;
    Valeur vkxy2 = (k_deux(0, 1)).va ;
    Valeur vkxz1 = (k_un(0, 2)).va ;
    Valeur vkxz2 = (k_deux(0, 2)).va ;
    Valeur vkyy1 = (k_un(1, 1)).va ;
    Valeur vkyy2 = (k_deux(1, 1)).va ;
    Valeur vkyz1 = (k_un(1, 2)).va ;
    Valeur vkyz2 = (k_deux(1, 2)).va ;
    Valeur vkzz1 = (k_un(2, 2)).va ;
    Valeur vkzz2 = (k_deux(2, 2)).va ;
    
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
    
    
    for (int i=0; i<np; i++) {
    
	double x0 = xx[i] * aa ;    // x in Lorene's unit
	double y0 = yy[i] * aa ;
	double z0 = zz[i] * aa ;
    
	// Values of (l1, xi1, theta1, phi1) (grid 1) 
	// corresponding to (x,y,z):
	// ------------------------------------------
	double r1, theta1, phi1 ;   // polar coordinates centered on b.h. 1
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
	
	if ( (r1 < aa) || (r2 < systeme(2).get_rayon()) ) {
	    // We are "inside" one of the throats:
	    nnn[i] = 0 ; 
	    beta_x[i] = 0 ; 
	    beta_y[i] = 0 ; 
	    beta_z[i] = 0 ; 
	    g_xx[i] = 0 ; 
	    g_xy[i] = 0 ; 
	    g_xz[i] = 0 ; 
	    g_yy[i] = 0 ; 
	    g_yz[i] = 0 ; 
	    g_zz[i] = 0 ; 
	    k_xx[i] = 0 ; 
	    k_xy[i] = 0 ; 
	    k_xz[i] = 0 ; 
	    k_yy[i] = 0 ; 
	    k_yz[i] = 0 ; 
	    k_zz[i] = 0 ; 
	}
	else {
	
	// Lapse function
	// --------------
	
	nnn[i] =    vnn1.c_cf->val_point_symy(l1, xi1, theta1, phi1) 
		 +  vnn2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ;
	
	// Shift vector
	// ------------
	
	beta_x[i] = vbetax1.c_cf->val_point_asymy(l1, xi1, theta1, phi1) 
		 -  vbetax2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) 
		 - omega * yy[i] ;

	beta_y[i] = vbetay1.c_cf->val_point_symy(l1, xi1, theta1, phi1) 
		 -  vbetay2.c_cf->val_point_symy(l2, xi2, theta2, phi2) 
		 + omega * xx[i] ;

	beta_z[i] = vbetaz1.c_cf->val_point_asymy(l1, xi1, theta1, phi1) 
		 +  vbetaz2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ;

	
	// 3-metric
	// --------
	
	double psi4 = pow( vpsi1.c_cf->val_point_symy(l1, xi1, theta1, phi1) 
		   +  vpsi2.c_cf->val_point_symy(l2, xi2, theta2, phi2), 4) ; 
	
	g_xx[i] = psi4 ; 

	g_yy[i] = psi4  ; 	
	g_zz[i] = psi4  ; 	

	g_xy[i] = 0 ; 
	g_xz[i] = 0 ; 
	g_yz[i] = 0 ; 
			
	// Extrinsic curvature
	// -------------------
	
	double pre = aa * psi4 ; 
	
	k_xx[i] = pre * ( vkxx1.c_cf->val_point_asymy(l1, xi1, theta1, phi1) 
		      +	  vkxx2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ) ;
		
	k_xy[i] = pre * ( vkxy1.c_cf->val_point_symy(l1, xi1, theta1, phi1) 
		      +   vkxy2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ) ;
		
	k_xz[i] = pre * ( vkxz1.c_cf->val_point_asymy(l1, xi1, theta1, phi1) 
		      -   vkxz2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ) ;
		
	k_yy[i] = pre * ( vkyy1.c_cf->val_point_asymy(l1, xi1, theta1, phi1) 
		      +   vkyy2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ) ;
		
	k_yz[i] = pre * ( vkyz1.c_cf->val_point_symy(l1, xi1, theta1, phi1) 
		      -   vkyz2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ) ;
		
	k_zz[i] = pre * ( vkzz1.c_cf->val_point_asymy(l1, xi1, theta1, phi1) 
		      +   vkzz2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ) ;
		      	
	}

    }	// End of loop on the points
    


}
		    


