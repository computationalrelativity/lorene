/*
 * Constructor of class Bin_BH (binary black hole exportation)
 * which depends explicitely on Lorene objects.
 *
 * (see file bin_bh.h for documentation).
 */

/*
 *   Copyright (c) 2001-2002  Eric Gourgoulhon
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
 * Revision 1.3  2002/03/20 08:24:56  e_gourgoulhon
 * Added the derivatives of Psi.
 *
 * Revision 1.2  2001/12/19 11:20:56  e_gourgoulhon
 * Initialisation of radius2 was missing !
 *
 * Revision 1.1  2001/12/18 22:27:04  e_gourgoulhon
 * Exportation of Lorene structures
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

// Local prototype:
double lagrange_parabol(double x, const double* xp, const double* yp) ;

		    //----------------------------------------//
		    //	    Constructor from LORENE data      //
		    //----------------------------------------//

Bin_BH::Bin_BH(int nbpoints, const double* xi, const double* yi,
	       const double* zi, int fill, const char* filename)
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
    double aasq = aa * aa ;
    double aa2 = systeme(2).get_rayon() ;
    radius2 = aa2 ;

    omega = systeme.get_omega() * aa ;
    dist = ( map_un.get_ori_x() - map_deux.get_ori_x() ) / aa ;

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

    const Cmp& cnn1 = systeme(1).get_n_auto()() ;
    const Cmp& cnn2 = systeme(2).get_n_auto()() ;
    const Valeur& vnn1 = cnn1.va ;
    const Valeur& vnn2 = cnn2.va ;
    vnn1.coef() ;		// The sprectral coefficients are required
    vnn2.coef() ;		

    const Cmp& cbetax1 = systeme(1).get_shift_auto()(0) ;
    const Cmp& cbetax2 = systeme(2).get_shift_auto()(0) ;
    const Cmp& cbetay1 = systeme(1).get_shift_auto()(1) ;
    const Cmp& cbetay2 = systeme(2).get_shift_auto()(1) ;
    const Cmp& cbetaz1 = systeme(1).get_shift_auto()(2) ;
    const Cmp& cbetaz2 = systeme(2).get_shift_auto()(2) ;
    const Valeur& vbetax1 = cbetax1.va ;
    const Valeur& vbetax2 = cbetax2.va ;
    const Valeur& vbetay1 = cbetay1.va ;
    const Valeur& vbetay2 = cbetay2.va ;
    const Valeur& vbetaz1 = cbetaz1.va ;
    const Valeur& vbetaz2 = cbetaz2.va ;
    vbetax1.coef() ;
    vbetax2.coef() ;
    vbetay1.coef() ; 
    vbetay2.coef() ; 
    vbetaz1.coef() ; 
    vbetaz2.coef() ; 

    const Cmp& cpsi1 = systeme(1).get_psi_auto()() ;
    const Cmp& cpsi2 = systeme(2).get_psi_auto()() ;
    const Valeur& vpsi1 = cpsi1.va ;
    const Valeur& vpsi2 = cpsi2.va ;
    vpsi1.coef() ;		
    vpsi2.coef() ;
    
    Tenseur_sym k_un (systeme(1).get_tkij_auto()) ;
    k_un.set_std_base() ;
    k_un.dec2_dzpuis() ;
    Tenseur_sym k_deux (systeme(2).get_tkij_auto()) ;
    k_deux.set_std_base() ;
    k_deux.dec2_dzpuis() ;

    const Cmp& ckxx1 = k_un(0, 0) ;
    const Cmp& ckxy1 = k_un(0, 1) ;
    const Cmp& ckxz1 = k_un(0, 2) ;
    const Cmp& ckyy1 = k_un(1, 1) ;
    const Cmp& ckyz1 = k_un(1, 2) ;
    const Cmp& ckzz1 = k_un(2, 2) ;
    const Cmp& ckxx2 = k_deux(0, 0) ;
    const Cmp& ckxy2 = k_deux(0, 1) ;
    const Cmp& ckxz2 = k_deux(0, 2) ;
    const Cmp& ckyy2 = k_deux(1, 1) ;
    const Cmp& ckyz2 = k_deux(1, 2) ;
    const Cmp& ckzz2 = k_deux(2, 2) ;

    const Valeur& vkxx1 = ckxx1.va ;
    const Valeur& vkxy1 = ckxy1.va ;
    const Valeur& vkxz1 = ckxz1.va ;
    const Valeur& vkyy1 = ckyy1.va ;
    const Valeur& vkyz1 = ckyz1.va ;
    const Valeur& vkzz1 = ckzz1.va ;
    const Valeur& vkxx2 = ckxx2.va ;
    const Valeur& vkxy2 = ckxy2.va ;
    const Valeur& vkxz2 = ckxz2.va ;
    const Valeur& vkyy2 = ckyy2.va ;
    const Valeur& vkyz2 = ckyz2.va ;
    const Valeur& vkzz2 = ckzz2.va ;

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

    // First derivatives of psi
    //-------------------------

    Tenseur dpsi1 = (systeme(1).get_psi_auto()).gradient() ;
    dpsi1.dec2_dzpuis() ;

    Tenseur dpsi2 = (systeme(2).get_psi_auto()).gradient() ;
    dpsi2.dec2_dzpuis() ;

    const Cmp& cdpsix1 = dpsi1(0) ;
    const Cmp& cdpsiy1 = dpsi1(1) ;
    const Cmp& cdpsiz1 = dpsi1(2) ;
    const Cmp& cdpsix2 = dpsi2(0) ;
    const Cmp& cdpsiy2 = dpsi2(1) ;
    const Cmp& cdpsiz2 = dpsi2(2) ;

    const Valeur& vdpsix1 = cdpsix1.va ;
    const Valeur& vdpsiy1 = cdpsiy1.va ;
    const Valeur& vdpsiz1 = cdpsiz1.va ;
    const Valeur& vdpsix2 = cdpsix2.va ;
    const Valeur& vdpsiy2 = cdpsiy2.va ;
    const Valeur& vdpsiz2 = cdpsiz2.va ;

    vdpsix1.coef() ;
    vdpsiy1.coef() ;
    vdpsiz1.coef() ;
    vdpsix2.coef() ;
    vdpsiy2.coef() ;
    vdpsiz2.coef() ;

    // Second derivatives of psi
    //---------------------------

    Tenseur_sym d2psi1( dpsi1.gradient() ) ;
    d2psi1.dec2_dzpuis() ;

    Tenseur_sym d2psi2( dpsi2.gradient() ) ;
    d2psi2.dec2_dzpuis() ;

    const Cmp& cd2psixx1 = d2psi1(0, 0) ;
    const Cmp& cd2psixy1 = d2psi1(0, 1) ;
    const Cmp& cd2psixz1 = d2psi1(0, 2) ;
    const Cmp& cd2psiyy1 = d2psi1(1, 1) ;
    const Cmp& cd2psiyz1 = d2psi1(1, 2) ;
    const Cmp& cd2psizz1 = d2psi1(2, 2) ;
    const Cmp& cd2psixx2 = d2psi2(0, 0) ;
    const Cmp& cd2psixy2 = d2psi2(0, 1) ;
    const Cmp& cd2psixz2 = d2psi2(0, 2) ;
    const Cmp& cd2psiyy2 = d2psi2(1, 1) ;
    const Cmp& cd2psiyz2 = d2psi2(1, 2) ;
    const Cmp& cd2psizz2 = d2psi2(2, 2) ;

    const Valeur& vd2psixx1 = cd2psixx1.va ;
    const Valeur& vd2psixy1 = cd2psixy1.va ;
    const Valeur& vd2psixz1 = cd2psixz1.va ;
    const Valeur& vd2psiyy1 = cd2psiyy1.va ;
    const Valeur& vd2psiyz1 = cd2psiyz1.va ;
    const Valeur& vd2psizz1 = cd2psizz1.va ;
    const Valeur& vd2psixx2 = cd2psixx2.va ;
    const Valeur& vd2psixy2 = cd2psixy2.va ;
    const Valeur& vd2psixz2 = cd2psixz2.va ;
    const Valeur& vd2psiyy2 = cd2psiyy2.va ;
    const Valeur& vd2psiyz2 = cd2psiyz2.va ;
    const Valeur& vd2psizz2 = cd2psizz2.va ;

    vd2psixx1.coef() ;
    vd2psixy1.coef() ;
    vd2psixz1.coef() ;
    vd2psiyy1.coef() ;
    vd2psiyz1.coef() ;
    vd2psizz1.coef() ;
    vd2psixx2.coef() ;
    vd2psixy2.coef() ;
    vd2psixz2.coef() ;
    vd2psiyy2.coef() ;
    vd2psiyz2.coef() ;
    vd2psizz2.coef() ;

    // Arrays describing the 3 points used for the parabolic extrapolation
    //  "inside" the throats when fill = 1
	double r1p[3], t1p[3], p1p[3] ;
	double r2p[3], t2p[3], p2p[3] ;
	double yp[3] ;

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
	
	//------------------------------------------------------------
	// 			"Inside" hole 1
	//------------------------------------------------------------	
		
	if (r1 < aa) {  	

		switch (fill) {
			case 0 : {
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
                                dpsi_x[i] = 0 ;
                                dpsi_y[i] = 0 ;
                                dpsi_z[i] = 0 ;
                                d2psi_xx[i] = 0 ;
                                d2psi_xy[i] = 0 ;
                                d2psi_xz[i] = 0 ;
                                d2psi_yy[i] = 0 ;
                                d2psi_yz[i] = 0 ;
                                d2psi_zz[i] = 0 ;
	    			break ;
	    		}
	    		
	    		case 1 : {
	    		
		r1p[0] = aa ;    	// 3 points outside the throat
		r1p[1] = 1.1 * aa ;     // for the parabolic extrapolation
		r1p[2] = 1.2 * aa ;     //
		
		double rot_phi = mp1.get_rot_phi() ;
		double orix = mp1.get_ori_x() ;
		double oriy = mp1.get_ori_y() ;
		double oriz = mp1.get_ori_z() ;
		
		for (int j=0; j<3; j++) {
			double phi = phi1 + rot_phi ;
			double xap = r1p[j] * sin(theta1) * cos(phi) + orix ;
			double yap = r1p[j] * sin(theta1) * sin(phi) + oriy ;
			double zap = r1p[j] * cos(theta1) + oriz ;
			mp2.convert_absolute(xap, yap, zap,
						r2p[j], t2p[j], p2p[j]) ;
		}

		// Lapse function
		for (int j=0; j<3; j++) {
			yp[j] =   cnn1.val_point(r1p[j], theta1, phi1)
				+ cnn2.val_point(r2p[j], t2p[j], p2p[j])  ;
                }			
                nnn[i] = lagrange_parabol(r1, r1p, yp) ;
                	    		
		// Shift vector
		for (int j=0; j<3; j++) {
			yp[j] =   cbetax1.val_point(r1p[j], theta1, phi1)
				- cbetax2.val_point(r2p[j], t2p[j], p2p[j]) ;
		}
  		beta_x[i] = lagrange_parabol(r1, r1p, yp) - omega * yy[i] ;

		for (int j=0; j<3; j++) {
			yp[j] =   cbetay1.val_point(r1p[j], theta1, phi1)
				- cbetay2.val_point(r2p[j], t2p[j], p2p[j]) ;
		}
  		beta_y[i] = lagrange_parabol(r1, r1p, yp) + omega * xx[i] ;

		for (int j=0; j<3; j++) {
			yp[j] =   cbetaz1.val_point(r1p[j], theta1, phi1)
				+ cbetaz2.val_point(r2p[j], t2p[j], p2p[j]) ;
		}
  		beta_z[i] = lagrange_parabol(r1, r1p, yp) ;

		// 3-metric
		for (int j=0; j<3; j++) {
			yp[j] =   cpsi1.val_point(r1p[j], theta1, phi1)
				+ cpsi2.val_point(r2p[j], t2p[j], p2p[j])  ;
                }			
                double psi4 = pow( lagrange_parabol(r1, r1p, yp), 4) ;
		g_xx[i] = psi4 ;
		g_yy[i] = psi4 ; 	
		g_zz[i] = psi4 ; 	
		g_xy[i] = 0 ;
		g_xz[i] = 0 ;
		g_yz[i] = 0 ;

 		// Extrinsic curvature
		double pre = aa * psi4 ;
		for (int j=0; j<3; j++) {
			yp[j] =   ckxx1.val_point(r1p[j], theta1, phi1)
				+ ckxx2.val_point(r2p[j], t2p[j], p2p[j])  ;
                }			
                k_xx[i] = pre * lagrange_parabol(r1, r1p, yp) ;

		for (int j=0; j<3; j++) {
			yp[j] =   ckxy1.val_point(r1p[j], theta1, phi1)
				+ ckxy2.val_point(r2p[j], t2p[j], p2p[j])  ;
                }
                k_xy[i] = pre * lagrange_parabol(r1, r1p, yp) ;

		for (int j=0; j<3; j++) {
			yp[j] =   ckxz1.val_point(r1p[j], theta1, phi1)
				- ckxz2.val_point(r2p[j], t2p[j], p2p[j])  ;
                }			
                k_xz[i] = pre * lagrange_parabol(r1, r1p, yp) ;

		for (int j=0; j<3; j++) {
			yp[j] =   ckyy1.val_point(r1p[j], theta1, phi1)
				+ ckyy2.val_point(r2p[j], t2p[j], p2p[j])  ;
                }			
                k_yy[i] = pre * lagrange_parabol(r1, r1p, yp) ;

		for (int j=0; j<3; j++) {
			yp[j] =   ckyz1.val_point(r1p[j], theta1, phi1)
				- ckyz2.val_point(r2p[j], t2p[j], p2p[j])  ;
                }			
                k_yz[i] = pre * lagrange_parabol(r1, r1p, yp) ;

		for (int j=0; j<3; j++) {
			yp[j] =   ckzz1.val_point(r1p[j], theta1, phi1)
				+ ckzz2.val_point(r2p[j], t2p[j], p2p[j])  ;
                }			
                k_zz[i] = pre * lagrange_parabol(r1, r1p, yp) ;

                // First derivatives of Psi
		for (int j=0; j<3; j++) {
			yp[j] = aa * ( cdpsix1.val_point(r1p[j], theta1, phi1)
				- cdpsix2.val_point(r2p[j], t2p[j], p2p[j]) ) ;
		}
  		dpsi_x[i] = lagrange_parabol(r1, r1p, yp)  ;

                for (int j=0; j<3; j++) {
			yp[j] = aa * ( cdpsiy1.val_point(r1p[j], theta1, phi1)
				- cdpsiy2.val_point(r2p[j], t2p[j], p2p[j]) ) ;
		}
  		dpsi_y[i] = lagrange_parabol(r1, r1p, yp)  ;

		for (int j=0; j<3; j++) {
			yp[j] = aa * ( cdpsiz1.val_point(r1p[j], theta1, phi1)
				+ cdpsiz2.val_point(r2p[j], t2p[j], p2p[j]) ) ;
		}
  		dpsi_z[i] = lagrange_parabol(r1, r1p, yp)  ;

               // Second derivatives of Psi
		for (int j=0; j<3; j++) {
			yp[j] = aasq * ( cd2psixx1.val_point(r1p[j], theta1, phi1)
				+ cd2psixx2.val_point(r2p[j], t2p[j], p2p[j]) ) ;
		}
  		d2psi_xx[i] = lagrange_parabol(r1, r1p, yp)  ;

		for (int j=0; j<3; j++) {
			yp[j] = aasq * ( cd2psixy1.val_point(r1p[j], theta1, phi1)
				+ cd2psixy2.val_point(r2p[j], t2p[j], p2p[j]) ) ;
		}
  		d2psi_xy[i] = lagrange_parabol(r1, r1p, yp)  ;

		for (int j=0; j<3; j++) {
			yp[j] = aasq * ( cd2psixz1.val_point(r1p[j], theta1, phi1)
				- cd2psixz2.val_point(r2p[j], t2p[j], p2p[j]) ) ;
		}
  		d2psi_xz[i] = lagrange_parabol(r1, r1p, yp)  ;

		for (int j=0; j<3; j++) {
			yp[j] = aasq * ( cd2psiyy1.val_point(r1p[j], theta1, phi1)
				+ cd2psiyy2.val_point(r2p[j], t2p[j], p2p[j]) ) ;
		}
  		d2psi_yy[i] = lagrange_parabol(r1, r1p, yp)  ;

		for (int j=0; j<3; j++) {
			yp[j] = aasq * ( cd2psiyz1.val_point(r1p[j], theta1, phi1)
				- cd2psiyz2.val_point(r2p[j], t2p[j], p2p[j]) ) ;
		}
  		d2psi_yz[i] = lagrange_parabol(r1, r1p, yp)  ;

		for (int j=0; j<3; j++) {
			yp[j] = aasq * ( cd2psizz1.val_point(r1p[j], theta1, phi1)
				+ cd2psizz2.val_point(r2p[j], t2p[j], p2p[j]) ) ;
		}
  		d2psi_zz[i] = lagrange_parabol(r1, r1p, yp)  ;


		break ;
	    		}	// end of case fill = 1
	    		
	    		default : {
	    		
	  cout << "Bin_BH::Bin_BH : the case fill = " << fill
	  	<< " is not known !" << endl ;
	  abort() ;
	  break ;
	    		}
		 }
	} // End of case r1 < a
	

	//------------------------------------------------------------
	// 			"Inside" hole 2
	//------------------------------------------------------------	
		
	if (r2 < aa2) {  	

		switch (fill) {
			case 0 : {
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
                                dpsi_x[i] = 0 ;
                                dpsi_y[i] = 0 ;
                                dpsi_z[i] = 0 ;
                                d2psi_xx[i] = 0 ;
                                d2psi_xy[i] = 0 ;
                                d2psi_xz[i] = 0 ;
                                d2psi_yy[i] = 0 ;
                                d2psi_yz[i] = 0 ;
                                d2psi_zz[i] = 0 ;
	    			break ;
	    		}
	    		
	    		case 1 : {
	    		
		r2p[0] = aa2 ;    	 // 3 points outside the throat
		r2p[1] = 1.1 * aa2 ;     // for the parabolic extrapolation
		r2p[2] = 1.2 * aa2 ;     //
		
		double rot_phi = mp2.get_rot_phi() ;
		double orix = mp2.get_ori_x() ;
		double oriy = mp2.get_ori_y() ;
		double oriz = mp2.get_ori_z() ;
		
		for (int j=0; j<3; j++) {
			double phi = phi2 + rot_phi ;
			double xap = r2p[j] * sin(theta2) * cos(phi) + orix ;
			double yap = r2p[j] * sin(theta2) * sin(phi) + oriy ;
			double zap = r2p[j] * cos(theta2) + oriz ;
			mp1.convert_absolute(xap, yap, zap,
						r1p[j], t1p[j], p1p[j]) ;
		}

		// Lapse function
		for (int j=0; j<3; j++) {
			yp[j] =   cnn1.val_point(r1p[j], t1p[j], p1p[j])
				+ cnn2.val_point(r2p[j], theta2, phi2)  ;
                }			
                nnn[i] = lagrange_parabol(r2, r2p, yp) ;
                	    		
		// Shift vector
		for (int j=0; j<3; j++) {
			yp[j] =   cbetax1.val_point(r1p[j], t1p[j], p1p[j])
				- cbetax2.val_point(r2p[j], theta2, phi2) ;
		}
  		beta_x[i] = lagrange_parabol(r2, r2p, yp) - omega * yy[i] ;

		for (int j=0; j<3; j++) {
			yp[j] =   cbetay1.val_point(r1p[j], t1p[j], p1p[j])
				- cbetay2.val_point(r2p[j], theta2, phi2) ;
		}
  		beta_y[i] = lagrange_parabol(r2, r2p, yp) + omega * xx[i] ;
		
		for (int j=0; j<3; j++) {
			yp[j] =   cbetaz1.val_point(r1p[j], t1p[j], p1p[j])
				+ cbetaz2.val_point(r2p[j], theta2, phi2) ;
		}
  		beta_z[i] = lagrange_parabol(r2, r2p, yp) ;

		// 3-metric
		for (int j=0; j<3; j++) {
			yp[j] =   cpsi1.val_point(r1p[j], t1p[j], p1p[j])
				+ cpsi2.val_point(r2p[j], theta2, phi2)  ;
                }			
                double psi4 = pow( lagrange_parabol(r2, r2p, yp), 4) ;
		g_xx[i] = psi4 ;
		g_yy[i] = psi4 ; 	
		g_zz[i] = psi4 ; 	
		g_xy[i] = 0 ;
		g_xz[i] = 0 ;
		g_yz[i] = 0 ;

 		// Extrinsic curvature
		double pre = aa * psi4 ;
		for (int j=0; j<3; j++) {
			yp[j] =   ckxx1.val_point(r1p[j], t1p[j], p1p[j])
				+ ckxx2.val_point(r2p[j], theta2, phi2)  ;
                }			
                k_xx[i] = pre * lagrange_parabol(r2, r2p, yp) ;

		for (int j=0; j<3; j++) {
			yp[j] =   ckxy1.val_point(r1p[j], t1p[j], p1p[j])
				+ ckxy2.val_point(r2p[j], theta2, phi2)  ;
                }			
                k_xy[i] = pre * lagrange_parabol(r2, r2p, yp) ;

		for (int j=0; j<3; j++) {
			yp[j] =   ckxz1.val_point(r1p[j], t1p[j], p1p[j])
				- ckxz2.val_point(r2p[j], theta2, phi2)  ;
                }
                k_xz[i] = pre * lagrange_parabol(r2, r2p, yp) ;

		for (int j=0; j<3; j++) {
			yp[j] =   ckyy1.val_point(r1p[j], t1p[j], p1p[j])
				+ ckyy2.val_point(r2p[j], theta2, phi2)  ;
                }			
                k_yy[i] = pre * lagrange_parabol(r2, r2p, yp) ;

		for (int j=0; j<3; j++) {
			yp[j] =   ckyz1.val_point(r1p[j], t1p[j], p1p[j])
				- ckyz2.val_point(r2p[j], theta2, phi2)  ;
                }			
                k_yz[i] = pre * lagrange_parabol(r2, r2p, yp) ;

		for (int j=0; j<3; j++) {
			yp[j] =   ckzz1.val_point(r1p[j], t1p[j], p1p[j])
				+ ckzz2.val_point(r2p[j], theta2, phi2)  ;
                }
                k_zz[i] = pre * lagrange_parabol(r2, r2p, yp) ;

                // First derivatives of Psi
                for (int j=0; j<3; j++) {
		        yp[j] = aa * (  cdpsix1.val_point(r1p[j], t1p[j], p1p[j])
				- cdpsix2.val_point(r2p[j], theta2, phi2) ) ;
		}
  		dpsi_x[i] = lagrange_parabol(r2, r2p, yp)  ;

                for (int j=0; j<3; j++) {
		        yp[j] = aa * (   cdpsiy1.val_point(r1p[j], t1p[j], p1p[j])
				- cdpsiy2.val_point(r2p[j], theta2, phi2) ) ;
		}
  		dpsi_y[i] = lagrange_parabol(r2, r2p, yp)  ;

                for (int j=0; j<3; j++) {
		        yp[j] = aa * (   cdpsiz1.val_point(r1p[j], t1p[j], p1p[j])
				+ cdpsiz2.val_point(r2p[j], theta2, phi2) ) ;
		}
  		dpsi_z[i] = lagrange_parabol(r2, r2p, yp)  ;

                // Second derivatives of Psi
                for (int j=0; j<3; j++) {
		        yp[j] = aasq * (   cd2psixx1.val_point(r1p[j], t1p[j], p1p[j])
				+ cd2psixx2.val_point(r2p[j], theta2, phi2) ) ;
		}
  		d2psi_xx[i] = lagrange_parabol(r2, r2p, yp)  ;

                for (int j=0; j<3; j++) {
		        yp[j] = aasq * (    cd2psixy1.val_point(r1p[j], t1p[j], p1p[j])
				+ cd2psixy2.val_point(r2p[j], theta2, phi2) ) ;
		}
  		d2psi_xy[i] = lagrange_parabol(r2, r2p, yp)  ;

                for (int j=0; j<3; j++) {
		        yp[j] =  aasq * (  cd2psixz1.val_point(r1p[j], t1p[j], p1p[j])
				- cd2psixz2.val_point(r2p[j], theta2, phi2) ) ;
		}
  		d2psi_xz[i] = lagrange_parabol(r2, r2p, yp)  ;

                for (int j=0; j<3; j++) {
		        yp[j] = aasq * (   cd2psiyy1.val_point(r1p[j], t1p[j], p1p[j])
				+ cd2psiyy2.val_point(r2p[j], theta2, phi2) ) ;
		}
  		d2psi_yy[i] = lagrange_parabol(r2, r2p, yp)  ;

                for (int j=0; j<3; j++) {
		        yp[j] = aasq * (   cd2psiyz1.val_point(r1p[j], t1p[j], p1p[j])
				- cd2psiyz2.val_point(r2p[j], theta2, phi2) ) ;
		}
  		d2psi_yz[i] = lagrange_parabol(r2, r2p, yp)  ;

                for (int j=0; j<3; j++) {
		        yp[j] = aasq * (   cd2psizz1.val_point(r1p[j], t1p[j], p1p[j])
				+ cd2psizz2.val_point(r2p[j], theta2, phi2) ) ;
		}
  		d2psi_zz[i] = lagrange_parabol(r2, r2p, yp)  ;


		break ;	    		
	    		}	// end of case fill = 1
	    		
	    		default : {
	    		
	  cout << "Bin_BH::Bin_BH : the case fill = " << fill
	  	<< " is not known !" << endl ;
	  abort() ;
	  break ;
	    		}
		 }
	} // End of case r2 < a
	
	
	
	//------------------------------------------------------------
	// 			Outside the two holes
	//------------------------------------------------------------	
		
	if ( (r1 >= aa) && (r2 >= aa2) ) {
	
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
	g_yy[i] = psi4 ; 	
	g_zz[i] = psi4 ; 	
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
		      	

        // First derviatives of psi
        // -----------------------

	dpsi_x[i] = vdpsix1.c_cf->val_point_symy(l1, xi1, theta1, phi1)
		 -  vdpsix2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ;

	dpsi_y[i] = vdpsiy1.c_cf->val_point_asymy(l1, xi1, theta1, phi1)
		 -  vdpsiy2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ;

	dpsi_z[i] = vdpsiz1.c_cf->val_point_symy(l1, xi1, theta1, phi1)
		 +  vdpsiz2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ;

        // Second derviatives of psi
        // -------------------------

	d2psi_xx[i] = vd2psixx1.c_cf->val_point_symy(l1, xi1, theta1, phi1)
		    + vd2psixx2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ;

	d2psi_xy[i] = vd2psixy1.c_cf->val_point_asymy(l1, xi1, theta1, phi1)
		    + vd2psixy2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ;

	d2psi_xz[i] = vd2psixz1.c_cf->val_point_symy(l1, xi1, theta1, phi1)
		    - vd2psixz2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ;

	d2psi_yy[i] = vd2psiyy1.c_cf->val_point_symy(l1, xi1, theta1, phi1)
		    + vd2psiyy2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ;

	d2psi_yz[i] = vd2psiyz1.c_cf->val_point_asymy(l1, xi1, theta1, phi1)
		    - vd2psiyz2.c_cf->val_point_asymy(l2, xi2, theta2, phi2) ;

	d2psi_zz[i] = vd2psizz1.c_cf->val_point_symy(l1, xi1, theta1, phi1)
		    + vd2psizz2.c_cf->val_point_symy(l2, xi2, theta2, phi2) ;



        }


    }	// End of loop on the points
    


}

		

//=========================================================================
// Interpolation/extrapolation by means of a parabola connecting
//  three points
//=========================================================================

double lagrange_parabol(double x, const double* xp, const double* yp) {

	double a0 = yp[0] / (xp[0] - xp[1]) / (xp[0] - xp[2]) ;
	double a1 = yp[1] / (xp[1] - xp[0]) / (xp[1] - xp[2]) ;
	double a2 = yp[2] / (xp[2] - xp[0]) / (xp[2] - xp[1]) ;

 	return 		(x - xp[1]) * (x - xp[2]) * a0
     	            +   (x - xp[0]) * (x - xp[2]) * a1
         	    +   (x - xp[0]) * (x - xp[1]) * a2 ;	
}


