/*
 * Method of class Star_bin to compute the velocity scalar potential $\psi$
 * by solving the continuity equation.
 *
 * (see file star.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004      Francois Limousin
 *                 2000-2001 Eric Gourgoulhon  (for precedent version)
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

char star_bin_vel_pot_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/06/22 12:53:09  f_limousin
 * Change qq, qq_auto and qq_comp to beta, beta_auto and beta_comp.
 *
 * Revision 1.2  2004/06/07 16:26:01  f_limousin
 * Many modif...
 *
 * Revision 1.1  2004/04/08 16:34:08  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "star.h"
#include "eos.h"
#include "param.h"
#include "cmp.h"
#include "tenseur.h"
#include "graphique.h"
#include "utilitaires.h"

// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ; 

double Star_bin::velocity_potential(int mermax, double precis, double relax) {
  
    int nzm1 = mp.get_mg()->get_nzone() - 1 ;

     //----------------------------------
    // Specific relativistic enthalpy		    ---> hhh
    //----------------------------------
    
    Scalar hhh = exp(ent) ;  // = 1 at the Newtonian limit
    hhh.std_spectral_base() ;

    //----------------------------------------------
    //  Computation of W^i = - A^2 h Gamma_n B^i/N
    // See Eq (62) from Gourgoulhon et al. (2001)
    //----------------------------------------------

    Vector www = - hhh * gam_euler * bsn * psi4 ;
	
    www.change_triad( mp.get_bvect_cart() ) ;	// components on the mapping
						// Cartesian basis
    
    //-------------------------------------------------
    // Constant value of W^i at the center of the star
    //-------------------------------------------------
    
    Vector v_orb(mp, CON, mp.get_bvect_cart()) ; 
    
    v_orb.set_etat_qcq() ;
    for (int i=1; i<=3; i++) {
	v_orb.set(i) = www(i).val_grid_point(0, 0, 0, 0) ;
    }

    v_orb.set_triad( *(www.get_triad()) ) ;  
    v_orb.std_spectral_base() ;
   
    //-------------------------------------------------
    // Source and coefficients a,b for poisson_compact (independent from psi0)
    //-------------------------------------------------
    
    Scalar dndh_log = eos.der_nbar_ent(ent, nzet) ; 

    // In order to avoid any division by zero in the computation of zeta_h
    //  the value of dndh_log is set to 1 in the external domains:
    for (int l=nzet; l <= nzm1; l++) {
	dndh_log.set_domain(l) = 1 ; 
    }
    
    Scalar zeta_h( ent / dndh_log ) ;
    zeta_h.std_spectral_base() ;

    Vector bb = (1 - zeta_h) * ent.derive_con(flat) + 
	zeta_h * beta.derive_con(flat) ;

    Scalar entmb = ent - beta ;  

    www.change_triad(mp.get_bvect_spher()) ;
    v_orb.change_triad(mp.get_bvect_spher()) ;

    // See Eq (63) from Gourgoulhon et al. (2001)
    Scalar source = contract(www - v_orb, 0, ent.derive_cov(flat), 0)
	+ zeta_h * ( contract(v_orb, 0, entmb.derive_cov(flat), 0) +
		     contract(www/gam_euler, 0, gam_euler.derive_cov(flat), 0)
		     + contract(hij, 0, 1, (psi0.derive_cov(flat)  
				+ v_orb.down(0, flat))*entmb.derive_cov(flat),
				0, 1) )
	- contract(hij, 0, 1, ent.derive_cov(flat) * (psi0.derive_cov(flat) 
						      + v_orb.down(0, flat)), 
						      0, 1) ;
  
/*
    des_meridian(zeta_h,0., 4., "zeta_h", 10) ; 
    arrete() ; 
    des_meridian(bb(1),0., 4., "bb(1)", 10) ; 
    arrete() ; 
    des_meridian(bb(2),0., 4., "bb(2)", 10) ; 
    arrete() ; 
    des_meridian(bb(3),0., 4., "bb(3)", 10) ; 
    arrete() ; 
    des_meridian(psi0,0., 4., "psi0", 10) ; 
    arrete() ; 
    des_meridian(source,0., 4., "source", 10) ; 
    arrete() ; 
*/  


    www.change_triad(mp.get_bvect_cart()) ;
    v_orb.change_triad(mp.get_bvect_cart()) ;
				 
    source.annule(nzet, nzm1) ; 
    			
    //---------------------------------------------------
    // Resolution by means of Map_radial::poisson_compact 
    //---------------------------------------------------

    Param par ; 
    int niter ; 
    par.add_int(mermax) ; 
    par.add_double(precis, 0) ; 
    par.add_double(relax, 1) ; 
    par.add_int_mod(niter) ; 
    
    
    if (psi0.get_etat() == ETATZERO) {
	psi0 = 0 ; 
    }

    Cmp source_cmp (source) ;
    Cmp zeta_h_cmp (zeta_h) ;
    Cmp psi0_cmp (psi0) ;
    
    Tenseur bb_cmp(mp, 1, CON, mp.get_bvect_spher()) ;
    bb_cmp.set_etat_qcq() ;
    Cmp bb_cmp1 (bb(1)) ;
    Cmp bb_cmp2 (bb(2)) ;
    Cmp bb_cmp3 (bb(3)) ;
    bb_cmp.set(0) = bb_cmp1 ;
    bb_cmp.set(1) = bb_cmp2 ;
    bb_cmp.set(2) = bb_cmp3 ;

    source_cmp.va.ylm() ;

    cout << "source" << endl << norme(source_cmp)<< endl ;
    cout << "zeta_h " << endl << norme(zeta_h_cmp) << endl ;
    cout << "bb(1)" << endl << norme(bb_cmp(0)) << endl ;
    cout << "bb(2)" << endl << norme(bb_cmp(1)) << endl ;
    cout << "bb(3)" << endl << norme(bb_cmp(2)) << endl ;
    cout << "psiO" << endl << norme(psi0_cmp) << endl ;

    mp.poisson_compact(source_cmp, zeta_h_cmp, bb_cmp, par, psi0_cmp ) ;
    
    psi0 = psi0_cmp ;

    cout << "psiO apres" << endl << norme(psi0) << endl ;

    //---------------------------------------------------
    // Check of the solution  
    //---------------------------------------------------
    
    Scalar bb_dpsi0 = contract(bb, 0, psi0.derive_cov(flat), 0) ;
    
    Scalar oper = zeta_h * psi0.laplacian() + bb_dpsi0 ; 
    
    source.set_spectral_va().ylm_i() ;

    double erreur = diffrel(oper, source)(0) ; 

    cout << "Check of the resolution of the continuity equation : " 
	 << endl ; 
    cout << "            norme(source) : " << norme(source)(0) 
         << "    diff oper/source : " << erreur << endl ; 
    
    
    //--------------------------------
    // Computation of grad(psi)
    //--------------------------------
    
    v_orb.change_triad(mp.get_bvect_spher()) ;
    d_psi.change_triad(mp.get_bvect_spher()) ;

    for (int i=1; i<=3; i++) 
	d_psi.set(i) = (psi0.derive_cov(flat))(i) + v_orb(i) ; 

    v_orb.change_triad(mp.get_bvect_cart()) ;    
    d_psi.change_triad(mp.get_bvect_cart()) ; 
   
    // C^1 continuation of d_psi outside the star
    //  (to ensure a smooth enthalpy field accross the stellar surface)
    // ----------------------------------------------------------------
    
    d_psi.annule(nzet, nzm1) ;	  
 
   
    for (int i=1; i<=3; i++) {
	Cmp d_psi_i (d_psi(i)) ;
	d_psi_i = raccord_c1(d_psi_i, nzet) ; 
	d_psi.set(i) = d_psi_i ; 
    }
    
    return erreur ; 

}
