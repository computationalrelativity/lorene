/*
 * Method of class Etoile_bin to compute the velocity scalar potential $\psi$
 * by solving the continuity equation.
 *
 * (see file etoile.h for documentation).
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char et_bin_vel_pot_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2003/10/24 11:43:57  e_gourgoulhon
 * beta is now computed as ln(AN) in the case beta_auto
 * is undefined (for instance, if the companion is black hole).
 *
 * Revision 1.4  2003/01/17 13:38:56  f_limousin
 * Add comments
 *
 * Revision 1.3  2003/01/13 15:31:50  e_gourgoulhon
 * Suppressed the desaliasing
 *  (did not worked due to missing basis in ylm).
 *
 * Revision 1.2  2002/12/10 14:44:21  k_taniguchi
 * Change the multiplication "*" to "%"
 *   and flat_scalar_prod to flat_scalar_prod_desal.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.9  2001/02/23  15:18:59  eric
 * Modification du calcul de zeta_h pour eviter division par zero
 *   dans les domaines externes a l'etoile.
 *
 * Revision 2.8  2001/02/07  09:47:42  eric
 * zeta_h est desormais donne par Eos::der_nbar_ent.
 *
 * Revision 2.7  2000/12/22  13:10:03  eric
 * Prolongement C^1 de dpsi en dehors de l'etoile.
 *
 * Revision 2.6  2000/03/22  12:56:44  eric
 * Nouveau prototype d'Etoile_bin::velocity_potential : l'erreur est
 * retournee en double.
 *
 * Revision 2.5  2000/02/25  17:35:29  eric
 * Annulation de la source dans les zones externes avant l'appel a
 * poisson_compact.
 *
 * Revision 2.4  2000/02/22  11:42:55  eric
 * Test resolution de l'equation.
 *
 * Revision 2.3  2000/02/22  10:42:25  eric
 * Correction erreur dans les termes sources: multiplication par unsurc2 de
 *  termes relativistes.
 *
 * Revision 2.2  2000/02/21  15:05:50  eric
 * Traitement du cas psi0 = 0 .
 *
 * Revision 2.1  2000/02/21  13:59:39  eric
 * Remplacement du membre psi par psi0.
 * Modif calcul de d_psi a la fin.
 *
 * Revision 2.0  2000/02/17  18:50:44  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "etoile.h"
#include "eos.h"
#include "param.h"

// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ; 

double Etoile_bin::velocity_potential(int mermax, double precis, double relax) {
    
    int nzm1 = mp.get_mg()->get_nzone() - 1 ;    

    //----------------------------------
    // Specific relativistic enthalpy		    ---> hhh
    //----------------------------------
    
    Tenseur hhh = exp(unsurc2 * ent) ;  // = 1 at the Newtonian limit
    hhh.set_std_base() ;

    //----------------------------------------------
    //  Computation of W^i = - A^2 h Gamma_n B^i/N
    // See Eq (62) from Gourgoulhon et al. (2001)
    //----------------------------------------------

    Tenseur www = - a_car * hhh * gam_euler * bsn ; 
    
    www.change_triad( mp.get_bvect_cart() ) ;	// components on the mapping
						// Cartesian basis
    
    //-------------------------------------------------
    // Constant value of W^i at the center of the star
    //-------------------------------------------------
    
    Tenseur v_orb(mp, 1, CON, mp.get_bvect_cart()) ; 
    
    v_orb.set_etat_qcq() ; 
    for (int i=0; i<3; i++) {
	v_orb.set(i) = www(i)(0, 0, 0, 0) ; 
    }

    v_orb.set_triad( *(www.get_triad()) ) ;  
    v_orb.set_std_base() ;
   
    //-------------------------------------------------
    // Source and coefficients a,b for poisson_compact (idenpendent from psi0)
    //-------------------------------------------------
    
    Cmp dndh_log = eos.der_nbar_ent(ent(), nzet) ; 

    // In order to avoid any division by zero in the computation of zeta_h
    //  the value of dndh_log is set to 1 in the external domains:
    for (int l=nzet; l <= nzm1; l++) {
	dndh_log.set(l) = 1 ; 
    }
    
    Tenseur zeta_h( ent() / dndh_log ) ;
    zeta_h.set_std_base() ;

	Tenseur beta(mp) ; 

	if (beta_auto.get_etat() == ETATNONDEF) {
		beta = log( sqrt(a_car) * nnn ) ; 
	}
	else {
		beta = beta_auto + beta_comp ; 
	}

    Tenseur tmp_zeta = 1 - unsurc2 * zeta_h ;
    tmp_zeta.set_std_base() ;

    Tenseur bb = tmp_zeta * ent.gradient_spher()
		    + unsurc2 * zeta_h * beta.gradient_spher() ;
		    
    Tenseur entmb = ent - beta ; 
    
    // See Eq (63) from Gourgoulhon et al. (2001)
    Tenseur source = flat_scalar_prod( www - v_orb, ent.gradient() )
		     + unsurc2 * zeta_h * (
			 flat_scalar_prod( v_orb, entmb.gradient() )
		       + flat_scalar_prod( www, gam_euler.gradient() )
		         / gam_euler ) ; 
				 
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
	psi0.set_etat_qcq() ; 
	psi0.set() = 0 ; 
    }

    source.set().va.ylm() ;

    mp.poisson_compact(source(), zeta_h(), bb, par, psi0.set() ) ;
    
    //---------------------------------------------------
    // Check of the solution  
    //---------------------------------------------------
    
    Tenseur bb_dpsi0 = flat_scalar_prod( bb, psi0.gradient_spher() ) ;
    
    Cmp oper = zeta_h() * psi0().laplacien() + bb_dpsi0() ; 
    
    source.set().va.ylm_i() ;

    double erreur = diffrel(oper, source())(0) ; 

    cout << "Check of the resolution of the continuity equation : " 
	 << endl ; 
    cout << "            norme(source) : " << norme(source())(0) 
         << "    diff oper/source : " << erreur << endl ; 
    
    
    //--------------------------------
    // Computation of grad(psi)
    //--------------------------------
    
    // The computation is done component by component because psi0.gradient()
    // is a covariant vector, whereas v_orb is a contravariant one. 
    
    d_psi.set_etat_qcq() ; 
    
    for (int i=0; i<3; i++) {
	d_psi.set(i) = (psi0.gradient())(i) + v_orb(i) ; 
    }

    d_psi.set_triad(  *(v_orb.get_triad()) ) ; 
   
    // C^1 continuation of d_psi outside the star
    //  (to ensure a smooth enthalpy field accross the stellar surface)
    // ----------------------------------------------------------------
    
    d_psi.annule(nzet, nzm1) ;	  
    for (int i=0; i<3; i++) {
	d_psi.set(i) = raccord_c1(d_psi(i), nzet) ; 
    }
    
    
    assert( d_psi.get_triad() == &(mp.get_bvect_cart()) ) ; 

    d_psi.change_triad(ref_triad) ; 
    
    return erreur ; 
 
}
