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
 * Revision 1.15  2014/10/13 08:52:56  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.14  2007/10/18 14:26:43  e_gourgoulhon
 * Changed the call to Eos::der_nbar_ent in order to allow for MEos type
 * of equation of state.
 *
 * Revision 1.13  2007/10/16 21:56:26  e_gourgoulhon
 * Can deal with more than one domain into the star,
 * thanks to the new function Map_radial::poisson_compact.
 *
 * Revision 1.12  2005/10/18 13:12:33  p_grandclement
 * update of the mixted binary codes
 *
 * Revision 1.11  2004/05/25 15:38:38  f_limousin
 * Minor modifs.
 *
 * Revision 1.10  2004/05/10 10:17:27  f_limousin
 * Add a new member ssjm1_psi of class Etoile for the resolution of the
 * oisson_interne equation
 *
 * Revision 1.9  2004/04/19 11:26:17  f_limousin
 * Add a new function Etoile_bin::velocity_potential( , , , ) for the
 * case of strange stars
 *
 * Revision 1.8  2004/04/08 17:02:00  f_limousin
 * Modif to avoid an error in the compilation
 *
 * Revision 1.7  2004/04/08 16:52:58  f_limousin
 * Minor change
 *
 * Revision 1.6  2004/04/08 16:36:36  f_limousin
 * Implement the resolution of the continuity equation for strange
 * stars.
 *
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
#include "scalar.h"
#include "metrique.h" 
#include "etoile.h"
#include "eos.h"
#include "param.h"
#include "et_bin_nsbh.h"
#include "utilitaires.h"

// Local prototype
namespace Lorene {
Cmp raccord_c1(const Cmp& uu, int l1) ; 

double Etoile_bin::velocity_potential(int mermax, double precis, double relax) {
  
  // Which star is that ?
  const Et_bin_nsbh* pnsbh = dynamic_cast<const Et_bin_nsbh*>(this) ;
  
  if (eos.identify() == 5 || eos.identify() == 4 || 
      eos.identify() == 3) {
    
    // Routine used for binary strange stars.

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

    v_orb.annule(nzm1, nzm1) ;	// set to zero in the ZEC


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
    
    double erreur ;				
   
    Tenseur zeta_h( ent() / dndh_log ) ;
    zeta_h.set_std_base() ;
    
    Scalar zeta_h_scalar (zeta_h()) ;
    zeta_h_scalar.set_outer_boundary(0, (ent() / dndh_log)(0,0,0,0)) ;
    for (int l=1; l<=nzm1; l++)
	zeta_h_scalar.set_domain(l) = 1 ;
	
    Cmp zeta_h_cmp (zeta_h_scalar) ;
    zeta_h.set() = zeta_h_cmp ;
    zeta_h.set_std_base() ;

    

    Tenseur beta(mp) ;
    
    if (pnsbh!=0x0) {
	beta = log( sqrt(a_car) * nnn ) ; 
	beta.set_std_base() ;
    }
    else {
	beta = beta_auto + beta_comp ; 
    }

    Tenseur tmp_zeta = 1 - unsurc2 * zeta_h ;
    tmp_zeta.set_std_base() ;
	
    Tenseur bb = tmp_zeta * ent.gradient_spher()
	+ unsurc2 * zeta_h * beta.gradient_spher() ;
	
    Tenseur entmb = ent - beta ; 
	
    Tenseur grad_ent (ent.gradient()) ;
    grad_ent.change_triad(mp.get_bvect_spher()) ;

    // Source for the poisson equation 
    // See Eq (63) from Gourgoulhon et al. (2001)
    Tenseur source = flat_scalar_prod( www - v_orb, ent.gradient() )
	+ unsurc2 * zeta_h * (
	    flat_scalar_prod( v_orb, entmb.gradient() )
	    + flat_scalar_prod( www, gam_euler.gradient() )
	    / gam_euler ) ; 

    for (int l=1; l<=nzm1; l++)
	source.set().annule(l) ;
   
    source = (source - flat_scalar_prod(bb, psi0.gradient_spher())) 
	/ zeta_h ;
    source.annule(nzet, nzm1) ; 

    Param par ; 
    int niter ; 
    par.add_int(mermax) ; 
    par.add_double(precis, 0) ; 
    par.add_double(relax, 1) ; 
    par.add_int_mod(niter) ; 

    par.add_cmp_mod(ssjm1_psi, 0) ;
	
    if (psi0.get_etat() == ETATZERO) {
	psi0.set_etat_qcq() ; 
	psi0.set() = 0 ; 
    }

    int nr = mp.get_mg()->get_nr(0);
    int nt = mp.get_mg()->get_nt(0);
    int np = mp.get_mg()->get_np(0);

    cout << "nr = " << nr << "   nt = " << nt << "   np = " << np << endl ;

    cout << "psi0" << endl << norme(psi0()/(nr*nt*np)) << endl ;
    cout << "d(psi)/dr" << endl << norme(psi0.set().dsdr()/(nr*nt*np)) << endl ;

    Valeur lim(mp.get_mg()->get_angu()) ;
    lim.annule_hard() ;

    Tenseur normal (mp, 1, CON, mp.get_bvect_cart()) ;
    Tenseur normal2 (mp, 1, COV, mp.get_bvect_cart()) ;
    normal.set_etat_qcq() ;
    normal2.set_etat_qcq() ;

    const Coord& rr0 = mp.r ;
    Tenseur rr(mp) ;
    rr.set_etat_qcq() ;
    rr.set() = rr0 ;
    rr.set_std_base() ;

    Tenseur_sym plat(mp, 2, COV, mp.get_bvect_cart() ) ;
    plat.set_etat_qcq() ;
    for (int i=0; i<3; i++) {
	for (int j=0; j<i; j++) {
	    plat.set(i,j) = 0 ;
	}
	plat.set(i,i) = 1 ;
    }
    plat.set_std_base() ;
	
    Metrique flat(plat, true) ; 
    Tenseur dcov_r = rr.derive_cov(flat) ;


    for (int i=0; i<3; i++) {
	normal.set(i) = dcov_r(i) ;
	normal2.set(i) = dcov_r(i) ;
    }

    normal.change_triad(mp.get_bvect_spher()) ;
    normal2.change_triad(mp.get_bvect_spher()) ;
 


    Tenseur bsn0 (bsn) ;
    bsn0.change_triad(mp.get_bvect_cart()) ;
    Tenseur aa (mp, 1, CON, mp.get_bvect_cart()) ;
    aa = - v_orb - a_car * gam_euler * hhh * bsn0 ;
    aa.change_triad(mp.get_bvect_spher()) ;
	

    Tenseur dcov_psi = psi0.derive_cov(flat) ;
    dcov_psi.change_triad(mp.get_bvect_spher()) ;

    Cmp limite (mp) ;
    limite =  ( - dcov_psi(1) * normal(1) - dcov_psi(2) * normal(2)
		+ contract(aa, 0, normal2, 0)()) 
	/normal(0) ;
	
    for (int j=0; j<nt; j++)
	for (int k=0; k<np; k++)
	    lim.set(0, k, j, 0) = limite(0, k, j, nr-1) ;
	
//	cout << "lim" << endl << lim << endl ;

    lim.std_base_scal() ;
    Cmp resu (psi0()) ;
    source().poisson_neumann_interne(lim, par, resu) ;
    psi0 = resu ;

/*
    resu.va.ylm() ;
    Scalar psi00(resu) ;
    psi00.spectral_display("psi00") ;

    cout << "value of d(psi)/dr at the surface after poisson" << endl ;
    for (int j=0; j<nt; j++)
    for (int k=0; k<np; k++)
    cout << "j = " << j << " ; k = " << k << " : " << 
    psi0.set().dsdr()(0, k, j, nr-1) << endl ;
*/
    for (int l=1; l<=nzm1; l++)
	psi0.set().annule(l) ;
    

    //---------------------------------------------------
    // Check of the solution  
    //---------------------------------------------------
	
    Cmp laplacien_psi0 =  psi0().laplacien() ; 
	
    erreur = diffrel(laplacien_psi0, source())(0) ; 

    cout << "Check of the resolution of the continuity equation for strange stars: " 
	 << endl ; 
    cout << "norme(source) : " << norme(source())(0) << endl 
	 << "Error in the solution : " << erreur << endl ; 
	
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
 

  }  // End of strange stars case

//=============================================================================

  else {
    
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
    
    Cmp dndh_log(mp) ; 
	dndh_log = 0 ; 

	for (int l=0; l<nzet; l++) {

        Param par ;       // Paramater for multi-domain equation of state
        par.add_int(l) ;

        dndh_log = dndh_log +  eos.der_nbar_ent(ent(), 1, l, &par) ;

    }
	
    // Cmp dndh_log = eos.der_nbar_ent(ent(), nzet) ; 

    // In order to avoid any division by zero in the computation of zeta_h
    //  the value of dndh_log is set to 1 in the external domains:
    for (int l=nzet; l <= nzm1; l++) {
	dndh_log.set(l) = 1 ; 
    }
    
    Tenseur zeta_h( ent() / dndh_log ) ;
    zeta_h.set_std_base() ;
    
    Tenseur beta(mp) ; 
    
    if (pnsbh!=0x0) {
	beta = log( sqrt(a_car) * nnn ) ;
	beta.set_std_base() ;
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
    
    mp.poisson_compact(nzet, source(), zeta_h(), bb, par, psi0.set() ) ;

    //---------------------------------------------------
    // Check of the solution  
    //---------------------------------------------------
	
    Tenseur bb_dpsi0 = flat_scalar_prod( bb, psi0.gradient_spher() ) ;
	
    Cmp oper = zeta_h() * psi0().laplacien() + bb_dpsi0() ; 
	
    source.set().va.ylm_i() ;
	
    cout << "Check of the resolution of the continuity equation : "  << endl ; 
    Tbl terr = diffrel(oper, source()) ;
	double erreur = 0 ; 
	for (int l=0; l<nzet; l++) { 
		double err = terr(l) ; 
    	cout << " domain " << l << " : norme(source) : " << norme(source())(l) 
	 		<< "    relative error : " << err << endl ;
		if (err > erreur) erreur = err ; 
	} 
	// arrete() ; 
    
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
}
}
