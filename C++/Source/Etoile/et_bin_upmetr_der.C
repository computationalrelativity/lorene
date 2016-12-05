/*
 * Methods Etoile_bin::update_metric_der_comp
 *
 * (see file etoile.h for documentation)
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.8  2016/12/05 16:17:53  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2014/10/13 08:52:56  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2003/10/24 12:27:16  k_taniguchi
 * Suppress the method of update metric for NS-BH
 *
 * Revision 1.5  2003/10/24 11:47:02  k_taniguchi
 * Change some notations
 *
 * Revision 1.4  2002/12/19 14:53:38  e_gourgoulhon
 * Added the new function
 * 	void update_metric_der_comp(const Bhole& comp)
 * to treat the case where the companion is a black hole
 *
 * Revision 1.3  2002/12/10 15:12:07  k_taniguchi
 * Change the multiplication "*" to "%".
 *
 * Revision 1.2  2001/12/14 15:15:30  k_taniguchi
 * Change of the method to calculate derivatives with respect to the companion star
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/03/13  14:03:38  eric
 * Modif commentaires.
 *
 * Revision 2.3  2000/03/07  14:54:54  eric
 * Ajout du calcul de akcar_comp.
 *
 * Revision 2.2  2000/03/07  08:34:04  eric
 * Appel de Cmp::import_sym / asym (pour tenir compte de la symetrie /
 *  plan y=0).
 *
 * Revision 2.1  2000/02/10  18:56:38  eric
 * Traitement du cas ETATZERO.
 *
 * Revision 2.0  2000/02/04  16:38:11  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "etoile.h"
#include "bhole.h"

namespace Lorene {
void Etoile_bin::update_metric_der_comp(const Etoile_bin& comp) {

    // Computation of d_logn_comp
    // --------------------------

    if ( (comp.d_logn_auto).get_etat() == ETATZERO ) {
	d_logn_comp.set_etat_zero() ;
    }
    else{
      d_logn_comp = logn_comp.gradient() ;
    }

    d_logn_comp.change_triad(ref_triad) ;

    // Computation of d_beta_comp
    // --------------------------

    if ( (comp.d_beta_auto).get_etat() == ETATZERO ) {
	d_beta_comp.set_etat_zero() ;
    }
    else {
      d_beta_comp = beta_comp.gradient() ;
    }

    d_beta_comp.change_triad(ref_triad) ;

    // Computation of tkij_comp
    // ------------------------

    if ( (comp.tkij_auto).get_etat() == ETATZERO ) {
	tkij_comp.set_etat_zero() ;
    }
    else{

      // Components of shift_comp with respect to the Cartesian triad
      //  (d/dx, d/dy, d/dz) of the mapping :
      Tenseur shift_comp_local = shift_comp ;
      shift_comp_local.change_triad( mp.get_bvect_cart() ) ;

      // Gradient (partial derivatives with respect to
      //           the Cartesian coordinates of the mapping)
      // D_j N^i

      Tenseur dn_comp = shift_comp_local.gradient() ;

      // Return to the absolute reference frame
      dn_comp.change_triad(ref_triad) ;

      // Trace of D_j N^i = divergence of N^i :
      Tenseur divn_comp = contract(dn_comp, 0, 1) ;

      // Computation of A^2 K^{ij}
      // -------------------------
      tkij_comp.set_etat_qcq() ;

      for (int i=0; i<3; i++) {
	for (int j=i; j<3; j++) {
	  tkij_comp.set(i, j) = dn_comp(i, j) + dn_comp(j, i)  ;
	}
	tkij_comp.set(i, i) -= double(2)/double(3) * divn_comp() ;
      }

      tkij_comp = - 0.5 * tkij_comp / nnn ;

    }

    tkij_comp.set_triad( *((comp.tkij_auto).get_triad()) ) ;
    tkij_comp.set_std_base() ;

    if (relativistic) {
	// Computation of akcar_comp
	// -------------------------
    
	akcar_comp.set_etat_qcq() ;
    
	akcar_comp.set() = 0 ;

	for (int i=0; i<3; i++) {
	    for (int j=0; j<3; j++) {

		akcar_comp.set() += tkij_auto(i, j) % tkij_comp(i, j) ;

	    }
	}

	akcar_comp.set_std_base() ;
	akcar_comp = a_car % akcar_comp ;

    }

    // The derived quantities are obsolete
    // -----------------------------------

    del_deriv() ;


    //-----------------------------------------------------
    // The previous way to calculate d_logn_comp and so on
    //  which we do not use
    //-----------------------------------------------------

    //#################################
    /*
    int nz = mp.get_mg()->get_nzone() ;
    int nzm1 = nz - 1 ;

    // Computation of d_logn_comp
    // --------------------------

    if ( (comp.d_logn_auto).get_etat() == ETATZERO ) {
	d_logn_comp.set_etat_zero() ;
    }
    else{

	// 1/ Division by r^2 of comp.d_logn_auto in the ZEC
	Tenseur vecttmp = comp.d_logn_auto ;
	vecttmp.dec2_dzpuis() ;

	// 2/ Interpolation of the result
	//## OUTSIDE THE ZEC

	d_logn_comp.set_etat_qcq() ;
	(d_logn_comp.set(0)).import_symy(nzm1, vecttmp(0) ) ;  // d/dx sym.
	(d_logn_comp.set(1)).import_asymy(nzm1, vecttmp(1) ) ; // d/dy antisym.
	(d_logn_comp.set(2)).import_symy(nzm1, vecttmp(2) ) ;  // d/dz sym.

    }
    
    d_logn_comp.set_triad( *((comp.d_logn_auto).get_triad()) ) ;


    // Computation of d_beta_comp
    // --------------------------

    if ( (comp.d_beta_auto).get_etat() == ETATZERO ) {
	d_beta_comp.set_etat_zero() ; 
    }
    else {
	// 1/ Division by r^2 of comp.d_logn_auto in the ZEC
	Tenseur vecttmp = comp.d_beta_auto ;
	vecttmp.dec2_dzpuis() ;

	// 2/ Interpolation of the result 
	//## OUTSIDE THE ZEC

	d_beta_comp.set_etat_qcq() ;

	(d_beta_comp.set(0)).import_symy(nzm1, vecttmp(0) ) ;  // d/dx sym.
	(d_beta_comp.set(1)).import_asymy(nzm1, vecttmp(1) ) ; // d/dy antisym.
	(d_beta_comp.set(2)).import_symy(nzm1, vecttmp(2) ) ;  // d/dz sym.

    }

    d_beta_comp.set_triad( *((comp.d_beta_auto).get_triad()) ) ;

    // Computation of tkij_comp
    // ------------------------
    
    if ( (comp.tkij_auto).get_etat() == ETATZERO ) {
	tkij_comp.set_etat_zero() ;
    }
    else{

	// 1/ Division by r^2 of comp.d_logn_auto in the ZEC
	Tenseur_sym tenstmp = comp.tkij_auto ;
	tenstmp.dec2_dzpuis() ;

	// 2/ Interpolation of the result
	//## OUTSIDE THE ZEC

	tkij_comp.set_etat_qcq() ;

	(tkij_comp.set(0, 0)).import_asymy(nzm1, tenstmp(0, 0) ) ; // K_xx antisym
	(tkij_comp.set(0, 1)).import_symy(nzm1, tenstmp(0, 1) ) ;  // K_xy sym.
	(tkij_comp.set(0, 2)).import_asymy(nzm1, tenstmp(0, 2) ) ; // K_xz antisym
	(tkij_comp.set(1, 1)).import_asymy(nzm1, tenstmp(1, 1) ) ; // K_yy antisym.
	(tkij_comp.set(1, 2)).import_symy(nzm1, tenstmp(1, 2) ) ;  // K_yz sym
	(tkij_comp.set(2, 2)).import_asymy(nzm1, tenstmp(2, 2) ) ; // K_zz antisym.

    }

    tkij_comp.set_triad( *((comp.tkij_auto).get_triad()) ) ;

    if (relativistic) {
	// Computation of akcar_comp
	// -------------------------

	akcar_comp.set_etat_qcq() ;

	akcar_comp.set() = 0 ;
    
	for (int i=0; i<3; i++) {
	    for (int j=0; j<3; j++) {

		akcar_comp.set() += tkij_auto(i, j) * tkij_comp(i, j) ; 

	    }
	}

	akcar_comp = a_car * akcar_comp ;
    }


    // The derived quantities are obsolete
    // -----------------------------------

    del_deriv() ;
    */
    //#################################

}
}
