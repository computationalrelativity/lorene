/*
 * Method Etoile_bin::update_metric_der_comp
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


char et_bin_upmetr_der_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
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

void Etoile_bin::update_metric_der_comp(const Etoile_bin& comp) {
    
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
    
    
}
