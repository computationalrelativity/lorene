/*
 * Methods Etoile_bin::update_metric
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
 * Revision 1.6  2016/12/05 16:17:53  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:52:56  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2003/10/24 12:26:38  k_taniguchi
 * Suppress the method of update metric for NS-BH
 *
 * Revision 1.3  2003/10/24 11:46:07  k_taniguchi
 * Change some notations
 *
 * Revision 1.2  2002/12/19 14:52:42  e_gourgoulhon
 * Added the new function
 *         void update_metric(const Bhole& comp)
 * to treat the case where the companion is a black hole
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.9  2000/09/27  12:49:57  keisuke
 * Utilisation de d_logn_auto_div dans le calcul de d_logn_auto dans
 * la version avec relaxation.
 *
 * Revision 2.8  2000/09/22  15:53:06  keisuke
 * Calcul de d_logn_auto prenant en compte d_logn_auto_div.
 *
 * Revision 2.7  2000/03/07  14:55:08  eric
 * Ajout de l'appel a extrinsic_curvature.
 *
 * Revision 2.6  2000/03/07  08:33:24  eric
 * Appel de Cmp::import_sym / asym (pour tenir compte de la symetrie /
 *  plan y=0).
 *
 * Revision 2.5  2000/02/12  18:38:11  eric
 * Ajout de la version avec relaxation.
 * Appel de set_std_base() sur nnn et a_car.
 *
 * Revision 2.4  2000/02/12  11:42:49  eric
 * Appel de Tenseur::set_std_base() sur les Tenseurs importes du
 * compagnon.
 *
 * Revision 2.3  2000/02/10  18:54:41  eric
 * Traitement du cas ETATZERO.
 *
 * Revision 2.2  2000/02/10  16:55:10  eric
 * Appel de change_triad sur d_logn_auto et d_beta_auto.
 *
 * Revision 2.1  2000/02/04  17:14:32  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/02/04  16:38:00  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "etoile.h"
#include "bhole.h"

		    //----------------------------------//
		    //	 Version without relaxation	//
		    //----------------------------------//

namespace Lorene {
void Etoile_bin::update_metric(const Etoile_bin& comp) {

    // Computation of quantities coming from the companion
    // ---------------------------------------------------

    if ( (comp.logn_auto).get_etat() == ETATZERO ) {
	logn_comp.set_etat_zero() ;
    }
    else{
	logn_comp.set_etat_qcq() ;
	(logn_comp.set()).import_symy( comp.logn_auto() ) ;
	logn_comp.set_std_base() ;   // set the bases for spectral expansions
    }


    if ( (comp.beta_auto).get_etat() == ETATZERO ) {
	beta_comp.set_etat_zero() ;
    }
    else{
	beta_comp.set_etat_qcq() ;
	(beta_comp.set()).import_symy( comp.beta_auto() ) ; 
	beta_comp.set_std_base() ;   // set the bases for spectral expansions
    }


    if ( (comp.shift_auto).get_etat() == ETATZERO ) {
	shift_comp.set_etat_zero() ;
    }
    else{
	shift_comp.set_etat_qcq() ; 

	(shift_comp.set(0)).import_asymy( comp.shift_auto(0) ) ;  // N^x antisym
	(shift_comp.set(1)).import_symy( comp.shift_auto(1) ) ;   // N^y sym.
	(shift_comp.set(2)).import_asymy( comp.shift_auto(2) ) ;  // N^z anisym

	shift_comp.set_std_base() ;   // set the bases for spectral expansions
    }
    shift_comp.set_triad( *((comp.shift_auto).get_triad()) ) ;
    


    // Lapse function N
    // ----------------

    Tenseur logn_total = logn_auto + logn_comp ; 

    nnn = exp( unsurc2 * logn_total ) ;

    nnn.set_std_base() ;   // set the bases for spectral expansions

    // Conformal factor A^2
    // ---------------------
    
    a_car = exp( 2*unsurc2*( beta_auto + beta_comp - logn_total ) ) ;

    a_car.set_std_base() ;   // set the bases for spectral expansions

    // Shift vector N^i
    // ----------------

    shift = shift_auto + shift_comp ;

    // Derivatives of metric coefficients
    // ----------------------------------

    // ... (d/dX,d/dY,d/dZ)(logn_auto) :
    d_logn_auto_regu = logn_auto_regu.gradient() ;    // (d/dx, d/dy, d/dz)
    d_logn_auto_regu.change_triad(ref_triad) ;   // -->  (d/dX, d/dY, d/dZ)

    if ( *(d_logn_auto_div.get_triad()) != ref_triad ) {

	// Change the basis from spherical coordinate to Cartesian one
	d_logn_auto_div.change_triad( mp.get_bvect_cart() ) ;

	// Change the basis from mapping coordinate to absolute one
	d_logn_auto_div.change_triad( ref_triad ) ;

    }

    d_logn_auto = d_logn_auto_regu + d_logn_auto_div ;

    // ... (d/dX,d/dY,d/dZ)(beta_auto) :
    d_beta_auto = beta_auto.gradient() ;    // (d/dx, d/dy, d/dz)
    d_beta_auto.change_triad(ref_triad) ;   // -->  (d/dX, d/dY, d/dZ)

    if (relativistic) {
	// ... extrinsic curvature (tkij_auto and akcar_auto)
	extrinsic_curvature() ;
    }

    // The derived quantities are obsolete
    // -----------------------------------

    del_deriv() ;


}



		    //----------------------------------//
		    //	  Version with relaxation       //
		    //----------------------------------//

void Etoile_bin::update_metric(const Etoile_bin& comp,
			       const Etoile_bin& star_jm1, double relax) {


    // Computation of quantities coming from the companion
    // ---------------------------------------------------

    if ( (comp.logn_auto).get_etat() == ETATZERO ) {
	logn_comp.set_etat_zero() ;
    }
    else{
	logn_comp.set_etat_qcq() ;
	(logn_comp.set()).import_symy( comp.logn_auto() ) ;
	logn_comp.set_std_base() ;   // set the bases for spectral expansions
    }


    if ( (comp.beta_auto).get_etat() == ETATZERO ) {
	beta_comp.set_etat_zero() ;
    }
    else{
	beta_comp.set_etat_qcq() ;
	(beta_comp.set()).import_symy( comp.beta_auto() ) ;
	beta_comp.set_std_base() ;   // set the bases for spectral expansions
    }

    
    if ( (comp.shift_auto).get_etat() == ETATZERO ) {
	shift_comp.set_etat_zero() ; 
    }
    else{  
	shift_comp.set_etat_qcq() ; 

	(shift_comp.set(0)).import_asymy( comp.shift_auto(0) ) ;  // N^x antisym
	(shift_comp.set(1)).import_symy( comp.shift_auto(1) ) ;   // N^y sym.
	(shift_comp.set(2)).import_asymy( comp.shift_auto(2) ) ;  // N^z anisym

	shift_comp.set_std_base() ;   // set the bases for spectral expansions
    }
    shift_comp.set_triad( *((comp.shift_auto).get_triad()) ) ;  
    
    // Relaxation on logn_comp, beta_comp, shift_comp
    // ----------------------------------------------
    double relaxjm1 = 1. - relax ; 
    
    logn_comp = relax * logn_comp + relaxjm1 * (star_jm1.get_logn_comp()) ; 
    
    beta_comp = relax * beta_comp + relaxjm1 * (star_jm1.get_beta_comp()) ; 
    
    shift_comp = relax * shift_comp + relaxjm1 * (star_jm1.get_shift_comp()) ; 
        
    // Lapse function N
    // ----------------
    
    Tenseur logn_total = logn_auto + logn_comp ; 
    
    nnn = exp( unsurc2 * logn_total ) ; 
    
    nnn.set_std_base() ;   // set the bases for spectral expansions
    
    // Conformal factor A^2
    // ---------------------
    
    a_car = exp( 2*unsurc2*( beta_auto + beta_comp - logn_total ) ) ; 

    a_car.set_std_base() ;   // set the bases for spectral expansions
    
    // Shift vector N^i
    // ----------------
    
    shift = shift_auto + shift_comp ; 

    // Derivatives of metric coefficients
    // ----------------------------------
    
    // ... (d/dX,d/dY,d/dZ)(logn_auto) : 
    d_logn_auto_regu = logn_auto_regu.gradient() ;    // (d/dx, d/dy, d/dz)
    d_logn_auto_regu.change_triad(ref_triad) ;   // -->  (d/dX, d/dY, d/dZ)  
    
    if ( *(d_logn_auto_div.get_triad()) != ref_triad ) {

	// Change the basis from spherical coordinate to Cartesian one
	d_logn_auto_div.change_triad( mp.get_bvect_cart() ) ;

	// Change the basis from mapping coordinate to absolute one
	d_logn_auto_div.change_triad( ref_triad ) ;

    }

    d_logn_auto = d_logn_auto_regu + d_logn_auto_div ;      
    
    // ... (d/dX,d/dY,d/dZ)(beta_auto) : 
    d_beta_auto = beta_auto.gradient() ;    // (d/dx, d/dy, d/dz)
    d_beta_auto.change_triad(ref_triad) ;   // -->  (d/dX, d/dY, d/dZ)        

    // ... extrinsic curvature (tkij_auto and akcar_auto)
    extrinsic_curvature() ; 
    
    // The derived quantities are obsolete
    // -----------------------------------
    
    del_deriv() ;                
    
			      
} 
}
