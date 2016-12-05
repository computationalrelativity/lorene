/*
 * Function Etoile_rot::hydro_euler
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
 * Revision 1.5  2016/12/05 16:17:54  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:52:58  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:09  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2006/09/14 07:37:47  j_novak
 * Removal of a test on u_euler.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2000/10/12  15:36:06  eric
 * Le test sur u_euler est desormais effectue sur la totalite de u_euler
 *  (ie comprenant le shift).
 *
 * Revision 2.5  2000/10/06  15:08:13  eric
 * Calcul 3-D de u_euler.
 * uuu s'en deduit.
 *
 * Revision 2.4  2000/08/31  15:38:34  eric
 * Appel du nouvel operateur Cmp::mult_rsint pour le calcul de uuu.
 *
 * Revision 2.3  2000/08/25  12:28:21  eric
 * *** empty log message ***
 *
 * Revision 2.2  2000/08/18  14:01:49  eric
 * *** empty log message ***
 *
 * Revision 2.1  2000/08/17  12:39:56  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/07/21  16:31:00  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>

// Headers Lorene
#include "etoile.h"
#include "utilitaires.h"

namespace Lorene {
void Etoile_rot::hydro_euler(){

    int nz = mp.get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 

    // Computation of u_euler
    // ----------------------
    
    const Coord& x = mp.x ; 
    const Coord& y = mp.y ; 
    
    u_euler.set_etat_qcq() ; 
    
    u_euler.set(0) = - omega * y ;	    // Cartesian components of solid rotation
    u_euler.set(1) =   omega * x ;
    u_euler.set(2) = 0 ;
    u_euler.annule(nzm1) ; 
    
    u_euler.set_triad( mp.get_bvect_cart() ) ;	// Triad = Cartesian triad
    
    u_euler.set_std_base() ;	// sets the standard bases for spectral expansions

    u_euler = ( u_euler - shift ) / nnn ; 

    u_euler.set_std_base() ;	// sets the standard bases for spectral expansions

    if ( (u_euler(0).get_etat() == ETATZERO) &&
	 (u_euler(1).get_etat() == ETATZERO) &&
	 (u_euler(2).get_etat() == ETATZERO) )    {
	
	u_euler = 0 ;    
    }


    // Computation of uuu (norme of u_euler)
    // ------------------

    // The scalar product is performed on the spherical components: 

    Tenseur us = u_euler ; 
    us.change_triad( mp.get_bvect_spher() ) ; 

    Cmp uuu2 =	a_car() * ( us(0) * us(0) + us(1) * us(1) ) 
	     +	b_car() * us(2) * us(2) ; 

    uuu = sqrt( uuu2 ) ; 
    
    if (uuu.get_etat() == ETATQCQ) {
	((uuu.set()).va).set_base( us(2).va.base ) ;   // Same basis as 
    }						   // (Omega -N^phi) r sin(theta)


    // Lorentz factor
    // --------------
    
    Tenseur u2(mp) ; 
    u2 = unsurc2 * uuu2 ; 
    
    Tenseur gam2 = 1 / (1 - u2) ; 
    
    gam_euler = sqrt(gam2) ; 

    gam_euler.set_std_base() ;  // sets the standard spectral bases for
				    // a scalar field

    //  Energy density E with respect to the Eulerian observer
    //------------------------------------
    
    ener_euler = gam2 * ( ener + press ) - press ; 

    ener_euler.set_std_base() ; 
    
    // Trace of the stress tensor with respect to the Eulerian observer
    //------------------------------------
    
    s_euler = 3 * press  +  ( ener_euler + press ) * u2  ;

    s_euler.set_std_base() ; 
    
    // The derived quantities are obsolete
    // -----------------------------------
    
    del_deriv() ;                
    

}
}
