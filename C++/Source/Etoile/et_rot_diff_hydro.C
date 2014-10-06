/*
 * Function Et_rot_diff::hydro_euler
 *
 * (see file et_rot_diff.h for documentation)
 *
 */

/*
 *   Copyright (c) 2001 Eric Gourgoulhon
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


char et_rot_diff_hydro_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:13:09  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2001/10/19  08:18:36  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */


// Headers C
#include <cstdlib>

// Headers Lorene
#include "et_rot_diff.h"
#include "utilitaires.h"

void Et_rot_diff::hydro_euler(){

    int nz = mp.get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 

    // Computation of u_euler
    // ----------------------
    
    Cmp x(mp) ; 
    Cmp y(mp) ; 
    x = mp.x ; 
    y = mp.y ; 
    
    u_euler.set_etat_qcq() ; 
    
    // Cartesian components of differential rotation:

    u_euler.set(0) = - omega_field() * y ;
    u_euler.set(1) =   omega_field() * x ;
    u_euler.set(2) = 0 ;
    u_euler.annule(nzm1) ; 
    
    u_euler.set_triad( mp.get_bvect_cart() ) ;	// Triad = Cartesian triad
    
    u_euler.set_std_base() ;	// sets the standard bases for spectral expansions

    u_euler = ( u_euler - shift ) / nnn ; 

    u_euler.set_std_base() ;	// sets the standard bases for spectral expansions

//## Test
    Tenseur utest(mp, 1, CON, mp.get_bvect_spher()) ; 
    utest.set_etat_qcq() ; 
    
    utest.set(0) = 0 ;	    // Spherical components of differential rotation
    utest.set(1) = 0 ;
    utest.set(2) = ( omega_field() - nphi() ) / nnn();

    utest.set(2).annule(nzm1) ; 
    utest.set(2).std_base_scal() ;
    utest.set(2).mult_rsint() ;	    //  Multiplication by r sin(theta)
    
    utest.set_triad( mp.get_bvect_spher() ) ; 

    utest.change_triad( mp.get_bvect_cart() ) ; 
    
    for (int i=0; i<3; i++) {
	Valeur& uu = u_euler.set(i).va ;
	Valeur& ut = utest.set(i).va ;
	
	if (uu.get_etat() != ETATZERO) {
	    uu.coef() ; 
	    
	    if (ut.get_etat() == ETATZERO) {
		ut.set_etat_cf_qcq() ; 
		*(ut.c_cf) = 0 ; 
		ut.c_cf->base = uu.c_cf->base ; 
	    }
	    else {
		ut.coef() ; 
	    }
	    
	    Mtbl_cf diff = *(uu.c_cf) - *(ut.c_cf) ;
	    cout << "Et_rot_diff::hydro_euler: test u_euler(" << i << ") : " 
		 << max( abs(diff) )(0) << endl ; 
	
	}
    }
//##

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
