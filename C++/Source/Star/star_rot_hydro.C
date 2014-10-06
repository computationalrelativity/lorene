/*
 * Method Star_rot::hydro_euler
 *
 * (see file star_rot.h for documentation)
 *
 */

/*
 *   Copyright (c) 2010 Eric Gourgoulhon
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


char star_rot_hydro_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:13:17  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2010/01/25 18:15:52  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>

// Headers Lorene
#include "star_rot.h"
#include "utilitaires.h"

void Star_rot::hydro_euler(){

    int nz = mp.get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 

    // Computation of u_euler
    // ----------------------
    
    const Coord& x = mp.x ; 
    const Coord& y = mp.y ; 
    
    u_euler.set_etat_qcq() ; 
    
    u_euler.set(1) = - omega * y ;	    // Cartesian components of solid rotation
    u_euler.set(2) =   omega * x ;
    u_euler.set(3) = 0 ;
    u_euler.annule_domain(nzm1) ; 
    
    u_euler.set_triad( mp.get_bvect_cart() ) ;	// Triad = Cartesian triad
    
    u_euler = ( u_euler + beta ) / nn ; 

    u_euler.std_spectral_base() ;	// sets the standard bases for spectral expansions

    if ( (u_euler(1).get_etat() == ETATZERO) &&
	 (u_euler(2).get_etat() == ETATZERO) &&
	 (u_euler(3).get_etat() == ETATZERO) )    {
	
	u_euler.set_etat_zero() ;    
    }


    // Computation of uuu (norme of u_euler)
    // ------------------

    // The scalar product is performed on the spherical components: 

    Vector us = u_euler ; 
    us.change_triad( mp.get_bvect_spher() ) ; 

    Scalar uuu2 = a_car * ( us(1)*us(1) + us(2)*us(2) ) + b_car * us(3)*us(3) ; 

    uuu = sqrt( uuu2 ) ; 
    
    if (uuu.get_etat() == ETATQCQ) {
	// Same basis as (Omega -N^phi) r sin(theta)
	(uuu.set_spectral_va()).set_base( us(3).get_spectral_va().get_base() ) ;   
    }						  

    // Lorentz factor
    // --------------
    
    Scalar u2 = unsurc2 * uuu2 ; 
    
    Scalar gam2 = double(1) / (double(1) - u2) ; 
    
    gam_euler = sqrt(gam2) ; 

    gam_euler.std_spectral_base() ; // sets the standard spectral bases for
				    // a scalar field

    //  Energy density E with respect to the Eulerian observer
    //--------------------------------------------------------
    
    ener_euler = gam2 * ( ener + press ) - press ; 

    // Trace of the stress tensor with respect to the Eulerian observer
    //------------------------------------
    
    s_euler = 3 * press  +  ( ener_euler + press ) * u2  ;
    
    // The derived quantities are obsolete
    // -----------------------------------
    
    del_deriv() ;                
    
}
