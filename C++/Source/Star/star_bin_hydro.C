/*
 * Methods of the class Star_bin for computing hydro quantities
 *
 * (see file star.h for documentation)
 */

/*
 *   Copyright (c) 2004 Francois Limousin
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


char star_bin_hydro_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2004/04/08 16:34:39  f_limousin
 * Changes for irrotationnal binaries.
 *
 * Revision 1.3  2004/02/27 09:54:48  f_limousin
 * Many minor changes.
 *
 * Revision 1.2  2004/01/20 15:18:31  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// Headers C

// Headers Lorene
#include "star.h"

void Star_bin::hydro_euler(){

    int nz = mp.get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 

    //----------------------------------
    // Specific relativistic enthalpy		    ---> hhh
    //----------------------------------
    
    Scalar hhh = exp(ent) ;  // = 1 at the Newtonian limit
    hhh.std_spectral_base() ;

    //---------------------------------------------------
    // Lorentz factor between the co-orbiting	          ---> gam00
    // observer and the Eulerian one
    // See Eq (23) and (24) from Gourgoulhon et al. (2001)
    //---------------------------------------------------

    Scalar gam0 = 1 / sqrt( 1 - sprod(bsn,bsn) ) ;
    gam0.std_spectral_base() ;

    //------------------------------------------
    // Lorentz factor and 3-velocity of the fluid 
    //  with respect to the Eulerian observer
    //------------------------------------------
    
    if (irrotational) {	

        d_psi.std_spectral_base() ;

	// See Eq (32) from Gourgoulhon et al. (2001)
	gam_euler = sqrt( 1 + sprod(d_psi, d_psi) / (hhh%hhh) ) ; 

	gam_euler.std_spectral_base() ; 
	
	d_psi.change_triad(mp.get_bvect_spher()) ;
	u_euler = contract(gamma.con(), 0, d_psi, 0)/( hhh % gam_euler ) ;
	u_euler.change_triad(mp.get_bvect_cart()) ;
	d_psi.change_triad(mp.get_bvect_cart()) ;

	u_euler.std_spectral_base() ; 

    }
    else {
	// Rigid rotation
	// --------------

	gam_euler = gam0 ; 
	gam_euler.std_spectral_base() ; 
	u_euler = - bsn ; 
	
    }
    
    //------------------------------------
    //  Energy density E with respect to the Eulerian observer
    // See Eq (53) from Gourgoulhon et al. (2001)  
    //------------------------------------
    
    ener_euler = gam_euler % gam_euler % ( ener + press ) - press ; 
    
    //------------------------------------
    // Trace of the stress tensor with respect to the Eulerian observer
    // See Eq (54) from Gourgoulhon et al. (2001)  
    //------------------------------------

    s_euler = 3 * press  +  ( ener_euler + press ) %
	sprod(u_euler, u_euler) ;
    s_euler.std_spectral_base() ; 


    //-------------------------------------------
    // Spatial part of the stress-energy tensor with respect
    // to the Eulerian observer. 
    //-------------------------------------------

    stress_euler.change_triad(mp.get_bvect_cart()) ;
    
    for(int i=1; i<=3; i++){
	for(int j=1; j<=3; j++){
	    stress_euler.set(i,j) = (ener_euler + press )*u_euler(i)
		*u_euler(j) + press*gtilde.con()(i,j) ;
	}
    }
    stress_euler.std_spectral_base() ;
    stress_euler.change_triad(mp.get_bvect_spher()) ;

    
    //-------------------------------------------
    //	Lorentz factor between the fluid and		---> gam
    //	co-orbiting observers
    // See Eq (58) from Gourgoulhon et al. (2001)  
    //--------------------------------------------
    
    if (irrotational) {	

	Scalar tmp = ( 1 + sprod(bsn,u_euler) ) ;
	tmp.std_spectral_base() ;
	Scalar gam = gam0 % gam_euler % tmp ;
	
	//-------------------------------------------
	//	Spatial projection of the fluid 3-velocity
	//  with respect to the co-orbiting observer
	//--------------------------------------------
	
	wit_w = gam_euler / gam * u_euler + gam0 * bsn ; 
	
	wit_w.std_spectral_base() ;  // set the bases for spectral expansions
	
	wit_w.annule_domain(nzm1) ;	// zero in the ZEC
	
	
	//-------------------------------------------
	//	Logarithm of the Lorentz factor between 
	//	the fluid and co-orbiting observers
	//--------------------------------------------
	
	loggam = log( gam ) ;
	
	loggam.std_spectral_base() ;   // set the bases for spectral expansions
	
	//-------------------------------------------------
	// Velocity fields set to zero in external domains
	//-------------------------------------------------
	
	loggam.annule_domain(nzm1) ;	    // zero in the ZEC only
	
	wit_w.annule_domain(nzm1) ;		// zero outside the star     
	
	u_euler.annule_domain(nzm1) ;	// zero outside the star     
	
	loggam.set_dzpuis(0) ; 
    }
    else {
	    
	loggam = 0 ; 
	wit_w.set_etat_zero() ; 
    }
    
    // The derived quantities are obsolete
    // -----------------------------------
    
    del_deriv() ;                
    
/*    cout << " comparaison " << endl ;
    cout << "a_car" << endl << norme(psi4) << endl ;
    cout << "logn_auto" << endl << norme(logn_auto) << endl ;
    
    shift_auto.change_triad(mp.get_bvect_cart()) ;
    cout << "shift_x" << endl << norme(shift_auto(1)) << endl ;
    cout << "shift_y" << endl << norme(shift_auto(2)) << endl ;
    cout << "shift_z" << endl << norme(shift_auto(3)) << endl ;
    shift_auto.change_triad(mp.get_bvect_spher()) ;
    
    cout << "qq" << endl << norme(qq) << endl ;
    cout << "hij" << endl << norme(hij(1,1)) << endl ;
  
    u_euler.change_triad(mp.get_bvect_cart()) ;
    cout << "u_euler_x" << endl << norme(u_euler(1)) << endl ;
    cout << "u_euler_y" << endl << norme(u_euler(2)) << endl ;
    cout << "u_euler_z" << endl << norme(u_euler(3)) << endl ;
    u_euler.change_triad(mp.get_bvect_spher()) ;
*/

}
