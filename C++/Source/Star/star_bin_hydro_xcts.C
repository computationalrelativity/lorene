/*
 * Methods of the class Star_bin_xcts for computing hydro quantities
 * (see file star.h for documentation)
 */
 
/*
 *   Copyright (c) 2010 Michal Bejger
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char star_bin_hydro_xcts_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2010/05/04 07:51:05  m_bejger
 * Initial version
 *
 * $Header$
 *
 */ 

// Headers Lorene
#include "star.h"

void Star_bin_xcts::hydro_euler(){

    int nz = mp.get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 

    Sym_tensor gamma_cov (gamma.cov()) ;
    Sym_tensor gamma_con (gamma.con()) ;
    gamma_cov.change_triad(mp.get_bvect_cart()) ;
    gamma_con.change_triad(mp.get_bvect_cart()) ;

    //----------------------------------
    // Specific relativistic enthalpy		    ---> hhh
    //----------------------------------
    
    Scalar hhh = exp(ent) ;  // = 1 at the Newtonian limit
    hhh.std_spectral_base() ;

    //---------------------------------------------------
    // Lorentz factor between the co-orbiting
    // observer and the Eulerian one
    // Eqs. 23 and 24 from Gourgoulhon et al. (2001)
    //---------------------------------------------------

    Scalar gam0 = 1 / sqrt( 1 - contract(gamma_cov, 0, 1, bsn * bsn, 0, 1)) ;
    gam0.std_spectral_base() ;
    
    //------------------------------------------
    // Lorentz factor and 3-velocity of the fluid 
    // with respect to the Eulerian observer
    //------------------------------------------
    
    if (irrotational) {	

        d_psi.std_spectral_base() ;

	// See Eq. 32 from Gourgoulhon et al. (2001)
	gam_euler = sqrt( 1 + contract(gamma_con, 0, 1, d_psi * d_psi, 0, 1) 
			  / (hhh%hhh) ) ; 

	gam_euler.std_spectral_base() ; 
	
	u_euler = contract(gamma_con, 0, d_psi, 0)/( hhh % gam_euler ) ;
	u_euler.std_spectral_base() ; 

    } else {
		
	// Rigid rotation
	// --------------

	gam_euler = gam0 ; 
	gam_euler.std_spectral_base() ; 
	u_euler = bsn ; 
	
    }
    
    //------------------------------------
    // Energy density E with respect to the Eulerian observer
    // Eq. 53 from Gourgoulhon et al. (2001)  
    //--------------------------------------
    
    ener_euler = gam_euler % gam_euler % ( ener + press ) - press ; 
    
    //-------------------------------------------
    // Trace of the stress tensor with respect to the Eulerian observer
    // See Eq (54) from Gourgoulhon et al. (2001)  
    //-------------------------------------------

    s_euler = 3 * press  +  ( ener_euler + press ) %
	contract(gamma_cov, 0, 1, u_euler * u_euler, 0 ,1) ;

    //-------------------------------------------
    // Spatial part of the stress-energy tensor with respect
    // to the Eulerian observer. 
    //-------------------------------------------

    for(int i=1; i<=3; i++) {
		for(int j=1; j<=3; j++){
		    stress_euler.set(i,j) = (ener_euler + press )*u_euler(i)
			*u_euler(j) + press * gamma_con(i,j) ;
		}
    }
    
    //-------------------------------------------
    // Lorentz factor between the fluid and		---> gam
    // co-orbiting observers
    // See Eq (58) from Gourgoulhon et al. (2001)  
    //--------------------------------------------
    
    if (irrotational) {	

	Scalar tmp = ( 1 - contract(gamma_cov, 0, 1, bsn * u_euler, 0, 1) ) ;
	tmp.std_spectral_base() ;
	Scalar gam = gam0 % gam_euler % tmp ;
	
	//-------------------------------------------
	//	Spatial projection of the fluid 3-velocity
	//  with respect to the co-orbiting observer
	//--------------------------------------------
	
	wit_w = - gam_euler / gam * u_euler + gam0 * bsn ; 
	
	wit_w.std_spectral_base() ;  // set the bases for spectral expansions
	wit_w.annule_domain(nzm1) ;	 // zero in the ZEC
	
	
	//-------------------------------------------
	//	Logarithm of the Lorentz factor between 
	//	the fluid and co-orbiting observers
	//--------------------------------------------
	
	loggam = log( gam ) ;
	
	loggam.std_spectral_base() ;   // set the bases for spectral expansions
	
	//------------------------------------------------
	// Velocity fields set to zero in external domains
	//------------------------------------------------
	
		loggam.annule_domain(nzm1) ;	// zero in the ZEC only
		wit_w.annule_domain(nzm1) ;		// zero outside the star     
		u_euler.annule_domain(nzm1) ;	// zero outside the star     
		loggam.set_dzpuis(0) ; 
    
    } else {
	    
		loggam = 0 ; 
		wit_w.set_etat_zero() ; 
    }
    
    // The derived quantities are obsolete
    // -----------------------------------
    
    del_deriv() ;                
    
}
