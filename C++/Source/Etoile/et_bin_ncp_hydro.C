/*
 * Methods of the class Et_bin_ncp for computing hydro quantities
 *
 * (see file etoile.h for documentation)
 */

/*
 *   Copyright (c) 2003 Francois Limousin
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


char et_bin_ncp_hydro_C[] = "$Header$" ;

/*
 * $Header$
 *
 */

// Headers C

// Headers Lorene
#include "et_bin_ncp.h"

void Et_bin_ncp::hydro_euler(){

    int nz = mp.get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 

    //----------------------------------
    // Specific relativistic enthalpy		    ---> hhh
    //----------------------------------
    
    Tenseur hhh = exp(unsurc2 * ent) ;  // = 1 at the Newtonian limit
    hhh.set_std_base() ;

    //---------------------------------------------------
    // Lorentz factor between the co-orbiting	          ---> gam0
    // observer and the Eulerian one
    // See Eq (23) and (24) from Gourgoulhon et al. (2001)
    //---------------------------------------------------

    Tenseur gam0 = 1/sqrt(1-unsurc2*sprod(bsn,bsn)) ;
    gam0.set_std_base() ;

    //------------------------------------------
    // Lorentz factor and 3-velocity of the fluid 
    //  with respect to the Eulerian observer
    //------------------------------------------
    
    if (irrotational) {	

//##	cout << "Etoile_bin::hydro_euler : ### warning : " << endl ; 
//	cout << "   d_psi.d_psi would be better computed in spher. coord. !"
//##	     << endl ; 

        d_psi.set_std_base() ;

	// See Eq (32) from Gourgoulhon et al. (2001)
	gam_euler = sqrt( 1 + unsurc2 * sprod(d_psi, d_psi)
				    / (hhh%hhh) ) ; 

	gam_euler.set_std_base() ;  // sets the standard spectral bases for
				    // a scalar field


	// See Eq (31) from Gourgoulhon et al. (2001)
	Tenseur vtmp = d_psi / ( hhh % gam_euler % a_car ) ; 

	// The assignment of u_euler is performed component by component 
	//  because u_euler is contravariant and d_psi is covariant
	if (vtmp.get_etat() == ETATZERO) {
	    u_euler.set_etat_zero() ; 
	}
	else{
	    assert(vtmp.get_etat() == ETATQCQ) ; 
	    u_euler.set_etat_qcq() ; 
	    for (int i=0; i<3; i++) {
		u_euler.set(i) = vtmp(i) ;
	    }
	    u_euler.set_triad( *(vtmp.get_triad()) ) ; 
	}
	
	u_euler.set_std_base() ; 

    }
    else {		// Rigid rotation
			// --------------

	gam_euler = gam0 ; 
	gam_euler.set_std_base() ;  // sets the standard spectral bases for
				    // a scalar field

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

    //Tenseur tmp00 = sprod(u_euler, u_euler) ; 
    //cout << hex << u_euler(0).va.base.b[0] << endl ; 
    //cout << hex << u_euler(1).va.base.b[0] << endl ; 
    //cout << hex << u_euler(2).va.base.b[0] << endl ; 
    //cout << hex << tmp00().va.base.b[0] << endl ; 

    s_euler = 3 * press  +  ( ener_euler + press ) *
      sprod(u_euler, u_euler) ;
    s_euler.set_std_base() ; 


    //-------------------------------------------
    // Spatial part of the stress-energy tensor with respect
    // to the Eulerian observer. 
    //-------------------------------------------

    stress.set_etat_qcq() ;
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
	stress.set(i,j) = (ener_euler() + press() )*u_euler(i)
	  *u_euler(j) + press()*gtilde.set_cov(i,j) ;
      }
    }
    stress.set_std_base() ;
    
    //-------------------------------------------
    //	Lorentz factor between the fluid and		---> gam
    //	co-orbiting observers
    // See Eq (58) from Gourgoulhon et al. (2001)  
    //--------------------------------------------
    
    Tenseur tmp = ( 1+unsurc2*sprod(bsn,u_euler) ) ;
    tmp.set_std_base() ;
    Tenseur gam = gam0 % gam_euler % tmp ;

    //-------------------------------------------
    //	Spatial projection of the fluid 3-velocity
    //  with respect to the co-orbiting observer
    //--------------------------------------------

    wit_w = gam_euler / gam % u_euler + gam0 % bsn ; 

    wit_w.set_std_base() ;  // set the bases for spectral expansions

    wit_w.annule(nzm1) ;	// zero in the ZEC


    //-------------------------------------------
    //	Logarithm of the Lorentz factor between 
    //	the fluid and co-orbiting observers
    //--------------------------------------------

    if (relativistic) {

	loggam = log( gam ) ;
    
        loggam.set_std_base() ;   // set the bases for spectral expansions
    }
    else {
  
	loggam = 0.5 * flat_scalar_prod_desal(wit_w, wit_w) ;

        loggam.set_std_base() ;   // set the bases for spectral expansions

//### Forcage a zero des termes en sin(m*phi) : 
//	loggam.coef() ; 
//	int np = mgrille->np[0] ; 
//	int nt = mgrille->nt[0] ;
//	int nr = mgrille->nr[0] ;
//	int ntnr = nt * nr ; 
//	for (int k = 1; k < np+2; k+=2) {
//	    for (int j = 0; j < nt; j++) {
//		for (int i = 0; i<nr; i++) {
//		    loggam.c_cf[0]->t[0]->t[ntnr*k + nr*j + i] = 0 ;
//		}
//	    }
//	}
//	loggam.c_ajx[0] = 0 ; 
//	loggam.c_aj = 0 ; 
//	loggam.coef_i() ; 
//###
    }


    //-------------------------------------------------
    // Velocity fields set to zero in external domains
    //-------------------------------------------------

    loggam.annule(nzm1) ;	    // zero in the ZEC only

    wit_w.annule(nzm1) ;		// zero outside the star     

    u_euler.annule(nzm1) ;	// zero outside the star     


    if (loggam.get_etat() != ETATZERO) {
	(loggam.set()).set_dzpuis(0) ; 
    }
    
    //################
    // verification: test on gam_euler

    // if (irrotational) {
    
	// Tenseur gam_test = 1. / sqrt( 1 - unsurc2 * sprod(u_euler, u_euler) ) ;
	
	// cout << "Etoile_bin::hydro_euler: test of gam_euler : " << endl ; 
	// cout << "  maximum error : " << endl ;
	// cout << max(gam_test() - gam_euler()) << endl ; 
	//cout << "  relative error : " << endl ; 
	// cout << diffrel(gam_test(), gam_euler()) << endl ; 
	
    // }
        
    //##################


    //### Test

    if (!irrotational) {

    //	Tenseur diff = gam - 1 ; 
    //	cout << "Etoile_bin::hydro_euler: rigid motion: difference gam <-> 1 : " 
    //	     << endl ; 
    //	cout << norme(diff()) / norme(gam()) << endl ; 
    //
    //	cout << "Etoile_bin::hydro_euler: rigid motion: difference wit_w <-> 0 : " 
    //	     << endl ; 
    //	cout << "Component x : " << endl << norme(wit_w(0)) << endl ; 
    //	cout << "Component y : " << endl << norme(wit_w(1)) << endl ; 
    //	cout << "Component z : " << endl << norme(wit_w(2)) << endl ; 

//####
	gam = 1 ; 
	loggam = 0 ; 
	wit_w = 0 ; 
    }
    
    // The derived quantities are obsolete
    // -----------------------------------
    
    del_deriv() ;                
    
    
}
