/*
 * Methods Map_radial::comp_x_from_spherical
 *	   Map_radial::comp_y_from_spherical
 *	   Map_radial::comp_z_from_spherical
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


char map_radial_comp_xyz_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2005/09/15 15:51:25  j_novak
 * The "rotation" (change of triad) methods take now Scalars as default
 * arguments.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/09/11  15:56:22  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <assert.h>

// Headers Lorene
#include "tensor.h"
#include "cmp.h"


		    //------------------------------------//
		    //		X  component		  //
		    //------------------------------------//
void Map_radial::comp_x_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
				       const Cmp& v_phi, Cmp& v_x) const {
    Scalar resu(v_x) ;
    comp_x_from_spherical(Scalar(v_r), Scalar(v_theta), Scalar(v_phi), 
			  resu ) ;
    v_x = resu ;
}

void Map_radial::comp_x_from_spherical(const Scalar& v_r, const Scalar& v_theta, 
				       const Scalar& v_phi, Scalar& v_x) const {
				       

    // Protections
    // -----------
    assert(v_r.get_etat() != ETATNONDEF) ; 
    assert(v_theta.get_etat() != ETATNONDEF) ; 
    assert(v_phi.get_etat() != ETATNONDEF) ; 

    assert(v_r.get_mp() == *this) ; 
    assert(v_theta.get_mp() == *this) ; 
    assert(v_phi.get_mp() == *this) ; 
    
    int dzp ;
    if ( v_r.dz_nonzero() ) {
	dzp = v_r.get_dzpuis() ; 
    }
    else{
	if ( v_theta.dz_nonzero() ) {
	    dzp = v_theta.get_dzpuis() ; 
	}
	else{
	    dzp = v_phi.get_dzpuis() ; 
	}
    }
     
    assert( v_r.check_dzpuis(dzp) ) ; 
    assert( v_theta.check_dzpuis(dzp) ) ; 
    assert( v_phi.check_dzpuis(dzp) ) ; 
    
    // Computation
    // -----------
    const Valeur& w_r = v_r.get_spectral_va() ; 
    const Valeur& w_t = v_theta.get_spectral_va() ; 
    const Valeur& w_p = v_phi.get_spectral_va() ; 
    
    Valeur tmp = w_r.mult_st() + w_t.mult_ct() ;

    v_x = tmp.mult_cp() - w_p.mult_sp() ; 
    
    v_x.set_dzpuis(dzp) ; 
	  
}
		    

		    //------------------------------------//
		    //		Y  component		  //
		    //------------------------------------//

void Map_radial::comp_y_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
				       const Cmp& v_phi, Cmp& v_y) const {

    Scalar resu(v_y) ;
    comp_y_from_spherical(Scalar(v_r), Scalar(v_theta), Scalar(v_phi), 
			  resu ) ;
    v_y = resu ;
}

void Map_radial::comp_y_from_spherical(const Scalar& v_r, const Scalar& v_theta, 
				       const Scalar& v_phi, Scalar& v_y) const {
				       

    // Protections
    // -----------
    assert(v_r.get_etat() != ETATNONDEF) ; 
    assert(v_theta.get_etat() != ETATNONDEF) ; 
    assert(v_phi.get_etat() != ETATNONDEF) ; 

    assert(v_r.get_mp() == *this) ; 
    assert(v_theta.get_mp() == *this) ; 
    assert(v_phi.get_mp() == *this) ; 
    
    int dzp ;
    if ( v_r.dz_nonzero() ) {
	dzp = v_r.get_dzpuis() ; 
    }
    else{
	if ( v_theta.dz_nonzero() ) {
	    dzp = v_theta.get_dzpuis() ; 
	}
	else{
	    dzp = v_phi.get_dzpuis() ; 
	}
    }
     
    assert( v_r.check_dzpuis(dzp) ) ; 
    assert( v_theta.check_dzpuis(dzp) ) ; 
    assert( v_phi.check_dzpuis(dzp) ) ; 
    
    // Computation
    // -----------
    const Valeur& w_r = v_r.get_spectral_va() ; 
    const Valeur& w_t = v_theta.get_spectral_va() ; 
    const Valeur& w_p = v_phi.get_spectral_va() ; 
    
    Valeur tmp = w_r.mult_st() + w_t.mult_ct() ;

    v_y = tmp.mult_sp() + w_p.mult_cp() ; 
    
    v_y.set_dzpuis(dzp) ; 
	  
}
		    
		    //------------------------------------//
		    //		Z  component		  //
		    //------------------------------------//

void Map_radial::comp_z_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
				       Cmp& v_z) const {
    Scalar resu = v_z ;
    comp_z_from_spherical(Scalar(v_r), Scalar(v_theta), resu) ;
    v_z = resu ;

}

void Map_radial::comp_z_from_spherical(const Scalar& v_r, const Scalar& v_theta, 
				       Scalar& v_z) const {
				       

    // Protections
    // -----------
    assert(v_r.get_etat() != ETATNONDEF) ; 
    assert(v_theta.get_etat() != ETATNONDEF) ; 

    assert(v_r.get_mp() == *this) ; 
    assert(v_theta.get_mp() == *this) ; 
    
    int dzp ;
    if ( v_r.dz_nonzero() ) {
	dzp = v_r.get_dzpuis() ; 
    }
    else{
	dzp = v_theta.get_dzpuis() ; 
    }
     
    assert( v_r.check_dzpuis(dzp) ) ; 
    assert( v_theta.check_dzpuis(dzp) ) ; 
    
    // Computation
    // -----------
    const Valeur& w_r = v_r.get_spectral_va() ; 
    const Valeur& w_t = v_theta.get_spectral_va() ; 
    
    v_z = w_r.mult_ct() - w_t.mult_st() ; 
    
    v_z.set_dzpuis(dzp) ; 
	  
}
		    
		    
