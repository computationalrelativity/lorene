/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


char map_radial_comp_rtp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.1  2000/09/19  15:25:50  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */


// Headers C
#include <assert.h>

// Headers Lorene
#include "map.h"
#include "cmp.h"


		    //------------------------------------//
		    //		r  component		  //
		    //------------------------------------//

void Map_radial::comp_r_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
				       const Cmp& v_z, Cmp& v_r) const {
				       

    // Protections
    // -----------
    assert(v_x.get_etat() != ETATNONDEF) ; 
    assert(v_y.get_etat() != ETATNONDEF) ; 
    assert(v_z.get_etat() != ETATNONDEF) ; 

    assert(v_x.get_mp() == this) ; 
    assert(v_y.get_mp() == this) ; 
    assert(v_z.get_mp() == this) ; 
    
    int dzp ;
    if ( v_x.dz_nonzero() ) {
	dzp = v_x.get_dzpuis() ; 
    }
    else{
	if ( v_y.dz_nonzero() ) {
	    dzp = v_y.get_dzpuis() ; 
	}
	else{
	    dzp = v_z.get_dzpuis() ; 
	}
    }
     
    assert( v_x.check_dzpuis(dzp) ) ; 
    assert( v_y.check_dzpuis(dzp) ) ; 
    assert( v_z.check_dzpuis(dzp) ) ; 
    
    // Computation
    // -----------
    const Valeur& w_x = v_x.va ; 
    const Valeur& w_y = v_y.va ; 
    const Valeur& w_z = v_z.va ; 
    
    Valeur tmp = w_x.mult_cp() + w_y.mult_sp() ;

    v_r = tmp.mult_st() + w_z.mult_ct() ; 
    
    v_r.set_dzpuis(dzp) ; 
	  
}
		    

		    //------------------------------------//
		    //		Theta  component	  //
		    //------------------------------------//

void Map_radial::comp_t_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
				       const Cmp& v_z, Cmp& v_t) const {
				       

    // Protections
    // -----------
    assert(v_x.get_etat() != ETATNONDEF) ; 
    assert(v_y.get_etat() != ETATNONDEF) ; 
    assert(v_z.get_etat() != ETATNONDEF) ; 

    assert(v_x.get_mp() == this) ; 
    assert(v_y.get_mp() == this) ; 
    assert(v_z.get_mp() == this) ; 
    
    int dzp ;
    if ( v_x.dz_nonzero() ) {
	dzp = v_x.get_dzpuis() ; 
    }
    else{
	if ( v_y.dz_nonzero() ) {
	    dzp = v_y.get_dzpuis() ; 
	}
	else{
	    dzp = v_z.get_dzpuis() ; 
	}
    }
     
    assert( v_x.check_dzpuis(dzp) ) ; 
    assert( v_y.check_dzpuis(dzp) ) ; 
    assert( v_z.check_dzpuis(dzp) ) ; 
    
    // Computation
    // -----------
    const Valeur& w_x = v_x.va ; 
    const Valeur& w_y = v_y.va ; 
    const Valeur& w_z = v_z.va ; 
    
    Valeur tmp = w_x.mult_cp() + w_y.mult_sp() ;

    v_t = tmp.mult_ct() - w_z.mult_st() ; 
    
    v_t.set_dzpuis(dzp) ; 
	  
}
		    
		    //------------------------------------//
		    //		Phi  component		  //
		    //------------------------------------//

void Map_radial::comp_p_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
				       Cmp& v_p) const {
				       

    // Protections
    // -----------
    assert(v_x.get_etat() != ETATNONDEF) ; 
    assert(v_y.get_etat() != ETATNONDEF) ; 

    assert(v_x.get_mp() == this) ; 
    assert(v_y.get_mp() == this) ; 
    
    int dzp ;
    if ( v_x.dz_nonzero() ) {
	dzp = v_x.get_dzpuis() ; 
    }
    else{
	dzp = v_y.get_dzpuis() ; 
    }
     
    assert( v_x.check_dzpuis(dzp) ) ; 
    assert( v_y.check_dzpuis(dzp) ) ; 
    
    // Computation
    // -----------
    const Valeur& w_x = v_x.va ; 
    const Valeur& w_y = v_y.va ; 
    
    v_p = - w_x.mult_sp() + w_y.mult_cp() ; 
    
    v_p.set_dzpuis(dzp) ; 
	  
}
		    
		    
