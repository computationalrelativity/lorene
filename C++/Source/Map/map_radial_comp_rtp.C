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
 * Revision 1.5  2014/10/13 08:53:06  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:13  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2005/09/15 15:51:25  j_novak
 * The "rotation" (change of triad) methods take now Scalars as default
 * arguments.
 *
 * Revision 1.2  2003/06/20 14:46:17  f_limousin
 * Les assert sur le mapping sont realise a partir du mapping meme et non a partir du pointeur sur ce mapping
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2000/09/19  15:25:50  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */


// Headers C
#include <cassert>

// Headers Lorene
#include "tensor.h"
#include "cmp.h"


		    //------------------------------------//
		    //		r  component		  //
		    //------------------------------------//
namespace Lorene {
void Map_radial::comp_r_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
				       const Cmp& v_z, Cmp& v_r) const {
    Scalar resu = v_r ;
    comp_r_from_cartesian(Scalar(v_x), Scalar(v_y), Scalar(v_z), resu) ;
    v_r = resu ;				       
}

void Map_radial::comp_r_from_cartesian(const Scalar& v_x, const Scalar& v_y, 
				       const Scalar& v_z, Scalar& v_r) const {
				       

    // Protections
    // -----------
    assert(v_x.get_etat() != ETATNONDEF) ; 
    assert(v_y.get_etat() != ETATNONDEF) ; 
    assert(v_z.get_etat() != ETATNONDEF) ; 

    assert(v_x.get_mp() == *this) ; 
    assert(v_y.get_mp() == *this) ; 
    assert(v_z.get_mp() == *this) ; 
    
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
    const Valeur& w_x = v_x.get_spectral_va() ; 
    const Valeur& w_y = v_y.get_spectral_va() ; 
    const Valeur& w_z = v_z.get_spectral_va() ; 
    
    Valeur tmp = w_x.mult_cp() + w_y.mult_sp() ;

    v_r = tmp.mult_st() + w_z.mult_ct() ; 
    
    v_r.set_dzpuis(dzp) ; 
	  
}
		    

		    //------------------------------------//
		    //		Theta  component	  //
		    //------------------------------------//
void Map_radial::comp_t_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
				       const Cmp& v_z, Cmp& v_t) const {
    Scalar resu = v_t ;
    comp_t_from_cartesian( Scalar(v_x), Scalar(v_y), Scalar(v_z), resu ) ;
    v_t = resu ;
}

void Map_radial::comp_t_from_cartesian(const Scalar& v_x, const Scalar& v_y, 
				       const Scalar& v_z, Scalar& v_t) const {
				       

    // Protections
    // -----------
    assert(v_x.get_etat() != ETATNONDEF) ; 
    assert(v_y.get_etat() != ETATNONDEF) ; 
    assert(v_z.get_etat() != ETATNONDEF) ; 

    assert(v_x.get_mp() == *this) ; 
    assert(v_y.get_mp() == *this) ; 
    assert(v_z.get_mp() == *this) ; 
    
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
    const Valeur& w_x = v_x.get_spectral_va() ; 
    const Valeur& w_y = v_y.get_spectral_va() ; 
    const Valeur& w_z = v_z.get_spectral_va() ; 
    
    Valeur tmp = w_x.mult_cp() + w_y.mult_sp() ;

    v_t = tmp.mult_ct() - w_z.mult_st() ; 
    
    v_t.set_dzpuis(dzp) ; 
	  
}
		    
		    //------------------------------------//
		    //		Phi  component		  //
		    //------------------------------------//
void Map_radial::comp_p_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
				       Cmp& v_p) const {
    Scalar resu = v_p ;
    comp_p_from_cartesian(Scalar(v_x), Scalar(v_y), resu) ;
    v_p = resu ;
}

void Map_radial::comp_p_from_cartesian(const Scalar& v_x, const Scalar& v_y, 
				       Scalar& v_p) const {
				       

    // Protections
    // -----------
    assert(v_x.get_etat() != ETATNONDEF) ; 
    assert(v_y.get_etat() != ETATNONDEF) ; 

    assert(v_x.get_mp() == *this) ; 
    assert(v_y.get_mp() == *this) ; 
    
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
    const Valeur& w_x = v_x.get_spectral_va() ; 
    const Valeur& w_y = v_y.get_spectral_va() ; 
    
    v_p = - w_x.mult_sp() + w_y.mult_cp() ; 
    
    v_p.set_dzpuis(dzp) ; 
	  
}
		    
		    
}
