/*
 *  Methods for changing the triad of a Vector
 *
 */

/*
 *   Copyright (c) 2003  Eric Gourgoulhon & Jerome Novak
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

char vector_change_triad_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/10/03 14:41:31  e_gourgoulhon
 * Changed some assert.
 *
 * Revision 1.1  2003/09/29 12:52:57  j_novak
 * Methods for changing the triad are implemented.
 *
 *
 * $Header$
 *
 */

// C headers
#include <assert.h>

// Lorene headers
#include "cmp.h"
#include "tensor.h"

void Vector::change_triad(const Base_vect& new_triad) {

  assert(triad != 0x0) ; 
  
  const Base_vect_cart* nbvc = dynamic_cast<const Base_vect_cart*>(&new_triad) ; 
  const Base_vect_spher* nbvs = dynamic_cast<const Base_vect_spher*>(&new_triad) ; 

  assert((nbvc != 0x0) || (nbvs != 0x0)) ;
   
  const Base_vect_cart* bvc = dynamic_cast<const Base_vect_cart*>(triad) ; 
  const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ; 
    
  assert((bvc != 0x0) || (bvs != 0x0)) ;

  // ---------------------------------------------
  // Case where the input triad is a Cartesian one
  // ---------------------------------------------
  if (nbvc != 0x0) {
    assert(nbvs == 0x0) ;
    
    // -----------------------------
    // Case cartesian -> cartesian
    // -----------------------------
    if (bvc != 0x0) {	// The old triad is a cartesian one
      assert(bvs == 0x0) ; 
      
      int ind = nbvc->get_align() * (bvc->get_align()) ; 
      
      switch (ind) {
	
      case 1 : {	// the two bases are aligned : nothing to do
			// -----------------------------------------
	break ; 		
      }
	
      case - 1 : {    // the two bases are anti-aligned 
	set(0) = - set(0) ;	 // V^x --> - V^x
	set(1) = - set(1) ; 	 // V^y --> - V^y
	                         // V^z unchanged
	break ; 
      }

      case 0 : {	// the two basis have not a special relative orientation
			// -----------------------------------------------------
	cout << 
	  "Vector::change_basis : general value of rot_phi "
	     << " not contemplated yet, sorry !" << endl ;
	abort() ; 
	break ; 		
      }
	
      default : {	    // error
	cout << 
	  "Vector::change_basis : unexpected value of ind !" << endl ;
	cout << "  ind = " << ind << endl ; 
	abort() ; 
	break ; 		
      }
      }
      
    }	// end of the cart -> cart basis case
    

    // -----------------------------
    // Case spherical -> cartesian
    // -----------------------------
    if (bvs != 0x0) {	// The old triad is a spherical one

      assert(bvc == 0x0) ; 
      
      // The triads should be the same as that associated 
      // with the mapping :
      assert( *nbvc == mp->get_bvect_cart() ) ; 
      assert( *bvs == mp->get_bvect_spher() ) ; 
      Cmp vr(*cmp[0]) ; 
      Cmp vt(*cmp[1]) ; 
      Cmp vp(*cmp[2]) ; 
      Cmp res1(*mp) ;
      Cmp res2(*mp) ;
      Cmp res3(*mp) ;
	    
      mp->comp_x_from_spherical(vr, vt, vp, res1) ; 
      mp->comp_y_from_spherical(vr, vt, vp, res2) ; 
      mp->comp_z_from_spherical(vr, vt, res3) ; 
	
      set(1) = res1 ;
      set(2) = res2 ;
      set(3) = res3 ;
      
    }// End of the spher -> cart case
  } // End of the case of cartesian new triad

  // ---------------------------------------------
  // Case where the new triad is a spherical one
  // ---------------------------------------------
  else {

    assert(nbvc == 0x0) ;

    // ---------------------------------
    //     Case cartesian -> spherical 
    // ---------------------------------
    if (bvc != 0x0) {	// The old triad is a cartesian one
      assert(bvs == 0x0) ; 
      
      // The triads should be the same as that associated 
      // with the mapping :
      assert( *nbvs == mp->get_bvect_spher() ) ; 
      assert( *bvc == mp->get_bvect_cart() ) ; 
      Cmp vx(*cmp[0]) ; 
      Cmp vy(*cmp[1]) ; 
      Cmp vz(*cmp[2]) ; 
      Cmp res1(*mp) ;
      Cmp res2(*mp) ;
      Cmp res3(*mp) ;
	    
      mp->comp_r_from_cartesian(vx, vy, vz, res1) ; 
      mp->comp_t_from_cartesian(vx, vy, vz, res2) ; 
      mp->comp_p_from_cartesian(vx, vy, res3) ; 

      set(1) = res1 ;
      set(2) = res2 ;
      set(3) = res3 ;
    }	// end of the  case cart -> spher

    
    // ------------------------------------
    //      Case spherical -> spherical
    // ------------------------------------
    if (bvs != 0x0) {	
      
      assert(bvc == 0x0) ; 
      
      cout << "Vector::change_triad : case not treated yet !" << endl ;
      abort() ; 
    }	// end of the spher->spher basis case

  } //  End of the case of spherical new triad

  triad = &new_triad ;

}
