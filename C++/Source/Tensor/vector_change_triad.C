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

 


/*
 * $Id$
 * $Log$
 * Revision 1.10  2016/12/05 16:18:18  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.9  2014/10/13 08:53:44  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2014/10/06 15:13:20  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2010/01/26 16:47:10  e_gourgoulhon
 * Suppressed the assert on np >= 4 (too strong ?).
 *
 * Revision 1.6  2005/09/15 15:51:26  j_novak
 * The "rotation" (change of triad) methods take now Scalars as default
 * arguments.
 *
 * Revision 1.5  2005/02/03 14:34:19  f_limousin
 * Improvement of the case Cartesian --> Cartesian to be consistent
 * with Tensor::change_triad(Base_vect&).
 *
 * Revision 1.4  2003/10/28 21:27:54  e_gourgoulhon
 * -- Read-only access to the components performed by operator()(int) instead of
 *     set(int ).
 * -- Corrected index range in the Cartesian->Cartesian case.
 *
 * Revision 1.3  2003/10/06 14:25:51  j_novak
 * Added a test #ifndef... to prevent a warning
 *
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
#include <cassert>

// Lorene headers
#include "tensor.h"

namespace Lorene {
void Vector::change_triad(const Base_vect& new_triad) {

  assert(triad != 0x0) ; 
  
  const Base_vect_cart* nbvc = dynamic_cast<const Base_vect_cart*>(&new_triad) ; 
#ifndef NDEBUG
  const Base_vect_spher* nbvs 
    = dynamic_cast<const Base_vect_spher*>(&new_triad) ; 
#endif

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
	  Vector copie (*this) ;
	  
	set(1) = - copie(1) ;	 // V^x --> - V^x
	set(2) = - copie(2) ; 	 // V^y --> - V^y
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
#ifndef NDEBUG
      int nz = mp->get_mg()->get_nzone() ;
      for (int i=0; i<nz; i++) {
//##	  assert( mp->get_mg()->get_np(i) >= 4) ;
	  assert( mp->get_mg()->get_nt(i) >= 5) ;
      }
#endif
      Scalar res1(*mp) ;
      Scalar res2(*mp) ;
	    
      mp->comp_x_from_spherical(*cmp[0], *cmp[1], *cmp[2], res1) ; 
      mp->comp_y_from_spherical(*cmp[0], *cmp[1], *cmp[2], res2) ; 
      mp->comp_z_from_spherical(*cmp[0], *cmp[1], set(3)) ; 
	
      set(1) = res1 ;
      set(2) = res2 ;
      
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
#ifndef NDEBUG
      int nz = mp->get_mg()->get_nzone() ;
      for (int i=0; i<nz; i++) {
//##	  assert( mp->get_mg()->get_np(i) >= 4) ;
	  assert( mp->get_mg()->get_nt(i) >= 5) ;
      }
#endif
      Scalar res1(*mp) ;
      Scalar res2(*mp) ;
	    
      mp->comp_r_from_cartesian(*cmp[0], *cmp[1], *cmp[2], res1) ; 
      mp->comp_t_from_cartesian(*cmp[0], *cmp[1], *cmp[2], res2) ; 
      mp->comp_p_from_cartesian(*cmp[0], *cmp[1], set(3)) ; 

      set(1) = res1 ;
      set(2) = res2 ;
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
}
