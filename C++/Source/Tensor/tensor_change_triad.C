/*
 *  Methods for changing the triad of a Tensor
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

char tensor_change_triad_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/10/27 10:50:24  e_gourgoulhon
 * Added the case of a twice contravariant tensor in the assert.
 *
 * Revision 1.3  2003/10/06 14:25:51  j_novak
 * Added a test #ifndef... to prevent a warning
 *
 * Revision 1.2  2003/10/05 21:12:19  e_gourgoulhon
 * - Modified some assert.
 * - Corrected bug on index range in line 200.
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

void Tensor::change_triad(const Base_vect& new_triad) {

  assert (valence == 2) ;
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
	// ------------------------------
	
	set(0, 2) = - set(0, 2) ; // {xz} --> - {xz}
	set(1, 2) = - set(1, 2) ; // {yz} --> - {yz}
	set(2, 0) = - set(2, 0) ; // {zx} --> - {zx}
	set(2, 1) = - set(2, 1) ; // {zy} --> - {zy}
	// all other components are unchanged
	break ; 
      }
      case 0 : {	// the two basis have not a special relative orientation
			// -----------------------------------------------------
	cout << 
		"Tensor::change_basis : general value of rot_phi "
	     << " not contemplated yet, sorry !" << endl ;
	abort() ; 
	break ; 		
      }
	    
      default : {	    // error
	cout << 
	  "Tensor::change_basis : unexpected value of ind !" << endl ;
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

      // Only for double-covariant tensors or double-contravariant tensors
		assert( ( (type_indice(0)==COV) && (type_indice(1)==COV) ) ||
		        ( (type_indice(0)==CON) && (type_indice(1)==CON) ) ) ;  
      
      // Temporary storage of the components
      // the Base_vect *this is not valid...
      Tensor tmp(*mp, 2, COV, *triad) ;
      Cmp res1(*mp) ;
      Cmp res2(*mp) ;
      Cmp res3(*mp) ;
      for (int i=1; i<=3; i++) { //## Pb: les comp_?_from... utilisent des Cmp!
	Cmp cp1(set(1,i)) ;
	Cmp cp2(set(2,i)) ;
	Cmp cp3(set(3,i)) ;
	mp->comp_x_from_spherical(cp1, cp2, cp3, res1) ; 
	mp->comp_y_from_spherical(cp1, cp2, cp3, res2) ; 
	mp->comp_z_from_spherical(cp1, cp2, res3 ) ;
	tmp.set(1,i) = res1 ;
	tmp.set(2,i) = res2 ;
	tmp.set(3,i) = res3 ;
      }
      for (int i=1; i<=3; i++) {
	Cmp cp1(tmp(i,1)) ;
	Cmp cp2(tmp(i,2)) ;
	Cmp cp3(tmp(i,3)) ;
	mp->comp_x_from_spherical(cp1, cp2, cp3, res1) ;
	mp->comp_y_from_spherical(cp1, cp2, cp3, res2) ;
	mp->comp_z_from_spherical(cp1, cp2, res3) ;
	
	set(i,1) = res1 ;
	set(i,2) = res2 ;
	set(i,3) = res3 ;
      }	
      
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

      // Only for double-covariant tensors or double-contravariant tensors
		assert( ( (type_indice(0)==COV) && (type_indice(1)==COV) ) ||
		        ( (type_indice(0)==CON) && (type_indice(1)==CON) ) ) ;  
      
      // Temporary storage of the components
      Tensor tmp(*mp, 2, COV, *triad) ;
      Cmp res1(*mp) ;
      Cmp res2(*mp) ;
      Cmp res3(*mp) ;
      for (int i=1; i<=3; i++) {
	Cmp cp1(set(1,i)) ;
	Cmp cp2(set(2,i)) ;
	Cmp cp3(set(3,i)) ;
	mp->comp_r_from_cartesian(cp1, cp2, cp3, res1) ; 
	mp->comp_t_from_cartesian(cp1, cp2, cp3, res2) ; 
	mp->comp_p_from_cartesian(cp1, cp2, res3) ;
	tmp.set(1,i) = res1 ;
	tmp.set(2,i) = res2 ;
	tmp.set(3,i) = res3 ;
      }
      for (int i=1; i<=3; i++) {
	Cmp cp1(tmp(i,1)) ;
	Cmp cp2(tmp(i,2)) ;
	Cmp cp3(tmp(i,3)) ;
	mp->comp_r_from_cartesian(cp1, cp2, cp3, res1) ; 
	mp->comp_t_from_cartesian(cp1, cp2, cp3, res2) ; 
	mp->comp_p_from_cartesian(cp1, cp2, res3) ;

 	set(i,1) = res1 ;
	set(i,2) = res2 ;
	set(i,3) = res3 ;
      }
    }	// end of the  case cart -> spher


    // ------------------------------------
    //      Case spherical -> spherical
    // ------------------------------------
    if (bvs != 0x0) {	
      
      assert(bvc == 0x0) ; 
      
      cout << "Tensor::change_triad : case not treated yet !" << endl ;
      abort() ; 
    }	// end of the spher->spher basis case

  } //  End of the case of spherical new triad

  triad = &new_triad ;

}
