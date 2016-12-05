/*
 *  Member functions of the class Scalar for various r manipulations
 *
 *    See file scalar.h for documentation. 
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Eric Gourgoulhon  (for a preceding Cmp version)
 *   Copyright (c) 1999-2001 Philippe Grandclement  (for a preceding Cmp version)
 *   Copyright (c) 2001 Jerome Novak (for a preceding Cmp version)
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.24  2016/12/05 16:18:19  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.23  2014/10/13 08:53:47  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.22  2005/05/25 16:11:05  j_novak
 * Better handling of the case with no compactified domain.
 *
 * Revision 1.21  2005/05/18 11:45:46  j_novak
 * Added del_deriv() calls at the end of inc/dec_dzpuis.
 *
 * Revision 1.20  2004/10/11 15:09:04  j_novak
 * The radial manipulation functions take Scalar as arguments, instead of Cmp.
 * Added a conversion operator from Scalar to Cmp.
 * The Cmp radial manipulation function make conversion to Scalar, call to the
 * Map_radial version with a Scalar argument and back.
 *
 * Revision 1.19  2004/10/08 13:34:37  j_novak
 * Scalar::div_r() does not need to pass through Cmp version anymore.
 *
 * Revision 1.18  2004/05/12 14:22:12  f_limousin
 * Treated the case ETATZERO in dec_dzpuis and inc_dzpuis (-> return).
 *
 * Revision 1.17  2004/01/29 09:31:44  j_novak
 * Better treatment of ETATUN
 *
 * Revision 1.16  2004/01/28 10:34:27  j_novak
 * Corrected some errors.
 *
 * Revision 1.15  2004/01/27 15:10:02  j_novak
 * New methods Scalar::div_r_dzpuis(int) and Scalar_mult_r_dzpuis(int)
 * which replace div_r_inc*. Tried to clean the dzpuis handling.
 * WARNING: no testing at this point!!
 *
 * Revision 1.14  2004/01/22 16:11:30  e_gourgoulhon
 * Added (provisory method) div_r_inc1().
 * Case inc = 3 treated in inc_dzpuis.
 *
 * Revision 1.13  2003/11/04 23:02:21  e_gourgoulhon
 * -- Method dec_dzpuis(int decrem) : the case decrem = 1 is now treated.
 * -- Method div_tant() is now defined in file scalar_th_manip.C.
 *
 * Revision 1.12  2003/11/03 22:36:37  e_gourgoulhon
 *  Method dec_dzpuis: changed the name of argument dec --> decrem
 *  (in order not to shadow some globally defined dec).
 *
 * Revision 1.11  2003/11/03 10:25:05  e_gourgoulhon
 * Scalar::dec_dzpuis(int dec)  : treated the case dec = 0
 * Scalar::inc_dzpuis(int inc)  : treated the case inc = 0
 *
 * Revision 1.10  2003/10/29 11:02:13  e_gourgoulhon
 * Functions dec_dzpuis and inc_dzpuis have now an integer argument to
 * specify by which amount dzpuis is to be increased.
 * Accordingly methods dec2_dzpuis and inc2_dzpuis have been suppressed
 *
 * Revision 1.9  2003/10/16 15:29:08  e_gourgoulhon
 * Name of method div_r_ced() changed to div_r_inc2().
 * Name of method div_rsint_ced() changed to div_rsint_inc2().
 *
 * Revision 1.8  2003/10/15 10:44:21  e_gourgoulhon
 * Added new method div_tant.
 *
 * Revision 1.7  2003/10/13 20:50:16  e_gourgoulhon
 * Replaced "va.base" by "va.get_base()".
 *
 * Revision 1.6  2003/10/10 15:57:29  j_novak
 * Added the state one (ETATUN) to the class Scalar
 *
 * Revision 1.5  2003/10/08 14:24:10  j_novak
 * replaced mult_r_zec with mult_r_ced
 *
 * Revision 1.4  2003/10/08 12:26:03  j_novak
 * Second part of the bug (sorry!)
 *
 * Revision 1.3  2003/10/08 12:19:12  j_novak
 * Bug corrected, thanks to purify
 *
 * Revision 1.2  2003/10/05 21:16:41  e_gourgoulhon
 * Added methods div_r_ced() and div_rsint_ced().
 *
 * Revision 1.1  2003/09/25 09:12:01  e_gourgoulhon
 * First version (uses Cmp as intermediate variable).
 *
 *
 * $Header$
 *
 */

#include "tensor.h" 
#include "cmp.h" 


			//-------------------//
			//	 div_r       //
			//-------------------//

namespace Lorene {
void Scalar::div_r() {
    
  mp->div_r(*this) ;   // Call of the appropriate routine of the mapping
    
  del_deriv() ;   // Delete the derived members

}


			//---------------------//
			//    div_r_dzpuis     //
			//---------------------//


void Scalar::div_r_dzpuis(int ced_mult_r) {
    
  assert (etat != ETATNONDEF) ;

  Base_val lbase = va.base ;
 
  if (etat != ETATZERO) {
	
    assert((etat == ETATQCQ) || (etat == ETATUN)) ;
    
    set_etat_qcq() ;

    int nzm1 = mp->get_mg()->get_nzone() - 1 ; // index of the CED
    bool ced = mp->get_mg()->get_type_r(nzm1) == UNSURR ;
	
    // Copy of the CED part of *this into uu_ext 
    Scalar uu_ext(*mp) ;
    if (ced) {
	uu_ext.allocate_all() ;
	uu_ext.annule(0,nzm1-1) ; // zero in all domains but the CED
	uu_ext.set_domain(nzm1) = domain(nzm1) ; 
	uu_ext.set_spectral_base(va.get_base()) ; 

	// Division by r in all domains but the CED
	annule(nzm1, nzm1) ; 	// zero in the CED
    }
    div_r() ; 
	
    if (ced) { // Add the CED part
	set_domain(nzm1) = uu_ext.domain(nzm1) ; 
	dzpuis += 1 ; 
    }
  }
	
  int diff_dzpuis = ced_mult_r - dzpuis ;
  (diff_dzpuis >= 0) ? inc_dzpuis(diff_dzpuis) : dec_dzpuis(-diff_dzpuis) ;
  
  lbase.sx() ;
  set_spectral_base(lbase) ;

  del_deriv() ;   // Delete the derived members
    
  return ;
}

			//---------------------------//
			//	    div_r_ced	     //
			//---------------------------//

void Scalar::div_r_ced() {

  mp->div_r_zec(*this) ;   // Call of the appropriate routine of the mapping
    
  del_deriv() ;   // Delete the derived members

}

			//---------------------------//
			//	    mult_r	     //
			//---------------------------//

void Scalar::mult_r() {
    
  mp->mult_r(*this) ;   // Call of the appropriate routine of the mapping
    
  del_deriv() ;   // Delete the derived members
    
}

			//---------------------//
			//    mult_r_dzpuis    //
			//---------------------//


void Scalar::mult_r_dzpuis(int ced_mult_r) {
    
  assert (etat != ETATNONDEF) ;

  Base_val lbase = va.base ;

  if (etat != ETATZERO) {
	
    assert((etat == ETATQCQ) || (etat == ETATUN)) ;
    
    set_etat_qcq() ;
	
    int nzm1 = mp->get_mg()->get_nzone() - 1 ; // index of the CED
    bool ced = mp->get_mg()->get_type_r(nzm1) == UNSURR ;
    Scalar uu_ext(*mp) ; 
	
    if (ced) {// Copy of the CED part of *this into uu_ext 
	uu_ext.allocate_all() ;
	uu_ext.annule(0,nzm1-1) ; // zero in all domains but the CED
	uu_ext.set_domain(nzm1) = domain(nzm1) ; 
	uu_ext.set_spectral_base(va.get_base()) ; 
	
	// Division by r in all domains but the CED
	annule(nzm1, nzm1) ; 	// zero in the CED
    }
    mult_r() ; 
	
    if (ced) {// Add the CED part
	set_domain(nzm1) = uu_ext.domain(nzm1) ; 
	dzpuis -= 1 ; 
    }
  }
	
  int diff_dzpuis = ced_mult_r - dzpuis ;
  (diff_dzpuis >= 0) ? inc_dzpuis(diff_dzpuis) : dec_dzpuis(-diff_dzpuis) ;
  
  lbase.mult_x() ;
  set_spectral_base(lbase) ;

  del_deriv() ;   // Delete the derived members
    
  return ;
}

			//---------------------------//
			//	    mult_r_ced	     //
			//---------------------------//

void Scalar::mult_r_ced() {
    
  mp->mult_r_zec(*this) ;   // Call of the appropriate routine of the mapping
    
  del_deriv() ;   // Delete the derived members

}

			//---------------------------//
			//	    mult_rsint	     //
			//---------------------------//

void Scalar::mult_rsint() {
    
  mp->mult_rsint(*this) ; // Call of the appropriate routine of the mapping 
    
  del_deriv() ;   // Delete the derived members

}

			//-------------------------//
			//    mult_rsint_dzpuis    //
			//-------------------------//


void Scalar::mult_rsint_dzpuis(int ced_mult_r) {
    
  assert (etat != ETATNONDEF) ;

  Base_val lbase = va.base ;
  int nzm1 = mp->get_mg()->get_nzone() - 1 ; // index of the CED
  bool ced = mp->get_mg()->get_type_r(nzm1) == UNSURR ;

  if (etat != ETATZERO) {

    assert((etat == ETATQCQ) || (etat == ETATUN)) ;
    
    set_etat_qcq() ;
    Scalar uu_ext(*mp) ; 
		
    if (ced) {// Copy of the CED part of *this into uu_ext and multiplication by r
	uu_ext.allocate_all() ;
	uu_ext.annule(0,nzm1-1) ; // zero in all domains but the CED
	uu_ext.set_domain(nzm1) = domain(nzm1) ; 
	uu_ext.set_spectral_base(va.get_base()) ; 
	int dzpuis_ext = dzpuis - 1 ;
	uu_ext.set_dzpuis(dzpuis_ext) ;
	int diff_dzpuis = ced_mult_r - dzpuis_ext ;
	(diff_dzpuis >= 0) ? uu_ext.inc_dzpuis(diff_dzpuis) 
	    : uu_ext.dec_dzpuis(-diff_dzpuis) ;

	// Multiplication by sin(theta) in the CED :
	// what follows does not apply if the mapping is not radial:
	assert( dynamic_cast<const Map_radial*>(mp) != 0x0 ) ; 
	uu_ext.mult_sint() ;
	
	// Multiplication by r sin(theta) in all domains but the CED
	annule(nzm1, nzm1) ; 	// zero in the CED
    }
    mult_rsint() ; 
	
    if (ced) // Add the CED part
	set_domain(nzm1) = uu_ext.domain(nzm1) ; 
  }
  
  if (ced) dzpuis = ced_mult_r ;
 
  lbase.mult_x() ;
  lbase.mult_sint() ; 
  set_spectral_base(lbase) ;

  del_deriv() ;   // Delete the derived members

}

			//---------------------------//
			//	    div_rsint	     //
			//---------------------------//

void Scalar::div_rsint() {

  mp->div_rsint(*this) ;   // Call of the appropriate routine of the mapping
    
  del_deriv() ;   // Delete the derived members

}


			//-------------------------//
			//     div_rsint_dzpuis    //
			//-------------------------//


void Scalar::div_rsint_dzpuis(int ced_mult_r) {
    
  assert (etat != ETATNONDEF) ;

  Base_val lbase = va.base ;
  int nzm1 = mp->get_mg()->get_nzone() - 1 ; // index of the CED
  bool ced = mp->get_mg()->get_type_r(nzm1) == UNSURR ;

  if (etat != ETATZERO) {

    assert((etat == ETATQCQ) || (etat == ETATUN)) ;
    
    set_etat_qcq() ;
    Scalar uu_ext(*mp) ; 
	
    if (ced) {// Copy of the CED part of *this into uu_ext and multiplication by r
	uu_ext.allocate_all() ;
	uu_ext.annule(0,nzm1-1) ; // zero in all domains but the CED
	uu_ext.set_domain(nzm1) = domain(nzm1) ; 
	uu_ext.set_spectral_base(va.get_base()) ; 
	int dzpuis_ext = dzpuis + 1 ;
	uu_ext.set_dzpuis(dzpuis_ext) ;
	int diff_dzpuis = ced_mult_r - dzpuis_ext ;
	(diff_dzpuis >= 0) ? uu_ext.inc_dzpuis(diff_dzpuis) 
	    : uu_ext.dec_dzpuis(-diff_dzpuis) ;
	
	// Division by sin(theta) in the CED :
	// what follows does not apply if the mapping is not radial:
	assert( dynamic_cast<const Map_radial*>(mp) != 0x0 ) ; 
	uu_ext.div_sint() ;

	// Division by r sin(theta) in all domains but the CED
	annule(nzm1, nzm1) ; 	// zero in the CED
    }
    div_rsint() ; 
	
    if (ced) // Add the CED part
	set_domain(nzm1) = uu_ext.domain(nzm1) ; 
	
  }
  if (ced) 
      dzpuis = ced_mult_r ;
	
  lbase.sx() ;
  lbase.ssint() ;
  set_spectral_base(lbase) ;

  del_deriv() ;   // Delete the derived members

}


			//-----------------------//
			//      dec_dzpuis       //
			//-----------------------//

void Scalar::dec_dzpuis(int decrem) {
    
	if (decrem == 0) return ; 
	
	if (etat == ETATZERO) {
//##	    dzpuis -= decrem ;   // a voir...
	    return ; 
	}

	switch (decrem) {
	
		case 0 : { 
   			break ; 
		}

		case 1 : { 
			mp->dec_dzpuis(*this) ;   
    		break ; 
		}

		case 2 : {
    		mp->dec2_dzpuis(*this) ;  
			break ; 
		}
		
		case 3 : {
    		mp->dec2_dzpuis(*this) ;  
		mp->dec_dzpuis(*this) ;   
			break ; 
		}
		
		case 4 : {
    		mp->dec2_dzpuis(*this) ;  
    		mp->dec2_dzpuis(*this) ;  
			break ; 
		}
		
		default : {
			cout << "Scalar::dec_dzpuis : unexpected value of decrem !"
				<< endl << "  decrem = " << decrem << endl ; 
			abort() ;
			break ; 
		}
	}
	del_deriv() ;
	
}

			//-----------------------//
			//      inc_dzpuis       //
			//-----------------------//

void Scalar::inc_dzpuis(int inc) {
    
	if (inc == 0) return ; 

	if (etat == ETATZERO) {
//##	    dzpuis += inc ;   // a voir ...
	    return ; 
	}

	switch (inc) {
	
		case 0 : { 
   		    break ; 
		}

		case 1 : { 
		    mp->inc_dzpuis(*this) ;   
    		    break ; 
		}

		case 2 : {
    		    mp->inc2_dzpuis(*this) ;  
		    break ; 
		}
		
		case 3 : {
		    mp->inc_dzpuis(*this) ;   
    		    mp->inc2_dzpuis(*this) ;  
		    break ; 
		}
		
		case 4 : {
    		    mp->inc2_dzpuis(*this) ;  
    		    mp->inc2_dzpuis(*this) ;  
		    break ; 
		}
		
		default : {
			cout << "Scalar::inc_dzpuis : unexpected value of inc !"
				<< endl << "  inc = " << inc << endl ; 
			abort() ;
			break ; 
		}
	}
	del_deriv() ;
	
}



}
