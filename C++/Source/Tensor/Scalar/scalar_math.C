/*
 *  Mathematical functions for the Scalar class.
 *
 *  These functions are not member functions of the Scalar class.
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
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

char scalar_math_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2003/10/10 15:57:29  j_novak
 * Added the state one (ETATUN) to the class Scalar
 *
 * Revision 1.2  2003/10/01 13:04:44  e_gourgoulhon
 * The method Tensor::get_mp() returns now a reference (and not
 * a pointer) onto a mapping.
 *
 * Revision 1.1  2003/09/25 08:06:56  e_gourgoulhon
 * First versions (use Cmp as intermediate quantities).
 *
 *
 * $Header$
 *
 */

// Headers C
#include <math.h>

// Headers Lorene
#include "tensor.h"

			    //-------//
			    // Sine  //
			    //-------//

Scalar sin(const Scalar& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
		return ci ;
    }
    
    // Cas ETATUN
    if (ci.get_etat() == ETATUN) {
      Scalar resu(ci.get_mp()) ;
      resu = sin(double(1)) ;
      return resu ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Scalar co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = sin( ci.va ) ; 

    return co ;
}

			    //---------//
			    // Cosine  //
			    //---------//

Scalar cos(const Scalar& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    Scalar co(ci.get_mp()) ;		// result

    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
		co.set_etat_qcq() ; 
		co.va = double(1) ;
    }
    else {
      // Cas ETATUN
      if (ci.get_etat() == ETATUN) {
	co = cos(double(1)) ;
      }
      else {
		// Cas general
	assert(ci.get_etat() == ETATQCQ) ;	// sinon...
		co.set_etat_qcq() ; 
		co.va = cos( ci.va ) ; 
      }
    }

    return co ;
}

			    //----------//
			    // Tangent  //
			    //----------//

Scalar tan(const Scalar& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
		return ci ;
    }
    
    // Cas ETATUN
    if (ci.get_etat() == ETATUN) {
      Scalar resu(ci.get_mp()) ;
      resu = tan(double(1)) ;
      return resu ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Scalar co(ci.get_mp()) ;		// result
    co.set_etat_qcq() ; 
    co.va = tan( ci.va ) ; 

    return co ;
}

			    //----------//
			    // Arcsine  //
			    //----------//

Scalar asin(const Scalar& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	return ci ;
    }
    
    // Cas ETATUN
    if (ci.get_etat() == ETATUN) {
      Scalar resu(ci.get_mp()) ;
      resu = M_PI_2 ;
      return resu ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Scalar co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = asin( ci.va ) ; 

    return co ;
}

			    //-------------//
			    // Arccosine  //
			    //-------------//

Scalar acos(const Scalar& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    Scalar co(ci.get_mp()) ;		// result

    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	co.set_etat_qcq() ; 
	co.va = double(0.5) * M_PI ;
    }
    else {
      // Cas ETATUN
      if (ci.get_etat() == ETATUN) {
	co.set_etat_zero() ;
      }
      else {
	// Cas general
	assert(ci.get_etat() == ETATQCQ) ;	// sinon...
	co.set_etat_qcq() ; 
	co.va = acos( ci.va ) ; 
      }
    }

    return co ;
}

			    //-------------//
			    // Arctangent  //
			    //-------------//

Scalar atan(const Scalar& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	return ci ;
    }
    
    // Cas ETATUN
    if (ci.get_etat() == ETATUN) {
      Scalar resu(ci.get_mp()) ;
      resu = 0.25*M_PI ;
      return resu ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Scalar co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = atan( ci.va ) ; 

    return co ;
}

			    //-------------//
			    // Square root //
			    //-------------//

Scalar sqrt(const Scalar& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	return ci ;
    }
    
    // Cas ETATUN
    if (ci.get_etat() == ETATUN) {
	return ci ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Scalar co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = sqrt( ci.va ) ; 

    return co ;
}

			    //-------------//
			    // Cube root   //
			    //-------------//

Scalar racine_cubique(const Scalar& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	return ci ;
    }
    
    // Cas ETATUN
    if (ci.get_etat() == ETATUN) {
	return ci ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Scalar co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = racine_cubique( ci.va ) ; 

    return co ;
}

			    //--------------//
			    // Exponential  //
			    //--------------//

Scalar exp(const Scalar& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    Scalar co(ci.get_mp()) ;		// result

    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	co.set_etat_one() ; 
    }
    else {
      // Cas ETATUN
      if (ci.get_etat() == ETATUN) {
	co.set_etat_qcq() ;
	co = M_E ;
      }
      else {
	// Cas general
	assert(ci.get_etat() == ETATQCQ) ;	// sinon...
	co.set_etat_qcq() ; 
	co.va = exp( ci.va ) ; 
      }
    }

    return co ;
}

			    //---------------------//
			    // Neperian logarithm  //
			    //---------------------//

Scalar log(const Scalar& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	cout << "Argument of log is ZERO in log(Scalar) !" << endl ; 
	abort() ; 
    }

    // Cas ETATUN
    if (ci.get_etat() == ETATUN) {
      Scalar resu(ci.get_mp()) ;
      resu.set_etat_zero() ;
      return resu ;
    }

    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Scalar co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = log( ci.va ) ; 

    return co ;
}

			    //---------------------//
			    // Decimal  logarithm  //
			    //---------------------//

Scalar log10(const Scalar& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	cout << "Argument of log10 is ZERO in log10(Scalar) !" << endl ; 
	abort() ; 
    }

    // Cas ETATUN
    if (ci.get_etat() == ETATUN) {
      Scalar resu(ci.get_mp()) ;
      resu.set_etat_zero() ;
      return resu ;
    }

    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Scalar co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = log10( ci.va ) ; 

    return co ;
}

			    //------------------//
			    // Power (integer)  //
			    //------------------//

Scalar pow(const Scalar& ci, int n)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	if (n > 0) {
	    return ci ; 
	}
	else {
	    cout << "pow(Scalar, int) : ETATZERO^n with n <= 0 !" << endl ; 
	    abort() ;
	} 
    }

    // Cas ETATUN
    if (ci.get_etat() == ETATUN) {
      return ci ;
    }

    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Scalar co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = pow(ci.va,  n) ; 

    return co ;
}

			    //-----------------//
			    // Power (double)  //
			    //-----------------//

Scalar pow(const Scalar& ci, double x)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	if (x > double(0)) {
	    return ci ; 
	}
	else {
	    cout << "pow(Scalar, double) : ETATZERO^x with x <= 0 !" << endl ; 
	    abort() ;
	} 
    }

    // Cas ETATUN
    if (ci.get_etat() == ETATUN) {
      return ci ;
    }

    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Scalar co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = pow(ci.va, x) ; 

    return co ;
}

			    //-----------------//
			    // Absolute value  //
			    //-----------------//

Scalar abs(const Scalar& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	return ci ;
    }
    
    // Cas ETATUN
    if (ci.get_etat() == ETATUN) {
      return ci ;
    }

    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Scalar co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = abs( ci.va ) ; 

    return co ;
}

		    //-------------------------------//
    	    	    //            max                //
		    //-------------------------------//

Tbl max(const Scalar& ci) {

    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    Tbl resu( ci.get_mp().get_mg()->get_nzone() ) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {
      if (ci.get_etat() == ETATUN) {
	resu = 1 ; 
      }
      else {
	assert(ci.get_etat() == ETATQCQ) ;	
	
	resu = max( ci.va ) ;		// max(Valeur) 
      }
    }
   
    return resu ; 
}

		    //-------------------------------//
    	    	    //            min                //
		    //-------------------------------//

Tbl min(const Scalar& ci) {

    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    Tbl resu( ci.get_mp().get_mg()->get_nzone() ) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {
      if (ci.get_etat() == ETATUN) {
	resu = 1 ; 
      }
      else {
	assert(ci.get_etat() == ETATQCQ) ;	
	
	resu = min( ci.va ) ;		// min(Valeur) 	
      }
    }
          
    return resu ; 
}

		    //-------------------------------//
    	    	    //            norme              //
		    //-------------------------------//

Tbl norme(const Scalar& ci) {

    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    Tbl resu( ci.get_mp().get_mg()->get_nzone() ) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {
      if (ci.get_etat() == ETATUN) {
	resu = 1 ; 
      }
      else {
	assert(ci.get_etat() == ETATQCQ) ;	
	
	resu = norme( ci.va ) ;		// norme(Valeur) 	
      }
    }
          
    return resu ; 
}

		    //-------------------------------//
    	    	    //          diffrel              //
		    //-------------------------------//

Tbl diffrel(const Scalar& c1, const Scalar& c2) {
    
    // Protections
    assert(c1.get_etat() != ETATNONDEF) ;
    assert(c2.get_etat() != ETATNONDEF) ;

    int nz = c1.get_mp().get_mg()->get_nzone() ;
    Tbl resu(nz) ; 
    
    Scalar diff = c1 - c2 ;     // la compatibilite dzpuis est testee a ce niveau

    Tbl normdiff = norme(diff) ; 
    Tbl norme2 = norme(c2) ;
    
    assert(normdiff.get_etat() == ETATQCQ) ;     
    assert(norme2.get_etat() == ETATQCQ) ; 

    resu.set_etat_qcq() ; 
    for (int l=0; l<nz; l++) {
	if ( norme2(l) == double(0) ) {
	    resu.set(l) = normdiff(l) ; 
	}
	else{
	    resu.set(l) = normdiff(l) / norme2(l) ; 		    
	}
    }
    
    return resu ; 
    
}

			//-------------------------------//
			//          diffrelmax           //
		    //-------------------------------//

Tbl diffrelmax(const Scalar& c1, const Scalar& c2) {
    
    // Protections
    assert(c1.get_etat() != ETATNONDEF) ;
    assert(c2.get_etat() != ETATNONDEF) ;
    
    int nz = c1.get_mp().get_mg()->get_nzone() ;
    Tbl resu(nz) ; 
    
    Tbl max2 = max(abs(c2)) ;
    
    Scalar diff = c1 - c2 ;    // la compatibilite dzpuis est testee a ce niveau

    Tbl maxdiff = max(abs(diff)) ; 
    
    assert(maxdiff.get_etat() == ETATQCQ) ;     
    assert(max2.get_etat() == ETATQCQ) ; 

    resu.set_etat_qcq() ; 
    for (int l=0; l<nz; l++) {
	if ( max2(l) == double(0) ) {
	    resu.set(l) = maxdiff(l) ; 
	}
	else{
	    resu.set(l) = maxdiff(l) / max2(l) ; 		    
	}
    }
    
    return resu ; 
    
}

