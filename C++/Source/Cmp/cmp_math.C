/*
 *  Mathematical functions for the Cmp class.
 *
 *  These functions are not member functions of the Cmp class.
 *
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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

char cmp_math_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.2  1999/12/02  17:59:16  phil
 * /
 *
 * Revision 2.1  1999/11/26  14:23:30  eric
 * Modif commentaires.
 *
 * Revision 2.0  1999/11/15  14:13:05  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <math.h>

// Headers Lorene
#include "cmp.h"

			    //-------//
			    // Sine  //
			    //-------//

Cmp sin(const Cmp& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	return ci ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Cmp co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = sin( ci.va ) ; 

    return co ;
}

			    //---------//
			    // Cosine  //
			    //---------//

Cmp cos(const Cmp& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    Cmp co(ci.get_mp()) ;		// result

    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	co.set_etat_qcq() ; 
	co.va = double(1) ;
    }
    else {
		// Cas general
	assert(ci.get_etat() == ETATQCQ) ;	// sinon...
	co.set_etat_qcq() ; 
	co.va = cos( ci.va ) ; 
    }

    return co ;
}

			    //----------//
			    // Tangent  //
			    //----------//

Cmp tan(const Cmp& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	return ci ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Cmp co(ci.get_mp()) ;		// result
    co.set_etat_qcq() ; 
    co.va = tan( ci.va ) ; 

    return co ;
}

			    //----------//
			    // Arcsine  //
			    //----------//

Cmp asin(const Cmp& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	return ci ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Cmp co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = asin( ci.va ) ; 

    return co ;
}

			    //-------------//
			    // Arccossine  //
			    //-------------//

Cmp acos(const Cmp& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    Cmp co(ci.get_mp()) ;		// result

    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	co.set_etat_qcq() ; 
	co.va = double(0.5) * M_PI ;
    }
    else {
		// Cas general
	assert(ci.get_etat() == ETATQCQ) ;	// sinon...
	co.set_etat_qcq() ; 
	co.va = acos( ci.va ) ; 
    }

    return co ;
}

			    //-------------//
			    // Arctangent  //
			    //-------------//

Cmp atan(const Cmp& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	return ci ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Cmp co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = atan( ci.va ) ; 

    return co ;
}

			    //-------------//
			    // Square root //
			    //-------------//

Cmp sqrt(const Cmp& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	return ci ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Cmp co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = sqrt( ci.va ) ; 

    return co ;
}

			    //-------------//
			    // Cube root   //
			    //-------------//

Cmp racine_cubique(const Cmp& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	return ci ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Cmp co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = racine_cubique( ci.va ) ; 

    return co ;
}

			    //--------------//
			    // Exponential  //
			    //--------------//

Cmp exp(const Cmp& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    Cmp co(ci.get_mp()) ;		// result

    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	co.set_etat_qcq() ; 
	co.va = double(1) ;
    }
    else {
		// Cas general
	assert(ci.get_etat() == ETATQCQ) ;	// sinon...
	co.set_etat_qcq() ; 
	co.va = exp( ci.va ) ; 
    }

    return co ;
}

			    //---------------------//
			    // Neperian logarithm  //
			    //---------------------//

Cmp log(const Cmp& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	cout << "Argument of log is ZERO in log(Cmp) !" << endl ; 
	abort() ; 
    }

    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Cmp co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = log( ci.va ) ; 

    return co ;
}

			    //---------------------//
			    // Decimal  logarithm  //
			    //---------------------//

Cmp log10(const Cmp& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	cout << "Argument of log10 is ZERO in log10(Cmp) !" << endl ; 
	abort() ; 
    }

    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Cmp co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = log10( ci.va ) ; 

    return co ;
}

			    //------------------//
			    // Power (integer)  //
			    //------------------//

Cmp pow(const Cmp& ci, int n)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	if (n > 0) {
	    return ci ; 
	}
	else {
	    cout << "pow(Cmp, int) : ETATZERO^n with n <= 0 !" << endl ; 
	    abort() ;
	} 
    }

    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Cmp co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = pow(ci.va,  n) ; 

    return co ;
}

			    //-----------------//
			    // Power (double)  //
			    //-----------------//

Cmp pow(const Cmp& ci, double x)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	if (x > double(0)) {
	    return ci ; 
	}
	else {
	    cout << "pow(Cmp, double) : ETATZERO^x with x <= 0 !" << endl ; 
	    abort() ;
	} 
    }

    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Cmp co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = pow(ci.va, x) ; 

    return co ;
}

			    //-----------------//
			    // Absolute value  //
			    //-----------------//

Cmp abs(const Cmp& ci)
{
    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ci.get_etat() == ETATZERO) {
	return ci ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...

    Cmp co(ci.get_mp()) ;		// result

    co.set_etat_qcq() ; 
    co.va = abs( ci.va ) ; 

    return co ;
}

		    //-------------------------------//
    	    	    //            max                //
		    //-------------------------------//

Tbl max(const Cmp& ci) {

    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    Tbl resu( ci.get_mp()->get_mg()->get_nzone() ) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {
	assert(ci.get_etat() == ETATQCQ) ;	

	resu = max( ci.va ) ;		// max(Valeur) 
    }
          
    return resu ; 
}

		    //-------------------------------//
    	    	    //            min                //
		    //-------------------------------//

Tbl min(const Cmp& ci) {

    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    Tbl resu( ci.get_mp()->get_mg()->get_nzone() ) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {
	assert(ci.get_etat() == ETATQCQ) ;	

	resu = min( ci.va ) ;		// min(Valeur) 	
    }
          
    return resu ; 
}

		    //-------------------------------//
    	    	    //            norme              //
		    //-------------------------------//

Tbl norme(const Cmp& ci) {

    // Protection
    assert(ci.get_etat() != ETATNONDEF) ;
    
    Tbl resu( ci.get_mp()->get_mg()->get_nzone() ) ; 
    
    if (ci.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {
	assert(ci.get_etat() == ETATQCQ) ;	

	resu = norme( ci.va ) ;		// norme(Valeur) 	
    }
          
    return resu ; 
}

		    //-------------------------------//
    	    	    //          diffrel              //
		    //-------------------------------//

Tbl diffrel(const Cmp& c1, const Cmp& c2) {
    
    // Protections
    assert(c1.get_etat() != ETATNONDEF) ;
    assert(c2.get_etat() != ETATNONDEF) ;

    int nz = c1.get_mp()->get_mg()->get_nzone() ;
    Tbl resu(nz) ; 
    
    Cmp diff = c1 - c2 ;     // la compatibilite dzpuis est testee a ce niveau

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

Tbl diffrelmax(const Cmp& c1, const Cmp& c2) {
    
    // Protections
    assert(c1.get_etat() != ETATNONDEF) ;
    assert(c2.get_etat() != ETATNONDEF) ;
    
    int nz = c1.get_mp()->get_mg()->get_nzone() ;
    Tbl resu(nz) ; 
    
    Tbl max2 = max(abs(c2)) ;
    
    Cmp diff = c1 - c2 ;    // la compatibilite dzpuis est testee a ce niveau

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

