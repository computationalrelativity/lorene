/*
 *  Mathematical functions for class Mtbl
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char mtbl_math_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.4  1999/12/02  17:55:03  phil
 * *** empty log message ***
 *
 * Revision 2.3  1999/10/29  15:46:35  eric
 * *** empty log message ***
 *
 * Revision 2.2  1999/10/29  15:06:38  eric
 * Ajout de fonctions mathematiques (abs, norme, etc...).
 *
 * Revision 2.1  1999/03/01  15:09:56  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/24  15:25:41  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
// ---------
#include <math.h>
#include <stdlib.h>

// Headers Lorene
// --------------
#include "mtbl.h"

			    //-------//
			    // Sinus //
			    //-------//

Mtbl sin(const Mtbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;		// sinon...
    Mtbl to(ti.get_mg()) ;			// Mtbl resultat
    to.set_etat_qcq() ;
    int nzone = ti.get_nzone() ;
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = sin( *(ti.t[i]) ) ;
    }
    return to ;
}

			    //---------//
			    // Cosinus //
			    //---------//

Mtbl cos(const Mtbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    Mtbl to(ti.get_mg()) ;			// Mtbl resultat
    to.set_etat_qcq() ;
    int nzone = ti.get_nzone() ;

    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	for (int i=0 ; i<nzone ; i++) {
	    *(to.t[i]) = 1. ;
	}
	return to ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = cos( *(ti.t[i]) ) ;
    }
    return to ;
}

			    //----------//
			    // Tangente //
			    //----------//

Mtbl tan(const Mtbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Mtbl to(ti.get_mg()) ;			// Mtbl resultat
    to.set_etat_qcq() ;
    int nzone = ti.get_nzone() ;
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = tan( *(ti.t[i]) ) ;
    }
    return to ;
}

			    //----------//
			    // ArcSinus //
			    //----------//

Mtbl asin(const Mtbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Mtbl to(ti.get_mg()) ;			// Mtbl resultat
    to.set_etat_qcq() ;
    int nzone = ti.get_nzone() ;
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = asin( *(ti.t[i]) ) ;
    }
    return to ;
}

			    //------------//
			    // ArcCosinus //
			    //------------//

Mtbl acos(const Mtbl& ti)
{
   // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    Mtbl to(ti.get_mg()) ;			// Mtbl resultat
    to.set_etat_qcq() ;
    int nzone = ti.get_nzone() ;

    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	for (int i=0 ; i<nzone ; i++) {
	    *(to.t[i]) = M_PI * .5 ;
	}
	return to ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = acos( *(ti.t[i]) ) ;
    }
    return to ;
}

			    //-------------//
			    // ArcTangente //
			    //-------------//

Mtbl atan(const Mtbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Mtbl to(ti.get_mg()) ;			// Mtbl resultat
    to.set_etat_qcq() ;
    int nzone = ti.get_nzone() ;
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = atan( *(ti.t[i]) ) ;
    }
    return to ;
}

			    //------//
			    // Sqrt //
			    //------//

Mtbl sqrt(const Mtbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Mtbl to(ti.get_mg()) ;			// Mtbl resultat
    to.set_etat_qcq() ;
    int nzone = ti.get_nzone() ;
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = sqrt( *(ti.t[i]) ) ;
    }
    return to ;
}


			    //------//
			    // Cubic //
			    //------//

Mtbl racine_cubique(const Mtbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Mtbl to(ti.get_mg()) ;			// Mtbl resultat
    to.set_etat_qcq() ;
    int nzone = ti.get_nzone() ;
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = racine_cubique( *(ti.t[i]) ) ;
    }
    return to ;
}
			    //---------------//
			    // Exponantielle //
			    //---------------//

Mtbl exp(const Mtbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    Mtbl to(ti.get_mg()) ;			// Mtbl resultat
    to.set_etat_qcq() ;
    int nzone = ti.get_nzone() ;

    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	for (int i=0 ; i<nzone ; i++) {
	    *(to.t[i]) = 1. ;
	}
	return to ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = exp( *(ti.t[i]) ) ;
    }
    return to ;
}

			    //-------------//
			    // Log naturel //
			    //-------------//

Mtbl log(const Mtbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	cout << "Mtbl log: log(ETATZERO) !" << endl  ;
	abort () ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Mtbl to(ti.get_mg()) ;		// Mtbl resultat
    to.set_etat_qcq() ;
    int nzone = ti.get_nzone() ;
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = log( *(ti.t[i]) ) ;
    }
    return to ;
}

			    //-------------//
			    // Log decimal //
			    //-------------//

Mtbl log10(const Mtbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	cout << "Mtbl log10: log10(ETATZERO) !" << endl ;
	abort () ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Mtbl to(ti.get_mg()) ;		// Mtbl resultat
    to.set_etat_qcq() ;
    int nzone = ti.get_nzone() ;
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = log10( *(ti.t[i]) ) ;
    }
    return to ;
}

			    //--------------//
			    // Power entier //
			    //--------------//

Mtbl pow(const Mtbl& ti, int n)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	if (n > 0) {
	    return ti ;
	}
	else {
	    cout << "Mtbl pow: ETATZERO^n avec n<=0 ! "<< endl  ;
	    abort () ;
	}
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Mtbl to(ti.get_mg()) ;		// Mtbl resultat
    to.set_etat_qcq() ;
    double x = n ;
    int nzone = ti.get_nzone() ;
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = pow( *(ti.t[i]), x ) ;
    }
    return to ;
}

			    //--------------//
			    // Power double //
			    //--------------//

Mtbl pow(const Mtbl& ti, double x)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	if (x > 0) {
	    return ti ;
	}
	else {
	    cout << "Mtbl pow: ETATZERO^x avec x<=0 !" << endl ;
	    abort () ;
	}
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Mtbl to(ti.get_mg()) ;			// Mtbl resultat
    to.set_etat_qcq() ;
    int nzone = ti.get_nzone() ;
    for (int i=0 ; i<nzone ; i++) {
	*(to.t[i]) = pow( *(ti.t[i]), x ) ;
    }
    return to ;
}

			    //----------------//
			    // Absolute value //
			    //----------------//

Mtbl abs(const Mtbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general

    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Mtbl to(ti.get_mg()) ;			// Mtbl resultat

    to.set_etat_qcq() ;

    int nzone = ti.get_nzone() ;

    for (int l=0 ; l<nzone ; l++) {
	*(to.t[l]) = abs( *(ti.t[l]) ) ;	
    }

    return to ;
}




		    //-------------------------------//
    	    	    //            max                //
		    //-------------------------------//

Tbl max(const Mtbl& mti) {

    // Protection
    assert(mti.get_etat() != ETATNONDEF) ;
    
    int nz = mti.get_nzone() ;
    
    Tbl resu(nz) ; 
    
    if (mti.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {  	// Cas general

	assert(mti.get_etat() == ETATQCQ) ;     // sinon....
    
	resu.set_etat_qcq() ; 
	for (int l=0 ; l<nz ; l++) {
	    resu.set(l) = max( *(mti.t[l]) ) ;
	}
    }
     
    return resu ; 
}

		    //-------------------------------//
    	    	    //            min                //
		    //-------------------------------//

Tbl min(const Mtbl& mti) {

    // Protection
    assert(mti.get_etat() != ETATNONDEF) ;
    
    int nz = mti.get_nzone() ;
    
    Tbl resu(nz) ; 
    
    if (mti.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {  	// Cas general

	assert(mti.get_etat() == ETATQCQ) ;     // sinon....
    
	resu.set_etat_qcq() ; 
	for (int l=0 ; l<nz ; l++) {
	    resu.set(l) = min( *(mti.t[l]) ) ;
	}
    }
     
    return resu ; 
}

		    //-------------------------------//
    	    	    //            norme              //
		    //-------------------------------//

Tbl norme(const Mtbl& mti) {

    // Protection
    assert(mti.get_etat() != ETATNONDEF) ;
    
    int nz = mti.get_nzone() ;
    
    Tbl resu(nz) ; 
    
    if (mti.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {  	// Cas general

	assert(mti.get_etat() == ETATQCQ) ;     // sinon....
    
	resu.set_etat_qcq() ; 
	for (int l=0 ; l<nz ; l++) {
	    resu.set(l) = norme( *(mti.t[l]) ) ;
	}
    }
     
    return resu ; 
}

		    //-------------------------------//
    	    	    //          diffrel              //
		    //-------------------------------//

Tbl diffrel(const Mtbl& mt1, const Mtbl& mt2) {
    
    // Protections
    assert(mt1.get_etat() != ETATNONDEF) ;
    assert(mt2.get_etat() != ETATNONDEF) ;
    
    int nz = mt1.get_nzone() ;
    Tbl resu(nz) ; 
    
    Tbl normdiff = norme(mt1 - mt2) ; 
    Tbl norme2 = norme(mt2) ;
    
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

Tbl diffrelmax(const Mtbl& mt1, const Mtbl& mt2) {
    
    // Protections
    assert(mt1.get_etat() != ETATNONDEF) ;
    assert(mt2.get_etat() != ETATNONDEF) ;
    
    int nz = mt1.get_nzone() ;
    Tbl resu(nz) ; 
    
    Tbl max2 = max(abs(mt2)) ;
    Tbl maxdiff = max(abs(mt1 - mt2)) ; 
    
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


