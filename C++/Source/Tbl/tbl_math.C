/*
 * Mathematical functions for class Tbl
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
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


char tbl_math_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2012/01/17 10:38:48  j_penner
 * added a Heaviside function
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  1999/12/02  17:47:48  phil
 * ajout de racine_cubique
 *
 * Revision 2.3  1999/10/29  15:05:32  eric
 * Ajout de fonctions mathematiques (abs, norme, etc...).
 *
 * Revision 2.2  1999/03/01  15:10:19  eric
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/24  15:26:13  hyc
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
#include "tbl.h"

			    //-------//
			    // Sinus //
			    //-------//

Tbl sin(const Tbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Tbl to(ti.dim) ;			// Tbl resultat
    to.set_etat_qcq() ;
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	to.t[i] = sin(ti.t[i]) ;
    }
    return to ;
}

			    //---------//
			    // Cosinus //
			    //---------//

Tbl cos(const Tbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;

    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	int taille = ti.get_taille() ;
	for (int i=0 ; i<taille ; i++) {
	    to.t[i] = 1 ;
	}
	return to ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	to.t[i] = cos(ti.t[i]) ;
    }
    return to ;
}

			    //----------//
			    // Tangente //
			    //----------//

Tbl tan(const Tbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	to.t[i] = tan(ti.t[i]) ;
    }
    return to ;
}

			    //----------//
			    // ArcSinus //
			    //----------//

Tbl asin(const Tbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	to.t[i] = asin(ti.t[i]) ;
    }
    return to ;
}

			    //------------//
			    // ArcCosinus //
			    //------------//

Tbl acos(const Tbl& ti)
{
   // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;

    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	int taille = ti.get_taille() ;
	for (int i=0 ; i<taille ; i++) {
	    to.t[i] = M_PI * .5 ;
	}
	return to ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	to.t[i] = acos(ti.t[i]) ;
    }
    return to ;
}

			    //-------------//
			    // ArcTangente //
			    //-------------//

Tbl atan(const Tbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	to.t[i] = atan(ti.t[i]) ;
    }
    return to ;
}

			    //------//
			    // Sqrt //
			    //------//

Tbl sqrt(const Tbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	to.t[i] = sqrt(ti.t[i]) ;
    }
    return to ;
}

			    //---------------//
			    // Exponantielle //
			    //---------------//

Tbl exp(const Tbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;

    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	int taille = ti.get_taille() ;
	for (int i=0 ; i<taille ; i++) {
	    to.t[i] = 1 ;
	}
	return to ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	to.t[i] = exp(ti.t[i]) ;
    }
    return to ;
}

			    //--------------------//
			    // Heaviside Function //
			    //--------------------//

Tbl Heaviside(const Tbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;

    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	int taille = ti.get_taille() ;
	for (int i=0 ; i<taille ; i++) {
	    to.t[i] = 0 ;
	}
	return to ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// Otherwise
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	if(ti.t[i] >= 0)
	to.t[i] = 1 ;
	else
	to.t[i] = 0 ;
    }
    return to ;
}

			    //-------------//
			    // Log naturel //
			    //-------------//

Tbl log(const Tbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	cout << "Tbl log: log(ETATZERO) !" << endl  ;
	abort () ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	to.t[i] = log(ti.t[i]) ;
    }
    return to ;
}

			    //-------------//
			    // Log decimal //
			    //-------------//

Tbl log10(const Tbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	cout << "Tbl log10: log10(ETATZERO) !" << endl ;
	abort () ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	to.t[i] = log10(ti.t[i]) ;
    }
    return to ;
}

			    //--------------//
			    // Power entier //
			    //--------------//

Tbl pow(const Tbl& ti, int n)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	if (n > 0) {
	    return ti ;
	}
	else {
	    cout << "Tbl pow: ETATZERO^n avec n<=0 ! "<< endl  ;
	    abort () ;
	}
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;
    double x = n ;
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	to.t[i] = pow(ti.t[i], x) ;		
    }
    return to ;
}

			    //--------------//
			    // Power double //
			    //--------------//

Tbl pow(const Tbl& ti, double x)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	if (x > 0) {
	    return ti ;
	}
	else {
	    cout << "Tbl pow: ETATZERO^x avec x<=0 !" << endl ;
	    abort () ;
	}
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;
    int taille = ti.get_taille() ;
    for (int i=0 ; i<taille ; i++) {
	to.t[i] = pow(ti.t[i], x) ;
    }
    return to ;
}

			    //----------------//
			    // Absolute value //
			    //----------------//

Tbl abs(const Tbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...

    Tbl to(ti.dim) ;		    // Tbl resultat
    to.set_etat_qcq() ;

    const double* xi = ti.t ; 
    double* xo = to.t ; 
    int taille = ti.get_taille() ;

    for (int i=0 ; i<taille ; i++) {
	xo[i] = fabs( xi[i] ) ;
    }
    
    return to ;
}
			    //----------------//
			    //	    Cubic     //
			    //----------------//

Tbl racine_cubique(const Tbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...

    Tbl absolute(abs(ti)) ;
    Tbl res (pow(absolute, 1./3.)) ;
    
    for (int i=0 ; i<ti.get_taille() ; i++)
	if (ti.t[i] < 0)
	    res.t[i] *= -1 ;
    
    return res ;
}

		    //-------------------------------//
    	    	    //            max                //
		    //-------------------------------//

double max(const Tbl& ti) {
    
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (ti.get_etat() == ETATZERO) {
    	return double(0) ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;     // sinon....
    
    const double* x = ti.t ; 
    double resu = x[0] ;
    for (int i=1; i<ti.get_taille(); i++) {
    	if ( x[i] > resu ) resu = x[i] ;
    }
	
    return resu ; 
}

		    //-------------------------------//
    	    	    //            min                //
		    //-------------------------------//

double min(const Tbl& ti) {
    
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (ti.get_etat() == ETATZERO) {
    	return double(0) ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;     // sinon....
    
    const double* x = ti.t ; 
    double resu = x[0] ;
    for (int i=1; i<ti.get_taille(); i++) {
    	if ( x[i] < resu ) resu = x[i] ;
    }
	
    return resu ; 
}

		    //-------------------------------//
    	    	    //            norme              //
		    //-------------------------------//

double norme(const Tbl& ti) {

    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    double resu = 0 ; 
    
    if (ti.get_etat() != ETATZERO) {  // on n'effectue la somme que si necessaire

	assert(ti.get_etat() == ETATQCQ) ;     // sinon....
	const double* x = ti.t ; 
	for (int i=0; i<ti.get_taille(); i++) {
	    resu += fabs( x[i] ) ;
	}

    }

    return resu ;
}

		    //-------------------------------//
    	    	    //          diffrel              //
		    //-------------------------------//

double diffrel(const Tbl& t1, const Tbl& t2) {
    
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    
    double norm2 = norme(t2) ;
    double normdiff = norme(t1-t2) ;
    double resu ;  
    if ( norm2 == double(0) ) {
	resu = normdiff ; 
    }
    else {
	resu =  normdiff / norm2 ; 
    }
    
    return resu ; 
    
}

		    //-------------------------------//
    	    	    //       diffrelmax              //
		    //-------------------------------//

double diffrelmax(const Tbl& t1, const Tbl& t2) {
    
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    
    double max2 = max(abs(t2)) ;
    double maxdiff = max(abs(t1-t2)) ;
    double resu ;  
    if ( max2 == double(0) ) {
	resu = maxdiff ; 
    }
    else {
	resu =  maxdiff / max2 ; 
    }
    
    return resu ; 
    
}
