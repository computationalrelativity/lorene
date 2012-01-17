/*
 *  Mathematical functions for class Valeur
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


char valeur_math_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2012/01/17 10:39:27  j_penner
 * added a Heaviside function
 *
 * Revision 1.2  2002/10/16 14:37:15  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  1999/12/02  17:57:42  phil
 * *** empty log message ***
 *
 * Revision 2.3  1999/11/09  16:18:02  phil
 * ajout de racine_cubique
 *
 * Revision 2.2  1999/10/29  15:14:50  eric
 * Ajout de fonctions mathematiques (abs, norme, etc...).
 *
 * Revision 2.1  1999/03/01  15:11:03  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/24  15:40:32  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Fichiers include
// ----------------
#include <math.h>
#include <stdlib.h>

#include "valeur.h"


			    //-------//
			    // Sinus //
			    //-------//

Valeur sin(const Valeur& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    if (ti.c == 0x0) {			// Il faut la valeur physique
	ti.coef_i() ;
    }
    Valeur to(ti.get_mg()) ;		// Valeur resultat
    to.set_etat_c_qcq() ;
    *(to.c) = sin( *(ti.c) ) ;
    return to ;
}

			    //---------//
			    // Cosinus //
			    //---------//

Valeur cos(const Valeur& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    Valeur to(ti.get_mg()) ;		// Valeur resultat
    to.set_etat_c_qcq() ;

    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	*(to.c) = 1. ;
	return to ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    if (ti.c == 0x0) {			// Il faut la valeur physique
	ti.coef_i() ;
    }
    *(to.c) = cos( *(ti.c) ) ;
    return to ;
}

			    //----------//
			    // Tangente //
			    //----------//

Valeur tan(const Valeur& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    if (ti.c == 0x0) {			// Il faut la valeur physique
	ti.coef_i() ;
    }
    Valeur to(ti.get_mg()) ;		// Valeur resultat
    to.set_etat_c_qcq() ;
    *(to.c) = tan( *(ti.c) ) ;
    return to ;
}

			    //----------//
			    // ArcSinus //
			    //----------//

Valeur asin(const Valeur& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    if (ti.c == 0x0) {			// Il faut la valeur physique
	ti.coef_i() ;
    }
    Valeur to(ti.get_mg()) ;		// Valeur resultat
    to.set_etat_c_qcq() ;
    *(to.c) = asin( *(ti.c) ) ;
    return to ;
}

			    //------------//
			    // ArcCosinus //
			    //------------//

Valeur acos(const Valeur& ti)
{
   // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    Valeur to(ti.get_mg()) ;			// Valeur resultat
    to.set_etat_c_qcq() ;

    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	*(to.c) = M_PI * .5 ;
	return to ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    if (ti.c == 0x0) {			// Il faut la valeur physique
	ti.coef_i() ;
    }
    *(to.c) = acos( *(ti.c) ) ;
    return to ;
}

			    //-------------//
			    // ArcTangente //
			    //-------------//

Valeur atan(const Valeur& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    if (ti.c == 0x0) {			// Il faut la valeur physique
	ti.coef_i() ;
    }
    Valeur to(ti.get_mg()) ;		// Valeur resultat
    to.set_etat_c_qcq() ;
    *(to.c) = atan( *(ti.c) ) ;
    return to ;
}

			    //------//
			    // Sqrt //
			    //------//

Valeur sqrt(const Valeur& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    if (ti.c == 0x0) {			// Il faut la valeur physique
	ti.coef_i() ;
    }
    Valeur to(ti.get_mg()) ;		// Valeur resultat
    to.set_etat_c_qcq() ;
    *(to.c) = sqrt( *(ti.c) ) ;
    return to ;
}

			    //---------------//
			    // Exponantielle //
			    //---------------//

Valeur exp(const Valeur& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    Valeur to(ti.get_mg()) ;			// Valeur resultat
    to.set_etat_c_qcq() ;

    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	*(to.c) = 1. ;
	return to ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    if (ti.c == 0x0) {			// Il faut la valeur physique
	ti.coef_i() ;
    }
    *(to.c) = exp( *(ti.c) ) ;
    return to ;
}

			    //--------------------//
			    // Heaviside Function //
			    //--------------------//

Valeur Heaviside(const Valeur& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    Valeur to(ti.get_mg()) ;			// Valeur resultat
    to.set_etat_c_qcq() ;

    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	*(to.c) = 0. ;
	return to ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// otherwise
    if (ti.c == 0x0) {			// Use the physical value << What?? (used Google translate)
	ti.coef_i() ;
    }

    *(to.c) = Heaviside( *(ti.c) ) ;

    return to ;
}

			    //-------------//
			    // Log naturel //
			    //-------------//

Valeur log(const Valeur& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	cout << "Valeur log: log(ETATZERO) !" << endl  ;
	abort () ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    if (ti.c == 0x0) {			// Il faut la valeur physique
	ti.coef_i() ;
    }
    Valeur to(ti.get_mg()) ;		// Valeur resultat
    to.set_etat_c_qcq() ;
    *(to.c) = log( *(ti.c) ) ;
    return to ;
}

			    //-------------//
			    // Log decimal //
			    //-------------//

Valeur log10(const Valeur& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	cout << "Valeur log10: log10(ETATZERO) !" << endl ;
	abort () ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    if (ti.c == 0x0) {			// Il faut la valeur physique
	ti.coef_i() ;
    }
    Valeur to(ti.get_mg()) ;		// Valeur resultat
    to.set_etat_c_qcq() ;
    *(to.c) = log10( *(ti.c) ) ;
    return to ;
}

			    //--------------//
			    // Power entier //
			    //--------------//

Valeur pow(const Valeur& ti, int n)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	if (n > 0) {
	    return ti ;
	}
	else {
	    cout << "Valeur pow: ETATZERO^n with n<=0 ! "<< endl  ;
	    abort () ;
	}
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    if (ti.c == 0x0) {			// Il faut la valeur physique
	ti.coef_i() ;
    }
    Valeur to(ti.get_mg()) ;		// Valeur resultat
    to.set_etat_c_qcq() ;
    double x = n ;
    *(to.c) = pow( *(ti.c), x ) ;
    return to ;
}

			    //--------------//
			    // Power double //
			    //--------------//

Valeur pow(const Valeur& ti, double x)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	if (x > 0) {
	    return ti ;
	}
	else {
	    cout << "Valeur pow: ETATZERO^x with x<=0 !" << endl ;
	    abort () ;
	}
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...
    if (ti.c == 0x0) {			// Il faut la valeur physique
	ti.coef_i() ;
    }
    Valeur to(ti.get_mg()) ;			// Valeur resultat
    to.set_etat_c_qcq() ;
    *(to.c) = pow( *(ti.c), x ) ;
    return to ;
}

			    //----------------//
			    // Absolute value //
			    //----------------//

Valeur abs(const Valeur& vi)
{
    // Protection
    assert(vi.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (vi.get_etat() == ETATZERO) {
	return vi ;
    }
    
    // Cas general

    assert(vi.get_etat() == ETATQCQ) ;	// sinon...
    if (vi.c == 0x0) {			// Il faut la valeur physique
	vi.coef_i() ;
    }

    Valeur vo(vi.get_mg()) ;			// Valeur resultat

    vo.set_etat_c_qcq() ;

    *(vo.c) = abs( *(vi.c) ) ;		// abs(Mtbl)

    return vo ;
}
			     //----------------//
			    //    Cube root   //
			    //----------------//

Valeur racine_cubique(const Valeur& vi)
{   
// Protection
    assert(vi.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (vi.get_etat() == ETATZERO) {
	return vi ;
    }
    
    // Cas general
    assert(vi.get_etat() == ETATQCQ) ;	// sinon...
    if (vi.c == 0x0) {			// Il faut la valeur physique
	vi.coef_i() ;
    }
    Valeur vo(vi.get_mg()) ;		// Valeur resultat
    vo.set_etat_c_qcq() ;
    *(vo.c) = racine_cubique( *(vi.c) ) ;
    return vo ;
}
		    //-------------------------------//
    	    	    //            totalmax           //
		    //-------------------------------//

double totalmax(const Valeur& vi) {

    // Protection
    assert(vi.get_etat() != ETATNONDEF) ;
    
//    Tbl resu(vi.get_mg()->get_nzone()) ; 
    double resu ; 
    
    if (vi.get_etat() == ETATZERO) {
	resu = 0 ; 
    }
    else {

	assert(vi.get_etat() == ETATQCQ) ;	
	if (vi.c == 0x0) {			// Il faut la valeur physique
	    vi.coef_i() ;
	}

	resu = totalmax( *(vi.c) ) ;		// max(Mtbl) 
	
    }

    return resu ;
}

		    //-------------------------------//
    	    	    //           totalmin            //
		    //-------------------------------//

double totalmin(const Valeur& vi) {

    // Protection
    assert(vi.get_etat() != ETATNONDEF) ;
    
    double resu ; 

    if (vi.get_etat() == ETATZERO) {
	resu = 0 ; 
    }
    else {

	assert(vi.get_etat() == ETATQCQ) ;	
	if (vi.c == 0x0) {			// Il faut la valeur physique
	    vi.coef_i() ;
	}

	resu = totalmin( *(vi.c) ) ;		// min(Mtbl) 
	
    }

    return resu ; 
}

		    //-------------------------------//
    	    	    //            max                //
		    //-------------------------------//

Tbl max(const Valeur& vi) {

    // Protection
    assert(vi.get_etat() != ETATNONDEF) ;
    
    Tbl resu(vi.get_mg()->get_nzone()) ; 
    
    if (vi.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {

	assert(vi.get_etat() == ETATQCQ) ;	
	if (vi.c == 0x0) {			// Il faut la valeur physique
	    vi.coef_i() ;
	}

	resu = max( *(vi.c) ) ;		// max(Mtbl) 
	
    }
          
    return resu ; 
}

		    //-------------------------------//
    	    	    //            min                //
		    //-------------------------------//

Tbl min(const Valeur& vi) {

    // Protection
    assert(vi.get_etat() != ETATNONDEF) ;
    
    Tbl resu(vi.get_mg()->get_nzone()) ; 
    
    if (vi.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {

	assert(vi.get_etat() == ETATQCQ) ;	
	if (vi.c == 0x0) {			// Il faut la valeur physique
	    vi.coef_i() ;
	}

	resu = min( *(vi.c) ) ;		// min(Mtbl) 
	
    }
          
    return resu ; 
}

		    //-------------------------------//
    	    	    //            norme              //
		    //-------------------------------//

Tbl norme(const Valeur& vi) {

    // Protection
    assert(vi.get_etat() != ETATNONDEF) ;
    
    Tbl resu(vi.get_mg()->get_nzone()) ; 
    
    if (vi.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {

	assert(vi.get_etat() == ETATQCQ) ;	
	if (vi.c == 0x0) {			// Il faut la valeur physique
	    vi.coef_i() ;
	}

	resu = norme( *(vi.c) ) ;		// norme(Mtbl) 
	
    }
          
    return resu ; 
}

		    //-------------------------------//
    	    	    //          diffrel              //
		    //-------------------------------//

Tbl diffrel(const Valeur& v1, const Valeur& v2) {
    
    // Protections
    assert(v1.get_etat() != ETATNONDEF) ;
    assert(v2.get_etat() != ETATNONDEF) ;

    int nz = v1.get_mg()->get_nzone() ;
    Tbl resu(nz) ; 
    
    Valeur diff = v1 - v2 ; 
    Tbl normdiff = norme(diff) ; 
    Tbl norme2 = norme(v2) ;
    
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

Tbl diffrelmax(const Valeur& v1, const Valeur& v2) {
    
    // Protections
    assert(v1.get_etat() != ETATNONDEF) ;
    assert(v2.get_etat() != ETATNONDEF) ;
    
    int nz = v1.get_mg()->get_nzone() ;
    Tbl resu(nz) ; 
    
    Tbl max2 = max(abs(v2)) ;
    Valeur diff = v1 - v2 ;
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

