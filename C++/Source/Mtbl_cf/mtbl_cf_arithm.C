/*
 * Arithmetical operations for class Mtbl_cf
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:17:59  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:08  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:15  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2002/10/16 14:36:43  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2000/09/27  14:25:15  eric
 * Multiplication par un double : on met le resultat a ETATZERO si x == 0.
 *
 * Revision 2.5  2000/08/16  10:43:18  eric
 * Suppression du membre dzpuis.
 *
 * Revision 2.4  1999/10/26  08:09:03  eric
 * Ajout de protection dzpuis dans +=, -=
 *
 * Revision 2.3  1999/10/18  15:07:34  eric
 * La fonction membre annule() est rebaptisee annule_hard().
 *
 * Revision 2.2  1999/10/13  15:50:49  eric
 * *** empty log message ***
 *
 * Revision 2.1  1999/10/01  14:50:14  eric
 * Ajout des operations manquantes.
 *
 * Revision 2.0  1999/06/23  12:36:43  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Fichiers include
// ----------------
#include <cmath>
#include <cassert>
#include <cstdlib>

#include "mtbl_cf.h"

			//********************//
			// OPERATEURS UNAIRES //
			//********************//

// + Mtbl_cf
// ---------
namespace Lorene {
Mtbl_cf operator+(const Mtbl_cf& t1)	    
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    
    return t1 ;
}


// - Mtbl_cf
// ---------
Mtbl_cf operator-(const Mtbl_cf& t1)	    
{

    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (t1.get_etat() == ETATZERO) {
	return t1 ;
    }
    
    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...
    Mtbl_cf r(t1) ;	// Mtbl_cf resultat

    for (int i=0 ; i<r.get_nzone() ; i++) {
	*(r.t)[i] = -(*(t1.t)[i]) ;
    }
    return r ;
}

			//**********//
			// ADDITION //
			//**********//

// Mtbl_cf + Mtbl_cf
// -----------------
Mtbl_cf operator+(const Mtbl_cf& t1, const Mtbl_cf& t2)	    
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mg() == t2.get_mg()) ;
    assert(t1.base == t2.base) ;
    
    // Cas particulier
    if (t1.get_etat() == ETATZERO) {
    	return t2 ;
    }
    if (t2.get_etat() == ETATZERO) {
    	return t1 ;
    }
    assert(t1.get_etat() == ETATQCQ) ;  // sinon...
    assert(t2.get_etat() == ETATQCQ) ;  // sinon...

    // Cas general
    int nz = t1.get_nzone() ;

    Mtbl_cf r(t1) ;	// Mtbl resultat

    for (int i=0 ; i<nz ; i++) {
	*(r.t)[i] += *(t2.t)[i] ;
    }

    return r ;
}



			//**************//
			// SOUSTRACTION //
			//**************//

// Mtbl_cf - Mtbl_cf
// -----------------
Mtbl_cf operator-(const Mtbl_cf& t1, const Mtbl_cf& t2)	 
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mg() == t2.get_mg()) ;
    assert(t1.base == t2.base) ;
    
    // Cas particulier
    if (t1.get_etat() == ETATZERO) {
    	return - t2 ;
    }
    if (t2.get_etat() == ETATZERO) {
    	return t1 ;
    }
    assert(t1.get_etat() == ETATQCQ) ;  // sinon...
    assert(t2.get_etat() == ETATQCQ) ;  // sinon...

    // Cas general
    int nz = t1.get_nzone() ;

    Mtbl_cf r(t1) ;	// Mtbl_cf resultat

    for (int i=0 ; i<nz ; i++) {
	*(r.t)[i] -= *(t2.t)[i] ;
    }

    return r ;
}


			//****************//
			// MULTIPLICATION //
			//****************//


// Mtbl_cf * double
// ----------------
Mtbl_cf operator*(const Mtbl_cf& t1, double x)	    // Mtbl_cf * double
{

    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if ((t1.get_etat() == ETATZERO) || ( x == double(1) )) {
	return t1 ;
    }

    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...

    Mtbl_cf r(t1) ;	// Mtbl_cf resultat

    if ( x == double(0) ) {
	r.set_etat_zero() ;
    }
    else{
	int nz = t1.get_nzone() ;
	for (int i=0 ; i<nz ; i++) {
	*(r.t)[i] *= x ;
	}
    }
    
    return r ;
}

// double * Mtbl_cf
// ----------------
Mtbl_cf operator*(double x, const Mtbl_cf& t1)	    
{
    return t1 * x ;
}

// Mtbl_cf * int
// -------------
Mtbl_cf operator*(const Mtbl_cf& t1, int m)
{
    return t1 * double(m) ;
}

// int * Mtbl_cf
// -------------
Mtbl_cf operator*(int m, const Mtbl_cf& t1)	    
{
    return t1 * double(m) ;
}


			//**********//
			// DIVISION //
			//**********//


// Mtbl_cf / double
// ----------------
Mtbl_cf operator/(const Mtbl_cf& t1, double x)	    
{

    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    if ( x == double(0) ) {
	cout << "Mtbl_cf division by 0 !" << endl ;
	abort() ;
    }
    
    // Cas particulier
    if ((t1.get_etat() == ETATZERO) || ( x == double(1) )) {
	return t1 ;
    }
    
    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...

    Mtbl_cf r(t1) ;	// Mtbl_cf resultat
    int nz = t1.get_nzone() ;
    for (int i=0 ; i<nz ; i++) {
	*(r.t)[i] /= x ;
    }

    return r ;
}

// Mtbl_cf / int
// -------------
Mtbl_cf operator/(const Mtbl_cf& t1, int n)	    
{
    return t1/double(n) ;
}



			//*******************//
			// operateurs +=,... //
			//*******************//

void Mtbl_cf::operator+=(const Mtbl_cf & mi) {
    
    // Protection
    assert(mg == mi.get_mg()) ;		    // meme grille
    assert(base == mi.base) ;		    // meme base
    assert(etat != ETATNONDEF) ;	    // etat defini
    assert(mi.get_etat() != ETATNONDEF) ;   // etat defini
    
    // Cas particulier
    if (mi.get_etat() == ETATZERO) {
	return ;
    }
    
    // Cas general

    if (etat == ETATZERO) {
	annule_hard() ;
    }
    for (int i=0 ; i<nzone ; i++) {
	*(t[i]) += *(mi.t[i]) ;
    }
}

void Mtbl_cf::operator-=(const Mtbl_cf & mi) {
    
    // Protection
    assert(mg == mi.get_mg()) ;		    // meme grille
    assert(base == mi.base) ;		    // meme base
    assert(etat != ETATNONDEF) ;	    // etat defini
    assert(mi.get_etat() != ETATNONDEF) ;   // etat defini
    
    // Cas particulier
    if (mi.get_etat() == ETATZERO) {
	return ;
    }
    
    // Cas general

    if (etat == ETATZERO) {
	annule_hard() ;
    }
    for (int i=0 ; i<nzone ; i++) {
	*(t[i]) -= *(mi.t[i]) ;
    }
}

}
