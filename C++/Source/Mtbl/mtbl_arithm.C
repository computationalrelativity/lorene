/*
 *  Arithmetical operations for class Mtbl
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


char mtbl_arithm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/13 08:53:08  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.8  2000/09/27  14:21:12  eric
 * Multiplication par un double : on met le resultat a ETATZERO si
 *   x == 0
 *
 * Revision 2.7  2000/08/16  10:30:14  eric
 * Suppression du membre dzpuis.
 *
 * Revision 2.6  1999/10/26  08:08:25  eric
 * Ajout de protection dzpuis dans +=, -=
 *
 * Revision 2.5  1999/10/18  15:07:22  eric
 * La fonction membre annule() est rebaptisee annule_hard().
 *
 * Revision 2.4  1999/10/15  13:58:04  eric
 * L'arithmetique liee aux Coord's se trouve desormais dans le fichier
 *   arithm_coord.C.
 *
 * Revision 2.3  1999/10/01  10:09:01  eric
 * Changement prototypes.
 * Correction erreurs lorsque etat=ETATZERO
 *
 * Revision 2.2  1999/04/26  17:12:26  phil
 * ajout de Coord * Mtbl
 *
 * Revision 2.1  1999/02/22  15:50:22  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/15  10:42:45  hyc
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/15  09:59:50  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1999/01/15  09:10:39  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "mtbl.h"
#include "coord.h"

			//********************//
			// OPERATEURS UNAIRES //
			//********************//

// + Mtbl
// ------
namespace Lorene {
Mtbl operator+(const Mtbl& t1)	    // + Mtbl
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    
    return t1 ;
}

// - Mtbl
// ------
Mtbl operator-(const Mtbl& t1)	    // - Mtbl
{

    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (t1.get_etat() == ETATZERO) {
	return t1 ;
    }
    
    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...
    Mtbl r(t1) ;	// Mtbl resultat

    for (int i=0 ; i<r.get_nzone() ; i++) {
	*(r.t)[i] = -(*(t1.t)[i]) ;
    }
    return r ;
}

			//**********//
			// ADDITION //
			//**********//

// Mtbl + Mtbl
// -----------
Mtbl operator+(const Mtbl& t1, const Mtbl& t2)	    // Mtbl + Mtbl
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mg() == t2.get_mg()) ;
    
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

    Mtbl r(t1) ;	// Mtbl resultat

    for (int i=0 ; i<nz ; i++) {
	*(r.t)[i] += *(t2.t)[i] ;
    }

    return r ;
}

// Mtbl + double
// -------------
Mtbl operator+(const Mtbl& t1, double x)	    // Mtbl + double
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if (x == double(0)) {
	return t1 ;
    }

    int nz = t1.get_nzone() ;


    Mtbl r(t1) ;	// Mtbl resultat

    if (r.get_etat() == ETATZERO) {
	r.set_etat_qcq() ;
	for (int i=0 ; i<nz ; i++) {
	    r.t[i]->set_etat_zero() ;
	    *(r.t)[i] += x ;
	}		
    }
    else{ 
	assert(r.get_etat() == ETATQCQ) ;

	for (int i=0 ; i<nz ; i++) {
	    *(r.t)[i] += x ;
	}
    }
    
    return r ;
}

// double + Mtbl
// -------------
Mtbl operator+(double x, const Mtbl& t1)	    // double + Mtbl
{
    return t1 + x ;
}

// Mtbl + int
// ----------
Mtbl operator+(const Mtbl& t1, int m)	    // Mtbl + int
{
    return t1 + double(m) ;
}

// int + Mtbl
// ----------
Mtbl operator+(int m, const Mtbl& t1)	    // int + Mtbl
{
    return t1 + double(m) ;
}


			//**************//
			// SOUSTRACTION //
			//**************//

// Mtbl - Mtbl
// -----------
Mtbl operator-(const Mtbl& t1, const Mtbl& t2)	    // Mtbl - Mtbl
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mg() == t2.get_mg()) ;
    
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

    Mtbl r(t1) ;	// Mtbl resultat

    for (int i=0 ; i<nz ; i++) {
	*(r.t)[i] -= *(t2.t)[i] ;
    }

    return r ;
}

// Mtbl - double
// -------------
Mtbl operator-(const Mtbl& t1, double x)	    // Mtbl - double
{

    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if (x == double(0)) {
	return t1 ;
    }
    
    // Cas general
    int nz = t1.get_nzone() ;

    Mtbl r(t1) ;	// Mtbl resultat

    if (r.get_etat() == ETATZERO) {
	r.set_etat_qcq() ;
	for (int i=0 ; i<nz ; i++) {
	    r.t[i]->set_etat_zero() ;
	    *(r.t)[i] -= x ;
	}		
    }
    else{ 
	assert(r.get_etat() == ETATQCQ) ;

	for (int i=0 ; i<nz ; i++) {
	    *(r.t)[i] -= x ;
	}
    }
    
    return r ;
}

// double - Mtbl
// -------------
Mtbl operator-(double x, const Mtbl& t1)	    // double - Mtbl
{
    return - (t1 -x) ;
}

// Mtbl - int
// ----------
Mtbl operator-(const Mtbl& t1, int m)	    // Mtbl - int
{
    return t1 - double(m) ;
}

// int - Mtbl
// ----------
Mtbl operator-(int m, const Mtbl& t1)	    // int - Mtbl
{
    return double(m) - t1 ;
}

			//****************//
			// MULTIPLICATION //
			//****************//

// Mtbl * Mtbl
// -----------
Mtbl operator*(const Mtbl& t1, const Mtbl& t2)	    // Mtbl * Mtbl
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mg() == t2.get_mg()) ;

    // Cas particulier
    if (t1.get_etat() == ETATZERO) {
	return t1 ;
    }
    if (t2.get_etat() == ETATZERO) {
	return t2 ;
    }
    
    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...
    assert(t2.get_etat() == ETATQCQ) ;	// sinon...

    Mtbl r(t1) ;	// Mtbl resultat

    int nz = t1.get_nzone() ;
    for (int i=0 ; i<nz ; i++) {
	*(r.t)[i] *= (*(t2.t)[i]) ;
    }

    return r ;
}

// Mtbl * double
// -------------
Mtbl operator*(const Mtbl& t1, double x)	    // Mtbl * double
{

    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if ((t1.get_etat() == ETATZERO) || ( x == double(1) )) {
	return t1 ;
    }

    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...

    Mtbl r(t1) ;	// Mtbl resultat

    if ( x == double(0) ) {
	r.set_etat_zero() ;
    }
    else {
	int nz = t1.get_nzone() ;
	for (int i=0 ; i<nz ; i++) {
	    *(r.t)[i] *= x ;
	}
    }
    
    return r ;
}

// double * Mtbl
// -------------
Mtbl operator*(double x, const Mtbl& t1)	    // double * Mtbl
{
    return t1 * x ;
}

// Mtbl * int
// ----------
Mtbl operator*(const Mtbl& t1, int m)	    // Mtbl * int
{
    return t1 * double(m) ;
}

// int * Mtbl
// ----------
Mtbl operator*(int m, const Mtbl& t1)	    // int * Mtbl
{
    return t1 * double(m) ;
}

			//**********//
			// DIVISION //
			//**********//

// Mtbl / Mtbl
// -----------
Mtbl operator/(const Mtbl& t1, const Mtbl& t2)	    // Mtbl / Mtbl
{

    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mg() == t2.get_mg()) ;

    // Cas particuliers
    if (t2.get_etat() == ETATZERO) {
	cout << "Mtbl division by 0 !" << endl ;
	abort() ; 
    }
    if (t1.get_etat() == ETATZERO) {
	return t1 ;
    }
    
    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...
    assert(t2.get_etat() == ETATQCQ) ;	// sinon...
    
    Mtbl r(t1) ;	// Mtbl resultat

    int nz = t1.get_nzone() ;
    for (int i=0 ; i<nz ; i++) {
	*(r.t)[i] /= (*(t2.t)[i]) ;
    }

    return r ;
}

// Mtbl / double
// -------------
Mtbl operator/(const Mtbl& t1, double x)	    // Mtbl / double
{

    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    if ( x == double(0) ) {
	cout << "Mtbl division by 0 !" << endl ;
	abort() ;
    }
    
    // Cas particulier
    if ((t1.get_etat() == ETATZERO) || ( x == double(1) )) {
	return t1 ;
    }
    
    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...

    Mtbl r(t1) ;	// Mtbl resultat
    int nz = t1.get_nzone() ;
    for (int i=0 ; i<nz ; i++) {
	*(r.t)[i] /= x ;
    }

    return r ;
}

// Mtbl / int
// ----------
Mtbl operator/(const Mtbl& t1, int n)	    // Mtbl / int
{
    return t1/double(n) ;
}

// double / Mtbl
// -------------
Mtbl operator/(double x, const Mtbl& t1)	    // double / Mtbl
{

    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particuliers
    if (t1.get_etat() == ETATZERO) {
	cout << "Division by 0 in double / Mtbl !" << endl ;
	abort() ; 
    }

    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...

    Mtbl r( *(t1.get_mg()) ) ;	// Mtbl resultat, a priori NONDEF

    if ( x == double(0) ) {
	r.set_etat_zero() ;
    }
    else {
    	r.set_etat_qcq() ;
    	int nz = t1.get_nzone() ;
    	for (int i=0 ; i<nz ; i++) {
	    *(r.t)[i] = x / (*(t1.t)[i]) ;
    	}
    }
    
    // Termine
    return r ;
}

// int / Mtbl
// ----------
Mtbl operator/(int m, const Mtbl& t1)	    // int / Mtbl
{
    return double(m)/t1 ;
}


			//*******************//
			// operateurs +=,... //
			//*******************//

void Mtbl::operator+=(const Mtbl & mi) {
    
    // Protection
    assert(mg == mi.get_mg()) ;		    // meme grille
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

void Mtbl::operator-=(const Mtbl & mi) {
    
    // Protection
    assert(mg == mi.get_mg()) ;		    // meme grille
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

void Mtbl::operator*=(const Mtbl & mi) {
    
    // Protection
    assert(mg == mi.get_mg()) ;		    // meme grille
    assert(etat != ETATNONDEF) ;	    // etat defini
    assert(mi.get_etat() != ETATNONDEF) ;   // etat defini
    
    // Cas particuliers
    if (etat == ETATZERO) {
	return ;
    }
    if (mi.get_etat() == ETATZERO) {
	set_etat_zero() ;
	return ;
    }
    
    // Cas general
    for (int i=0 ; i<nzone ; i++) {
	*(t[i]) *= *(mi.t[i]) ;
    }

}
}
