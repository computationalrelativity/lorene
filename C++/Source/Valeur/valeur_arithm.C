/*
 *  Arithmetical operations for class Valeur
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2005 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 2004 Jerome Novak
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
 * Revision 1.8  2016/12/05 16:18:20  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2014/10/13 08:53:49  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:13:22  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2008/08/27 08:52:55  jl_cornou
 * Added Jacobi(0,2) polynomials case
 *
 * Revision 1.4  2005/11/17 15:19:23  e_gourgoulhon
 * Added Valeur + Mtbl and Valeur - Mtbl.
 *
 * Revision 1.3  2004/07/06 13:36:30  j_novak
 * Added methods for desaliased product (operator |) only in r direction.
 *
 * Revision 1.2  2002/10/16 14:37:15  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.20  2001/08/31  15:23:48  eri * operator% : traitement du cas ou le Tbl est zero dans une zone.
 *
 * Revision 2.19  2001/05/28  12:42:27  eric
 * Passage en ylm pour le desaliasing dans operator%.
 *
 * Revision 2.18  2001/05/26  14:50:21  eric
 * Ajout de l'operator% : produit de deux Valeur avec desaliasage.
 *
 * Revision 2.17  2000/01/14  14:42:47  eric
 * Valeur * double : operation effectue dans l'espace des coefficients si
 *                   la Valeur n'est pas a jour ds l'espace des config.
 * Valeur / double : idem.
 * += Valeur : idem.
 * -= Valeur : idem.
 *
 * Pour += Valeur et -= Valeur les schemas sont desormais calques sur
 * Valeur + Valeur et Valeur - Valeur.
 *
 * Revision 2.16  1999/12/10  16:33:36  eric
 * Dans l'arithmetique membre (+=, -=, *=), on n'appelle desormais
 *  del_deriv() que tout a la fin.
 *
 * Revision 2.15  1999/11/30  14:12:54  phil
 * gestion de base dans operator/ (double,Vlaeur)
 *
 * Revision 2.14  1999/11/30  12:42:10  eric
 * Le membre base est desormais un objet de type Base_val et non plus
 *  un pointeur vers une Base_val.
 *
 * Revision 2.13  1999/11/29  13:28:06  eric
 * *** empty log message ***
 *
 * Revision 2.12  1999/11/29  10:20:50  eric
 * Ajout de Valeur/Mtbl et Mtbl / Valeur.
 *
 * Revision 2.11  1999/11/29  10:06:05  eric
 * Ajout de Valeur*Mtbl et Mtbl*Valeur
 *
 * Revision 2.10  1999/10/26  14:40:29  phil
 * On gere les bases pour *, *=, /= et /
 *
 * Revision 2.9  1999/10/20  15:31:28  eric
 * Ajout de la plupart des fonctions arithmetiques.
 *
 * Revision 2.8  1999/10/18  15:07:47  eric
 * La fonction membre annule() est rebaptisee annule_hard().
 *
 * Revision 2.7  1999/10/13  15:50:56  eric
 * Depoussierage.
 *
 * Revision 2.6  1999/09/15  10:02:26  phil
 * gestion de la base dans Valeur operator (double, const Valeur &)
 *
 * Revision 2.5  1999/09/14  17:18:36  phil
 * aout de Valeur operator* (double, const Valeur&)
 *
 * Revision 2.4  1999/09/13  14:52:55  phil
 * *** empty log message ***
 *
 * Revision 2.3  1999/04/09  12:38:05  phil
 * *** empty log message ***
 *
 * Revision 2.2  1999/04/09  12:26:59  phil
 * Ajout de valeur * coord
 *
 * Revision 2.1  1999/02/22  15:49:28  hyc
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

#include "mtbl.h"
#include "valeur.h"
#include "coord.h"
			//********************//
			// OPERATEURS UNAIRES //
			//********************//

namespace Lorene {
Valeur operator+(const Valeur & vi) {
    
    // Protection
    assert(vi.get_etat() != ETATNONDEF) ;
    
    return vi ;
}

Valeur operator-(const Valeur & vi) {
    
    // Protection
    assert(vi.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (vi.get_etat() == ETATZERO) {
	return vi ;
    }
    
    // Cas general
    assert(vi.get_etat() == ETATQCQ) ;	// sinon...

    Valeur resu(vi.get_mg()) ;	// Valeur resultat

    if (vi.c != 0x0) {
	resu = - *(vi.c) ; 
	resu.base = vi.base ;  // N'oublions pas la base...
    }
    else{
	assert(vi.c_cf != 0x0) ;
	resu = - *(vi.c_cf) ;  // Dans ce cas la base est prise en 
			           // charge par l'operator=(const Mtbl_cf&).
    }    
    
    // Termine
    return resu ;
}

			//**********//
			// ADDITION //
			//**********//

// Valeur + Valeur
// ---------------
Valeur operator+(const Valeur& t1, const Valeur& t2)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mg() == t2.get_mg()) ;
    
    // Cas particuliers
    if (t1.get_etat() == ETATZERO) {
    	return t2 ;
    }
    if (t2.get_etat() == ETATZERO) {
    	return t1 ;
    }
    assert(t1.get_etat() == ETATQCQ) ;  // sinon...
    assert(t2.get_etat() == ETATQCQ) ;  // sinon...

    Valeur resu(t1.get_mg()) ; 
    
    // On privelegie l'addition dans l'espace des configurations: 
    // ----------------------------------------------------------
    if (t1.c != 0x0) {
	if (t2.c != 0x0) {
	    resu = *(t1.c) + *(t2.c) ;	    // Addition des Mtbl
	    resu.base = t1.base ;  
	}
	else {
	    assert(t2.c_cf != 0x0) ;
	    if (t1.c_cf != 0x0) {
		resu = *(t1.c_cf) + *(t2.c_cf) ;    // Addition des Mtbl_cf
	    }
	    else {
		t2.coef_i() ; 
		resu = *(t1.c) + *(t2.c) ;	    // Addition des Mtbl
		resu.base = t1.base ;  
	    }
	}
    }
    else{	// Cas ou t1.c n'est pas a jour
	assert(t1.c_cf != 0x0) ;
	if (t2.c_cf != 0x0) {
	    resu = *(t1.c_cf) + *(t2.c_cf) ;    // Addition des Mtbl_cf
	}
	else {
	    assert(t2.c != 0x0) ; 
	    t1.coef_i() ; 
	    resu = *(t1.c) + *(t2.c) ;	    // Addition des Mtbl
	    resu.base = t1.base ;  
	}	    
    }

    return resu ;
}

// Valeur + Mtbl
// -------------
Valeur operator+(const Valeur& t1, const Mtbl& mi)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(mi.get_etat() != ETATNONDEF) ;
    
    // Cas particuliers
    if (mi.get_etat() == ETATZERO) {
	return t1 ;
    }

    Valeur resu(t1) ;
    
    if (t1.get_etat() == ETATZERO) {
	resu.set_etat_c_qcq() ;
	*(resu.c) = mi ; 
    }
    else{ 
	assert(resu.get_etat() == ETATQCQ) ; // sinon ...
	resu.coef_i() ;	// l'addition se fait dans l'espace des configurations
	resu.set_etat_c_qcq() ;
	*(resu.c) += mi ;
    }
        
    return resu ;
}

// Mtbl + Valeur 
// -------------
Valeur operator+(const Mtbl& mi, const Valeur& t1) {
    return t1 + mi ; 
}


// Valeur + double
// ---------------
Valeur operator+(const Valeur& t1, double x)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    
    // Cas particuliers
    if (x == double(0)) {
	return t1 ;
    }

    Valeur resu(t1) ;
    
    if (t1.get_etat() == ETATZERO) {
	resu.set_etat_c_qcq() ;
	*(resu.c) = x ; 
    }
    else{ 
	assert(resu.get_etat() == ETATQCQ) ; // sinon ...
	resu.coef_i() ;	// l'addition se fait dans l'espace des configurations
	resu.set_etat_c_qcq() ;
	*(resu.c) = *(resu.c) + x ;
    }
        
    return resu ;
}

// double + Valeur
// ---------------
Valeur operator+(double x, const Valeur& t1)	    // double + Valeur
{
    return t1 + x ;
}

// Valeur + int
// ------------
Valeur operator+(const Valeur& t1, int m)	    // Valeur + int
{
    return t1 + double(m) ;
}

// int + Valeur
// -------------
Valeur operator+(int m, const Valeur& t1)	    // int + Valeur
{
    return t1 + double(m) ;
}



			//**************//
			// SOUSTRACTION //
			//**************//

// Valeur - Valeur
// ---------------
Valeur operator-(const Valeur& t1, const Valeur& t2)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mg() == t2.get_mg()) ;
    
    // Cas particuliers
    if (t1.get_etat() == ETATZERO) {
    	return - t2 ;
    }
    if (t2.get_etat() == ETATZERO) {
    	return t1 ;
    }
    assert(t1.get_etat() == ETATQCQ) ;  // sinon...
    assert(t2.get_etat() == ETATQCQ) ;  // sinon...

    Valeur resu(t1.get_mg()) ; 
    
    // On privelegie la soustraction dans l'espace des configurations: 
    // ---------------------------------------------------------------
    if (t1.c != 0x0) {
	if (t2.c != 0x0) {
	    resu = *(t1.c) - *(t2.c) ;	    // Soustraction des Mtbl
	    resu.base = t1.base ;  
	}
	else {
	    assert(t2.c_cf != 0x0) ;
	    if (t1.c_cf != 0x0) {
		resu = *(t1.c_cf) - *(t2.c_cf) ;    // Soustraction des Mtbl_cf
	    }
	    else {
		t2.coef_i() ; 
		resu = *(t1.c) - *(t2.c) ;	    // Soustraction des Mtbl
		resu.base = t1.base ;  
	    }
	}
    }
    else{	// Cas ou t1.c n'est pas a jour
	assert(t1.c_cf != 0x0) ;
	if (t2.c_cf != 0x0) {
	    resu = *(t1.c_cf) - *(t2.c_cf) ;    // Soustraction des Mtbl_cf
	}
	else {
	    assert(t2.c != 0x0) ; 
	    t1.coef_i() ; 
	    resu = *(t1.c) - *(t2.c) ;	    // Soustraction des Mtbl
	    resu.base = t1.base ;  
	}	    
    }

    return resu ;
}


// Valeur - Mtbl
// -------------
Valeur operator-(const Valeur& t1, const Mtbl& mi)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(mi.get_etat() != ETATNONDEF) ;
    
    // Cas particuliers
    if (mi.get_etat() == ETATZERO) {
	return t1 ;
    }

    Valeur resu(t1) ;
    
    if (t1.get_etat() == ETATZERO) {
	resu.set_etat_c_qcq() ;
	*(resu.c) = - mi ; 
    }
    else{ 
	assert(resu.get_etat() == ETATQCQ) ; // sinon ...
	resu.coef_i() ;	// substraction in configuration space
	resu.set_etat_c_qcq() ;
	*(resu.c) -= mi ;
    }
        
    return resu ;
}

// Mtbl - Valeur
// -------------
Valeur operator-(const Mtbl& mi, const Valeur& t1) {
    return - (t1 - mi) ; 
}

// Valeur - double
// ---------------
Valeur operator-(const Valeur& t1, double x)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    
    // Cas particuliers
    if (x == double(0)) {
	return t1 ;
    }

    Valeur resu(t1) ;
    
    if (t1.get_etat() == ETATZERO) {
	resu.set_etat_c_qcq() ;
	*(resu.c) = - x ; 
    }
    else{ 
	assert(resu.get_etat() == ETATQCQ) ; // sinon ...
	resu.coef_i() ;	// l'addition se fait dans l'espace des configurations
	resu.set_etat_c_qcq() ;
	*(resu.c) = *(resu.c) - x ;
    }
        
    return resu ;
}

// double - Valeur
// ---------------
Valeur operator-(double x, const Valeur& t1)	    // double - Valeur
{
    return - (t1 - x) ;
}

// Valeur - int
// ------------
Valeur operator-(const Valeur& t1, int m)	    // Valeur - int
{
    return t1 - double(m) ;
}

// int - Valeur
// -------------
Valeur operator-(int m, const Valeur& t1)	    // int - Valeur
{
    return double(m) - t1 ;
}





			//****************//
			// MULTIPLICATION //
			//****************//
			
// Valeur * Valeur
// ---------------
Valeur operator*(const Valeur& t1, const Valeur& t2)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mg() == t2.get_mg()) ;
    
    // Cas particuliers
    if (t1.get_etat() == ETATZERO) {
    	return t1 ;
    }
    if (t2.get_etat() == ETATZERO) {
    	return t2 ;
    }
    assert(t1.get_etat() == ETATQCQ) ;  // sinon...
    assert(t2.get_etat() == ETATQCQ) ;  // sinon...

    Valeur resu(t1.get_mg()) ; 

    // La multiplication est faite dans l'espace des configurations:    
    if (t1.c == 0x0) {
	t1.coef_i() ; 
    }
    if (t2.c == 0x0) {
	t2.coef_i() ; 
    }
    
    resu = (*(t1.c)) * (*(t2.c)) ;	// Multiplication des Mtbl 
    
    // affectation de la base :
    resu.base = t1.base * t2.base;
    
    return resu ;
}


// Valeur * double
// ---------------
Valeur operator*(const Valeur& c1, double a) {
    
    // Protection
    assert(c1.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if ((c1.get_etat() == ETATZERO) || ( a == double(1) )) {
	return c1 ;
    }

    // Cas general
    assert(c1.get_etat() == ETATQCQ) ;	// sinon...

    Valeur result(c1.get_mg()) ;

    if (c1.c != 0x0) {
	result = *(c1.c) * a ;	   // Mtbl * double
	result.base = c1.base ;    // in this case, result.base must be set
				   // by hand
    }
    else {
	assert(c1.c_cf != 0x0) ; 
	result = *(c1.c_cf) * a ;   // Mtbl_cf * double
    }

    return result ;
}

// double * Valeur
// ---------------
Valeur operator*(double x, const Valeur& t1)	    // double * Valeur
{
    return t1 * x ;
}

// Valeur * int
// ------------
Valeur operator*(const Valeur& t1, int m)	    // Valeur * int
{
    return t1 * double(m) ;
}

// int * Valeur
// ------------
Valeur operator*(int m, const Valeur& t1)	    // int * Valeur
{
    return t1 * double(m) ;
}


// Valeur * Mtbl
// --------------
Valeur operator*(const Valeur & c1, const Mtbl& c2) {
    
    // Protection
    assert(c1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if (c1.get_etat() == ETATZERO) {
	return c1 ;
    }

    // Cas general

    assert(c1.get_etat() == ETATQCQ) ;

    Valeur result(c1.get_mg()) ;

    // La multiplication  est faite dans l'espace des configurations:    
    if (c1.c == 0x0) {
	c1.coef_i() ; 
    }
    result = *(c1.c) * c2 ;	// Mtbl * Mtbl 

    return result ;
}

// Mtbl * Valeur
// --------------
Valeur operator*(const Mtbl& c1, const Valeur& t1) {
    
    return t1 * c1 ; 
    
}


// Valeur * Coord
// --------------
Valeur operator*(const Valeur & c1, const Coord& c2) {
    
    // Protection
    assert(c1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if (c1.get_etat() == ETATZERO) {
	return c1 ;
    }

    // Cas general

    assert(c1.get_etat() == ETATQCQ) ;

    Valeur result(c1.get_mg()) ;

    // La multiplication  est faite dans l'espace des configurations:    
    if (c1.c == 0x0) {
	c1.coef_i() ; 
    }
    result = *(c1.c) * c2 ;	// Mtbl * Coord

    return result ;
}

// Coord * Valeur
// --------------
Valeur operator*(const Coord& c1, const Valeur& t1) {
    
    return t1 * c1 ; 
    
}

			//**********//
			// DIVISION //
			//**********//

// Valeur / Valeur
// ---------------
Valeur operator/(const Valeur& t1, const Valeur& t2)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mg() == t2.get_mg()) ;
    
    // Cas particuliers
    if (t2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Valeur / Valeur !" << endl ;
	abort() ; 
    }
    if (t1.get_etat() == ETATZERO) {
    	return t1 ;
    }

    assert(t1.get_etat() == ETATQCQ) ;  // sinon...
    assert(t2.get_etat() == ETATQCQ) ;  // sinon...

    Valeur resu(t1.get_mg()) ; 

    // La division est faite dans l'espace des configurations:    
    if (t1.c == 0x0) {
	t1.coef_i() ; 
    }
    if (t2.c == 0x0) {
	t2.coef_i() ; 
    }
    
    resu = (*(t1.c)) / (*(t2.c)) ;	// Division des Mtbl 
    
    // affectation de la base :
    resu.base = t1.base * t2.base;
    
    return resu ;
}

// Valeur / double 
// ---------------
Valeur operator/(const Valeur& t1, double x)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    
    // Cas particuliers
    if ( x == double(0) ) {
	cout << "Division by 0 in Valeur / double !" << endl ;
	abort() ;
    }
    if ((t1.get_etat() == ETATZERO) || ( x == double(1) )) {
	return t1 ;
    }

    assert(t1.get_etat() == ETATQCQ) ;  // sinon...

    Valeur resu(t1.get_mg()) ; 

    if (t1.c != 0x0) {
	resu = *(t1.c) / x ;	 // Mtbl / double
	resu.base = t1.base ;    // in this case, resu.base must be set by hand
    }
    else {
	assert(t1.c_cf != 0x0) ; 
	resu = *(t1.c_cf) / x ;   // Mtbl_cf * double
    }
    
    return resu ;
}

// double / Valeur
// ---------------
Valeur operator/(double x, const Valeur & c2) {
    
    // Protection
    assert(c2.get_etat() != ETATNONDEF) ;

    // Cas particuliers
    if (c2.get_etat() == ETATZERO) {
	cout << "Division by 0 in double / Valeur !" << endl ;
	abort() ; 
    }
    
    // Cas general
    assert(c2.get_etat() == ETATQCQ) ;	// sinon...

    // Il faut les valeurs physiques de c2
    if (c2.c == 0x0) {
	c2.coef_i() ;
    }
    
    // Le resultat
    Valeur r(c2.get_mg()) ;
    r.set_etat_c_qcq() ;
    *(r.c) = x / *(c2.c) ;
    
    // affectation de la base :
    r.base = c2.get_mg()->std_base_scal() * c2.base ; 
    
    // Termine
    return r ;
}

// Valeur / int
// ------------
Valeur operator/(const Valeur& t1,  int m) {
    return t1 / double(m) ;
}

// int / Valeur
// ------------
Valeur operator/(int m, const Valeur& t1) {
    return double(m) / t1 ;
}

// Valeur / Mtbl 
// ---------------
Valeur operator/(const Valeur& t1, const Mtbl& m2)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(m2.get_etat() != ETATNONDEF) ;
    
    // Cas particuliers
    if ( m2.get_etat() == ETATZERO ) {
	cout << "Division by 0 in Valeur / Mtbl !" << endl ;
	abort() ;
    }
    if (t1.get_etat() == ETATZERO) {
	return t1 ;
    }

    assert(t1.get_etat() == ETATQCQ) ;  // sinon...

    Valeur resu(t1.get_mg()) ; 

    // La division est faite dans l'espace des configurations:    
    if (t1.c == 0x0) {
	t1.coef_i() ; 
    }
    
    resu = (*(t1.c)) / m2 ;	// Division Mtbl / Mtbl 
        
    return resu ;
}

// Mtbl / Valeur
// ---------------
Valeur operator/(const Mtbl& m1, const Valeur & c2) {
    
    // Protection
    assert(m1.get_etat() != ETATNONDEF) ;
    assert(c2.get_etat() != ETATNONDEF) ;

    // Cas particuliers
    if (c2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Mtbl / Valeur !" << endl ;
	abort() ; 
    }
    
    // Cas general
    assert(c2.get_etat() == ETATQCQ) ;	// sinon...

    // Le resultat
    Valeur resu(c2.get_mg()) ;

    // Il faut les valeurs physiques de c2
    if (c2.c == 0x0) {
	c2.coef_i() ;
    }
    
    resu = m1 / (*(c2.c))  ;	// Division Mtbl / Mtbl
    
    // Termine
    return resu ;
}






			//*******************//
			// operateurs +=,... //
			//*******************//

void Valeur::operator+=(const Valeur & vi) {
    
    // Protection
    assert(mg == vi.get_mg()) ;		    // meme grille
    assert(etat != ETATNONDEF) ;	    // etat defini
    assert(vi.get_etat() != ETATNONDEF) ;   // etat defini
    
    // Cas particulier
    if (vi.get_etat() == ETATZERO) {
	return ;
    }
    
    // Cas general
    
    // Cas de l'etat ZERO
    if (etat == ETATZERO) {
	annule_hard() ;
    }
    

    if (c != 0x0) {
	if (vi.c != 0x0) {
	    *c += *(vi.c) ;	    // += Mtbl
	    delete c_cf ; 
	    c_cf = 0x0 ; 
	}
	else {
	    assert(vi.c_cf != 0x0) ;
	    if (c_cf != 0x0) {
		*c_cf += *(vi.c_cf) ;    // += Mtbl_cf
		delete c ; 
		c = 0x0 ; 
	    }
	    else {
		vi.coef_i() ; 
		*c += *(vi.c) ;		 // += Mtbl
		delete c_cf ; 
		c_cf = 0x0 ; 
	    }
	}
    }
    else{	// Case where c is not up to date
	assert(c_cf != 0x0) ;
	if (vi.c_cf != 0x0) {
	    *c_cf += *(vi.c_cf) ;    // += Mtbl_cf
	    delete c ; 
	    c = 0x0 ; 
	}
	else {
	    assert(vi.c != 0x0) ; 
	    coef_i() ; 
	    *c += *(vi.c) ;		 // += Mtbl
	    delete c_cf ; 
	    c_cf = 0x0 ; 
	}	    
    }
    
    // Menage (a ne faire qu'a la fin seulement)
    del_deriv() ;

    
}

void Valeur::operator-=(const Valeur & vi) {
    
    // Protection
    assert(mg == vi.get_mg()) ;		    // meme grille
    assert(etat != ETATNONDEF) ;	    // etat defini
    assert(vi.get_etat() != ETATNONDEF) ;   // etat defini
    
    // Cas particulier
    if (vi.get_etat() == ETATZERO) {
	return ;
    }
    
    // Cas general
    
    // Cas de l'etat ZERO
    if (etat == ETATZERO) {
	annule_hard() ;
    }
    
    if (c != 0x0) {
	if (vi.c != 0x0) {
	    *c -= *(vi.c) ;	    // -= Mtbl
	    delete c_cf ; 
	    c_cf = 0x0 ; 
	}
	else {
	    assert(vi.c_cf != 0x0) ;
	    if (c_cf != 0x0) {
		*c_cf -= *(vi.c_cf) ;    // -= Mtbl_cf
		delete c ; 
		c = 0x0 ; 
	    }
	    else {
		vi.coef_i() ; 
		*c -= *(vi.c) ;		 // -= Mtbl
		delete c_cf ; 
		c_cf = 0x0 ; 
	    }
	}
    }
    else{			    // Case where c is not up to date
	assert(c_cf != 0x0) ;
	if (vi.c_cf != 0x0) {
	    *c_cf -= *(vi.c_cf) ;    // -= Mtbl_cf
	    delete c ; 
	    c = 0x0 ; 
	}
	else {
	    assert(vi.c != 0x0) ; 
	    coef_i() ; 
	    *c -= *(vi.c) ;		 // -= Mtbl
	    delete c_cf ; 
	    c_cf = 0x0 ; 
	}	    
    }

    // Menage (a ne faire qu'a la fin seulement)
    del_deriv() ;

}

void Valeur::operator*=(const Valeur & vi) {
    
    // Protection
    assert(mg == vi.get_mg()) ;		    // meme grille
    assert(etat != ETATNONDEF) ;	    // etat defini
    assert(vi.get_etat() != ETATNONDEF) ;   // etat defini
    
    // Cas particuliers
    if (etat == ETATZERO) {
	return ;
    }
    if (vi.get_etat() == ETATZERO) {
	set_etat_zero() ;
	return ;
    }
    
    // Cas general
    
    // Calcul dans l'espace physique
    if (c == 0x0) {
	coef_i() ;
    }
    if (vi.c == 0x0) {
	vi.coef_i() ;
    }
    
    // Calcul
    *c *= *(vi.c) ;
    
    // Affectation de la base :
    base = base * vi.base ;
    
    // Coefficients
    delete c_cf ;
    c_cf = 0x0 ;

    // Menage (a ne faire qu'a la fin seulement)
    del_deriv() ;

}

		    //-----------------------------------//
		    //	Multiplication without aliasing  //
		    //-----------------------------------//

Valeur operator%(const Valeur& t1, const Valeur& t2)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mg() == t2.get_mg()) ;
    
    // Cas particuliers
    if (t1.get_etat() == ETATZERO) {
    	return t1 ;
    }
    if (t2.get_etat() == ETATZERO) {
    	return t2 ;
    }
    assert(t1.get_etat() == ETATQCQ) ;  // sinon...
    assert(t2.get_etat() == ETATQCQ) ;  // sinon...

    // Grid
    const Mg3d& mg = *(t1.get_mg()) ; 

    // Grid with twice the number of points in each dimension:
    const Mg3d& mg2 = *(mg.get_twice()) ;  
    
    // The coefficients are required
    if (t1.c_cf == 0x0) {
	t1.coef() ; 
    }
    if (t2.c_cf == 0x0) {
	t2.coef() ; 
    }

    const Mtbl_cf& c1 = *(t1.c_cf) ; 
    const Mtbl_cf& c2 = *(t2.c_cf) ; 

    assert( c1.get_etat() == ETATQCQ ) ; 
    assert( c2.get_etat() == ETATQCQ ) ; 

    // The number of coefficients is multiplied by 2 and the additionnal
    //  coefficients are set to zero
    // -----------------------------------------------------------------
    
    Mtbl_cf cc1( mg2, c1.base ) ;
    Mtbl_cf cc2( mg2, c2.base ) ; 
    
    cc1.set_etat_qcq() ; 
    cc2.set_etat_qcq() ; 
    
    for (int l=0; l<mg.get_nzone(); l++) {

	int nr = mg.get_nr(l) ; 
	int nt = mg.get_nt(l) ; 
	int np = mg.get_np(l) ; 
	int nr2 = mg2.get_nr(l) ; 
	int nt2 = mg2.get_nt(l) ; 
	int np2 = mg2.get_np(l) ; 

	// Copy of the coefficients of t1
	// ------------------------------

	if ( c1.t[l]->get_etat() == ETATZERO ) {
	    cc1.t[l]->set_etat_zero() ; 
	}
	else {

	    assert( c1.t[l]->get_etat() == ETATQCQ ) ; 
	    cc1.t[l]->set_etat_qcq() ; 

	    // Copy of the coefficients of t1 
	    for (int k=0; k<np+1; k++) {
		for (int j=0; j<nt; j++) {
		    for (int i=0; i<nr; i++) {
			cc1.t[l]->set(k, j, i) = (*(c1.t[l]))(k, j, i) ; 
		    }
		}
	    }

	    // The extra phi coefficients are set to zero
	    for (int k=np+1; k<np2+2; k++) {
		for (int j=0; j<nt2; j++) {
		    for (int i=0; i<nr2; i++) {
			cc1.t[l]->set(k, j, i) = 0 ; 
		    }
		}
	    }

	    // The extra theta coefficients are set to zero
	    for (int k=0; k<np+1; k++) {
		for (int j=nt; j<nt2; j++) {
		    for (int i=0; i<nr2; i++) {
			cc1.t[l]->set(k, j, i) = 0 ; 
		    }
		}
	    }

	    // The extra r coefficients are set to zero
	    for (int k=0; k<np+1; k++) {
		for (int j=0; j<nt; j++) {
		    for (int i=nr; i<nr2; i++) {
			cc1.t[l]->set(k, j, i) = 0 ; 
		    }
		}
	    }
	    
	}

	// Copy of the coefficients of t2
	// ------------------------------

	if ( c2.t[l]->get_etat() == ETATZERO ) {
	    cc2.t[l]->set_etat_zero() ; 
	}
	else {
	
	    assert( c2.t[l]->get_etat() == ETATQCQ ) ; 
	    cc2.t[l]->set_etat_qcq() ; 

	    // Copy of the coefficients of t2 
	    for (int k=0; k<np+1; k++) {
		for (int j=0; j<nt; j++) {
		    for (int i=0; i<nr; i++) {
			cc2.t[l]->set(k, j, i) = (*(c2.t[l]))(k, j, i) ; 
		    }
		}
	    }

	    // The extra phi coefficients are set to zero
	    for (int k=np+1; k<np2+2; k++) {
		for (int j=0; j<nt2; j++) {
		    for (int i=0; i<nr2; i++) {
			cc2.t[l]->set(k, j, i) = 0 ; 
		    }
		}
	    }

	    // The extra theta coefficients are set to zero
	    for (int k=0; k<np+1; k++) {
		for (int j=nt; j<nt2; j++) {
		    for (int i=0; i<nr2; i++) {
			cc2.t[l]->set(k, j, i) = 0 ; 
		    }
		}
	    }

	    // The extra r coefficients are set to zero
	    for (int k=0; k<np+1; k++) {
		for (int j=0; j<nt; j++) {
		    for (int i=nr; i<nr2; i++) {
			cc2.t[l]->set(k, j, i) = 0 ; 
		    }
		}
	    }
	    
	}

    }  // End of loop on the domains
    
    
    Valeur tt1( mg2 ) ; 
    Valeur tt2( mg2 ) ; 
    
    tt1 = cc1 ; 
    tt2 = cc2 ; 
    

    // Multiplication (in the configuration space) on the large grids
    // --------------------------------------------------------------
    
    tt1 = tt1 * tt2 ; 
    
    // Coefficients of the result
    // --------------------------
    
    tt1.coef() ; 
    
    tt1.ylm() ; 
    
    const Mtbl_cf& cr2 = *(tt1.c_cf) ; 
    
    Mtbl_cf cr(mg, cr2.base ) ; 
    
    cr.set_etat_qcq() ; 
    
    for (int l=0; l<mg.get_nzone(); l++) {

	if ( cr2.t[l]->get_etat() == ETATZERO ) {
	    
	    cr.t[l]->set_etat_zero() ; 
	    
	}
	else {
    
	    assert( cr2.t[l]->get_etat() == ETATQCQ ) ; 
	    
	    cr.t[l]->set_etat_qcq() ; 

	    int nr = mg.get_nr(l) ; 
	    int nt = mg.get_nt(l) ; 
	    int np = mg.get_np(l) ; 

	    // Copy of the coefficients of cr2
	    for (int k=0; k<np+1; k++) {
		for (int j=0; j<nt; j++) {
		    for (int i=0; i<nr; i++) {
			cr.t[l]->set(k, j, i) = (*(cr2.t[l]))(k, j, i) ; 
		    }
		}
	    }
    
	    // The last coefficient in phi is set to zero
	    for (int j=0; j<nt; j++) {
		for (int i=0; i<nr; i++) {
		    cr.t[l]->set(np+1, j, i) = 0 ; 
		}
	    }
	
	}
    
    }  // End of loop on the domains
        
    Valeur resu( mg ) ; 

    resu = cr ; 
    
    resu.ylm_i() ; 
    
    return resu ;
}

		    //---------------------------------------//
		    //	Multiplication with de-aliasing in r //
		    //---------------------------------------//

Valeur operator|(const Valeur& t1, const Valeur& t2)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mg() == t2.get_mg()) ;
    
    // Cas particuliers
    if (t1.get_etat() == ETATZERO) {
    	return t1 ;
    }
    if (t2.get_etat() == ETATZERO) {
    	return t2 ;
    }
    assert(t1.get_etat() == ETATQCQ) ;  // sinon...
    assert(t2.get_etat() == ETATQCQ) ;  // sinon...

    // Grid
    const Mg3d& mg = *(t1.get_mg()) ; 

    // Grid with twice the number of points in each dimension:
    const Mg3d& mg2 = *(mg.plus_half()) ; 
    
    // The coefficients are required
    if (t1.c_cf == 0x0) {
	t1.coef() ; 
    }
    if (t2.c_cf == 0x0) {
	t2.coef() ; 
    }

    const Mtbl_cf& c1 = *(t1.c_cf) ; 
    const Mtbl_cf& c2 = *(t2.c_cf) ; 

    assert( c1.get_etat() == ETATQCQ ) ; 
    assert( c2.get_etat() == ETATQCQ ) ; 

    // The number of coefficients is increased by 50% in r 
    //  and the additionnal coefficients are set to zero
    // ---------------------------------------------------
    
    Mtbl_cf cc1( mg2, c1.base ) ;
    Mtbl_cf cc2( mg2, c2.base ) ; 
    
    cc1.set_etat_qcq() ; 
    cc2.set_etat_qcq() ; 
    
    for (int l=0; l<mg.get_nzone(); l++) {

	int nr = mg.get_nr(l) ; 
	int nt = mg.get_nt(l) ; 
	int np = mg.get_np(l) ; 
	int nr2 = mg2.get_nr(l) ; 

	// Copy of the coefficients of t1
	// ------------------------------

	if ( c1.t[l]->get_etat() == ETATZERO ) {
	    cc1.t[l]->set_etat_zero() ; 
	}
	else {

	    assert( c1.t[l]->get_etat() == ETATQCQ ) ; 
	    cc1.t[l]->set_etat_qcq() ; 

	    // Copy of the coefficients of t1 
	    for (int k=0; k<np+1; k++) {
		for (int j=0; j<nt; j++) {
		    for (int i=0; i<nr; i++) {
			cc1.t[l]->set(k, j, i) = (*(c1.t[l]))(k, j, i) ; 
		    }
		}
	    }

	    // The extra phi coefficient is set to zero
	    for (int j=0; j<nt; j++) {
	      for (int i=0; i<nr2; i++) {
		cc1.t[l]->set(np+1, j, i) = 0 ; 
	      }
	    }


	    // The extra r coefficients are set to zero
	    for (int k=0; k<np+1; k++) {
		for (int j=0; j<nt; j++) {
		    for (int i=nr; i<nr2; i++) {
			cc1.t[l]->set(k, j, i) = 0 ; 
		    }
		}
	    }
	    
	}

	// Copy of the coefficients of t2
	// ------------------------------

	if ( c2.t[l]->get_etat() == ETATZERO ) {
	    cc2.t[l]->set_etat_zero() ; 
	}
	else {
	
	    assert( c2.t[l]->get_etat() == ETATQCQ ) ; 
	    cc2.t[l]->set_etat_qcq() ; 

	    // Copy of the coefficients of t2 
	    for (int k=0; k<np+1; k++) {
		for (int j=0; j<nt; j++) {
		    for (int i=0; i<nr; i++) {
			cc2.t[l]->set(k, j, i) = (*(c2.t[l]))(k, j, i) ; 
		    }
		}
	    }

	    // The extra phi coefficient is set to zero
	    for (int j=0; j<nt; j++) {
	      for (int i=0; i<nr2; i++) {
		cc2.t[l]->set(np+1, j, i) = 0 ; 
	      }
	    }

	    // The extra r coefficients are set to zero
	    for (int k=0; k<np+1; k++) {
		for (int j=0; j<nt; j++) {
		    for (int i=nr; i<nr2; i++) {
			cc2.t[l]->set(k, j, i) = 0 ; 
		    }
		}
	    }
	    
	}

    }  // End of loop on the domains
    
    
    Valeur tt1( mg2 ) ; 
    Valeur tt2( mg2 ) ; 
    
    tt1 = cc1 ; 
    tt2 = cc2 ; 
    
    // Multiplication (in the configuration space) on the large grids
    // --------------------------------------------------------------
    
    tt1 = tt1 * tt2 ; 
    
    // Coefficients of the result
    // --------------------------
    
    tt1.coef() ; 
    
    const Mtbl_cf& cr2 = *(tt1.c_cf) ; 
    
    Mtbl_cf cr(mg, cr2.base ) ; 
    
    cr.set_etat_qcq() ; 
    
    for (int l=0; l<mg.get_nzone(); l++) {

	if ( cr2.t[l]->get_etat() == ETATZERO ) {
	    
	    cr.t[l]->set_etat_zero() ; 
	    
	}
	else {
    
	    assert( cr2.t[l]->get_etat() == ETATQCQ ) ; 
	    
	    cr.t[l]->set_etat_qcq() ; 

	    int nr = mg.get_nr(l) ; 
	    int nt = mg.get_nt(l) ; 
	    int np = mg.get_np(l) ; 

	    // Copy of the coefficients of cr2
	    for (int k=0; k<np+1; k++) {
		for (int j=0; j<nt; j++) {
		    for (int i=0; i<nr; i++) {
			cr.t[l]->set(k, j, i) = (*(cr2.t[l]))(k, j, i) ; 
		    }
		}
	    }
    
	    // The last coefficient in phi is set to zero
	    for (int j=0; j<nt; j++) {
		for (int i=0; i<nr; i++) {
		    cr.t[l]->set(np+1, j, i) = 0 ; 
		}
	    }
	
	}
    
    }  // End of loop on the domains
        
    Valeur resu( mg ) ; 

    resu = cr ; 
    
    return resu ;
}

}
