/*
 * Methods for making calculations with Godunov-type arrays.
 *
 * See the file tbl_val.h for documentation
 *
 */

/*
 *   Copyright (c) 2001 Jerome Novak
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


char TBL_VAL_MATH_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:53:48  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:22  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2002/11/12 10:03:54  j_novak
 * The method "Tbl_val::get_gval" has been changed to "get_grid".
 *
 * Revision 1.1  2001/11/22 13:41:54  j_novak
 * Added all source files for manipulating Valencia type objects and making
 * interpolations to and from Meudon grids.
 *
 *
 * $Header$
 *
 */

// Headers C
// ---------
#include <cmath>
#include <cstdlib>

// Headers Lorene
// --------------
#include "tbl_val.h"

//-------//
// Sinus //
//-------//

namespace Lorene {
Tbl_val sin(const Tbl_val& ti)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    return ti ;
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...
  Tbl_val to(ti.get_grille()) ;			// Tbl_val resultat
  to.set_etat_qcq() ;
  int taille = ti.get_taille() ;
  for (int i=0 ; i<taille ; i++) {
    to.t[i] = sin(ti.t[i]) ;
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = sin(ti.tzri[i]) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = sin(ti.txti[i]) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = sin(ti.typi[i]) ;
  return to ;
}

//---------//
// Cosinus //
//---------//

Tbl_val cos(const Tbl_val& ti)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
  Tbl_val to(ti.get_grille()) ;		    // Tbl_val resultat
  to.set_etat_qcq() ;

    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    to = 1 ;
    return to ;
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...
  int taille = ti.get_taille() ;
  for (int i=0 ; i<taille ; i++) {
    to.t[i] = cos(ti.t[i]) ;
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = cos(ti.tzri[i]) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = cos(ti.txti[i]) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = cos(ti.typi[i]) ;
  return to ;
}

//----------//
// Tangente //
//----------//

Tbl_val tan(const Tbl_val& ti)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    return ti ;
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...
  Tbl_val to(ti.get_grille()) ;		    // Tbl_val resultat
  to.set_etat_qcq() ;
  int taille = ti.get_taille() ;
  for (int i=0 ; i<taille ; i++) {
    to.t[i] = tan(ti.t[i]) ;
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = tan(ti.tzri[i]) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = tan(ti.txti[i]) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = tan(ti.typi[i]) ;
  return to ;
}

//----------//
// ArcSinus //
//----------//

Tbl_val asin(const Tbl_val& ti)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    return ti ;
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...
  Tbl_val to(ti.get_grille()) ;		    // Tbl_val resultat
  to.set_etat_qcq() ;
  int taille = ti.get_taille() ;
  for (int i=0 ; i<taille ; i++) {
    to.t[i] = asin(ti.t[i]) ;
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = asin(ti.tzri[i]) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = asin(ti.txti[i]) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = asin(ti.typi[i]) ;

  return to ;
}

//------------//
// ArcCosinus //
//------------//

Tbl_val acos(const Tbl_val& ti)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
  Tbl_val to(ti.get_grille()) ;		    // Tbl_val resultat
  to.set_etat_qcq() ;

    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    to = M_PI * .5 ;
    return to ;
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...
  int taille = ti.get_taille() ;
  for (int i=0 ; i<taille ; i++) {
    to.t[i] = acos(ti.t[i]) ;
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = acos(ti.tzri[i]) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = acos(ti.txti[i]) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = acos(ti.typi[i]) ;
  return to ;
}

//-------------//
// ArcTangente //
//-------------//

Tbl_val atan(const Tbl_val& ti)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    return ti ;
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...
  Tbl_val to(ti.get_grille()) ;		    // Tbl_val resultat
  to.set_etat_qcq() ;
  int taille = ti.get_taille() ;
  for (int i=0 ; i<taille ; i++) {
    to.t[i] = atan(ti.t[i]) ;
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = atan(ti.tzri[i]) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = atan(ti.txti[i]) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = atan(ti.typi[i]) ;
  return to ;
}

//------//
// Sqrt //
//------//

Tbl_val sqrt(const Tbl_val& ti)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    return ti ;
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...
  Tbl_val to(ti.get_grille()) ;		    // Tbl_val resultat
  to.set_etat_qcq() ;
  int taille = ti.get_taille() ;
  for (int i=0 ; i<taille ; i++) {
    to.t[i] = sqrt(ti.t[i]) ;
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = sqrt(ti.tzri[i]) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = sqrt(ti.txti[i]) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = sqrt(ti.typi[i]) ;
  return to ;
}

//---------------//
// Exponentielle //
//---------------//

Tbl_val exp(const Tbl_val& ti)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
  Tbl_val to(ti.get_grille()) ;		    // Tbl_val resultat
  to.set_etat_qcq() ;

    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    to = 1 ;
    return to ;
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...
  int taille = ti.get_taille() ;
  for (int i=0 ; i<taille ; i++) {
    to.t[i] = exp(ti.t[i]) ;
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = exp(ti.tzri[i]) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = exp(ti.txti[i]) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = exp(ti.typi[i]) ;
  return to ;
}

//-------------//
// Log naturel //
//-------------//

Tbl_val log(const Tbl_val& ti)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    cout << "Tbl_val log: log(ETATZERO) !" << endl  ;
    abort () ;
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...
  Tbl_val to(ti.get_grille()) ;		    // Tbl_val resultat
  to.set_etat_qcq() ;
  int taille = ti.get_taille() ;
  for (int i=0 ; i<taille ; i++) {
    to.t[i] = log(ti.t[i]) ;
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = log(ti.tzri[i]) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = log(ti.txti[i]) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = log(ti.typi[i]) ;
  return to ;
}

//-------------//
// Log decimal //
//-------------//

Tbl_val log10(const Tbl_val& ti)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    cout << "Tbl_val log10: log10(ETATZERO) !" << endl ;
    abort () ;
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...
  Tbl_val to(ti.get_grille()) ;		    // Tbl_val resultat
  to.set_etat_qcq() ;
  int taille = ti.get_taille() ;
  for (int i=0 ; i<taille ; i++) {
    to.t[i] = log10(ti.t[i]) ;
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = log(ti.tzri[i]) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = log(ti.txti[i]) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = log(ti.typi[i]) ;
  return to ;
}

//--------------//
// Power entier //
//--------------//

Tbl_val pow(const Tbl_val& ti, int n)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    if (n > 0) {
      return ti ;
    }
    else {
      cout << "Tbl_val pow: ETATZERO^n avec n<=0 ! "<< endl  ;
      abort () ;
    }
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...
  Tbl_val to(ti.get_grille()) ;		    // Tbl_val resultat
  to.set_etat_qcq() ;
  double x = n ;
  int taille = ti.get_taille() ;
  for (int i=0 ; i<taille ; i++) {
    to.t[i] = pow(ti.t[i], x) ;		
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = pow(ti.tzri[i], x) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = pow(ti.txti[i], x) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = pow(ti.typi[i], x) ;
  return to ;
}

//--------------//
// Power double //
//--------------//

Tbl_val pow(const Tbl_val& ti, double x)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    if (x > 0) {
      return ti ;
    }
    else {
      cout << "Tbl_val pow: ETATZERO^x avec x<=0 !" << endl ;
      abort () ;
    }
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...
  Tbl_val to(ti.get_grille()) ;		    // Tbl_val resultat
  to.set_etat_qcq() ;
  int taille = ti.get_taille() ;
  for (int i=0 ; i<taille ; i++) {
    to.t[i] = pow(ti.t[i], x) ;
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = pow(ti.tzri[i], x) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = pow(ti.txti[i], x) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = pow(ti.typi[i], x) ;
  return to ;
}

//----------------//
// Absolute value //
//----------------//

Tbl_val abs(const Tbl_val& ti)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    return ti ;
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...

  Tbl_val to(ti.get_grille()) ;		    // Tbl_val resultat
  to.set_etat_qcq() ;

  const double* xi = ti.t ; 
  double* xo = to.t ; 
  int taille = ti.get_taille() ;

  for (int i=0 ; i<taille ; i++) {
    xo[i] = fabs( xi[i] ) ;
  }
  for (int i=0; i<ti.get_taille_i(0); i++)
    to.tzri[i] = fabs(ti.tzri[i]) ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    to.txti[i] = fabs(ti.txti[i]) ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    to.typi[i] = fabs(ti.typi[i]) ;
    
  return to ;
}
//----------------//
//	    Cubic     //
//----------------//

Tbl_val racine_cubique(const Tbl_val& ti)
{
  // Protection
  assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
  if (ti.get_etat() == ETATZERO) {
    return ti ;
  }
    
  // Cas general
  assert(ti.get_etat() == ETATQCQ) ;	// sinon...

  Tbl_val absolute(abs(ti)) ;
  Tbl_val res (pow(absolute, 1./3.)) ;
    
  for (int i=0 ; i<ti.get_taille() ; i++)
    if (ti.t[i] < 0)
      res.t[i] *= -1 ;
  for (int i=0; i<ti.get_taille_i(0); i++)
    if (ti.tzri[i] < 0) res.tzri[i] *= -1 ;
  if (ti.txti != 0x0) for (int i=0; i<ti.get_taille_i(1); i++)
    if (ti.txti[i] < 0) res.txti[i] *= -1 ;
  if (ti.typi != 0x0) for (int i=0; i<ti.get_taille_i(2); i++)
    if (ti.typi[i] < 0) res.typi[i] *= -1 ;
    
  return res ;
}

//-------------------------------//
//            max                //
//-------------------------------//

double max(const Tbl_val& ti) {
    
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

double min(const Tbl_val& ti) {
    
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

double norme(const Tbl_val& ti) {

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

double diffrel(const Tbl_val& t1, const Tbl_val& t2) {
    
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

double diffrelmax(const Tbl_val& t1, const Tbl_val& t2) {
    
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
}
