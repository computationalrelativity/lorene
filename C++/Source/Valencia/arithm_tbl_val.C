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


char ARITHM_TBL_VAL_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/11/12 10:03:54  j_novak
 * The method "Tbl_val::get_gval" has been changed to "get_grid".
 *
 * Revision 1.2  2002/10/16 14:37:15  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1  2001/11/22 13:41:54  j_novak
 * Added all source files for manipulating Valencia type objects and making
 * interpolations to and from Meudon grids.
 *
 *
 * $Header$
 *
 */

// headers Lorene
#include "tbl_val.h"

			//********************//
			// OPERATEURS UNAIRES //
			//********************//

// + Tbl_val
// -----
Tbl_val operator+(const Tbl_val& t1)
{
  // Protection
  assert(t1.get_etat() != ETATNONDEF) ;
  
  return t1 ;
}

// - Tbl_val
// -----
Tbl_val operator-(const Tbl_val& t1)
{
  // Protection
  assert(t1.get_etat() != ETATNONDEF) ;
  
  // Cas particulier
  if (t1.get_etat() == ETATZERO) {
    return t1 ;
  }
  
  // Cas general
  Tbl_val r(t1.get_grille()) ;		    // Tbl_val resultat
  r.set_etat_qcq() ;
  for (int i=0 ; i<r.get_taille() ; i++) 
    (r.t)[i] = - (t1.t)[i] ;
  
  for (int i=0 ; i<r.get_taille_i(0) ; i++) 
    (r.tzri)[i] = - (t1.tzri)[i] ;
  
  if (t1.txti != 0x0) for (int i=0 ; i<r.get_taille_i(1) ; i++) 
    (r.txti)[i] = - (t1.txti)[i] ;
  
  if (t1.typi != 0x0) for (int i=0 ; i<r.get_taille_i(2) ; i++) 
    (r.typi)[i] = - (t1.typi)[i] ;
  
  return r ;
}

                        //**********//
			// ADDITION //
			//**********//

// Tbl_val + Tbl_val
// ---------
Tbl_val operator+(const Tbl_val& t1, const Tbl_val& t2)
{

  // Protection
  assert(t1.get_etat() != ETATNONDEF) ;
  assert(t2.get_etat() != ETATNONDEF) ;
  assert(t1.get_ndim() == t2.get_ndim()) ;
  for (int i=0 ; i<t1.get_ndim() ; i++) {
    assert( t1.get_dim(i) == t2.get_dim(i) ) ;
  }
  
  // Traitement des cas particuliers
  if (t1.get_etat() == ETATZERO) {
    return t2 ;
  }
  if (t2.get_etat() == ETATZERO) {
    return t1 ;
  }
  
  // Cas general
  assert(t1.get_etat() == ETATQCQ) ;	// sinon...
  assert(t2.get_etat() == ETATQCQ) ;	// sinon...
  
  Tbl_val r(t1) ;	    // Tbl_val resultat
  for (int i=0 ; i<r.get_taille() ; i++) {
    (r.t)[i] += (t2.t)[i] ;
  }
  for (int i=0 ; i<r.get_taille_i(0) ; i++) 
    (r.tzri)[i] += (t2.tzri)[i] ;
  
  if (t1.txti != 0x0) for (int i=0 ; i<r.get_taille_i(1) ; i++) 
    (r.txti)[i] *= (t2.txti)[i] ;
  
  if (t1.typi != 0x0) for (int i=0 ; i<r.get_taille_i(2) ; i++) 
    (r.typi)[i] *= (t2.typi)[i] ;
  
  
  // Termine
  return r ;
}

// Tbl_val + double
// ------------
Tbl_val operator+(const Tbl_val& t1, double x)
{
  // Protection
  assert(t1.get_etat() != ETATNONDEF) ;
  
  // Cas particulier
  if ( x == double(0) ) {
    return t1 ;
  }
  
  // Cas general
  Tbl_val r(t1) ;		// Tbl_val resultat
  r.set_etat_qcq() ;
  for (int i=0 ; i<r.get_taille() ; i++) {
    (r.t)[i] += x ;
  }
  for (int i=0 ; i<r.get_taille_i(0) ; i++) 
    (r.tzri)[i] += x ;
  
  if (t1.txti != 0x0) for (int i=0 ; i<r.get_taille_i(1) ; i++) 
    (r.txti)[i] += x ;
  
  if (t1.typi != 0x0) for (int i=0 ; i<r.get_taille_i(2) ; i++) 
    (r.typi)[i] += x;
  
  return r ;
}

// double + Tbl_val
// ------------
Tbl_val operator+(double x, const Tbl_val& t1)
{
  return t1 + x ;
}

// Tbl_val + int
// ---------
Tbl_val operator+(const Tbl_val& t1, int n)
{
  return t1 + double(n) ;
}

// int + Tbl_val
// ---------
Tbl_val operator+(int n, const Tbl_val& t1)
{
  return t1 + double(n) ;
}


			//**************//
			// SOUSTRACTION //
			//**************//

// Tbl_val - Tbl_val
// ---------
Tbl_val operator-(const Tbl_val& t1, const Tbl_val& t2)
{

  // Protection
  assert(t1.get_etat() != ETATNONDEF) ;
  assert(t2.get_etat() != ETATNONDEF) ;
  assert(t1.get_ndim() == t2.get_ndim()) ;
  for (int i=0 ; i<t1.get_ndim() ; i++) {
    assert( t1.get_dim(i) == t2.get_dim(i) ) ;
  }
  
  // Traitement des cas particuliers
  if (t1.get_etat() == ETATZERO) {
    return -t2 ;
  }
  if (t2.get_etat() == ETATZERO) {
    return t1 ;
  }
  
  // Cas general
  assert(t1.get_etat() == ETATQCQ) ;	// sinon...
  assert(t2.get_etat() == ETATQCQ) ;	// sinon...

  Tbl_val r(t1) ;	    // Tbl_val resultat
  for (int i=0 ; i<r.get_taille() ; i++) {
    (r.t)[i] -= (t2.t)[i] ;
  }
  for (int i=0 ; i<r.get_taille_i(0) ; i++) 
    (r.tzri)[i] -= (t2.tzri)[i] ;
    
  if (t1.txti != 0x0) for (int i=0 ; i<r.get_taille_i(1) ; i++) 
    (r.txti)[i] -= (t2.txti)[i] ;
    
  if (t1.typi != 0x0) for (int i=0 ; i<r.get_taille_i(2) ; i++) 
    (r.typi)[i] -= (t2.typi)[i] ;
    

  // Termine
  return r ;
}


// Tbl_val - double
// ------------
Tbl_val operator-(const Tbl_val& t1, double x)
{
  // Protection
  assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
  if ( x == double(0) ) {
    return t1 ;
  }
    
  // Cas general
  Tbl_val r(t1) ;		// Tbl_val resultat
  r.set_etat_qcq() ;
  for (int i=0 ; i<r.get_taille() ; i++) {
    (r.t)[i] -= x ;
  }
  for (int i=0 ; i<r.get_taille_i(0) ; i++) 
    (r.tzri)[i] -= x ;
    
  if (t1.txti != 0x0) for (int i=0 ; i<r.get_taille_i(1) ; i++) 
    (r.txti)[i] -= x ;
    
  if (t1.typi != 0x0) for (int i=0 ; i<r.get_taille_i(2) ; i++) 
    (r.typi)[i] -= x ;
    

  return r ;
}

// Tbl_val - int
// ---------
Tbl_val operator-(const Tbl_val& t1, int n)
{
  return t1 - double(n) ;
}

// double - Tbl_val
// ------------
Tbl_val operator-(double x, const Tbl_val& t1)
{
  // Protection
  assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
  if ( x == double(0) ) {
    return -t1 ;
  }
    
  // Cas general
  Tbl_val r(t1) ;		// Tbl_val resultat
  r.set_etat_qcq() ;
  for (int i=0 ; i<r.get_taille() ; i++) {
    (r.t)[i] -= x ;
  }
  for (int i=0 ; i<r.get_taille_i(0) ; i++) 
    (r.tzri)[i] -= x ;
    
  if (t1.txti != 0x0) for (int i=0 ; i<r.get_taille_i(1) ; i++) 
    (r.txti)[i] -= x ;
    
  if (t1.typi != 0x0) for (int i=0 ; i<r.get_taille_i(2) ; i++) 
    (r.typi)[i] -= x ;
    
  return -r ;
}

// int - Tbl_val
// ---------
Tbl_val operator-(int n, const Tbl_val& t1)
{
  return double(n) - t1 ;
}

//****************//
// MULTIPLICATION //
//****************//

// Tbl_val * Tbl_val
// ---------
Tbl_val operator*(const Tbl_val& t1, const Tbl_val& t2)
{
  // Protection
  assert(t1.get_etat() != ETATNONDEF) ;
  assert(t2.get_etat() != ETATNONDEF) ;
  assert(t1.get_ndim() == t2.get_ndim()) ;
  for (int i=0 ; i<t1.get_ndim() ; i++) {
    assert( t1.get_dim(i) == t2.get_dim(i) ) ;
  }

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

  Tbl_val r(t1) ;
  for (int i=0 ; i<r.get_taille() ; i++) {
    (r.t)[i] *= (t2.t)[i] ;
  }
  for (int i=0 ; i<r.get_taille_i(0) ; i++) 
    (r.tzri)[i] *= (t2.tzri)[i] ;
    
  if (t1.txti != 0x0) for (int i=0 ; i<r.get_taille_i(1) ; i++) 
    (r.txti)[i] *= (t2.txti)[i] ;
    
  if (t1.typi != 0x0) for (int i=0 ; i<r.get_taille_i(2) ; i++) 
    (r.typi)[i] *= (t2.typi)[i] ;
    
    // Termine
  return r ;
}

// Tbl_val * double
// ------------
Tbl_val operator*(const Tbl_val& t1, double x)
{
  // Protection
  assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
  if ((t1.get_etat() == ETATZERO) || ( x == double(1) )) {
    return t1 ;
  }

  // Cas general
  assert(t1.get_etat() == ETATQCQ) ;	// sinon...

  Tbl_val r(t1) ;		    // Tbl_val resultat

  if (x == double(0)) {
    r.set_etat_zero() ;
  }
  else {
    for (int i=0 ; i<r.get_taille() ; i++) {
      (r.t)[i] *= x ;
    }
    for (int i=0 ; i<r.get_taille_i(0) ; i++) 
      (r.tzri)[i] *= x ;
    
    if (t1.txti != 0x0) for (int i=0 ; i<r.get_taille_i(1) ; i++) 
      (r.txti)[i] *= x ;
    
    if (t1.typi != 0x0) for (int i=0 ; i<r.get_taille_i(2) ; i++) 
      (r.typi)[i] *= x ;
    

  }
    
  // Termine
  return r ;
}

// double * Tbl_val
// ------------
Tbl_val operator*(double x, const Tbl_val& t1)
{
  return t1 * x ;
}

// Tbl_val * int
// ---------
Tbl_val operator*(const Tbl_val& t1, int n)
{
  return t1 * double(n) ;
}

// int * Tbl_val
// ---------
Tbl_val operator*(int n, const Tbl_val& t1)
{
  return t1 * double(n) ;
}

//**********//
// DIVISION //
//**********//

// Tbl_val / Tbl_val
// ---------
Tbl_val operator/(const Tbl_val& t1, const Tbl_val& t2)
{
  // Protection
  assert(t1.get_etat() != ETATNONDEF) ;
  assert(t2.get_etat() != ETATNONDEF) ;
  assert(t1.get_ndim() == t2.get_ndim()) ;
  for (int i=0 ; i<t1.get_ndim() ; i++) {
    assert( t1.get_dim(i) == t2.get_dim(i) ) ;
  }

  // Cas particuliers
  if (t2.get_etat() == ETATZERO) {
    cout << "Division by 0 in Tbl_val/Tbl_val !" << endl ;
    abort() ; 
  }
  if (t1.get_etat() == ETATZERO) {
    return t1 ;
  }
    
  // Cas general
  assert(t1.get_etat() == ETATQCQ) ;	// sinon...
  assert(t2.get_etat() == ETATQCQ) ;	// sinon...

  Tbl_val r(t1) ;		    // Tbl_val resultat
  for (int i=0 ; i<r.get_taille() ; i++) {
    (r.t)[i] /= (t2.t)[i] ;
  }
  for (int i=0 ; i<r.get_taille_i(0) ; i++) 
    (r.tzri)[i] /= (t2.tzri)[i] ;
    
  if (t1.txti != 0x0) for (int i=0 ; i<r.get_taille_i(1) ; i++) 
    (r.txti)[i] /= (t2.txti)[i] ;
    
  if (t1.typi != 0x0) for (int i=0 ; i<r.get_taille_i(2) ; i++) 
    (r.typi)[i] /= (t2.typi)[i] ;
        

    // Termine
  return r ;
}

// Tbl_val / double
// ------------
Tbl_val operator/(const Tbl_val& t1, double x)
{
  // Protection
  assert(t1.get_etat() != ETATNONDEF) ;
  if ( x == double(0) ) {
    cout << "Division by 0 in Tbl_val/double !" << endl ;
    abort() ;
  }
    
  // Cas particulier
  if ((t1.get_etat() == ETATZERO) || ( x == double(1) )) {
    return t1 ;
  }
    
  // Cas general
  assert(t1.get_etat() == ETATQCQ) ;	// sinon...

  Tbl_val r(t1) ;		    // Tbl_val resultat
  for (int i=0 ; i<r.get_taille() ; i++) {
    (r.t)[i] /= x ;
  }
  for (int i=0 ; i<r.get_taille_i(0) ; i++) 
    (r.tzri)[i] /= x ;
    
  if (t1.txti != 0x0) for (int i=0 ; i<r.get_taille_i(1) ; i++) 
    (r.txti)[i] /= x ;
    
  if (t1.typi != 0x0) for (int i=0 ; i<r.get_taille_i(2) ; i++) 
    (r.typi)[i] /= x ;
    


  return r ;
}

// Tbl_val / int
// ---------
Tbl_val operator/(const Tbl_val& t1, int n)
{
  return t1 / double(n) ;
}

// double / Tbl_val
// ------------
Tbl_val operator/(double x, const Tbl_val& t1)
{
  // Protection
  assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particuliers
  if (t1.get_etat() == ETATZERO) {
    cout << "Division by 0 in double/Tbl_val !" << endl ;
    abort() ; 
  }

  // Cas general
  assert(t1.get_etat() == ETATQCQ) ;	// sinon...

  Tbl_val r(t1.get_grille()) ;		// Tbl_val resultat, a priori NONDEF
    
  if ( x == double(0) ) {
    r.set_etat_zero() ;
  }
  else {
    r.set_etat_qcq() ;
    for (int i=0 ; i<r.get_taille() ; i++) {
      (r.t)[i] = x / (t1.t)[i] ;
    }
    for (int i=0 ; i<r.get_taille_i(0) ; i++) 
      (r.tzri)[i] = x / (t1.tzri)[i] ;
    
    if (t1.txti != 0x0) for (int i=0 ; i<r.get_taille_i(1) ; i++) 
      (r.txti)[i] = x / (t1.txti)[i] ;
    
    if (t1.typi != 0x0) for (int i=0 ; i<r.get_taille_i(2) ; i++) 
      (r.typi)[i] = x / (t1.typi)[i] ;
    

  }
    
  // Termine
  return r ;
}

// int / Tbl_val
// ---------
Tbl_val operator/(int n, const Tbl_val& t1)
{
  return double(n) / t1 ;
}

//*******************//
// operateurs +=,... //
//*******************//

void Tbl_val::operator+=(const Tbl_val & ti) {

  // Protection
  assert(gval = ti.gval) ;
  assert(etat != ETATNONDEF) ;
  assert(ti.get_etat() != ETATNONDEF) ;
    
  // Cas particulier
  if (ti.get_etat() == ETATZERO) {
    return ;
  }
    
  // Cas general
  int n = get_taille() ;
  switch(etat) {
  case ETATZERO:
    set_etat_qcq() ;
    for (int i=0 ; i<n ; i++) {
      t[i] = ti.t[i] ;
    }
    for (int i=0 ; i < get_taille_i(0) ; i++) 
      tzri[i] = ti.tzri[i] ;
    
    if (ti.txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
      txti[i] = ti.txti[i] ;
    
    if (ti.typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
      typi[i] = ti.typi[i] ;

    break ;
	
  case ETATQCQ:
    for (int i=0 ; i<n ; i++) {
      t[i] += ti.t[i] ;
    }
    for (int i=0 ; i < get_taille_i(0) ; i++) 
      tzri[i] += ti.tzri[i] ;
    
    if (ti.txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
      txti[i] += ti.txti[i] ;
    
    if (ti.typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
      typi[i] += ti.typi[i] ;
    break ;
	
  default:
    cout << "etat inconnu " << __FILE__ << endl ;
    abort() ;
    break ;
  }
    
  // Termine
}

void Tbl_val::operator+=(double x) {

  // Protection
  assert(etat != ETATNONDEF) ;

  // Cas particulier
  if ( x == double(0) ) {
    return ;
  }
    
  // Cas general
  int n = get_taille() ;
  switch(etat) {
  case ETATZERO:
    set_etat_qcq() ;
    for (int i=0 ; i<n ; i++) {
      t[i] = x ;
    }
    for (int i=0 ; i < get_taille_i(0) ; i++) 
      tzri[i] = x ;
    
    if (txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
      txti[i] = x ;
    
    if (typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
      typi[i] = x ;
    break ;
	
  case ETATQCQ:
    for (int i=0 ; i<n ; i++) {
      t[i] += x ;
    }
    for (int i=0 ; i < get_taille_i(0) ; i++) 
      tzri[i] += x ;
    
    if (txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
      txti[i] += x ;
    
    if (typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
      typi[i] += x ;

    break ;
	
  default:
    cout << "etat inconnu " << __FILE__ << endl ;
    abort() ;
    break ;
  }
    
  // Termine
}

void Tbl_val::operator-=(const Tbl_val & ti) {

  // Protection
  assert(gval = ti.gval) ;
  assert(etat != ETATNONDEF) ;
  assert(ti.get_etat() != ETATNONDEF) ;
    
  // Cas particulier
  if (ti.get_etat() == ETATZERO) {
    return ;
  }
    
  // Cas general
  int n = get_taille() ;
  switch(etat) {
  case ETATZERO:
    set_etat_qcq() ;
    for (int i=0 ; i<n ; i++) {
      t[i] = - ti.t[i] ;
    }
    for (int i=0 ; i < get_taille_i(0) ; i++) 
      tzri[i] = - ti.tzri[i] ;
    
    if (ti.txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
      txti[i] = - ti.txti[i] ;
    
    if (ti.typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
      typi[i] = - ti.typi[i] ;

    break ;
	
  case ETATQCQ:
    for (int i=0 ; i<n ; i++) {
      t[i] -= ti.t[i] ;
    }
    for (int i=0 ; i < get_taille_i(0) ; i++) 
      tzri[i] -= ti.tzri[i] ;
    
    if (ti.txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
      txti[i] -= ti.txti[i] ;
    
    if (ti.typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
      typi[i] -= ti.typi[i] ;

    break ;
	
  default:
    cout << "etat inconnu " << __FILE__ << endl ;
    abort() ;
    break ;
  }
    
  // Termine
}

void Tbl_val::operator-=(double x) {

  // Protection
  assert(etat != ETATNONDEF) ;

  // Cas particulier
  if ( x == double(0) ) {
    return ;
  }
    
  // Cas general
  int n = get_taille() ;
  switch(etat) {
  case ETATZERO:
    set_etat_qcq() ;
    for (int i=0 ; i<n ; i++) {
      t[i] = - x ;
    }
    for (int i=0 ; i < get_taille_i(0) ; i++) 
      tzri[i] = - x ;
    
    if (txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
      txti[i] = - x ;
    
    if (typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
      typi[i] = - x ;

    break ;
	
  case ETATQCQ:
    for (int i=0 ; i<n ; i++) {
      t[i] -= x ;
    }
    for (int i=0 ; i < get_taille_i(0) ; i++) 
      tzri[i] -= x ;
    
    if (txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
      txti[i] -= x ;
    
    if (typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
      typi[i] -= x ;

    break ;
	
  default:
    cout << "etat inconnu " << __FILE__ << endl ;
    abort() ;
    break ;
  }
    
  // Termine
}

void Tbl_val::operator*=(const Tbl_val & ti) {

  // Protection
  assert(gval = ti.gval) ;
  assert(etat != ETATNONDEF) ;
  assert(ti.get_etat() != ETATNONDEF) ;
    
  // Cas particulier
  if (etat == ETATZERO) {
    return ;
  }
  if (ti.get_etat() == ETATZERO) {
    set_etat_zero() ;
    return ;
  }
    
  // Cas general
  assert(etat == ETATQCQ) ;
  for (int i=0 ; i<get_taille() ; i++) {
    t[i] *= ti.t[i] ;
  }
  for (int i=0 ; i < get_taille_i(0) ; i++) 
    tzri[i] *= ti.tzri[i] ;
    
  if (ti.txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
    txti[i] *= ti.txti[i] ;
    
  if (ti.typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
    typi[i] *= ti.typi[i] ;

  // Termine
}

void Tbl_val::operator*=(double x) {

  // Protection
  assert(etat != ETATNONDEF) ;

  // Cas particulier
  if ( x == double(0) ) {
    set_etat_zero() ;
    return ;
  }
  if (etat == ETATZERO) {
    return ;
  }
    
  // Cas general
  assert(etat == ETATQCQ) ;
  for (int i=0 ; i<get_taille() ; i++) {
    t[i] *= x ;
  }
  for (int i=0 ; i < get_taille_i(0) ; i++) 
    tzri[i] *= x ;
    
  if (txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
    txti[i] *= x ;
    
  if (typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
    typi[i] *= x ;
    
  // Termine
}

void Tbl_val::operator/=(const Tbl_val & ti) {

  // Protection
  assert(gval = ti.gval) ;
  assert(etat != ETATNONDEF) ;
  assert(ti.get_etat() != ETATNONDEF) ;
    
  // Cas particulier
  if (ti.get_etat() == ETATZERO) {
    cout << "Division by 0 in Tbl_val::operator/=(const Tbl_val &) !" << endl ;
    abort() ;
  }
  if (etat == ETATZERO) {
    return ;
  }
    
  // Cas general
  assert(etat == ETATQCQ) ;
  assert(ti.get_etat() == ETATQCQ) ;
  int n = get_taille() ;
  for (int i=0 ; i<n ; i++) {
    t[i] /= ti.t[i] ;
  }
  for (int i=0 ; i < get_taille_i(0) ; i++) 
    tzri[i] /= ti.tzri[i] ;
    
  if (ti.txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
    txti[i] /= ti.txti[i] ;
    
  if (ti.typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
    typi[i] /= ti.typi[i] ;

  // Termine
}

void Tbl_val::operator/=(double x) {

  // Protection
  assert(etat != ETATNONDEF) ;

  // Cas particulier
  if ( x == double(0) ) {
    cout << "Division by 0 in Tbl_val::operator/=(double ) !" << endl ;
    abort() ;
  }
  if (etat == ETATZERO) {
    return ;
  }
    
  // Cas general
  assert(etat == ETATQCQ) ;
  int n = get_taille() ;
  for (int i=0 ; i<n ; i++) {
    t[i] /= x ;
  }
  for (int i=0 ; i < get_taille_i(0) ; i++) 
    tzri[i] /= x ;
    
  if (txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
    txti[i] /= x ;
    
  if (typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
    typi[i] /= x ;

  // Termine
}

