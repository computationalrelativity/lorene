/*
 * Methods for the class Tbl_val.
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


char TBL_VAL_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2014/10/06 15:13:22  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2008/02/18 13:53:48  j_novak
 * Removal of special indentation instructions.
 *
 * Revision 1.4  2007/11/02 15:45:58  j_novak
 * Added an ugly method "append_array", which substitutes the argument to the
 * main array t.
 *
 * Revision 1.3  2002/10/16 14:37:15  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2001/12/04 21:27:54  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1  2001/11/22 13:41:54  j_novak
 * Added all source files for manipulating Valencia type objects and making
 * interpolations to and from Meudon grids.
 *
 *
 * $Header$
 *
 */

// headers C
#include <cmath>

// headers Lorene
#include "headcpp.h"
#include "tbl_val.h"
#include "utilitaires.h"


			//---------------//
			// Constructeurs //
			//---------------//


// Constructeur a partir d'une grille de Valence
Tbl_val::Tbl_val(const Grille_val* g) : etat(ETATNONDEF), 
  dim(g->get_dim_tbl()), gval(g), t(0x0), tzri(0x0), txti(0x0), typi(0x0) {}

// Copie
Tbl_val::Tbl_val(const Tbl_val& tc) : etat(tc.etat), dim(tc.dim), 
  gval(tc.gval)  {

  // La valeur eventuelle
  if (tc.etat == ETATQCQ) {
    t = new double[get_taille()] ;
    for (int i=0 ; i<get_taille() ; i++) {
      t[i] = tc.t[i] ;
    }
    
    tzri = new double[get_taille_i(0)] ;
    for (int i=0 ; i<get_taille_i(0) ; i++) {
      tzri[i] = tc.tzri[i] ;
    }
    
    if (get_ndim() > 1) {
      txti= new double[get_taille_i(1)] ;
      for (int i=0 ; i<get_taille_i(1) ; i++) txti[i] = tc.txti[i] ;
    }
    else txti = 0x0 ;
    if (get_ndim() > 2) {
      typi = new double[get_taille_i(2)] ;
      for (int i=0; i<get_taille_i(2); i++) typi[i] = tc.typi[i] ;
    } 
    else typi = 0x0 ;
  }
  else{
    t = 0x0 ; 
    tzri = 0x0 ;
    txti = 0x0 ;
    typi = 0x0 ;
  }   
}

// From file
Tbl_val::Tbl_val(const Grille_val* g, FILE* fd) : dim(g->get_dim_tbl()), 
  gval(g) {
  
  fread_be(&etat, sizeof(int), 1, fd) ;		// etat
  
  // Le tableau
  if (etat == ETATQCQ) {
    t = new double[get_taille()] ;
    fread_be(t, sizeof(double), get_taille(), fd) ;	    // le tableau
    tzri = new double[get_taille_i(0)] ;
    fread_be(tzri, sizeof(double), get_taille_i(0), fd) ;
    if (get_ndim() > 1) {
      txti = new double[get_taille_i(1)] ;
      fread_be(txti, sizeof(double), get_taille_i(1), fd) ; }
    else txti = 0x0 ;
    if (get_ndim() > 2) {
      typi = new double[get_taille_i(2)] ;
      fread_be(typi, sizeof(double), get_taille_i(2), fd) ; }
    else typi = 0x0 ;
  }
  else{
    t = 0x0 ; 
    tzri = 0x0 ;
    txti = 0x0 ;
    typi = 0x0 ;
  }
}

			//-------------//
			// Destructeur //
			//-------------//

Tbl_val::~Tbl_val() {
  del_t() ;
}

			//-------------//
			// Affectation //
			//-------------//

// From Tbl_val
void Tbl_val::operator=(const Tbl_val& tx)
{
  // Protection
  assert( gval == tx.gval ) ;
  assert(tx.get_etat() != ETATNONDEF) ;
  
  int n = get_taille() ;
  int ndim = get_ndim() ;
  switch (tx.etat) {
  case ETATZERO:
    set_etat_zero() ;
    break ;
    
  case ETATQCQ:
    set_etat_qcq() ;
    for (int i=0 ; i<n ; i++) {
      t[i] = tx.t[i] ;
    }
    for (int i=0; i<get_taille_i(0); i++) tzri[i] = tx.tzri[i] ;
    if (ndim > 1) for(int i=0; i < get_taille_i(1); i++)
      txti[i] = tx.txti[i] ;
    if (ndim > 2) for(int i=0; i < get_taille_i(2); i++)
      typi[i] = tx.typi[i] ; 
    break ;
    
  default:
    cout << "Erreur bizarre !" << endl ;
    abort() ;
    break ;
  }
}

// From double
void Tbl_val::operator=(double a)
{
  if ( a == double(0) ) {
    set_etat_zero() ;
  }
  else {
    int n = get_taille() ;
    set_etat_qcq() ;
    for (int i=0 ; i<n ; i++) {
      t[i] = a ;
    }
    for (int i=0 ; i < get_taille_i(0) ; i++) 
      tzri[i] = a ;
    
    if (txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
      txti[i] = a ;
    
    if (typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
      typi[i] = a ;
  }
}

// From int
void Tbl_val::operator=(int m)
{
  if (m == 0) {
    set_etat_zero() ;
  }
  else {
    int n = get_taille() ;
    set_etat_qcq() ;
    for (int i=0 ; i<n ; i++) {
      t[i] = m ;
    }
    for (int i=0 ; i < get_taille_i(0) ; i++) 
      tzri[i] = m ;
    
    if (txti != 0x0) for (int i=0 ; i < get_taille_i(1) ; i++) 
      txti[i] = m ;
    
    if (typi != 0x0) for (int i=0 ; i < get_taille_i(2) ; i++) 
      typi[i] = m ;
    
  }
}

    
			//------------//
			// Sauvegarde //
			//------------//

// save in a file

void Tbl_val::sauve(FILE* fd) const {
  
  fwrite_be(&etat, sizeof(int), 1, fd) ;		    // etat
  if (etat == ETATQCQ) {
    fwrite_be(t, sizeof(double), get_taille(), fd) ;	    // le tableau
    fwrite_be(tzri, sizeof(double), get_taille_i(0), fd) ;
    if (get_ndim() > 1) 
      fwrite_be(txti, sizeof(double), get_taille_i(1), fd) ;
    if (get_ndim() > 2) 
      fwrite_be(typi, sizeof(double), get_taille_i(2), fd) ;
  }
}

//-----------------//
// Gestion memoire //
//-----------------//

// Destructeur logique
void Tbl_val::del_t() {
  if (t != 0x0) delete [] t ;
  t = 0x0 ;
  if (tzri != 0x0) delete [] tzri ;
  tzri = 0x0 ;
  if (txti != 0x0) delete [] txti ;
  txti = 0x0 ;
  if (typi != 0x0) delete [] typi ;
  typi = 0x0 ;
  etat = ETATNONDEF ;
}

// ETATZERO
void Tbl_val::set_etat_zero() {
  if (etat == ETATZERO) return ;
  del_t() ;
  etat = ETATZERO ;
}

// ETATNONDEF
void Tbl_val::set_etat_nondef() {
  if (etat == ETATNONDEF) return ;
  del_t() ;
  etat = ETATNONDEF ;
}

// ETATQCQ
void Tbl_val::set_etat_qcq() {
  if (etat == ETATQCQ) return ;
  
  // Protection
  assert( (etat == ETATZERO) || (etat == ETATNONDEF) ) ; // sinon...
  
  t = new double[get_taille()] ;
  tzri = new double[get_taille_i(0)] ;
  int ndim = get_ndim() ;
  if (ndim > 1) {txti = new double[get_taille_i(1)] ;}
  else txti = 0x0 ;    
  if (ndim > 2) {typi = new double[get_taille_i(2)] ;}
  else typi = 0x0 ;
  etat = ETATQCQ ;
}

// ZERO hard
void Tbl_val::annule_hard() {
  if (t == 0x0) {
    t = new double[get_taille()] ;
  }
  for (int i=0 ; i<get_taille() ; i++) {
    t[i] = 0. ;
  }
  tzri = new double[get_taille_i(0)] ;
  for (int i=0 ; i < get_taille_i(0) ; i++) 
    tzri[i] = 0 ;
  int ndim = get_ndim() ;
  if (ndim > 1) {
    txti = new double[get_taille_i(1)] ;
    for (int i=0 ; i < get_taille_i(1) ; i++) txti[i] = 0 ;
  }
  else txti = 0x0 ;    
  if (ndim > 2) {
    typi = new double[get_taille_i(2)] ;
    for (int i=0 ; i < get_taille_i(2) ; i++) typi[i] = 0 ;
  }
  else typi = 0x0 ;
  
  etat = ETATQCQ ;
}

void Tbl_val::append_array(double* t_in) {
    assert (t_in != 0x0) ;
    del_t() ;
    t = t_in ;
    etat = ETATQCQ ;
}

			//------------------------//
			//	Display		  //
			//------------------------//
			
//-----------			
// Operator<<
//-----------			

ostream& operator<<(ostream& o, const Tbl_val& t) {
    
  int ndim = t.get_ndim() ;
  o.precision(4);
  o.setf(ios::showpoint);
  o << "*** Tbl_val " << ndim << "D" << "   size: " ; 
  for (int i = 0; i<ndim-1; i++) {
    o << t.get_dim(ndim-1-i) ;
    if (ndim-i == 3) o << "(Y)" << " x " ;
    if (ndim-i == 2) o << "(X)" << " x " ;
  } 
  o << t.get_dim(0) << "(Z)" << " + " << t.gval->get_fantome() << 
    " hidden cells on each side =  " << t.get_taille() << endl ;
  
  if (t.get_etat() == ETATZERO) {
    o << "Identically ZERO" << endl ;
    return o ;
  }
  
  if (t.get_etat() == ETATNONDEF) {
    o << "UNDEFINED STATE" << endl ;
    return o ;
  }
  
  assert(t.etat == ETATQCQ) ;
  switch (ndim) {
    
  case 1 : {
    for (int i=0 ; i<t.get_dim(0) ; i++) {
      o << " " << t(i)  ;
    }
    o << endl ;
    break ;
  }
  
  
  case 2 : {
    for (int j=0 ; j<t.get_dim(1) ; j++) {
      o << " J_x " << j << " : " << endl ;
      for (int i=0 ; i<t.get_dim(0) ; i++) {
	o << " " << t(j, i)  ;
      }
      o << endl ;
    }
    o << endl ;
    break ;
  }
  
  case 3 : {
    for (int k=0 ; k<t.get_dim(2) ; k++) {
      o << " K_y = " << k << " : " << endl ;
      for (int j=0 ; j<t.get_dim(1) ; j++) {
	o << " J_x = " << j << " : "  ;
	for (int i=0 ; i<t.get_dim(0) ; i++) {
	  o << " " << t(k, j, i)  ;
	}
	o << endl ;
      }
      o << endl ;
    }
    o << endl ;
    break ;
  }
  
  default : {
    cout << "operator<< Tbl_val : unexpected dimension !" << endl ;
    cout << " ndim = " << ndim << endl ; 	
    abort() ;
    break ;
  }
  }
  return o ;
}

//---------------
// Affiche_seuil
//---------------

void Tbl_val::affiche_seuil(ostream& ost, int precis,  double seuil) const {

  int ndim = get_ndim() ;
  ost << "*** Tbl_val " << ndim << "D" << "   size: " ; 
  for (int i = 0; i<ndim-1; i++) {
    ost << get_dim(i) << " x " ;
  } 
  ost << get_dim(ndim-1) << " + " << gval->get_fantome() <<
    " hidden cells on each side  =  " << get_taille() << endl ; 
  
  // Cas particuliers
  //-----------------
  
  if (etat == ETATNONDEF) {
    ost << "    state: UNDEFINED" << endl ;
    return ;
  }
  
  if (etat == ETATZERO) {
    ost << "    state: ZERO" << endl ;
    return ;
  }
  
  // Affichage des elements du tableau 
  //----------------------------------
  
  ost << "    threshold for display : " << seuil << endl ; 
  ost.precision(precis);
  ost.setf(ios::showpoint);
  
  ost << "   Values on the nodes, without hidden cells:" << endl ;
  switch (get_ndim()) {
  case 1 : {			    	// cas 1-D
    
    for (int i=0; i<get_dim(0); i++) {
      ost <<  " " << setw(precis) << (*this)(i)  ;
    }
    ost << endl ;
    break ;
  }
  
  case 2 : {				// cas 2-D
    
    for (int j=0; j<get_dim(1); j++) {
      ost <<  " #j=" << j << " : "  ;
      for (int i=0; i<get_dim(0); i++){
	ost <<  " " << setw(precis) << (*this)(j, i)  ;
      }
      ost << endl;
    }
    ost << endl;
    break;
	}
  
  case 3 : {				// cas 3-D
    for (int k=0; k<get_dim(2); k++) {
      for (int j=0; j<get_dim(1); j++){
	int test_imp = 0 ;
	for (int i=0; i<get_dim(0); i++){
	  if ( fabs( (*this)(k, j, i) ) >= seuil ) 
	    test_imp = 1 ; 
	}
	if (test_imp == 1 ) {
	  ost <<  " #k=" << k <<",j=" << j << " : "  ;
	  for (int i=0; i<get_dim(0); i++){
	    ost <<  " " << setw(precis) << (*this)(k, j, i) ;
	  }
	  ost << endl ;
		    }
      }
    }
    ost << endl;
    break;
  }
  
  default : {
    cout << "Tbl_val:affiche_seuil : unexpected dimension !" << endl ;
    cout << " get_ndim() = " << get_ndim() << endl ; 	
    abort() ; 
    break;
  }               
  
  }     // fin du switch sur le nombre de dimensions
  
  // On restaure l'etat du flot ost a ses valeurs standards:
  ost.precision(6);
  ost.unsetf(ios::showpoint);
}



