/*
 *  Methods of class Qtenseur. 
 *
 *   (see file qtenseur.h for documentation)
 *
 */

/*
 *   Copyright (c) 2002 Jerome Novak
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

char qtenseur_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/09/19 09:52:43  j_novak
 * Added objects Qtenseur and Qmetrique for 4D tensor and metric handling.
 *
 *
 * $Header$
 *
 */

// Headers C++
#include <iostream.h>

// Headers C
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// Headers Lorene
#include "qtenseur.h"
#include "qmetrique.h"
#include "utilitaires.h"

//--------------//
// Constructors //
//--------------//

// Constructor for a scalar field
// ------------------------------
Qtenseur::Qtenseur (const Map& map) : 
  mp(&map), valence(0), triad(0x0),
  type_indice(0), n_comp(1), etat(ETATNONDEF), met_depend(0x0) {
    
  c = new Cmp*[n_comp] ;
  c[0] = 0x0 ;
  set_der_0x0() ;
}



// Constructor for a scalar field and from a {\tt Cmp} 
// ---------------------------------------------------
Qtenseur::Qtenseur (const Cmp& ci) : 
  mp(ci.get_mp()), valence(0), triad(0x0),
  type_indice(0), n_comp(1), etat(ci.get_etat()), 
  met_depend(0x0) {
    
  assert(ci.get_etat() != ETATNONDEF) ; 
    
  c = new Cmp*[n_comp] ;
  set_der_0x0() ;

  if ( ci.get_etat() != ETATZERO ) {
    assert( ci.get_etat() == ETATQCQ ) ; 
    c[0] = new Cmp(ci) ;
  }
  else {
    c[0] = 0x0 ;
  }
}

// Standard constructor 
// --------------------
Qtenseur::Qtenseur(const Map& map, int val, const Itbl& tipe, 
		   const Base_vect& triad_i) 
  : mp(&map), valence(val), triad(&triad_i), type_indice(tipe), 
  n_comp(int(pow(4., val))), etat(ETATNONDEF) {
		
  // Des verifs :
  assert (valence >= 0) ;
  assert (tipe.get_ndim() == 1) ;
  assert (valence == tipe.get_dim(0)) ;
  for (int i=0 ; i<valence ; i++)
    assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
    
  c = new Cmp*[n_comp] ;
  for (int i=0 ; i<n_comp ; i++)
    c[i] = 0x0 ;
  set_der_0x0() ;
}

// Standard constructor with the triad passed as a pointer
// -------------------------------------------------------
Qtenseur::Qtenseur(const Map& map, int val, const Itbl& tipe, 
		   const Base_vect* triad_i) 
  : mp(&map), valence(val), triad(triad_i), type_indice(tipe), 
  n_comp(int(pow(4., val))), etat(ETATNONDEF),
  met_depend(0x0) {
		
  // Des verifs :
  assert (valence >= 0) ;
  assert (tipe.get_ndim() == 1) ;
  assert (valence == tipe.get_dim(0)) ;
  for (int i=0 ; i<valence ; i++)
    assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
    
  if (valence == 0) {	    // particular case of a scalar 
    triad = 0x0 ; 
  }   
    
  c = new Cmp*[n_comp] ;
  for (int i=0 ; i<n_comp ; i++)
    c[i] = 0x0 ;
  set_der_0x0() ;
}




// Standard constructor when all the indices are of the same type
// --------------------------------------------------------------
Qtenseur::Qtenseur(const Map& map, int val, int tipe, const Base_vect& triad_i)
  : mp(&map), valence(val), triad(&triad_i), type_indice(val), 
  n_comp(int(pow(4., val))), etat (ETATNONDEF), met_depend(0x0) {
    
  // Des verifs :
  assert (valence >= 0) ;
  assert ((tipe == COV) || (tipe == CON)) ;
  type_indice.set_etat_qcq() ;
  type_indice = tipe ;
    
  c = new Cmp*[n_comp] ;
  for (int i=0 ; i<n_comp ; i++)
    c[i] = 0x0 ;
  set_der_0x0() ;
}	
	
// Copy constructor
// ----------------
Qtenseur::Qtenseur (const Qtenseur& source) : 
  mp(source.mp), valence(source.valence), triad(source.triad), 
  type_indice(source.type_indice), etat (source.etat) {
  
  n_comp = int(pow(4., valence)) ;
        
  c = new Cmp*[n_comp] ;
  for (int i=0 ; i<n_comp ; i++) {
    int place_source = source.donne_place(donne_indices(i)) ;
    if (source.c[place_source] == 0x0)
      c[i] = 0x0 ;
    else
      c[i] = new Cmp(*source.c[place_source]) ;
  }
  set_der_0x0() ;
    
  if (source.p_gradient != 0x0)
    p_gradient = new Qtenseur (*source.p_gradient) ;
    
  if (source.p_gradient_spher != 0x0)
    p_gradient_spher = new Qtenseur (*source.p_gradient_spher) ;

  met_depend = source.met_depend ;
  if (met_depend != 0x0) {

    set_dependance (*met_depend) ;

    if (source.p_derive_cov != 0x0)
      p_derive_cov = new Qtenseur (*source.p_derive_cov) ;
    if (source.p_derive_con != 0x0)
      p_derive_con = new Qtenseur (*source.p_derive_con) ;
    if (source.p_carre_scal != 0x0)
      p_carre_scal = new Qtenseur (*source.p_carre_scal) ;
  }
}   

// Constructor from a symmetric tensor
// -----------------------------------
Qtenseur::Qtenseur (const Qtenseur_sym& source) :
  mp(source.mp), valence(source.valence), triad(source.triad), 
  type_indice(source.type_indice), etat(source.etat) {
    
  n_comp = int(pow(4., valence)) ;
        
  c = new Cmp*[n_comp] ;
  for (int i=0 ; i<n_comp ; i++) {
    int place_source = source.donne_place(donne_indices(i)) ;
    if (source.c[place_source] == 0x0)
      c[i] = 0x0 ;
    else
      c[i] = new Cmp(*source.c[place_source]) ;
  }
  set_der_0x0() ;
    
  if (source.p_gradient != 0x0)
    p_gradient = new Qtenseur_sym (*source.p_gradient) ;

  met_depend = source.met_depend ;
  if (met_depend != 0x0) {
	
    set_dependance (*met_depend) ;
      
    if (source.p_derive_cov != 0x0)
      p_derive_cov = new Qtenseur (*source.p_derive_cov) ;
    if (source.p_derive_con != 0x0)
      p_derive_con = new Qtenseur (*source.p_derive_con) ;
    if (source.p_carre_scal != 0x0)
      p_carre_scal = new Qtenseur (*source.p_carre_scal) ;
      
  }
}

// Constructor from a file
// -----------------------
Qtenseur::Qtenseur(const Map& mapping, const Base_vect& triad_i, FILE* fd)
  : mp(&mapping), triad(&triad_i), type_indice(fd) {
   
  fread_be(&valence, sizeof(int), 1, fd) ;

  if (valence != 0) {
    Base_vect* triad_fich = Base_vect::bvect_from_file(fd) ; 
    assert( *triad_fich == *triad) ; 
    delete triad_fich ; 
  }
  else{
    triad = 0x0 ; 
  }
    
  fread_be(&n_comp, sizeof(int), 1, fd) ;
  fread_be(&etat, sizeof(int), 1, fd) ;
    
  c = new Cmp*[n_comp] ;
  for (int i=0 ; i<n_comp ; i++)
    c[i] = 0x0 ;
  if (etat == ETATQCQ)
    for (int i=0 ; i<n_comp ; i++)
      c[i] = new Cmp(*mp, *mp->get_mg(), fd) ;

  set_der_0x0() ;

  met_depend = 0x0 ;
    
}


// Constructor from a file for a scalar field
// ------------------------------------------
Qtenseur::Qtenseur (const Map& mapping, FILE* fd) 
  : mp(&mapping), type_indice(fd) {
   
  fread_be(&valence, sizeof(int), 1, fd) ;

  assert(valence == 0) ; 
    
  triad = 0x0 ; 
    
  fread_be(&n_comp, sizeof(int), 1, fd) ;
    
  assert(n_comp == 1) ; 

  fread_be(&etat, sizeof(int), 1, fd) ;
    
  c = new Cmp*[n_comp] ;

  if (etat == ETATQCQ) {
    c[0] = new Cmp(*mp, *mp->get_mg(), fd) ;
  }
  else{
    c[0] = 0x0 ; 
  }
    
  set_der_0x0() ;

  met_depend = 0x0 ;
}




// Constructor used by the derived classes
// ---------------------------------------
Qtenseur::Qtenseur (const Map& map, int val, const Itbl& tipe, int compo, 
		    const Base_vect& triad_i) :
  mp(&map), valence(val), triad(&triad_i), type_indice(tipe), n_comp(compo),
  etat (ETATNONDEF), met_depend(0x0) {
     
  // Des verifs :
  assert (valence >= 0) ;
  assert (tipe.get_ndim() == 1) ;
  assert (n_comp > 0) ;   
  assert (valence == tipe.get_dim(0)) ;
  for (int i=0 ; i<valence ; i++)
    assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
    
  c = new Cmp*[n_comp] ;
  for (int i=0 ; i<n_comp ; i++)
    c[i] = 0x0 ;
  set_der_0x0() ;
}

// Constructor used by the derived classes when all the indices are of 
// the same type.
// -------------------------------------------------------------------
Qtenseur::Qtenseur (const Map& map, int val, int tipe, int compo, 
		    const Base_vect& triad_i) :
  mp(&map), valence(val), triad(&triad_i), type_indice(val), n_comp(compo), 
  etat (ETATNONDEF), met_depend(0x0) {
  // Des verifs :
  assert (valence >= 0) ;
  assert (n_comp >= 0) ;
  assert ((tipe == COV) || (tipe == CON)) ;
  type_indice.set_etat_qcq() ;
  type_indice = tipe ;
    
  c = new Cmp*[n_comp] ;
  for (int i=0 ; i<n_comp ; i++)
    c[i] = 0x0 ;
  set_der_0x0() ;
}

//--------------//
//  Destructor  //
//--------------//


Qtenseur::~Qtenseur () {
    
  del_t() ;
  delete [] c ;
}



void Qtenseur::del_t() {
  del_derive() ;
  for (int i=0 ; i<n_comp ; i++)
    if (c[i] != 0x0) {
      delete c[i] ;
      c[i] = 0x0 ;
    }
}

void Qtenseur::del_derive_met() const {
  // On gere la metrique ...
  if (met_depend != 0x0) {
    for (int i=0 ; i<N_DEPEND ; i++)
      if (met_depend->dependances[i] == this)
	met_depend->dependances[i] = 0x0 ;
    if (p_derive_cov != 0x0)
      delete p_derive_cov ;
    if (p_derive_con != 0x0)
      delete p_derive_con ;
    if (p_carre_scal != 0x0)
      delete p_carre_scal ;
    set_der_met_0x0() ;
  }
}


void Qtenseur::del_derive () const {
  del_derive_met() ;
  if (p_gradient != 0x0)
    delete p_gradient ;
  if (p_gradient_spher != 0x0)
    delete p_gradient_spher ;
  set_der_0x0() ;
}

void Qtenseur::set_der_met_0x0() const {
  p_derive_cov = 0x0 ;
  p_derive_con = 0x0 ;
  p_carre_scal = 0x0 ;
}


void Qtenseur::set_der_0x0() const {
  set_der_met_0x0() ;
  p_gradient = 0x0 ;   
  p_gradient_spher = 0x0 ;   
}

void Qtenseur::set_dependance (const Qmetrique& met) const {
        
  // Cas ou on a calcule des trucs avec une autre metrique ...
  if ((met_depend != 0x0) && (&met != met_depend)) {
    del_derive_met() ;
    met_depend = 0x0 ;
  }
  
  if (met_depend == 0x0) {
    int conte = 0 ;
    while ((conte < N_DEPEND) && (met.dependances[conte] != 0x0))
      conte ++ ; 
    
    if (conte == N_DEPEND) {
      cout << "Too many dependancies in Qtenseur::set_dependances " << endl ;
      abort() ;
    }
    else {
      met.dependances[conte] = this ;
      met_depend = &met ;
    }
  }
}

void Qtenseur::set_etat_qcq() { 
    
  del_derive() ;
  for (int i=0 ; i<n_comp ; i++)
    if (c[i] == 0x0)
      c[i] = new Cmp(mp) ;
  etat = ETATQCQ ;
}

void Qtenseur::set_etat_zero() { 
  del_t() ;
  etat = ETATZERO ;
}

void Qtenseur::set_etat_nondef() { 
  del_t() ;
  etat = ETATNONDEF ;
}

// Allocates everything
// --------------------
void Qtenseur::allocate_all() {
    
  set_etat_qcq() ; 
  for (int i=0 ; i<n_comp ; i++) {
    c[i]->allocate_all() ; 
  }
	
} 



//  void Qtenseur::change_triad(const Base_vect& bi) {
    
//    bi.change_basis(*this) ; 
    
//  }

void Qtenseur::set_triad(const Base_vect& bi) {
    
  triad = &bi ; 
    
}

int Qtenseur::donne_place (const Itbl& idx) const {
    
  assert (idx.get_ndim() == 1) ;
  assert (idx.get_dim(0) == valence) ;
    
  for (int i=0 ; i<valence ; i++)
    assert ((idx(i)>=0) && (idx(i)<4)) ;
  int res = 0 ;
  for (int i=0 ; i<valence ; i++)
    res = 4*res+idx(i) ;
    
  return res;
}

Itbl Qtenseur::donne_indices (int place) const {
    
  assert ((place >= 0) && (place < n_comp)) ;

  Itbl res(valence) ;
  res.set_etat_qcq() ;
    	    
  for (int i=valence-1 ; i>=0 ; i--) {
    res.set(i) = div(place, 4).rem ;
    place = int((place-res(i))/4) ;
  }
  return res ;
}

void Qtenseur::operator=(const Qtenseur& t) {
    
  assert (valence == t.valence) ;

  triad = t.triad ; 

  for (int i=0 ; i<valence ; i++)
    assert (t.type_indice(i) == type_indice(i)) ;
	
  switch (t.etat) {
  case ETATNONDEF: {
    set_etat_nondef() ;
    break ;
  }
	
  case ETATZERO: {
    set_etat_zero() ;
    break ;
  }
	
  case ETATQCQ: {
    set_etat_qcq() ;
    for (int i=0 ; i<n_comp ; i++) {
      int place_t = t.donne_place(donne_indices(i)) ;
      if (t.c[place_t] == 0x0)
	c[i] = 0x0 ;
      else 
	*c[i] = *t.c[place_t] ;
    }
    break ;
  }
	
  default: {
    cout << "Unknown state in Qtenseur::operator= " << endl ;
    abort() ;
    break ;
  }
  }
}


void Qtenseur::operator=(const Cmp& ci) {
    
  assert (valence == 0) ;

  switch (ci.get_etat()) {
  case ETATNONDEF: {
    set_etat_nondef() ;
    break ;
  }
	
  case ETATZERO: {
    set_etat_zero() ;
    break ;
  }
	
  case ETATQCQ: {
    set_etat_qcq() ;
    *(c[0]) = ci ; 
    break ;
  }
	
  default: {
    cout << "Unknown state in Qtenseur::operator= " << endl ;
    abort() ;
    break ;
  }
  }
}

void Qtenseur::operator=(double x) {

  if (x == double(0)) {
    set_etat_zero() ;
  }
  else {
    assert(valence == 0) ; 
    set_etat_qcq() ;
    *(c[0]) = x ; 
  }

}

void Qtenseur::operator=(int x) {

  if (x == 0) {
    set_etat_zero() ;
  }
  else {
    assert(valence == 0) ; 
    set_etat_qcq() ;
    *(c[0]) = x ; 
  }

}


// Affectation d'un scalaire ...
Cmp& Qtenseur::set () {
    
  del_derive() ;
  assert(etat == ETATQCQ) ;
  assert (valence == 0) ;  
  return *c[0] ;
}


// Affectation d'un vecteur :
Cmp& Qtenseur::set (int ind) {
    
  del_derive() ;
  assert(valence == 1) ;
  assert (etat == ETATQCQ) ;
  assert ((ind >= 0) && (ind < 4)) ;
    
  return *c[ind] ;
}

// Affectation d'un tenseur d'ordre 2 :
Cmp& Qtenseur::set (int ind1, int ind2) {
    
  del_derive() ;
  assert (valence == 2) ;
  assert (etat == ETATQCQ) ;
  assert ((ind1 >= 0) && (ind1 < 4)) ;
  assert ((ind2 >= 0) && (ind2 < 4)) ;
    
  Itbl ind (valence) ;
  ind.set_etat_qcq() ;
  ind.set(0) = ind1 ;
  ind.set(1) = ind2 ;
    
  int place = donne_place(ind) ;
    
  return *c[place] ;
}

// Affectation d'un tenseur d'ordre 3 :
Cmp& Qtenseur::set (int ind1, int ind2, int ind3) {
    
  del_derive() ;
  assert (valence == 3) ;
  assert (etat == ETATQCQ) ;
  assert ((ind1 >= 0) && (ind1 < 4)) ;
  assert ((ind2 >= 0) && (ind2 < 4)) ;
  assert ((ind3 >= 0) && (ind3 < 4)) ;
    
  Itbl indices(valence) ;
  indices.set_etat_qcq() ;
  indices.set(0) = ind1 ;
  indices.set(1) = ind2 ;
  indices.set(2) = ind3 ;
  int place = donne_place(indices) ;
 
  return *c[place] ;
}

// Affectation cas general
Cmp& Qtenseur::set(const Itbl& indices) {
    
  assert (indices.get_ndim() == 1) ;
  assert (indices.get_dim(0) == valence) ;
    
  del_derive() ;
  assert (etat == ETATQCQ) ;
  for (int i=0 ; i<valence ; i++)
    assert ((indices(i)>=0) && (indices(i)<4)) ;
	
  int place = donne_place(indices) ;
    
  return *c[place] ;
}

// Annulation dans des domaines
void Qtenseur::annule(int l) {
    
  annule(l, l) ;     
}

void Qtenseur::annule(int l_min, int l_max) {
    
  // Cas particulier: annulation globale : 
  if ( (l_min == 0) && (l_max == mp->get_mg()->get_nzone()-1) ) {
    set_etat_zero() ;
    return ; 
  }
    
  assert( etat != ETATNONDEF ) ; 
    
  if ( etat == ETATZERO ) {
    return ;		// rien n'a faire si c'est deja zero
  }
  else {
    assert( etat == ETATQCQ ) ;	// sinon...
	
    // Annulation des composantes:
    for (int i=0 ; i<n_comp ; i++) {
      c[i]->annule(l_min, l_max) ; 
    }
	
    // Annulation des membres derives
    if (p_gradient != 0x0) p_gradient->annule(l_min, l_max) ;
    if (p_gradient_spher != 0x0) p_gradient_spher->annule(l_min, l_max) ;
    if (p_derive_cov != 0x0) p_derive_cov->annule(l_min, l_max) ;
    if (p_derive_con != 0x0) p_derive_con->annule(l_min, l_max) ;
    if (p_carre_scal != 0x0) p_carre_scal->annule(l_min, l_max) ;

  }
    
}




// Exctraction :
const Cmp& Qtenseur::operator()() const {
    
  assert(valence == 0) ;
    
  if (etat == ETATQCQ) return *c[0] ;	    // pour la performance,
  // ce cas est traite en premier,
  // en dehors du switch
  switch (etat) {
	
  case ETATNONDEF : {
    cout << "Undefined 4-Tensor in Qtenseur::operator() ..." << endl ;
    abort() ;
    return *c[0] ;  // bidon pour satisfaire le compilateur
  }
	    
  case ETATZERO : {
    return mp->cmp_zero() ;
  }
	    
	    
  default : {
    cout <<"Unknown state in Qtenseur::operator()" << endl ;
    abort() ;
    return *c[0] ;  // bidon pour satisfaire le compilateur
  }
  }
}


const Cmp& Qtenseur::operator() (int indice) const {
    
  assert ((indice>=0) && (indice<4)) ;
  assert(valence == 1) ;
    
  if (etat == ETATQCQ) return *c[indice] ;	 // pour la performance,
  // ce cas est traite en premier,
  // en dehors du switch
  switch (etat) {
	
  case ETATNONDEF : {
    cout << "Undefined 4-Tensor in Qtenseur::operator(int) ..." << endl ;
    abort() ;
    return *c[0] ;  // bidon pour satisfaire le compilateur
  }
	    
  case ETATZERO : {
    return mp->cmp_zero() ;
  }
	    	    
  default : {
    cout <<"Unknown state in Qtenseur::operator(int)" << endl ;
    abort() ;
    return *c[0] ;  // bidon pour satisfaire le compilateur
  }
  }
}

const Cmp& Qtenseur::operator() (int indice1, int indice2) const {
    
  assert ((indice1>=0) && (indice1<4)) ;
  assert ((indice2>=0) && (indice2<4)) ;
  assert(valence == 2) ;
    
  if (etat == ETATQCQ) {		// pour la performance,
    Itbl idx(2) ;		// ce cas est traite en premier,
    idx.set_etat_qcq() ;	// en dehors du switch
    idx.set(0) = indice1 ;
    idx.set(1) = indice2 ;
    return *c[donne_place(idx)] ;	
  }

  switch (etat) {
	
  case ETATNONDEF : {
    cout << "Undefined 4-Tensor in Qtenseur::operator(int, int) ..." << endl ;
    abort() ;
    return *c[0] ;  // bidon pour satisfaire le compilateur
  }
	    
  case ETATZERO : {
    return mp->cmp_zero() ;
  }
	   	    
  default : {
    cout <<"Unknown state in Qtenseur::operator(int, int)" << endl ;
    abort() ;
    return *c[0] ;  // bidon pour satisfaire le compilateur
  }
  } 
}

const Cmp& Qtenseur::operator() (int indice1, int indice2, int indice3) const {
    
  assert ((indice1>=0) && (indice1<4)) ;
  assert ((indice2>=0) && (indice2<4)) ;
  assert ((indice3>=0) && (indice3<4)) ;
  assert(valence == 3) ;
    
  if (etat == ETATQCQ) {		// pour la performance,
    Itbl idx(3) ;		// ce cas est traite en premier,
    idx.set_etat_qcq() ;	// en dehors du switch
    idx.set(0) = indice1 ;
    idx.set(1) = indice2 ;
    idx.set(2) = indice3 ;
    return *c[donne_place(idx)] ;
  }

  switch (etat) {
	
  case ETATNONDEF : {
    cout << "Undefined 4-Tensor in Qtenseur::operator(int, int, int) ..." << endl ;
    abort() ;
    return *c[0] ;  // bidon pour satisfaire le compilateur
  }
	    
  case ETATZERO : {
    return mp->cmp_zero() ;
  }
	    	    
  default : {
    cout <<"Unknown state in Qtenseur::operator(int, int, int)" << endl ;
    abort() ;
    return *c[0] ;  // bidon pour satisfaire le compilateur
  }
  }
}


const Cmp& Qtenseur::operator() (const Itbl& indices) const {
    
  assert (indices.get_ndim() == 1) ;
  assert (indices.get_dim(0) == valence) ;
  for (int i=0 ; i<valence ; i++)
    assert ((indices(i)>=0) && (indices(i)<4)) ;
    
  if (etat == ETATQCQ) {		    // pour la performance,
    return *c[donne_place(indices)]	;   // ce cas est traite en premier,
  }					    // en dehors du switch

  switch (etat) {
	
  case ETATNONDEF : {
    cout << "Undefined 4-Tensor in Qtenseur::operator(const Itbl&) ..." << endl ;
    abort() ;
    return *c[0] ;  // bidon pour satisfaire le compilateur
  }
	    
  case ETATZERO : {
    return mp->cmp_zero() ;
  }
	    
  default : {
    cout <<"Unknown state in Qtenseur::operator(const Itbl& )" << endl ;
    abort() ;
    return *c[0] ;  // bidon pour satisfaire le compilateur
  }
  }
    
}

// Gestion de la ZEC :
void Qtenseur::dec_dzpuis() {
    
  if (etat == ETATZERO) {
    return ; 
  }
    
  assert(etat == ETATQCQ) ;
   
  for (int i=0 ; i<n_comp ; i++)
    if (c[i] != 0x0)
      c[i]->dec_dzpuis() ;
}

void Qtenseur::inc_dzpuis() {
    
  if (etat == ETATZERO) {
    return ; 
  }

  assert(etat == ETATQCQ) ;
    
  for (int i=0 ; i<n_comp ; i++)
    if (c[i] != 0x0)
      c[i]->inc_dzpuis() ;
}

void Qtenseur::dec2_dzpuis() {
    
  if (etat == ETATZERO) {
    return ; 
  }
    
  assert(etat == ETATQCQ) ;
   
  for (int i=0 ; i<n_comp ; i++)
    if (c[i] != 0x0)
      c[i]->dec2_dzpuis() ;
}

void Qtenseur::inc2_dzpuis() {
    
  if (etat == ETATZERO) {
    return ; 
  }

  assert(etat == ETATQCQ) ;
    
  for (int i=0 ; i<n_comp ; i++)
    if (c[i] != 0x0)
      c[i]->inc2_dzpuis() ;
}

void Qtenseur::mult_r_zec() {
    
  if (etat == ETATZERO) {
    return ; 
  }

  assert(etat == ETATQCQ) ;
    
  for (int i=0 ; i<n_comp ; i++) 
    if (c[i] != 0x0)
      c[i]->mult_r_zec() ;
}

// Gestion des bases spectrales (valence <= 2)
void Qtenseur::set_std_base() {
    
  if (etat == ETATZERO) {
    return ; 
  }
    
  assert(etat == ETATQCQ) ;
  switch (valence) {
	
  case 0 : {
    c[0]->std_base_scal() ;
    break ;
  }
	    
  case 1 : {

    c[0]->std_base_scal() ;

    if ( triad->identify() == (mp->get_bvect_cart()).identify() ) {

      Base_val** bases = mp->get_mg()->std_base_vect_cart() ;

      for (int i=1 ; i<4 ; i++)
	(c[i]->va).set_base( *bases[i-1] ) ;
      for (int i=0 ; i<3 ; i++)
	delete bases[i] ;
      delete [] bases ;
    }
    else {
      assert( triad->identify() == (mp->get_bvect_spher()).identify()) ;
      Base_val** bases = mp->get_mg()->std_base_vect_spher() ;

      for (int i=1 ; i<4 ; i++)
	(c[i]->va).set_base( *bases[i-1] ) ;
      for (int i=0 ; i<3 ; i++)
	delete bases[i] ;
      delete [] bases ;
    }
    break ;
	    
  }
	    
  case 2 : { 

    if( triad->identify() == (mp->get_bvect_cart()).identify() ) {

      Base_val** bases = new Base_val*[4] ;
      Base_val** base_spa = mp->get_mg()->std_base_vect_cart() ;
      Base_val scal(mp->get_mg()->std_base_scal()) ;
      bases[0] = &scal ;
      for (int i=0; i<3; i++) 
	bases[i+1] = base_spa[i] ;
	    
      Itbl indices (2) ;
      indices.set_etat_qcq() ;
      for (int i=0 ; i<n_comp ; i++) {   
	indices = donne_indices(i) ;
	(c[i]->va).set_base( (*bases[indices(0)]) * 
			     (*bases[indices(1)]) ) ;
      }
      for (int i=0 ; i<3 ; i++) 
	delete bases[i+1] ;

      delete [] bases ;
      delete [] base_spa ;
    }
    else {
      assert( triad->identify() == (mp->get_bvect_spher()).identify()) ;
      Base_val** bases = new Base_val*[4] ;
      Base_val** base_spa = mp->get_mg()->std_base_vect_spher() ;
      Base_val scal(mp->get_mg()->std_base_scal()) ;
      bases[0] = &scal ;
      for (int i=0; i<3; i++) 
	bases[i+1] = base_spa[i] ;
	    
      Itbl indices (2) ;
      indices.set_etat_qcq() ;
      for (int i=0 ; i<n_comp ; i++) {   
	indices = donne_indices(i) ;
	(c[i]->va).set_base( (*bases[indices(0)]) * 
			     (*bases[indices(1)]) ) ;
      }
      for (int i=0 ; i<3 ; i++) 
	delete bases[i+1] ;

      delete [] bases ;
      delete [] base_spa ;
    }
    break ;
  }
	   
  default : {
    cout << 
      "Qtenseur::set_std_base() : the case valence = " << valence
	 << " is not treated !" << endl ;
    abort() ;
    break ;
  }
  }
}

// Le cout :
ostream& operator<<(ostream& flux, const Qtenseur &source ) {
    
  flux << "4D Tenseur" << endl ;
  flux << "Valence : " << source.valence << endl ;
  if (source.get_triad() != 0x0) {
    flux << "Vectorial basis (triad) on which the components are defined :" 
	 << endl ; 
    flux << *(source.get_triad()) << endl ;
  }
    
  if (source.valence != 0)
    flux << "Type of the indices : " << endl ;
  for (int i=0 ; i<source.valence ; i++) {
    flux << "Index " << i << " : " ;
    if (source.type_indice(i) == CON)
      flux << " contravariant." << endl ;
    else
      flux << " covariant." << endl ;
  }
    
  switch (source.etat) {
	
  case ETATZERO : {
    flux << "Null Qtenseur. " << endl ;
    break ;
  }
	
  case ETATNONDEF : {
    flux << "Undefined Qtenseur. " << endl ;
    break ;
  }
	
  case ETATQCQ : {
    for (int i=0 ; i<source.n_comp ; i++) {

      Itbl num_indices (source.donne_indices(i)) ;
      flux << "Component " ;
		
      if (source.valence != 0) {
	for (int j=0 ; j<source.valence ; j++)
	  flux << "  " << num_indices(j) ;
      }
      else
	flux << "  " << 0 ;
      flux << " : " << endl ;
      flux << "-------------" << endl ; 


      if (source.c[i] != 0x0)
	flux << *source.c[i] << endl ;
      else
	flux << "Unknown component ... " << endl ;
		    
    }
    break ;
  }
  default : {
    cout << "Unknown case in operator<< (ostream&, const Qtenseur&)" << endl ;
    abort() ;
    break ;
  }
  }
    
  flux << " -----------------------------------------------------" << endl ;
  return flux ;
}

void Qtenseur::sauve(FILE* fd) const {
    
  type_indice.sauve(fd) ;	// type des composantes
  fwrite_be(&valence, sizeof(int), 1, fd) ;    // la valence
    
  if (valence != 0) {
    triad->sauve(fd) ;	    // Vectorial basis
  }
    
  fwrite_be(&n_comp, sizeof(int), 1, fd) ; // nbre composantes
  fwrite_be(&etat, sizeof(int), 1, fd) ; // etat
   
  if (etat == ETATQCQ)
    for (int i=0 ; i<n_comp ; i++)
      c[i]->sauve(fd) ;
}


void Qtenseur::fait_gradient (const Param& para) const {
    
  assert (etat != ETATNONDEF) ;
    
  if (p_gradient != 0x0)
    return ;
  else {
 
    // Construction du resultat :
    Itbl tipe (valence+1) ;
    tipe.set_etat_qcq() ;
    tipe.set(0) = COV ;
    for (int i=0 ; i<valence ; i++)
      tipe.set(i+1) = type_indice(i) ;
	
    // Vectorial basis
    // ---------------
    if (valence != 0) {
      assert (*triad == mp->get_bvect_cart()) ;
    }

    p_gradient = new Qtenseur(*mp, valence+1, tipe, mp->get_bvect_cart()) ;
    p_gradient->set_etat_qcq() ;
    assert (para.get_n_double() == 1) ;
    double dt = para.get_double() ;
    Qtenseur dtemp(*mp, valence, type_indice, triad) ;
    switch (para.get_n_qtenseur_mod()) {
    case 0 : {
      dtemp.set_etat_zero() ;
      break ;
    }
    case 1 : {
      dtemp.set_etat_qcq() ;
      dtemp = (*this - para.get_qtenseur_mod(0))/dt ;
      break ;
    }
    case 2 : {
      dtemp.set_etat_qcq() ;
      dtemp = (3*(*this) - 4*para.get_qtenseur_mod(0) 
	       + para.get_qtenseur_mod(1))/(2*dt) ;
      break ;
    }
    default : {
      cout << "Problem in Qtenseur::gradient!" << endl ;
      abort() ;
      break ;
    }
    }

    // Boucle sur les indices :
	
    Itbl indices_source(valence) ;
    indices_source.set_etat_qcq() ;
    
    for (int j=0 ; j<p_gradient->n_comp ; j++) {    
      
      Itbl indices_res(p_gradient->donne_indices(j)) ;
      int i0 = indices_res(0) ;
      for (int m=0 ; m<valence ; m++)
	indices_source.set(m) = indices_res(m+1) ;
      if (i0 == 0) 
	p_gradient->set(indices_res) = dtemp(indices_source) ;
      else {
	p_gradient->set(indices_res) = 
	  (*this)(indices_source).deriv(indices_res(0)-1) ;
	p_gradient->set(indices_res).dec2_dzpuis() ;
      }
    }
  }
}

void Qtenseur::fait_gradient_spher (const Param& para) const {
    
  assert (etat != ETATNONDEF) ;
    
  if (p_gradient_spher != 0x0)
    return ;
  else {
    
    // Construction du resultat :
    
    if (valence != 0) {
      cout << 
	"Qtenseur::fait_gradient_spher : the valence must be zero !" 
	   << endl ; 
      abort() ; 
    }
    
    p_gradient_spher = new Qtenseur(*mp, 1, COV, mp->get_bvect_spher() ) ;
    
    p_gradient_spher->set_etat_qcq() ;
    
    p_gradient_spher->set(1) = c[0]->dsdr() ;	    // d/dr 
    p_gradient_spher->set(1).dec2_dzpuis() ;
    p_gradient_spher->set(2) = c[0]->srdsdt() ;	    // 1/r d/dtheta
    p_gradient_spher->set(2).dec2_dzpuis() ;
    p_gradient_spher->set(3) = c[0]->srstdsdp() ;   // 1/(r sin(theta))d/dphi
    p_gradient_spher->set(3).dec2_dzpuis() ;
    assert (para.get_n_double() == 1) ;
    double dt = para.get_double() ;
    Qtenseur dtemp(*mp) ;
    switch (para.get_n_qtenseur_mod()) {
    case 0 : {
      dtemp.set_etat_zero() ;
      break ;
    }
    case 1 : {
      dtemp.set_etat_qcq() ;
      dtemp = (*this - para.get_qtenseur_mod(0))/dt ;
      break ;
    }
    case 2 : {
      dtemp.set_etat_qcq() ;
      dtemp = (3*(*this) - 4*para.get_qtenseur_mod(0) 
	       + para.get_qtenseur_mod(1))/(2*dt) ;
      break ;
    }
    default : {
      cout << "Problem in Qtenseur::gradient_spher!" << endl ;
      abort() ;
      break ;
    }
    }
    p_gradient_spher->set(0) = dtemp() ;
    
  }
}


void Qtenseur::fait_derive_cov (const Param& para, const Qmetrique& metre) const {
    
  assert (etat != ETATNONDEF) ;
  assert (valence != 0) ;
  
  if (p_derive_cov != 0x0)
    return ;
  else {
    assert (para.get_n_double() == 1) ;
    double dt = para.get_double() ;
    p_derive_cov = new Qtenseur (gradient(para)) ;
    
    if ((valence != 0) && (etat != ETATZERO)) {
      
      assert( metre.gamma(dt).get_triad() == triad ) ; 
      Qtenseur* auxi ;
      for (int i=0 ; i<valence ; i++) {
	
	if (type_indice(i) == COV) {
	  auxi = new Qtenseur(contract(metre.gamma(dt), 0,(*this), i)) ;
	  
	  Itbl indices_gamma(p_derive_cov->valence) ;
	  indices_gamma.set_etat_qcq() ;
	  //On range comme il faut :
	  for (int j=0 ; j<p_derive_cov->n_comp ; j++) {
	    
	    Itbl indices (p_derive_cov->donne_indices(j)) ;
	    indices_gamma.set(0) = indices(0) ;
	    indices_gamma.set(1) = indices(i+1) ;
	    for (int idx=2 ; idx<p_derive_cov->valence ; idx++)
	      if (idx<=i+1)
		indices_gamma.set(idx) = indices(idx-1) ;
	      else
		indices_gamma.set(idx) = indices(idx) ;
	    
	    p_derive_cov->set(indices) -= (*auxi)(indices_gamma) ;
	  }
	}   
	else {
	  auxi = new Qtenseur(contract(metre.gamma(dt), 1, (*this), i)) ;
	  
	  Itbl indices_gamma(p_derive_cov->valence) ;
	  indices_gamma.set_etat_qcq() ;
	  
	  //On range comme il faut :
	  for (int j=0 ; j<p_derive_cov->n_comp ; j++) {
	    
	    Itbl indices (p_derive_cov->donne_indices(j)) ;
	    indices_gamma.set(0) = indices(i+1) ;
	    indices_gamma.set(1) = indices(0) ;
	    for (int idx=2 ; idx<p_derive_cov->valence ; idx++)
	      if (idx<=i+1)
		indices_gamma.set(idx) = indices(idx-1) ;
	      else
		indices_gamma.set(idx) = indices(idx) ;
	    p_derive_cov->set(indices) += (*auxi)(indices_gamma) ;
	  }
	}
	delete auxi ;
      }
    }
  }
}



void Qtenseur::fait_derive_con (const Param& para, const Qmetrique& metre) const {
    
  if (p_derive_con != 0x0)
    return ;
  else {
    // On calcul la derivee covariante :
    if (valence != 0)
      p_derive_con = new Qtenseur
	(contract(metre.con(), 1, derive_cov(para, metre), 0)) ;
	
    else
      p_derive_con = new Qtenseur
	(contract(metre.con(), 1, gradient(para), 0)) ;
  }
}

void Qtenseur::fait_carre_scal (const Qmetrique& met) const {
    
  if (p_carre_scal != 0x0)
    return ;
  else {
    assert (valence != 0) ;   // A ne pas appeler sur un scalaire ;
       
    // On bouge tous les indices :
    Qtenseur op_t(manipule(*this, met)) ;
   
    Qtenseur* auxi = new Qtenseur(contract(*this, 0, op_t, 0)) ;
    Qtenseur* auxi_old ;
    
    // On contracte tous les indices restant :
    for (int indice=1 ; indice<valence ; indice++) {
      auxi_old = new Qtenseur(contract(*auxi, 0, valence-indice)) ;
      delete auxi ;
      auxi = new Qtenseur(*auxi_old) ;
      delete auxi_old ;
    }
    p_carre_scal = new Qtenseur (*auxi) ;
    delete auxi ;
  }
}
    
const Qtenseur& Qtenseur::gradient (const Param& para) const {
  if (p_gradient == 0x0)
    fait_gradient(para) ;
  return *p_gradient ;
}

const Qtenseur& Qtenseur::gradient_spher(const Param& para) const {
  if (p_gradient_spher == 0x0)
    fait_gradient_spher(para) ;
  return *p_gradient_spher ;
}

const Qtenseur& Qtenseur::derive_cov (const Param& para, const Qmetrique& metre) const {
    
  if (valence == 0)
    return gradient(para) ;
  else {
    set_dependance(metre) ;
    if (p_derive_cov == 0x0)
      fait_derive_cov (para, metre) ;
    return *p_derive_cov ;
  }
}

const Qtenseur& Qtenseur::derive_con (const Param& para, const Qmetrique& metre) const {
  set_dependance(metre) ;
  if (p_derive_con == 0x0)
    fait_derive_con (para, metre) ;
  return *p_derive_con ;
}

const Qtenseur& Qtenseur::carre_scal (const Qmetrique& metre) const {
  set_dependance(metre) ;
  if (p_carre_scal == 0x0)
    fait_carre_scal (metre) ;
  return *p_carre_scal ;
}
