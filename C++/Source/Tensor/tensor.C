/*
 *  Methods of class Tensor
 *
 *   (see file tensor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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


char tensor_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2003/09/24 15:10:54  j_novak
 * Suppression of the etat flag in class Tensor (still present in Scalar)
 *
 * Revision 1.2  2003/09/23 08:52:17  e_gourgoulhon
 * new version
 *
 * Revision 1.1  2003/09/22 12:52:51  e_gourgoulhon
 * First version: not ready yet !
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// Headers Lorene
#include "tensor.h"
#include "utilitaires.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
Tensor::Tensor(const Map& map, int val, const Itbl& tipe, 
		 const Base_vect& triad_i) 
		: mp(&map), valence(val), triad(&triad_i), type_indice(tipe), 
		   n_comp(int(pow(3., val))){
		
    // Des verifs :
    assert (valence >= 0) ;
    assert (tipe.get_ndim() == 1) ;
    assert (valence == tipe.get_dim(0)) ;
    for (int i=0 ; i<valence ; i++)
	assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
    
    cmp = new Scalar*[n_comp] ;

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(map) ;
	
}

// Standard constructor with the triad passed as a pointer
// -------------------------------------------------------
Tensor::Tensor(const Map& map, int val, const Itbl& tipe, 
		 const Base_vect* triad_i) 
		: mp(&map), valence(val), triad(triad_i), type_indice(tipe), 
		   n_comp(int(pow(3., val))){
		
    // Des verifs :
    assert (valence >= 0) ;
    assert (tipe.get_ndim() == 1) ;
    assert (valence == tipe.get_dim(0)) ;
    for (int i=0 ; i<valence ; i++)
	assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
    
    cmp = new Scalar*[n_comp] ;

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(map) ;
}




// Standard constructor when all the indices are of the same type
// --------------------------------------------------------------
Tensor::Tensor(const Map& map, int val, int tipe, const Base_vect& triad_i) 
		: mp(&map), valence(val), triad(&triad_i), type_indice(val), 
                  n_comp(int(pow(3., val))){
    
    // Des verifs :
    assert (valence >= 0) ;
    assert ((tipe == COV) || (tipe == CON)) ;
    type_indice.set_etat_qcq() ;
    type_indice = tipe ;
    
    cmp = new Scalar*[n_comp] ;

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(map) ;

}	
	
// Copy constructor
// ----------------
Tensor::Tensor (const Tensor& source) : 
    mp(source.mp), valence(source.valence), triad(source.triad), 
    type_indice(source.type_indice) {
  
    n_comp = int(pow(3., valence)) ;
        
    cmp = new Scalar*[n_comp] ;
    for (int i=0 ; i<n_comp ; i++) {
	int place_source = source.position(indices(i)) ;
	cmp[i] = new Scalar(*source.cmp[place_source]) ;
    }


}   


// Constructor from a file
// -----------------------
Tensor::Tensor(const Map& mapping, const Base_vect& triad_i, FILE* fd)
		 : mp(&mapping), triad(&triad_i), type_indice(fd){
   
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
    
    cmp = new Scalar*[n_comp] ;
    for (int i=0 ; i<n_comp ; i++)
      cmp[i] = new Scalar(*mp, *(mp->get_mg()), fd) ;

}


//  Constructor for a scalar field: to be used by the derived
//  class {\tt Scalar}
//-----------------------------------------------------------
Tensor::Tensor(const Map& map) : mp(&map), valence(0), triad(0x0),
		type_indice(0), n_comp(1) {
		
		cmp = new Scalar*[n_comp] ; 
		cmp[0] = 0x0 ; 
		
}


// Constructor used by the derived classes
// ---------------------------------------
Tensor::Tensor (const Map& map, int val, const Itbl& tipe, int compo, 
		const Base_vect& triad_i) :
     mp(&map), valence(val), triad(&triad_i), type_indice(tipe), n_comp(compo)
{
     
    // Des verifs :
    assert (valence >= 0) ;
    assert (tipe.get_ndim() == 1) ;
    assert (n_comp > 0) ;   
    assert (valence == tipe.get_dim(0)) ;
    for (int i=0 ; i<valence ; i++)
	assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
    
    cmp = new Scalar*[n_comp] ;

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(map) ;
    
}

// Constructor used by the derived classes when all the indices are of 
// the same type.
// -------------------------------------------------------------------
Tensor::Tensor (const Map& map, int val, int tipe, int compo, 
		const Base_vect& triad_i) :
     mp(&map), valence(val), triad(&triad_i), type_indice(val), n_comp(compo)
{

    // Des verifs :
    assert (valence >= 0) ;
    assert (n_comp >= 0) ;
    assert ((tipe == COV) || (tipe == CON)) ;
    type_indice.set_etat_qcq() ;
    type_indice = tipe ;
    
    cmp = new Scalar*[n_comp] ;

    for (int i=0 ; i<n_comp ; i++)
      cmp[i] = new Scalar(map) ;

}

			//--------------//
			//  Destructor  //
			//--------------//


Tensor::~Tensor () {
    
    del_deriv() ;

    for (int i=0 ; i<n_comp ; i++)
      delete cmp[i] ;
    delete [] cmp ;
}



void Tensor::del_deriv() {

}


void Tensor::set_etat_qcq() { 
    
    del_deriv() ;
    for (int i=0 ; i<n_comp ; i++) {
      cmp[i]->set_etat_qcq() ; 
    }
}

void Tensor::set_etat_nondef() { 
    
    del_deriv() ;
    for (int i=0 ; i<n_comp ; i++) {
      cmp[i]->set_etat_nondef() ; 
    }
}

void Tensor::set_etat_zero() { 
    
    del_deriv() ;
    for (int i=0 ; i<n_comp ; i++) {
      cmp[i]->set_etat_zero() ; 
    }
}


// Allocates everything
// --------------------
void Tensor::allocate_all() {
    
  del_deriv() ;
  for (int i=0 ; i<n_comp ; i++) {
    cmp[i]->allocate_all() ; 
  }
	
} 



void Tensor::change_triad(const Base_vect& bi) {
    
    // bi.change_basis(*this) ; 
    cout << "Tensor::change_triad not ready yet !" << endl ; 
    abort(); 
	
}

void Tensor::set_triad(const Base_vect& bi) {
    
    triad = &bi ; 
    
}

int Tensor::position (const Itbl& idx) const {
    
    assert (idx.get_ndim() == 1) ;
    assert (idx.get_dim(0) == valence) ;
    
    for (int i=0 ; i<valence ; i++)
	assert ((idx(i)>=1) && (idx(i)<=3)) ;
    int res = 0 ;
    for (int i=0 ; i<valence ; i++)
        res = 3*res+(idx(i)-1) ;
    
    return res;
}

Itbl Tensor::indices (int place) const {
    
    assert ((place >= 0) && (place < n_comp)) ;

    Itbl res(valence) ;
    res.set_etat_qcq() ;
    	    
    for (int i=valence-1 ; i>=0 ; i--) {
		res.set(i) = div(place, 3).rem ;
		place = int((place-res(i))/3) ;
		res.set(i)++ ; 
	}
    return res ;
}

void Tensor::operator=(const Tensor& t) {
    
    assert (valence == t.valence) ;

    del_deriv() ;

    triad = t.triad ; 

    for (int i=0 ; i<valence ; i++)
      assert(t.type_indice(i) == type_indice(i)) ;
	
    for (int i=0 ; i<n_comp ; i++) {
      int place_t = t.position(indices(i)) ;
      *cmp[i] = *t.cmp[place_t] ;
    }
}



// Affectation d'un tenseur d'ordre 2 :
Scalar& Tensor::set(int ind1, int ind2) {
    
    del_deriv() ;
    assert (valence == 2) ;
    
    Itbl ind (valence) ;
    ind.set_etat_qcq() ;
    ind.set(0) = ind1 ;
    ind.set(1) = ind2 ;
    
    int place = position(ind) ;
    
    return *cmp[place] ;
}

// Affectation d'un tenseur d'ordre 3 :
Scalar& Tensor::set(int ind1, int ind2, int ind3) {
    
    del_deriv() ;
    assert (valence == 3) ;
    
    Itbl indices(valence) ;
    indices.set_etat_qcq() ;
    indices.set(0) = ind1 ;
    indices.set(1) = ind2 ;
    indices.set(2) = ind3 ;
    int place = position(indices) ;
 
    return *cmp[place] ;
}

// Affectation cas general
Scalar& Tensor::set(const Itbl& indices) {
    
    assert (indices.get_ndim() == 1) ;
    assert (indices.get_dim(0) == valence) ;
    
    del_deriv() ;
	
    int place = position(indices) ;
    
    return *cmp[place] ;
}

// Annulation dans des domaines
void Tensor::annule(int l) {
    
    annule(l, l) ;     
}

void Tensor::annule(int l_min, int l_max) {
    
    // Cas particulier: annulation globale : 
    if ( (l_min == 0) && (l_max == mp->get_mg()->get_nzone()-1) ) {
      set_etat_zero() ;
      return ; 
    }
    
    // Annulation des composantes:
    for (int i=0 ; i<n_comp ; i++) {
      cmp[i]->annule(l_min, l_max) ; 
    }
	
    //## Annulation des membres derives
    //## pas avec un del_deriv() ;
    
}



const Scalar& Tensor::operator()(int indice1, int indice2) const {
    
    assert(valence == 2) ;
    
    Itbl idx(2) ;		
    idx.set_etat_qcq() ;	
    idx.set(0) = indice1 ;
    idx.set(1) = indice2 ;
    return *cmp[position(idx)] ;

}

const Scalar& Tensor::operator()(int indice1, int indice2, int indice3) const {
    
    assert(valence == 3) ;
    
    Itbl idx(3) ;		
    idx.set_etat_qcq() ;	
    idx.set(0) = indice1 ;
    idx.set(1) = indice2 ;
    idx.set(2) = indice3 ;
    return *cmp[position(idx)] ;
}


const Scalar& Tensor::operator()(const Itbl& ind) const {
    
    assert (ind.get_ndim() == 1) ;
    assert (ind.get_dim(0) == valence) ;
    return *cmp[position(ind)] ;
    
}


// Gestion de la CED :
void Tensor::dec_dzpuis() {
    
  for (int i=0 ; i<n_comp ; i++)
    cmp[i]->dec_dzpuis() ;
}

void Tensor::inc_dzpuis() {

    for (int i=0 ; i<n_comp ; i++)
      cmp[i]->inc_dzpuis() ;
}

void Tensor::dec2_dzpuis() {
    
  for (int i=0 ; i<n_comp ; i++)
    cmp[i]->dec2_dzpuis() ;
}

void Tensor::inc2_dzpuis() {
    
  for (int i=0 ; i<n_comp ; i++)
    cmp[i]->inc2_dzpuis() ;
}

void Tensor::mult_r_zec() {
    
  for (int i=0 ; i<n_comp ; i++) 
    cmp[i]->mult_r_zec() ;
}


// Le cout :
ostream& operator<<(ostream& flux, const Tensor &source ) {
    
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
    
    for (int i=0 ; i<source.n_comp ; i++) {

      Itbl num_indices (source.indices(i)) ;
      flux << "Component " ;
		
      if (source.valence != 0) {
	for (int j=0 ; j<source.valence ; j++)
	  flux << "  " << num_indices(j) ;
      }
      else
	flux << "  " << 0 ;
      flux << " : " << endl ;
      flux << "-------------" << endl ; 
      
      flux << *source.cmp[i] << endl ;
    }
    
    flux << " -----------------------------------------------------" << endl ;
    return flux ;
}


void Tensor::sauve(FILE* fd) const {
    
    type_indice.sauve(fd) ;	// type des composantes
    fwrite_be(&valence, sizeof(int), 1, fd) ;    // la valence
    
    if (valence != 0) {
		triad->sauve(fd) ;	    // Vectorial basis
    }
    
    fwrite_be(&n_comp, sizeof(int), 1, fd) ; // nbre composantes
    for (int i=0 ; i<n_comp ; i++)
      cmp[i]->sauve(fd) ;

}


