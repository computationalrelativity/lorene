/*
 *  Methods of class Tensor
 *
 *   (see file tensor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Philippe Grandclement (for preceding class Tenseur)
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
 * Revision 1.31  2004/02/05 15:03:47  e_gourgoulhon
 * Corrected bug in method derive_con().
 *
 * Revision 1.30  2004/01/27 13:05:11  j_novak
 * Removed the method Tensor::mult_r_ced()
 *
 * Revision 1.29  2004/01/19 16:32:13  e_gourgoulhon
 * Added operator()(int, int, int, int) and set(int, int, int, int)
 * for direct access to components of valence 4 tensors.
 *
 * Revision 1.28  2004/01/04 20:55:23  e_gourgoulhon
 * Method spectral_display(): added printing of type of class (through typeid).
 *
 * Revision 1.27  2003/12/30 23:09:47  e_gourgoulhon
 * Change in methods derive_cov() and divergence() to take into account
 *  the change of name: Metric::get_connect() --> Metric::connect().
 *
 * Revision 1.26  2003/12/27 15:01:50  e_gourgoulhon
 * Method derive_con(): the position of the derivation index has
 * been changed from the first one to the last one.
 *
 * Revision 1.25  2003/11/03 22:34:41  e_gourgoulhon
 * Method dec_dzpuis: changed the name of argument dec --> decrem
 * (in order not to shadow some globally defined dec).
 *
 * Revision 1.24  2003/10/29 11:02:54  e_gourgoulhon
 * Functions dec_dzpuis and inc_dzpuis have now an integer argument to
 * specify by which amount dzpuis is to be increased.
 * Accordingly methods dec2_dzpuis and inc2_dzpuis have been suppressed
 *
 * Revision 1.23  2003/10/27 10:49:48  e_gourgoulhon
 * Slightly modified operator<<.
 *
 * Revision 1.22  2003/10/19 19:57:00  e_gourgoulhon
 * -- Added new method spectral_display
 * -- slightly rearranged the operator<<
 *
 * Revision 1.21  2003/10/16 15:27:46  e_gourgoulhon
 * Name of method annule(int ) changed to annule_domain(int ).
 *
 * Revision 1.20  2003/10/16 14:21:36  j_novak
 * The calculation of the divergence of a Tensor is now possible.
 *
 * Revision 1.19  2003/10/13 13:52:40  j_novak
 * Better managment of derived quantities.
 *
 * Revision 1.18  2003/10/11 16:47:10  e_gourgoulhon
 * Suppressed the call to Ibtl::set_etat_qcq() after the construction
 * of the Itbl's, thanks to the new property of the Itbl class.
 *
 * Revision 1.17  2003/10/11 14:48:40  e_gourgoulhon
 * Line 344: suppressed assert(resu == -1)
 *           and added a break command to quit the loop.
 *
 * Revision 1.16  2003/10/08 14:24:09  j_novak
 * replaced mult_r_zec with mult_r_ced
 *
 * Revision 1.15  2003/10/07 15:25:38  e_gourgoulhon
 * Added call to del_derive_met in del_deriv().
 *
 * Revision 1.14  2003/10/07 09:10:00  j_novak
 * Use of ::contract instead of up()
 *
 * Revision 1.13  2003/10/06 20:51:43  e_gourgoulhon
 * In methods set: changed name "indices" to "idx" to avoid shadowing
 *  of class member.
 *
 * Revision 1.12  2003/10/06 16:17:31  j_novak
 * Calculation of contravariant derivative and Ricci scalar.
 *
 * Revision 1.11  2003/10/06 13:58:48  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.10  2003/10/05 21:11:22  e_gourgoulhon
 * - Added method std_spectral_base().
 * - Removed method change_triad() from this file.
 *
 * Revision 1.9  2003/10/03 15:09:38  j_novak
 * Improved display
 *
 * Revision 1.8  2003/10/03 11:21:48  j_novak
 * More methods for the class Metric
 *
 * Revision 1.7  2003/10/01 11:56:31  e_gourgoulhon
 * Corrected error: '=' replaced by '==' in two assert tests.
 *
 * Revision 1.6  2003/09/30 08:38:23  j_novak
 * added a header typeinfo
 *
 * Revision 1.5  2003/09/29 12:52:57  j_novak
 * Methods for changing the triad are implemented.
 *
 * Revision 1.4  2003/09/26 14:33:53  j_novak
 * Arithmetic functions for the class Tensor
 *
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
#include "metric.h"
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

    set_der_0x0() ;
	
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

    set_der_0x0() ;

}




// Standard constructor when all the indices are of the same type
// --------------------------------------------------------------
Tensor::Tensor(const Map& map, int val, int tipe, const Base_vect& triad_i) 
		: mp(&map), valence(val), triad(&triad_i), type_indice(val), 
                  n_comp(int(pow(3., val))){
    
    // Des verifs :
    assert (valence >= 0) ;
    assert ((tipe == COV) || (tipe == CON)) ;

    type_indice = tipe ;
    
    cmp = new Scalar*[n_comp] ;

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(map) ;

    set_der_0x0() ;
}	
	
// Copy constructor
// ----------------
Tensor::Tensor (const Tensor& source) : 
    mp(source.mp), valence(source.valence), triad(source.triad), 
    type_indice(source.type_indice) {
  
    n_comp = int(pow(3., valence)) ;
        
    cmp = new Scalar*[n_comp] ;
    for (int i=0 ; i<n_comp ; i++) {

        // the following writing takes into account the case where
        // source belongs to a derived class of Tensor with a different
        // storage of components :

	int place_source = source.position(indices(i)) ;
	cmp[i] = new Scalar(*source.cmp[place_source]) ;
    }

    set_der_0x0() ;

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

    set_der_0x0() ;
}


//  Constructor for a scalar field: to be used by the derived
//  class {\tt Scalar}
//-----------------------------------------------------------
Tensor::Tensor(const Map& map) : mp(&map), valence(0), triad(0x0),
		type_indice(0), n_comp(1) {
		
  cmp = new Scalar*[n_comp] ; 
  cmp[0] = 0x0 ; 
  
  set_der_0x0() ;
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

    set_der_0x0() ;
  
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

    type_indice = tipe ;
    
    cmp = new Scalar*[n_comp] ;

    for (int i=0 ; i<n_comp ; i++)
      cmp[i] = new Scalar(map) ;

    set_der_0x0() ;

}

			//--------------//
			//  Destructor  //
			//--------------//


Tensor::~Tensor () {
    
    Tensor::del_deriv() ;

    for (int i=0 ; i<n_comp ; i++) {
      if (cmp[i] != 0x0)
	delete cmp[i] ;
    }
    delete [] cmp ;
}



void Tensor::del_deriv() const {

  for (int i=0; i<N_MET_MAX; i++) 
    del_derive_met(i) ;

  set_der_0x0() ;

}

void Tensor::set_der_0x0() const {

  for (int i=0; i<N_MET_MAX; i++) 
    set_der_met_0x0(i) ;

}

void Tensor::del_derive_met(int j) const {

  assert( (j>=0) && (j<N_MET_MAX) ) ;

  if (met_depend[j] != 0x0) {
    for (int i=0 ; i<N_TENSOR_DEPEND ; i++)
      if (met_depend[j]->tensor_depend[i] == this)
		met_depend[j]->tensor_depend[i] = 0x0 ;
    if (p_derive_cov[j] != 0x0)
      delete p_derive_cov[j] ;
    if (p_derive_con[j] != 0x0)
      delete p_derive_con[j] ;
    if (p_divergence[j] != 0x0)
      delete p_divergence[j] ;

    set_der_met_0x0(j) ;
  }
}

void Tensor::set_der_met_0x0(int i) const {

  assert( (i>=0) && (i<N_MET_MAX) ) ;
  met_depend[i] = 0x0 ;
  p_derive_cov[i] = 0x0 ;
  p_derive_con[i] = 0x0 ;
  p_divergence[i] = 0x0 ;

}

int Tensor::get_place_met(const Metric& metre) const {
  int resu = -1 ;
  for (int i=0; i<N_MET_MAX; i++) 
    if (met_depend[i] == &metre) {
      	resu = i ;
		break ; 
    }
  return resu ;
}

void Tensor::set_dependance (const Metric& met) const {
    
  int nmet = 0 ;
  bool deja = false ;
  for (int i=0; i<N_MET_MAX; i++) {
    if (met_depend[i] == &met) deja = true ;
    if ((!deja) && (met_depend[i] != 0x0)) nmet++ ;
  }
  if (nmet == N_MET_MAX) {
    cout << "Too many metrics in Tensor::set_dependances" << endl ;
    abort() ;
  }
  if (!deja) { 
    int conte = 0 ;
    while ((conte < N_TENSOR_DEPEND) && (met.tensor_depend[conte] != 0x0))
      conte ++ ;
    
    if (conte == N_TENSOR_DEPEND) {
      cout << "Too many dependancies in Tensor::set_dependances " << endl ;
      abort() ;
    }
    else {
      met.tensor_depend[conte] = this ;
      met_depend[nmet] = &met ;
    }
  }
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
    	    
    for (int i=valence-1 ; i>=0 ; i--) {
		res.set(i) = div(place, 3).rem ;
		place = int((place-res(i))/3) ;
		res.set(i)++ ; 
	}
    return res ;
}

void Tensor::operator=(const Tensor& t) {
    
    assert (valence == t.valence) ;

    triad = t.triad ; 

    for (int i=0 ; i<valence ; i++)
      assert(t.type_indice(i) == type_indice(i)) ;
	
    for (int i=0 ; i<n_comp ; i++) {
      int place_t = t.position(indices(i)) ;
      *cmp[i] = *t.cmp[place_t] ;
    }

    del_deriv() ;

}

void Tensor::operator+=(const Tensor& t) {
    
    assert (valence == t.valence) ;
    assert (triad == t.triad) ; 
    for (int i=0 ; i<valence ; i++)
      assert(t.type_indice(i) == type_indice(i)) ;
	
    for (int i=0 ; i<n_comp ; i++) {
      int place_t = t.position(indices(i)) ;
      *cmp[i] += *t.cmp[place_t] ;
    }

    del_deriv() ;

}

void Tensor::operator-=(const Tensor& t) {
    
    assert (valence == t.valence) ;
    assert (triad == t.triad) ; 
    for (int i=0 ; i<valence ; i++)
      assert(t.type_indice(i) == type_indice(i)) ;
	
    for (int i=0 ; i<n_comp ; i++) {
      int place_t = t.position(indices(i)) ;
      *cmp[i] -= *t.cmp[place_t] ;
    }

    del_deriv() ;

}



// Affectation d'un tenseur d'ordre 2 :
Scalar& Tensor::set(int ind1, int ind2) {
    
    assert (valence == 2) ;
    
    Itbl ind (valence) ;
    ind.set(0) = ind1 ;
    ind.set(1) = ind2 ;
    
    int place = position(ind) ;
    
    del_deriv() ;
    return *cmp[place] ;
}

// Affectation d'un tenseur d'ordre 3 :
Scalar& Tensor::set(int ind1, int ind2, int ind3) {
    
    assert (valence == 3) ;
    
    Itbl idx(valence) ;
    idx.set(0) = ind1 ;
    idx.set(1) = ind2 ;
    idx.set(2) = ind3 ;
    int place = position(idx) ;
    del_deriv() ;
 
    return *cmp[place] ;
}


// Affectation d'un tenseur d'ordre 4 :
Scalar& Tensor::set(int ind1, int ind2, int ind3, int ind4) {
    
    assert (valence == 4) ;
    
    Itbl idx(valence) ;
    idx.set(0) = ind1 ;
    idx.set(1) = ind2 ;
    idx.set(2) = ind3 ;
    idx.set(3) = ind4 ;
    int place = position(idx) ;
    del_deriv() ;
 
    return *cmp[place] ;
}


// Affectation cas general
Scalar& Tensor::set(const Itbl& idx) {
    
    assert (idx.get_ndim() == 1) ;
    assert (idx.get_dim(0) == valence) ;
    	
    int place = position(idx) ;
    
    del_deriv() ;
    return *cmp[place] ;
}

// Annulation dans des domaines
void Tensor::annule_domain(int l) {
    
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
    idx.set(0) = indice1 ;
    idx.set(1) = indice2 ;
    return *cmp[position(idx)] ;

}

const Scalar& Tensor::operator()(int indice1, int indice2, int indice3) const {
    
    assert(valence == 3) ;
    
    Itbl idx(3) ;		
    idx.set(0) = indice1 ;
    idx.set(1) = indice2 ;
    idx.set(2) = indice3 ;
    return *cmp[position(idx)] ;
}


const Scalar& Tensor::operator()(int indice1, int indice2, int indice3,
                                 int indice4) const {
    
    assert(valence == 4) ;
    
    Itbl idx(4) ;		
    idx.set(0) = indice1 ;
    idx.set(1) = indice2 ;
    idx.set(2) = indice3 ;
    idx.set(3) = indice4 ;
    return *cmp[position(idx)] ;
}



const Scalar& Tensor::operator()(const Itbl& ind) const {
    
    assert (ind.get_ndim() == 1) ;
    assert (ind.get_dim(0) == valence) ;
    return *cmp[position(ind)] ;
    
}


// Gestion de la CED :
void Tensor::dec_dzpuis(int decrem) {
    
  for (int i=0 ; i<n_comp ; i++)
    cmp[i]->dec_dzpuis(decrem) ;
}

void Tensor::inc_dzpuis(int inc) {

    for (int i=0 ; i<n_comp ; i++)
      cmp[i]->inc_dzpuis(inc) ;
}


// Le cout :
ostream& operator<<(ostream& flux, const Tensor &source ) {

  	flux << '\n' ;
  	flux << "Lorene class : " << typeid(source).name() 
        << "           Valence : " << source.valence << '\n' ;

    if (source.get_triad() != 0x0) {
		flux << "Vectorial basis (triad) on which the components are defined :" 
	     << '\n' ; 
		flux << *(source.get_triad()) << '\n' ;
    }
    
    if (source.valence != 0) {
		flux <<    "Type of the indices : " ;
		for (int i=0 ; i<source.valence ; i++) {
			flux << "index " << i << " : " ;
			if (source.type_indice(i) == CON)
	    		flux << " contravariant." << '\n' ;
			else
	    		flux << " covariant." << '\n' ;
			if ( i < source.valence-1 ) flux << "                      " ;
		}
		flux << '\n' ; 
	}
    
    for (int i=0 ; i<source.n_comp ; i++) {

		if (source.valence == 0) {
			flux << 
			"===================== Scalar field ========================= \n" ;
		}
		else { 
      		flux << "================ Component " ;
			Itbl num_indices = source.indices(i) ;
			for (int j=0 ; j<source.valence ; j++) {
	  			flux << " " << num_indices(j) ;
      		}
			flux << " ================ \n" ; 
		}
		flux << '\n' ; 
      
      	flux << *source.cmp[i] << '\n' ;
    }
	
	return flux ;
}


void Tensor::spectral_display(double thres, int precis, ostream& ost) const {

  	ost << "Lorene class : " << typeid(*this).name() 
            << "           Valence : " << valence << '\n' ;

	for (int ic=0; ic<n_comp; ic++) {
		
		if (valence == 0) {
			ost << 
			"===================== Scalar field ========================= \n" ;
		}
		else { 
      		ost << "================ Component " ;
			Itbl num_indices = indices(ic) ;
			for (int j=0 ; j<valence ; j++) {
	  			ost << " " << num_indices(j) ;
      		}
			ost << " ================ \n" ; 
		}
		ost << '\n' ; 
		
		cmp[ic]->spectral_display(thres, precis, ost) ; 
		ost << '\n' ; 			
	}
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





// Sets the standard spectal bases of decomposition for each component

void Tensor::std_spectral_base() {

	switch (valence) {

		case 0 : {
			cmp[0]->std_spectral_base() ; 
			break ; 
		}	
		
		case 1 : {
			cout << 
			"Tensor::std_spectral_base: should not be called on a Tensor"
			<< " of valence 1 but on a Vector !" << endl ;  
			abort() ; 
			break ; 
		}
	
		case 2 : {
		
			Base_val** bases = 0x0 ; 
			if( triad->identify() == (mp->get_bvect_cart()).identify() ) {
				bases = mp->get_mg()->std_base_vect_cart() ;
			}
			else {
				assert( triad->identify() == (mp->get_bvect_spher()).identify()) ;
				bases = mp->get_mg()->std_base_vect_spher() ;
			}

	    	Itbl ind(2) ;
	    	for (int i=0 ; i<n_comp ; i++) {   
				ind = indices(i) ;
				cmp[i]->set_spectral_base( (*bases[ind(0)-1]) * 
				     (*bases[ind(1)-1]) ) ;
	    	}
	    
			for (int i=0 ; i<3 ; i++) {
				delete bases[i] ;
			}
			delete [] bases ;
			break ; 

		}
	    
	   
	    default : {

			cout << "Tensor::std_spectral_base: the case valence = " << valence
		 	<< " is not treated yet !" << endl ;
			abort() ;
			break ;
		}
	}
}


const Tensor& Tensor::derive_cov(const Metric& metre) const {
  
  set_dependance(metre) ;
  int j = get_place_met(metre) ;
  assert ((j>=0) && (j<N_MET_MAX)) ;
  if (p_derive_cov[j] == 0x0) {
    p_derive_cov[j] = metre.connect().p_derive_cov(*this) ;
  }
  return *p_derive_cov[j] ;
}

const Tensor& Tensor::derive_con(const Metric& metre) const {
  
  set_dependance(metre) ;
  int j = get_place_met(metre) ;
  assert ((j>=0) && (j<N_MET_MAX)) ;
  if (p_derive_con[j] == 0x0) {
    p_derive_con[j] = 
      new Tensor( ::contract(derive_cov(metre), valence, metre.con(), 0) ) ;
            // valence is the number of the last index of derive_cov 
            //  (the "derivation" index)
  }
  
  return *p_derive_con[j] ;

}

const Tensor& Tensor::divergence(const Metric& metre) const {
  
  set_dependance(metre) ;
  int j = get_place_met(metre) ;
  assert ((j>=0) && (j<N_MET_MAX)) ;
  if (p_divergence[j] == 0x0) {
    p_divergence[j] = metre.connect().p_divergence(*this) ;
  }
  return *p_divergence[j] ;
}















