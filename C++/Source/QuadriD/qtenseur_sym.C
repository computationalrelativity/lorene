/*
 *  Methods of class Qtenseur_sym
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

char qtenseur_sym_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:37:13  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1  2002/09/19 09:52:43  j_novak
 * Added objects Qtenseur and Qmetrique for 4D tensor and metric handling.
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
#include "qtenseur.h"
#include "qmetrique.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
Qtenseur_sym::Qtenseur_sym(const Map& map, int val, const Itbl& tipe, 
			 const Base_vect& triad_i) 
		: Qtenseur(map, val, tipe, int(pow(4., val-2)) * 10, triad_i) {

	assert (val >= 2) ;
}

// Standard constructor when all the indices are of the same type
// --------------------------------------------------------------
Qtenseur_sym::Qtenseur_sym(const Map& map, int val, int tipe, 
			 const Base_vect& triad_i)  
		: Qtenseur(map, val, tipe, int(pow(4., val-2)) * 10, triad_i) {

	assert (val >= 2) ;
}

// Copy constructor
// ----------------
Qtenseur_sym::Qtenseur_sym (const Qtenseur_sym& source) : 
    Qtenseur (*source.mp, source.valence, source.type_indice, 
	     int(pow(4., source.valence-2)*10), *(source.triad)) {
    
    assert (valence >= 2) ;   
    for (int i=0 ; i<n_comp ; i++) {
	int place_source = source.donne_place(donne_indices(i)) ;
	if (source.c[place_source] == 0x0)
	    c[i] = 0x0 ;
	else
	    c[i] = new Cmp (*source.c[place_source]) ;
    }
    etat = source.etat ;
    if (source.p_gradient != 0x0)
	p_gradient = new Qtenseur_sym (*source.p_gradient) ;
    
    met_depend = source.met_depend ;
    if (met_depend != 0x0) {
	set_dependance(*met_depend) ;
   
	if (source.p_derive_cov != 0x0)
	    p_derive_cov = new Qtenseur_sym (*source.p_derive_cov) ;
	if (source.p_derive_con != 0x0)
	    p_derive_con = new Qtenseur_sym (*source.p_derive_con) ;
	if (source.p_carre_scal != 0x0)
	    p_carre_scal = new Qtenseur (*source.p_carre_scal) ;
    }
}

// Constructor from a Qtenseur
// --------------------------
Qtenseur_sym::Qtenseur_sym (const Qtenseur& source) :
   Qtenseur (*source.mp, source.valence, source.type_indice, 
	    int(pow(4., source.valence-2)*10), *(source.triad)) {
	
    assert (valence >= 2) ;

    for (int i=0 ; i<n_comp ; i++) {
	int place_source = source.donne_place(donne_indices(i)) ;
	if (source.c[place_source] == 0x0)
	    c[i] = 0x0 ;
	else
    	    c[i] = new Cmp (*source.c[place_source]) ;
    }
	
    etat = source.etat ;
    
    if (source.p_gradient != 0x0)
	p_gradient = new Qtenseur (*source.p_gradient) ;
    met_depend = 0x0 ;
    set_der_met_0x0() ;
}
	
// Constructor from a file
// -----------------------
Qtenseur_sym::Qtenseur_sym(const Map& map, const Base_vect& triad_i, FILE* fd)
			: Qtenseur(map, triad_i, fd) {
	
	assert (valence >= 2) ;
	assert (n_comp == int(pow(4., valence-2))*10) ;
}

			//--------------//
			//  Destructor  //
			//--------------//

Qtenseur_sym::~Qtenseur_sym() {}




	
int Qtenseur_sym::donne_place (const Itbl& idx) const {
    
    assert (idx.get_ndim() == 1) ;
    assert (idx.get_dim(0) == valence) ;
    for (int i=0 ; i<valence ; i++)
	assert ((idx(i) >= 0) && (idx(i) < 4)) ;
	
   
     // Gestion des deux derniers indices :
    int last = idx(valence-1) ;
    int lastm1 = idx(valence-2) ;
    if (last < lastm1) {
	int auxi = last ;
	last = lastm1 ;
	lastm1 = auxi ;
    }
    
    int place_fin ;
    switch (lastm1) {
			case 0 : {
			    place_fin = last ;
			    break ;
			    }
			case 1 : {
			    place_fin = 3+last ;
			    break ;
			    }
			case 2 : {
			    place_fin = 5+last ;
			    break ;
			    }
                       case 3 : {
			    place_fin = 9 ;
			    break ;
		            }
			default : {
			    abort() ;
			    }
		    }
    
    int res = 0 ;
    for (int i=0 ; i<valence-2 ; i++)
	res = 4*res+idx(i) ;
    
    res = 10*res + place_fin ;
    
    return res ;
}

Itbl Qtenseur_sym::donne_indices (int place) const {
    Itbl res(valence) ;
    res.set_etat_qcq() ;
    assert ((place>=0) && (place<n_comp)) ;
    
    int reste = div(place, 10).rem ;
    place = int((place-reste)/10) ;
    
    for (int i=valence-3 ; i>=0 ; i--) {
	res.set(i) = div(place, 4).rem ;
	place = int((place-res(i))/4) ;
	}
	
    if (reste<4) {
	res.set(valence-2) = 0 ;
	res.set(valence-1) = reste ;
	}
    
    if ((reste>3) && (reste<7)) {
	res.set(valence-2) = 1 ;
	res.set(valence-1) = reste - 3 ;
	}
    
    if ((reste>6) && (reste<9)) {
	res.set(valence-2) = 2 ;
	res.set(valence-1) = reste - 5 ;
	}
    
    if (reste == 9) {
	res.set(valence-2) = 3 ;
	res.set(valence-1) = 3 ;
	}
 
    return res ;
}
	
void Qtenseur_sym::operator= (const Qtenseur& t) {
    
    assert (valence == t.get_valence()) ;
    
    triad = t.triad ; 
    
    for (int i=0 ; i<valence ; i++)
	assert (type_indice(i) == t.type_indice(i)) ;
    
    switch (t.get_etat()) {
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
	    cout << "Unknown state in Qtenseur_sym::operator= " << endl ;
	    abort() ;
	    break ;
	    }
    }
}

void Qtenseur_sym::fait_gradient (const Param& para) const {
    
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
    assert (*triad == mp->get_bvect_cart()) ;

    p_gradient = new Qtenseur_sym(*mp, valence+1, tipe, mp->get_bvect_cart()) ;
    p_gradient->set_etat_qcq() ;
    assert (para.get_n_double() == 1) ;
    double dt = para.get_double() ;
    Qtenseur_sym dtemp(*mp, valence, type_indice, *triad) ;
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
      cout << "Problem in Qtenseur_sym::gradient!" << endl ;
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

void Qtenseur_sym::fait_derive_cov (const Param& para, const Qmetrique& metre) const {
    
    assert (etat != ETATNONDEF) ;
    assert (valence != 0) ;
    
    if (p_derive_cov != 0x0)
	return ;
    else {
	p_derive_cov = new Qtenseur_sym (gradient(para)) ;
	assert(para.get_n_double() == 1) ;
	double dt = para.get_double() ;
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



void Qtenseur_sym::fait_derive_con (const Param& para, const Qmetrique& metre) const {
    
    if (p_derive_con != 0x0)
	return ;
    else {
	// On calcul la derivee covariante :
	if (valence != 0)
	    p_derive_con = new Qtenseur_sym
		(contract(metre.con(), 1, derive_cov(para, metre), 0)) ;
	
    else
	p_derive_con = new Qtenseur_sym
		(contract(metre.con(), 1, gradient(para), 0)) ;
    }
}

Qtenseur_sym operator*(const Qtenseur& t1, const Qtenseur_sym& t2) {
   
    assert ((t1.get_etat() != ETATNONDEF) && (t2.etat != ETATNONDEF)) ;
    assert (t1.get_mp() == t2.mp) ;
    
    int val_res = t1.get_valence() + t2.valence ;
    
   
    Itbl tipe (val_res) ;
    tipe.set_etat_qcq() ;
    for (int i=0 ; i<t1.get_valence() ; i++)
	tipe.set(i) = t1.get_type_indice(i) ;
    for (int i=0 ; i<t2.valence ; i++)
	tipe.set(i+t1.get_valence()) = t2.type_indice(i) ;
	

    if ( t1.get_valence() != 0 ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }

    Qtenseur_sym res(*t1.get_mp(), val_res, tipe, *(t2.get_triad()) ) ;
	

    if ((t1.get_etat() == ETATZERO) || (t2.etat == ETATZERO))
	res.set_etat_zero() ;
    else {
	res.set_etat_qcq() ;
	Itbl jeux_indice_t1 (t1.get_valence()) ;
	jeux_indice_t1.set_etat_qcq() ;
	Itbl jeux_indice_t2 (t2.valence) ;
	jeux_indice_t2.set_etat_qcq() ;
	    
	for (int i=0 ; i<res.n_comp ; i++) {
	    Itbl jeux_indice_res(res.donne_indices(i)) ;
	    for (int j=0 ; j<t1.get_valence() ; j++)
		jeux_indice_t1.set(j) = jeux_indice_res(j) ;
	    for (int j=0 ; j<t2.valence ; j++)
		jeux_indice_t2.set(j) = jeux_indice_res(j+t1.get_valence()) ;
		
	    res.set(jeux_indice_res) = t1(jeux_indice_t1)*t2(jeux_indice_t2) ;
	}
    }
    return res ;
}

