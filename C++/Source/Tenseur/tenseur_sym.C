/*
 *  Methods of class Tenseur_sym
 *
 *   (see file tenseur.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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
 * Revision 1.10  2021/04/13 11:24:53  j_novak
 * Corrected a bug in Etoile_rot::fait_shift() which was missing the np=1 case.
 *
 * Revision 1.9  2016/12/05 16:18:17  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.8  2014/10/13 08:53:42  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.7  2014/10/06 15:13:19  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2003/03/03 19:39:58  f_limousin
 * Modification of an assert to have a check on a triad and not only on a pointer.
 *
 * Revision 1.5  2002/10/16 14:37:14  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.4  2002/09/06 14:49:25  j_novak
 * Added method lie_derive for Tenseur and Tenseur_sym.
 * Corrected various errors for derive_cov and arithmetic.
 *
 * Revision 1.3  2002/08/14 13:46:15  j_novak
 * Derived quantities of a Tenseur can now depend on several Metrique's
 *
 * Revision 1.2  2002/08/07 16:14:11  j_novak
 * class Tenseur can now also handle tensor densities, this should be transparent to older codes
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2001/10/10  13:55:23  eric
 * Modif Joachim: pow(3, *) --> pow(3., *)
 *
 * Revision 2.2  2000/02/09  19:33:38  eric
 * MODIF IMPORTANTE: la triade de decomposition est desormais passee en
 * argument des constructeurs.
 *
 * Revision 2.1  2000/01/11  11:14:41  eric
 * Gestion de la base vectorielle (triad).
 *
 * Revision 2.0  1999/12/02  17:18:37  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cassert>
#include <cmath>

// Headers Lorene
#include "tenseur.h"
#include "metrique.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
namespace Lorene {
Tenseur_sym::Tenseur_sym(const Map& map, int val, const Itbl& tipe, 
			 const Base_vect& triad_i, const Metrique* met, 
			 double weight) 
		: Tenseur(map, val, tipe, int(pow(3., val-2)) * 6, triad_i, 
			  met, weight) {

	assert (val >= 2) ;
}

// Standard constructor when all the indices are of the same type
// --------------------------------------------------------------
Tenseur_sym::Tenseur_sym(const Map& map, int val, int tipe, 
			 const Base_vect& triad_i, const Metrique* met,
			 double weight)  
		: Tenseur(map, val, tipe, int(pow(3., val-2)) * 6, triad_i,
			  met, weight) {

	assert (val >= 2) ;
}

// Copy constructor
// ----------------
Tenseur_sym::Tenseur_sym (const Tenseur_sym& source) : 
    Tenseur (*source.mp, source.valence, source.type_indice, 
	     int(pow(3., source.valence-2)*6), *(source.triad), source.metric,
	     source.poids) {
    
    assert (valence >= 2) ;   
    for (int i=0 ; i<n_comp ; i++) {
	int place_source = source.donne_place(donne_indices(i)) ;
	if (source.c[place_source] == 0x0)
	    c[i] = 0x0 ;
	else
	    c[i] = new Cmp (*source.c[place_source]) ;
    }
    etat = source.etat ;
    assert(source.met_depend != 0x0) ;
    assert(source.p_derive_cov != 0x0) ;
    assert(source.p_derive_con != 0x0) ;
    assert(source.p_carre_scal != 0x0) ;
    
    if (source.p_gradient != 0x0)
	    p_gradient = new Tenseur_sym (*source.p_gradient) ;
    
    for (int i=0; i<N_MET_MAX; i++) {
      met_depend[i] = source.met_depend[i] ;
      if (met_depend[i] != 0x0) {
	
	set_dependance (*met_depend[i]) ;
	
	if (source.p_derive_cov[i] != 0x0)
	  p_derive_cov[i] = new Tenseur (*source.p_derive_cov[i]) ;
	if (source.p_derive_con[i] != 0x0)
	  p_derive_con[i] = new Tenseur (*source.p_derive_con[i]) ;
	if (source.p_carre_scal[i] != 0x0)
	    p_carre_scal[i] = new Tenseur (*source.p_carre_scal[i]) ;
      }
    }
}   


// Constructor from a Tenseur
// --------------------------
Tenseur_sym::Tenseur_sym (const Tenseur& source) :
   Tenseur (*source.mp, source.valence, source.type_indice, 
	    int(pow(3., source.valence-2)*6), *(source.triad), source.metric,
	    source.poids) {
	
    assert (valence >= 2) ;

    for (int i=0 ; i<n_comp ; i++) {
	int place_source = source.donne_place(donne_indices(i)) ;
	if (source.c[place_source] == 0x0)
	    c[i] = 0x0 ;
	else
    	    c[i] = new Cmp (*source.c[place_source]) ;
    }
	
    etat = source.etat ;
    
    assert(source.met_depend != 0x0) ;
    assert(source.p_derive_cov != 0x0) ;
    assert(source.p_derive_con != 0x0) ;
    assert(source.p_carre_scal != 0x0) ;
    
    if (source.p_gradient != 0x0)
	    p_gradient = new Tenseur (*source.p_gradient) ;
    
}   

	
// Constructor from a file
// -----------------------
Tenseur_sym::Tenseur_sym(const Map& map, const Base_vect& triad_i, FILE* fd,
			 const Metrique* met)
			: Tenseur(map, triad_i, fd, met) {
	
	assert (valence >= 2) ;
	assert (n_comp == int(pow(3., valence-2))*6) ;
}

			//--------------//
			//  Destructor  //
			//--------------//

Tenseur_sym::~Tenseur_sym() {}




	
int Tenseur_sym::donne_place (const Itbl& idx) const {
    
    assert (idx.get_ndim() == 1) ;
    assert (idx.get_dim(0) == valence) ;
    for (int i=0 ; i<valence ; i++)
	assert ((idx(i) >= 0) && (idx(i) < 3)) ;
	
   
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
			    place_fin = 2+last ;
			    break ;
			    }
			case 2 : {
			    place_fin = 5 ;
			    break ;
			    }
			default : {
			    abort() ;
			    }
		    }
    
    int res = 0 ;
    for (int i=0 ; i<valence-2 ; i++)
	res = 3*res+idx(i) ;
    
    res = 6*res + place_fin ;
    
    return res ;
}

Itbl Tenseur_sym::donne_indices (int place) const {
    Itbl res(valence) ;
    res.set_etat_qcq() ;
    assert ((place>=0) && (place<n_comp)) ;
    
    int reste = div(place, 6).rem ;
    place = int((place-reste)/6) ;
    
    for (int i=valence-3 ; i>=0 ; i--) {
	res.set(i) = div(place, 3).rem ;
	place = int((place-res(i))/3) ;
	}
	
    if (reste<3) {
	res.set(valence-2) = 0 ;
	res.set(valence-1) = reste ;
	}
    
    if ((reste>2) && (reste<5)) {
	res.set(valence-2) = 1 ;
	res.set(valence-1) = reste - 2 ;
	}
    
    if (reste == 5) {
	res.set(valence-2) = 2 ;
	res.set(valence-1) = 2 ;
	}
 
    return res ;
}
	
void Tenseur_sym::operator= (const Tenseur& t) {
    
    assert (valence == t.get_valence()) ;
    
    triad = t.triad ; 
    poids = t.poids ;
    metric = t.metric ;
    
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
	    cout << "Unknown state in Tenseur_sym::operator= " << endl ;
	    abort() ;
	    break ;
	    }
    }
}

void Tenseur_sym::operator= (const Tenseur_sym& t) {
  operator=(t) ;
}

void Tenseur_sym::fait_gradient () const {
    
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
      	assert(*triad == mp->get_bvect_cart()) ;

	p_gradient = new Tenseur_sym(*mp, valence+1, tipe, 
				     mp->get_bvect_cart(), metric, poids) ;
    

	if (etat == ETATZERO)
	    p_gradient->set_etat_zero() ;
	else {
	    p_gradient->set_etat_qcq() ;
	    // Boucle sur les indices :
	
	    Itbl indices_source(valence) ;
	    indices_source.set_etat_qcq() ;
	
	    for (int j=0 ; j<p_gradient->n_comp ; j++) {    
		
		Itbl indices_res(p_gradient->donne_indices(j)) ;

		for (int m=0 ; m<valence ; m++)
		    indices_source.set(m) = indices_res(m+1) ;
		 
		p_gradient->set(indices_res) =  
		    (*this)(indices_source).deriv(indices_res(0)) ;
		}
	  }
    }
}


void Tenseur_sym::fait_derive_cov (const Metrique& metre, int ind) const {
    
  assert (etat != ETATNONDEF) ;
  assert (valence != 0) ;
    
  if (p_derive_cov[ind] != 0x0)
    return ;
  else {
    p_derive_cov[ind] = new Tenseur_sym (gradient()) ;
    
    if ((valence != 0) && (etat != ETATZERO)) {

      assert( *(metre.gamma().get_triad()) == *triad ) ; 

      Tenseur* auxi ;
      for (int i=0 ; i<valence ; i++) {
	
	if (type_indice(i) == COV) {
	  auxi = new Tenseur(contract(metre.gamma(), 0,(*this), i)) ;

	  Itbl indices_gamma(p_derive_cov[ind]->valence) ;
	  indices_gamma.set_etat_qcq() ;
	  //On range comme il faut :
	  for (int j=0 ; j<p_derive_cov[ind]->n_comp ; j++) {
		    
	    Itbl indices (p_derive_cov[ind]->donne_indices(j)) ;
	    indices_gamma.set(0) = indices(0) ;
	    indices_gamma.set(1) = indices(i+1) ;
	    for (int idx=2 ; idx<p_derive_cov[ind]->valence ; idx++)
	      if (idx<=i+1)
		indices_gamma.set(idx) = indices(idx-1) ;
	      else
		indices_gamma.set(idx) = indices(idx) ;
		    
	    p_derive_cov[ind]->set(indices) -= (*auxi)(indices_gamma) ;
	  }
	}   
	else {
	  auxi = new Tenseur(contract(metre.gamma(), 1, (*this), i)) ;

	  Itbl indices_gamma(p_derive_cov[ind]->valence) ;
	  indices_gamma.set_etat_qcq() ;
		
	  //On range comme il faut :
	  for (int j=0 ; j<p_derive_cov[ind]->n_comp ; j++) {
		    
	    Itbl indices (p_derive_cov[ind]->donne_indices(j)) ;
	    indices_gamma.set(0) = indices(i+1) ;
	    indices_gamma.set(1) = indices(0) ;
	    for (int idx=2 ; idx<p_derive_cov[ind]->valence ; idx++)
	      if (idx<=i+1)
		indices_gamma.set(idx) = indices(idx-1) ;
	      else
		indices_gamma.set(idx) = indices(idx) ;
	    p_derive_cov[ind]->set(indices) += (*auxi)(indices_gamma) ;
	  }
	}
	delete auxi ;
      }
    }
    if ((poids != 0.)&&(etat != ETATZERO)) 
      *p_derive_cov[ind] = *p_derive_cov[ind] - 
	poids*contract(metre.gamma(), 0, 2) * (*this) ;
    
  }
}



void Tenseur_sym::fait_derive_con (const Metrique& metre, int ind) const {
    
    if (p_derive_con[ind] != 0x0)
	return ;
    else {
	// On calcul la derivee covariante :
	if (valence != 0)
	    p_derive_con[ind] = new Tenseur_sym
		(contract(metre.con(), 1, derive_cov(metre), 0)) ;
	
    else
	p_derive_con[ind] = new Tenseur_sym
		(contract(metre.con(), 1, gradient(), 0)) ;
    }
}
}
