/*
 *  Methods of class Tenseur
 *
 *   (see file tenseur.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char tenseur_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.10  2002/08/30 13:21:38  j_novak
 * Corrected error in constructor
 *
 * Revision 1.9  2002/08/14 13:46:15  j_novak
 * Derived quantities of a Tenseur can now depend on several Metrique's
 *
 * Revision 1.8  2002/08/13 08:02:45  j_novak
 * Handling of spherical vector/tensor components added in the classes
 * Mg3d and Tenseur. Minor corrections for the class Metconf.
 *
 * Revision 1.7  2002/08/08 15:10:45  j_novak
 * The flag "plat" has been added to the class Metrique to show flat metrics.
 *
 * Revision 1.6  2002/08/07 16:14:11  j_novak
 * class Tenseur can now also handle tensor densities, this should be transparent to older codes
 *
 * Revision 1.5  2002/08/02 15:07:41  j_novak
 * Member function determinant has been added to the class Metrique.
 * A better handling of spectral bases is now implemented for the class Tenseur.
 *
 * Revision 1.4  2002/05/07 07:36:03  e_gourgoulhon
 * Compatibilty with xlC compiler on IBM SP2:
 *    suppressed the parentheses around argument of instruction new:
 * 	e.g.   t = new (Tbl *[nzone])  -->   t = new Tbl*[nzone]
 *
 * Revision 1.3  2002/05/02 15:16:22  j_novak
 * Added functions for more general bi-fluid EOS
 *
 * Revision 1.2  2001/12/04 21:27:54  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.21  2001/10/10  13:54:40  eric
 * Modif Joachim: pow(3, *) --> pow(3., *)
 *
 * Revision 2.20  2000/12/20  09:50:08  eric
 * Correction erreur dans operator<< : la sortie doit etre flux et non cout !
 *
 * Revision 2.19  2000/10/12  13:11:23  eric
 * Methode set_std_base(): traitement du cas etat = ETATZERO (return).
 *
 * Revision 2.18  2000/09/13  12:11:40  eric
 * Ajout de la fonction allocate_all().
 *
 * Revision 2.17  2000/05/22  14:40:09  phil
 * ajout de inc_dzpuis et dec_dzpuis
 *
 * Revision 2.16  2000/03/22  09:18:57  eric
 * Traitement du cas ETATZERO dans dec2_dzpuis, inc2_dzpuis et mult_r_zec.
 *
 * Revision 2.15  2000/02/12  11:35:58  eric
 * Modif de la fonction set_std_base : appel de Valeur::set_base plutot
 * que l'affectation directe du membre Valeur::base.
 *
 * Revision 2.14  2000/02/10  18:30:47  eric
 * La fonction set_triad ne fait plus que l'affectation du membre triad.
 *
 * Revision 2.13  2000/02/10  16:11:07  eric
 * Ajout de la fonction change_triad.
 *
 * Revision 2.12  2000/02/09  19:32:39  eric
 * MODIF IMPORTANTE: la triade de decomposition est desormais passee en
 * argument des constructeurs.
 *
 * Revision 2.11  2000/01/24  13:02:39  eric
 * Traitement du cas triad = 0x0 dans la sauvegarde/lecture fichier
 * Constructeur par lecture de fichier: met_depend est desormais initialise
 *  a 0x0.
 *
 * Revision 2.10  2000/01/20  16:02:57  eric
 * Ajout des operator=(double ) et operator=(int ).
 *
 * Revision 2.9  2000/01/20  15:34:39  phil
 * changement traid dans fait_gradient()
 *
 * Revision 2.8  2000/01/14  14:07:26  eric
 * Ajout de la fonction annule.
 *
 * Revision 2.7  2000/01/13  14:10:53  eric
 * Ajout du constructeur par copie d'un Cmp (pour un scalaire)
 * ainsi que l'affectation a un Cmp.
 *
 * Revision 2.6  2000/01/13  13:46:38  eric
 * Ajout du membre p_gradient_spher et des fonctions fait_gradient_spher(),
 *  gradient_spher() pour le calcul du gradient d'un scalaire en
 *  coordonnees spheriques sur la triade spherique associee.
 *
 * Revision 2.5  2000/01/12  13:19:04  eric
 * Les operator::(...) renvoient desormais une reference const sur le c[...]
 * correspondant et non plus un Cmp copie de c[...].
 * (ceci grace a la nouvelle fonction Map::cmp_zero()).
 *
 * Revision 2.4  2000/01/11  11:13:59  eric
 * Changement de nom pour la base vectorielle : base --> triad
 *
 * Revision 2.3  2000/01/10  17:23:07  eric
 * Modif affichage.
 * Methode fait_derive_cov : ajout de
 *   assert( metre.gamma().get_base() == base )
 * Methode set_std_base : ajout de
 *    assert( *base == mp->get_bvect_cart() )
 *
 * Revision 2.2  2000/01/10  15:15:26  eric
 * Ajout du membre base (base vectorielle sur laquelle sont definies
 *   les composantes).
 *
 * Revision 2.1  1999/12/09  12:39:23  phil
 * changement prototypage des derivees
 *
 * Revision 2.0  1999/12/02  17:18:31  phil
 * *** empty log message ***
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
#include "tenseur.h"
#include "metrique.h"
#include "utilitaires.h"

			//--------------//
			// Constructors //
			//--------------//
// Consistency check for tensor densities
//---------------------------------------
bool Tenseur::verif() const {
  return ( (poids == 0.) || (metric != 0x0) ) ;
}

void Tenseur::new_der_met() {
    met_depend = new const Metrique*[N_MET_MAX] ;
    p_derive_cov = new Tenseur*[N_MET_MAX] ;
    p_derive_con = new Tenseur*[N_MET_MAX] ;
    p_carre_scal = new Tenseur*[N_MET_MAX] ;
    for (int i=0; i<N_MET_MAX; i++) {
      met_depend[i] = 0x0 ;
    }
    set_der_0x0() ;
}

// Constructor for a scalar field
// ------------------------------
Tenseur::Tenseur (const Map& map, const Metrique* met, double weight) : 
		mp(&map), valence(0), triad(0x0),
		type_indice(0), n_comp(1), etat(ETATNONDEF), poids(weight),
		metric(met) {
    
    assert(verif()) ;
    c = new Cmp*[n_comp] ;
    c[0] = 0x0 ;
    new_der_met() ;
}



// Constructor for a scalar field and from a {\tt Cmp} 
// ---------------------------------------------------
Tenseur::Tenseur (const Cmp& ci, const Metrique* met, double weight) : 
		mp(ci.get_mp()), valence(0), triad(0x0),
		type_indice(0), n_comp(1), etat(ci.get_etat()), poids(weight),
		metric(met){
    
    assert(ci.get_etat() != ETATNONDEF) ; 
    assert(verif()) ;
    
    c = new Cmp*[n_comp] ;

    if ( ci.get_etat() != ETATZERO ) {
      assert( ci.get_etat() == ETATQCQ ) ; 
      c[0] = new Cmp(ci) ;
    }
    else {
      c[0] = 0x0 ;
    }
    new_der_met() ;
}

// Standard constructor 
// --------------------
Tenseur::Tenseur(const Map& map, int val, const Itbl& tipe, 
		 const Base_vect& triad_i, const Metrique* met, double weight) 
		: mp(&map), valence(val), triad(&triad_i), type_indice(tipe), 
		   n_comp(int(pow(3., val))), etat(ETATNONDEF), poids(weight),
		metric(met){
		
    // Des verifs :
    assert (valence >= 0) ;
    assert (tipe.get_ndim() == 1) ;
    assert (valence == tipe.get_dim(0)) ;
    for (int i=0 ; i<valence ; i++)
	assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
    assert(verif()) ;
    
    c = new Cmp*[n_comp] ;
    for (int i=0 ; i<n_comp ; i++)
	c[i] = 0x0 ;
    new_der_met() ;
}

// Standard constructor with the triad passed as a pointer
// -------------------------------------------------------
Tenseur::Tenseur(const Map& map, int val, const Itbl& tipe, 
		 const Base_vect* triad_i, const Metrique* met, double weight) 
		: mp(&map), valence(val), triad(triad_i), type_indice(tipe), 
		   n_comp(int(pow(3., val))), etat(ETATNONDEF), poids(weight),
		metric(met){
		
    // Des verifs :
    assert (valence >= 0) ;
    assert (tipe.get_ndim() == 1) ;
    assert (valence == tipe.get_dim(0)) ;
    for (int i=0 ; i<valence ; i++)
	assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
    assert(verif()) ;
    
    if (valence == 0) {	    // particular case of a scalar 
	triad = 0x0 ; 
    }   
    
    c = new Cmp*[n_comp] ;
    for (int i=0 ; i<n_comp ; i++)
	c[i] = 0x0 ;
    new_der_met() ;
}




// Standard constructor when all the indices are of the same type
// --------------------------------------------------------------
Tenseur::Tenseur(const Map& map, int val, int tipe, const Base_vect& triad_i,
		 const Metrique* met, double weight) 
		: mp(&map), valence(val), triad(&triad_i), type_indice(val), 
                  n_comp(int(pow(3., val))), etat (ETATNONDEF), poids(weight), 
		  metric(met){
    
    // Des verifs :
    assert (valence >= 0) ;
    assert ((tipe == COV) || (tipe == CON)) ;
    assert(verif()) ;
    type_indice.set_etat_qcq() ;
    type_indice = tipe ;
    
    c = new Cmp*[n_comp] ;
    for (int i=0 ; i<n_comp ; i++)
	c[i] = 0x0 ;
    new_der_met() ;
}	
	
// Copy constructor
// ----------------
Tenseur::Tenseur (const Tenseur& source) : 
    mp(source.mp), valence(source.valence), triad(source.triad), 
    type_indice(source.type_indice), etat (source.etat), poids(source.poids),
    metric(source.metric) {
  
    assert(verif()) ;
    
    n_comp = int(pow(3., valence)) ;
        
    c = new Cmp*[n_comp] ;
    for (int i=0 ; i<n_comp ; i++) {
	int place_source = source.donne_place(donne_indices(i)) ;
	if (source.c[place_source] == 0x0)
	    c[i] = 0x0 ;
	else
	    c[i] = new Cmp(*source.c[place_source]) ;
    }

    assert(source.met_depend != 0x0) ;
    assert(source.p_derive_cov != 0x0) ;
    assert(source.p_derive_con != 0x0) ;
    assert(source.p_carre_scal != 0x0) ;
    new_der_met() ;
    
    if (source.p_gradient != 0x0)
	    p_gradient = new Tenseur (*source.p_gradient) ;
    
    if (source.p_gradient_spher != 0x0)
	    p_gradient_spher = new Tenseur (*source.p_gradient_spher) ;

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

// Constructor from a symmetric tensor
// -----------------------------------
Tenseur::Tenseur (const Tenseur_sym& source) :
	mp(source.mp), valence(source.valence), triad(source.triad), 
	type_indice(source.type_indice), etat(source.etat), 
	poids(source.poids), metric(source.metric) {
    
    assert(verif()) ;
    n_comp = int(pow(3., valence)) ;
        
    c = new Cmp*[n_comp] ;
    for (int i=0 ; i<n_comp ; i++) {
	int place_source = source.donne_place(donne_indices(i)) ;
	if (source.c[place_source] == 0x0)
	    c[i] = 0x0 ;
	else
	    c[i] = new Cmp(*source.c[place_source]) ;
    }

    assert(source.met_depend != 0x0) ;
    assert(source.p_derive_cov != 0x0) ;
    assert(source.p_derive_con != 0x0) ;
    assert(source.p_carre_scal != 0x0) ;
    new_der_met() ;
    
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

// Constructor from a file
// -----------------------
Tenseur::Tenseur(const Map& mapping, const Base_vect& triad_i, FILE* fd, 
		 const Metrique* met)
		 : mp(&mapping), triad(&triad_i), type_indice(fd), 
                   metric(met) {
   
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

    new_der_met() ;
    
    if (met == 0x0) 
      poids = 0. ;
    else
      fread_be(&poids, sizeof(double), 1, fd) ;
}


// Constructor from a file for a scalar field
// ------------------------------------------
Tenseur::Tenseur (const Map& mapping, FILE* fd, const Metrique* met) 
		 : mp(&mapping), type_indice(fd), metric(met){
   
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
    
    new_der_met() ;

    if (met == 0x0) 
      poids = 0. ;
    else
      fread_be(&poids, sizeof(double), 1, fd) ;
}




// Constructor used by the derived classes
// ---------------------------------------
Tenseur::Tenseur (const Map& map, int val, const Itbl& tipe, int compo, 
		const Base_vect& triad_i, const Metrique* met, double weight) :
     mp(&map), valence(val), triad(&triad_i), type_indice(tipe), n_comp(compo), 
	    etat (ETATNONDEF), poids(weight), metric(met) {
     
    // Des verifs :
    assert (valence >= 0) ;
    assert (tipe.get_ndim() == 1) ;
    assert (n_comp > 0) ;   
    assert (valence == tipe.get_dim(0)) ;
    for (int i=0 ; i<valence ; i++)
	assert ((tipe(i) == COV) || (tipe(i) == CON)) ;
    assert(verif()) ;
    
    c = new Cmp*[n_comp] ;
    for (int i=0 ; i<n_comp ; i++)
	c[i] = 0x0 ;
    
    new_der_met() ;
}

// Constructor used by the derived classes when all the indices are of 
// the same type.
// -------------------------------------------------------------------
Tenseur::Tenseur (const Map& map, int val, int tipe, int compo, 
		const Base_vect& triad_i, const Metrique* met, double weight) :
     mp(&map), valence(val), triad(&triad_i), type_indice(val), n_comp(compo), 
	    etat (ETATNONDEF), poids(weight), metric(met) {
    // Des verifs :
    assert (valence >= 0) ;
    assert (n_comp >= 0) ;
    assert ((tipe == COV) || (tipe == CON)) ;
    assert(verif()) ;
    type_indice.set_etat_qcq() ;
    type_indice = tipe ;
    
    c = new Cmp*[n_comp] ;
    for (int i=0 ; i<n_comp ; i++)
	c[i] = 0x0 ;

    new_der_met() ;
}

			//--------------//
			//  Destructor  //
			//--------------//


Tenseur::~Tenseur () {
    
    del_t() ;
    delete [] met_depend ;
    delete [] p_derive_cov ;
    delete [] p_derive_con ;
    delete [] p_carre_scal ;
    delete [] c ;
}



void Tenseur::del_t() {
    del_derive() ;
    for (int i=0 ; i<n_comp ; i++)
	if (c[i] != 0x0) {
	    delete c[i] ;
	    c[i] = 0x0 ;
	    }
}

void Tenseur::del_derive_met(int j) const {

  assert( (j>=0) && (j<N_MET_MAX) ) ;
  // On gere la metrique ...
  if (met_depend[j] != 0x0) {
    for (int i=0 ; i<N_DEPEND ; i++)
      if (met_depend[j]->dependances[i] == this)
	met_depend[j]->dependances[i] = 0x0 ;
    if (p_derive_cov[j] != 0x0)
      delete p_derive_cov[j] ;
    if (p_derive_con[j] != 0x0)
      delete p_derive_con[j] ;
    if (p_carre_scal[j] != 0x0)
      delete p_carre_scal[j] ;
    set_der_met_0x0(j) ;
  }
}


void Tenseur::del_derive () const {
  for (int i=0; i<N_MET_MAX; i++) 
    del_derive_met(i) ;
  if (p_gradient != 0x0)
    delete p_gradient ;
  if (p_gradient_spher != 0x0)
    delete p_gradient_spher ;
  set_der_0x0() ;
}

void Tenseur::set_der_met_0x0(int i) const {
  met_depend[i] = 0x0 ;
    p_derive_cov[i] = 0x0 ;
    p_derive_con[i] = 0x0 ;
    p_carre_scal[i] = 0x0 ;
}


void Tenseur::set_der_0x0() const {
  for (int i=0; i<N_MET_MAX; i++) 
    set_der_met_0x0(i) ;
  p_gradient = 0x0 ;   
  p_gradient_spher = 0x0 ;   
}

int Tenseur::get_place_met(const Metrique& metre) const {
  int resu = -1 ;
  for (int i=0; i<N_MET_MAX; i++) 
    if (met_depend[i] == &metre) {
      assert(resu == -1) ;
      resu = i ;
    }
  return resu ;
}

void Tenseur::set_dependance (const Metrique& met) const {
    
  int nmet = 0 ;
  bool deja = false ;
  for (int i=0; i<N_MET_MAX; i++) {
    if (met_depend[i] == &met) deja = true ;
    if ((!deja) && (met_depend[i] != 0x0)) nmet++ ;
  }
  if (nmet == N_MET_MAX) {
    cout << "Too many metrics in Tenseur::set_dependances" << endl ;
    abort() ;
  }
  if (!deja) { 
    int conte = 0 ;
    while ((conte < N_DEPEND) && (met.dependances[conte] != 0x0))
      conte ++ ;
    
    if (conte == N_DEPEND) {
      cout << "Too many dependancies in Tenseur::set_dependances " << endl ;
      abort() ;
    }
    else {
      met.dependances[conte] = this ;
      met_depend[nmet] = &met ;
    }
  }
}

void Tenseur::set_etat_qcq() { 
    
    del_derive() ;
    for (int i=0 ; i<n_comp ; i++)
	if (c[i] == 0x0)
	    c[i] = new Cmp(mp) ;
    etat = ETATQCQ ;
}

void Tenseur::set_etat_zero() { 
    del_t() ;
    etat = ETATZERO ;
}

void Tenseur::set_etat_nondef() { 
    del_t() ;
    etat = ETATNONDEF ;
}

// Allocates everything
// --------------------
void Tenseur::allocate_all() {
    
	set_etat_qcq() ; 
	for (int i=0 ; i<n_comp ; i++) {
	    c[i]->allocate_all() ; 
	}
	
} 



void Tenseur::change_triad(const Base_vect& bi) {
    
    bi.change_basis(*this) ; 
    
}

void Tenseur::set_triad(const Base_vect& bi) {
    
    triad = &bi ; 
    
}

void Tenseur::set_poids(double weight) {

    poids = weight ;
}

void Tenseur::set_metric(const Metrique& met) {

    metric = &met ;
}

int Tenseur::donne_place (const Itbl& idx) const {
    
    assert (idx.get_ndim() == 1) ;
    assert (idx.get_dim(0) == valence) ;
    
    for (int i=0 ; i<valence ; i++)
	assert ((idx(i)>=0) && (idx(i)<3)) ;
    int res = 0 ;
    for (int i=0 ; i<valence ; i++)
        res = 3*res+idx(i) ;
    
    return res;
}

Itbl Tenseur::donne_indices (int place) const {
    
    assert ((place >= 0) && (place < n_comp)) ;

    Itbl res(valence) ;
    res.set_etat_qcq() ;
    	    
    for (int i=valence-1 ; i>=0 ; i--) {
	res.set(i) = div(place, 3).rem ;
	place = int((place-res(i))/3) ;
	}
    return res ;
}

void Tenseur::operator=(const Tenseur& t) {
    
    assert (valence == t.valence) ;

    triad = t.triad ; 
    poids = t.poids ;
    metric = t.metric ;

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
	    cout << "Unknown state in Tenseur::operator= " << endl ;
	    abort() ;
	    break ;
	    }
    }
}


void Tenseur::operator=(const Cmp& ci) {
    
    assert (valence == 0) ;
    poids = 0. ;
    metric = 0x0 ;

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
	    cout << "Unknown state in Tenseur::operator= " << endl ;
	    abort() ;
	    break ;
	}
    }
}

void Tenseur::operator=(double x) {

    poids = 0. ;
    metric = 0x0 ;
    if (x == double(0)) {
	set_etat_zero() ;
    }
    else {
	assert(valence == 0) ; 
	set_etat_qcq() ;
	*(c[0]) = x ; 
    }

}

void Tenseur::operator=(int x) {

    poids = 0. ;
    metric = 0x0 ;
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
Cmp& Tenseur::set () {
    
    del_derive() ;
    assert(etat == ETATQCQ) ;
    assert (valence == 0) ;  
    return *c[0] ;
}


// Affectation d'un vecteur :
Cmp& Tenseur::set (int ind) {
    
    del_derive() ;
    assert(valence == 1) ;
    assert (etat == ETATQCQ) ;
    assert ((ind >= 0) && (ind < 3)) ;
    
    return *c[ind] ;
}

// Affectation d'un tenseur d'ordre 2 :
Cmp& Tenseur::set (int ind1, int ind2) {
    
    del_derive() ;
    assert (valence == 2) ;
    assert (etat == ETATQCQ) ;
    assert ((ind1 >= 0) && (ind1 < 3)) ;
    assert ((ind2 >= 0) && (ind2 < 3)) ;
    
    Itbl ind (valence) ;
    ind.set_etat_qcq() ;
    ind.set(0) = ind1 ;
    ind.set(1) = ind2 ;
    
    int place = donne_place(ind) ;
    
    return *c[place] ;
}

// Affectation d'un tenseur d'ordre 3 :
Cmp& Tenseur::set (int ind1, int ind2, int ind3) {
    
    del_derive() ;
    assert (valence == 3) ;
    assert (etat == ETATQCQ) ;
    assert ((ind1 >= 0) && (ind1 < 3)) ;
    assert ((ind2 >= 0) && (ind2 < 3)) ;
    assert ((ind3 >= 0) && (ind3 < 3)) ;
    
    Itbl indices(valence) ;
    indices.set_etat_qcq() ;
    indices.set(0) = ind1 ;
    indices.set(1) = ind2 ;
    indices.set(2) = ind3 ;
    int place = donne_place(indices) ;
 
    return *c[place] ;
}

// Affectation cas general
Cmp& Tenseur::set(const Itbl& indices) {
    
    assert (indices.get_ndim() == 1) ;
    assert (indices.get_dim(0) == valence) ;
    
    del_derive() ;
    assert (etat == ETATQCQ) ;
    for (int i=0 ; i<valence ; i++)
	assert ((indices(i)>=0) && (indices(i)<3)) ;
	
    int place = donne_place(indices) ;
    
    return *c[place] ;
}

// Annulation dans des domaines
void Tenseur::annule(int l) {
    
    annule(l, l) ;     
}

void Tenseur::annule(int l_min, int l_max) {
    
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
	for (int j=0; j<N_MET_MAX; j++) {
	  if (p_derive_cov[j] != 0x0) p_derive_cov[j]->annule(l_min, l_max) ;
	  if (p_derive_con[j] != 0x0) p_derive_con[j]->annule(l_min, l_max) ;
	  if (p_carre_scal[j] != 0x0) p_carre_scal[j]->annule(l_min, l_max) ;
	}

    }
    
}




// Exctraction :
const Cmp& Tenseur::operator()() const {
    
    assert(valence == 0) ;
    
    if (etat == ETATQCQ) return *c[0] ;	    // pour la performance,
					    // ce cas est traite en premier,
					    // en dehors du switch
    switch (etat) {
	
	case ETATNONDEF : {
	    cout << "Undefined Tensor in Tenseur::operator() ..." << endl ;
	    abort() ;
	    return *c[0] ;  // bidon pour satisfaire le compilateur
	    }
	    
	case ETATZERO : {
	    return mp->cmp_zero() ;
	    }
	    
	    
	default : {
	    cout <<"Unknown state in Tenseur::operator()" << endl ;
	    abort() ;
	    return *c[0] ;  // bidon pour satisfaire le compilateur
	    }
	}
}


const Cmp& Tenseur::operator() (int indice) const {
    
    assert ((indice>=0) && (indice<3)) ;
    assert(valence == 1) ;
    
    if (etat == ETATQCQ) return *c[indice] ;	 // pour la performance,
					    // ce cas est traite en premier,
					    // en dehors du switch
    switch (etat) {
	
	case ETATNONDEF : {
	    cout << "Undefined Tensor in Tenseur::operator(int) ..." << endl ;
	    abort() ;
	    return *c[0] ;  // bidon pour satisfaire le compilateur
	    }
	    
	case ETATZERO : {
	    return mp->cmp_zero() ;
	    }
	    	    
	default : {
	    cout <<"Unknown state in Tenseur::operator(int)" << endl ;
	    abort() ;
	    return *c[0] ;  // bidon pour satisfaire le compilateur
	    }
	}
}

const Cmp& Tenseur::operator() (int indice1, int indice2) const {
    
    assert ((indice1>=0) && (indice1<3)) ;
    assert ((indice2>=0) && (indice2<3)) ;
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
	    cout << "Undefined Tensor in Tenseur::operator(int, int) ..." << endl ;
	    abort() ;
	    return *c[0] ;  // bidon pour satisfaire le compilateur
	    }
	    
	case ETATZERO : {
	    return mp->cmp_zero() ;
	    }
	   	    
	default : {
	    cout <<"Unknown state in Tenseur::operator(int, int)" << endl ;
	    abort() ;
	    return *c[0] ;  // bidon pour satisfaire le compilateur
	    }
	} 
}

const Cmp& Tenseur::operator() (int indice1, int indice2, int indice3) const {
    
    assert ((indice1>=0) && (indice1<3)) ;
    assert ((indice2>=0) && (indice2<3)) ;
    assert ((indice3>=0) && (indice3<3)) ;
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
	    cout << "Undefined Tensor in Tenseur::operator(int, int, int) ..." << endl ;
	    abort() ;
	    return *c[0] ;  // bidon pour satisfaire le compilateur
	    }
	    
	case ETATZERO : {
	    return mp->cmp_zero() ;
	    }
	    	    
	default : {
	    cout <<"Unknown state in Tenseur::operator(int, int, int)" << endl ;
	    abort() ;
	    return *c[0] ;  // bidon pour satisfaire le compilateur
	    }
	}
}


const Cmp& Tenseur::operator() (const Itbl& indices) const {
    
    assert (indices.get_ndim() == 1) ;
    assert (indices.get_dim(0) == valence) ;
    for (int i=0 ; i<valence ; i++)
	assert ((indices(i)>=0) && (indices(i)<3)) ;
    
    if (etat == ETATQCQ) {		    // pour la performance,
	return *c[donne_place(indices)]	;   // ce cas est traite en premier,
    }					    // en dehors du switch

    switch (etat) {
	
	case ETATNONDEF : {
	    cout << "Undefined Tensor in Tenseur::operator(const Itbl&) ..." << endl ;
	    abort() ;
	    return *c[0] ;  // bidon pour satisfaire le compilateur
	    }
	    
	case ETATZERO : {
	    return mp->cmp_zero() ;
	    }
	    
	default : {
	    cout <<"Unknown state in Tenseur::operator(const Itbl& )" << endl ;
	    abort() ;
	    return *c[0] ;  // bidon pour satisfaire le compilateur
	    }
	}
    
}

// Gestion de la ZEC :
void Tenseur::dec_dzpuis() {
    
    if (etat == ETATZERO) {
	return ; 
    }
    
    assert(etat == ETATQCQ) ;
   
    for (int i=0 ; i<n_comp ; i++)
	if (c[i] != 0x0)
	    c[i]->dec_dzpuis() ;
}

void Tenseur::inc_dzpuis() {
    
    if (etat == ETATZERO) {
	return ; 
    }

    assert(etat == ETATQCQ) ;
    
    for (int i=0 ; i<n_comp ; i++)
	if (c[i] != 0x0)
	    c[i]->inc_dzpuis() ;
}

void Tenseur::dec2_dzpuis() {
    
    if (etat == ETATZERO) {
	return ; 
    }
    
    assert(etat == ETATQCQ) ;
   
    for (int i=0 ; i<n_comp ; i++)
	if (c[i] != 0x0)
	    c[i]->dec2_dzpuis() ;
}

void Tenseur::inc2_dzpuis() {
    
    if (etat == ETATZERO) {
	return ; 
    }

    assert(etat == ETATQCQ) ;
    
    for (int i=0 ; i<n_comp ; i++)
	if (c[i] != 0x0)
	    c[i]->inc2_dzpuis() ;
}

void Tenseur::mult_r_zec() {
    
    if (etat == ETATZERO) {
	return ; 
    }

    assert(etat == ETATQCQ) ;
    
    for (int i=0 ; i<n_comp ; i++) 
    	if (c[i] != 0x0)
	    c[i]->mult_r_zec() ;
}

// Gestion des bases spectrales (valence <= 2)
void Tenseur::set_std_base() {
    
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

	  if ( triad->identify() == (mp->get_bvect_cart()).identify() ) {

	    Base_val** bases = mp->get_mg()->std_base_vect_cart() ;

	    for (int i=0 ; i<3 ; i++)
		(c[i]->va).set_base( *bases[i] ) ;
	    for (int i=0 ; i<3 ; i++)
		delete bases[i] ;
	    delete [] bases ;
	  }
	  else {
	    assert( triad->identify() == (mp->get_bvect_spher()).identify()) ;
	    Base_val** bases = mp->get_mg()->std_base_vect_spher() ;

	    for (int i=0 ; i<3 ; i++)
		(c[i]->va).set_base( *bases[i] ) ;
	    for (int i=0 ; i<3 ; i++)
		delete bases[i] ;
	    delete [] bases ;
	  }
	  break ;
	    
	}
	    
	case 2 : {

	  if( triad->identify() == (mp->get_bvect_cart()).identify() ) {

	    Base_val** bases = mp->get_mg()->std_base_vect_cart() ;
	    
	    Itbl indices (2) ;
	    indices.set_etat_qcq() ;
	    for (int i=0 ; i<n_comp ; i++) {   
		indices = donne_indices(i) ;
		(c[i]->va).set_base( (*bases[indices(0)]) * 
				     (*bases[indices(1)]) ) ;
	    }
	    for (int i=0 ; i<3 ; i++)
		delete bases[i] ;
	    delete [] bases ;
	  }
	  else {
	    assert( triad->identify() == (mp->get_bvect_spher()).identify()) ;
	    Base_val** bases = mp->get_mg()->std_base_vect_spher() ;
	    
	    Itbl indices (2) ;
	    indices.set_etat_qcq() ;
	    for (int i=0 ; i<n_comp ; i++) {   
		indices = donne_indices(i) ;
		(c[i]->va).set_base( (*bases[indices(0)]) * 
				     (*bases[indices(1)]) ) ;
	    }
	    for (int i=0 ; i<3 ; i++)
		delete bases[i] ;
	    delete [] bases ;
	  }
	  break ;
	}
	   
	    default : {
		cout << 
		"Tenseur::set_std_base() : the case valence = " << valence
		 << " is not treated !" << endl ;
		abort() ;
		break ;
	    }
	}
}

 // Le cout :
ostream& operator<<(ostream& flux, const Tenseur &source ) {
    
    flux << "Valence : " << source.valence << endl ;
    if (source.get_poids() != 0.) 
      flux << "Tensor density of weight " << source.poids << endl ;

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
	    flux << "Null Tenseur. " << endl ;
	    break ;
	    }
	
	case ETATNONDEF : {
	    flux << "Undefined Tenseur. " << endl ;
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
	    cout << "Unknown case in operator<< (ostream&, const Tenseur&)" << endl ;
	    abort() ;
	    break ;
	}
    }
    
    flux << " -----------------------------------------------------" << endl ;
    return flux ;
}

void Tenseur::sauve(FILE* fd) const {
    
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
    if (poids != 0.)
      fwrite_be(&poids, sizeof(double), 1, fd) ; //poids, si pertinent
}


void Tenseur::fait_gradient () const {
    
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

	p_gradient = new Tenseur(*mp, valence+1, tipe, mp->get_bvect_cart(), 
				 metric, poids) ;
    
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

void Tenseur::fait_gradient_spher () const {
    
    assert (etat != ETATNONDEF) ;
    
    if (p_gradient_spher != 0x0)
	return ;
    else {
 
	// Construction du resultat :
	
	if (valence != 0) {
	    cout << 
	    "Tenseur::fait_gradient_spher : the valence must be zero !" 
	    << endl ; 
	    abort() ; 
	}

	p_gradient_spher = new Tenseur(*mp, 1, COV, mp->get_bvect_spher(),
				       metric, poids) ;
    
	if (etat == ETATZERO) {
	    p_gradient_spher->set_etat_zero() ;
	}
	else {
	    assert( etat == ETATQCQ ) ; 
	    p_gradient_spher->set_etat_qcq() ;

	    p_gradient_spher->set(0) = c[0]->dsdr() ;	    // d/dr 
	    p_gradient_spher->set(1) = c[0]->srdsdt() ;	    // 1/r d/dtheta
	    p_gradient_spher->set(2) = c[0]->srstdsdp() ;   // 1/(r sin(theta))
							    //		  d/dphi

	}
    }
}


void Tenseur::fait_derive_cov (const Metrique& metre, int ind) const {
    
  assert (etat != ETATNONDEF) ;
  assert (valence != 0) ;
  
  if (p_derive_cov[ind] != 0x0)
    return ;
  else {
    
    p_derive_cov[ind] = new Tenseur (gradient()) ;
    
    if ((valence != 0) && (etat != ETATZERO)) {
      
      assert( metre.gamma().get_triad() == triad ) ; 
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
	    indices.set_etat_qcq() ;
	    indices_gamma.set(0) = indices(0) ;
	    indices_gamma.set(1) = indices(i+1) ;
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



void Tenseur::fait_derive_con (const Metrique& metre, int ind) const {
    
    if (p_derive_con[ind] != 0x0)
	return ;
    else {
	// On calcul la derivee covariante :
	if (valence != 0)
	    p_derive_con[ind] = new Tenseur
		(contract(metre.con(), 1, derive_cov(metre), 0)) ;
	
    else
	p_derive_con[ind] = new Tenseur
		(contract(metre.con(), 1, gradient(), 0)) ;
    }
}

void Tenseur::fait_carre_scal (const Metrique& met, int ind) const {
    
    if (p_carre_scal[ind] != 0x0)
	return ;
    else {
	assert (valence != 0) ;   // A ne pas appeler sur un scalaire ;
       
	// On bouge tous les indices :
	Tenseur op_t(manipule(*this, met)) ;
   
	Tenseur* auxi = new Tenseur(contract(*this, 0, op_t, 0)) ;
	Tenseur* auxi_old ;
    
	// On contracte tous les indices restant :
	for (int indice=1 ; indice<valence ; indice++) {
	    auxi_old = new Tenseur(contract(*auxi, 0, valence-indice)) ;
	    delete auxi ;
	    auxi = new Tenseur(*auxi_old) ;
	    delete auxi_old ;
	}
	p_carre_scal[ind] = new Tenseur (*auxi) ;
	delete auxi ;
    }
}
    
const Tenseur& Tenseur::gradient () const {
    if (p_gradient == 0x0)
	fait_gradient() ;
    return *p_gradient ;
}

const Tenseur& Tenseur::gradient_spher() const {
    if (p_gradient_spher == 0x0)
	fait_gradient_spher() ;
    return *p_gradient_spher ;
}

const Tenseur& Tenseur::derive_cov (const Metrique& metre) const {
    
    if (valence == 0)
	return gradient() ;
    else {
	set_dependance(metre) ;
	int j = get_place_met(metre) ;
	assert(j!=-1) ;
	if (p_derive_cov[j] == 0x0)
	    fait_derive_cov (metre,j) ;
	return *p_derive_cov[j] ;
    }
}

const Tenseur& Tenseur::derive_con (const Metrique& metre) const {
    set_dependance(metre) ;
    int j = get_place_met(metre) ;
    assert(j!=-1) ;
    if (p_derive_con[j] == 0x0)
	fait_derive_con (metre, j) ;
    return *p_derive_con[j] ;
}

const Tenseur& Tenseur::carre_scal (const Metrique& metre) const {
    set_dependance(metre) ;
    int j = get_place_met(metre) ;
    assert(j!=-1) ;
    if (p_carre_scal[j] == 0x0)
	fait_carre_scal (metre, j) ;
    return *p_carre_scal[j] ;
}
