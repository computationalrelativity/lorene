/*
 *  Methods of class Qmetrique
 *
 *   (see file qmetrique.h for documentation)
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

char qmetrique_C[] = "$Header$" ;

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
#include "qmetrique.h"
#include "utilitaires.h"

//Constructeur standard (ne fait pas grand chose) :

Qmetrique::Qmetrique (const Map& mapping, bool plate) : 
  mp(&mapping), lapse(0x0), shift(0x0), gamij(0x0), plat(plate),
  p_met_con(0x0), p_met_cov(0x0), p_met_cov_jm1(0x0), p_met_cov_jm2(0x0){
    
    etat = ETATNONDEF ;
    set_der_0x0() ;
    dependances = new (const Qtenseur* [N_DEPEND]) ;
    for (int i=0 ; i<N_DEPEND ; i++)
	dependances[i] = 0x0 ;
}

// COPY :
Qmetrique::Qmetrique (const Qmetrique& source) : mp(source.mp),
					      plat(source.plat) {
	
    if (source.lapse != 0x0)
	lapse = new Tenseur(*source.lapse) ;
    else 
	lapse = 0x0 ;
    
     if (source.shift != 0x0)
	shift = new Tenseur(*source.shift) ;
    else 
	shift = 0x0 ;
    
    if (source.gamij != 0x0)
	gamij = new Metrique(*source.gamij) ;
    else 
	gamij = 0x0 ;
    
   if (source.p_met_con != 0x0)
	p_met_con = new Qtenseur_sym(*source.p_met_con) ;
    else 
	p_met_con = 0x0 ;
    
    if (source.p_met_cov != 0x0)
	p_met_cov = new Qtenseur_sym(*source.p_met_cov) ;
    else 
	p_met_cov = 0x0 ;	

    if (source.p_met_cov_jm1 != 0x0)
	p_met_cov_jm1 = new Qtenseur_sym(*source.p_met_cov_jm1) ;
    else 
	p_met_cov_jm1 = 0x0 ;	

    if (source.p_met_cov_jm2 != 0x0)
	p_met_cov_jm2 = new Qtenseur_sym(*source.p_met_cov_jm2) ;
    else 
	p_met_cov_jm2 = 0x0 ;	

    if (source.p_gamma != 0x0)
	p_gamma = new Qtenseur_sym(*source.p_gamma) ;
    else 
	p_gamma = 0x0 ;
	
    if (source.p_gamma_jm1 != 0x0)
	p_gamma_jm1 = new Qtenseur_sym(*source.p_gamma_jm1) ;
    else 
	p_gamma_jm1 = 0x0 ;
	
    if (source.p_gamma_jm2 != 0x0)
	p_gamma_jm2 = new Qtenseur_sym(*source.p_gamma_jm2) ;
    else 
	p_gamma_jm2 = 0x0 ;	
    
    if (source.p_ricci != 0x0)
	p_ricci = new Qtenseur_sym(*source.p_ricci) ;
    else
	p_ricci = 0x0 ;
	
    if (source.p_ricci_scal != 0x0)
	p_ricci_scal = new Qtenseur(*source.p_ricci_scal) ;
    else
	p_ricci_scal = 0x0 ;
    
    if (source.p_determinant != 0x0)
	p_determinant = new Qtenseur(*source.p_determinant) ;
    else
	p_determinant = 0x0 ;
    
    dependances = new (const Qtenseur* [N_DEPEND]) ;
    for (int i=0 ; i<N_DEPEND ; i++)
	dependances[i] = 0x0 ;
    etat = source.etat ;
}

// Constructeur a partir des champs definis en formalisme 3+1
Qmetrique::Qmetrique(const Tenseur& alpha, const Tenseur& beta,
		     const Metrique& gamma, bool plate) :
  mp(alpha.get_mp()), etat(ETATQCQ), plat(plate) {
  assert(alpha.get_valence() == 0) ;  assert(alpha.get_etat() != ETATNONDEF) ;
  assert(beta.get_valence() == 1) ;   assert(beta.get_etat() != ETATNONDEF) ;
  assert(gamma.get_etat() != ETATNONDEF) ;
  assert(beta.get_type_indice(0) == CON) ;
  p_met_cov = 0x0 ;
  p_met_cov_jm1 = 0x0 ;
  p_met_cov_jm2 = 0x0 ;
  p_met_con = 0x0 ;
  lapse = new Tenseur(alpha) ;
  shift = new Tenseur(beta) ;
  gamij = new Metrique(gamma) ;
  set_der_0x0() ;
  dependances = new (const Qtenseur* [N_DEPEND]) ;
  for (int i=0 ; i<N_DEPEND ; i++)
    dependances[i] = 0x0 ;
}

// Constructeur depuis un 4-Tenseur d'ordre 2 symetrique  :
Qmetrique::Qmetrique (const Qtenseur_sym &source, bool plate) : 
  mp(source.get_mp()), etat(source.get_etat()), plat(plate) {
      
    assert (source.get_etat() != ETATNONDEF) ;
    assert (source.get_valence() == 2) ;
    
    // On regarde si on est en covariant ou contravariant ;
    int tipe = source.get_type_indice(0) ;
    assert (source.get_type_indice(1) == tipe) ;
    if (etat == ETATQCQ) {
      lapse = new Tenseur(*mp) ; 
      shift = new Tenseur(*mp, 1, CON, *source.get_triad()) ;
    }
    else {
      lapse = 0x0 ;
      shift = 0x0 ;
      gamij = 0x0 ;
    }
    if (tipe == CON) {
	p_met_con = new Qtenseur_sym (source) ;
	p_met_cov = 0x0 ;
	p_met_cov_jm1 = 0x0 ;
	p_met_cov_jm2 = 0x0 ;
	if (etat == ETATQCQ) {
	  Tenseur_sym gtmp(*mp, 2, tipe, *source.get_triad()) ;
	  lapse->set_etat_qcq() ;
	  shift->set_etat_qcq() ;
	  gtmp.set_etat_qcq() ;
	  lapse->set() = sqrt(-1/source(0,0)) ;
	  for (int i=0; i<3; i++) 
	    shift->set(i) = source(0,i+1)/source(0,0) ;
	  for (int i=0; i<3; i++) 
	    for(int j=i; j<3; j++)
	      gtmp.set(i,j) = source(i+1,j+1) 
		- source(0,i+1)*source(0,j+1)/source(0,0) ;
	  gamij = new Metrique(gtmp) ;
	}
    }
    else {
	p_met_cov = new Qtenseur_sym (source) ;
	p_met_cov_jm1 = 0x0 ;
	p_met_cov_jm2 = 0x0 ;
	p_met_con = 0x0 ;
	if (etat == ETATQCQ) {
	  Tenseur_sym gtmp(*mp, 2, tipe, *source.get_triad()) ;
	  gtmp.set_etat_qcq() ;
	  for (int i=0; i<3; i++) 
	    for(int j=i; j<3; j++)
	      gtmp.set(i,j) = source(i+1,j+1) ;
	  gamij = new Metrique(gtmp) ;
	  Tenseur tshift(*mp, 1, COV, *source.get_triad()) ;
	  tshift.set_etat_qcq() ;
	  for (int i=0; i<3; i++)
	    tshift.set(i) = source(0,i+1) ;
	  *shift = contract(gamij->con(), 1, tshift, 0) ;
	  lapse->set_etat_qcq() ;
	  lapse->set() = sqrt(contract(tshift, 0, *shift, 0)() - source(0,0)) ;
	}
    }
    if (etat == ETATQCQ) {
      shift->set_std_base() ;
      lapse->set_std_base() ;
      gamij->set_std_base() ;
    }
    set_der_0x0() ;
    dependances = new (const Qtenseur* [N_DEPEND]) ;
    for (int i=0 ; i<N_DEPEND ; i++)
	dependances[i] = 0x0 ;
}

// From a file and a mapping :
Qmetrique::Qmetrique (const Map& mapping, const Base_vect& triad,
		    FILE* fd) : mp(&mapping){
    
    fread_be (&etat, sizeof(int), 1, fd) ;
    int plate ;
    fread_be (&plate, sizeof(int), 1, fd) ;
    plat = plate ;
    
    p_met_cov = 0x0 ;
    p_met_cov_jm1 = 0x0 ;
    p_met_cov_jm2 = 0x0 ;
    p_met_con = 0x0 ;
    lapse = 0x0 ;
    shift = 0x0 ;
    gamij = 0x0 ;
    
    set_der_0x0() ;
    
    if (etat == ETATQCQ) {
      lapse = new Tenseur(mapping, fd) ;
      shift = new Tenseur(mapping, triad, fd) ;
      gamij = new Metrique(mapping, triad, fd) ;
    }
   dependances = new (const Qtenseur* [N_DEPEND]) ;
   for (int i=0 ; i<N_DEPEND ; i++)
	dependances[i] = 0x0 ; 
}

//Destructor :
Qmetrique::~Qmetrique() {
    del_dependances() ;
    delete [] dependances ;
    del_t() ;
}

//Destructeur logique
void Qmetrique::del_t() {
    
    if (p_met_con != 0x0) {
	delete p_met_con ;
	p_met_con = 0x0 ;
	}
    
    if (p_met_cov != 0x0) {
	delete p_met_cov ;
	p_met_cov = 0x0 ;
	}
    if (p_met_cov_jm1 != 0x0) {
	delete p_met_cov_jm1 ;
	p_met_cov_jm1 = 0x0 ;
	}
    if (p_met_cov_jm2 != 0x0) {
	delete p_met_cov_jm2 ;
	p_met_cov_jm2 = 0x0 ;
	}
    if (lapse != 0x0) {
      delete lapse ;
      lapse = 0x0 ;
    }
    if (shift != 0x0) {
      delete shift ;
      shift = 0x0 ;
    }
    if (gamij != 0x0) {
      delete gamij ;
      gamij = 0x0 ;
    }
    
    del_deriv() ;
    etat = ETATNONDEF ;
}

void Qmetrique::del_deriv() {
    if (p_gamma != 0x0)
	delete p_gamma ;
	
    if (p_gamma_jm1 != 0x0)
	delete p_gamma_jm1 ;
	
    if (p_gamma_jm2 != 0x0)
	delete p_gamma_jm2 ;
	
    if (p_ricci != 0x0)
	delete p_ricci ;
    
    if (p_ricci_scal != 0x0)
	delete p_ricci_scal ;
    
    if (p_determinant != 0x0)
	delete p_determinant ;
    
    set_der_0x0() ;
}

void Qmetrique::set_der_0x0() {
    p_gamma = 0x0 ;
    p_gamma_jm1 = 0x0 ;
    p_gamma_jm2 = 0x0 ;
    p_ricci = 0x0 ;
    p_ricci_scal = 0x0 ;
    p_determinant = 0x0 ;
}

void Qmetrique::del_dependances() {
    for (int i=0 ; i<N_DEPEND ; i++)
	if (dependances[i] != 0x0) {
	  dependances[i]->del_derive_met() ;
	  dependances[i] = 0x0 ;
	}
}

void Qmetrique::set_etat_nondef() {
    
    del_t() ;
    del_dependances() ;
    etat = ETATNONDEF ;   
}

void Qmetrique::set_etat_zero() {
    
    del_t() ;
    del_dependances() ;
    etat = ETATZERO ;
}


void Qmetrique::set_etat_con_qcq() {
    
    del_dependances() ;
    if (p_met_con == 0x0)
	p_met_con = new Qtenseur_sym (*mp, 2, CON, mp->get_bvect_cart()) ;
    
    if (p_met_cov != 0x0) {
	delete p_met_cov ;
	p_met_cov = 0x0 ;
    }
    
    del_deriv() ;
    etat = ETATQCQ ;
}


void Qmetrique::set_etat_cov_qcq() {
    
    del_dependances() ;
    if (p_met_cov == 0x0)
	p_met_cov = new Qtenseur_sym(*mp, 2, COV, mp->get_bvect_cart()) ;
    
    if (p_met_con != 0x0) {
	delete p_met_con ;
	p_met_con = 0x0 ;
    }
    
    del_deriv() ;
    etat = ETATQCQ ;
}


//AFFECTATIONS :

void Qmetrique::operator= (const Qmetrique& source) {
    
    assert(source.etat != ETATNONDEF) ;
    // Sur la meme grille ?
    assert (source.mp == mp) ;
    
    del_dependances() ;
    del_t() ;
   
    plat = source.plat ;
    if (source.etat == ETATZERO)
	set_etat_zero() ;
    else {
      if (source.lapse != 0x0)
	lapse = new Tenseur(*source.lapse) ;
      if (source.shift != 0x0)
	shift = new Tenseur(*source.shift) ;
      if (source.gamij != 0x0)
	gamij = new Metrique(*source.gamij) ;
      if (source.p_met_con != 0x0)
	p_met_con = new Qtenseur_sym (*source.p_met_con) ;
      if (source.p_met_cov != 0x0)
	p_met_cov = new Qtenseur_sym (*source.p_met_cov) ;
      if (source.p_met_cov_jm1 != 0x0)
	p_met_cov_jm1 = new Qtenseur_sym (*source.p_met_cov_jm1) ;
      if (source.p_met_cov_jm2 != 0x0)
	p_met_cov_jm2 = new Qtenseur_sym (*source.p_met_cov_jm2) ;
      if (source.p_gamma != 0x0)
	    p_gamma = new Qtenseur_sym (*source.p_gamma) ;
      if (source.p_gamma_jm1 != 0x0)
	    p_gamma_jm1 = new Qtenseur_sym (*source.p_gamma_jm1) ;
      if (source.p_gamma_jm2 != 0x0)
	    p_gamma_jm2 = new Qtenseur_sym (*source.p_gamma_jm2) ;
      if (source.p_ricci != 0x0)
	p_ricci = new Qtenseur_sym(*source.p_ricci) ;
      if (source.p_ricci_scal != 0x0)
	p_ricci_scal = new Qtenseur(*source.p_ricci_scal) ;
      if (source.p_determinant != 0x0)
	p_determinant = new Qtenseur(*source.p_determinant) ;
      etat = ETATQCQ ;
    }
}


void Qmetrique::operator= (const Qtenseur_sym& source) {
    
    assert (source.get_etat() != ETATNONDEF) ;
    // Sur la meme grille ?
    assert (source.get_mp() == mp) ;
    
    etat = source.get_etat() ;
    del_dependances() ;
    del_t() ;
    
    assert (source.get_valence() == 2) ;
    int tipe = source.get_type_indice(0) ;
    assert (source.get_type_indice(1) == tipe) ;
    
    if (etat == ETATQCQ) {
      lapse = new Tenseur(*mp) ; 
      shift = new Tenseur(*mp, 1, CON, *source.get_triad()) ;
    }
    if (tipe == COV)
	p_met_cov = new Qtenseur_sym (source) ;
	if (etat == ETATQCQ) {
	  lapse->set_etat_qcq() ;
	  Tenseur_sym gtmp(*mp, 2, tipe, *source.get_triad()) ;
	  gtmp.set_etat_qcq() ;
	  for (int i=0; i<3; i++) 
	    for(int j=i; j<3; j++)
	      gtmp.set(i,j) = source(i,j) ;
	  gamij = new Metrique(gtmp) ;
	  Tenseur tshift(*mp, 1, COV, *source.get_triad()) ;
	  tshift.set_etat_qcq() ;
	  for (int i=0; i<3; i++)
	    tshift.set(i) = source(0,i) ;
	  *shift = contract(gamij->con(), 1, tshift, 0) ;
	  lapse->set() = sqrt(contract(tshift, 0, *shift, 0)() - source(0,0)) ;
	}
    else
	p_met_con = new Qtenseur_sym (source) ;
	if (etat == ETATQCQ) {
	  Tenseur_sym gtmp(*mp, 2, tipe, *source.get_triad()) ;
	  lapse->set_etat_qcq() ;
	  shift->set_etat_qcq() ;
	  gtmp.set_etat_qcq() ;
	  lapse->set() = sqrt(-1/source(0,0)) ;
	  for (int i=0; i<3; i++) 
	    shift->set(i) = source(0,i)/source(0,0) ;
	  for (int i=0; i<3; i++) 
	    for(int j=i; j<3; j++)
	      gtmp.set(i,j) = source(i,j) 
		- source(0,i)*source(0,j)/source(0,0) ;
	  gamij = new Metrique(gtmp) ;
	}
}


// Le cout :
ostream& operator<<(ostream& flux, const Qmetrique & source) {
    
    switch (source.etat) {
	
	case ETATNONDEF : {
	    cout << "Undefined 4-metric in operator << ." << endl ;
	    break ;
	    }
	case ETATZERO : {
	    cout << "Null 4-metric." << endl ;
	    break ;
	   }
	   
	case ETATQCQ : {
	  cout << "4D metric" << endl ;
	if (source.plat) cout << "Flat 4-metric" << endl ; 
	if (source.p_met_con != 0x0) {
	    cout << "CONTRA-variant representation : " << endl ;
	    cout << *source.p_met_con << endl ;
	    cout << "-------------------------------------------" << endl ;
	    }
	else 
	    cout << "CONTRA-variant representation unknown : " << endl ;
	
	if (source.p_met_cov != 0x0) {
	    cout << "CO-variant representation : " << endl ;
	    cout << *source.p_met_cov << endl ;
	    cout << "-------------------------------------------" << endl ;
	}
	else 
	    cout << "CO-variant representation unknown : " << endl ;
	
	if (source.p_gamma == 0x0)
	    cout << "Christoffel unknown." << endl ;
	else
	    cout << "Christoffel known." << endl ;
	
	if (source.p_ricci == 0x0)
	    cout << "Ricci unknown." << endl ;
	else
	    cout << "Ricci known." << endl ;
	
	if (source.p_ricci_scal == 0x0)
	    cout << "Ricci scalar unknown." << endl ;
	else
	    cout << "Ricci scalar known." << endl ;
	
	if (source.p_determinant == 0x0)
	    cout << "determinant unknown." << endl ;
	else
	    cout << "determinant known." << endl ;
	break ;
	}
	default : {
	    abort() ;
	    break ;
	    }
    }
    return flux ;
}

void Qmetrique::sauve(FILE* fd) const {

    fwrite_be(&etat, sizeof(int), 1, fd) ;
    int plate = plat ;
    fwrite_be(&plate, sizeof(int), 1, fd) ;
    
    if (etat == ETATQCQ) {
      lapse->sauve(fd) ;
      shift->sauve(fd) ;
      gamij->sauve(fd) ;
    }
      
}


// Gestion des bases spectrales :
void Qmetrique:: set_std_base() {
    
    del_deriv() ;   
    if (p_met_con != 0x0)
	p_met_con->set_std_base() ;
    
    if (p_met_cov != 0x0)
	p_met_cov->set_std_base() ;
}

//Avancee dans le temps :
void Qmetrique::avance_temps(const Tenseur& alpha, const Tenseur& beta,
			     const Metrique& gama, double dt) {
  assert(alpha.get_valence() == 0) ;  assert(alpha.get_etat() != ETATNONDEF) ;
  assert(beta.get_valence() == 1) ;   assert(beta.get_etat() != ETATNONDEF) ;
  assert(gama.get_etat() != ETATNONDEF) ;
  assert(beta.get_type_indice(0) == CON) ;
  if (p_gamma_jm2 != 0x0) 
    delete p_gamma_jm2 ;
  p_gamma_jm2 = p_gamma_jm1 ;
  if (p_gamma == 0x0) fait_gamma(dt) ;
  p_gamma_jm1 = p_gamma ;
  p_gamma = 0x0 ;

  if (p_met_cov_jm2 != 0x0) 
    delete p_met_cov_jm2 ;
  p_met_cov_jm2 = p_met_cov_jm1 ;
  if (p_met_cov == 0x0) fait_cov() ;
  p_met_cov_jm1 = p_met_cov ;
  p_met_cov = 0x0 ;
  if (p_met_con != 0x0) {
    delete p_met_con ;
    p_met_con = 0x0 ;
  }
  *lapse = alpha ;
  *shift = beta ;
  *gamij = gama ;
  del_dependances() ;
}
  


// LES ROUTINES D'INVERSION ...
void Qmetrique::fait_con() const {
    
    if (p_met_con != 0x0)
	return ;
    else {
      assert(etat == ETATQCQ) ;
      if ((lapse == 0x0)||(shift == 0x0)||(gamij == 0x0)) {
	cout << "Covariant representation unknown. " << endl ;
	abort() ;
      }
      
      else {
	p_met_con = new Qtenseur_sym (*mp, 2, CON, *shift->get_triad()) ; 
	p_met_con->set_etat_qcq() ;
	p_met_con->set(0,0) = -1/((*lapse)()*(*lapse)()) ;
	for (int i=1; i<4; i++)
	  p_met_con->set(0,i) = (*shift)(i-1)*(*p_met_con)(0,0) ;
	for (int i=1; i<4; i++) 
	  for (int j=i; j<4; j++) 
	    p_met_con->set(i,j) = gamij->con()(i-1,j-1) 
	      + (*shift)(i-1)*(*shift)(j-1)* (*p_met_con)(0,0) ;
	p_met_con->set_std_base() ;
      }
    }
    
}

void Qmetrique::fait_cov() const {
    
    if (p_met_cov != 0x0)
	return ;
    else {
      if (etat == ETATZERO) { //##!! Defined on cartesian components
	p_met_cov = new Qtenseur_sym (*mp, 2, COV, mp->get_bvect_cart()) ;
	p_met_cov->set_etat_zero() ;
      }
      else 
	assert (etat == ETATQCQ) ;
	if ((lapse == 0x0)||(shift == 0x0)||(gamij == 0x0)) {
	  cout << "Contravariant representation unknown. " << endl ;
	  abort() ;
	}
      
	else {
	  p_met_cov = new Qtenseur_sym(*mp, 2, COV, *shift->get_triad()) ;
	  p_met_cov->set_etat_qcq() ;
	  Tenseur betad(contract(*shift, 0, gamij->cov(), 0)) ;
	  p_met_cov->set(0,0) = contract(*shift, 0, betad, 0)() 
	    - (*lapse)()*(*lapse)() ;
	  for (int i=1; i<4; i++) {
	    p_met_cov->set(0,i) = - betad(i-1) ;
	    for (int j=i; j<4; j++) 
	      p_met_cov->set(i,j) = gamij->cov()(i-1,j-1) ;
	  }
	  p_met_cov->set_std_base() ;
	}
    }
}




// Le calcul des Christoffel, cas general :
void Qmetrique::fait_gamma(double dt) const {
    
  assert (etat != ETATNONDEF) ;
  if (p_gamma != 0x0)
    return ;
	
  else    // Calcul a faire :
    {
      Itbl tipe (3) ;
      tipe.set_etat_qcq() ;
      tipe.set(0) = CON ; tipe.set(1) = COV ; tipe.set(2) = COV ; 
      p_gamma = new Qtenseur_sym (*mp, 3, tipe, mp->get_bvect_cart() ) ;
      bool cart = cov().get_triad()->identify() ==
	(mp->get_bvect_cart()).identify() ;

      if ( (etat == ETATZERO) || (plat && cart) )
	p_gamma->set_etat_zero() ;
      else {
	p_gamma->set_etat_qcq() ;
	Param para ;
	para.add_double(dt) ;
	if (p_met_cov_jm1 != 0x0) {
	  para.add_qtenseur_mod(*p_met_cov_jm1) ;
	  if (p_met_cov_jm2 != 0x0)
	    para.add_qtenseur_mod(*p_met_cov_jm2, 1) ;
	}
	Qtenseur t1 (contract(con(), 1, cov().gradient(para), 2)) ;
	Qtenseur t2 (contract(con(), 1, cov().gradient(para), 0)) ;

	Cmp auxi(mp) ;
	    
	// Boucle sur les composantes :
	for (int i=0 ; i<4 ; i++)
	  for (int j=0 ; j<4 ; j++)
	    for (int k=j ; k<4 ; k++) {
	      auxi = 0.5*(t1(i, j, k)+t1(i, k, j)-t2(i, j, k)) ;
	      p_gamma->set(i, j, k) = auxi ;
	    }
      }
    }
}


// Calcul de ricci :
void Qmetrique::fait_ricci(double dt) const {
    
    assert(etat != ETATNONDEF) ;
    if (p_ricci != 0x0)
	return ;
	
    else {
	p_ricci = new Qtenseur_sym (*mp, 2, COV, mp->get_bvect_cart() ) ;
	if ( (etat == ETATZERO) || (plat) )	     
	    p_ricci->set_etat_zero() ;
	else {
	    Param para ; 
	    p_ricci->set_etat_qcq() ;
	    para.add_double(dt) ;
	    if (p_gamma_jm1 != 0x0) {
	      para.add_qtenseur_mod(*p_gamma_jm1) ;
	      if (p_gamma_jm2 != 0x0)
		para.add_qtenseur_mod(*p_gamma_jm2, 1) ;
	    }

	    Qtenseur_sym grad(gamma(dt).gradient(para)) ;
	    
	    Qtenseur_sym T1 (contract(grad, 0, 1)) ;
	    
	    Qtenseur_sym T2 (contract(grad, 1, 2)) ;
	    
	    Qtenseur auxi_un (contract(gamma(dt), 0, 2)) ;
	    Qtenseur_sym T3 (contract(auxi_un, 0, gamma(dt), 0)) ;
	
	    Qtenseur auxi_deux (contract(gamma(dt), 1, gamma(dt), 0)) ;
	    Qtenseur_sym T4 (contract(auxi_deux, 0, 3)) ;
	    
	    *p_ricci = T1-T2+T3-T4 ;
	}
    }
}

// Calcul du scalaire de ricci :
void Qmetrique::fait_ricci_scal(double dt) const {
    
    assert(etat != ETATNONDEF) ;
    if (p_ricci_scal != 0x0)
	return ;
	
    else {
	
	p_ricci_scal = new Qtenseur(*mp) ;	    // Il s'agit d'un scalaire ...
	if ( (etat == ETATZERO) || (plat) )
	    p_ricci_scal->set_etat_zero() ;
	else {
	    
	    p_ricci_scal->set_etat_qcq() ;
	    Qtenseur auxi(contract(con(), 1, ricci(dt), 1)) ;
	    *p_ricci_scal = contract(auxi, 0, 1) ;
	}
    }
}


// Calcul du determinant :
void Qmetrique::fait_determinant() const {
    
  assert(etat != ETATNONDEF) ;
  if (p_determinant != 0x0)
    return ;
  
  else {
    
    p_determinant = new Qtenseur(*mp) ;	   
    if (etat == ETATZERO)
      p_determinant->set_etat_zero() ;
    else {
      assert( (lapse != 0x0) && (gamij != 0x0) );
      p_determinant->set_etat_qcq() ;
      
      p_determinant->set() = -(*lapse)()*(*lapse)()*gamij->determinant()() ;
    }
  }
}

const Qtenseur_sym& Qmetrique::con() const{
    if (p_met_con == 0x0)
	fait_con() ;
    return *p_met_con ;
}

const Qtenseur_sym& Qmetrique::cov() const{
    if (p_met_cov == 0x0)
	fait_cov() ;
    return *p_met_cov ;
}

const Qtenseur_sym& Qmetrique::gamma(double dt) const{
    if (p_gamma == 0x0)
	fait_gamma(dt) ;
    return *p_gamma ;
}

const Qtenseur_sym& Qmetrique::ricci(double dt) const{
    if (p_ricci == 0x0)
	fait_ricci(dt) ;
    return *p_ricci ;
}

const Qtenseur& Qmetrique::ricci_scal(double dt) const{
    if (p_ricci_scal == 0x0)
	fait_ricci_scal(dt) ;
    return *p_ricci_scal ;
}

const Qtenseur& Qmetrique::determinant() const{
    if (p_determinant == 0x0)
	fait_determinant() ;
    return *p_determinant ;
}
