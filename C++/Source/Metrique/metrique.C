/*
 *  Methods of class Metrique
 *
 *   (see file metrique.h for documentation)
 *
 */

/*
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


char metrique_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.2  2000/02/09  19:31:42  eric
 * Ajout de l'argument triad dans le constructeur par lecture de fichier.
 * La triade de decomposition doit desormais figurer en argument des
 * constructeurs de Tenseur.
 *
 * Revision 2.1  2000/01/10  17:21:30  eric
 * Modif des #include
 *
 * Revision 2.0  1999/12/02  17:19:15  phil
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
#include "cmp.h"
#include "metrique.h"

//Constructeur standard (ne fait pas grand chose) :

Metrique::Metrique (const Map& mapping) : mp(&mapping) {
    
    p_met_con = 0x0 ;
    p_met_cov = 0x0 ;
    set_der_0x0() ;
    etat = ETATNONDEF ;
    
    dependances = new (const Tenseur* [N_DEPEND]) ;
    for (int i=0 ; i<N_DEPEND ; i++)
	dependances[i] = 0x0 ;
}

// COPY :
Metrique::Metrique (const Metrique& source) : mp(source.mp) {
	
    if (source.p_met_con != 0x0)
	p_met_con = new Tenseur_sym(*source.p_met_con) ;
    else 
	p_met_con = 0x0 ;
    
    if (source.p_met_cov != 0x0)
	p_met_cov = new Tenseur_sym(*source.p_met_cov) ;
    else 
	p_met_cov = 0x0 ;	

    if (source.p_gamma != 0x0)
	p_gamma = new Tenseur_sym(*source.p_gamma) ;
    else 
	p_gamma = 0x0 ;
	
    
    if (source.p_ricci != 0x0)
	p_ricci = new Tenseur_sym(*source.p_ricci) ;
    else
	p_ricci = 0x0 ;
	
    if (source.p_ricci_scal != 0x0)
	p_ricci_scal = new Tenseur(*source.p_ricci_scal) ;
    else
	p_ricci_scal = 0x0 ;
    
    dependances = new (const Tenseur* [N_DEPEND]) ;
    for (int i=0 ; i<N_DEPEND ; i++)
	dependances[i] = 0x0 ;
    etat = source.etat ;
}


// Constructeur from Tensor d'ordre 2 symetrique  :
Metrique::Metrique (const Tenseur_sym &source) : mp(source.get_mp()) {
      
    assert (source.get_etat() != ETATNONDEF) ;
    assert (source.get_valence() == 2) ;
    
    // On regarde si on est en covariant ou contravariant ;
    int tipe = source.get_type_indice(0) ;
    assert (source.get_type_indice(1) == tipe) ;
    
    if (tipe == CON) {
	p_met_con = new Tenseur_sym (source) ;
	p_met_cov = 0x0 ;
	}
    else {
	p_met_cov = new Tenseur_sym (source) ;
	p_met_con = 0x0 ;
	}
	
    set_der_0x0() ;
    dependances = new (const Tenseur* [N_DEPEND]) ;
    for (int i=0 ; i<N_DEPEND ; i++)
	dependances[i] = 0x0 ;
    etat = source.get_etat() ;
}

// From a file and a mapping :
Metrique::Metrique (const Map& mapping, const Base_vect& triad,
		    FILE* fd) : mp(&mapping){
    
    fread (&etat, sizeof(int), 1, fd) ;
    
    p_met_cov = 0x0 ;
    p_met_con = 0x0 ;
    set_der_0x0() ;
    
    if (etat == ETATQCQ) {
	int indic ;
	fread (&indic, sizeof(int), 1, fd) ;
	switch (indic) {
	    case COV : {
		p_met_cov = new Tenseur_sym(mapping, triad, fd) ;
		break ;
		}
	    case CON : {
		p_met_con = new Tenseur_sym(mapping, triad, fd) ;
		break ;
		}
	    default :
		break ;
	}
    }
   dependances = new (const Tenseur* [N_DEPEND]) ;
   for (int i=0 ; i<N_DEPEND ; i++)
	dependances[i] = 0x0 ; 
}

//Destructor :
Metrique::~Metrique() {
    del_t() ;
    delete [] dependances ;
}

//Destructeur logique
void Metrique::del_t() {
    
    if (p_met_con != 0x0) {
	delete p_met_con ;
	p_met_con = 0x0 ;
	}
    
    if (p_met_cov != 0x0) {
	delete p_met_cov ;
	p_met_cov = 0x0 ;
	}
    
    del_gamma() ;
    etat = ETATNONDEF ;
}

void Metrique::del_gamma() {
    if (p_gamma != 0x0)
	delete p_gamma ;
	
    if (p_ricci != 0x0)
	delete p_ricci ;
    
    if (p_ricci_scal != 0x0)
	delete p_ricci_scal ;
    
    set_der_0x0() ;
}

void Metrique::set_der_0x0() {
    p_gamma = 0x0 ;
    p_ricci = 0x0 ;
    p_ricci_scal = 0x0 ;
}

void Metrique::del_dependances() {
    for (int i=0 ; i<N_DEPEND ; i++)
	if (dependances[i] != 0x0) {
	    dependances[i]->del_derive_met() ;
	    dependances[i] = 0x0 ;
	}
}

void Metrique::set_etat_nondef() {
    
    del_t() ;
    del_dependances() ;
    etat = ETATNONDEF ;   
}

void Metrique::set_etat_zero() {
    
    del_t() ;
    del_dependances() ;
    etat = ETATZERO ;
}


void Metrique::set_etat_con_qcq() {
    
    del_dependances() ;
    if (p_met_con == 0x0)
	p_met_con = new Tenseur_sym (*mp, 2, CON, mp->get_bvect_cart()) ;
    
    if (p_met_cov != 0x0) {
	delete p_met_cov ;
	p_met_cov = 0x0 ;
    }
    
    del_gamma() ;
    etat = ETATQCQ ;
}


void Metrique::set_etat_cov_qcq() {
    
    del_dependances() ;
    if (p_met_cov == 0x0)
	p_met_cov = new Tenseur_sym(*mp, 2, COV, mp->get_bvect_cart()) ;
    
    if (p_met_con != 0x0) {
	delete p_met_con ;
	p_met_con = 0x0 ;
    }
    
    del_gamma() ;
    etat = ETATQCQ ;
}


//AFFECTATIONS :

void Metrique::operator= (const Metrique& source) {
    
    assert(source.etat != ETATNONDEF) ;
    // Sur la meme grille ?
    assert (source.mp == mp) ;
    
    del_dependances() ;
    del_t() ;
   
    if (source.etat == ETATZERO)
	set_etat_zero() ;
    else {
	if (source.p_met_con != 0x0)
	    p_met_con = new Tenseur_sym (*source.p_met_con) ;
	if (source.p_met_cov != 0x0)
	    p_met_cov = new Tenseur_sym (*source.p_met_cov) ;
	if (source.p_gamma != 0x0)
	    p_gamma = new Tenseur_sym (*source.p_gamma) ;
	if (source.p_ricci != 0x0)
	    p_ricci = new Tenseur_sym(*source.p_ricci) ;
	if (source.p_ricci_scal != 0x0)
	    p_ricci_scal = new Tenseur(*source.p_ricci_scal) ;
	etat = ETATQCQ ;
    }
}


void Metrique::operator= (const Tenseur_sym& source) {
    
    assert (source.get_etat() != ETATNONDEF) ;
    // Sur la meme grille ?
    assert (source.get_mp() == mp) ;
    
    del_dependances() ;
    del_t() ;
    
    assert (source.get_valence() == 2) ;
    int tipe = source.get_type_indice(0) ;
    assert (source.get_type_indice(1) == tipe) ;
    
    if (tipe == COV)
	p_met_cov = new Tenseur_sym (source) ;
    else
	p_met_con = new Tenseur_sym (source) ;
    
    etat = source.get_etat() ;
}


// Le cout :
ostream& operator<<(ostream& flux, const Metrique & source) {
    
    switch (source.etat) {
	
	case ETATNONDEF : {
	    cout << "Undefined metric in operator << ." << endl ;
	    break ;
	    }
	case ETATZERO : {
	    cout << "Nul metric." << endl ;
	    break ;
	   }
	   
	case ETATQCQ : {
	
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
	break ;
	}
	default : {
	    abort() ;
	    break ;
	    }
    }
    return flux ;
}

void Metrique::sauve(FILE* fd) const {

    fwrite(&etat, sizeof(int), 1, fd) ;
    
    if (etat == ETATQCQ) {
    
    // Dis ce que l'on a sauve ...
    int indic ;
    if (p_met_cov != 0x0)
	indic = COV ;
    else if (p_met_con != 0x0)
	    indic = CON ;
	 else indic = 0 ;
    fwrite(&indic, sizeof(int), 1, fd) ;
    switch (indic) {
	case COV : {
	    p_met_cov->sauve(fd) ;
	    break ;
	    }
	case CON : {
	    p_met_con->sauve(fd) ;
	    break ;
	    }
	default : {
	    break ;
	    }
    }
  }
}


// Gestion des bases spectrales :
void Metrique:: set_std_base() {
    
    del_gamma() ;   
    if (p_met_con != 0x0)
	p_met_con->set_std_base() ;
    
    if (p_met_cov != 0x0)
	p_met_cov->set_std_base() ;
}



// LES ROUTINES D'INVERSION ...
void Metrique::fait_con() const {
    
    if (p_met_con != 0x0)
	return ;
    else
	if (p_met_cov == 0x0) {
	    cout << "Covariant representation unknown. " << endl ;
	    abort() ;
	    }
	    
	 else
	    p_met_con = new Tenseur_sym (fait_inverse(*p_met_cov)) ;    
}

void Metrique::fait_cov() const {
    
    if (p_met_cov != 0x0)
	return ;
    else
	if (p_met_con == 0x0) {
	    cout << "Contravariant representation unknown. " << endl ;
	    abort() ;
	    }
	    
	 else
	    p_met_cov = new Tenseur_sym(fait_inverse(*p_met_con)) ;    
    
}




// Le calcul des Christoffel, cas general :
void Metrique::fait_gamma() const {
    
    assert (etat != ETATNONDEF) ;
    if (p_gamma != 0x0)
	return ;
	
    else    // Calcul a faire :
    {
	Itbl tipe (3) ;
	tipe.set_etat_qcq() ;
	tipe.set(0) = CON ; tipe.set(1) = COV ; tipe.set(2) = COV ; 
	p_gamma = new Tenseur_sym (*mp, 3, tipe, mp->get_bvect_cart() ) ;
	
	if (etat == ETATZERO)
	    p_gamma->set_etat_zero() ;
	else {
	    p_gamma->set_etat_qcq() ;
	   
	    Tenseur t1 (contract(con(), 1, cov().gradient(), 2)) ;
	    Tenseur t2 (contract(con(), 1, cov().gradient(), 0)) ;

	    Cmp auxi(mp) ;
	    
	    // Boucle sur les composantes :
	    for (int i=0 ; i<3 ; i++)
		for (int j=0 ; j<3 ; j++)
		    for (int k=j ; k<3 ; k++) {
			auxi = 0.5*(t1(i, j, k)+t1(i, k, j)-t2(i, j, k)) ;
			p_gamma->set(i, j, k) = auxi ;
			}
	   }
    }
}


// Calcul de ricci :
void Metrique::fait_ricci() const {
    
    assert(etat != ETATNONDEF) ;
    if (p_ricci != 0x0)
	return ;
	
    else {
	p_ricci = new Tenseur_sym (*mp, 2, COV, mp->get_bvect_cart() ) ;
	if (etat == ETATZERO)
	    p_ricci->set_etat_zero() ;
	else {
	    
	    p_ricci->set_etat_qcq() ;
	    
	    Tenseur_sym copie_gamma(gamma()) ;
	    copie_gamma.dec2_dzpuis() ;
	    Tenseur_sym grad(copie_gamma.gradient()) ;
	    grad.inc2_dzpuis() ;
	    
	    Tenseur_sym T1 (contract(grad, 0, 1)) ;
	    
	    Tenseur_sym T2 (contract(grad, 1, 2)) ;
	    
	    Tenseur auxi_un (contract(gamma(), 0, 2)) ;
	    Tenseur_sym T3 (contract(auxi_un, 0, gamma(), 0)) ;
	
	    Tenseur auxi_deux (contract(gamma(), 1, gamma(), 0)) ;
	    Tenseur_sym T4 (contract(auxi_deux, 0, 3)) ;
	    
	    *p_ricci = T1-T2+T3-T4 ;
	}
    }
}

// Calcul du scalaire de ricci :
void Metrique::fait_ricci_scal() const {
    
    assert(etat != ETATNONDEF) ;
    if (p_ricci_scal != 0x0)
	return ;
	
    else {
	
	p_ricci_scal = new Tenseur(*mp) ;	    // Il s'agit d'un scalaire ...
	if (etat == ETATZERO)
	    p_ricci_scal->set_etat_zero() ;
	else {
	    
	    p_ricci_scal->set_etat_qcq() ;
	    Tenseur auxi(contract(con(), 1, ricci(), 1)) ;
	    *p_ricci_scal = contract(auxi, 0, 1) ;
	}
    }
}


Tenseur_sym fait_inverse (const Tenseur_sym& s) {
    
    assert (s.get_etat() == ETATQCQ) ;
    assert (s.get_valence() == 2) ;
    assert (s.get_type_indice(0) == s.get_type_indice(1)) ;
  
    //Le resultat :
    Tenseur_sym res (*s.get_mp(), 2, -s.get_type_indice(0), *(s.get_triad()) ) ;
    res.set_etat_qcq() ;
    
    // le determinant :
    Cmp determ1(s.get_mp()) ;
    determ1 = double(1)/(s(0, 0)*s(1, 1)*s(2, 2)+s(0, 1)*s(1, 2)*s(0, 2)
	    +s(0, 2)*s(0, 1)*s(1, 2)-s(0, 2)*s(1, 1)*s(0, 2)
	    -s(1, 2)*s(1, 2)*s(0, 0)-s(2, 2)*s(0, 1)*s(0, 1)) ;
    
    int sgn ;	// Le signe du co-facteur ...
    int l_up, l_down, c_left, c_right ;	    // Coordonnees du cofacteur :
    
    Cmp cofacteur(s.get_mp()) ;
    
    for (int i=0 ; i<3 ; i++) {
	sgn = 1 ;
	for (int j=i ; j<3 ; j++) {
	    
	    switch (j) {
		
		case 0 : {
		    c_left = 1 ;
		    c_right = 2 ;
		    break ;
		    }
		case 1 : {
		    c_left = 0 ;
		    c_right = 2 ;
		    break ;
		    }
		default : {
		    c_left = 0 ;
		    c_right = 1 ;
		    break ;
		    }
	    }
	    
	    switch (i) {
		
		case 0 : {
		    l_up = 1 ;
		    l_down = 2 ;
		    break ;
		    }
		case 1 : {
		    l_up = 0 ;
		    l_down = 2 ;
		    break ;
		    }
		default : {
		    l_up = 0 ;
		    l_down = 1 ;
		    break ;
		    } 
	    }
	    
	    cofacteur = sgn*(s(l_up, c_left)*s(l_down, c_right)-
			    s(l_up, c_right)*s(l_down, c_left))*determ1 ;
	    
	    res.set(i, j) = cofacteur ;
	    sgn *= -1 ;
	    }
	}
   return res ;
}

const Tenseur_sym& Metrique::con() const{
    if (p_met_con == 0x0)
	fait_con() ;
    return *p_met_con ;
}

const Tenseur_sym& Metrique::cov() const{
    if (p_met_cov == 0x0)
	fait_cov() ;
    return *p_met_cov ;
}

const Tenseur_sym& Metrique::gamma() const{
    if (p_gamma == 0x0)
	fait_gamma() ;
    return *p_gamma ;
}

const Tenseur_sym& Metrique::ricci() const{
    if (p_ricci == 0x0)
	fait_ricci() ;
    return *p_ricci ;
}

const Tenseur& Metrique::ricci_scal() const{
    if (p_ricci_scal == 0x0)
	fait_ricci_scal() ;
    return *p_ricci_scal ;
}
