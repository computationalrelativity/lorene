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

char qtenseur_operateur_C[] = "$Header$" ;

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


Qtenseur operator*(const Qtenseur& t1, const Qtenseur& t2) {
   
    assert ((t1.etat != ETATNONDEF) && (t2.etat != ETATNONDEF)) ;
    assert (t1.mp == t2.mp) ;
    
    int val_res = t1.valence + t2.valence ;
    
   // cas scalaire :
    if (val_res == 0) {
	Qtenseur scal(*t1.mp) ;
	// cas ou un des deux est nul :
	if ((t1.etat == ETATZERO) || (t2.etat == ETATZERO))
	    scal.set_etat_zero() ;
	else {
	    scal.set_etat_qcq() ;
	    scal.set() = t1() * t2() ;
	}
    return scal ;
   }
    
    else {
	
	Itbl tipe (val_res) ;
	tipe.set_etat_qcq() ;
	for (int i=0 ; i<t1.valence ; i++)
	    tipe.set(i) = t1.type_indice(i) ;
	for (int i=0 ; i<t2.valence ; i++)
	    tipe.set(i+t1.valence) = t2.type_indice(i) ;
	

	if ( (t1.valence != 0) && (t2.valence != 0) ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
	}

	const Base_vect* triad_res ; 
	if (t1.valence != 0) {
	    triad_res = t1.get_triad() ; 
	}
	else{
	    triad_res = t2.get_triad() ; 
	}

	Qtenseur res(*t1.mp, val_res, tipe, triad_res) ;


	
	if ((t1.etat == ETATZERO) || (t2.etat == ETATZERO))
	    res.set_etat_zero() ;
	else {
	    res.set_etat_qcq() ;
	    Itbl jeux_indice_t1 (t1.valence) ;
	    jeux_indice_t1.set_etat_qcq() ;
	    Itbl jeux_indice_t2 (t2.valence) ;
	    jeux_indice_t2.set_etat_qcq() ;
	    
	    for (int i=0 ; i<res.n_comp ; i++) {
		Itbl jeux_indice_res(res.donne_indices(i)) ;
		for (int j=0 ; j<t1.valence ; j++)
		    jeux_indice_t1.set(j) = jeux_indice_res(j) ;
		for (int j=0 ; j<t2.valence ; j++)
		    jeux_indice_t2.set(j) = jeux_indice_res(j+t1.valence) ;
		
		res.set(jeux_indice_res) = t1(jeux_indice_t1)*t2(jeux_indice_t2) ;
	    }
	}
	return res ;
    }
}

    //------------------------------------//
    // Tensorial product wiht desaliasing //
    //------------------------------------//
    

Qtenseur operator%(const Qtenseur& t1, const Qtenseur& t2) {
   
    assert ((t1.etat != ETATNONDEF) && (t2.etat != ETATNONDEF)) ;
    assert (t1.mp == t2.mp) ;
    
    int val_res = t1.valence + t2.valence ;
    
   // cas scalaire :
    if (val_res == 0) {
	Qtenseur scal(*t1.mp) ;
	// cas ou un des deux est nul :
	if ((t1.etat == ETATZERO) || (t2.etat == ETATZERO))
	    scal.set_etat_zero() ;
	else {
	    scal.set_etat_qcq() ;
	    scal.set() = t1() % t2() ;
	}
    return scal ;
   }
    
    else {
	
	Itbl tipe (val_res) ;
	tipe.set_etat_qcq() ;
	for (int i=0 ; i<t1.valence ; i++)
	    tipe.set(i) = t1.type_indice(i) ;
	for (int i=0 ; i<t2.valence ; i++)
	    tipe.set(i+t1.valence) = t2.type_indice(i) ;
	

	if ( (t1.valence != 0) && (t2.valence != 0) ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
	}

	const Base_vect* triad_res ; 
	if (t1.valence != 0) {
	    triad_res = t1.get_triad() ; 
	}
	else{
	    triad_res = t2.get_triad() ; 
	}

	Qtenseur res(*t1.mp, val_res, tipe, triad_res) ;


	
	if ((t1.etat == ETATZERO) || (t2.etat == ETATZERO))
	    res.set_etat_zero() ;
	else {
	    res.set_etat_qcq() ;
	    Itbl jeux_indice_t1 (t1.valence) ;
	    jeux_indice_t1.set_etat_qcq() ;
	    Itbl jeux_indice_t2 (t2.valence) ;
	    jeux_indice_t2.set_etat_qcq() ;
	    
	    for (int i=0 ; i<res.n_comp ; i++) {
		Itbl jeux_indice_res(res.donne_indices(i)) ;
		for (int j=0 ; j<t1.valence ; j++)
		    jeux_indice_t1.set(j) = jeux_indice_res(j) ;
		for (int j=0 ; j<t2.valence ; j++)
		    jeux_indice_t2.set(j) = jeux_indice_res(j+t1.valence) ;
		
		res.set(jeux_indice_res) = t1(jeux_indice_t1) % 
						    t2(jeux_indice_t2) ;
	    }
	}
	return res ;
    }
}



Qtenseur contract(const Qtenseur& source, int ind_1, int ind_2)  {
    
    
    // Les verifications :
    assert (source.etat != ETATNONDEF) ;
    assert ((ind_1 >= 0) && (ind_1 < source.valence)) ;
    assert ((ind_2 >= 0) && (ind_2 < source.valence)) ;
    assert (source.type_indice(ind_1) != source.type_indice(ind_2))  ;
 
    // On veut ind_1 < ind_2 :
    if (ind_1 > ind_2) {
	int auxi = ind_2 ;
	ind_2 = ind_1 ;
	ind_1 = auxi ;
    }
    
    // On construit le resultat :
    int val_res = source.valence - 2 ;
   
    Itbl tipe (val_res) ;
    tipe.set_etat_qcq() ;
    for (int i=0 ; i<ind_1 ; i++)
	tipe.set(i) = source.type_indice(i) ;
    for (int i=ind_1 ; i<ind_2-1 ; i++)
	tipe.set(i) = source.type_indice(i+1) ;
    for (int i = ind_2-1 ; i<val_res ; i++)
	tipe.set(i) = source.type_indice(i+2) ;
	
    Qtenseur res(*source.mp, val_res, tipe, source.triad) ;

    // Cas particulier d'une source nulle
    if (source.etat == ETATZERO) {
	res.set_etat_zero() ; 
	return res ; 
    }

    res.set_etat_qcq() ;
	
    Cmp work(source.mp) ;
    
    // Boucle sur les composantes de res :
	
    Itbl jeux_indice_source(source.valence) ;
    jeux_indice_source.set_etat_qcq() ;
	
    for (int i=0 ; i<res.n_comp ; i++) {
	Itbl jeux_indice_res (res.donne_indices(i)) ;
	for (int j=0 ; j<ind_1 ; j++)
	    jeux_indice_source.set(j) = jeux_indice_res(j) ;
	for (int j=ind_1+1 ; j<ind_2 ; j++)
	    jeux_indice_source.set(j) = jeux_indice_res(j-1) ;
	for (int j=ind_2+1 ; j<source.valence ; j++)
	    jeux_indice_source.set(j) = jeux_indice_res(j-2) ;
	    
	    
	work.set_etat_zero() ;
	for (int j=0 ; j<4 ; j++) {
	    jeux_indice_source.set(ind_1) = j ;
	    jeux_indice_source.set(ind_2) = j ;
	    work = work + source(jeux_indice_source) ;
	    }
	    
	res.set(jeux_indice_res) = work ;
	}
    return res ;
}


Qtenseur contract (const Qtenseur& t1, int ind1, const Qtenseur& t2, int ind2) {
    
    assert ((t1.etat != ETATNONDEF) && (t2.etat != ETATNONDEF)) ;
    // Verifs :
    assert ((ind1>=0) && (ind1<t1.valence)) ;
    assert ((ind2>=0) && (ind2<t2.valence)) ;
    assert (t1.mp == t2.mp) ;
    
    // Contraction possible ?
    if ( (t1.valence != 0) && (t2.valence != 0) ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    assert (t1.type_indice(ind1) != t2.type_indice(ind2)) ;
    
    int val_res = t1.valence + t2.valence - 2;
    Itbl tipe(val_res) ;
    tipe.set_etat_qcq() ;
    for (int i=0 ; i<ind1 ; i++)
	tipe.set(i) = t1.type_indice(i) ;
    for (int i=ind1 ; i<t1.valence-1 ; i++)
	tipe.set(i) = t1.type_indice(i+1) ;
    for (int i=t1.valence-1 ; i<t1.valence+ind2-1 ; i++)
	tipe.set(i) = t2.type_indice(i-t1.valence+1) ;
    for (int i = t1.valence+ind2-1 ; i<val_res ; i++)
	tipe.set(i) = t2.type_indice(i-t1.valence+2) ;
	
    const Base_vect* triad_res = (val_res == 0) ? 0x0 : t1.get_triad() ; 

    Qtenseur res(*t1.mp, val_res, tipe, triad_res) ;

    // Cas particulier ou l'un des deux tenseurs est nul
    if ( (t1.etat == ETATZERO) || (t2.etat == ETATZERO) ) {
	res.set_etat_zero() ; 
	return res ; 
    }

    res.set_etat_qcq() ;
	
    Cmp work(t1.mp) ;
    
    // Boucle sur les composantes de res :
	
    Itbl jeux_indice_t1(t1.valence) ;
    Itbl jeux_indice_t2(t2.valence) ;
    jeux_indice_t1.set_etat_qcq() ;
    jeux_indice_t2.set_etat_qcq() ;
    
    for (int i=0 ; i<res.n_comp ; i++) {
	Itbl jeux_indice_res (res.donne_indices(i)) ;
	for (int i=0 ; i<ind1 ; i++)
	    jeux_indice_t1.set(i) = jeux_indice_res(i) ;
	for (int i=ind1+1 ; i<t1.valence ; i++)
	    jeux_indice_t1.set(i) = jeux_indice_res(i-1) ;
	for (int i=0 ; i<ind2 ; i++)
	    jeux_indice_t2.set(i) = jeux_indice_res(t1.valence+i-1) ;
	for (int i=ind2+1 ; i<t2.valence ; i++)
	    jeux_indice_t2.set(i) = jeux_indice_res(t1.valence+i-2) ;
	
	    
	    
	work.set_etat_zero() ;
	for (int j=0 ; j<4 ; j++) {
	    jeux_indice_t1.set(ind1) = j ;
	    jeux_indice_t2.set(ind2) = j ;
	    work = work + t1(jeux_indice_t1)*t2(jeux_indice_t2) ;
	    }
	    
	res.set(jeux_indice_res) = work ;
	}
    return res ;
}


Qtenseur manipule(const Qtenseur& t1, const Qmetrique& met, int place) {
    
    assert (t1.etat != ETATNONDEF) ;
    assert (met.get_etat() != ETATNONDEF) ;
    
    int valen = t1.valence ;
    assert (valen != 0) ;	    // Aucun interet pour un scalaire...
    assert ((place >=0) && (place < valen)) ;
    
    Itbl tipe (valen) ;
    tipe.set_etat_qcq() ;
    tipe.set(0) = -t1.type_indice(place) ;
    for (int i=1 ; i<place+1 ; i++)
	tipe.set(i) = t1.type_indice(i-1) ;
    for (int i=place+1 ; i<valen ; i++)
	tipe.set(i) = t1.type_indice(i) ;
    
    Qtenseur auxi(*t1.mp, valen, tipe, t1.triad) ;
    
    if (t1.type_indice(place) == COV)
	auxi = contract (met.con(), 1, t1, place) ;
    else
	auxi = contract (met.cov(), 1, t1, place) ;
    
    // On doit remettre les indices a la bonne place ...
    
    for (int i=0 ; i<valen ; i++)
	tipe.set(i) = t1.type_indice(i) ;
    tipe.set(place) *= -1 ;
    
    Qtenseur res(*t1.mp, valen, tipe, t1.triad) ;
    res.set_etat_qcq() ;
    
    Itbl place_auxi(valen) ;
    place_auxi.set_etat_qcq() ;
    
    for (int i=0 ; i<res.n_comp ; i++) {
	
	Itbl place_res (res.donne_indices(i)) ;
	
	place_auxi.set(0) = place_res(place) ;
	for (int j=1 ; j<place+1 ; j++)
	    place_auxi.set(j) = place_res(j-1)  ;
	place_res.set(place) = place_auxi(0) ;
	for (int j=place+1 ; j<valen ; j++)
	     place_auxi.set(j) = place_res(j);
	
	
	res.set(place_res) = auxi(place_auxi) ;
    }
    return res ;
}

Qtenseur manipule (const Qtenseur& t1, const Qmetrique& met) {
    
    Qtenseur* auxi ;
    Qtenseur* auxi_old = new Qtenseur(t1) ;
    
    for (int i=0 ; i<t1.valence ; i++) {
	auxi = new Qtenseur(manipule(*auxi_old, met, i)) ;
	delete auxi_old ;
	auxi_old = new Qtenseur(*auxi) ;
	delete auxi ;
    }
    
    Qtenseur result(*auxi_old) ;
    delete auxi_old ;
    return result ;
}


Qtenseur flat_scalar_prod(const Qtenseur& t1, const Qtenseur& t2) {
    
    assert ((t1.etat != ETATNONDEF) && (t2.etat != ETATNONDEF)) ;
    // Verifs :
    assert (t1.mp == t2.mp) ;
    
    // Contraction possible ?
    if ( (t1.valence != 0) && (t2.valence != 0) ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    
    int val_res = t1.valence + t2.valence - 2;
    Itbl tipe(val_res) ;
    tipe.set_etat_qcq() ;

    for (int i=0 ; i<t1.valence - 1 ; i++)
	tipe.set(i) = t1.type_indice(i) ;
    for (int i = t1.valence-1 ; i<val_res ; i++)
	tipe.set(i) = t2.type_indice(i-t1.valence+2) ;
	
    Qtenseur res(*t1.mp, val_res, tipe, t1.triad) ;

    // Cas particulier ou l'un des deux Qtenseurs est nul
    if ( (t1.etat == ETATZERO) || (t2.etat == ETATZERO) ) {
	res.set_etat_zero() ; 
	return res ; 
    }

    res.set_etat_qcq() ;
	
    Cmp work(t1.mp) ;
    
    // Boucle sur les composantes de res :
	
    Itbl jeux_indice_t1(t1.valence) ;
    Itbl jeux_indice_t2(t2.valence) ;
    jeux_indice_t1.set_etat_qcq() ;
    jeux_indice_t2.set_etat_qcq() ;
    
    for (int ir=0 ; ir<res.n_comp ; ir++) {    // Boucle sur les composantes
					       // du resultat 

	// Indices du resultat correspondant a la position ir : 
	Itbl jeux_indice_res = res.donne_indices(ir) ;

	// Premiers indices correspondant dans t1 : 
	for (int i=0 ; i<t1.valence - 1 ; i++) {
	    jeux_indice_t1.set(i) = jeux_indice_res(i) ;
	}
	
	// Derniers indices correspondant dans t2 : 
	for (int i=1 ; i<t2.valence ; i++) {
	    jeux_indice_t2.set(i) = jeux_indice_res(t1.valence+i-2) ;
	}
	
	work.set_etat_zero() ;

	// Sommation sur le dernier indice de t1 et le premier de t2 : 
	
	for (int j=0 ; j<4 ; j++) {
	    jeux_indice_t1.set(t1.valence - 1) = j ;
	    jeux_indice_t2.set(0) = j ;
	    work = work + t1(jeux_indice_t1)*t2(jeux_indice_t2) ;
	    }
	    
	res.set(jeux_indice_res) = work ;
    
    }	// fin de la boucle sur les composantes du resultat
    
    return res ;
}



Qtenseur flat_scalar_prod_desal(const Qtenseur& t1, const Qtenseur& t2) {
    
    assert ((t1.etat != ETATNONDEF) && (t2.etat != ETATNONDEF)) ;
    // Verifs :
    assert (t1.mp == t2.mp) ;
    
    // Contraction possible ?
    if ( (t1.valence != 0) && (t2.valence != 0) ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    
    int val_res = t1.valence + t2.valence - 2;
    Itbl tipe(val_res) ;
    tipe.set_etat_qcq() ;

    for (int i=0 ; i<t1.valence - 1 ; i++)
	tipe.set(i) = t1.type_indice(i) ;
    for (int i = t1.valence-1 ; i<val_res ; i++)
	tipe.set(i) = t2.type_indice(i-t1.valence+2) ;
	
    Qtenseur res(*t1.mp, val_res, tipe, t1.triad) ;

    // Cas particulier ou l'un des deux Qtenseurs est nul
    if ( (t1.etat == ETATZERO) || (t2.etat == ETATZERO) ) {
	res.set_etat_zero() ; 
	return res ; 
    }

    res.set_etat_qcq() ;
	
    Cmp work(t1.mp) ;
    
    // Boucle sur les composantes de res :
	
    Itbl jeux_indice_t1(t1.valence) ;
    Itbl jeux_indice_t2(t2.valence) ;
    jeux_indice_t1.set_etat_qcq() ;
    jeux_indice_t2.set_etat_qcq() ;
    
    for (int ir=0 ; ir<res.n_comp ; ir++) {    // Boucle sur les composantes
					       // du resultat 

	// Indices du resultat correspondant a la position ir : 
	Itbl jeux_indice_res = res.donne_indices(ir) ;

	// Premiers indices correspondant dans t1 : 
	for (int i=0 ; i<t1.valence - 1 ; i++) {
	    jeux_indice_t1.set(i) = jeux_indice_res(i) ;
	}
	
	// Derniers indices correspondant dans t2 : 
	for (int i=1 ; i<t2.valence ; i++) {
	    jeux_indice_t2.set(i) = jeux_indice_res(t1.valence+i-2) ;
	}
	
	work.set_etat_zero() ;

	// Sommation sur le dernier indice de t1 et le premier de t2 : 
	
	for (int j=0 ; j<4 ; j++) {
	    jeux_indice_t1.set(t1.valence - 1) = j ;
	    jeux_indice_t2.set(0) = j ;
	    work = work + t1(jeux_indice_t1) % t2(jeux_indice_t2) ;
	    }
	    
	res.set(jeux_indice_res) = work ;
    
    }	// fin de la boucle sur les composantes du resultat
    
    return res ;
}

