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


char tenseur_operateur_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.9  2014/10/06 15:13:19  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2004/05/27 07:17:19  p_grandclement
 * Correction of some shadowed variables
 *
 * Revision 1.7  2003/06/20 14:53:38  f_limousin
 * Add the function contract_desal()
 *
 * Revision 1.6  2003/03/03 19:38:41  f_limousin
 * Suppression of an assert on a metric associated with a tensor.
 *
 * Revision 1.5  2002/10/16 14:37:14  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.4  2002/09/10 13:44:17  j_novak
 * The method "manipule" of one indice has been removed for Tenseur_sym objects
 * (the result cannot be a Tenseur_sym).
 * The method "sans_trace" now computes the traceless part of a Tenseur (or
 * Tenseur_sym) of valence 2.
 *
 * Revision 1.3  2002/09/06 14:49:25  j_novak
 * Added method lie_derive for Tenseur and Tenseur_sym.
 * Corrected various errors for derive_cov and arithmetic.
 *
 * Revision 1.2  2002/08/07 16:14:11  j_novak
 * class Tenseur can now also handle tensor densities, this should be transparent to older codes
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.11  2001/08/27  10:04:21  eric
 * Ajout de l'operator% (produit tensoriel avec desaliasing)
 *
 * Revision 2.10  2001/05/26  15:43:17  eric
 * Ajout de la fonction flat_scalar_prod_desal (desaliasage)
 *
 * Revision 2.9  2000/10/06  15:08:40  eric
 * Traitement des cas ETATZERO dans contract et flat_scal_prod
 *
 * Revision 2.8  2000/09/13  09:43:29  eric
 * Modif skxk : appel des nouvelles fonctions Valeur::mult_cp() et
 *  Valeur::mult_sp() pour la multiplication par cos(phi) et sin(phi).
 *
 * Revision 2.7  2000/02/09  19:32:11  eric
 * MODIF IMPORTANTE: la triade de decomposition est desormais passee en
 * argument des constructeurs.
 *
 * Revision 2.6  2000/02/01  14:14:25  eric
 * Ajout de la fonction amie flat_scalar_prod.
 *
 * Revision 2.5  2000/01/21  12:59:18  phil
 * ajout de skxk
 *
 * Revision 2.4  2000/01/14  09:29:43  eric
 * *** empty log message ***
 *
 * Revision 2.3  2000/01/13  17:22:37  phil
 * la fonction contraction de deux tenseurs ne passe plus par produit tensoriel
 *
 * Revision 2.2  2000/01/11  11:14:29  eric
 * Changement de nom pour la base vectorielle : base --> triad
 *
 * Revision 2.1  2000/01/10  17:25:15  eric
 * Gestion des bases vectorielles (triades de decomposition).
 *
 * Revision 2.0  1999/12/02  17:19:06  phil
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


Tenseur operator*(const Tenseur& t1, const Tenseur& t2) {
   
    assert ((t1.etat != ETATNONDEF) && (t2.etat != ETATNONDEF)) ;
    assert (t1.mp == t2.mp) ;
    
    int val_res = t1.valence + t2.valence ;
    double poids_res = t1.poids + t2.poids ;
    poids_res = (fabs(poids_res) < 1.e-10 ? 0. : poids_res) ;
    const Metrique* met_res = 0x0 ;
    if (poids_res != 0.) {
      //      assert((t1.metric != 0x0) || (t2.metric != 0x0)) ;
      if (t1.metric != 0x0) met_res = t1.metric ;
      else met_res = t2.metric ;
    }
    
   // cas scalaire :
    if (val_res == 0) {
	Tenseur scal(*t1.mp, met_res, poids_res) ;
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

	Tenseur res(*t1.mp, val_res, tipe, triad_res, met_res, poids_res) ;
	
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
    

Tenseur operator%(const Tenseur& t1, const Tenseur& t2) {
   
    assert ((t1.etat != ETATNONDEF) && (t2.etat != ETATNONDEF)) ;
    assert (t1.mp == t2.mp) ;
    
    int val_res = t1.valence + t2.valence ;
    double poids_res = t1.poids + t2.poids ;
    poids_res = (fabs(poids_res) < 1.e-10 ? 0. : poids_res) ;
    const Metrique* met_res = 0x0 ;
    if (poids_res != 0.) {
      // assert((t1.metric != 0x0) || (t2.metric != 0x0)) ;
      if (t1.metric != 0x0) met_res = t1.metric ;
      else met_res = t2.metric ;
    }
    
   // cas scalaire :
    if (val_res == 0) {
	Tenseur scal(*t1.mp, met_res, poids_res) ;
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

	Tenseur res(*t1.mp, val_res, tipe, triad_res, met_res, poids_res) ;


	
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



Tenseur contract(const Tenseur& source, int ind_1, int ind_2)  {
    
    
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
	
    Tenseur res(*source.mp, val_res, tipe, source.triad, source.metric, 
		source.poids) ;

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
	for (int j=0 ; j<3 ; j++) {
	    jeux_indice_source.set(ind_1) = j ;
	    jeux_indice_source.set(ind_2) = j ;
	    work = work + source(jeux_indice_source) ;
	    }
	    
	res.set(jeux_indice_res) = work ;
	}
    return res ;
}


Tenseur contract (const Tenseur& t1, int ind1, const Tenseur& t2, int ind2) {
    
    assert ((t1.etat != ETATNONDEF) && (t2.etat != ETATNONDEF)) ;
    // Verifs :
    assert ((ind1>=0) && (ind1<t1.valence)) ;
    assert ((ind2>=0) && (ind2<t2.valence)) ;
    assert (*(t1.mp) == *(t2.mp)) ;
    
    // Contraction possible ?
    if ( (t1.valence != 0) && (t2.valence != 0) ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    assert (t1.type_indice(ind1) != t2.type_indice(ind2)) ;
    
    int val_res = t1.valence + t2.valence - 2;
    double poids_res = t1.poids + t2.poids ;
    poids_res = (fabs(poids_res) < 1.e-10 ? 0. : poids_res) ;
    const Metrique* met_res = 0x0 ;
    if (poids_res != 0.) {
      //  assert((t1.metric != 0x0) || (t2.metric != 0x0)) ;
      if (t1.metric != 0x0) met_res = t1.metric ;
      else met_res = t2.metric ;
    }
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

    Tenseur res(*t1.mp, val_res, tipe, triad_res, met_res, poids_res) ;

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
    
    for (int comp=0 ; comp<res.n_comp ; comp++) {
	Itbl jeux_indice_res (res.donne_indices(comp)) ;
	for (int i=0 ; i<ind1 ; i++)
	    jeux_indice_t1.set(i) = jeux_indice_res(i) ;
	for (int i=ind1+1 ; i<t1.valence ; i++)
	    jeux_indice_t1.set(i) = jeux_indice_res(i-1) ;
	for (int i=0 ; i<ind2 ; i++)
	    jeux_indice_t2.set(i) = jeux_indice_res(t1.valence+i-1) ;
	for (int i=ind2+1 ; i<t2.valence ; i++)
	    jeux_indice_t2.set(i) = jeux_indice_res(t1.valence+i-2) ;
	
	    
	    
	work.set_etat_zero() ;
	for (int j=0 ; j<3 ; j++) {
	    jeux_indice_t1.set(ind1) = j ;
	    jeux_indice_t2.set(ind2) = j ;
	    work = work + t1(jeux_indice_t1)*t2(jeux_indice_t2) ;
	    }
	    
	res.set(jeux_indice_res) = work ;
	}
    return res ;
}

Tenseur contract_desal (const Tenseur& t1, int ind1, const Tenseur& t2, int ind2) {
    
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
    double poids_res = t1.poids + t2.poids ;
    poids_res = (fabs(poids_res) < 1.e-10 ? 0. : poids_res) ;
    const Metrique* met_res = 0x0 ;
    if (poids_res != 0.) {
      //  assert((t1.metric != 0x0) || (t2.metric != 0x0)) ;
      if (t1.metric != 0x0) met_res = t1.metric ;
      else met_res = t2.metric ;
    }
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

    Tenseur res(*t1.mp, val_res, tipe, triad_res, met_res, poids_res) ;

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
    
    for (int comp=0 ; comp<res.n_comp ; comp++) {
	Itbl jeux_indice_res (res.donne_indices(comp)) ;
	for (int i=0 ; i<ind1 ; i++)
	    jeux_indice_t1.set(i) = jeux_indice_res(i) ;
	for (int i=ind1+1 ; i<t1.valence ; i++)
	    jeux_indice_t1.set(i) = jeux_indice_res(i-1) ;
	for (int i=0 ; i<ind2 ; i++)
	    jeux_indice_t2.set(i) = jeux_indice_res(t1.valence+i-1) ;
	for (int i=ind2+1 ; i<t2.valence ; i++)
	    jeux_indice_t2.set(i) = jeux_indice_res(t1.valence+i-2) ;
	
	    
	    
	work.set_etat_zero() ;
	for (int j=0 ; j<3 ; j++) {
	    jeux_indice_t1.set(ind1) = j ;
	    jeux_indice_t2.set(ind2) = j ;
	    work = work + t1(jeux_indice_t1)%t2(jeux_indice_t2) ;
	    }
	    
	res.set(jeux_indice_res) = work ;
	}
    return res ;
}


Tenseur manipule(const Tenseur& t1, const Metrique& met, int place) {
    
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
    
    Tenseur auxi(*t1.mp, valen, tipe, t1.triad) ;
    
    if (t1.type_indice(place) == COV)
	auxi = contract (met.con(), 1, t1, place) ;
    else
	auxi = contract (met.cov(), 1, t1, place) ;
    
    // On doit remettre les indices a la bonne place ...
    
    for (int i=0 ; i<valen ; i++)
	tipe.set(i) = t1.type_indice(i) ;
    tipe.set(place) *= -1 ;
    
    Tenseur res(*t1.mp, valen, tipe, t1.triad, auxi.metric, auxi.poids) ;
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

Tenseur manipule (const Tenseur& t1, const Metrique& met) {
    
    Tenseur* auxi ;
    Tenseur* auxi_old = new Tenseur(t1) ;
    
    for (int i=0 ; i<t1.valence ; i++) {
	auxi = new Tenseur(manipule(*auxi_old, met, i)) ;
	delete auxi_old ;
	auxi_old = new Tenseur(*auxi) ;
	delete auxi ;
    }
    
    Tenseur result(*auxi_old) ;
    delete auxi_old ;
    return result ;
}


Tenseur skxk(const Tenseur& source) {
    
    // Verification
    assert (source.valence > 0) ;
    assert (source.etat != ETATNONDEF) ;
    assert (*source.triad == source.mp->get_bvect_cart()) ;
    
    // Le resultat :
    int val_res = source.valence-1 ;
    Itbl tipe (val_res) ;
    tipe.set_etat_qcq() ;
    for (int i=0 ; i<val_res ; i++)
	tipe.set(i) = source.type_indice(i) ;
    
    
    Tenseur res (*source.mp, val_res, tipe, source.triad, source.metric,
		 source.poids) ;
    
    if (source.etat == ETATZERO)
	res.set_etat_zero() ;
    else {
	res.set_etat_qcq() ;
	
	for (int i=0 ; i<res.n_comp ; i++) {
	    Itbl indices_res (res.donne_indices(i)) ;
	    Itbl indices_so (val_res+1) ;
	    indices_so.set_etat_qcq() ;
	    for (int j=0 ; j<val_res ; j++)
		indices_so.set(j) = indices_res(j) ;
	    // x S_x
	    // -----
	    indices_so.set(val_res) = 0 ;
	    Cmp resu(source(indices_so)) ;
        
	    resu.mult_r() ;			    // Multipl. by r

	    // What follows is valid only for a mapping of class Map_radial : 
	    assert( dynamic_cast<const Map_radial*>(source.get_mp()) != 0x0) ; 

	    resu.va = (resu.va).mult_st() ;	    // Multipl. by sin(theta)
	    resu.va = (resu.va).mult_cp() ;	    // Multipl. by cos(phi)

	// y S_y
	// -----
	    indices_so.set(val_res) = 1 ;
	    Cmp auxiliaire (source(indices_so)) ;
        
	    auxiliaire.mult_r() ;			    // Multipl. by r

	    auxiliaire.va = (auxiliaire.va).mult_st() ;  // Multipl. by sin(theta)
	    auxiliaire.va = (auxiliaire.va).mult_sp() ;  // Multipl. by sin(phi)
    
	    resu = resu + auxiliaire ; 
    
	    // z S_z
	    // -----
	    indices_so.set(val_res) = 2 ;
	    auxiliaire = source(indices_so) ;
        
	    auxiliaire.mult_r() ;			    // Multipl. by r

	    auxiliaire.va = (auxiliaire.va).mult_ct() ;     // Multipl. by cos(theta)
   
	    resu = resu + auxiliaire ; 
    
	    res.set(indices_res) = resu ;
	    // The End 
	    // -------
	    }
    }
    return res ;
}

Tenseur flat_scalar_prod(const Tenseur& t1, const Tenseur& t2) {
    
    assert ((t1.etat != ETATNONDEF) && (t2.etat != ETATNONDEF)) ;
    // Verifs :
    assert (t1.mp == t2.mp) ;
    
    // Contraction possible ?
    if ( (t1.valence != 0) && (t2.valence != 0) ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    
    int val_res = t1.valence + t2.valence - 2;
    double poids_res = t1.poids + t2.poids ;
    poids_res = (fabs(poids_res) < 1.e-10 ? 0. : poids_res) ;
    const Metrique* met_res = 0x0 ;
    if (poids_res != 0.) {
      assert((t1.metric != 0x0) || (t2.metric != 0x0)) ;
      if (t1.metric != 0x0) met_res = t1.metric ;
      else met_res = t2.metric ;
    }
    Itbl tipe(val_res) ;
    tipe.set_etat_qcq() ;

    for (int i=0 ; i<t1.valence - 1 ; i++)
	tipe.set(i) = t1.type_indice(i) ;
    for (int i = t1.valence-1 ; i<val_res ; i++)
	tipe.set(i) = t2.type_indice(i-t1.valence+2) ;
	
    Tenseur res(*t1.mp, val_res, tipe, t1.triad, met_res, poids_res) ;

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
	
	for (int j=0 ; j<3 ; j++) {
	    jeux_indice_t1.set(t1.valence - 1) = j ;
	    jeux_indice_t2.set(0) = j ;
	    work = work + t1(jeux_indice_t1)*t2(jeux_indice_t2) ;
	    }
	    
	res.set(jeux_indice_res) = work ;
    
    }	// fin de la boucle sur les composantes du resultat
    
    return res ;
}



Tenseur flat_scalar_prod_desal(const Tenseur& t1, const Tenseur& t2) {
    
    assert ((t1.etat != ETATNONDEF) && (t2.etat != ETATNONDEF)) ;
    // Verifs :
    assert (t1.mp == t2.mp) ;
    
    // Contraction possible ?
    if ( (t1.valence != 0) && (t2.valence != 0) ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    
    int val_res = t1.valence + t2.valence - 2;
    double poids_res = t1.poids + t2.poids ;
    poids_res = (fabs(poids_res) < 1.e-10 ? 0. : poids_res) ;
    const Metrique* met_res = 0x0 ;
    if (poids_res != 0.) {
      assert((t1.metric != 0x0) || (t2.metric != 0x0)) ;
      if (t1.metric != 0x0) met_res = t1.metric ;
      else met_res = t2.metric ;
    }
    Itbl tipe(val_res) ;
    tipe.set_etat_qcq() ;

    for (int i=0 ; i<t1.valence - 1 ; i++)
	tipe.set(i) = t1.type_indice(i) ;
    for (int i = t1.valence-1 ; i<val_res ; i++)
	tipe.set(i) = t2.type_indice(i-t1.valence+2) ;
	
    Tenseur res(*t1.mp, val_res, tipe, t1.triad, met_res, poids_res) ;

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
	
	for (int j=0 ; j<3 ; j++) {
	    jeux_indice_t1.set(t1.valence - 1) = j ;
	    jeux_indice_t2.set(0) = j ;
	    work = work + t1(jeux_indice_t1) % t2(jeux_indice_t2) ;
	    }
	    
	res.set(jeux_indice_res) = work ;
    
    }	// fin de la boucle sur les composantes du resultat
    
    return res ;
}


Tenseur lie_derive (const Tenseur& t, const Tenseur& x, const Metrique* met)
{
  assert ( (t.get_etat() != ETATNONDEF) && (x.get_etat() != ETATNONDEF) ) ;
  assert(x.get_valence() == 1) ;
  assert(x.get_type_indice(0) == CON) ;
  assert(x.get_poids() == 0.) ;
  assert(t.get_mp() == x.get_mp()) ;
 
  int val = t.get_valence() ;
  double poids = t.get_poids() ;
  Itbl tipe (val+1) ;
  tipe.set_etat_qcq() ;
  tipe.set(0) = COV ;
  Itbl tipx(2) ;
  tipx.set_etat_qcq() ;
  tipx.set(0) = COV ;
  tipx.set(1) = CON ;
  for (int i=0 ; i<val ; i++)
    tipe.set(i+1) = t.get_type_indice(i) ;
  Tenseur dt(*t.get_mp(), val+1, tipe, t.get_triad(), t.get_metric(), poids) ;
  Tenseur dx(*x.get_mp(), 2, tipx, x.get_triad()) ; 
  if (met == 0x0) {
    dx = x.gradient() ;
    dt = t.gradient() ;
  }
  else {
    dx = x.derive_cov(*met) ;
    dt = t.derive_cov(*met) ;
  }
  Tenseur resu(contract(x,0,dt,0)) ;
  Tenseur* auxi ;
  if ((val!=0)&&(t.get_etat()!=ETATZERO)&&(x.get_etat()!=ETATZERO)) {
    assert(t.get_triad()->identify() == x.get_triad()->identify()) ;

    for (int i=0 ; i<val ; i++) {
      if (t.get_type_indice(i) == COV) {
	auxi = new Tenseur(contract(t,i,dx,1)) ;
	
	Itbl indices_aux(val) ;
	indices_aux.set_etat_qcq() ;
	for (int j=0 ; j<resu.get_n_comp() ; j++) {
	  
	  Itbl indices (resu.donne_indices(j)) ;
	  indices_aux.set(val-1) = indices(i) ;
	  for (int idx=0 ; idx<val-1 ; idx++)
	    if (idx<i)
	      indices_aux.set(idx) = indices(idx) ;
	    else
	      indices_aux.set(idx) = indices(idx+1) ;
	  
	  resu.set(indices) += (*auxi)(indices_aux) ;
	}
      }   
      else {
	auxi = new Tenseur(contract(t,i,dx,0)) ;
	
	Itbl indices_aux(val) ;
	indices_aux.set_etat_qcq() ;
	
	//On range comme il faut :
	for (int j=0 ; j<resu.get_n_comp() ; j++) {
	  
	  Itbl indices (resu.donne_indices(j)) ;
	  indices_aux.set(val-1) = indices(i) ;
	  for (int idx=0 ; idx<val-1 ; idx++)
	    if (idx<i)
	      indices_aux.set(idx) = indices(idx) ;
	    else
	      indices_aux.set(idx) = indices(idx+1) ;
	  resu.set(indices) -= (*auxi)(indices_aux) ;
	}
      }
      delete auxi ;
    }
  }
  if ((poids != 0.)&&(t.get_etat()!=ETATZERO)&&(x.get_etat()!=ETATZERO)) 
    resu = resu + poids*contract(dx,0,1)*t ;

  return resu ;
}

Tenseur sans_trace(const Tenseur& t, const Metrique& metre) 
{
  assert(t.get_etat() != ETATNONDEF) ;
  assert(metre.get_etat() != ETATNONDEF) ;
  assert(t.get_valence() == 2) ;

  Tenseur resu(t) ;
  if (resu.get_etat() == ETATZERO) return resu ;
  assert(resu.get_etat() == ETATQCQ) ;

  int t0 = t.get_type_indice(0) ;
  int t1 = t.get_type_indice(1) ;
  Itbl mix(2) ;
  mix.set_etat_qcq() ;
  mix.set(0) = (t0 == t1 ? -t0 : t0) ;
  mix.set(1) = t1 ;

  Tenseur tmp(*t.get_mp(), 2, mix, *t.get_triad(), t.get_metric(), 
	      t.get_poids()) ;
  if (t0 == t1)
    tmp = manipule(t, metre, 0) ;
  else
    tmp = t ;

  Tenseur trace(contract(tmp, 0, 1)) ;

  if (t0 == t1) {
	switch (t0) {
	case COV : {
	  resu = resu - 1./3.*trace * metre.cov() ;
	  break ;
	}
	case CON : {
	  resu = resu - 1./3.*trace * metre.con() ;	
	  break ;
	}
	default :
	  cout << "Erreur bizarre dans sans_trace!" << endl ;
	  abort() ;
	  break ;
	}
  }
  else {
    Tenseur_sym delta(*t.get_mp(), 2, mix, *t.get_triad(), 
		      t.get_metric(), t.get_poids()) ;
    delta.set_etat_qcq() ;
    for (int i=0; i<3; i++) 
      for (int j=i; j<3; j++)
	delta.set(i,j) = (i==j ? 1 : 0) ;
    resu = resu - trace/3. * delta ;
  }
  resu.set_std_base() ;
  return resu ;
}
    

  
