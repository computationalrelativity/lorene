/*
 *  Methods of class Mtbl_cf
 *
 *   (see file mtbl_cf.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char mtbl_cf_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/05/07 07:36:03  e_gourgoulhon
 * Compatibilty with xlC compiler on IBM SP2:
 *    suppressed the parentheses around argument of instruction new:
 * 	e.g.   t = new (Tbl *[nzone])  -->   t = new Tbl*[nzone]
 *
 * Revision 1.2  2001/12/04 21:27:54  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.13  2000/08/16  10:43:09  eric
 * Suppression du membre dzpuis.
 *
 * Revision 2.12  2000/02/25  10:27:44  eric
 * Suppression des appels a nettoie() dans l'affichage.
 *
 * Revision 2.11  1999/10/29  15:07:27  eric
 * Suppression des fonctions membres min() et max():
 * elles deviennent des fonctions externes.
 *
 * Revision 2.10  1999/10/21  13:42:05  eric
 * *** empty log message ***
 *
 * Revision 2.9  1999/10/18  15:16:12  eric
 * *** empty log message ***
 *
 * Revision 2.8  1999/10/18  15:08:44  eric
 * La fonction membre annule() est rebaptisee annule_hard().
 * Introduction de la fonction membre annule(int, int).
 *
 * Revision 2.7  1999/10/13  15:52:03  eric
 * Depoussierage.
 * Ajout du membre base.
 *
 * Revision 2.6  1999/10/01  14:49:38  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.5  1999/03/03  10:35:37  hyc
 * *** empty log message ***
 *
 * Revision 2.4  1999/03/02  16:26:39  eric
 * Modif des indentations dans <<
 *
 * Revision 2.3  1999/03/02  15:34:30  eric
 * Anglicisation des messages...
 *
 *
 * $Header$
 *
 */

// headers C++
#include <iostream.h>

// headers C
#include <assert.h>
#include <stdlib.h>
#include <math.h>

// headers Lorene
#include "mtbl_cf.h"
#include "utilitaires.h"

// Prototypage
void c_est_pas_fait(char * ) ;

// Constructeur
// ------------

Mtbl_cf::Mtbl_cf(const Mg3d& g, const Base_val& ba) : mg(&g), 
						      etat(ETATNONDEF),
						      base(ba),  
						      t(0x0) {
    nzone = g.get_nzone() ;
    assert(base.nzone == nzone) ; 
}

Mtbl_cf::Mtbl_cf(const Mg3d* g, const Base_val& ba) : mg(g), 
						      etat(ETATNONDEF), 
						      base(ba),  
						      t(0x0) {
    nzone = g->get_nzone() ;
    assert(base.nzone == nzone) ; 
}


// Destructeur
// -----------
Mtbl_cf::~Mtbl_cf() {

    del_t() ;
}

// Copie
// -----
Mtbl_cf::Mtbl_cf(const Mtbl_cf& mtc) : mg(mtc.mg), 
				       nzone(mtc.nzone),
				       base(mtc.base) {

    // Protection
    assert(mtc.get_etat() != ETATNONDEF) ;
    
    t = 0x0 ;
    etat = ETATNONDEF ;
    if (mtc.etat == ETATQCQ) {
	set_etat_qcq() ;
	for (int i=0 ; i<nzone ; i++) {
	    *t[i] = *mtc.t[i] ;
	}
    }
    else {
	assert(mtc.etat == ETATZERO) ;	// sinon...
    }
    etat = mtc.etat ;
}

// Constructeur a partir d'une grille et d'un fichier
// --------------------------------------------------
Mtbl_cf::Mtbl_cf(const Mg3d& g, FILE* fd) : mg(&g), 
					    base(fd) {
    
    // La multi-grille
    Mg3d* mg_tmp = new Mg3d(fd) ;	// la multi-grille d'origine
    if (*mg != *mg_tmp) {
	cout << "Mtbl_cf::Mtbl_cf(const Mg3d& , FILE*): grid not consistent !" 
	     << endl ;
	abort() ;
    }
    delete mg_tmp ;
    
    // Lecture
    nzone = mg->get_nzone() ;
    assert(base.nzone == nzone) ; 
    fread_be(&etat, sizeof(int), 1, fd) ;		// etat
    
    // Le tableau
    t = 0x0 ;
    if (etat == ETATQCQ) {
	t = new Tbl*[nzone] ;
	for (int i=0 ; i<nzone ; i++) {
	    t[i] = new Tbl(fd) ;
	}
    }
    int dzpuis_vieux ;
    fread_be(&dzpuis_vieux, sizeof(int), 1, fd) ;	    // le vieux dzpuis !
}

// Sauvegarde sur un fichier
// -------------------------
void Mtbl_cf::sauve(FILE* fd) const {

    base.sauve(fd) ;			    // la base
    mg->sauve(fd) ;			    // la multi-grille
    fwrite_be(&etat, sizeof(int), 1, fd) ;		    // etat
    if (etat == ETATQCQ) {
	for (int i=0 ; i<nzone ; i++) {
	    t[i]->sauve(fd) ;
	}
    }
    int dzpuis_vieux = 0 ;
    fwrite_be(&dzpuis_vieux, sizeof(int), 1, fd) ;	    // le vieux dzpuis !
}

// Affectations
// ------------
void Mtbl_cf::operator=(const Mtbl_cf& mtc)
{
    // Protection
    assert (mg == mtc.mg) ;
    assert(mtc.get_etat() != ETATNONDEF) ;

    // Les choses fixes
    base = mtc.base ; 
    
    // Gestion des donnees
    if (mtc.get_etat() == ETATZERO) {
	set_etat_zero() ;
    }
    else {
    	assert(mtc.get_etat() == ETATQCQ) ; // sinon...
	set_etat_qcq() ;
	for (int i=0 ; i<nzone ; i++) {
	    *t[i] = *mtc.t[i] ;
	}
    }
}

void Mtbl_cf::operator=(double x)
{
    if ( x == double(0) ) {
	set_etat_zero() ;
    }
    else {
	set_etat_qcq() ;
	for (int i=0 ; i<nzone ; i++) {
	    *t[i] = x ;
	}
    }

}

void Mtbl_cf::operator=(int m)
{
    if (m == 0) {
	set_etat_zero() ;
    }
    else {
	set_etat_qcq() ;
	for (int i=0 ; i<nzone ; i++) {
	    *t[i] = m ;
	}
    }

}


		    //-----------------//
    	    	    // Gestion memoire //
		    //-----------------//

// Destructeur logique
void Mtbl_cf::del_t() {
    // Detruit le mtbl_cf
    if (t != 0x0) {
    for (int l=0 ; l<nzone ; l++) {
	delete t[l] ;
    }
    delete [] t ;
    t = 0x0 ;
    }
    etat = ETATNONDEF ;
}
// ETATZERO
void Mtbl_cf::set_etat_zero() {
    if (etat == ETATZERO) return ;
    del_t() ;
    etat = ETATZERO ;
}
// ETATNONDEF
void Mtbl_cf::set_etat_nondef() {
    if (etat == ETATNONDEF) return ;
    del_t() ;
    etat = ETATNONDEF ;
}
// ETATQCQ
void Mtbl_cf::set_etat_qcq() {
    if (etat == ETATQCQ) return ;
    t = new Tbl*[nzone] ;
    for (int i=0 ; i<nzone ; i++) {
	int nbr = mg->get_nr(i) ;
	int nbt = mg->get_nt(i) ;
	int nbp = mg->get_np(i) ;
	t[i] = new Tbl(nbp+2, nbt, nbr) ;
    }
    etat = ETATQCQ ;
}
// ZERO hard
void Mtbl_cf::annule_hard() {
    if (t == 0x0) {
	t = new Tbl*[nzone] ;
	for (int i=0 ; i<nzone ; i++) {
	    int nbr = mg->get_nr(i) ;
	    int nbt = mg->get_nt(i) ;
	    int nbp = mg->get_np(i) ;
	    t[i] = new Tbl(nbp+2, nbt, nbr) ;
	}
    }
    
    for (int i=0 ; i<nzone ; i++) {
	t[i]->annule_hard() ;
    }
    etat = ETATQCQ ;
}

// Sets the {\tt Mtbl_cf} to zero in some domains
// ----------------------------------------------

void Mtbl_cf::annule(int l_min, int l_max) {

    assert( (l_min >= 0) && (l_min < nzone) ) ; 
    assert( (l_max >= 0) && (l_max < nzone) ) ; 
    
    // Cas particulier: annulation globale : 
    if ( (l_min == 0) && (l_max == nzone-1) ) {
	set_etat_zero() ;
	return ; 
    }
    
    assert( etat != ETATNONDEF ) ; 
    
    if ( etat == ETATZERO ) {
	return ;		// rien n'a faire si c'est deja zero
    }
    else {
	assert( etat == ETATQCQ ) ;	// sinon...
	for (int l=l_min; l<=l_max; l++) {
	    t[l]->set_etat_zero() ; 
	}
	 
    }
    
}

			//------------------------//
			//	Display		  //
			//------------------------//
			
//-----------			
// Operator<<
//-----------			
			
ostream& operator<<(ostream& o, const Mtbl_cf& mt) { 
    // Protection
    assert(mt.get_etat() != ETATNONDEF) ;
    
    int nzone = mt.get_nzone() ;
    o.precision(4);
    o.setf(ios::showpoint);
    o << "*** MTBL_CF " << nzone << " domains" << incindent << iendl ;

    o << mt.base << iendl ; 

    o << "Values of the coefficients : " << iendl ; 
    if (mt.get_etat() == ETATZERO) {
	o << "Logically NUL" << iendl ;
    }
    else {
	for (int l=0 ; l<nzone ; l++) {
	    o << " Domain #" << l << incindent << iendl ;
	    o << *(mt.t[l]) ;
	    o << decindent << iendl ;
	}
    }

    o << decindent << iendl ;
    return o ;
}

//---------------
// Affiche_seuil
//---------------

void Mtbl_cf::affiche_seuil(ostream& ost, int precis,  double seuil) const {
    ost << "*** Mtbl_cf " << nzone << " domains" << endl ;
    ost << base << endl ; 
    ost << "Values of the coefficients : " << endl ; 
   
    // Cas particuliers
    //-----------------

    if (etat == ETATNONDEF) {
	ost << "    state: UNDEFINED" << endl ;
	return ;
    }

    if (etat == ETATZERO) {
	ost << "    state: ZERO" << endl ;
	return ;
    }

    // Affichage des Tbl
    //------------------
    assert( t != 0x0 ) ; 
        
    for (int l=0; l < nzone; l++) {
	t[l]->affiche_seuil( ost , precis, seuil ) ;
    }

}


// To be done
//-----------

void Mtbl_cf::operator*=(double ) {
    char* f = __FILE__ ;
    c_est_pas_fait(f) ;
}

void Mtbl_cf::operator/=(double ) {
    char* f = __FILE__ ;
    c_est_pas_fait(f) ;
}


