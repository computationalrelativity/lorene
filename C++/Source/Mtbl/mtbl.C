/*
 *  Methods of class Mtbl
 *
 *   (see file mtbl.h for documentation)
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


char mtbl_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.10  2000/08/16  10:30:04  eric
 * Suppression du membre dzpuis.
 *
 * Revision 2.9  1999/11/23  13:33:04  eric
 *  Le constructeur Tbl::Tbl(const Grille3d* ) est devenu Tbl::Tbl(const Grille3d& ).
 *
 * Revision 2.8  1999/10/29  15:06:24  eric
 * Suppression des fonctions membres min() et max():
 * elles deviennent des fonctions externes.
 *
 * Revision 2.7  1999/10/18  15:16:05  eric
 * *** empty log message ***
 *
 * Revision 2.6  1999/10/18  15:08:22  eric
 * La fonction membre annule() est rebaptisee annule_hard().
 * Introduction de la fonction membre annule(int, int).
 *
 * Revision 2.5  1999/10/01  12:35:54  eric
 * Ajout de la fonction affiche_seuil.
 *
 * Revision 2.4  1999/10/01  10:08:41  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.3  1999/03/02  16:26:26  eric
 * Modif des indentations dans <<
 *
 * Revision 2.2  1999/03/02  15:34:08  eric
 * Anglicisation des commentaires...
 *
 * Revision 2.1  1999/02/22  15:24:13  hyc
 * *** empty log message ***
 *
 *
 * Revision 2.0  1999/01/15  09:10:39  hyc
 * *** empty log message ***
 *
 * $Header$
 *
 */

// headers C++
#include <iostream.h>

// headers C
#include <assert.h>

// headers Lorene
#include "mtbl.h"
#include "coord.h"
#include "type_parite.h"

// Prototypage
void c_est_pas_fait(char * ) ;

// Constructeurs
// -------------
Mtbl::Mtbl(const Mg3d& g) : mg(&g), etat(ETATNONDEF), t(0x0) {

    nzone = g.get_nzone() ;
    
}

Mtbl::Mtbl(const Mg3d* g) : mg(g), etat(ETATNONDEF), t(0x0) {

    nzone = g->get_nzone() ;

}

Mtbl::Mtbl(const Coord& c1) {
    
    // La coordonnee est-elle a jour ?
    if (c1.c == 0x0) c1.fait() ;

    // Les donnees fixes
    mg = c1.c->get_mg() ;
    nzone = mg->get_nzone() ;

    // L'etat
    t = 0x0 ;
    etat = ETATNONDEF ;
        
    // La transformation
    *this = *(c1.c) ;
}


// Destructeur
// -----------
Mtbl::~Mtbl() {
    del_t() ;
}

// Copie
// -----
Mtbl::Mtbl(const Mtbl& mtc) : mg(mtc.mg), nzone(mtc.nzone) {

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
Mtbl::Mtbl(const Mg3d & g, FILE* fd) : mg(&g) {
    
    // La multi-grille
    Mg3d* mg_tmp = new Mg3d(fd) ;	// la multi-grille d'origine
    if (*mg != *mg_tmp) {
	cout << "Mtbl::Mtbl(Mg3d & , FILE*): grid not consistent !" << endl ;
	abort() ;
    }
    delete mg_tmp ;
    
    // Lecture
    nzone = mg->get_nzone() ;
    fread(&etat, sizeof(int), 1, fd) ;		// etat
    
    // Le tableau
    t = 0x0 ;
    if (etat == ETATQCQ) {
	t = new (Tbl *[nzone]) ;
	for (int i=0 ; i<nzone ; i++) {
	    t[i] = new Tbl(fd) ;
	}
    }
    int dzpuis_vieux ; 
    fread(&dzpuis_vieux, sizeof(int), 1, fd) ;	    // le vieux dzpuis
}

// Sauvegarde sur un fichier
void Mtbl::sauve(FILE* fd) const {

    mg->sauve(fd) ;			    // la multi-grille
    fwrite(&etat, sizeof(int), 1, fd) ;		    // etat
    if (etat == ETATQCQ) {
	for (int i=0 ; i<nzone ; i++) {
	    t[i]->sauve(fd) ;
	}
    }
    int dzpuis_vieux = 0 ; 
    fwrite(&dzpuis_vieux, sizeof(int), 1, fd) ;	    // le vieux dzpuis
}

// Affectations
// ------------
void Mtbl::operator=(const Mtbl& mtc)
{
    // Protection
    assert (mg == mtc.mg) ;
    assert(mtc.get_etat() != ETATNONDEF) ;

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

void Mtbl::operator=(double x)
{
    if (x == double(0)) {
	set_etat_zero() ;
    }
    else {
	set_etat_qcq() ;
	for (int i=0 ; i<nzone ; i++) {
	    *t[i] = x ;
	}
    }

}

void Mtbl::operator=(int m)
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
void Mtbl::del_t() {
    if (t != 0x0) {
    for (int l=0 ; l<nzone ; l++) {
	delete t[l] ;
    }
    delete [] t ;
    t = 0x0 ;
    }
}
// ETATZERO
void Mtbl::set_etat_zero() {
    if (etat == ETATZERO) return ;
    del_t() ;
    etat = ETATZERO ;
}
// ETATNONDEF
void Mtbl::set_etat_nondef() {
    if (etat == ETATNONDEF) return ;
    del_t() ;
    etat = ETATNONDEF ;
}
// ETATQCQ
void Mtbl::set_etat_qcq() {
    if (etat == ETATQCQ) return ;

    // Protection
    assert( (etat == ETATZERO) || (etat == ETATNONDEF) ) ; // sinon...
    
    t = new (Tbl *[nzone]) ;
    for (int i=0 ; i<nzone ; i++) {
	t[i] = new Tbl( *(mg->get_grille3d(i)) ) ;
    }
    etat = ETATQCQ ;
}
// ZERO hard
void Mtbl::annule_hard() {
    if (t == 0x0) {
	t = new (Tbl *[nzone]) ;
	for (int i=0 ; i<nzone ; i++) {
	    t[i] = new Tbl( *(mg->get_grille3d(i)) ) ;
	}
    }
    
    for (int i=0 ; i<nzone ; i++) {
	t[i]->annule_hard() ;
    }
    etat = ETATQCQ ;
}

// Sets the {\tt Mtbl} to zero in some domains
// -------------------------------------------

void Mtbl::annule(int l_min, int l_max) {

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

ostream& operator<<(ostream& o, const Mtbl& mt) {
    // Protection
    assert(mt.get_etat() != ETATNONDEF) ;
    
    int nzone = mt.get_nzone() ;
    o.precision(4);
    o.setf(ios::showpoint);
    o << "*** Mtbl " << nzone << " domains" << incindent << iendl ;

    if (mt.get_etat() == ETATZERO) {
	o << "Logically NULL" << iendl ;
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

void Mtbl::affiche_seuil(ostream& ost, int precis,  double seuil) const {
    ost << "*** Mtbl " << nzone << " domains" << endl ;
    
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
    
    for (int l=0; l < nzone; l++) {
	t[l]->affiche_seuil( ost , precis, seuil ) ;
    }


}


// To be done
//-----------

void Mtbl::operator+=(double ) {
    char* f = __FILE__ ;
    c_est_pas_fait(f) ;
}

void Mtbl::operator-=(double ) {
    char* f = __FILE__ ;
    c_est_pas_fait(f) ;
}

void Mtbl::operator*=(double ) {
    char* f = __FILE__ ;
    c_est_pas_fait(f) ;
}

void Mtbl::operator/=(double ) {
    char* f = __FILE__ ;
    c_est_pas_fait(f) ;
}


