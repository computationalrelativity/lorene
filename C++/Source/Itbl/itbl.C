/*
 * Methods of class Itbl
 *
 *  (see file itbl.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char itbl_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/10/16 14:36:37  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2001/12/04 21:27:53  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  1999/11/23  13:17:09  eric
 * Le constructeur Itbl::Itbl(const Dim_tbl ) devient desormais
 *   tbl::Itbl(const Dim_tbl& ).
 * La taille zero est autorisee par le constructeur 1D.
 * Modif affichage.
 *
 * Revision 2.0  1999/11/17  16:04:38  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// headers C
#include <math.h>

// headers Lorene
#include "itbl.h"
#include "indent.h"
#include "utilitaires.h"


			//---------------//
			// Constructeurs //
			//---------------//


// Constructeur 1D
Itbl::Itbl(int n1) : etat(ETATNONDEF), dim(n1), t(0x0) {
    if (n1 == 0) {
	set_etat_zero() ; 
    }
}

// Constructeur 2D
Itbl::Itbl(int n1, int n0) : etat(ETATNONDEF), dim(n1, n0), t(0x0) {}

// Constructeur 3D
Itbl::Itbl(int n2, int n1, int n0) : etat(ETATNONDEF), dim(n2, n1, n0), t(0x0) {}

// Constructeur a partir d'un Dim_tbl
Itbl::Itbl(const Dim_tbl& dt) : etat(ETATNONDEF), dim(dt), t(0x0) {
    if (get_taille() == 0) {
	set_etat_zero() ; 
    }
}

// Copie
Itbl::Itbl(const Itbl& tc) : etat(tc.etat), dim(tc.dim) {

    // La valeur eventuelle
    if (tc.etat == ETATQCQ) {
	t = new int[get_taille()] ;
	for (int i=0 ; i<get_taille() ; i++) {
	    t[i] = tc.t[i] ;
	}
    }
    else{
	t = 0x0 ; 
    }
    
}

// From file
Itbl::Itbl(FILE* fd) : dim(fd) {

    fread_be(&etat, sizeof(int), 1, fd) ;		// etat
    
    // Le tableau
    if (etat == ETATQCQ) {
	t = new int[get_taille()] ;
	fread_be(t, sizeof(int), get_taille(), fd) ;	    // le tableau
    }
    else{
	t = 0x0 ; 
    }
}

			//-------------//
			// Destructeur //
			//-------------//

Itbl::~Itbl() {
    delete [] t ;
}

			//-------------//
			// Affectation //
			//-------------//

// From Itbl
void Itbl::operator=(const Itbl& tx)
{
    // Protection
    assert( dim == tx.dim ) ;
    assert(tx.get_etat() != ETATNONDEF) ;

    int n = get_taille() ;
    switch (tx.etat) {
	case ETATZERO:
	set_etat_zero() ;
	break ;
	
	case ETATQCQ:
	set_etat_qcq() ;
	for (int i=0 ; i<n ; i++) {
	    t[i] = tx.t[i] ;
	}
	break ;
	
	default:
	cout << "Erreur bizarre !" << endl ;
	abort() ;
	break ;
    }
}

// From int
void Itbl::operator=(int a)
{
    if ( a == 0 ) {
	set_etat_zero() ;
    }
    else {
	int n = get_taille() ;
	if (n > 0) {
	    set_etat_qcq() ;
	    for (int i=0 ; i<n ; i++) {
		t[i] = a ;
	    }
	}
    }
}

   
			//------------//
			// Sauvegarde //
			//------------//

// save in a file

void Itbl::sauve(FILE* fd) const {

    dim.sauve(fd) ;	    	    	    	    // dim
    fwrite_be(&etat, sizeof(int), 1, fd) ;		    // etat
    if (etat == ETATQCQ) {
	fwrite_be(t, sizeof(int), get_taille(), fd) ;	    // le tableau
    }
}
    
		    //-----------------//
    	    	    // Gestion memoire //
		    //-----------------//

// Destructeur logique
void Itbl::del_t() {
    delete [] t ;
    t = 0x0 ;
    etat = ETATNONDEF ;
}

// ETATZERO
void Itbl::set_etat_zero() {
    if (etat == ETATZERO) return ;
    del_t() ;
    etat = ETATZERO ;
}

// ETATNONDEF
void Itbl::set_etat_nondef() {
    if (etat == ETATNONDEF) return ;
    del_t() ;
    etat = ETATNONDEF ;
}

// ETATQCQ
void Itbl::set_etat_qcq() {
    if (etat == ETATQCQ) return ;

    // Protection
    assert( (etat == ETATZERO) || (etat == ETATNONDEF) ) ; // sinon...

    t = new int[get_taille()] ;
    etat = ETATQCQ ;
}

// ZERO hard
void Itbl::annule_hard() {
    if (t == 0x0) {
	t = new int[get_taille()] ;
    }
    for (int i=0 ; i<get_taille() ; i++) {
	t[i] = 0 ;
    }
    etat = ETATQCQ ;
}


			//------------------------//
			//	Display		  //
			//------------------------//
			
//-----------			
// Operator<<
//-----------			

ostream& operator<<(ostream& o, const Itbl& t) {
    
    int ndim = t.get_ndim() ;
    o.precision(4);
    o.setf(ios::showpoint);
    o << "*** Itbl " << ndim << "D" << "   size: " ; 
    for (int i = 0; i<ndim-1; i++) {
	o << t.get_dim(i) << " x " ;
    } 
    o << t.get_dim(ndim-1) << incindent << iendl ;

    if (t.get_etat() == ETATZERO) {
	o << "Identically ZERO" << decindent << iendl ;
	return o ;
    }

    if (t.get_etat() == ETATNONDEF) {
	o << "UNDEFINED STATE" << decindent << iendl ;
	return o ;
    }

    assert(t.etat == ETATQCQ) ;
    switch (ndim) {

	case 1 : {
	    for (int i=0 ; i<t.get_dim(0) ; i++) {
		o << " " << t(i)  ;
	    }
	    o << decindent << endl ;
	    break ;
	}


	case 2 : {
	    for (int j=0 ; j<t.get_dim(1) ; j++) {
		o << " J = " << j << " : " << incindent << iendl ;
		for (int i=0 ; i<t.get_dim(0) ; i++) {
		    o << " " << t(j, i)  ;
		}
		o << decindent << iendl ;
	    }
	    o << decindent << endl ;
	    break ;
	}
		
	case 3 : {
	    for (int k=0 ; k<t.get_dim(2) ; k++) {
		o << " K = " << k << " : " << incindent << iendl ;
		for (int j=0 ; j<t.get_dim(1) ; j++) {
		    o << " J = " << j << " : " << incindent ;
		    for (int i=0 ; i<t.get_dim(0) ; i++) {
			o << " " << t(k, j, i)  ;
		    }
		    o << decindent << iendl ;
		}
		o << decindent << iendl ;
	    }
	    o << decindent << endl ;
	    break ;
	}
		
	default : {
	    cout << "operator<< Itbl : unexpected dimension !" << endl ;
	    cout << " ndim = " << ndim << endl ; 	
	    abort() ;
	    break ;
	}
    }
    return o ;
}
