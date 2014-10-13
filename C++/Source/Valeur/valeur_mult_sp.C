/*
 * Methods for multiplication by sin(phi) for classes
 *   - Valeur
 *   - Mtbl_cf
 */

/*
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


char valeur_mult_sp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:53:50  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:23  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2012/01/17 15:08:32  j_penner
 * using MAX_BASE_2 for the phi coordinate
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2000/09/18  10:15:20  eric
 * Ajout des bases P_COSSIN_P et P_COSSIN_I
 *
 * Revision 2.1  2000/09/11  15:04:04  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/09/11  13:54:00  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cassert>

// Headers Lorene
#include "mtbl_cf.h"
#include "valeur.h"

// Local prototypes
namespace Lorene {
void _mult_sp_pas_prevu(Tbl*, int&) ;
void _mult_sp_p_cossin(Tbl*, int&) ;
void _mult_sp_p_cossin_p(Tbl*, int&) ;
void _mult_sp_p_cossin_i(Tbl*, int&) ;

// Version membre d'un Valeur
// --------------------------

const Valeur& Valeur::mult_sp() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_mult_sp != 0x0) {
	return *p_mult_sp ;
    }
    
    // ... si, il faut bosser

    p_mult_sp = new Valeur(mg) ;

    if (etat == ETATZERO) {
	p_mult_sp->set_etat_zero() ; 
    }
    else {
	assert(etat == ETATQCQ) ; 
	p_mult_sp->set_etat_cf_qcq() ;
	Mtbl_cf* cfp = p_mult_sp->c_cf ; // Pointeur sur le Mtbl_cf qui vient d'etre
					// cree par le set_etat_cf_qcq()

	// Initialisation de *cfp : recopie des coef. de la fonction
	if (c_cf == 0x0) {
	    coef() ;
	}
	*cfp = *c_cf ;	
 
	cfp->mult_sp() ;	// calcul 
    
	p_mult_sp->base = cfp->base ; // On remonte la base de sortie au niveau Valeur
    }

    // Termine
    return *p_mult_sp ;
}



// Version membre d'un Mtbl_cf
// ---------------------------

void Mtbl_cf::mult_sp() {

// Routines de derivation
static void (*_mult_sp[MAX_BASE_2])(Tbl *, int &) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE_2 ; i++) {
	    _mult_sp[i] = _mult_sp_pas_prevu ;
	}
	// Les routines existantes
	_mult_sp[P_COSSIN >> TRA_P] = _mult_sp_p_cossin ;
	_mult_sp[P_COSSIN_P >> TRA_P] = _mult_sp_p_cossin_p ;
	_mult_sp[P_COSSIN_I >> TRA_P] = _mult_sp_p_cossin_i ;
    }

    // Debut de la routine 

    // Protection
    assert(etat == ETATQCQ) ;
    
    // Boucle sur les zones
    for (int l=0 ; l<nzone ; l++) {
	int base_p = (base.b[l] & MSQ_P) >> TRA_P ;
	assert(t[l] != 0x0) ;
	_mult_sp[base_p](t[l], base.b[l]) ;
    }
}
}
