/*
 * Computation of d/dphi
 *
 * for:
 *   - Valeur
 *   - Mtbl_cf
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


char valeur_dsdp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2012/01/17 15:08:16  j_penner
 * using MAX_BASE_2 for the phi coordinate
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.10  2000/10/04  12:50:53  eric
 * Ajout de la base P_COSSIN_I.
 *
 * Revision 2.9  1999/11/30  12:43:52  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.8  1999/11/23  16:16:44  eric
 * Reorganisation du calcul dans le cas ETATZERO.
 *
 * Revision 2.7  1999/11/19  09:30:46  eric
 * La valeur de retour est desormais const Valeur &.
 *
 * Revision 2.6  1999/10/28  07:59:37  eric
 * Modif commentaires.
 *
 * Revision 2.5  1999/10/18  13:41:24  eric
 * Suppression de l'argument base dans les routines de derivation des mtbl_cf.
 *
 * Revision 2.4  1999/09/07  15:46:22  phil
 * correction des bases
 *
 * Revision 2.3  1999/03/01  15:09:01  eric
 * *** empty log message ***
 *
 * Revision 2.2  1999/02/23  11:26:27  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <assert.h>

// Headers Lorene
#include "mtbl_cf.h"
#include "valeur.h"

// Prototypage 
void _dsdphi_pas_prevu(Tbl *, int &) ;
void _dsdphi_p_cossin(Tbl *, int &) ;
void _dsdphi_p_cossin_p(Tbl *, int &) ;
void _dsdphi_p_cossin_i(Tbl *, int &) ;

// Version membre d'un Valeur
// --------------------------

const Valeur& Valeur::dsdp() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_dsdp != 0x0) {
	return *p_dsdp ;
    }
    
    // ... si, il faut bosser

    p_dsdp = new Valeur(mg) ;

    if (etat == ETATZERO) {
	p_dsdp->set_etat_zero() ; 
    }
    else {
	assert(etat == ETATQCQ) ; 
	p_dsdp->set_etat_cf_qcq() ;
	Mtbl_cf* cfp = p_dsdp->c_cf ; // Pointeur sur le Mtbl_cf qui vient d'etre
					// cree par le set_etat_cf_qcq()

	// Initialisation de *cfp : recopie des coef. de la fonction
	if (c_cf == 0x0) {
	    coef() ;
	}
	*cfp = *c_cf ;	
 
	cfp->dsdp() ;	// calcul 
    
	p_dsdp->base = cfp->base ; // On remonte la base de sortie au niveau Valeur
    }

    // Termine
    return *p_dsdp ;
}

// Version membre d'un Mtbl_cf
// ---------------------------

void Mtbl_cf::dsdp() {

// Routines de derivation
static void (*_dsdphi[MAX_BASE_2])(Tbl *, int &) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE_2 ; i++) {
	    _dsdphi[i] = _dsdphi_pas_prevu ;
	}
	// Les routines existantes
	    _dsdphi[P_COSSIN >> TRA_P] = _dsdphi_p_cossin ;
	    _dsdphi[P_COSSIN_P >> TRA_P] = _dsdphi_p_cossin_p ;
	    _dsdphi[P_COSSIN_I >> TRA_P] = _dsdphi_p_cossin_i ;
    }

    //- Debut de la routine -

    // Protection
    assert(etat == ETATQCQ) ;
    
    // Boucle sur les zones
    for (int l=0 ; l<nzone ; l++) {
	int base_p = (base.b[l] & MSQ_P) >> TRA_P ;
	assert(t[l] != 0x0) ;
	_dsdphi[base_p](t[l], base.b[l]) ;
    }
}
