/*
 * Computation of d/dx
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


char valeur_dsdx_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.7  1999/11/30  12:44:13  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.6  1999/11/23  16:17:06  eric
 * Reorganisation du calcul dans le cas ETATZERO.
 *
 * Revision 2.5  1999/11/19  09:31:06  eric
 * La valeur de retour est desormais const Valeur &.
 *
 * Revision 2.4  1999/10/28  08:00:20  eric
 * Modif commentaires.
 *
 * Revision 2.3  1999/10/18  13:41:46  eric
 * Suppression de l'argument base dans les routines de derivation des mtbl_cf.
 *
 * Revision 2.2  1999/04/09  14:03:58  phil
 * Correction erreur base dans Valeur Valeur::dsdx.
 *
 * Revision 2.1  1999/03/01  14:56:40  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/22  15:39:13  hyc
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

// Prototypage temporaire
void _dsdx_pas_prevu(Tbl *, int &) ;
void _dsdx_r_cheb(Tbl *, int &) ;
void _dsdx_r_chebu(Tbl *, int &) ;
void _dsdx_r_chebp(Tbl *, int &) ;
void _dsdx_r_chebi(Tbl *, int &) ;
void _dsdx_r_chebpim_p(Tbl *, int &) ;
void _dsdx_r_chebpim_i(Tbl *, int &) ;

// Version membre d'un Valeur
// --------------------------

const Valeur& Valeur::dsdx() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_dsdx != 0x0) {
	return *p_dsdx ;
    }
    
    // ... si, il faut bosser

    p_dsdx = new Valeur(mg) ;

    if (etat == ETATZERO) {
	p_dsdx->set_etat_zero() ; 
    }
    else {
	assert(etat == ETATQCQ) ; 
	p_dsdx->set_etat_cf_qcq() ;
	Mtbl_cf* cfp = p_dsdx->c_cf ; // Pointeur sur le Mtbl_cf qui vient d'etre
					// cree par le set_etat_cf_qcq()

	// Initialisation de *cfp : recopie des coef. de la fonction
	if (c_cf == 0x0) {
	    coef() ;
	}
	*cfp = *c_cf ;	
 
	cfp->dsdx() ;	// calcul 
    
	p_dsdx->base = cfp->base ; // On remonte la base de sortie au niveau Valeur
    }
    
    // Termine
    return *p_dsdx ;
}


// Version membre d'un Mtbl_cf
// ---------------------------

void Mtbl_cf::dsdx() {

// Routines de derivation
static void (*_dsdx[MAX_BASE])(Tbl *, int &) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    _dsdx[i] = _dsdx_pas_prevu ;
	}
	// Les routines existantes
	_dsdx[R_CHEB >> TRA_R] = _dsdx_r_cheb ;
	_dsdx[R_CHEBU >> TRA_R] = _dsdx_r_chebu ;
	_dsdx[R_CHEBP >> TRA_R] = _dsdx_r_chebp ;
	_dsdx[R_CHEBI >> TRA_R] = _dsdx_r_chebi ;
	_dsdx[R_CHEBPIM_P >> TRA_R] = _dsdx_r_chebpim_p ;
	_dsdx[R_CHEBPIM_I >> TRA_R] = _dsdx_r_chebpim_i ;
    }

    //- Debut de la routine -

    // Protection
    assert(etat == ETATQCQ) ;
    
    // Boucle sur les zones
    for (int l=0 ; l<nzone ; l++) {
	int base_r = (base.b[l] & MSQ_R) >> TRA_R ;
	assert(t[l] != 0x0) ;
	_dsdx[base_r](t[l], base.b[l]) ;
    }
}
