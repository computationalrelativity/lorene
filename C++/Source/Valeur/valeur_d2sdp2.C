/*
 * Computation of d^2/dphi^2
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:18:20  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:49  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:23  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.8  1999/11/30  12:43:12  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.7  1999/11/23  16:15:57  eric
 * Reorganisation du calcul dans le cas ETATZERO.
 *
 * Revision 2.6  1999/11/19  09:29:56  eric
 * La valeur de retour est desormais const Valeur &.
 *
 * Revision 2.5  1999/10/28  07:58:56  eric
 * Modif commentaires.
 *
 * Revision 2.4  1999/10/18  13:40:38  eric
 * Suppression de l'argument base dans les routines de derivation des mtbl_cf.
 *
 * Revision 2.3  1999/09/07  15:47:03  phil
 * correction de la gestion des bases
 *
 * Revision 2.2  1999/03/01  15:07:51  eric
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/23  11:43:58  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/23  11:32:52  hyc
 * *** empty log message ***
 *
 * $Header$
 *
 */

// Headers C
#include <cassert>

// Headers Lorene
#include "mtbl_cf.h"
#include "valeur.h"

// Prototypage temporaire
namespace Lorene {
void _d2sdphi2_pas_prevu(Tbl *, int &) ;
void _d2sdphi2_p_cossin(Tbl *, int &) ;

// Version membre d'un Valeur
// --------------------------

const Valeur& Valeur::d2sdp2() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_d2sdp2 != 0x0) {
	return *p_d2sdp2 ;
    }
    
    // ... si, il faut bosser

    p_d2sdp2 = new Valeur(mg) ;

    if (etat == ETATZERO) {
	p_d2sdp2->set_etat_zero() ; 
    }
    else {
	assert(etat == ETATQCQ) ; 
	p_d2sdp2->set_etat_cf_qcq() ;
	Mtbl_cf* cfp = p_d2sdp2->c_cf ; // Pointeur sur le Mtbl_cf qui vient d'etre
					// cree par le set_etat_cf_qcq()
	// Initialisation de *cfp : recopie des coef. de la fonction

	if (c_cf == 0x0) {
	    coef() ;
	}
	*cfp = *c_cf ;	
 
	cfp->d2sdp2() ;	// calcul 
    
	p_d2sdp2->base = cfp->base ; // On remonte la base de sortie au niveau Valeur
    }
    
    // Termine
    return *p_d2sdp2 ;
}

// Version membre d'un Mtbl_cf
// ---------------------------
void Mtbl_cf::d2sdp2() {

// Routines de derivation
static void (*_d2sdphi2[MAX_BASE])(Tbl *, int &) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    _d2sdphi2[i] = _d2sdphi2_pas_prevu ;
	}
	// Les routines existantes
	    _d2sdphi2[P_COSSIN >> TRA_P] = _d2sdphi2_p_cossin ;
    }

    //- Debut de la routine -

    // Protection
    assert(etat == ETATQCQ) ;
    
    // Boucle sur les zones
    for (int l=0 ; l<nzone ; l++) {
	int base_p = (base.b[l] & MSQ_P) >> TRA_P ;
	assert(t[l] != 0x0) ;
	_d2sdphi2[base_p](t[l], base.b[l]) ;
    }
}
}
