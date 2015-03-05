/*
 * Computation of 1/x^2
 *
 * for:
 *   - Valeur
 *   - Mtbl_cf
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char valeur_sx2_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2015/03/05 08:49:33  j_novak
 * Implemented operators with Legendre bases.
 *
 * Revision 1.4  2014/10/13 08:53:51  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:24  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2004/11/23 15:17:20  m_forot
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.8  1999/11/30  12:46:01  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.7  1999/11/23  16:19:22  eric
 * Reorganisation du calcul dans le cas ETATZERO.
 *
 * Revision 2.6  1999/11/19  09:32:02  eric
 * La valeur de retour est desormais const Valeur &.
 *
 * Revision 2.5  1999/10/28  08:00:57  eric
 * Modif commentaires.
 *
 * Revision 2.4  1999/10/18  13:43:07  eric
 * Suppression de l'argument base dans les routines de derivation des mtbl_cf.
 *
 * Revision 2.3  1999/04/28  10:10:59  phil
 * *** empty log message ***
 *
 * Revision 2.2  1999/04/26  17:15:57  phil
 * *** empty log message ***
 *
 * Revision 2.1  1999/04/26  17:14:20  phil
 * *** empty log message ***
 *
 * Revision 2.0  1999/04/26  17:13:09  phil
 * *** empty log message ***
 *
 * Revision 2.0  1999/04/26  14:59:31  phil
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

// Local prototypes: 
namespace Lorene {
void _sx2_pas_prevu(Tbl *, int &) ;
void _sx2_r_chebu(Tbl *, int &) ;
void _sx2_r_chebp(Tbl *, int &) ;
void _sx2_r_chebi(Tbl *, int &) ;
void _sx2_r_chebpim_p(Tbl *, int &) ;
void _sx2_r_chebpim_i(Tbl *, int &) ;
void _sx2_identite (Tbl *, int &) ;
void _sx2_r_chebpi_p(Tbl *, int &) ;
void _sx2_r_chebpi_i(Tbl *, int &) ;
void _sx2_r_legp(Tbl *, int &) ;
void _sx2_r_legi(Tbl *, int &) ;

// Version membre d'un Valeur
// --------------------------

const Valeur& Valeur::sx2() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_sx2 != 0x0) {
	return *p_sx2 ;
    }
    
    // ... si, il faut bosser

    p_sx2 = new Valeur(mg) ;

    if (etat == ETATZERO) {
	p_sx2->set_etat_zero() ; 
    }
    else {
	assert(etat == ETATQCQ) ; 
	p_sx2->set_etat_cf_qcq() ;
	Mtbl_cf* cfp = p_sx2->c_cf ; // Pointeur sur le Mtbl_cf qui vient d'etre
					// cree par le set_etat_cf_qcq()

	// Initialisation de *cfp : recopie des coef. de la fonction
	if (c_cf == 0x0) {
	    coef() ;
	}
	*cfp = *c_cf ;	
 
	cfp->sx2() ;	// calcul 
    
	p_sx2->base = cfp->base ; // On remonte la base de sortie au niveau Valeur
    }

    // Termine
    return *p_sx2 ;
}



// Version membre d'un Mtbl_cf
// ---------------------------
\
void Mtbl_cf::sx2() {

// Routines de derivation
static void (*_sx2[MAX_BASE])(Tbl *, int &) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    _sx2[i] = _sx2_pas_prevu ;
	}
	// Les routines existantes
	_sx2[R_CHEB >> TRA_R] = _sx2_identite ;
	_sx2[R_CHEBU >> TRA_R] = _sx2_r_chebu ;
	_sx2[R_CHEBP >> TRA_R] = _sx2_r_chebp ;
	_sx2[R_CHEBI >> TRA_R] = _sx2_r_chebi ;
	_sx2[R_CHEBPIM_P >> TRA_R] = _sx2_r_chebpim_p ;
	_sx2[R_CHEBPIM_I >> TRA_R] = _sx2_r_chebpim_i ;
	_sx2[R_CHEBPI_P >> TRA_R] = _sx2_r_chebpi_p ;
	_sx2[R_CHEBPI_I >> TRA_R] = _sx2_r_chebpi_i ;
	_sx2[R_LEG >> TRA_R] = _sx2_identite ;
	_sx2[R_LEGP >> TRA_R] = _sx2_r_legp ;
	_sx2[R_LEGI >> TRA_R] = _sx2_r_legi ;
    }

    //- Debut de la routine -

    // Protection
    assert(etat == ETATQCQ) ;
    
    // Boucle sur les zones
    for (int l=0 ; l<nzone ; l++) {
	int base_r = (base.b[l] & MSQ_R) >> TRA_R ;
	assert(t[l] != 0x0) ;
	_sx2[base_r](t[l], base.b[l]) ;
    }
}
}
