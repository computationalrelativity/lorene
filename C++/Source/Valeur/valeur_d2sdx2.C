/*
 * Computation of d^2/dx^2
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


char valeur_d2sdx2_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2008/08/27 08:52:55  jl_cornou
 * Added Jacobi(0,2) polynomials case
 *
 * Revision 1.2  2004/11/23 15:17:19  m_forot
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.9  2000/08/16  10:46:03  eric
 * Suppression de Mtbl_cf::dzpuis.
 *
 * Revision 2.8  1999/11/30  12:43:41  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.7  1999/11/23  16:16:32  eric
 * Reorganisation du calcul dans le cas ETATZERO.
 *
 * Revision 2.6  1999/11/19  09:30:38  eric
 * La valeur de retour est desormais const Valeur &.
 *
 * Revision 2.5  1999/10/28  07:59:29  eric
 * Modif commentaires.
 *
 * Revision 2.4  1999/10/18  13:41:11  eric
 * Suppression de l'argument base dans les routines de derivation des mtbl_cf.
 *
 * Revision 2.3  1999/09/07  15:43:08  phil
 * coreection bases
 *
 * Revision 2.2  1999/03/01  15:08:40  eric
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/23  10:39:02  hyc
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
void _d2sdx2_pas_prevu(Tbl *, int &) ;
void _d2sdx2_r_cheb(Tbl *, int &) ;
void _d2sdx2_r_chebu_0(Tbl *, int &) ;
void _d2sdx2_r_chebp(Tbl *, int &) ;
void _d2sdx2_r_chebi(Tbl *, int &) ;
void _d2sdx2_r_chebpim_p(Tbl *, int &) ;
void _d2sdx2_r_chebpim_i(Tbl *, int &) ;
void _d2sdx2_r_chebpi_p(Tbl *, int &) ;
void _d2sdx2_r_chebpi_i(Tbl *, int &) ;
void _d2sdx2_r_jaco02(Tbl *, int &) ;

// Version membre d'un Valeur
// --------------------------

const Valeur& Valeur::d2sdx2() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_d2sdx2 != 0x0) {
	return *p_d2sdx2 ;
    }
    
    // ... si, il faut bosser

    p_d2sdx2 = new Valeur(mg) ;

    if (etat == ETATZERO) {
	p_d2sdx2->set_etat_zero() ; 
    }
    else {
	assert(etat == ETATQCQ) ; 
	p_d2sdx2->set_etat_cf_qcq() ;
	Mtbl_cf* cfp = p_d2sdx2->c_cf ; // Pointeur sur le Mtbl_cf qui vient d'etre
					// cree par le set_etat_cf_qcq()

	// Initialisation de *cfp : recopie des coef. de la fonction
	if (c_cf == 0x0) {
	    coef() ;
	}
	*cfp = *c_cf ;	
 
	cfp->d2sdx2() ;	// calcul 
    
	p_d2sdx2->base = cfp->base ; // On remonte la base de sortie au niveau Valeur
    }

    // Termine
    return *p_d2sdx2 ;
}

// Version membre d'un Mtbl_cf
// ---------------------------

void Mtbl_cf::d2sdx2() {

// Routines de derivation
static void (*_d2sdx2[MAX_BASE])(Tbl *, int &) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    _d2sdx2[i] = _d2sdx2_pas_prevu ;
	}
	// Les routines existantes
	_d2sdx2[R_CHEB >> TRA_R] = _d2sdx2_r_cheb ;
	_d2sdx2[R_CHEBU >> TRA_R] = _d2sdx2_r_chebu_0 ;
	_d2sdx2[R_CHEBP >> TRA_R] = _d2sdx2_r_chebp ;
	_d2sdx2[R_CHEBI >> TRA_R] = _d2sdx2_r_chebi ;
	_d2sdx2[R_CHEBPIM_P >> TRA_R] = _d2sdx2_r_chebpim_p ;
	_d2sdx2[R_CHEBPIM_I >> TRA_R] = _d2sdx2_r_chebpim_i ;
	_d2sdx2[R_CHEBPI_P >> TRA_R] = _d2sdx2_r_chebpi_p ;
	_d2sdx2[R_CHEBPI_I >> TRA_R] = _d2sdx2_r_chebpi_i ;
	_d2sdx2[R_JACO02 >> TRA_R] = _d2sdx2_r_jaco02 ;
    }

    //- Debut de la routine -

    // Protection
    assert(etat == ETATQCQ) ;

    // Boucle sur les zones
    for (int l=0 ; l<nzone ; l++) {
	int base_r = (base.b[l] & MSQ_R) >> TRA_R ;
	assert(t[l] != 0x0) ;
	_d2sdx2[base_r](t[l], base.b[l]) ;
    }
}
