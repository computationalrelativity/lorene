/*
 * Computation of d^2/dtheta^2
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


char valeur_d2sdt2_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:53:49  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:23  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2009/10/09 14:01:06  j_novak
 * New bases T_cos and T_SIN.
 *
 * Revision 1.2  2004/11/23 15:17:19  m_forot
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.7  1999/11/30  12:43:25  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.6  1999/11/23  16:16:23  eric
 * Reorganisation du calcul dans le cas ETATZERO.
 *
 * Revision 2.5  1999/11/19  09:30:26  eric
 * La valeur de retour est desormais const Valeur &.
 *
 * Revision 2.4  1999/10/28  07:59:16  eric
 * Modif commentaires.
 *
 * Revision 2.3  1999/10/18  13:41:00  eric
 * Suppression de l'argument base dans les routines de derivation des mtbl_cf.
 *
 * Revision 2.2  1999/09/07  15:45:25  phil
 * correction gestion des bases
 *
 * Revision 2.1  1999/03/01  15:08:18  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/23  11:23:27  hyc
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

// Prototypage temporaire
namespace Lorene {
void _d2sdtet2_pas_prevu(Tbl *, int &) ;
void _d2sdtet2_t_cos(Tbl *, int &) ;
void _d2sdtet2_t_sin(Tbl *, int &) ;
void _d2sdtet2_t_cos_p(Tbl *, int &) ;
void _d2sdtet2_t_sin_p(Tbl *, int &) ;
void _d2sdtet2_t_sin_i(Tbl *, int &) ;
void _d2sdtet2_t_cos_i(Tbl *, int &) ;
void _d2sdtet2_t_cossin_cp(Tbl *, int &) ;
void _d2sdtet2_t_cossin_sp(Tbl *, int &) ;
void _d2sdtet2_t_cossin_si(Tbl *, int &) ;
void _d2sdtet2_t_cossin_c(Tbl *, int &) ;
void _d2sdtet2_t_cossin_s(Tbl *, int &) ;

// Version membre d'un Valeur
// --------------------------

const Valeur& Valeur::d2sdt2() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_d2sdt2 != 0x0) {
	return *p_d2sdt2 ;
    }
    
    // ... si, il faut bosser

    p_d2sdt2 = new Valeur(mg) ;

    if (etat == ETATZERO) {
	p_d2sdt2->set_etat_zero() ; 
    }
    else {
	assert(etat == ETATQCQ) ; 
	p_d2sdt2->set_etat_cf_qcq() ;
	Mtbl_cf* cfp = p_d2sdt2->c_cf ; // Pointeur sur le Mtbl_cf qui vient d'etre
					// cree par le set_etat_cf_qcq()

	// Initialisation de *cfp : recopie des coef. de la fonction
	if (c_cf == 0x0) {
	    coef() ;
	}
	*cfp = *c_cf ;	
 
	cfp->d2sdt2() ;	// calcul 
    
	p_d2sdt2->base = cfp->base ; // On remonte la base de sortie au niveau Valeur
    }

    // Termine
    return *p_d2sdt2 ;
}

// Version membre d'un Mtbl_cf
// ---------------------------

void Mtbl_cf::d2sdt2() {

// Routines de derivation
static void (*_d2sdtet2[MAX_BASE])(Tbl *, int &) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    _d2sdtet2[i] = _d2sdtet2_pas_prevu ;
	}
	// Les routines existantes
	_d2sdtet2[T_COS >> TRA_T] = _d2sdtet2_t_cos ;
	_d2sdtet2[T_SIN >> TRA_T] = _d2sdtet2_t_sin ;
	_d2sdtet2[T_COS_P >> TRA_T] = _d2sdtet2_t_cos_p ;
	_d2sdtet2[T_SIN_P >> TRA_T] = _d2sdtet2_t_sin_p ;
	_d2sdtet2[T_SIN_I >> TRA_T] = _d2sdtet2_t_sin_i ;
	_d2sdtet2[T_COS_I >> TRA_T] = _d2sdtet2_t_cos_i ;
	_d2sdtet2[T_COSSIN_CP >> TRA_T] = _d2sdtet2_t_cossin_cp ;
	_d2sdtet2[T_COSSIN_SP >> TRA_T] = _d2sdtet2_t_cossin_sp ;
	_d2sdtet2[T_COSSIN_SI >> TRA_T] = _d2sdtet2_t_cossin_si ;
	_d2sdtet2[T_COSSIN_S >> TRA_T] = _d2sdtet2_t_cossin_s ;
	_d2sdtet2[T_COSSIN_C >> TRA_T] = _d2sdtet2_t_cossin_c ;
    }

    //- Debut de la routine -

    // Protection
    assert(etat == ETATQCQ) ;
    
    // Boucle sur les zones
    for (int l=0 ; l<nzone ; l++) {
	int base_t = (base.b[l] & MSQ_T) >> TRA_T ;
	assert(t[l] != 0x0) ;
	_d2sdtet2[base_t](t[l], base.b[l]) ;
    }
}
}
