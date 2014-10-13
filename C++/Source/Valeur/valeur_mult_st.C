/*
 * Multiplication by sin(theta)
 *
 * for:
 *   - Valeur
 *   - Mtbl_cf
 */

/*
 *   Copyright (c) 1999-2001 Jerome Novak
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


char valeur_mult_st_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:53:50  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:24  j_novak
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
 * Revision 1.4  1999/11/30  12:44:54  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 1.3  1999/11/24  09:33:41  eric
 * Modif commentaires.
 *
 * Revision 1.2  1999/11/23  16:17:55  eric
 * Reorganisation du calcul dans le cas ETATZERO.
 *
 * Revision 1.1  1999/11/23  14:31:41  novak
 * Initial revision
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
void _mult_st_pas_prevu (Tbl*, int&) ;
void _mult_st_t_cos (Tbl*, int&) ;
void _mult_st_t_sin (Tbl*, int&) ;
void _mult_st_t_sin_p (Tbl*, int&) ;
void _mult_st_t_sin_i (Tbl*, int&) ;
void _mult_st_t_cos_i (Tbl*, int&) ;
void _mult_st_t_cos_p (Tbl*, int&) ;
void _mult_st_t_cossin_si (Tbl*, int&) ;
void _mult_st_t_cossin_ci (Tbl*, int&) ;
void _mult_st_t_cossin_cp (Tbl*, int&) ;
void _mult_st_t_cossin_sp (Tbl*, int&) ;
void _mult_st_t_cossin_c (Tbl*, int&) ;
void _mult_st_t_cossin_s (Tbl*, int&) ;

// Version membre d'un Valeur
// --------------------------

const Valeur& Valeur::mult_st() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_mult_st != 0x0) {
	return *p_mult_st ;
    }
    
    // ... si, il faut bosser

    p_mult_st = new Valeur(mg) ;

    if (etat == ETATZERO) {
	p_mult_st->set_etat_zero() ; 
    }
    else {
	assert(etat == ETATQCQ) ; 
	p_mult_st->set_etat_cf_qcq() ;
	Mtbl_cf* cfp = p_mult_st->c_cf ; // Pointeur sur le Mtbl_cf qui vient d'etre
					// cree par le set_etat_cf_qcq()

	// Initialisation de *cfp : recopie des coef. de la fonction
	if (c_cf == 0x0) {
	    coef() ;
	}
	*cfp = *c_cf ;	
 
	cfp->mult_st() ;	// calcul 
    
	p_mult_st->base = cfp->base ; // On remonte la base de sortie au niveau Valeur
    }

    // Termine
    return *p_mult_st ;
}



// Version membre d'un Mtbl_cf
// ---------------------------

void Mtbl_cf::mult_st() {

// Routines de derivation
static void (*_mult_st[MAX_BASE])(Tbl *, int &) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    _mult_st[i] = _mult_st_pas_prevu ;
	}
	// Les routines existantes
	_mult_st[T_COS >> TRA_T] = _mult_st_t_cos ;
	_mult_st[T_SIN >> TRA_T] = _mult_st_t_sin ;
	_mult_st[T_COS_P >> TRA_T] = _mult_st_t_cos_p ;
	_mult_st[T_COS_I >> TRA_T] = _mult_st_t_cos_i ;
	_mult_st[T_SIN_P >> TRA_T] = _mult_st_t_sin_p ;
	_mult_st[T_SIN_I >> TRA_T] = _mult_st_t_sin_i ;
	_mult_st[T_COSSIN_SI >> TRA_T] = _mult_st_t_cossin_si ;
	_mult_st[T_COSSIN_CI >> TRA_T] = _mult_st_t_cossin_ci ;
	_mult_st[T_COSSIN_CP >> TRA_T] = _mult_st_t_cossin_cp ;
	_mult_st[T_COSSIN_SP >> TRA_T] = _mult_st_t_cossin_sp ;
	_mult_st[T_COSSIN_C >> TRA_T] = _mult_st_t_cossin_c ;
	_mult_st[T_COSSIN_S >> TRA_T] = _mult_st_t_cossin_s ;
    }

    // Debut de la routine 

    // Protection
    assert(etat == ETATQCQ) ;
    
    // Boucle sur les zones
    int base_t ;
    for (int l=0 ; l<nzone ; l++) {
	base_t = (base.b[l] & MSQ_T) >> TRA_T ;
	assert(t[l] != 0x0) ;
	_mult_st[base_t](t[l], base.b[l]) ;
    }
}
}
