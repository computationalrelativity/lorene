/*
 * Computaiton of d/dtheta
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
 * Revision 1.6  2016/12/05 16:18:20  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:53:49  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:23  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2009/10/08 16:23:14  j_novak
 * Addition of new bases T_COS and T_SIN.
 *
 * Revision 1.2  2004/11/23 15:17:19  m_forot
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.8  1999/11/30  12:44:02  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.7  1999/11/23  16:16:55  eric
 * Reorganisation du calcul dans le cas ETATZERO.
 *
 * Revision 2.6  1999/11/19  09:30:55  eric
 * La valeur de retour est desormais const Valeur &.
 *
 * Revision 2.5  1999/10/28  08:00:06  eric
 * Modif commentaires.
 *
 * Revision 2.4  1999/10/18  13:41:37  eric
 * Suppression de l'argument base dans les routines de derivation des mtbl_cf.
 *
 * Revision 2.3  1999/09/07  15:41:40  phil
 * correction Valeur Valeur::dsdx
 *
 * Revision 2.2  1999/03/01  15:09:32  eric
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/23  10:52:11  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/22  15:39:23  hyc
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
void _dsdtet_pas_prevu(Tbl *, int &) ;
void _dsdtet_t_cos(Tbl *, int &) ;
void _dsdtet_t_sin(Tbl *, int &) ;
void _dsdtet_t_cos_p(Tbl *, int &) ;
void _dsdtet_t_sin_p(Tbl *, int &) ;
void _dsdtet_t_sin_i(Tbl *, int &) ;
void _dsdtet_t_cos_i(Tbl *, int &) ;
void _dsdtet_t_cossin_cp(Tbl *, int &) ;
void _dsdtet_t_cossin_sp(Tbl *, int &) ;
void _dsdtet_t_cossin_ci(Tbl *, int &) ;
void _dsdtet_t_cossin_si(Tbl *, int &) ;
void _dsdtet_t_cossin_s(Tbl *, int &) ;
void _dsdtet_t_cossin_c(Tbl *, int &) ;


// Version membre d'un Valeur
// --------------------------

const Valeur& Valeur::dsdt() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_dsdt != 0x0) {
	return *p_dsdt ;
    }
    
    // ... si, il faut bosser

    p_dsdt = new Valeur(mg) ;

    if (etat == ETATZERO) {
	p_dsdt->set_etat_zero() ; 
    }
    else {
	assert(etat == ETATQCQ) ; 
	p_dsdt->set_etat_cf_qcq() ;
	Mtbl_cf* cfp = p_dsdt->c_cf ; // Pointeur sur le Mtbl_cf qui vient d'etre
					// cree par le set_etat_cf_qcq()

	// Initialisation de *cfp : recopie des coef. de la fonction
	if (c_cf == 0x0) {
	    coef() ;
	}
	*cfp = *c_cf ;	
 
	cfp->dsdt() ;	// calcul 
    
	p_dsdt->base = cfp->base ; // On remonte la base de sortie au niveau Valeur
    }

    // Termine
    return *p_dsdt ;
}

// Version membre d'un Mtbl_cf
// ---------------------------
void Mtbl_cf::dsdt() {

// Routines de derivation
static void (*_dsdtet[MAX_BASE])(Tbl *, int &) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    _dsdtet[i] = _dsdtet_pas_prevu ;
	}
	// Les routines existantes
	_dsdtet[T_COS >> TRA_T] = _dsdtet_t_cos ;
	_dsdtet[T_SIN >> TRA_T] = _dsdtet_t_sin ;
	_dsdtet[T_COS_P >> TRA_T] = _dsdtet_t_cos_p ;
	_dsdtet[T_SIN_P >> TRA_T] = _dsdtet_t_sin_p ;
	_dsdtet[T_SIN_I >> TRA_T] = _dsdtet_t_sin_i ;
	_dsdtet[T_COS_I >> TRA_T] = _dsdtet_t_cos_i ;
	_dsdtet[T_COSSIN_CP >> TRA_T] = _dsdtet_t_cossin_cp ;
	_dsdtet[T_COSSIN_SP >> TRA_T] = _dsdtet_t_cossin_sp ;
	_dsdtet[T_COSSIN_CI >> TRA_T] = _dsdtet_t_cossin_ci ;
	_dsdtet[T_COSSIN_SI >> TRA_T] = _dsdtet_t_cossin_si ;
	_dsdtet[T_COSSIN_C >> TRA_T] = _dsdtet_t_cossin_c ;
	_dsdtet[T_COSSIN_S >> TRA_T] = _dsdtet_t_cossin_s ;
    }

    //- Debut de la routine -

    // Protection
    assert(etat == ETATQCQ) ;
    
    // Boucle sur les zones
    for (int l=0 ; l<nzone ; l++) {
	int base_t = (base.b[l] & MSQ_T) >> TRA_T ;
	assert(t[l] != 0x0) ;
	_dsdtet[base_t](t[l], base.b[l]) ;
    }
}
}
