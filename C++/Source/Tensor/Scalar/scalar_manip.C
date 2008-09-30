/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *   
 *   Copyright (c) 2000-2001 Philippe Grandclement (Cmp version)
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


char scalar_manip_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.16  2008/09/30 08:35:18  j_novak
 * Correction of forgotten call to coef()
 *
 * Revision 1.15  2008/09/29 13:23:51  j_novak
 * Implementation of the angular mapping associated with an affine
 * mapping. Things must be improved to take into account the domain index.
 *
 * Revision 1.14  2008/09/22 19:08:01  j_novak
 * New methods to deal with boundary conditions
 *
 * Revision 1.13  2007/06/21 20:00:00  k_taniguchi
 * Addition of the method filtre_r (int n, int nz)
 *
 * Revision 1.12  2006/06/28 07:46:41  j_novak
 * Better treatment in the case of a domain set to zero.
 *
 * Revision 1.11  2005/09/07 13:39:10  j_novak
 * *** empty log message ***
 *
 * Revision 1.10  2005/09/07 13:10:48  j_novak
 * Added a filter setting to zero all mulitpoles in a given range.
 *
 * Revision 1.9  2004/11/23 12:47:44  f_limousin
 * Add function filtre_tp(int nn, int nz1, int nz2).
 *
 * Revision 1.8  2004/05/07 11:24:54  f_limousin
 * Implement new method filtre_r (int* nn)
 *
 * Revision 1.7  2004/02/27 09:47:26  f_limousin
 * New methods filtre_phi(int) and filtre_theta(int).
 *
 * Revision 1.6  2004/01/23 13:26:28  e_gourgoulhon
 *  Added methods set_inner_boundary and set_outer_boundary.
 *  Methods set_val_inf and set_val_hor, which are particular cases of
 *  the above, have been suppressed.
 *
 * Revision 1.5  2003/11/13 13:43:55  p_grandclement
 * Addition of things needed for Bhole::update_metric (const Etoile_bin&, double, double)
 *
 * Revision 1.4  2003/10/11 14:47:17  e_gourgoulhon
 * Lines 112 and 145 : replaced "0" by "double(0)" in the comparison
 * statement.
 *
 * Revision 1.3  2003/10/10 15:57:29  j_novak
 * Added the state one (ETATUN) to the class Scalar
 *
 * Revision 1.2  2003/10/08 14:24:09  j_novak
 * replaced mult_r_zec with mult_r_ced
 *
 * Revision 1.1  2003/09/25 09:33:36  j_novak
 * Added methods for integral calculation and various manipulations
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>
#include <math.h>

// Lorene
#include "tensor.h"
#include "proto.h"
#include "utilitaires.h"

/*
 * Annule tous les l entre l_min et l_max (compris)
 */

void Scalar::annule_l (int l_min, int l_max, bool ylm_output) {

    assert (etat != ETATNONDEF) ;
    assert (l_min <= l_max) ;
    assert (l_min >= 0) ;
    if (etat == ETATZERO )
	return ;

    if (etat == ETATUN) {
	if (l_min == 0) set_etat_zero() ;
	else return ;
    }

    va.ylm() ;
    Mtbl_cf& m_coef = *va.c_cf ;
    const Base_val& base = va.base ; 
    int l_q, m_q, base_r ;
    for (int lz=0; lz<mp->get_mg()->get_nzone(); lz++)
	if (m_coef(lz).get_etat() != ETATZERO)
	    for (int k=0; k<mp->get_mg()->get_np(lz)+1; k++)
		for (int j=0; j<mp->get_mg()->get_nt(lz); j++)
		    for (int i=0; i<mp->get_mg()->get_nr(lz); i++) {
			base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
			if ((l_min <= l_q) && (l_q<= l_max)) 
			    m_coef.set(lz, k, j, i) = 0 ;
		    }
    if (va.c != 0x0) {
	delete va.c ;
	va.c = 0x0 ;
    }
    if (!ylm_output)
	va.ylm_i() ;
    
    return ;
}

/*
 * Annule les n derniers coefficients en r dans la derniere zone
 */
 
void Scalar::filtre (int n) {
    
    assert (etat != ETATNONDEF) ;
    if ( (etat == ETATZERO) || (etat == ETATUN) )
	return ;
    
    int nz = mp->get_mg()->get_nzone() ;
    int np = mp->get_mg()->get_np(nz-1) ;
    int nt = mp->get_mg()->get_nt(nz-1) ;
    int nr = mp->get_mg()->get_nr(nz-1) ;
    
    del_deriv() ;
    
    va.coef() ;
    va.set_etat_cf_qcq() ;

    
    for (int k=0 ; k<np+1 ; k++)
	if (k!=1)
	    for (int j=0 ; j<nt ; j++)
		for (int i=nr-1 ; i>nr-1-n ; i--)
		    va.c_cf->set(nz-1, k, j, i) = 0 ;
}


/*
 * Annule les n derniers coefficients en r dans toutes les zones
 */
 
void Scalar::filtre_r (int* nn) {
    assert (etat != ETATNONDEF) ;
    if ( (etat == ETATZERO) || (etat == ETATUN) )
	return ;
    
    del_deriv() ;
    
    va.coef() ;
    va.set_etat_cf_qcq() ;
    int nz = mp->get_mg()->get_nzone() ;
    int* nr = new int[nz];
    int* nt = new int[nz];
    int* np = new int[nz];
    for (int l=0; l<=nz-1; l++) {
	nr[l] = mp->get_mg()->get_nr(l) ; 
	nt[l] = mp->get_mg()->get_nt(l) ; 
	np[l] = mp->get_mg()->get_np(l) ; 
    }
 
    for (int l=0; l<=nz-1; l++)
	for (int k=0 ; k<np[l]+1 ; k++)
	    if (k!=1)
		for (int j=0 ; j<nt[l] ; j++)
		    for (int i=nr[l]-1; i>nr[l]-1-nn[l] ; i--)
			va.c_cf->set(l, k, j, i) = 0 ;

	if (va.c != 0x0) {
	    delete va.c ;
	    va.c = 0x0 ;
	}
    
}


/*
 * Annule les n derniers coefficients en r dans zone nz
 */

void Scalar::filtre_r (int n, int nz) {
    assert (etat != ETATNONDEF) ;
    if ( (etat == ETATZERO) || (etat == ETATUN) )
	return ;
    
    del_deriv() ;
    
    va.coef() ;
    va.set_etat_cf_qcq() ;
    int nr = mp->get_mg()->get_nr(nz) ; 
    int nt = mp->get_mg()->get_nt(nz) ; 
    int np = mp->get_mg()->get_np(nz) ; 

    for (int k=0 ; k<np+1 ; k++)
        if (k!=1)
	    for (int j=0 ; j<nt ; j++)
	        for (int i=nr-1; i>nr-1-n ; i--)
		    va.c_cf->set(nz, k, j, i) = 0 ;

    if (va.c != 0x0) {
        delete va.c ;
	va.c = 0x0 ;
    }

}


/*
 * Annule les n derniers coefficients en phi dans zone nz
 */
 
void Scalar::filtre_phi (int n, int nz) {
    assert (etat != ETATNONDEF) ;
    if ( (etat == ETATZERO) || (etat == ETATUN) )
	return ;
    
    del_deriv() ;
    
    va.coef() ;
    va.set_etat_cf_qcq() ;
    int np = mp->get_mg()->get_np(nz) ;
    int nt = mp->get_mg()->get_nt(nz) ;
    int nr = mp->get_mg()->get_nr(nz) ;
    
    for (int k=np+1-n ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++)
	    for (int i=0 ; i<nr ; i++)
		va.c_cf->set(nz, k, j, i) = 0 ;

}


void Scalar::filtre_tp(int nn, int nz1, int nz2) {

    va.filtre_tp(nn, nz1, nz2) ;

}


  /* Sets the value of the {\tt Scalar} at the inner boundary of a given 
   * domain. 
   * @param l [input] domain index
   * @param x [input] (constant) value at the inner boundary of domain no. {\tt l}
   */

void Scalar::set_inner_boundary(int l0, double x0) {
    
    assert (etat != ETATNONDEF) ;
    if (etat == ETATZERO) {
	    if (x0 == double(0)) return ;
	    else annule_hard() ;
    }
    
    if (etat == ETATUN) {
	    if (x0 == double(1)) return ;
	    else etat = ETATQCQ ;
    }
    
    del_deriv() ;
    
    int nt = mp->get_mg()->get_nt(l0) ;
    int np = mp->get_mg()->get_np(l0) ;
    
    va.coef_i() ;
    va.set_etat_c_qcq() ;
    
    for (int k=0 ; k<np ; k++)
	for (int j=0 ; j<nt ; j++)
	    va.set(l0, k, j, 0) = x0 ;
}

  /* Sets the value of the {\tt Scalar} at the outer boundary of a given 
   * domain. 
   * @param l [input] domain index
   * @param x [input] (constant) value at the outer boundary of domain no. {\tt l}
   */

void Scalar::set_outer_boundary(int l0, double x0) {
    
    assert (etat != ETATNONDEF) ;
    if (etat == ETATZERO) {
	    if (x0 == double(0)) return ;
	    else annule_hard() ;
    }
    
    if (etat == ETATUN) {
	    if (x0 == double(1)) return ;
	    else etat = ETATQCQ ;
    }
    
    del_deriv() ;
    
    int nrm1 = mp->get_mg()->get_nr(l0) - 1 ;
    int nt = mp->get_mg()->get_nt(l0) ;
    int np = mp->get_mg()->get_np(l0) ;
    
    va.coef_i() ;
    va.set_etat_c_qcq() ;
    
    for (int k=0 ; k<np ; k++)
	for (int j=0 ; j<nt ; j++)
	    va.set(l0, k, j, nrm1) = x0 ;
}

/*
 * Permet de fixer la decroissance du cmp a l infini en viurant les 
 * termes en 1/r^n
 */
void Scalar::fixe_decroissance (int puis) {
    
    if (puis<dzpuis)
	return ;
    else {
	
	int nbre = puis-dzpuis ;
	
	// le confort avant tout ! (c'est bien le confort ...)
	int nz = mp->get_mg()->get_nzone() ;
	int np = mp->get_mg()->get_np(nz-1) ;
	int nt = mp->get_mg()->get_nt(nz-1) ;
	int nr = mp->get_mg()->get_nr(nz-1) ;
	
	const Map_af* map  = dynamic_cast<const Map_af*>(mp) ;
	if (map == 0x0) {
	    cout << "Le mapping doit etre affine" << endl ;
	    abort() ;
	}
	
	double alpha = map->get_alpha()[nz-1] ;
	
	Scalar courant (*this) ;
	
	va.coef() ;
	va.set_etat_cf_qcq() ;
	
	for (int conte=0 ; conte<nbre ; conte++) {
	    
	    int base_r = courant.va.base.get_base_r(nz-1) ;
	    
	    courant.va.coef() ;
	    
	    // On calcul les coefficients de 1/r^conte
	    double* coloc = new double [nr] ;
	    int * deg = new int[3] ;
	    deg[0] = 1 ; 
	    deg[1] = 1 ;
	    deg[2] = nr ;
		    
	    for (int i=0 ; i<nr ; i++)
		coloc[i] =pow(alpha, double(conte))*
		    pow(-1-cos(M_PI*i/(nr-1)), double(conte)) ;
		    
	    cfrcheb(deg, deg, coloc, deg, coloc) ;
	    
	    for (int k=0 ; k<np+1 ; k++)
		if (k != 1)
		for (int j=0 ; j<nt ; j++) {
		    
		    // On doit determiner le coefficient du truc courant :
		    double* coef = new double [nr] ;
		    double* auxi = new double[1] ;
		    for (int i=0 ; i<nr ; i++)
			coef[i] = (*courant.va.c_cf)(nz-1, k, j, i) ;
		    switch (base_r) {
			case R_CHEBU :
			som_r_chebu (coef, nr, 1, 1, 1, auxi) ;
			break ;
		    default :
			som_r_pas_prevu (coef, nr, 1, 1, 1, auxi) ;
			break ;
		    }
		    
		    // On modifie le cmp courant :
		    courant.va.coef() ;
		    courant.va.set_etat_cf_qcq() ;
		    courant.va.c_cf->set(nz-1, k, j, 0) -= *auxi ;  
			
		    for (int i=0 ; i<nr ; i++)
		    	this->va.c_cf->set(nz-1, k, j, i) -= *auxi * coloc[i] ;

			  
		    delete [] coef ;
		    delete [] auxi ;
		}
	    delete [] coloc ;
	    delete [] deg ;
	    
	    courant.mult_r_ced() ;
	}
    }
}

Tbl Scalar::tbl_out_bound(int l_zone, bool output_ylm) {
    va.coef() ;
    if (output_ylm) va.ylm() ;

    int np = mp->get_mg()->get_np(l_zone) ;
    int nt = mp->get_mg()->get_nt(l_zone) ;
    Tbl resu(np+2, nt) ;
    if (etat == ETATZERO) resu.set_etat_zero() ;
    else {
	assert(etat == ETATQCQ) ;
	resu.set_etat_qcq() ;
	for (int k=0; k<np+2; k++)
	    for (int j=0; j<nt; j++)
		resu.set(k, j) = va.c_cf->val_out_bound_jk(l_zone, j, k) ;
    }
    return resu ;
}

Tbl Scalar::tbl_in_bound(int l_zone, bool output_ylm) {
    assert(mp->get_mg()->get_type_r(l_zone) != RARE) ;
    va.coef() ;
    if (output_ylm) va.ylm() ;

    int np = mp->get_mg()->get_np(l_zone) ;
    int nt = mp->get_mg()->get_nt(l_zone) ;
    Tbl resu(np+2, nt) ;
    if (etat == ETATZERO) resu.set_etat_zero() ;
    else {
	assert(etat == ETATQCQ) ;
	resu.set_etat_qcq() ;
	for (int k=0; k<np+2; k++)
	    for (int j=0; j<nt; j++)
		resu.set(k, j) = va.c_cf->val_in_bound_jk(l_zone, j, k) ;
    }
    return resu ;
}

Scalar Scalar::scalar_out_bound(int l_zone, bool output_ylm) {
    va.coef() ;
    if (output_ylm) va.ylm() ;

    Scalar resu(mp->mp_angu(l_zone)) ;
    if (etat == ETATZERO) resu.set_etat_zero() ;
    else {
	assert(etat == ETATQCQ) ;
	resu.allocate_all() ;
	int np = mp->get_mg()->get_np(l_zone) ;
	int nt = mp->get_mg()->get_nt(l_zone) ;
	for (int k=0; k<np; k++)
	    for (int j=0; j<nt; j++)
		resu.set_grid_point(0, k, j, 0) 
		    = va.c_cf->val_out_bound_jk(l_zone, j, k) ;
    }
    resu.std_spectral_base() ;
    Base_val base = resu.get_spectral_base() ;
    base.set_base_t(va.base.get_base_t(l_zone)) ;
    resu.set_spectral_base(base) ;
    return resu ;
}
